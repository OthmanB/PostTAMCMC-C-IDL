@read_ascii_kepler
@derivatives_discret
; CLEANED FOR MANY ECHELLE DIAGRAM FUNCTION THAT ARE NOT USEFUL FOR MAIN SEQUENCE STARS

; This is an echelle diagram where identical frequency are on an horizontal line
pro echellecorot_v4,freq,spec_reg,delta,nbegin,uu,freq_table=freq_t,freq_scd2=freq_t2,smooth_coef=smooth_coef,ps=ps,overplot=overplot,$
	color=color, scd_color=scd_color,show_ylabel=show_ylabel, charsize=charsize, $
	n_number=n_number,print=print, doublecolor=doublecolor, doublepos=doublepos, bar=bar, $
	tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, tag_color=tag_color, $
	extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, int_coef=int_coef, $
	vert_smooth=vert_smooth, recenter=recenter,$
	fbegin=fbegin, fend=fend, set_plot=set_plot ;, dots_extra=dots_extra

A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A) ;, /FILL

printer='' ; relicat of the appourchaux function
index=''

if n_elements(set_plot) eq 0 then set_plot='X' ; Default, we plot on a linux system 

if n_elements(n_number) eq 0 then n_number= 21
if n_elements(smooth_coef) eq 0 then smooth_coef=1
;if n_elements(tag_label) eq 0 then tag_label=''
if n_elements(tag_pos) eq 0 then tag_pos=[5.1, 90]
if n_elements(tag_size) eq 0 then tag_size=1.
if n_elements(int_coef) eq 0 then int_coef=1
if n_elements(vert_smooth) eq 0 then vert_smooth=1
;if n_elements(dots_extra) eq 0 then dots_extra=-1

if n_elements(show_ylabel) eq 0 then show_ylabel=1
if n_elements(charsize) eq 0 then charsize=2.

freq0=freq
spec_reg0=spec_reg
Delf=max(freq0) - min(freq0)
df0=freq[2] - freq[0]
N0=n_elements(freq0)

N1=int_coef*N0
df1=Delf/N1

freq1=findgen(N1)*df1 + min(freq0)
spec_reg1=interpol(spec_reg0, freq0, freq1, /quadratic)

freq2=freq1
spec_reg2=spec_reg1

a=smooth(spec_reg,smooth_coef*int_coef)


if n_elements(extralines) eq 0 then begin
	extralines=-1
	extralines_color=-1
endif
if n_elements(extralines) ne n_elements(pextralines_col) then pextralines_col=strarr(n_elements(extralines)) + 'Black'


n_numb=n_number
nbeg=nbegin

if n_elements(freq_t) ne 0  then if  freq_t[0] ne -1 then begin
if overplot ne 4 AND overplot ne 5 then freq_table=freq_t
if overplot eq 4 then begin ; if overplot = 4 we assume that the table is COMPATIBLE with Draw_BoxAndWiskers_horizontal
	stat_synthese_freq=freq_t
	freq_table=dblarr(n_elements(stat_synthese_freq[0,*,0]),n_elements(stat_synthese_freq[0,0,*]))
	freq_table[*,*]=stat_synthese_freq[3,*,*] ; freq_table contain ONLY THE MEDIAN POSITION
endif
endif

if n_elements(freq_t2) ne 0  then if  freq_t2[0] ne -1 then begin
if overplot ne 4 AND overplot ne 5 then freq_scd2=freq_t2
if overplot eq 4 then begin ; if overplot = 4 we assume that the table is COMPATIBLE with Draw_BoxAndWiskers_horizontal
	stat_synthese_freq=freq_t2
	freq_scd2=dblarr(n_elements(stat_synthese_freq[0,*,0]),n_elements(stat_synthese_freq[0,0,*]))
	freq_scd2[*,*]=stat_synthese_freq[3,*,*] ; freq_table contain ONLY THE MEDIAN POSITION
endif
endif

if n_elements(overplot) eq 0 then overplot=-1

	nn=N_params()
	s={echelleplot,windowi:0.,delta:0.,echelletitle:''}
	if (nn LT N_tags(s)) then begin
		nn=nn-1
		read_structure,s,nn
		names=tag_names(s)
		for i=nn,N_tags(s)-1 do begin
			zzz=names(i)+'=s.'+names(i)
			;print,zzz
			r=execute(names(i)+'=s.'+names(i))
		endfor
	endif

	if (printer ne '') then begin
		set_plot,'ps',/interpolate
		file='corotvg'+index+'.ps'
		device,/landscape,/color,bits=8,filename=file
	endif
	if (n_elements(ps) ne 0) then begin
		set_plot,'ps',/interpolate
		device,/times
	endif

	del=1d*delta

	windowi=del

	window=windowi/2.

	Na=N_elements(a) ; initial size of the power spectrum
	resol=1d*(freq2[10]-freq2[9]) ; intial resolution

	Nz=1d*del/resol+1 ; intial number of bins
	Nz0=Nz ; save the initial number of bins
	if Nz gt 500 then begin
		Nz=500d ; to avoid incredibly long loops
		Dsmooth=ceil(1d*Nz0/Nz) ; new smooth coefficient
		if Dsmooth mod 2 ne 0 then Dsmooth=Dsmooth + 1 ; even numbers only... mayh be can be discarded...
		a=smooth(a,Dsmooth) ; take the power spectrum and smooth it over Dsmooth bins
		f_out=freq2[0]
		s_out=a[0]
		for i=long(0), Na-1 do begin
			if i mod Dsmooth eq 0 then begin
				f_out=[f_out, freq2[i]]
				s_out=[s_out, a[i]]
			endif
		endfor
		freq2=f_out
		a=s_out

		Na=N_elements(a) ; the new smoothed power spectrum size
		resol=1d*(freq2[10]-freq2[9]) ; the new resolution

		Nz=1d*del/resol+1 ; the final number of bins
	endif

	Ny=floor(n_numb*Nz) +1
	b=fltarr(Nz,Ny) ; each slice has Nz subdivisions


	if overplot eq 5 then b2=fltarr(Nz,n_numb*Nz)
	middlefreq=1d*del/2d/resol

	nmax=floor(max(freq2)/del)-1

	nn=min([nmax,nbeg+n_numb-1])

	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
		nu=dblarr(n_elements(freq_table[*,0]),n_elements(freq_table[0,*])*3)
		freq_y=dblarr(n_elements(freq_table[*,0]),n_elements(freq_table[0,*])*3)
	endif
	if n_elements(freq_scd2) ne 0  then if  freq_t2[0] ne -1 then begin
		nu2=dblarr(n_elements(freq_scd2[*,0]),n_elements(freq_scd2[0,*])*3)
		freq_y2=dblarr(n_elements(freq_scd2[*,0]),n_elements(freq_scd2[0,*])*3)
	endif
	if nbeg le 0 then begin
		nbeg=15
		nn=30
	endif

	iend_max=round(1d*(nn+1)*del/resol)
	if iend_max gt n_elements(a) then nn=fix(n_elements(a)*resol/del)-1 ; if we are at upper edge of the spectrum!

	fbegin=freq2[round(1d*(nbeg-1)*del/resol)] ;+ del/2
	fend=freq2[round(1d*(nn)*del/resol)]

	v=min(abs(freq2 - fbegin),i1)
	v=min(abs(freq2 - fend),i2)
	i1=i1-Nz
	i2=i2+Nz
	i2 = i2-1

	FOR j=i1,i2 DO BEGIN
  		;x1=round((freq2(j) mod del)/resol)
		;y1=round((freq2(j)-fbegin)/resol)
		x1=floor((freq2(j) mod del)/resol)
		y1=floor((freq2(j)-fbegin)/resol)

 		FOR k=CEIL(-Nz/2.),FLOOR(Nz/2.) DO BEGIN
   			y2=y1+k
   			IF (y2 gt 0) AND (y2 lt ny) THEN b(x1,y2)= a(j)
 		ENDFOR
	ENDFOR
	b0=b
	if vert_smooth gt 1 then begin
		for i=0, n_elements(b[*, 0])-1 do b[i, *]= smooth(b[i,*], vert_smooth)
	endif

	for i=nbeg,nn do begin
		ibegin=round(1d*i*del/resol)
		iend=round(1d*(i+1)*del/resol)

		if overplot eq 5 then b2(0:(iend-ibegin),i-nbeg)=a2(ibegin:iend)

		;if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
		;	f1=freq[floor(ibegin)] & f2=freq[ceil(ibegin)]
		;	f= (f2 - f1) * (ibegin - floor(ibegin)) + f1
		;	corr_t=corr_t + f/(resol*ibegin)
		;endif
	endfor

; ------------------ TABLE 1 --------------
	if n_elements(freq_table) ne 0 then lmax=n_elements(freq_table[*,0])-1

	kk=0d & phase=0d
	for i=nbeg,nn do begin
		loop=0d
		ibegin=1d*i*del/resol
		iend=1d*(i+1)*del/resol ;-1

		l_tmpQ=dblarr(lmax+1)
		tmpQ=dblarr(lmax+1, 20)-1
		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
			for l=0, lmax do begin
				ta=where(freq_table[l,*] ge 1d*freq2[ibegin] AND freq_table[l,*] le 1d*freq2[iend])
				if ta[0] ne -1 then begin
					tmpQ[l, 0:n_elements(ta)-1]=ta
					l_tmpQ[l]=n_elements(ta)
				endif
			endfor
			dup_max=max(l_tmpQ)
			for j=0, dup_max-1 do begin
				for l=0, lmax do begin
					if tmpQ[l,j] ne -1 then begin
						nu[l, kk]=freq_table[l, tmpQ[l,j]]  mod del
						freq_y[l, kk]=freq_table[l, tmpQ[l,j]]
					endif
				endfor
				kk=kk+1d
			endfor
		endif
	endfor

; ------------------ TABLE 2 --------------
	if n_elements(freq_scd2) ne 0 then lmax=n_elements(freq_scd2[*,0])-1

	kk=0d & phase=0d
	for i=nbeg,nn do begin
		loop=0d
		ibegin=1d*i*del/resol
		iend=1d*(i+1)*del/resol ;-1

		l_tmpQ=dblarr(lmax+1)
		tmpQ=dblarr(lmax+1, 20)-1
		if n_elements(freq_scd2) ne 0  then if  freq_t2[0] ne -1 then begin
			for l=0, lmax do begin
				ta=where(freq_scd2[l,*] ge 1d*freq2[ibegin] AND freq_scd2[l,*] le 1d*freq2[iend])
				if ta[0] ne -1 then begin
					tmpQ[l, 0:n_elements(ta)-1]=ta
					l_tmpQ[l]=n_elements(ta)
				endif
			endfor
			dup_max=max(l_tmpQ)
			for j=0, dup_max-1 do begin
				for l=0, lmax do begin
					if tmpQ[l,j] ne -1 then begin
						nu2[l, kk]=freq_scd2[l, tmpQ[l,j]]  mod del
						freq_y2[l, kk]=freq_scd2[l, tmpQ[l,j]]
					endif
				endfor
				kk=kk+1d
			endfor
		endif
	endfor

	zz=b
	uu=zz(middlefreq-window/resol:middlefreq+window/resol,*)

	if overplot eq 5 then begin
		zz2=congrid(b2,Nz,100)
		uu2=zz2(middlefreq-window/resol:middlefreq+window/resol,*)
	endif

	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin

		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then nu_scale=nu-middlefreq*resol
		pos_dummy=where(nu_scale eq -middlefreq*resol) ; these positions had 0... we need to let the 0
		nu_scale[pos_dummy]=0d

		if overplot eq 4 then begin
			param_stat=dblarr(n_elements(nu_scale[0,*]),n_elements(nu_scale[*,0]),n_elements(stat_synthese_freq[*,0,0]))

			for i=0,n_elements(nu_scale[0,*])-1 do begin
				for j=0,n_elements(stat_synthese_freq[*,0,0])-1 do begin
					for l=0,n_elements(nu_scale[*,0])-1 do $
						param_stat[i,l,j]=nu_scale[l,i]+(stat_synthese_freq[j,l,i]-stat_synthese_freq[3,l,i]) ; using nu_scale, we compute the error boxes
				endfor
			endfor
		endif

		if n_elements(overplot) eq 0 then overplot=2 ; definee the defaut value of overplot : we plot lines + diamonds !
		if n_elements(color) eq 0 then color=['Red','Red','Red','Red','Red']
		if n_elements(color) eq 1 then color=[color,color,color] ; 1 color code by degree

		x0=-delta/2
		no_zero=where(freq_y ne 0)

		if n_elements(doublepos) ne 0 then $
			posy_dble=where(freq_y eq doublepos) ; find the position of the specially tagged color !!

	endif

; ---------------- TABLE 2 --------------
	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, COS(A), SIN(A) ;, /FILL
	if n_elements(freq_scd2) ne 0  then if  freq_t2[0] ne -1 then begin

		if n_elements(freq_scd2) ne 0  then if  freq_t2[0] ne -1 then nu_scale2=nu2-middlefreq*resol
		pos_dummy=where(nu_scale2 eq -middlefreq*resol) ; these positions had 0... we need to let the 0
		nu_scale2[pos_dummy]=0d

		if overplot eq 4 then begin
			param_stat=dblarr(n_elements(nu_scale2[0,*]),n_elements(nu_scale2[*,0]),n_elements(stat_synthese_freq[*,0,0]))

			for i=0,n_elements(nu_scale2[0,*])-1 do begin
				for j=0,n_elements(stat_synthese_freq[*,0,0])-1 do begin
					for l=0,n_elements(nu_scale2[*,0])-1 do $
						param_stat[i,l,j]=nu_scale2[l,i]+(stat_synthese_freq[j,l,i]-stat_synthese_freq[3,l,i]) ; using nu_scale, we compute the error boxes
				endfor
			endfor
		endif

		if n_elements(overplot) eq 0 then overplot=2 ; definee the defaut value of overplot : we plot lines + diamonds !
		if n_elements(color) eq 0 then color=['Red','Red','Red','Red','Red']
		if n_elements(color) eq 1 then color=[color,color,color] ; 1 color code by degree

		x0=-delta/2
		no_zero=where(freq_y2 ne 0)

		if n_elements(doublepos) ne 0 then $
			posy_dble=where(freq_y2 eq doublepos) ; find the position of the specially tagged color !!

	endif

	help,uu
	freq_u=TeXtoIDL('\mu')+'Hz'
	print,max(uu)



if n_elements(ps) eq 0 then begin
	set_plot,set_plot
	device, decomposed=0
endif
if n_elements(ps) eq 1 then if ps eq 0 then begin
	set_plot,set_plot
	device, decomposed=0
endif
!p.font=0 ;& device,set_font='times' ; police de caractere times
;loadct,39  ; met une table de couleur avec 255 = Noir et 1 = Blanc
loadct,39
if n_elements(print) eq 0 then print=0
if print eq 1 then loadct, 0 ; table with a white background

!EDIT_INPUT=50

uu0=uu
freq_y0=freq_y
nu_scale0=nu_scale

if n_elements(freq_scd2) ne 0 then begin
	freq_y02=freq_y2
	nu_scale02=nu_scale2
endif

xlabel='(!3Frequency mod '+textoidl('\Delta') +') - ' +textoidl('\Delta/2') + ' ('+freq_u+')'
if show_ylabel eq 1 then begin
	ylabel='!3Frequency ('+freq_u+')'
endif else begin
	ylabel=''
endelse

if n_elements(bar) eq 0 then bar=1
	if print eq 0 then begin
		if bar eq 1 then tvframe,uu,xrange=[-window,window],yrange=[fbegin,fend],charsize=1.5,$
			ytitle=ylabel,xtitle=xlabel,/bar else $ ;,/center
		tvframe,uu,xrange=[-window,window],yrange=[fbegin,fend],charsize=1.5,$
			ytitle=ylabel,xtitle=xlabel;,/center
	endif
	if print ne 0 then begin
		if bar eq 1 then tvframe,-uu,xrange=[-window,window],yrange=[fbegin,fend],charsize=1.5,$
			ytitle=ylabel,xtitle=xlabel,/bar else $ ;,/center
		tvframe,-uu,xrange=[-window,window],yrange=[fbegin,fend],charsize=1.5,$
			ytitle=ylabel,xtitle=xlabel ;,/center
	endif
	if extralines[0] ne -1 then begin
		for i=0, n_elements(extralines)-1 do $
			if extralines[i]-window ge fbegin AND extralines[i]+window le fend then begin
				plots, [-window , window], [extralines[i], extralines[i]], color=fsc_color(pextralines_col[i]), linestyle=style_extralines[i]
			endif
	endif

;	if dots_extra[0] ne -1 then begin
;		for i=0, n_elements(dots_extra)-1 do plots, dots_extra[i] mod del, dots_extra[i], color=fsc_color(color[5]), psym=5, symsize=1.5
;	endif


; -------------- TABLE 1 ------------
	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
		if n_elements(nu_scale[*,0]) ge 6 then begin

			maxi=3
			test=where(nu_scale[5,*] ne 0) ; the last coloumn contains extra_dots
			if test[0] ne -1 then $
				oplot, nu_scale[5,where(nu_scale[5,*] ne 0)],freq_y[5,where(nu_scale[5,*] ne 0)],color=fsc_color(color[5]),psym=6,thick=3.,symsize=1.5

		endif else maxi=n_elements(nu_scale[*,0])-1


		if overplot eq 1 then for l=0,maxi do begin
				;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color)
				test=where(nu_scale[l,*] ne 0)
				if test[0] ne -1 then $
					oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
				endfor

		if overplot eq 2 then for l=0,maxi do begin
					test=where(nu_scale[l,*] ne 0)
					if test[0] ne -1 then begin
						tab=reform(nu_scale[l,where(nu_scale[l,*] ne 0)])
						freq_tab=reform(freq_y[l,where(nu_scale[l,*] ne 0)])
					endif
					if test[0] ne -1 then begin
						secu=0d & jj=0d & secu2=0d & ii=0d
						while (ii lt n_elements(tab)-1) AND secu lt 10 do begin
							d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
							d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
							secu2=0d
							while (d1_p le d2_p) AND (d1_m le d2_m) AND ii lt n_elements(tab)-2 AND secu2 lt 100 do begin
								;if ii eq 7 then stop
								plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(color[l])
								ii=ii+1
								d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
								d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
								secu2=secu2+1
							endwhile
							secu=secu+1
							ii=ii+1
						endwhile
						ii=ii-1
						if d1_p le d2_p then if n_elements(tab) ge 2 then plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(color[l])
						;oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
						oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l]),psym=4,thick=3.,symsize=1.5
					endif
				endfor
		if overplot eq 3 then for l=0,maxi do begin
					;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color),psym=4,symsize=1.5,thick=3
					test=where(nu_scale[l,*] ne 0)
					if test[0] ne -1 then $
						oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l]),psym=4,thick=3.,symsize=1.5
				endfor
		if overplot eq 4 then for l=0,n_elements(nu_scale[*,0])-1 do begin

					if l eq 0 then begin
						c='Orange' & t=0
					endif
					if l eq 1 then begin
						c='Yellow' & t=0
					endif
					if l eq 2 then begin
						c='red' & t=0
					endif
					if l eq 3 then begin
						c='brown' & t=0
					endif
					x_axis=findgen(nn)+nbeg+0.5
					;width=delta/4
					width=1d/4
					for i=0,n_elements(param_stat[*,0,0])-1 do $
						;Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,c], WIDTH=width, XLOC=stat_synthese_freq[3,l,i],1,t
						Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,'white'], WIDTH=width, XLOC=x_axis[i],1,t
					endfor
	endif
	if overplot eq 5 then stop

	if n_elements(doublepos) ne 0 then begin
		A = FINDGEN(17) * (!PI*2/16.)
		USERSYM, COS(A), SIN(A), /FILL
		sym_perso=dblarr(2, n_elements(A)) ; for the legends
		sym_perso[0,*]=cos(A) ; for the legends
		sym_perso[1,*]=sin(A) ; for the legends
		if posy_dble[0] ne -1 then $
			plots, nu_scale[posy_dble], freq_y[posy_dble], color=fsc_color(doublecolor), psym=8, symsize=2; else $
				;stop
	endif


; ------------ TABLE 2 ----------
	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, COS(A), SIN(A), thick=3 ;, /FILL
	if n_elements(freq_scd2) ne 0  then if  freq_t2[0] ne -1 then begin
		if n_elements(nu_scale2[*,0]) ge 6 then begin

			maxi=3
			test=where(nu_scale2[5,*] ne 0) ; the last coloumn contains extra_dots
			if test[0] ne -1 then $
				oplot, nu_scale2[5,where(nu_scale[5,*] ne 0)],freq_y2[5,where(nu_scale[5,*] ne 0)],color=fsc_color(scd_color[5]),psym=8,thick=8.,symsize=2

		endif else maxi=n_elements(nu_scale2[*,0])-1


		if overplot eq 1 then for l=0,maxi do begin
				;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color)
				test=where(nu_scale2[l,*] ne 0)
				if test[0] ne -1 then $
					oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale[l,*] ne 0)],color=fsc_color(scd_color[l])
				endfor

		if overplot eq 2 then for l=0,maxi do begin
					test=where(nu_scale2[l,*] ne 0)
					if test[0] ne -1 then begin
						tab=reform(nu_scale2[l,where(nu_scale2[l,*] ne 0)])
						freq_tab=reform(freq_y2[l,where(nu_scale2[l,*] ne 0)])
					endif
					if test[0] ne -1 then begin
						secu=0d & jj=0d & secu2=0d & ii=0d
						while (ii lt n_elements(tab)-1) AND secu lt 10 do begin
							d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
							d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
							secu2=0d
							while (d1_p le d2_p) AND (d1_m le d2_m) AND ii lt n_elements(tab)-2 AND secu2 lt 100 do begin
								;if ii eq 7 then stop
								plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(scd_color[l]), thick=4 ;, psym=8
								ii=ii+1
								d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
								d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
								secu2=secu2+1
							endwhile
							secu=secu+1
							ii=ii+1
						endwhile
						ii=ii-1
						if d1_p le d2_p then if n_elements(tab) ge 2 then plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(scd_color[l])
						;oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
						oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(scd_color[l]),psym=8,thick=10.,symsize=2
					endif
				endfor
		if overplot eq 3 then for l=0,maxi do begin
					;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color),psym=4,symsize=1.5,thick=3
					test=where(nu_scale2[l,*] ne 0)
					if test[0] ne -1 then $
						oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(scd_color[l]),psym=8,thick=8.,symsize=2
				endfor
		if overplot eq 4 then for l=0,n_elements(nu_scale2[*,0])-1 do begin

					if l eq 0 then begin
						c='Orange' & t=0
					endif
					if l eq 1 then begin
						c='Yellow' & t=0
					endif
					if l eq 2 then begin
						c='red' & t=0
					endif
					if l eq 3 then begin
						c='brown' & t=0
					endif
					x_axis=findgen(nn)+nbeg+0.5
					;width=delta/4
					width=1d/4
					for i=0,n_elements(param_stat[*,0,0])-1 do $
						;Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,c], WIDTH=width, XLOC=stat_synthese_freq[3,l,i],1,t
						Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,'white'], WIDTH=width, XLOC=x_axis[i],1,t
					endfor
	endif

if n_elements(tag_label) ne 0 then xyouts, tag_pos[0]*window/100 - window, tag_pos[1]*fend/100, tag_label, $
										charsize=tag_size, color=fsc_color(tag_color)

end



; This echelle diagram is for plotting 2 tables and showing the difference between these table on the left side of the echelle diagam
; This is very usefull if one want to show differences between model and measured frequencies, all together with the data
pro special_ED_v3, freq, spec_reg, delta, nbegin , n_number, freq_t, freq_m, err_freq_t, smooth_coef, ps, overplot, $
	color=color, scd_color=scd_color, extra_text=extra_text, x_text=x_text, y_text=y_text, size_text=size_text, clor_text=clor_text
!p.multi=[0,1,2]

size_char=0.75
size_char_axes=1.4
bar=0

if n_elements(extra_text) ne 0 then begin
	if n_elements(x_text) eq 0 then x_text=0.8
	if n_elements(y_text) eq 0 then y_text=0.85
endif
if n_elements(extra_text) eq 0 then begin
	extra_text=''
	x_text=0.
	y_text=0.
endif
if n_elements(size_text) eq 0 then size_text=size_char_axes
if n_elements(clor_text) eq 0 then clor_text=fsc_color('Black')


dif=freq_t-freq_m
err_dif=err_freq_t; - dif

nmax=floor(max(freq)/delta)-1
resol=1d*(freq[10]-freq[9])
nn=min([nmax,nbegin+n_number])
fmin=freq[round(1d*(nbegin)*delta/resol)] ;+delta/2
fmax=freq[round(1d*nn*delta/resol)] ;+ delta/2


xrange=[0.90*min(dif[where(freq_t ne 0)]-err_dif[where(freq_t ne 0)]), $
		1.05*max(dif[where(freq_t ne 0)]+ err_dif[where(freq_t ne 0)])]

col=color

!p.position=[0.35, 0.175, 0.9, 0.9]

!y.tickname=replicate(' ', 24)

recenter=0
echellecorot_v3,freq,spec_reg,delta,nbegin,uu,freq_table=freq_t, freq_scd2=freq_m,smooth_coef=smooth_coef,ps=ps,overplot=overplot,$
	color=color, scd_color=scd_color,show_ylabel=0, charsize=size_char_axes,$
	n_number=n_number-1,print=print, doublecolor=doublecolor, doublepos=doublepos, bar=bar, $
	;tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, $
	tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, tag_color=tag_color, $
	extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, recenter=recenter ;, rota=rota, congridon=congridon


!y.tickname=''
lmax=n_elements(freq_t[*,0])-1
for el=0, lmax do begin
	test=where(freq_t[el, *] ne 0)
	if el eq 0 then begin
		plot, dif[el,test], freq_t[el, test], color=fsc_color('Black'), background=fsc_color('white'),charsize=size_char_axes, $
			xtitle=textoidl('\nu_{obs} - \nu_{m} (\mu' + 'Hz)'), ytitle='Frequency (' + textoidl('\mu') + 'Hz)', $
			psym=4, yr=[fmin, fmax], xr=xrange, /xst, /yst, /nodata, $; ytickname=replicate(' ', 24), $
			pos=[0.125, 0.175, 0.33, 0.9]
	endif
	oplot, dif[el,test], freq_t[el, test], color=fsc_color(col[el]), psym=-4
	oploterr_astro, reform(dif[el,test]), reform(freq_t[el, test]), reform(err_dif[el, test]), reform(err_freq_t[el,test]), $
		color=fsc_color(col[el]), errcolor=fsc_color(col[el])
endfor
plots, [0, 0], [fmin, fmax], color=fsc_color('Dark Gray'), linestyle=2

for i=0, n_elements(extra_text)-1 do $
	xyouts, x_text[i], y_text[i], extra_text[i], /normal, charsize=size_text[i], color=clor_text[i]
end


; This echelle diagram is for plotting 2 tables and showing the difference between these table on the left side of the echelle diagam
; This is very usefull if one want to show differences between model and measured frequencies, all together with the data
; THIS IS THE V4 THAT ENSURE THAT FREQUENCIES ARE IN SAME LINES... SHOULD BE BETTER THAN V3 FOR PLOTING IN THE SIDE FREQUENCY DIFFERENCES
pro special_ED_v4, freq, spec_reg, delta, nbegin , n_number, freq_t, freq_m, err_freq_t, smooth_coef, vert_smooth, ps, overplot, $
	color=color, scd_color=scd_color, extra_text=extra_text, x_text=x_text, y_text=y_text, size_text=size_text, clor_text=clor_text,$
	correction_factor=correction_factor
!p.multi=[0,1,2]

size_char=0.75
size_char_axes=1.4
bar=0

if n_elements(extra_text) ne 0 then begin
	if n_elements(x_text) eq 0 then x_text=0.8
	if n_elements(y_text) eq 0 then y_text=0.85
endif
if n_elements(extra_text) eq 0 then begin
	extra_text=''
	x_text=0.
	y_text=0.
endif
if n_elements(size_text) eq 0 then size_text=size_char_axes
if n_elements(clor_text) eq 0 then clor_text=fsc_color('Black')


dif=freq_t-freq_m
err_dif=err_freq_t; - dif
;
nmax=floor(max(freq)/delta)-1
resol=1d*(freq[10]-freq[9])
nn=min([nmax,nbegin+n_number])
fmin=freq[round(1d*(nbegin)*delta/resol)] ;+delta/2
fmax=freq[round(1d*nn*delta/resol)] ;+ delta/2


xrange=[0.90*min(dif[where(freq_t ne 0)]-err_dif[where(freq_t ne 0)]), $
		1.05*max(dif[where(freq_t ne 0)]+ err_dif[where(freq_t ne 0)])]

col=color

!p.position=[0.35, 0.175, 0.9, 0.9]

!y.tickname=replicate(' ', 24)

spec_reg[where(freq eq fmin)+1]=max(spec_reg[where(freq eq fmin AND freq le fmax)])*correction_factor ; we put one point with high value to allow a good scaling of colors

recenter=0
echellecorot_v4,freq,spec_reg,delta,nbegin,uu,freq_table=freq_t, freq_scd2=freq_m,smooth_coef=smooth_coef,ps=ps,overplot=overplot,$
	color=color, scd_color=scd_color,show_ylabel=0, charsize=size_char_axes,$
	n_number=n_number,print=print, doublecolor=doublecolor, doublepos=doublepos, bar=bar, $
	vert_smooth=vert_smooth, $
	tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, tag_color=tag_color, $
	extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, recenter=recenter, $ ;, rota=rota, congridon=congridon
	fbegin=fbegin, fend=fend

fmin=fbegin
fmax=fend

!y.tickname=''
lmax=n_elements(freq_t[*,0])-1
for el=0, lmax do begin
	test=where(freq_t[el, *] ne 0)
	if el eq 0 then begin
		plot, dif[el,test], freq_t[el, test], color=fsc_color('Black'), background=fsc_color('white'),charsize=size_char_axes, $
			xtitle=textoidl('\nu_{obs} - \nu_{m} (\mu' + 'Hz)'), ytitle='Frequency (' + textoidl('\mu') + 'Hz)', $
			psym=4, yr=[fmin, fmax], xr=xrange, /xst, /yst, /nodata, $; ytickname=replicate(' ', 24), $
			pos=[0.125, 0.175, 0.33, 0.9]
	endif
	oplot, dif[el,test], freq_t[el, test], color=fsc_color(col[el]), psym=-4
	oploterr_astro, reform(dif[el,test]), reform(freq_t[el, test]), reform(err_dif[el, test]), reform(err_freq_t[el,test]), $
		color=fsc_color(col[el]), errcolor=fsc_color(col[el])
endfor
plots, [0, 0], [fmin, fmax], color=fsc_color('Dark Gray'), linestyle=2

for i=0, n_elements(extra_text)-1 do $
	xyouts, x_text[i], y_text[i], extra_text[i], /normal, charsize=size_text[i], color=clor_text[i]
end



; procedure that overplot frequencies on the echelle diagram of the data
; freq : the frequencies
; spec_reg : the power spectrum
; delta : the large separation
; nbegin : first n index
; uu : an output of the echelle diagram
; freq_table : input table of frequencies to overplot
; smooth_coef : coeficient of smoothing of the power spectrum (default = 1, ie no smoothing !)
; ps : if ps = 1 then we switch into device, decomposed=0 ... resolve problems when ploting on a ps file
; overplot : nature of the overplot (default overplot=2),
	; - overplot=1 then we show only lines for the input frequency table. Table syntax: [0:lmax, 0:Nmax]
	; - overplot=2 then show lines + diamonds for the input frequency table. Table syntax: [0:lmax, 0:Nmax]
	; - overplot=3 then show only diamonds for the input frequency table. Table syntax: [0:lmax, 0:Nmax]
	; - overpolot=4 then use whisker boxes ... frequency_table MUST BE a table COMPATIBLE WITH THE FUNCTION : Draw_BoxAndWiskers_horizontal
	; - overplot=5 then we overplot 2 echelle diagram. For example, one containing true data and one containing a model of the data (variable spec2)
; color : color code to be used (default color='red')... used by the function fsc_color.

;ADDED 01/04/2014 : freq_scd2... usefull to overplot two tables : eg. model and obs
pro echellecorot_v3,freq,spec_reg,delta,nbegin,uu,freq_table=freq_t,freq_scd2=freq_t2, smooth_coef=smooth_coef,ps=ps,overplot=overplot,$
	color=color, scd_color=scd_color, show_ylabel=show_ylabel, charsize=charsize,$
	n_number=n_number,print=print, doublecolor=doublecolor, doublepos=doublepos, bar=bar, $
	;tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, $
	tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, tag_color=tag_color, set_plot=set_plot,$
	extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, recenter=recenter ;, rota=rota, congridon=congridon


A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A);, /FILL
sym_perso=dblarr(2, n_elements(A)) ; for the legends
sym_perso[0,*]=cos(A) ; for the legends
sym_perso[1,*]=sin(A) ; for the legends

printer='' ; relicat of the appourchaux function
index=''
;bar=1

if n_elements(set_plot) eq 0 then set_plot='X' ; Default, we plot on a linux system
if n_elements(recenter) eq 0 then recenter=0
if n_elements(n_number) eq 0 then n_number= 21
if n_elements(smooth_coef) eq 0 then smooth_coef=1
;if n_elements(tag_label) eq 0 then tag_label=''
if n_elements(tag_pos) eq 0 then tag_pos=[5.1, 90]
if n_elements(tag_size) eq 0 then tag_size=1.
;if n_elements(rota) eq 0 then rota=0
;if n_elements(congridon) eq 0 then congridon=0

if n_elements(show_ylabel) eq 0 then show_ylabel=1
if n_elements(charsize) eq 0 then charsize=2.
a=smooth(spec_reg,smooth_coef)


if n_elements(extralines) eq 0 then begin
	extralines=-1
	extralines_color=-1
endif
if n_elements(extralines) ne n_elements(pextralines_col) then pextralines_col=strarr(n_elements(extralines)) + 'Black'


n_numb=n_number
if nbegin lt 0 then nbegin=0
nbeg=nbegin

if n_elements(freq_t) ne 0  then if  freq_t[0] ne -1 then begin
if overplot ne 4 AND overplot ne 5 then freq_table=freq_t
if overplot eq 4 then begin ; if overplot = 4 we assume that the table is COMPATIBLE with Draw_BoxAndWiskers_horizontal
	stat_synthese_freq=freq_t
	freq_table=dblarr(n_elements(stat_synthese_freq[0,*,0]),n_elements(stat_synthese_freq[0,0,*]))
	freq_table[*,*]=stat_synthese_freq[3,*,*] ; freq_table contain ONLY THE MEDIAN POSITION
endif
endif

if n_elements(freq_t2) ne 0  then if  freq_t2[0] ne -1 then begin
if overplot ne 4 AND overplot ne 5 then freq_table2=freq_t2
if overplot eq 4 then begin ; if overplot = 4 we assume that the table is COMPATIBLE with Draw_BoxAndWiskers_horizontal
	stat_synthese_freq2=freq_t2
	freq_table2=dblarr(n_elements(stat_synthese_freq2[0,*,0]),n_elements(stat_synthese_freq2[0,0,*]))
	freq_table2[*,*]=stat_synthese_freq2[3,*,*] ; freq_table contain ONLY THE MEDIAN POSITION
endif
endif

if n_elements(overplot) eq 0 then overplot=-1

	nn=N_params()
	s={echelleplot,windowi:0.,delta:0.,echelletitle:''}
	if (nn LT N_tags(s)) then begin
		nn=nn-1
		read_structure,s,nn
		names=tag_names(s)
		for i=nn,N_tags(s)-1 do begin
			zzz=names(i)+'=s.'+names(i)
			;print,zzz
			r=execute(names(i)+'=s.'+names(i))
		endfor
	endif

	if (printer ne '') then begin
		set_plot,'ps',/interpolate
		file='corotvg'+index+'.ps'
		device,/landscape,/color,bits=8,filename=file
	endif
	if (n_elements(ps) ne 0) then begin
		set_plot,'ps',/interpolate
		device,/times
	endif

	del=1d*delta

	windowi=del

	window=windowi/2.

	Na=N_elements(a)
	resol=1d*(freq[10]-freq[9])

	;print,'Na=',Na
	;print,'resolution',resol
	Nz=1d*del/resol+2
	b=fltarr(Nz,n_numb)
	if overplot eq 5 then b2=fltarr(Nz,n_numb)
	middlefreq=1d*del/2d/resol

	nmax=floor(max(freq)/del)-1
	;print,"nmax=",nmax
;	nn=min([nmax,nbeg+20])
	nn=min([nmax,nbeg+n_numb-1])

	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
		nu=dblarr(n_elements(freq_table[*,0]),n_elements(freq_table[0,*])*3)
		freq_y=dblarr(n_elements(freq_table[*,0]),n_elements(freq_table[0,*])*3)
	endif
	if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
		nu2=dblarr(n_elements(freq_table2[*,0]),n_elements(freq_table2[0,*])*3)
		freq_y2=dblarr(n_elements(freq_table2[*,0]),n_elements(freq_table2[0,*])*3)
	endif

	if nbeg lt 0 then begin
		nbeg=0
		nn=10
	endif
	iend_max=round(1d*(nn+1)*del/resol)
	if iend_max gt n_elements(a) then nn=fix(n_elements(a)*resol/del)-1 ; if we are at upper edge of the spectrum!
	corr_t=0d
	corr_t2=0d
	for i=nbeg,nn do begin
		ibegin=round(1d*i*del/resol)
		iend=round(1d*(i+1)*del/resol)
		b(0:(iend-ibegin),i-nbeg)=a(ibegin:iend)
		if overplot eq 5 then b2(0:(iend-ibegin),i-nbeg)=a2(ibegin:iend)

		;print,i,min(b(*,i-nbeg)),max(b(*,i-nbeg))
		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
			f1=freq[floor(ibegin)] & f2=freq[ceil(ibegin)]
			f= (f2 - f1) * (ibegin - floor(ibegin)) + f1
			corr_t=corr_t+ f/(resol*ibegin)
		endif
		if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
			f1=freq[floor(ibegin)] & f2=freq[ceil(ibegin)]
			f= (f2 - f1) * (ibegin - floor(ibegin)) + f1
			corr_t2=corr_t2+ f/(resol*ibegin)
		endif
	endfor

	if n_elements(freq_table) ne 0 then begin

		lmax=n_elements(freq_table[*,0])-1

	corr_t=corr_t/(nn-nbeg +1)
	kk=0d & phase=0d
	for i=nbeg,nn do begin
		loop=0d
		ibegin=1d*i*del/resol
		iend=1d*(i+1)*del/resol ;-1

		l_tmpQ=dblarr(lmax+1)
		tmpQ=dblarr(lmax+1, 20)-1
		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
			for l=0, lmax do begin
				ta=where(freq_table[l,*] ge 1d*freq[ibegin] AND freq_table[l,*] le 1d*freq[iend])
				if ta[0] ne -1 then begin
					tmpQ[l, 0:n_elements(ta)-1]=ta
					l_tmpQ[l]=n_elements(ta)
				endif
			endfor
			dup_max=max(l_tmpQ)
			for j=0, dup_max-1 do begin
				for l=0, lmax do begin
					if tmpQ[l,j] ne -1 then begin
						nu[l, kk]=freq_table[l, tmpQ[l,j]]- corr_t*resol*ibegin
						freq_y[l, kk]=freq_table[l, tmpQ[l,j]]
					endif
				endfor
				kk=kk+1d
			endfor
		endif
	endfor

	endif
; --------- table 2 ---------
if n_elements(freq_table2) ne 0 then begin

	lmax=n_elements(freq_table2[*,0])-1

	corr_t2=corr_t2/(nn-nbeg +1)
	kk=0d & phase=0d
	for i=nbeg,nn do begin
		loop=0d
		ibegin=1d*i*del/resol
		iend=1d*(i+1)*del/resol ;-1

		l_tmpQ2=dblarr(lmax+1)
		tmpQ2=dblarr(lmax+1, 20)-1
		if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
			for l=0, lmax do begin
				ta=where(freq_table2[l,*] ge 1d*freq[ibegin] AND freq_table2[l,*] le 1d*freq[iend])
				if ta[0] ne -1 then begin
					tmpQ2[l, 0:n_elements(ta)-1]=ta
					l_tmpQ2[l]=n_elements(ta)
				endif
			endfor
			dup_max=max(l_tmpQ2)
			for j=0, dup_max-1 do begin
				for l=0, lmax do begin
					if tmpQ2[l,j] ne -1 then begin
						nu2[l, kk]=freq_table2[l, tmpQ2[l,j]]- corr_t2*resol*ibegin
						freq_y2[l, kk]=freq_table2[l, tmpQ2[l,j]]
					endif
				endfor
				kk=kk+1d
			endfor
		endif
	endfor

	endif

;	help,b
	;if congridon eq 0 then zz=b else zz=congrid(b, Nz, Nz)
	zz=b
	if recenter eq 1 then uu=zz(middlefreq-window/resol:middlefreq+window/resol,*)
	if recenter eq 0 then begin
		uu=zz ;(middlefreq-window/resol:middlefreq+window/resol,*)
		;middlefreq=0
	endif

	if overplot eq 5 AND recenter eq 1 then begin
		zz2=congrid(b2,Nz,100)
		uu2=zz2(middlefreq-window/resol:middlefreq+window/resol,*)
	endif
	if overplot eq 5 AND recenter eq 0 then begin
		zz2=congrid(b2,Nz,100)
		uu2=zz2 ;(middlefreq-window/resol:middlefreq+window/resol,*)
	endif

	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin

		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then nu_scale=nu-middlefreq*resol*(1d +(corr_t-1)*2)

		pos_dummy=where(nu_scale eq -middlefreq*resol*(1d +(corr_t-1)*2)) ; these positions had 0... we need to let the 0
		nu_scale[pos_dummy]=0d

		if overplot eq 4 then begin
			param_stat=dblarr(n_elements(nu_scale[0,*]),n_elements(nu_scale[*,0]),n_elements(stat_synthese_freq[*,0,0]))

			for i=0,n_elements(nu_scale[0,*])-1 do begin
				for j=0,n_elements(stat_synthese_freq[*,0,0])-1 do begin
					for l=0,n_elements(nu_scale[*,0])-1 do $
						param_stat[i,l,j]=nu_scale[l,i]+(stat_synthese_freq[j,l,i]-stat_synthese_freq[3,l,i]) ; using nu_scale, we compute the error boxes
				endfor
			endfor
		endif

		if n_elements(overplot) eq 0 then overplot=2 ; definee the defaut value of overplot : we plot lines + diamonds !
		if n_elements(color) eq 0 then color=['Red','Red','Red','Red','Red']
		if n_elements(color) eq 1 then color=[color,color,color] ; 1 color code by degree
		if n_elements(scd_color) eq 0 then color_scd=color
		if n_elements(scd_color) eq 1 then color_scd=[scd_color, scd_color,scd_color]
		if n_elements(scd_color) gt 1 then color_scd=scd_color

		x0=-delta/2
		no_zero=where(freq_y ne 0)

		if n_elements(doublepos) ne 0 then $
			posy_dble=where(freq_y eq doublepos) ; find the position of the specially tagged color !!

		freq_y[no_zero]=freq_y[no_zero] - (nu_scale[no_zero]- x0) ;+ window ; the last term is to center on the y axis the frequency plot

	endif

; ----- Table 2 ------
	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, COS(A), SIN(A) ;, /FILL
	if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin

		if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then nu_scale2=nu2-middlefreq*resol*(1d +(corr_t2-1)*2)

		pos_dummy2=where(nu_scale2 eq -middlefreq*resol*(1d +(corr_t2-1)*2)) ; these positions had 0... we need to let the 0
		nu_scale2[pos_dummy2]=0d

		if overplot eq 4 then begin
			param_stat2=dblarr(n_elements(nu_scale2[0,*]),n_elements(nu_scale2[*,0]),n_elements(stat_synthese_freq2[*,0,0]))

			for i=0,n_elements(nu_scale2[0,*])-1 do begin
				for j=0,n_elements(stat_synthese_freq2[*,0,0])-1 do begin
					for l=0,n_elements(nu_scale2[*,0])-1 do $
						param_stat2[i,l,j]=nu_scale2[l,i]+(stat_synthese_freq2[j,l,i]-stat_synthese_freq2[3,l,i]) ; using nu_scale, we compute the error boxes
				endfor
			endfor
		endif

		if n_elements(overplot) eq 0 then overplot=2 ; definee the defaut value of overplot : we plot lines + diamonds !
		if n_elements(color) eq 0 then color=['Red','Red','Red','Red','Red']
		if n_elements(color) eq 1 then color=[color,color,color] ; 1 color code by degree

		x0=-delta/2
		no_zero=where(freq_y2 ne 0)

		if n_elements(doublepos2) ne 0 then $
			posy_dble=where(freq_y2 eq doublepos2) ; find the position of the specially tagged color !!

		freq_y2[no_zero]=freq_y2[no_zero] - (nu_scale2[no_zero]- x0) ;+ window ; the last term is to center on the y axis the frequency plot

	endif

	help,uu
	freq_u=TeXtoIDL('\mu')+'Hz'
	print,max(uu)


	fbegin=freq[round(1d*(nbeg)*del/resol)] ;+delta/2
	fend=freq[round(1d*nn*del/resol)] ;+ delta/2
	print, 'Frequency range:', fbegin, fend

if n_elements(ps) eq 0 then begin
	set_plot,set_plot
	device, decomposed=0
endif
if n_elements(ps) eq 1 then if ps eq 0 then begin
	set_plot,set_plot
	device, decomposed=0
endif
!p.font=0 ;& device,set_font='times' ; police de caractere times
;loadct,39  ; met une table de couleur avec 255 = Noir et 1 = Blanc
loadct,39
if n_elements(print) eq 0 then print=0
if print eq 1 then loadct, 0 ; table with a white background

!EDIT_INPUT=50

uu0=uu

xlabel='!3Frequency ('+freq_u+')'
if show_ylabel eq 1 then begin
	ylabel='!3Frequency ('+freq_u+')'
endif else begin
	ylabel=''
endelse

;if recenter eq 0 then begin
;	bmin=0
;	bmax=delta
;endif
if recenter eq 1 OR recenter eq 0 then begin
	bmin=-window
	bmax=window
endif

if n_elements(bar) eq 0 then bar=1
	if print eq 0 then begin
		if bar eq 1 then tvframe,uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel,/bar else $ ;,/center
		tvframe,uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel;,/center
	endif
	if print ne 0 then begin
		if bar eq 1 then tvframe,-uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel,/bar else $ ;,/center
		tvframe,-uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel ;,/center
	endif
	if extralines[0] ne -1 then begin
		for i=0, n_elements(extralines)-1 do $
			if extralines[i]+ bmin ge fbegin AND extralines[i]+bmax le fend then begin
;				norm0=sqrt( ( -window + window)^2 + (extralines[i] + window - nu_center)^2 )
;				norm1=sqrt( ( window + window)^2 + (extralines[i] - window - nu_center)^2 )
;				x0= -window + norm0 * cos( rota*!pi/180d + atan( (extralines[i] + window - nu_center)/(-window + window)) )
;				x1= -window + norm1 * cos( rota*!pi/180d + atan( (extralines[i] - window - nu_center)/( window + window)) )
;				y0=nu_center + norm0 * sin( rota*!pi/180d + atan( (extralines[i] + window - nu_center)/(-window + window)) )
;				y1=nu_center + norm1 * sin( rota*!pi/180d + atan( (extralines[i] - window - nu_center)/( window + window)) )
;				plots, [x0 , x1], [y0, y1], color=fsc_color(pextralines_col[i]), linestyle=style_extralines[i]
				plots, [bmin , bmax], [extralines[i]+ bmax, extralines[i]+bmin], color=fsc_color(pextralines_col[i]), linestyle=style_extralines[i]
				;plots, [-window , window], [extralines[i], extralines[i]], color=fsc_color(pextralines_col[i]), linestyle=style_extralines[i]
			endif
	endif


	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
		if n_elements(nu_scale[*,0]) ge 6 then begin

			maxi=3
			test=where(nu_scale[5,*] ne 0) ; the last coloumn contains extra_dots
			if test[0] ne -1 then $
				oplot, nu_scale[5,where(nu_scale[5,*] ne 0)],freq_y[5,where(nu_scale[5,*] ne 0)],color=fsc_color(color[5]),psym=4,thick=3.,symsize=1.5

		endif else maxi=n_elements(nu_scale[*,0])-1

		if overplot eq 1 then for l=0,maxi do begin
				;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color)
				test=where(nu_scale[l,*] ne 0)
				if test[0] ne -1 then $
					oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
				endfor

		if overplot eq 2 then for l=0,maxi do begin
					test=where(nu_scale[l,*] ne 0)
					if test[0] ne -1 then begin
						tab=reform(nu_scale[l,where(nu_scale[l,*] ne 0)])
						freq_tab=reform(freq_y[l,where(nu_scale[l,*] ne 0)])
					endif
					if test[0] ne -1 then begin
						secu=0d & jj=0d & secu2=0d & ii=0d
						while (ii lt n_elements(tab)-1) AND secu lt 10 do begin
							d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
							d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
							secu2=0d
							while (d1_p le d2_p) AND (d1_m le d2_m) AND ii lt n_elements(tab)-2 AND secu2 lt 100 do begin
								;if ii eq 7 then stop
								plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(color[l])
								ii=ii+1
								d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
								d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
								secu2=secu2+1
							endwhile
							secu=secu+1
							ii=ii+1
						endwhile
						ii=ii-1
						if d1_p le d2_p then if n_elements(tab) ge 2 then plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(color[l])
						;oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
						oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l]),psym=4,thick=3.,symsize=1.5
					endif
				endfor
		if overplot eq 3 then for l=0,maxi do begin
					;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color),psym=4,symsize=1.5,thick=3
					test=where(nu_scale[l,*] ne 0)
					if test[0] ne -1 then $
						oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l]),psym=4,thick=3.,symsize=1.5
				endfor
		if overplot eq 4 then for l=0,n_elements(nu_scale[*,0])-1 do begin

					if l eq 0 then begin
						c='Orange' & t=0
					endif
					if l eq 1 then begin
						c='Yellow' & t=0
					endif
					if l eq 2 then begin
						c='red' & t=0
					endif
					if l eq 3 then begin
						c='brown' & t=0
					endif
					x_axis=findgen(nn)+nbeg+0.5
					;width=delta/4
					width=1d/4
					for i=0,n_elements(param_stat[*,0,0])-1 do $
						;Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,c], WIDTH=width, XLOC=stat_synthese_freq[3,l,i],1,t
						Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,'white'], WIDTH=width, XLOC=x_axis[i],1,t
					endfor
	endif

; ------- Table 2 -------
	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, COS(A), SIN(A) ;, /FILL
	if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
		if n_elements(nu_scale2[*,0]) ge 6 then begin

			maxi=3
			test=where(nu_scale2[5,*] ne 0) ; the last coloumn contains extra_dots
			if test[0] ne -1 then $
				oplot, nu_scale2[5,where(nu_scale2[5,*] ne 0)],freq_y2[5,where(nu_scale2[5,*] ne 0)],color=fsc_color(color_scd[5]),psym=8,thick=5.,symsize=1.5

		endif else maxi=n_elements(nu_scale2[*,0])-1

		if overplot eq 1 then for l=0,maxi do begin
				;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color)
				test=where(nu_scale2[l,*] ne 0)
				if test[0] ne -1 then $
					oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(color_scd[l]),psym=8
				endfor

		if overplot eq 2 then for l=0,maxi do begin
					test=where(nu_scale2[l,*] ne 0)
					if test[0] ne -1 then begin
						tab2=reform(nu_scale2[l,where(nu_scale2[l,*] ne 0)])
						freq_tab2=reform(freq_y2[l,where(nu_scale2[l,*] ne 0)])
					endif
					if test[0] ne -1 then begin
						secu=0d & jj=0d & secu2=0d & ii=0d
						while (ii lt n_elements(tab)-1) AND secu lt 10 do begin
							d1_p=abs(tab2[ii+1] -tab2[ii]) & d2_p=abs(tab2[ii+1]+2*window -tab2[ii])
							d1_m=abs(tab2[ii+1] -tab2[ii]) & d2_m=abs(tab2[ii+1]-2*window -tab2[ii])
							secu2=0d
							while (d1_p le d2_p) AND (d1_m le d2_m) AND ii lt n_elements(tab2)-2 AND secu2 lt 100 do begin
								;if ii eq 7 then stop
								plots, [tab2[ii], tab2[ii+1]], [freq_tab2[ii], freq_tab2[ii+1]], color=fsc_color(color_scd[l])
								ii=ii+1
								d1_p=abs(tab2[ii+1] -tab2[ii]) & d2_p=abs(tab2[ii+1]+2*window -tab2[ii])
								d1_m=abs(tab2[ii+1] -tab2[ii]) & d2_m=abs(tab2[ii+1]-2*window -tab2[ii])
								secu2=secu2+1
							endwhile
							secu=secu+1
							ii=ii+1
						endwhile
						ii=ii-1
						if d1_p le d2_p then if n_elements(tab2) ge 2 then plots, [tab2[ii], tab2[ii+1]], [freq_tab2[ii], freq_tab2[ii+1]], color=fsc_color(color_scd[l])
						;oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
						oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(color_scd[l]),psym=8,thick=3.,symsize=1.5
					endif
				endfor
		if overplot eq 3 then for l=0,maxi do begin
					;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color),psym=8,symsize=1.5,thick=3
					test=where(nu_scale2[l,*] ne 0)
					if test[0] ne -1 then $
						oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(color_scd[l]),psym=8,thick=3.,symsize=1.5
				endfor
		if overplot eq 4 then for l=0,n_elements(nu_scale2[*,0])-1 do begin

					if l eq 0 then begin
						c='Orange' & t=0
					endif
					if l eq 1 then begin
						c='Yellow' & t=0
					endif
					if l eq 2 then begin
						c='red' & t=0
					endif
					if l eq 3 then begin
						c='brown' & t=0
					endif
					x_axis=findgen(nn)+nbeg+0.5
					;width=delta/4
					width=1d/4
					for i=0,n_elements(param_stat[*,0,0])-1 do $
						;Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,c], WIDTH=width, XLOC=stat_synthese_freq[3,l,i],1,t
						Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,'white'], WIDTH=width, XLOC=x_axis[i],1,t
					endfor
	endif

	if overplot eq 5 then stop

	if n_elements(doublepos) ne 0 then begin
		A = FINDGEN(17) * (!PI*2/16.)
		USERSYM, COS(A), SIN(A), /FILL
		sym_perso=dblarr(2, n_elements(A)) ; for the legends
		sym_perso[0,*]=cos(A) ; for the legends
		sym_perso[1,*]=sin(A) ; for the legends
		if posy_dble[0] ne -1 then $
			plots, nu_scale[posy_dble], freq_y[posy_dble], color=fsc_color(doublecolor), psym=8, symsize=2; else $
				;stop
	endif

if n_elements(tag_label) ne 0 then xyouts, tag_pos[0]*window/100 - window, tag_pos[1]*fend/100, tag_label, $
										charsize=tag_size, color=fsc_color(tag_color)


;	if n_elements(ps) eq 0 then device,/close
;	if n_elements(ps) eq 1 then if ps eq 0 then device,/close
end


pro show_ech_diag_v3, dir_fit, file_spec, ps=ps, extended=extended,smooth_coef=smooth_coef, $
	bckgrd_cor=bckgrd_cor,shifts=shifts, file_out=file_out, print=print, bar=bar, noise_file=noise_file, $
	tag_m=tag_m, tag_pos=tag_pos, tag_color=tag_color, tag_size=tag_size, trunc_spec=trunc_spec, a_law=a_law, show_spec=show_spec

if n_elements(ps) eq 0 then ps=0
if n_elements(extended) eq 0 then extended=0 else cmax=extended
if n_elements(smooth_coef) eq 0 then smooth_coef=12
if n_elements(bckgrd_cor) eq 0 then cor=1 else cor=bckgrd_cor
if n_elements(shifts) eq 0 then shifts=0 ; no shift by defaut
if n_elements(file_out) eq 0 then file_out2=dir_fit else file_out2=file_out
if n_elements(print) eq 0 then print=1 ; To have Blackk/white data
if n_elements(bar) eq 0 then bar=0 ; by defaut we put the graduation bar
if n_elements(tag_color) eq 0 then tag_color='Black'
if n_elements(show_spec) eq 0 then show_spec=0

;restore, 'C:\Work_dir\Analysis_Results\TimWTwin-B\2.1M\Syntheses\10124866-StarB_Mtmp_data_enhanced.sav'
;tab_B=reform(stat_synthese_freq[3,*,*])
;restore, 'C:\Work_dir\Analysis_Results\TimWTwin-B\2.1M\Syntheses\post_spect_param.sav'
;model_spectrum2=model_spectrum

;tab_B[0,0]=0 & tab_B[0,2:3]=0
;tab_B[0,11:*]=0
;
;tab_B[1,0]=0  & tab_B[1,2]=0
;
;tab_B[2,0:1]=0 & tab_B[2,3:5]=0
;tab_B[2,11:*]=0
;

file_fit=file_search(dir_fit+'*_data_enhanced.sav')
file_postfit=file_search(dir_fit+'post_spect_param.sav')
if file_postfit eq '' then file_postfit=file_search(dir_fit+'*postspec.sav') ; try the second syntax
if file_postfit eq '' then stop
restore, file_fit
restore, file_postfit
restore, file_spec

freq_tab_t0=reform(stat_synthese_freq[3,*,*])
;Nf_max=n_elements(tab_B[0,*])

l=1
freq_tab_t0[l,*]=freq_tab_t0[l,sort(freq_tab_t0[l,*])]

;freq_tab_t0[0,0:4]=0
;freq_tab_t0[0,15]=0
;
;freq_tab_t0[1,0:2]=0
;freq_tab_t0[1,15]=0
;
;freq_tab_t0[2,0:4]=0
;freq_tab_t0[2,14:15]=0

freq_tab_t00=freq_tab_t0

;freq_tab=dblarr(6, 16)
;freq_tab[0:2,*]=freq_tab_t0
;freq_tab[3:*,0:Nf_max-1]=tab_B

;freq_tab_t0=freq_tab

Nf_l0=n_elements(freq_tab_t0[0,where(freq_tab_t0[0,*] ne 0)])
Nf_l1=n_elements(freq_tab_t0[1,where(freq_tab_t0[1,*] ne 0)])

test=where(freq_tab_t0[0,*] ne 0)
;residu=continuity_condition(reform(freq_tab_t0[0,0:Nf_l0-1]),  n_elements(freq_tab_t0[0,0:Nf_l0-1]), n_elements(freq_tab_t0[0,0:Nf_l0-1]))
residu=continuity_condition(reform(freq_tab_t0[0,test]),  n_elements(freq_tab_t0[0,test]), n_elements(freq_tab_t0[0,test]))

residu_l1=continuity_condition(reform(freq_tab_t0[1,0:Nf_l1-1]),  n_elements(freq_tab_t0[1,0:Nf_l1-1]), n_elements(freq_tab_t0[1,0:Nf_l1-1]))

Dnu=mean(residu.residual_st)
;Dnu=155.963384

r=min(freq_tab_t0[1,where(freq_tab_t0[1,*] ne 0)])/Dnu
n0=floor(r)
epsilon=r-n0
n0=n0-1
;n0=11

overplot=3 ; show diamonds only
color=['Orange', 'Green','Red', 'Brown', 'Yellow','Blue']
;color=['Orange', 'Red','Dark Gray', 'Orange', 'Red','Dark Gray']

n_number=n_elements(freq_tab_t0[0,*]) +5
n0=n0-4

if ps eq 0 then window, 1, xsize=1400, ysize=800 else $
	nimp,name=file_out2+'Ech_diag.eps',/paper,/eps ; Savita's code
	;e=write_on_ps_on(file_fit+'_EchDiag')

if n_elements(freq) eq 0 then begin
	print, 'BEWARE: We didn t find the original power spectrum (files freq and spec_reg)'
	;print, 'Trying to use the file_spec option...'
	;if n_elements(file_spec) ne 0 then restore, file_spec else begin
	print, 'file_spec has not been specified properly ! Emmergency stop!'
	stop
	;endelse
endif
if max(freq) lt 10 then x0=freq*1d6 else x0=freq
if n_elements(noise_file) ne 0 then begin
	restore, noise_file
	backgrd=N
	N_local=interpol(N, freq, x)
	model_spectrum=model_spectrum+N_local
endif
if cor eq 1 then begin
	val_des_median=reform(stat_synthese_unsorted[3,*])
	 noise_param=val_des_median[total(parameters_length[0:4]):total(parameters_length[0:5])-1]
	 N=0d & cpt=0d
	 for i=0,1 do begin
		N=N + noise_param[cpt]/(1d  + (noise_param[cpt +1]*x0*1d-3)^noise_param[cpt+2])
		cpt=cpt+3
	 endfor
	 N=N+noise_param[6] ; we add the white noise
	 N=N+dblarr(n_elements(freq)) ;+ 0.22

	 if n_elements(noise_file) ne 0 then N=backgrd ; use the restored file

	 s0=spec_reg/N
endif else s0=spec_reg

if n_elements(trunc_spec) eq 0 then limit=max(s0) else limit=trunc_spec
test=where(s0 ge limit)
if test[0] ne -1 then s0[test]=limit
if extended eq 0 OR extended eq 1 then $
	echellecorot_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,overplot=overplot,color=color,n_number=n_number, print=print, bar=bar
if extended ge 2 then $
	echellecorot_UNFOLD_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,$
		overplot=overplot,color=color, n_number=n_number,cmax=cmax, print=print, bar=bar, $
		tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color, a_law=a_law
if ps ne 0 then $
	fimp

if show_spec eq 1 then begin
if ps eq 0 then window, 0, xsize=1400, ysize=800

x_min=3500
if ps ne 0 then $
	nimp,name=dir_fit+'SpecFit.eps',/paper,/eps ; Savita's code
	;e=write_on_ps_on(dir_fit + 'SpecFit')
	if n0*Dnu lt min(x) then x_min=min(x) else x_min=n0*Dnu
	if (n0+n_number)*Dnu gt max(x) then x_max=max(x) else x_max=(n0+n_number)*Dnu
plot, x,smooth(s, smooth_coef),xr=[x_min, x_max],/xst, thick=1
oplot, x, model_spectrum, color=fsc_color('blue')
oplot, freq, N , color=fsc_color('red')
;Nloc=N[118815:271579]
;oplot, x, model_spectrum2+Nloc, color=fsc_color('Red'), thick=5
; + N[where(freq eq x)]

level=0.999
if smooth_coef ge 30 then smooth_coef2=30 else smooth_coef2=smooth_coef
threshold=proba_calc_smooth_signal(1d + 0.05, smooth_coef, level)

plot, x,smooth(s/model_spectrum, smooth_coef2),xr=[x_min, x_max],/xst, thick=1, title='Residual spectrum'
oplot, x,  dblarr(n_elements(x)) + threshold, color=fsc_color('Red'), linestyle=2, thick=2
oplot, x,  dblarr(n_elements(x)) + 1, color=fsc_color('Red'), linestyle=2, thick=2

if ps ne 0 then fimp

if ps ne 0 then begin
	x_max=0
	dn=3 & i=0
	while x_max lt max(x) do begin
		nimp,name=dir_fit+'SpecFit-'+strtrim(i,1)+'.eps',/paper,/eps ; Savita's code
		x_min=min(x) + i*dn*Dnu & x_max=min(x) + (i+1)*dn*Dnu
		if x_max gt max(x) then x_max=max(x)
		plot, x,smooth(s, smooth_coef2),xr=[x_min, x_max],/xst, thick=1
		oplot, x, model_spectrum, color=fsc_color('blue'), thick=5
		;oplot, x, model_spectrum2+Nloc, color=fsc_color('Red'), thick=5
		oplot, freq, N , color=fsc_color('red')
		i=i+1
		wait, 0.1
		fimp
	endwhile
;	e=write_on_ps_off('')
endif
endif
print, 'yo'
end


pro show_ech_diag_v3_CPP, synthese_file, data, ps=ps, shifts=shifts, trunc_spec=trunc_spec

restore, synthese_file

freq=reform(data[0,*])
spec_reg=reform(data[1,*])
noise=reform(data[2,*])

fit_freq=linfit(findgen(n_elements(freq)), freq) ; This fit is way much more precise than the difference of the first points (subject to round-off errors)
resol=fit_freq[1] ; fit_freq[0] is by definition min(freq0) and fit_freq[1] is the exact resolution

smooth_coef=median(stat_synthese_width[3,0,*])/resol/2.

s0=spec_reg/noise

if n_elements(ps) eq 0 then ps=0
if n_elements(shifts) eq 0 then shifts=0 ; no shift by defaut
if n_elements(trunc_spec) eq 0 then limit=max(s0) else limit=trunc_spec
cor=1
print=0 ; To have Black/white data or Color
bar=0; by defaut we put the graduation bar
extended=0

freq_tab_t0=reform(stat_synthese_freq[3,*,*])

l=1
freq_tab_t0[l,*]=freq_tab_t0[l,sort(freq_tab_t0[l,*])]

Nf_l0=n_elements(freq_tab_t0[0,where(freq_tab_t0[0,*] gt 0)])
Nf_l1=n_elements(freq_tab_t0[1,where(freq_tab_t0[1,*] gt 0)])

residu=continuity_condition(reform(freq_tab_t0[0,0:Nf_l0-1]),  n_elements(freq_tab_t0[0,0:Nf_l0-1]), n_elements(freq_tab_t0[0,0:Nf_l0-1]))

residu_l1=continuity_condition(reform(freq_tab_t0[1,0:Nf_l1-1]),  n_elements(freq_tab_t0[1,0:Nf_l1-1]), n_elements(freq_tab_t0[1,0:Nf_l1-1]))

c=linfit(findgen(Nf_l0),reform(freq_tab_t0[0,0:Nf_l0-1]) ) 
Dnu=c[1]
;if n_elements(Dnu) eq 0 then Dnu=mean(residu.residual_st)

r=min(freq_tab_t0[1,where(freq_tab_t0[1,*] ne 0)])/Dnu
n0=floor(r)
epsilon=r-n0
n0=n0-1

overplot=2 ; show diamonds only
color=['Orange', 'Red','Dark Gray', 'Yellow']

n_number=n_elements(freq_tab_t0[0,*])+2
n0=n0-2

if n_elements(freq) eq 0 then begin
	print, 'BEWARE: We didn t find the original power spectrum (files freq and spec_reg)'
	print, 'Trying to use the file_spec option...'
	if n_elements(file_spec) ne 0 then restore, file_spec else begin
		print, 'file_spec has not been specified ! Emmergency stop!'
		stop
	endelse
endif

if max(freq) lt 10 then x0=freq*1d6 else x0=freq
if min(x0) gt resol then begin ; The echelle diagram works only if freq and spec_reg are vectors such that freq[0]=0... ensure this is not the case

	Nadd=x0[0]/resol ; number of points that need to be added
	fadd=findgen(Nadd-1)*resol
	x0=[fadd[1:*], x0]
	mean_at0=mean(s0[0:5./resol]) ; average over 5 microHz of the spectrum values
	s0=[replicate(mean_at0, n_elements(fadd)), s0]
endif

if ps eq 0 then window, 0, xsize=1400, ysize=800

test=where(s0 ge limit)
if test[0] ne -1 then s0[test]=limit

;stop
echellecorot_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,overplot=overplot,$
			color=color,n_number=n_number, print=print, bar=bar


;stop
;level=0.999
;if smooth_coef ge 30 then smooth_coef=30
;threshold=proba_calc_smooth_signal(1d + 0.05, smooth_coef, level)


end


pro show_ech_diag_v4iter, dir_fit, KIC, ps=ps, extended=extended,smooth_coef=smooth_coef, $
	bckgrd_cor=bckgrd_cor,shifts=shifts, file_out=file_out, dir_spec=dir_spec, print=print, bar=bar,$
	trunc_spec=trunc_spec, pmodes=pmodes,  a_law=a_law, bound=bound, overplot=overplot, $
	tag_label=tag_m, tag_pos=tag_pos, tag_color=tag_color, $
	extralines=extralines, pextralines_col=pextralines_col , style_extralines=style_extralines,  int_coef=int_coef,$
	extra_dots=extra_dots, vert_smooth=vert_smooth, colors_extradots=colors_extradots

if n_elements(ps) eq 0 then ps=0
if n_elements(extended) eq 0 then extended=0 else cmax=extended
if n_elements(smooth_coef) eq 0 then smooth_coef=12
if n_elements(bckgrd_cor) eq 0 then cor=1 else cor=bckgrd_cor
if n_elements(shifts) eq 0 then shifts=0 ; no shift by defaut
if n_elements(file_out) eq 0 then file_out=dir_fit
if n_elements(print) eq 0 then print=1 ; To have Blackk/white data
if n_elements(bar) eq 0 then bar=0 ; by defaut we put the graduation bar
if n_elements(pmodes) eq 0 then pmodes=0

file_fit=file_search(dir_fit+KIC+'_*_data_enhanced.sav')
file_postfit=file_search(dir_fit+KIC+'_postspec.sav')
if file_postfit[0] eq '' then file_postfit=file_search(dir_fit+'post_spect_param.sav')
file_spec=file_search(dir_spec+'*'+KIC+'.sav')

restore, file_fit
restore, file_postfit

freq_tab_t0=reform(stat_synthese_freq[3,*,*])

if pmodes eq 1 then freq_tab_f=dblarr(5, n_elements(stat_synthese_freq[3,0,*])) else freq_tab_f=dblarr(4, n_elements(stat_synthese_freq[3,0,*]))
lmax=n_elements(stat_synthese_freq[3,*,0])-1

if n_elements(overplot) eq 0 then overplot=2 ; show diamonds + lines
l3col='Yellow'
lpcol='Blue'
color=['Orange', 'Red','Dark Gray', l3col]


if pmodes eq 1 then color=[color, lpcol]
if lmax gt 3 then color=[color, replicate(l3col, lmax-2)]

if n_elements(tag_color) eq 0 then tag_color='Black'
if n_elements(extralines) eq 0 then extralines=-1
if n_elements(extra_dots) eq 0 then begin
	extra_dots=-1
	dim_tab=1
endif else 	begin
	dim_tab=size(extra_dots) ; give the dimensionality of the table
	if dim_tab[0] eq 1 then dim_tab=1 else dim_tab=dim_tab[1]
endelse

if n_elements(colors_extradots) eq 0 then begin
	colors_extradots=replicate('Yellow', dim_tab)
endif
color=[color, colors_extradots]

;freq_tab_f[0:lmax, *]=reform(stat_synthese_freq[3,*,*])

for i=0, lmax do freq_tab_f[i, *]=reform(stat_synthese_freq[3,i,*])

if pmodes eq 1 then freq_tab_f[4,0:n_elements(stat_synthese_freq_p[3,1,*])-1]=reform(stat_synthese_freq_p[3,1,*])

freq_tab_tp=freq_tab_f

if extra_dots[0] ne -1 then begin
	if dim_tab eq 1 then begin
		tmp=extra_dots
		extra_dots=dblarr(1, n_elements(tmp))
		extra_dots[0,*]=tmp
	endif

	n_f=n_elements(stat_synthese_freq[3,0,*])
	n_extra=n_elements(extra_dots[0,*])
	ntot=max([n_extra, n_f])

    freq_tab_new=dblarr(5 + dim_tab, ntot) ; we add the column of the extra dots

	freq_tab_new[0:n_elements(freq_tab_tp[*,0])-1, 0:n_elements(freq_tab_tp[0,*])-1]=freq_tab_tp ; we put the old table
	freq_tab_new[5:5+dim_tab-1 , 0:n_extra-1]=extra_dots ; we add the extra dots
	freq_tab_tp=freq_tab_new
endif

if n_elements(bound) eq 0 then begin
	bound=dblarr(lmax+1, 2)
	bound[*,1]=n_elements(freq_tab_tp[0,*])
endif
for el=0, lmax do begin
	Nf_lx=n_elements(where(freq_tab_tp[el,*] ne 0))
	KEPT=freq_tab_tp[el,bound[el,0]: bound[el,1]-1]
	freq_tab_tp[el,*]=0
	freq_tab_tp[el,0: abs(bound[el,1] - bound[el,0])-1]=KEPT
	;if bound[el,0] ne 0 then freq_tab_tp[el,0:bound[el,0]-1]=0
	;if bound[el,1] lt Nf_lx then freq_tab_tp[el,bound[el,1]:*]=0
endfor

l=1
freq_tab_t0[l,*]=freq_tab_t0[l,sort(freq_tab_t0[l,*])]

Nf_l0=n_elements(freq_tab_t0[0,where(freq_tab_t0[0,*] ne 0)])

residu=continuity_condition(reform(freq_tab_t0[0,0:Nf_l0-1]),  n_elements(freq_tab_t0[0,0:Nf_l0-1]), n_elements(freq_tab_t0[0,0:Nf_l0-1]))
Dnu=mean(residu.residual_st)
r=min(freq_tab_t0[0,where(freq_tab_t0[0,*] ne 0)])/Dnu
n0=floor(r)
epsilon=r-n0

n0_l1=floor(min(freq_tab_t0[1,where(freq_tab_t0[1,*] ne 0)])/Dnu)

n0=n0_l1 ; -3
;n0=11d
n_number=ceil(max(freq_tab_t0[1,*])/Dnu) - n0_l1-2 +6 ;+6  ;6

if n_elements(freq) eq 0 then begin
	print, 'BEWARE: We didn t find the original power spectrum (files freq and spec_reg)'
	print, 'Trying to use the file_spec option...'
	if n_elements(file_spec) ne 0 then if file_spec ne '' then restore, file_spec else begin
		print, 'file_spec has not been specified ! Emmergency stop!'
		stop
	endelse
endif
if max(freq) lt 10 then x0=freq*1d6 else x0=freq
if cor eq 1 then begin
	val_des_median=reform(stat_synthese_unsorted[3,*])
	if total(parameters_length[0:7]) lt total(parameters_length[0:8]) then 	noise_param=val_des_median[total(parameters_length[0:7]):total(parameters_length[0:8])-1]
	if parameters_length[7] eq 0 AND parameters_length[8] eq 0 then $
		noise_param=val_des_median[total(parameters_length[0:4]):total(parameters_length[0:5])-1]
	 N=0d & cpt=0d
	 for i=0,1 do begin
		N=N + noise_param[cpt]/(1d  + (noise_param[cpt +1]*x0*1d-3)^noise_param[cpt+2])
		cpt=cpt+3
	 endfor
	 N=N+noise_param[6] ; we add the white noise
	 s0=spec_reg/N
endif else s0=spec_reg

if n_elements(trunc_spec) eq 0 then limit=max(s0) else limit=trunc_spec
test=where(s0 ge limit)
if test[0] ne -1 then s0[test]=limit
if extended eq 0 OR extended eq 1 then begin
if ps eq 0 then window, 1, xsize=1400, ysize=800 else $
	nimp,name=file_out+KIC+'_Ech_diag.eps',/paper,/eps ; Savita's code
;	echellecorot_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,overplot=overplot,color=color,$
;		n_number=n_number, print=print, bar=bar, $
;		tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color, $
;		extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines;, rota=rota
	xtmp=x0
	stmp=s0
	echellecorot_v4,x0,s0,Dnu+shifts,n0,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,overplot=overplot,color=color,$
		n_number=n_number, print=print, bar=bar, $
		tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color, $
		extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, int_coef=int_coef, $
		vert_smooth=vert_smooth;, dots_extra=dots_extra

if ps eq 0 then window, 1, xsize=1400, ysize=800 else $
	nimp,name=file_out+KIC+'_Ech_diag_p.eps',/paper,/eps ; Savita's code
;	echellecorot_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_tp,smooth_coef=smooth_coef,ps=ps,overplot=overplot,color=color,$
;		n_number=n_number, print=print, bar=bar, $
;		tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color, $
;		extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines;, rota=rota
	x1=xtmp
	s1=stmp
	echellecorot_v4,x1,s1,Dnu,n0,uu,freq_table=freq_tab_tp,smooth_coef=smooth_coef,ps=ps,overplot=overplot,color=color,$
		n_number=n_number, print=print, bar=bar, $
		tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color, $
		extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, int_coef=int_coef, $
		vert_smooth=vert_smooth;, dots_extra=dots_extra
endif
if extended ge 2 then begin
	if ps eq 0 then window, 1, xsize=1400, ysize=800 else $
	nimp,name=file_out+KIC+'_Ech_diag_p.eps',/paper,/eps ; Savita's code
		echellecorot_UNFOLD_v4,x0,s0,Dnu+shifts,n0+shifts,uu,freq_table=freq_tab_tp,smooth_coef=smooth_coef,ps=ps,$
			overplot=overplot,color=color, n_number=n_number,cmax=cmax, print=print, bar=bar, a_law=a_law, $
			tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color, $
			extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, vert_smooth=vert_smooth
	if ps ne 0 then fimp
	if ps eq 0 then window, 1, xsize=1400, ysize=800 else $
	nimp,name=file_out+KIC+'_Ech_diag.eps',/paper,/eps ; Savita's code
		echellecorot_UNFOLD_v4,x0,s0,Dnu+shifts,n0,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,$
			overplot=overplot,color=color, n_number=n_number,cmax=cmax, print=print, bar=bar,  a_law=a_law, $
			tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color, $
			extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, vert_smooth=vert_smooth
	if ps ne 0 then fimp
endif


if ps eq 0 then window, 0, xsize=1400, ysize=800

if ps ne 0 then $
	nimp,name=dir_fit+KIC+'_SpecFit.eps',/paper,/eps ; Savita's code
	;e=write_on_ps_on(dir_fit + 'SpecFit')
	if n0*Dnu lt min(x) then x_min=min(x) else x_min=n0*Dnu
	if (n0+n_number)*Dnu gt max(x) then x_max=max(x) else x_max=(n0+n_number)*Dnu
plot, x,smooth(s, smooth_coef),xr=[x_min, x_max],/xst, thick=1
oplot, x, model_spectrum, color=fsc_color('blue')

if ps ne 0 then fimp

if ps ne 0 then begin
	x_max=0
	dn=3 & i=0
	while x_max lt max(x) do begin
		nimp,name=dir_fit+KIC+'_SpecFit-'+strtrim(i,1)+'.eps',/paper,/eps ; Savita's code
		x_min=min(x) + i*dn*Dnu & x_max=min(x) + (i+1)*dn*Dnu
		if x_max gt max(x) then x_max=max(x)
		plot, x,smooth(s, smooth_coef),xr=[x_min, x_max],/xst, thick=1
		oplot, x, model_spectrum, color=fsc_color('blue'), thick=5
		i=i+1
		wait, 0.1
		fimp
	endwhile
	;e=write_on_ps_off('')
endif
end




pro show_ech_diag_v3_Binary, dir_fit1, dir_fit2, ps=ps, extended=extended,smooth_coef=smooth_coef, $
	bckgrd_cor=bckgrd_cor,shifts=shifts, file_out=file_out, file_spec=file_spec, print=print, bar=bar, $
	tag_m=tag_m, tag_pos=tag_pos, tag_color=tag_color, tag_size=tag_size

if n_elements(ps) eq 0 then ps=0
if n_elements(extended) eq 0 then extended=0 else cmax=extended
if n_elements(smooth_coef) eq 0 then smooth_coef=12
if n_elements(bckgrd_cor) eq 0 then cor=1 else cor=bckgrd_cor
if n_elements(shifts) eq 0 then shifts=0 ; no shift by defaut
if n_elements(file_out) eq 0 then file_out2=dir_fit1 else file_out2=file_out
if n_elements(print) eq 0 then print=1 ; To have Blackk/white data
if n_elements(bar) eq 0 then bar=0 ; by defaut we put the graduation bar
if n_elements(tag_color) eq 0 then tag_color='Black'

file_postfit2=file_search(dir_fit2+'post_spect_param.sav')
restore, file_postfit2
x2=x & s2=s & model_spectrum2=model_spectrum ;& freq2=freq  & spec_reg2=spec_reg
file_fit=file_search(dir_fit1+'*_data_enhanced.sav')
file_postfit=file_search(dir_fit1+'post_spect_param.sav')
restore, file_fit
restore, file_postfit


freq_tab_t0=reform(stat_synthese_freq[3,*,*])

;freq_tab_t0[0,0]=0 & freq_tab_t0[0,12]=0 & freq_tab_t0[0,13]=0
;freq_tab_t0[2,*]=0 ; ADHOC
;
;freq_tab_t0=freq_tab_t0[0:1,*] & l0=freq_tab_t0[0,*]
;freq_tab_t0[0,*]=0
;freq_tab_t0[0,0: n_elements(where(l0 ne 0))-1]=l0[where(l0 ne 0)]

;for l=0, n_elements(freq_tab_t0[*,0])-1 do
l=1
freq_tab_t0[l,*]=freq_tab_t0[l,sort(freq_tab_t0[l,*])]


Nf_l0=n_elements(freq_tab_t0[0,where(freq_tab_t0[0,*] ne 0)])
Nf_l1=n_elements(freq_tab_t0[1,where(freq_tab_t0[1,*] ne 0)])

residu=continuity_condition(reform(freq_tab_t0[0,0:Nf_l0-1]),  n_elements(freq_tab_t0[0,0:Nf_l0-1]), n_elements(freq_tab_t0[0,0:Nf_l0-1]))

residu_l1=continuity_condition(reform(freq_tab_t0[1,0:Nf_l1-1]),  n_elements(freq_tab_t0[1,0:Nf_l1-1]), n_elements(freq_tab_t0[1,0:Nf_l1-1]))

Dnu=mean(residu.residual_st)

r=min(freq_tab_t0[0,where(freq_tab_t0[0,*] ne 0)])/Dnu
n0=floor(r)
epsilon=r-n0

overplot=2 ; show diamonds only
color=['Orange', 'Red','Gray', 'Yellow']
;color=['Orange', 'Red','Blue', 'Yellow']

n_number=n_elements(freq_tab_t0[0,*])+6
n0=n0-2

if ps eq 0 then window, 1, xsize=1400, ysize=800 else $
	nimp,name=file_out2+'Ech_diag.eps',/paper,/eps ; Savita's code
	;e=write_on_ps_on(file_fit+'_EchDiag')

if n_elements(freq) eq 0 then begin
	print, 'BEWARE: We didn t find the original power spectrum (files freq and spec_reg)'
	print, 'Trying to use the file_spec option...'
	if n_elements(file_spec) ne 0 then restore, file_spec else begin
		print, 'file_spec has not been specified ! Emmergency stop!'
		stop
	endelse
endif
if max(freq) lt 10 then x0=freq*1d6 else x0=freq


if cor eq 1 then begin
	val_des_median=reform(stat_synthese_unsorted[3,*])
	 noise_param=val_des_median[total(parameters_length[0:4]):total(parameters_length[0:5])-1]
	 N=0d & cpt=0d
	 for i=0,1 do begin
		N=N + noise_param[cpt]/(1d  + (noise_param[cpt +1]*x0*1d-3)^noise_param[cpt+2])
		cpt=cpt+3
	 endfor
	 N=N+noise_param[6] ; we add the white noise
	 N=N+dblarr(n_elements(freq)) ;+ 0.22

	 if n_elements(noise_file) ne 0 then N=backgrd ; use the restored file

	 s0=spec_reg/N
endif else s0=spec_reg

N_local=interpol(N, freq, x)
model_spectrum2_cor=model_spectrum2+N_local


s0[where(s0 ge 4)]=4
if extended eq 0 OR extended eq 1 then $
	echellecorot_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,overplot=overplot,color=color,n_number=n_number, print=print, bar=bar
if extended ge 2 then $
	echellecorot_UNFOLD_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,$
		overplot=overplot,color=color, n_number=n_number,cmax=cmax, print=print, bar=bar, $
		tag_label=tag_m, tag_pos=tag_pos, tag_size=2, tag_color=tag_color
if ps ne 0 then $
	fimp

if ps eq 0 then window, 0, xsize=1400, ysize=800

x_min=3500
if ps ne 0 then $
	nimp,name=dir_fit1+'SpecFit.eps',/paper,/eps ; Savita's code
	;e=write_on_ps_on(dir_fit + 'SpecFit')
	if n0*Dnu lt min(x) then x_min=min(x) else x_min=n0*Dnu
	if (n0+n_number)*Dnu gt max(x) then x_max=max(x) else x_max=(n0+n_number)*Dnu
	bmax=max([model_spectrum, model_spectrum2])*1.1
	bmin=min([model_spectrum])*0.7
plot, x,smooth(s, smooth_coef),xr=[x_min, x_max],/xst, thick=1, /yst, yr=[bmin, bmax]
oplot, x, model_spectrum, color=fsc_color('blue')
oplot, x, model_spectrum2_cor, color=fsc_color('Red') ; with noise


level=0.98
if smooth_coef ge 34 then smooth_coef2=34 else smooth_coef2=smooth_coef
threshold=proba_calc_smooth_signal(1d + 0.05, smooth_coef2, level)

s_red=s/(model_spectrum + model_spectrum2) ; normalized spectra
bmax=max(smooth(s_red[where(x ge x_min AND x le x_max)], smooth_coef, /EDGE_TRUNCATE))
plot, x,smooth(s_red, smooth_coef),xr=[x_min, x_max],/xst, thick=1, title='Residual spectrum', /yst, yr=[0, bmax]
oplot, x,  dblarr(n_elements(x)) + threshold, color=fsc_color('Red'), linestyle=2, thick=2
oplot, x,  dblarr(n_elements(x)) + 1, color=fsc_color('Red'), linestyle=2, thick=2

if ps ne 0 then fimp

if ps ne 0 then begin
	x_max=0
	dn=3 & i=0
	while x_max lt max(x) do begin
		nimp,name=dir_fit1+'SpecFit-'+strtrim(i,1)+'.eps',/paper,/eps ; Savita's code
		x_min=min(x) + i*dn*Dnu & x_max=min(x) + (i+1)*dn*Dnu
		if x_max gt max(x) then x_max=max(x)
		plot, x,smooth(s, smooth_coef2),xr=[x_min, x_max],/xst, thick=1
		oplot, x, model_spectrum, color=fsc_color('blue'), thick=5
		oplot, x, model_spectrum2_cor, color=fsc_color('Red') ; with noise
		i=i+1
		wait, 0.1
		fimp
	endwhile
	;e=write_on_ps_off('')
endif
print, 'yo'
end
