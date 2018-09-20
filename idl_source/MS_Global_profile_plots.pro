
;***************************************************************************************
;******************************* DATA PLOT VERSION 2 ***********************************
;***************************************************************************************
; fonction qui affiche les amplitudes, les largeurs, les hauteurs et les rapport signal sur bruit
;structure de param  : param=[lmax,noise_model_type,first_window,save_image_status]
;noise_model_type : 'Harvey', 'Harvey Free'
pro MS_Global_profile_plots,Stat_Synthese_freq,Stat_Synthese_Height,Stat_synthese_largeur,Stat_Synthese_Amplitude,$
	bruit_local, file

;e=0

ps_writing_plot=1
param=[2,1,30,ps_writing_plot]

j=param[2]

Nmax=n_elements(stat_synthese_Height[0,0,*])

;---- Choix d'une ecriture sur fichier ou sur ecran -----
ind=0
if param[3] eq 1 then begin
	!p.multi=0
	!p.thick=7.
	;e=write_on_ps_on(file+'Heights')
	nimp,name=file +'Heights.eps',/paper,/eps
	ind=ind+1
	print, ' Generating profiles for Amplitudes, Heights, Widths...'
	endif else begin
		!p.multi=[0,2,2]
		window, j, xsize=1024,ysize=600
		j=j+1
endelse

frequency_unit=TeXtoIDL('\mu')+'Hz'
height_unit=TeXtoIDL('ppm^2/\mu')+'Hz'
amplitudes_unit=TeXtoIDL('ppm^2')
angle_unit=TextoIDL('^\circ')
;;;****************************************************************

borne_sup=max(stat_Synthese_Height[4,0,*]+Bruit_local[0,*])*1.2; permet de tronquer la barre d'erreur à 2 sigma si elle depasse de la zone de graphique (3.5 pour HD141820)
print, 'Upper bound for graphs of the Heights: '+strtrim(borne_sup,1)
plotsym, 0, 1., /fill
plot, Stat_Synthese_freq[3,0,*],/NoData,$;/ylog,$
	psym=8,background=255,color=1,yr=[0.0,borne_sup],$
	xtitle='Frequency ('+frequency_unit+')',ytitle='Height ('+height_unit+')',$
	xr=[0.9*min(Stat_Synthese_freq[3,*,*]),1.05*max(Stat_Synthese_freq[3,*,*])],charsize=2.,/yst,/xst
oplot,  Stat_Synthese_freq[3,0,*],Bruit_local[0,*],linestyle=2,thick=5,color=1

for i=0, Nmax-1 do begin
	if (Stat_Synthese_Height[5,0,i]+Bruit_local[0,i]) gt borne_sup then begin
		Stat_Synthese_Height[5,0,i]=borne_sup-Bruit_local[0,i]
		Draw_BoxAndWiskers_no_upper_bound, Stat_Synthese_Height[*,0,i]+Bruit_local[0,i], COLOR='blue',$
			XLOC=Stat_Synthese_freq[3,0,i],1,0
	endif else Draw_BoxAndWiskers, Stat_Synthese_Height[*,0,i]+Bruit_local[0,i], COLOR='blue', XLOC=Stat_Synthese_freq[3,0,i],1,0
endfor

if param[3] eq 1 then begin
	;e=write_on_ps_on(file+'Widths')
	nimp,name= + file +'Widths.eps',/paper,/eps
	ind=ind+1
endif

borne_sup=max(stat_synthese_largeur[4,0,*])*1.2
print, 'Upper bound for graphs of the Widths: '+strtrim(borne_sup,1)
plot, Stat_Synthese_freq[3,0,*],/NoData,$;/ylog,$
	psym=8,background=255,color=1,yr=[0.0,borne_sup],$
	xtitle='Frequency ('+frequency_unit+')',ytitle='FWHM ('+frequency_unit+')',$
	xr=[0.9*min(Stat_Synthese_freq[3,*,*]),1.1*max(Stat_Synthese_freq[3,*,*])],charsize=2.,/yst,/xst

for i=0, Nmax-1 do begin
	if Stat_Synthese_largeur[5,0,i] gt borne_sup then begin
		Stat_Synthese_largeur[5,0,i]=borne_sup
		Draw_BoxAndWiskers_no_upper_bound, Stat_Synthese_largeur[*,0,i], COLOR='blue', XLOC=Stat_Synthese_freq[3,0,i],1,0
	endif else Draw_BoxAndWiskers, Stat_Synthese_largeur[*,0,i], COLOR='blue', XLOC=Stat_Synthese_freq[3,0,i],1,0
endfor

if param[3] eq 1 then begin
	;e=write_on_ps_on(file+'Amplitudes')
	nimp,name= + file +'Amplitudes.eps',/paper,/eps
	ind=ind+1
endif

borne_sup=max(Stat_synthese_Amplitude[4,0,*])*1.3
print, 'Upper bound for graphs of the Amplitudes: '+strtrim(borne_sup,1)
plot, Stat_Synthese_freq[3,0,*],/Nodata,$;/ylog,$
	psym=8,background=255,color=1,yr=[0.0,borne_sup],$
	xtitle='Frequency ('+frequency_unit+')',ytitle='Amplitude (ppm)',$
	xr=[0.9*min(Stat_Synthese_freq[3,*,*]),1.1*max(Stat_Synthese_freq[3,*,*])],charsize=2.,/xst,/yst
;oplot, Stat_Synthese_freq[7,0,*],Stat_Synthese_Amplitude[7,0,*],color=FSC_Color('Green'),psym=8
for i=0, Nmax-1 do begin
	Draw_BoxAndWiskers, Stat_Synthese_Amplitude[*,0,i], COLOR='blue', XLOC=Stat_Synthese_freq[3,0,i],1,0
endfor

fimp
;if e eq -1 then e=write_on_ps_off(file) else $
	print, ' Profiles for Amplitudes, Heights, Widths successfully generated'

end
