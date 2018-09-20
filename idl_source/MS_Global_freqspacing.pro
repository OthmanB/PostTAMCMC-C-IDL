; gere les fichiers de synthese contenant les calculs
; sur la base de cumulated probability density function
; Affiche la d01, d02, la grande separation
; Les variables generés suivent la meme syntaxe que celle d'entrée
function MS_Global_freqspacing,file,rep_out,Nmax, lmax, show=show

;e=write_on_ps_on('a')
nimp,name='dummy.eps',/paper,/eps
frequency_unit=TeXtoIDL('\mu')+'Hz'
height_unit=TeXtoIDL('ppm^2/\mu')+'Hz'
amplitudes_unit=TeXtoIDL('ppm^2')
angle_unit=TextoIDL('^\circ')
fimp
;e=write_on_ps_off('a')

tab_critere=[2.25,16,50,84,97.75]

restore,file

Dnu=dblarr(8,Lmax+1,Nmax-1)  ;
d01=dblarr(8,2,Nmax-1)
dd01=dblarr(8,2,Nmax-2) ; attention : dd01 est centré en i+1 par rapport à la table de Dnu
dd01_red=dblarr(8,2,Nmax-2) ; valeurs reduite ie : dd01(n)/Dnu(n)=dd01(i)/Dnu(i+1)
d02=dblarr(8,Nmax)

; ******* Calcul de la d02 est des erreurs associés ********
if Lmax ge 2 then begin
d02[3,*]=stat_synthese_freq[3,0,0:Nmax-1]-stat_synthese_freq[3,2,0:Nmax-1]
d02[7,*]=stat_synthese_freq[7,0,0:Nmax-1]-stat_synthese_freq[7,2,0:Nmax-1]
for j=0,2 do begin	; ------ les erreurs -----
	d02[j,*]=sqrt((stat_synthese_freq[j,0,0:Nmax-1]-stat_synthese_freq[3,0,0:Nmax-1])^2+$
		(stat_synthese_freq[j,2,0:Nmax-1]-stat_synthese_freq[3,2,0:Nmax-1])^2)
	d02[j+4,*]=sqrt((stat_synthese_freq[j+4,0,*]-stat_synthese_freq[3,0,0:Nmax-1])^2+$
		(stat_synthese_freq[j+4,2,0:Nmax-1]-stat_synthese_freq[3,2,0:Nmax-1])^2)
	d02[j,*]=d02[3,0:Nmax-1]-d02[j,0:Nmax-1] ; bornes inferieurs
	d02[j+4,*]=d02[3,0:Nmax-1]+d02[j+4,0:Nmax-1] ; bornes superieurs
endfor
endif

; ********* Calcul de la Grande separation *************
for i=0,Nmax-2 do begin
	for l=0,Lmax do begin ; ----- pur chaque degrée l -------
	Dnu[3,l,i]=stat_synthese_freq[3,l,i+1]-stat_synthese_freq[3,l,i]
	Dnu[7,l,i]=stat_synthese_freq[7,l,i+1]-stat_synthese_freq[7,l,i]
		for j=0,2 do begin ; ------ les errerus --------
			Dnu[j,l,i]=sqrt((stat_synthese_freq[j,l,i+1]-stat_synthese_freq[3,l,i+1])^2+$
				(stat_synthese_freq[j,l,i]-stat_synthese_freq[3,l,i])^2)
			Dnu[j+4,l,i]=sqrt((stat_synthese_freq[j+4,l,i+1]-stat_synthese_freq[3,l,i+1])^2+$
				(stat_synthese_freq[j+4,l,i]-stat_synthese_freq[3,l,i])^2)
			Dnu[j,l,i]=Dnu[3,l,i]-Dnu[j,l,i] ; bornes inferieurs
			Dnu[j+4,l,i]=Dnu[3,l,i]+Dnu[j+4,l,i] ; bornes superieurs
		endfor
	;	print, Dnu[*,l,i]
	endfor
endfor

; ******* Calculs des d01 ********
if n_elements(stat_synthese_freq[3,0,*]) eq Nmax then begin
	for i=0,Nmax-2 do begin
	 ; ---- les deux manieres de calculer la d01 ----
		tmp=0.5*(stat_synthese_freq[3,0,i]+stat_synthese_freq[3,0,i+1])-stat_synthese_freq[3,1,i]
		tmp_2=0.5*(stat_synthese_freq[3,0,i]+stat_synthese_freq[3,0,i+1])-stat_synthese_freq[3,1,i+1]

		if tmp lt mean(Dnu[3,*,*])/2 AND tmp gt -mean(Dnu[3,*,*])/2 then begin
			d01[3,0,i]=tmp
			d01[7,0,i]=0.5*(stat_synthese_freq[7,0,i]+stat_synthese_freq[7,0,i+1])-stat_synthese_freq[7,1,i]
			for j=0,2 do begin ; ------ les erreurs --------
				d01[j,0,i]=sqrt(0.25*((stat_synthese_freq[j,0,i]-stat_synthese_freq[3,0,i])^2+$
					(stat_synthese_freq[j,0,i+1]-stat_synthese_freq[3,0,i+1])^2)+$
					(stat_synthese_freq[j,1,i]-stat_synthese_freq[3,1,i])^2)
				d01[j+4,0,i]=sqrt(0.25*((stat_synthese_freq[j+4,0,i]-stat_synthese_freq[3,0,i])^2+$
					(stat_synthese_freq[j+4,0,i+1]-stat_synthese_freq[3,0,i+1])^2)+$
					(stat_synthese_freq[j+4,1,i]-stat_synthese_freq[3,1,i])^2)

				d01[j,0,i]=d01[3,0,i]-d01[j,0,i] ; bornes inferieurs
				d01[j+4,0,i]=d01[3,0,i]+d01[j+4,0,i] ; bornes superieurs
			endfor
		endif
		if tmp_2 lt mean(Dnu[3,*,*])/2  AND tmp_2 gt -mean(Dnu[3,*,*])/2 then begin
			d01[3,0,i]=tmp_2
			d01[7,0,i]=0.5*(stat_synthese_freq[7,0,i]+stat_synthese_freq[7,0,i+1])-stat_synthese_freq[7,1,i+1]
			for j=0,2 do begin ; ------ les erreurs --------
				d01[j,0,i]=sqrt(0.25*((stat_synthese_freq[j,0,i]-stat_synthese_freq[3,0,i])^2+$
					(stat_synthese_freq[j,0,i+1]-stat_synthese_freq[3,0,i+1])^2)+$
					(stat_synthese_freq[j,1,i]-stat_synthese_freq[3,1,i])^2)
				d01[j+4,0,i]=sqrt(0.25*((stat_synthese_freq[j+4,0,i]-stat_synthese_freq[3,0,i])^2+$
					(stat_synthese_freq[j+4,0,i+1]-stat_synthese_freq[3,0,i+1])^2)+$
					(stat_synthese_freq[j+4,1,i+1]-stat_synthese_freq[3,1,i+1])^2)

				d01[j,0,i]=d01[3,0,i]-d01[j,0,i] ; bornes inferieurs
				d01[j+4,0,i]=d01[3,0,i]+d01[j+4,0,i] ; bornes superieurs
			endfor
		endif
		if tmp ge mean(Dnu[3,*,*])/2 AND tmp_2 ge mean(Dnu[3,*,*])/2 then print, 'Probleme avec d01 !!'

		; * - * - * - * - * - * - * - * - * - * - * - * - *
		tmp=-0.5*(stat_synthese_freq[3,1,i]+stat_synthese_freq[3,1,i+1])+stat_synthese_freq[3,0,i]
		tmp_2=-0.5*(stat_synthese_freq[3,1,i]+stat_synthese_freq[3,1,i+1])+stat_synthese_freq[3,0,i+1]

		if tmp lt mean(Dnu[3,*,*])/2  AND tmp gt -mean(Dnu[3,*,*])/2 then begin
			d01[3,1,i]=tmp
			d01[7,1,i]=-0.5*(stat_synthese_freq[7,1,i]+stat_synthese_freq[7,1,i+1])+stat_synthese_freq[7,0,i]
			for j=0,2 do begin ; ------ les erreurs --------
				d01[j,1,i]=sqrt(0.25*((stat_synthese_freq[j,1,i]-stat_synthese_freq[3,1,i])^2+$
					(stat_synthese_freq[j,1,i+1]-stat_synthese_freq[3,1,i+1])^2)+$
					(stat_synthese_freq[j,0,i]-stat_synthese_freq[3,0,i])^2)
				d01[j+4,1,i]=sqrt(0.25*((stat_synthese_freq[j+4,1,i]-stat_synthese_freq[3,1,i])^2+$
					(stat_synthese_freq[j+4,1,i+1]-stat_synthese_freq[3,1,i+1])^2)+$
					(stat_synthese_freq[j+4,0,i]-stat_synthese_freq[3,0,i])^2)

				d01[j,1,i]=d01[3,1,i]-d01[j,1,i] ; bornes inferieurs
				d01[j+4,1,i]=d01[3,1,i]+d01[j+4,1,i] ; bornes superieurs
			endfor
		endif

		if tmp_2 lt mean(Dnu[3,*,*])/2  AND tmp_2 gt -mean(Dnu[3,*,*])/2 then begin
			d01[3,1,i]=tmp_2
			d01[7,1,i]=-0.5*(stat_synthese_freq[7,1,i]+stat_synthese_freq[7,1,i+1])+stat_synthese_freq[7,0,i+1]
			for j=0,2 do begin ; ------ les errerus --------
				d01[j,1,i]=sqrt(0.25*((stat_synthese_freq[j,1,i]-stat_synthese_freq[3,1,i])^2+$
					(stat_synthese_freq[j,1,i+1]-stat_synthese_freq[3,1,i+1])^2)+$
					(stat_synthese_freq[j,0,i+1]-stat_synthese_freq[3,0,i+1])^2)
				d01[j+4,1,i]=sqrt(0.25*((stat_synthese_freq[j+4,1,i]-stat_synthese_freq[3,1,i])^2+$
					(stat_synthese_freq[j+4,1,i+1]-stat_synthese_freq[3,1,i+1])^2)+$
					(stat_synthese_freq[j+4,0,i+1]-stat_synthese_freq[3,0,i+1])^2)
	
				d01[j,1,i]=d01[3,1,i]-d01[j,1,i] ; bornes inferieurs
				d01[j+4,1,i]=d01[3,1,i]+d01[j+4,1,i] ; bornes superieurs
			endfor
		endif

		if tmp ge mean(Dnu[3,*,*])/2 AND tmp_2 ge mean(Dnu[3,*,*])/2 then print, 'Probleme avec d01 !!'

	endfor

;********************************************************************
;********** dd01 (Roxburgh & Vorontsov 2003 A&A 411, 215) ***********
;********************************************************************

	for i=0,Nmax-3 do begin
		 ; ---- les deux manieres de calculer la dd01 selon l'identification ----
		tmp=(stat_synthese_freq[3,0,i]-4*stat_synthese_freq[3,1,i]+6*stat_synthese_freq[3,0,i+1]-$
			4*stat_synthese_freq[3,1,i+1]+stat_synthese_freq[3,0,i+2])/8
		tmp_2=(stat_synthese_freq[3,0,i]-4*stat_synthese_freq[3,1,i+1]+6*stat_synthese_freq[3,0,i+1]-$
			4*stat_synthese_freq[3,1,i+2]+stat_synthese_freq[3,0,i+2])/8

		if tmp lt mean(Dnu[3,*,*])/2 AND tmp gt -mean(Dnu[3,*,*])/2 then begin
			dd01[3,0,i]=tmp	
			dd01[7,0,i]=(stat_synthese_freq[7,0,i]-4*stat_synthese_freq[7,1,i]+6*stat_synthese_freq[7,0,i+1]-$
				4*stat_synthese_freq[7,1,i+1]+stat_synthese_freq[7,0,i+2])/8

			for j=0,2 do begin ; ------ les erreurs --------
				dd01[j,0,i]=sqrt((stat_synthese_freq[j,0,i]-stat_synthese_freq[3,0,i])^2+$
					16*(stat_synthese_freq[j,1,i]-stat_synthese_freq[3,1,i])^2+$
					36*(stat_synthese_freq[j,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
					16*(stat_synthese_freq[j,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
					(stat_synthese_freq[j,0,i+2]-stat_synthese_freq[3,0,i+2])^2)/8

				dd01[j+4,0,i]=sqrt((stat_synthese_freq[j+4,0,i]-stat_synthese_freq[3,0,i])^2+$
					16*(stat_synthese_freq[j+4,1,i]-stat_synthese_freq[3,1,i])^2+$
					36*(stat_synthese_freq[j+4,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
					16*(stat_synthese_freq[j+4,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
					(stat_synthese_freq[j+4,0,i+2]-stat_synthese_freq[3,0,i+2])^2)/8
	
				dd01[j,0,i]=d01[3,0,i]-d01[j,0,i] ; bornes inferieurs
				dd01[j+4,0,i]=d01[3,0,i]+d01[j+4,0,i] ; bornes superieurs
			endfor
		endif
		if tmp_2 lt mean(Dnu[3,*,*])/2  AND tmp_2 gt -mean(Dnu[3,*,*])/2 then begin
			dd01[3,0,i]=tmp_2
			dd01[7,0,i]=(stat_synthese_freq[7,0,i]-4*stat_synthese_freq[7,1,i+1]+6*stat_synthese_freq[7,0,i+1]-$
				4*stat_synthese_freq[7,1,i+2]+stat_synthese_freq[7,0,i+2])/8
			for j=0,2 do begin ; ------ les erreurs --------
				dd01[j,0,i]=sqrt((stat_synthese_freq[j,0,i]-stat_synthese_freq[3,0,i])^2+$
					16*(stat_synthese_freq[j,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
					36*(stat_synthese_freq[j,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
					16*(stat_synthese_freq[j,1,i+2]-stat_synthese_freq[3,1,i+2])^2+$
					(stat_synthese_freq[j,0,i+2]-stat_synthese_freq[3,0,i+2])^2)/8
				dd01[j+4,0,i]=sqrt((stat_synthese_freq[j+4,0,i]-stat_synthese_freq[3,0,i])^2+$
					16*(stat_synthese_freq[j+4,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
					36*(stat_synthese_freq[j+4,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
					16*(stat_synthese_freq[j+4,1,i+2]-stat_synthese_freq[3,1,i+2])^2+$
					(stat_synthese_freq[j+4,0,i+2]-stat_synthese_freq[3,0,i+2])^2)/8
		
				dd01[j,0,i]=dd01[3,0,i]-dd01[j,0,i] ; bornes inferieurs
				dd01[j+4,0,i]=dd01[3,0,i]+dd01[j+4,0,i] ; bornes superieurs
			endfor
		endif

		if tmp ge mean(Dnu[3,*,*])/2 AND tmp_2 ge mean(Dnu[3,*,*])/2 then print, 'Probleme avec d01 !!'

	endfor
; -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

	for i=0,Nmax-3 do begin
	 ; ---- les deux manieres de calculer la d01 selon l'identification----
		tmp=(stat_synthese_freq[3,1,i]-4*stat_synthese_freq[3,0,i+1]+6*stat_synthese_freq[3,1,i+1]-$
			4*stat_synthese_freq[3,0,i+2]+stat_synthese_freq[3,1,i+2])/8
		tmp_2=-(stat_synthese_freq[3,1,i]-4*stat_synthese_freq[3,0,i]+6*stat_synthese_freq[3,1,i+1]-$
			4*stat_synthese_freq[3,0,i+1]+stat_synthese_freq[3,1,i+2])/8

		if tmp lt mean(Dnu[3,*,*])/2 AND tmp gt -mean(Dnu[3,*,*])/2 then begin
			dd01[3,1,i]=tmp
			dd01[7,1,i]=(stat_synthese_freq[7,1,i]-4*stat_synthese_freq[7,0,i+1]+6*stat_synthese_freq[7,1,i+1]-$
				4*stat_synthese_freq[7,0,i+2]+stat_synthese_freq[7,1,i+2])/8

			for j=0,2 do begin ; ------ les erreurs --------
				dd01[j,1,i]=sqrt((stat_synthese_freq[j,1,i]-stat_synthese_freq[3,1,i])^2+$
					16*(stat_synthese_freq[j,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
					36*(stat_synthese_freq[j,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
					16*(stat_synthese_freq[j,0,i+2]-stat_synthese_freq[3,0,i+2])^2+$
					(stat_synthese_freq[j,1,i+2]-stat_synthese_freq[3,1,i+2])^2)/8

				dd01[j+4,1,i]=sqrt((stat_synthese_freq[j+4,1,i]-stat_synthese_freq[3,1,i])^2+$
					16*(stat_synthese_freq[j+4,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
					36*(stat_synthese_freq[j+4,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
					16*(stat_synthese_freq[j+4,0,i+2]-stat_synthese_freq[3,0,i+2])^2+$
					(stat_synthese_freq[j+4,1,i+2]-stat_synthese_freq[3,1,i+2])^2)/8

				dd01[j,1,i]=dd01[3,1,i]-dd01[j,1,i] ; bornes inferieurs
				dd01[j+4,1,i]=dd01[3,1,i]+dd01[j+4,1,i] ; bornes superieurs
			endfor
		endif
		if tmp_2 lt mean(Dnu[3,*,*])/2  AND tmp_2 gt -mean(Dnu[3,*,*])/2 then begin
			dd01[3,1,i]=tmp_2
			dd01[7,1,i]=-(stat_synthese_freq[7,1,i]-4*stat_synthese_freq[7,0,i]+6*stat_synthese_freq[7,1,i+1]-$
				4*stat_synthese_freq[7,0,i+1]+stat_synthese_freq[7,1,i+2])/8
				for j=0,2 do begin ; ------ les erreurs --------
					dd01[j,1,i]=sqrt((stat_synthese_freq[j,1,i]-stat_synthese_freq[3,1,i])^2+$
						16*(stat_synthese_freq[j,0,i]-stat_synthese_freq[3,0,i])^2+$
						36*(stat_synthese_freq[j,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
						16*(stat_synthese_freq[j,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
						(stat_synthese_freq[j,1,i+2]-stat_synthese_freq[3,1,i+2])^2)/8
					dd01[j+4,1,i]=sqrt((stat_synthese_freq[j+4,1,i]-stat_synthese_freq[3,1,i])^2+$
						16*(stat_synthese_freq[j+4,0,i]-stat_synthese_freq[3,0,i])^2+$
						36*(stat_synthese_freq[j+4,1,i+1]-stat_synthese_freq[3,1,i+1])^2+$
						16*(stat_synthese_freq[j+4,0,i+1]-stat_synthese_freq[3,0,i+1])^2+$
						(stat_synthese_freq[j+4,1,i+2]-stat_synthese_freq[3,1,i+2])^2)/8

					dd01[j,1,i]=dd01[3,1,i]-dd01[j,1,i] ; bornes inferieurs
					dd01[j+4,1,i]=dd01[3,1,i]+dd01[j+4,1,i] ; bornes superieurs
				endfor
		endif

		if tmp ge mean(Dnu[3,*,*])/2 AND tmp_2 ge mean(Dnu[3,*,*])/2 then print, 'Probleme avec d01 !!'

	endfor

	;*******************************************************
	;******** Petite separation / grande separation ********
	;*******************************************************

	for i=0,Nmax-3 do for l=0,1 do begin
		dd01_red[3,l,i]=dd01[3,l,i]/Dnu[3,0,i+1]
	
		for j=0,2 do begin ; ------ les erreurs --------
			dd01_red[j,l,i]=sqrt((dd01[3,l,i]^2/Dnu[3,l,i+1]^4)*(Dnu[j,l,i+1]-Dnu[3,l,i+1])^2+$
				(dd01[j,l,i]-dd01[3,l,i])^2/Dnu[3,l,i+1]^2)
			dd01_red[j+4,l,i]=sqrt((dd01[3,l,i]^2/Dnu[3,l,i+1]^4)*(Dnu[j+4,l,i+1]-Dnu[3,l,i+1])^2+$
				(dd01[j+4,l,i]-dd01[3,l,i])^2/Dnu[3,l,i+1]^2)

			dd01_red[j,l,i]=dd01_red[3,l,i]-dd01_red[j,l,i] ; bornes inferieurs
			dd01_red[j+4,l,i]=dd01_red[3,l,i]+dd01_red[j+4,l,i] ; bornes superieurs
		endfor
	endfor
endif else begin
	print, 'd01 not calculated: The number of l=1 modes is not matching the number of l=0 modes. This might be normal if you fit a star with mixed modes'
endelse

; *********************** ECRITURE SUR FICHIER ET SUR ECRAN ***************************

Dnu_m=dblarr(Lmax+1)
Err_Dnu=dblarr(Lmax+1)

d01_m=dblarr(2)
Err_d01=dblarr(2)

dd01_m=dblarr(2)
Err_dd01=dblarr(2)

d02_m=0.
err_d02=0.

print, '*** Mean Large separation and epsilon from linear fit ***'
if n_elements(stat_synthese_freq[3,0,*]) eq Nmax then begin
	for l=0,Lmax do begin
		errors=sqrt( (stat_synthese_freq[3,l,*]-stat_synthese_freq[2,l,*])^2 + (stat_synthese_freq[4,l,*]-stat_synthese_freq[3,l,*])^2 ) / sqrt(2)
		if total(errors) ge 1d-3 then c=linfit(findgen(n_elements(stat_synthese_freq[3,l,*])), reform(stat_synthese_freq[3,l,*]), MEASURE_ERRORS=errors, sigma=err_asslaw)
		if total(errors) lt 1d-3 then c=linfit(findgen(n_elements(stat_synthese_freq[3,l,*])), reform(stat_synthese_freq[3,l,*]))
	
		;errors=sqrt( (Dnu[3,l,*]-Dnu[2,l,*])^2 + (Dnu[4,l,*]-Dnu[3,l,*])^2 ) / sqrt(2)
		;c=linfit(reform(Dnu[3,l,*]), findgen(n_elements(Dnu[3,l,*])), MEASURE_ERRORS=errors, sigma=err_asslaw)
		Dnu_m[l]=c[1] ;mean(Dnu[3,l,*])
		Err_Dnu[l]=err_asslaw[1]*sqrt(n_elements(stat_synthese_freq[3,l,*])) ; come back to the error for each point... not the asymtpotic one ;sqrt(variance(Dnu[3,l,*]))/sqrt(n_elements(Dnu[0,0,*]))
		if l eq 0 then begin
			epsilon0=c[0]/c[1] - fix(c[0]/c[1])
			err_epsilon0=sqrt(err_asslaw[0]^2/c[1]^2 + err_asslaw[1]^2/c[0]^4)*sqrt(n_elements(stat_synthese_freq[3,l,*]))
		endif
		print, 'Dnu l='+strtrim(l,2), Dnu_m[l],'+/-', Err_Dnu[l]
	endfor
		print, 'epsilon l=0 ', epsilon0,'+/-', Err_epsilon0

	;stop
	print, '*** Average ***'
	Dnu_mm=mean(Dnu_m) ;mean(Dnu[3,0:1,*])
	Err_Dnu_mm=total(Err_Dnu^2)/sqrt(n_elements(err_dnu))  ;sqrt(variance(Dnu[3,0:1,*]))/sqrt(n_elements(Dnu[0,0:1,*]))
	print, Dnu_mm,'+/-',Err_Dnu_mm

	print, '*** Mean small separation d01 ***'
	for p=0,1 do begin
		d01_m[p]=mean(d01[3,p,*])
		Err_d01[p]=sqrt(variance(d01[3,p,*]))/sqrt(n_elements(d01[0,0,*]))
		print, 'l='+strtrim(p,2), d01_m[p], '+/-',Err_d01[p]
	endfor
	print, '*** Average d01 ***'
	d01_mm=mean(d01[3,*,*])
	Err_d01_mm=sqrt(variance(d01[3,*,*]))/sqrt(n_elements(d01[0,*,*]))
	print, d01_mm,'+/-',Err_d01_mm

	;print, '*** petites separations moyennes dd01 (Roxburgh)***'
	;for p=0,1 do begin
	;	dd01_m[p]=mean(dd01[3,p,*])
	;	Err_dd01[p]=sqrt(variance(dd01[3,p,*]))/sqrt(n_elements(dd01[0,0,*]))
	;	print, p, dd01_m[p], '+/-',Err_dd01[p]
	;endfor
	;print, '*** moyenne dd01 global (Roxburgh) ***'
	;dd01_mm=mean(dd01[3,*,*])
	;Err_dd01_mm=sqrt(variance(dd01[3,*,*]))/sqrt(n_elements(dd01[0,*,*]))
	;print, dd01_mm,'+/-',Err_dd01_mm
endif else begin
l=0
		errors=sqrt( (stat_synthese_freq[3,l,*]-stat_synthese_freq[2,l,*])^2 + (stat_synthese_freq[4,l,*]-stat_synthese_freq[3,l,*])^2 ) / sqrt(2)
		if total(errors) ge 1d-3 then c=linfit(findgen(n_elements(stat_synthese_freq[3,l,*])), reform(stat_synthese_freq[3,l,*]), MEASURE_ERRORS=errors, sigma=err_asslaw)
		if total(errors) lt 1d-3 then c=linfit(findgen(n_elements(stat_synthese_freq[3,l,*])), reform(stat_synthese_freq[3,l,*]))
	
		;errors=sqrt( (Dnu[3,l,*]-Dnu[2,l,*])^2 + (Dnu[4,l,*]-Dnu[3,l,*])^2 ) / sqrt(2)
		;c=linfit(reform(Dnu[3,l,*]), findgen(n_elements(Dnu[3,l,*])), MEASURE_ERRORS=errors, sigma=err_asslaw)
		Dnu_m[l]=c[1] ;mean(Dnu[3,l,*])
		Err_Dnu[l]=err_asslaw[1]*sqrt(n_elements(stat_synthese_freq[3,l,*])) ; come back to the error for each point... not the asymtpotic one ;sqrt(variance(Dnu[3,l,*]))/sqrt(n_elements(Dnu[0,0,*]))
		if l eq 0 then begin
			epsilon0=c[0]/c[1] - fix(c[0]/c[1])
			err_epsilon0=sqrt(err_asslaw[0]^2/c[1]^2 + err_asslaw[1]^2/c[0]^4)*sqrt(n_elements(stat_synthese_freq[3,l,*]))
		endif
		print, 'Dnu l='+strtrim(l,2), Dnu_m[l],'+/-', Err_Dnu[l]
		print, 'epsilon l=0 ', epsilon0,'+/-', Err_epsilon0
endelse
if Lmax ge 2 then begin
	print, '*** Mean small separation d02 ***'
	d02_m=mean(d02[3,0:Nmax-1])
	err_d02=sqrt(variance(d02[3,0:Nmax-1]))/n_elements(d02[3,0:Nmax-1])
endif
print, d02_m,'+/-',Err_d02

;**************** Ecriture sur fichier Texte de Sortie ********************
save,tab_critere,d01,dd01,dd01_red,Dnu,d02,filename=rep_out+'freq_spacings.sav'

title=['****** Dnu ****** ','****** d01 ******']
if Lmax ge 2 then title=[title,'****** d02 ******']

under_title=['l',strtrim(tab_critere[0],1)+'%',strtrim(tab_critere[1],1)+'%',$
        strtrim(tab_critere[2],1)+'%',strtrim(tab_critere[3],1)+'%',strtrim(tab_critere[4],1)+'%']

openw,5,rep_out+'freq_spacings.txt'
printf,5,format='(A38)',''

printf,5,format='(A38)',''
printf,5,format='(A38)',title[0]
printf,5, format='(A38)', ' ---- Individual parameters ----'


if n_elements(stat_synthese_freq[3,0,*]) eq Nmax then begin
	lmax0=lmax
endif else begin
	lmax0=0
endelse

; on ecrit dans l'ordre l'ensemble decile/quartile/mediane/quartile/decile/maximum
for ll=0,n_elements(Dnu[3,0:lmax0,0])-1 do for i=0,n_elements(Dnu[0,0,*])-1 do begin
        printf,5,format='(I4,5F14.6)',ll,Dnu[1:5,ll,i]
endfor
printf,5, format='(A38)', ' ---- Mean of Dnu using l=0 OR l=1 ----'
for ll=0,lmax0 do printf,5,format='(I4,F14.6,A4,F14.6)',ll,Dnu_m[ll],'+/-',Err_Dnu[ll]
printf,5, format='(A38)', ' ---- Mean of Dnu using l=0 AND l=1 ----'
printf,5,format='(F14.6,A4,F14.6)',Dnu_mm,'+/-',Err_Dnu_mm

printf,5,format='(A38)',''
printf,5,format='(A38)',title[1]
printf,5, format='(A38)', ' ---- Individual parameters ----'
if lmax0 ne 0 then begin
	for p=0,n_elements(d01[3,*,0])-1 do for i=0,n_elements(d01[0,0,*])-1 do begin
	        printf,5,format='(I4,5F14.6)',p,d01[1:5,p,i]
	endfor
	printf,5, format='(A60)', '---- Mean on d01 ----'
	for p=0,1 do printf,5,format='(F14.6,A4,F14.6)',d01_m[p],'+/-',Err_d01[p]
	printf,5, format='(A60)', '---- Moyenne global ----'
	printf,5,format='(F14.6,A4,F14.6)',d01_mm,'+/-',Err_d01_mm
endif else begin
	printf, 5, 'd01 not calculated because the number of l=0 and l=1 does not match'
endelse
if Lmax ge 2 then begin
	printf,5,format='(A38)',''
	printf,5,format='(A38)',title[2]
	printf,5, format='(A38)', ' ---- Parametres individuels ----'
	for i=0,n_elements(d02[0,*])-1 do begin
  	      printf,5,format='(I4,5F14.6)',p,d02[1:5,i]
	endfor
printf,5, format='(A60)', '---- d02 global mean ----'
printf,5,format='(F14.6,A4,F14.6)',d02_m,'+/-',Err_d02
endif

close,5


;*********** AFFICHAGE *******
if n_elements(show) eq 0 then show=1
if show eq 1 then begin
;----------- GRANDE SEPARATION --------
file_out=rep_out+'Dnu'
;e=write_on_ps_on(file_out)
nimp,name=file_out,/paper,/eps
plot, Stat_Synthese_freq[3,0,*],/NoData,$;/ylog,$
	psym=8,background=fsc_color('White'),color=1,$
	xtitle='Frequency ('+frequency_unit+')',ytitle='Large separation ('+frequency_unit+')',$
	title='Large Separation : Black => l=0, Blue => l=1' ,$
	xr=[0.9*min(Stat_Synthese_freq[*,0,*]),1.1*max(Stat_Synthese_freq[*,0,*])],charsize=2.,/yst,/xst,$
	yr=[min(Dnu[1:5,*,*]),max(Dnu[1:5,*,*])]

col_table=['BLACK','BLUE','GREEN']
if lmax0 ne 0 then begin
	for i=0, Nmax-2 do begin
		for l=0,Lmax-1 do begin
			Draw_BoxAndWiskers, Dnu[*,l,i], COLOR=col_table[l], XLOC=Stat_Synthese_freq[3,l,i],1
		endfor
	endfor
endif else begin
	l=0
	for i=0, Nmax-2 do begin
			Draw_BoxAndWiskers, Dnu[*,l,i], COLOR=col_table[l], XLOC=Stat_Synthese_freq[3,l,i],1
	endfor
endelse

fimp
;e=write_on_ps_off(file_out)
;------------ PETITE SEPARATION d02 --------
file_out=rep_out+'d02'
;e=write_on_ps_on(file_out)
nimp,name=file_out,/paper,/eps
plot, Stat_Synthese_freq[3,0,0:Nmax-1],/NoData,$;/ylog,$
	psym=8,background=255,color=1,$
	xtitle='Frequency ('+frequency_unit+')',ytitle='Large separation ('+frequency_unit+')',$
	title='Small separation d02' ,$
	xr=[0.9*min(Stat_Synthese_freq[*,0,0:Nmax-1]),1.1*max(Stat_Synthese_freq[*,0,0:Nmax-1])],charsize=2.,/yst,/xst,$
	yr=[min(d02[1:5,*,0:Nmax-1]),max(d02[1:5,*,0:Nmax-1])]

col_table=['BLACK','BLUE','GREEN']
for i=0, Nmax-2 do begin
		Draw_BoxAndWiskers, d02[*,i], COLOR=col_table[1], XLOC=Stat_Synthese_freq[3,0,i],1
endfor
fimp
;e=write_on_ps_off(file_out)

;---------- PETITE SEPARATION d01 ---------
if lmax0 ne 0 then begin
	file_out=rep_out+'d01'
	;e=write_on_ps_on(file_out)
	nimp,name=file_out,/paper,/eps
	plot, Stat_Synthese_freq[3,0,*],/NoData,$;/ylog,$
		psym=8,background=255,color=1,$
		xtitle='Frequency ('+frequency_unit+')',ytitle='Small separation ('+frequency_unit+')',$
		title='Small separation d01: Green => Def 1 (2*l=0 & 1*l=1), Blue => Def 2 (2*l=1 & 1*l=0)' ,$
		xr=[0.9*min(Stat_Synthese_freq[*,0,*]),1.1*max(Stat_Synthese_freq[*,0,*])],charsize=1.,/yst,/xst,$
		yr=[min(d01[1:5,*,*]),max(d01[1:5,*,*])]

	col_table=['BLUE','GREEN']
	for i=0, Nmax-2 do begin
		for l=0,1 do begin
			Draw_BoxAndWiskers, d01[*,l,i], COLOR=col_table[l], XLOC=Stat_Synthese_freq[3,l,i],1
		endfor
	endfor
	fimp
;e=write_on_ps_off(file_out)
endif
endif

arr_out=dblarr(2, n_elements(Dnu_m))
arr_out[0,*]=Dnu_m
arr_out[1,*]=err_Dnu
return, arr_out
end
