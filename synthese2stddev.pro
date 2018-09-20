; procedure lisant les fichiers de synthese (.sav)
; de la loi assymptotique et des tables de frequences
; et construit des tables digestes pour les theoriciens

pro synthese2stddev,file1,file2,file_out

restore, file1

err_m_f=stat_synthese_freq[3,*,*]-stat_synthese_freq[2,*,*]
err_p_f=stat_synthese_freq[4,*,*]-stat_synthese_freq[3,*,*]

err_2m_f=stat_synthese_freq[3,*,*]-stat_synthese_freq[1,*,*]
err_2p_f=stat_synthese_freq[5,*,*]-stat_synthese_freq[3,*,*]

err_m_l=stat_synthese_width[3,*,*]-stat_synthese_width[2,*,*]
err_p_l=stat_synthese_width[4,*,*]-stat_synthese_width[3,*,*]

err_2m_l=stat_synthese_width[3,*,*]-stat_synthese_width[1,*,*]
err_2p_l=stat_synthese_width[5,*,*]-stat_synthese_width[3,*,*]

err_m_h=stat_synthese_height[3,*,*]-stat_synthese_height[2,*,*]
err_p_h=stat_synthese_height[4,*,*]-stat_synthese_height[3,*,*]

err_2m_h=stat_synthese_height[3,*,*]-stat_synthese_height[1,*,*]
err_2p_h=stat_synthese_height[5,*,*]-stat_synthese_height[3,*,*]

err_m_a=stat_synthese_amplitude[3,*,*]-stat_synthese_amplitude[2,*,*]
err_p_a=stat_synthese_amplitude[4,*,*]-stat_synthese_amplitude[3,*,*]

err_2m_a=stat_synthese_amplitude[3,*,*]-stat_synthese_amplitude[1,*,*]
err_2p_a=stat_synthese_amplitude[5,*,*]-stat_synthese_amplitude[3,*,*]


Lmax=n_elements(stat_synthese_freq[0,*,0])
Nmax=n_elements(stat_synthese_freq[0,0,*])

std_dev_synth_freq=dblarr(5,Lmax,Nmax)
std_dev_synth_height=dblarr(5,Lmax,Nmax)
std_dev_synth_width=dblarr(5,Lmax,Nmax)
std_dev_synth_amplitude=dblarr(5,Lmax,Nmax)

;for l=0, Lmax-1 do for n=0,Nmax-1 do begin
	std_dev_synth_freq[0,*,*]=err_2m_f
	std_dev_synth_freq[1,*,*]=err_m_f
	std_dev_synth_freq[2,*,*]=stat_synthese_freq[3,*,*]
	std_dev_synth_freq[3,*,*]=err_p_f
	std_dev_synth_freq[4,*,*]=err_2p_f

	std_dev_synth_height[0,*,*]=err_2m_h
	std_dev_synth_height[1,*,*]=err_m_h
	std_dev_synth_height[2,*,*]=stat_synthese_height[3,*,*]
	std_dev_synth_height[3,*,*]=err_p_h
	std_dev_synth_height[4,*,*]=err_2p_h

	std_dev_synth_width[0,*,*]=err_2m_l
	std_dev_synth_width[1,*,*]=err_m_l
	std_dev_synth_width[2,*,*]=stat_synthese_width[3,*,*]
	std_dev_synth_width[3,*,*]=err_p_l
	std_dev_synth_width[4,*,*]=err_2p_l

	std_dev_synth_amplitude[0,*,*]=err_2m_a
	std_dev_synth_amplitude[1,*,*]=err_m_a
	std_dev_synth_amplitude[2,*,*]=stat_synthese_amplitude[3,*,*]
	std_dev_synth_amplitude[3,*,*]=err_p_a
	std_dev_synth_amplitude[4,*,*]=err_2p_a

;endfor

; **** petites et grande separations ****
if file2 ne '' then begin
	restore, file2

	err_m_d01=d01[3,*,*]-d01[2,*,*]
	err_2m_d01=d01[3,*,*]-d01[1,*,*]
	err_p_d01=d01[4,*,*]-d01[3,*,*]
	err_2p_d01=d01[5,*,*]-d01[3,*,*]

	err_m_dd01=dd01[3,*,*]-dd01[2,*,*]
	err_2m_dd01=dd01[3,*,*]-dd01[1,*,*]
	err_p_dd01=dd01[4,*,*]-dd01[3,*,*]
	err_2p_dd01=dd01[5,*,*]-dd01[3,*,*]

	err_m_Dnu=Dnu[3,*,*]-Dnu[2,*,*]
	err_2m_Dnu=Dnu[3,*,*]-Dnu[1,*,*]
	err_p_Dnu=Dnu[4,*,*]-Dnu[3,*,*]
	err_2p_Dnu=Dnu[5,*,*]-Dnu[3,*,*]

	err_m_d02=d02[3,*]-d02[2,*]
	err_2m_d02=d02[3,*]-d02[1,*]
	err_p_d02=d02[4,*]-d02[3,*]
	err_2p_d02=d02[5,*]-d02[3,*]

	std_dev_synth_Dnu=dblarr(5,Lmax,Nmax-1)
	std_dev_synth_dd01=dblarr(5,2,Nmax-2)
	std_dev_synth_d01=dblarr(5,2,Nmax-1)
	std_dev_synth_d02=dblarr(5,Nmax)

	std_dev_synth_Dnu[0,*,*]=err_2m_Dnu
	std_dev_synth_Dnu[1,*,*]=err_m_Dnu
	std_dev_synth_Dnu[2,*,*]=Dnu[3,*,*]
	std_dev_synth_Dnu[3,*,*]=err_p_Dnu
	std_dev_synth_Dnu[4,*,*]=err_2p_Dnu

	std_dev_synth_d01[0,*,*]=err_2m_d01
	std_dev_synth_d01[1,*,*]=err_m_d01
	std_dev_synth_d01[2,*,*]=d01[3,*,*]
	std_dev_synth_d01[3,*,*]=err_p_d01
	std_dev_synth_d01[4,*,*]=err_2p_d01

	std_dev_synth_dd01[0,*,*]=err_2m_dd01
	std_dev_synth_dd01[1,*,*]=err_m_dd01
	std_dev_synth_dd01[2,*,*]=dd01[3,*,*]
	std_dev_synth_dd01[3,*,*]=err_p_dd01
	std_dev_synth_dd01[4,*,*]=err_2p_dd01

	std_dev_synth_d02[0,*]=err_2m_d02
	std_dev_synth_d02[1,*]=err_m_d02
	std_dev_synth_d02[2,*]=d02[3,*]
	std_dev_synth_d02[3,*]=err_p_d02
	std_dev_synth_d02[4,*]=err_2p_d02
endif
;**************************************************
; ******* Phase d'ecriture sur fichier texte *****
;**************************************************

save,std_dev_synth_freq,std_dev_synth_height,std_dev_synth_width,std_dev_synth_amplitude,$
	std_dev_synth_Dnu,std_dev_synth_dd01,std_dev_synth_d01,std_dev_synth_d02,$
		filename=file_out+'.sav'

;****** Ecriture des variables frequences uniquement ********
openw,5,file_out+'_freq.txt'
printf,5,format='(A38)','** parametres individuel des modes **'

printf,5,format='(A38)','** frequences **'
legend=['l','-2sigma','-1sigma','median','+1sigma','+2sigma']
printf,5,format='(A3,5A14)',legend
for l=0,Lmax-1 do for n=0,Nmax-1 do begin
	printf,5,format='(I3,5F14.6)',l,std_dev_synth_freq[0,l,n],std_dev_synth_freq[1,l,n],$
		std_dev_synth_freq[2,l,n],std_dev_synth_freq[3,l,n],std_dev_synth_freq[4,l,n]
endfor

if file2 ne '' then begin

printf,5,format='(A38)','** Grande separation **'
legend=['l','freq','-2sigma','-1sigma','median','+1sigma','+2sigma']
printf,5,format='(A3,6A14)',legend
for l=0,Lmax-1 do for n=0,Nmax-2 do begin
	printf,5,format='(I3,6F14.6)',l,std_dev_synth_freq[2,l,n],std_dev_synth_Dnu[0,l,n],std_dev_synth_Dnu[1,l,n],$
		std_dev_synth_Dnu[2,l,n],std_dev_synth_Dnu[3,l,n],std_dev_synth_Dnu[4,l,n]
endfor

printf,5,format='(A38)','** separation d01 **'
legend=['meth','freq','-2sigma','-1sigma','median','+1sigma','+2sigma']
printf,5,format='(A6,6A14)',legend
for l=0,1 do begin
	for n=0,Nmax-3 do begin ; attention, ici l est un indice
	printf,5,format='(I6,6F14.6)',l,std_dev_synth_freq[2,l,n],std_dev_synth_d01[0,l,n],std_dev_synth_d01[1,l,n],$
		std_dev_synth_d01[2,l,n],std_dev_synth_d01[3,l,n],std_dev_synth_d01[4,l,n]
endfor
endfor

;printf,5,format='(A38)','** separation dd01 (Roxburgh) **'
;legend=['meth','freq','-2sigma','-1sigma','median','+1sigma','+2sigma']
;printf,5,format='(A6,6A14)',legend
;for l=0,1 do begin
;	if l eq 0 then printf,5,format='(A65)','formule : (nu(n-1,0)-4nu(n-1,1)+6nu(n,0)-4nu(n,1)+nu(n+1,0))/8'
;	if l eq 1 then printf,5,format='(A65)','formule : (nu(n-1,1)-4nu(n,0)+6nu(n,1)-4nu(n+1,0)+nu(n+1,1))/8'
;	for n=0,Nmax-3 do begin ; attention, ici l est un indice
;	printf,5,format='(I6,6F14.6)',l,std_dev_synth_freq[2,l,n+1],std_dev_synth_dd01[0,l,n],std_dev_synth_dd01[1,l,n],$
;		std_dev_synth_dd01[2,l,n],std_dev_synth_dd01[3,l,n],std_dev_synth_dd01[4,l,n]
;endfor
;endfor

printf,5,format='(A38)','** separation d02 **'
legend=['','freq','-2sigma','-1sigma','median','+1sigma','+2sigma']
printf,5,format='(A3,6A14)',legend
for n=0,Nmax-1 do begin
	printf,5,format='(I3,6F14.6)',l,std_dev_synth_freq[2,0,n],std_dev_synth_d02[0,n],std_dev_synth_d02[1,n],$
		std_dev_synth_d02[2,n],std_dev_synth_d02[3,n],std_dev_synth_d02[4,n]
endfor
endif

close,5
;*********************************************************************************************************

;********** Ecriture en colonne de tous les parametres (sauf les separations) ***********

openw,2,file_out+'_all_param.txt'

legend=['frequency','error','height','+error','-error','linewidth','+error','-error', 'amplitude','+error','-error']
printf, 2, '# Synthese of the parameters in a simple format.' 
printf, 2, '# CAUTION: For Amplitudes, please check whether your spectrum is double or single sided!'
printf, 2, '# Here, it is assumed single side (Kepler data convention). If it is not the case, you will need to multiply by two all amplitudes (and heights)!'
printf,2,format='(11A14)',legend
for l=0,Lmax-1 do begin
	printf,2,format='(A5)','l='+strtrim(l,1)

	for n=0,Nmax-1 do begin
		freq_err_mean=(std_dev_synth_freq[1,l,n]+std_dev_synth_freq[3,l,n])/2

		printf,2,format='(11F14.6)',std_dev_synth_freq[2,l,n],freq_err_mean,$
			std_dev_synth_height[2,l,n],std_dev_synth_height[3,l,n],std_dev_synth_height[1,l,n],$
			std_dev_synth_width[2,l,n],std_dev_synth_width[3,l,n],std_dev_synth_width[1,l,n], $
			std_dev_synth_amplitude[2,l,n],std_dev_synth_amplitude[3,l,n],std_dev_synth_amplitude[1,l,n]
	endfor

endfor


close,2
end
