; differs from oridinal version by the fact that it does not include 'param_' in the name file
; It also incorporate the new organization of parameters_length (with Nf_ls of fix size = 4)
; File number begins at 0 now
;.comp /home/obenomar/Dropbox/IDL/IDL_library/astro/pro/legend.pro
;.comp /home/obenomar/Difrot-Data/CPP_PostMCMC/build_histogram_v2.pro 
;.comp /home/obenomar/Difrot-Data/CPP_PostMCMC/MS_Global_rotinc_correlations.pro 
@MS_Global_rotinc_correlations
pro  local_rotinc_correlations, parameters_length, dir_files_samples, dir_out, identifier, idl_format

G=6.667d-8
Teff_sun= 5777d ; same values as in the function: seismic_vsini
Dnu_sun=135.1d
numax_sun=3150d
;numax_sun=3090d
R_sun=6.96342d5 ; in km
M_sun=1.98855d30 ; in kg
rho_sun=M_sun*1d3/(4d*!pi*(R_sun*1d5)^3d/3d) ; in g.cm-3

Nmax=parameters_length[0]
lmax=parameters_length[1] ; number of visibilities
if n_elements(idl_format) eq 0 then idl_format=0
if idl_format eq 0 then begin
	Nf=total(parameters_length[2:5]) ; The total number of parameter is total(Nf_ls)
	k=parameters_length[6]
	l=parameters_length[7]
	mm=parameters_length[8]
endif else begin
	Nf=parameters_length[2] ; The total number of parameter is total(Nf_ls)
	k=parameters_length[3]
	l=parameters_length[4]
	mm=parameters_length[5]
endelse
;; --------- MODIF AD HOC TO REMOVE WHAT WE THINK IS A SPURIOUS SOLUTION ------
;i_a3=Nmax+Nf+lmax+2
;restore, dir_files_samples + format_filename(i_a3) +'.sav'
;pos_OK=where(param gt -0.1 AND param le 0.1)
;a3p=param[pos_OK]
;; -------------------------------------------------

; Finding the first array of samples that contains the correct number of 
; samples. This because some parameters can be fixed and will lead to pos_OK=-1
; if we don't get the correct file with relevant (non-fixed) samples.
i=0
while n_elements(pos_OK) lt 2 do begin
	restore, dir_files_samples  + format_filename(i, idl_format) +'.sav'
	pos_OK=where(param ge 0)
	i=i+1
endwhile

i_a1=Nmax+Nf+lmax
restore, dir_files_samples  + format_filename(i_a1, idl_format) +'.sav'
;if n_elements(pos_OK) eq 0 then pos_OK=where(param ge 0) ; Basically all values by default

a1p=param[pos_OK]

i_a2=Nmax+Nf+lmax+1
restore, dir_files_samples  + format_filename(i_a2, idl_format) +'.sav'
a2p=param[pos_OK]

i_a3=Nmax+Nf+lmax+2
restore, dir_files_samples + format_filename(i_a3, idl_format) +'.sav'
a3p=param

i_a1cosi=Nmax+Nf+lmax+3
restore, dir_files_samples + format_filename(i_a1cosi, idl_format) +'.sav'
a1cosi=param[pos_OK]

i_a1sini=Nmax+Nf+lmax+4
restore, dir_files_samples + format_filename(i_a1sini, idl_format) +'.sav'
a1sini=param[pos_OK]

i_asym=Nmax+Nf+lmax+5
restore, dir_files_samples + format_filename(i_asym, idl_format) +'.sav'
asymp=param[pos_OK]

i_inc=Nmax+Nf+k+l+mm+lmax
restore, dir_files_samples  + format_filename(i_inc, idl_format) +'.sav'
incp=param[pos_OK]

if n_elements(a1p) gt 1 then begin
	if variance(a1p) eq 0 then begin
		;i_a1cosi=Nmax+Nf+lmax+6
		;restore, dir_files_samples + format_filename(i_a1cosi, idl_format) +'.sav'
		;a1cosi=param[pos_OK]
		;i_a1sini=Nmax+Nf+lmax+7
		;restore, dir_files_samples + format_filename(i_a1sini, idl_format) +'.sav'
		;a1sini=param[pos_OK]
		if variance(a1cosi) eq 0 then begin
			print, 'Warning: Could not find any sqrt(a1).cosi or sqrt(a1).sini values or even a1 values'
		endif else begin
			a1p=a1cosi^2 + a1sini^2
			incp=atan(a1sini/a1cosi)
			incp=incp*180./!pi
			Nb_classes=n_elements(a1p)/150.
			hist_a1p=build_histogram(a1p,Nb_classes, normalize=1)
			tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma
			stats_a1p=estimate_1_sigma_error(hist_a1p[0,*],hist_a1p[1,*],68.3,2,tab_critere)
		
			hist_incp=build_histogram(incp,Nb_classes, normalize=1)
			tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma
			stats_incp=estimate_1_sigma_error(hist_incp[0,*],hist_incp[1,*],68.3,2,tab_critere)
			print, '--------------------------'
			print, 'ID = ' + identifier 
			print, 'index_a1cosi =' + strtrim(i_a1cosi,2) + '  /  index_a1sini =' + strtrim(i_a1sini,2)
			print, '--------------------------'
			;stop
			file_split_inc=dir_out + identifier + '_rotinc.txt'
			openw, 3, file_split_inc
				printf, 3, '# Made by local_rotinc_correlations from sqrt(a1).cos(i) and sqrt(a1).sin(i)'
				printf, 3, '# Short summary files with the splitting values and the stellar inclination'
				printf, 3, '# index_a1cosi =' + strtrim(i_a1cosi,2) + '  /  index_a1sini =' + strtrim(i_a1sini,2)
				printf, 3, '# Format of the outputs'
				printf, 3, '# a1   0%     2.25%   16%     50%    84%    97.75%   100%'
				str=''
				for j=0, n_elements(stats_a1p)-1 do begin
					str=str + string(stats_a1p[j], format='(f20.9)')
				endfor
				printf, 3, str
				printf, 3, '# inc   0%     2.25%   16%     50%    84%    97.75%   100%'			
				str=''
				for j=0, n_elements(stats_incp)-1 do begin
					str=str + string(stats_incp[j], format='(f20.9)')
				endfor
				printf, 3, str
			close,3
		endelse
	endif
endif
;stop

if n_elements(asymp) gt 1 then begin 
	if variance(asymp) eq 0 then asymp=0
endif else begin
	asymp=0
endelse
if n_elements(a2p) gt 1 then begin 
	if variance(a2p) lt 1d-10 or finite(variance(a2p)) eq 0 then a2p=0
endif else begin
	a2p=0
endelse
if n_elements(a3p) gt 1 then begin
	if variance(a3p) eq 0 then a3p=0
endif else begin
	a3p=0
endelse
if n_elements(mag_bp) gt 1 then begin 
	if variance(mag_bp) eq 0 then mag_bp=0
endif else begin
	mag_bp=0
endelse
if n_elements(mag_alfap) gt 1 then begin 
	if variance(mag_alfap) eq 0 then mag_alfap=0
endif else begin
	mag_alfap=0
endelse

starID=identifier
dataID=''
if identifier eq 'kplr006116048_45_COR_PSD_filt_inp' then begin
	starID='HD 181207'
	dataID='Inpainted'
endif
if identifier eq 'kplr006116048_kasoc-psd_slc_v1' then begin
	starID='HD 181207'
	dataID='KASOC'
endif
if identifier eq 'kplr008379927_91_COR_PSD_filt_inp' then begin
	starID='HD 187160'
	dataID='Inpainted'
endif
if identifier eq 'kplr008379927_kasoc-psd_slc_v2' then begin
	starID='HD 187160'
	dataID='KASOC'
endif
if identifier eq '12069424' then begin
	starID='16 Cyg A'
	dataID=''
endif
if identifier eq '12069449' then begin
	starID='16 Cyg B'
	dataID=''
endif
if identifier eq '' then stop

if n_elements(a3p) ne 1 then begin
	Nb_classes=n_elements(a3p)/150.
	hist_a3=build_histogram(a3p,Nb_classes, normalize=1)
    cdf_a3=dblarr(2, n_elements(hist_a3[1,*]))
    for i=0, n_elements(hist_a3[0,*])-1 do cdf_a3[1,i]=total(hist_a3[1,0:i])/total(hist_a3[1,*])
	cdf_a3[0,*]=hist_a3[0,*]

        ; --- calculate the probability to be below 0 ----
        x0=0
        if x0 ge min(hist_a3[0,*]) AND x0 le max(hist_a3[0,*]) then begin ; case where p(x<0) > 0 AND p(x>0) >0
                proba_a3_neg=interpol(cdf_a3[1,*], hist_a3[0,*], 0., /quadratic) * 100.
                proba_a3_pos= 100. - proba_a3_neg
        endif
        if min(hist_a3[0,*]) ge x0 then begin ; case where p(x<0) = 0% and p(x>0) = 100%
                proba_a3_neg=0.
                proba_a3_pos=100.
        endif
        if max(hist_a3[0,*]) le x0 then begin ; case where p(x<0) = 100% and p(x>0) = 0%
                proba_a3_neg=100.
                proba_a3_pos=0.
        endif
	openw, 3, dir_out + 'a3_proba.txt'
		printf, 3, 'Probability for a3>0 :', proba_a3_pos
		printf, 3, 'Probability for a3<0 :', proba_a3_neg
		printf, 3, '#cdf'
	        for i=0, n_elements(hist_a3[0,*])-1 do printf, 3, cdf_a3[0,i], cdf_a3[1,i]
  
	close,3 
endif

if n_elements(a2p) ne 1 then begin
	print, 'WARNING STOP: Cannot deal with a2p in the current version of the local fit'
	print, '              Please update the code before dealing with this'
	print, '              The program will stop now'
	stop
	; --------------------- TO BE UPDATED ---------------------
	; ---- Legacy code from MS_Global_rotinc_correlations -----
	; ---------------------------------------------------------
	;beta0=a2p
	;p=linfit(findgen(n_elements(stat_synthese_freq[3,0,*])), stat_synthese_freq[3,0,*], /double)
	;Dnu=p[1]
	;eta0=4.*!pi*(Dnu_sun/Dnu)^2 / (3. * G)
	;eta0a1=eta0*(a1p*1d-6)^2
	;
	;Nb_classes=n_elements(beta0)/150.
	;hist_beta0=build_histogram(beta0,Nb_classes, normalize=1)
	;hist_eta0a1=build_histogram(eta0a1,Nb_classes, normalize=1)
	;hist_eta0a1[1,*]=hist_eta0a1[1,*]*max(hist_beta0[1,*])/max(hist_eta0a1[1,*])

	;hist_DR1=build_histogram(beta0*3./(8. * !pi),Nb_classes, normalize=1)
	;hist_DR1[0,*]=hist_DR1[0,*]*100. ; convert into %
	;hist_DR0=build_histogram(eta0a1*3./(8. * !pi),Nb_classes, normalize=1)
	;hist_DR0[0,*]=hist_DR0[0,*]*100. ; convert into %
	;hist_DR0[1,*]=hist_DR0[1,*]*max(hist_DR1[1,*])/max(hist_DR0[1,*])
    ;
	;cdf_DR1=dblarr(n_elements(hist_DR1[0,*]))
	;for i=0, n_elements(hist_DR1[0,*])-1 do cdf_DR1[i]=total(hist_DR1[1,0:i])/total(hist_DR1[1,*])
	;cdf_DR0=dblarr(n_elements(hist_DR0[0,*]))
	;for i=0, n_elements(hist_DR0[0,*])-1 do cdf_DR0[i]=total(hist_DR0[1,0:i])/total(hist_DR0[1,*])

	; --- calculate the probability to be below 0 ----
	;x0=0
	;if x0 ge min(hist_DR1[0,*]) AND x0 le max(hist_DR1[0,*]) then begin ; case where p(x<0) > 0 AND p(x>0) >0
	;	proba_DR1_prolate=interpol(cdf_DR1, hist_DR1[0,*], 0., /quadratic) * 100.
	;	proba_DR1_oblate= 100. - proba_DR1_prolate
	;endif 
	;if min(hist_DR1[0,*]) ge x0 then begin ; case where p(x<0) = 0% and p(x>0) = 100%
	;	proba_DR1_prolate=0.
	;	proba_DR1_oblate=100.
	;endif
	;if max(hist_DR1[0,*]) le x0 then begin ; case where p(x<0) = 100% and p(x>0) = 0%
	;	proba_DR1_prolate=100.
	;	proba_DR1_oblate=0.
	;endif

   	;file_out=dir_out + 'PDF_beta0_vs_centforce.eps'
	;nimp,name=file_out,/paper,/eps
;	;out=nice_hist1D_compare(hist_beta0, hist_eta0a1, xr=[min([hist_beta0[0,*], hist_eta0a1[0,*]]), max([hist_beta0[0,*], hist_eta0a1[0,*]])], $
;	;	y_max=max( [max(hist_beta0[1,*]), max(hist_eta0a1[1,*])] ), $
;	;	title=title, xtitle=textoidl('\beta_0') + ' vs. ' + textoidl('\eta_0 a^2_1') + ' (no unit)', show_stats=1, $
;	;	legendsize=legendsize, legend_precision=legend_precision, $
;	;	col_hists=col_hists, no_number=1)
	;out=nice_hist1D_compare(hist_DR1, hist_DR0, xr=[min([hist_DR0[0,*], hist_DR1[0,*]]), max([hist_DR0[0,*], hist_DR1[0,*]])], y_max=1.1*max( [max(hist_DR0[1,*]), max(hist_DR1[1,*])] ), $
	;	title=title, xtitle=textoidl('\Delta R / R') + ' (%)', show_stats=1, $
	;	legendsize=legendsize, legend_precision=legend_precision, col_hists=col_hists, no_number=1)
	;labels=[textoidl('\beta_0'), textoidl('a^2_1 \eta_0')]
	;label_code=[0,0]
	;labels_colors=fsc_color(['Black', 'Red'])
	;;stop
	;legend, labels, psym=label_code, colors=labels_colors, textcolors=labels_colors, charsize=1.4, /right, box=0
	;deltax=(max(hist_DR1[0,*]) - min(hist_DR1[0,*]))
	;deltay=1.1*max( [max(hist_DR1[1,*]), max(hist_DR0[1,*])] )
	;xyouts, min(hist_DR1[0,*]) + deltax*0.05, 0.930*deltay, starID, charsize=1.5, color=fsc_color('Black')
	;xyouts, min(hist_DR1[0,*]) + deltax*0.05, 0.895*deltay, dataID, charsize=1.4, color=fsc_color('Black')
	;xyouts, min(hist_DR1[0,*]) + deltax*0.05, 0.840*deltay, textoidl('P_{oblate} = ') + strtrim(string(proba_DR1_oblate, format='(f5.1)'),2) + '%', charsize=1.4, color=fsc_color('Black')
	;fimp

   	;print, 'Prolate Probability: ', proba_DR1_prolate
	;print, 'Oblate Probability: ', proba_DR1_oblate
    ;    openw, 3, dir_out + 'asphericity_proba.txt'
    ;            printf, 3, 'Prolate Probability: ', proba_DR1_prolate
    ;    printf ,3, 'Oblate Probability: ', proba_DR1_oblate
    ;    close, 3	
	;save, eta0a1, beta0, hist_beta0, hist_eta0a1, hist_DR0, hist_DR1, cdf_DR1, cdf_DR0, filename=dir_out + 'beta0_etaa1_samples_pdfs_cdfs.sav'
	;;stop
endif else begin
	print, 'a2 fucntionnality disabled'
	;posn0=where(stat_synthese_freq[3,0,*] gt 0)
	;p=linfit(findgen(n_elements(stat_synthese_freq[3,0,posn0])), stat_synthese_freq[3,0,posn0], /double)
	;Dnu=p[1]
	;eta0=4.*!pi*(Dnu_sun/Dnu)^2 / (3. * G)
	;eta0a1=eta0*(a1p*1d-6)^2
	;beta0=eta0a1
endelse

name_out=identifier + '_'
;stop
if (stddev(a1p) ne 0 AND stddev(incp) ne 0) AND n_elements(a1p) gt 1 AND n_elements(incp) gt 1 then begin 
	;show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=1d5 *beta0*3./(8. * !pi), a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=0, extra=1d5 *eta0a1*3./(8. * !pi), dir_out, name_out+'rotonly_'  ; all except asymetry
	;show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=1d5 *beta0*3./(8. * !pi), a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=asymp, extra=1d5 *eta0a1*3./(8. * !pi), dir_out, name_out ; all
	show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=a2p, a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=0, extra=0, dir_out, name_out+'rotonly_'
endif else print, 'No correlation diagram to show: a1 and/or inclination are fixed'

end
