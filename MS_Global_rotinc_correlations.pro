; differs from oridinal version by the fact that it does not include 'param_' in the name file
; It also incorporate the new organization of parameters_length (with Nf_ls of fix size = 4)
; File number begins at 0 now
;.comp /home/obenomar/Dropbox/IDL/IDL_library/astro/pro/legend.pro
;.comp /home/obenomar/Difrot-Data/CPP_PostMCMC/build_histogram_v2.pro 
;.comp /home/obenomar/Difrot-Data/CPP_PostMCMC/MS_Global_rotinc_correlations.pro 
pro  MS_Global_rotinc_correlations, file_synthese, dir_files_samples, dir_out, identifier, idl_format

G=6.667d-8
Teff_sun= 5777d ; same values as in the function: seismic_vsini
Dnu_sun=135.1d
numax_sun=3150d
;numax_sun=3090d
R_sun=6.96342d5 ; in km
M_sun=1.98855d30 ; in kg
rho_sun=M_sun*1d3/(4d*!pi*(R_sun*1d5)^3d/3d) ; in g.cm-3

restore, file_synthese ; this will give us the index number for a1,a2,a3,inc
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

i_a1=Nmax+Nf+lmax
restore, dir_files_samples  + format_filename(i_a1, idl_format) +'.sav'
if n_elements(pos_OK) eq 0 then pos_OK=where(param ge 0) ; Basically all values by default
a1p=param[pos_OK]

i_a2=Nmax+Nf+lmax+1
restore, dir_files_samples  + format_filename(i_a2, idl_format) +'.sav'
a2p=param[pos_OK]

i_a3=Nmax+Nf+lmax+2
restore, dir_files_samples + format_filename(i_a3, idl_format) +'.sav'
a3p=param

i_mag_b=Nmax+Nf+lmax+3
restore, dir_files_samples + format_filename(i_mag_b, idl_format) +'.sav'
mag_bp=param[pos_OK]

i_mag_alfa=Nmax+Nf+lmax+4
restore, dir_files_samples + format_filename(i_mag_alfa, idl_format) +'.sav'
mag_alfap=param[pos_OK]

i_asym=Nmax+Nf+lmax+5
restore, dir_files_samples + format_filename(i_asym, idl_format) +'.sav'
asymp=param[pos_OK]

i_inc=Nmax+Nf+k+l+mm+lmax
restore, dir_files_samples  + format_filename(i_inc, idl_format) +'.sav'
incp=param[pos_OK]

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
	;stop
endif

if n_elements(a2p) ne 1 then begin
	beta0=a2p
	p=linfit(findgen(n_elements(stat_synthese_freq[3,0,*])), stat_synthese_freq[3,0,*], /double)
	Dnu=p[1]
	eta0=4.*!pi*(Dnu_sun/Dnu)^2 / (3. * G)
	eta0a1=eta0*(a1p*1d-6)^2
	
	Nb_classes=n_elements(beta0)/150.
	hist_beta0=build_histogram(beta0,Nb_classes, normalize=1)
	hist_eta0a1=build_histogram(eta0a1,Nb_classes, normalize=1)
	hist_eta0a1[1,*]=hist_eta0a1[1,*]*max(hist_beta0[1,*])/max(hist_eta0a1[1,*])

	hist_DR1=build_histogram(beta0*3./(8. * !pi),Nb_classes, normalize=1)
	hist_DR1[0,*]=hist_DR1[0,*]*100. ; convert into %
	hist_DR0=build_histogram(eta0a1*3./(8. * !pi),Nb_classes, normalize=1)
	hist_DR0[0,*]=hist_DR0[0,*]*100. ; convert into %
	hist_DR0[1,*]=hist_DR0[1,*]*max(hist_DR1[1,*])/max(hist_DR0[1,*])

	cdf_DR1=dblarr(n_elements(hist_DR1[0,*]))
	for i=0, n_elements(hist_DR1[0,*])-1 do cdf_DR1[i]=total(hist_DR1[1,0:i])/total(hist_DR1[1,*])
	cdf_DR0=dblarr(n_elements(hist_DR0[0,*]))
	for i=0, n_elements(hist_DR0[0,*])-1 do cdf_DR0[i]=total(hist_DR0[1,0:i])/total(hist_DR0[1,*])

	; --- calculate the probability to be below 0 ----
	x0=0
	if x0 ge min(hist_DR1[0,*]) AND x0 le max(hist_DR1[0,*]) then begin ; case where p(x<0) > 0 AND p(x>0) >0
		proba_DR1_prolate=interpol(cdf_DR1, hist_DR1[0,*], 0., /quadratic) * 100.
		proba_DR1_oblate= 100. - proba_DR1_prolate
	endif 
	if min(hist_DR1[0,*]) ge x0 then begin ; case where p(x<0) = 0% and p(x>0) = 100%
		proba_DR1_prolate=0.
		proba_DR1_oblate=100.
	endif
	if max(hist_DR1[0,*]) le x0 then begin ; case where p(x<0) = 100% and p(x>0) = 0%
		proba_DR1_prolate=100.
		proba_DR1_oblate=0.
	endif

	file_out=dir_out + 'PDF_beta0_vs_centforce.eps'
	nimp,name=file_out,/paper,/eps
;	out=nice_hist1D_compare(hist_beta0, hist_eta0a1, xr=[min([hist_beta0[0,*], hist_eta0a1[0,*]]), max([hist_beta0[0,*], hist_eta0a1[0,*]])], $
;		y_max=max( [max(hist_beta0[1,*]), max(hist_eta0a1[1,*])] ), $
;		title=title, xtitle=textoidl('\beta_0') + ' vs. ' + textoidl('\eta_0 a^2_1') + ' (no unit)', show_stats=1, $
;		legendsize=legendsize, legend_precision=legend_precision, $
;		col_hists=col_hists, no_number=1)
	out=nice_hist1D_compare(hist_DR1, hist_DR0, xr=[min([hist_DR0[0,*], hist_DR1[0,*]]), max([hist_DR0[0,*], hist_DR1[0,*]])], y_max=1.1*max( [max(hist_DR0[1,*]), max(hist_DR1[1,*])] ), $
		title=title, xtitle=textoidl('\Delta R / R') + ' (%)', show_stats=1, $
		legendsize=legendsize, legend_precision=legend_precision, col_hists=col_hists, no_number=1)
	labels=[textoidl('\beta_0'), textoidl('a^2_1 \eta_0')]
	label_code=[0,0]
	labels_colors=fsc_color(['Black', 'Red'])
	;stop
	legend, labels, psym=label_code, colors=labels_colors, textcolors=labels_colors, charsize=1.4, /right, box=0
	deltax=(max(hist_DR1[0,*]) - min(hist_DR1[0,*]))
	deltay=1.1*max( [max(hist_DR1[1,*]), max(hist_DR0[1,*])] )
	xyouts, min(hist_DR1[0,*]) + deltax*0.05, 0.930*deltay, starID, charsize=1.5, color=fsc_color('Black')
	xyouts, min(hist_DR1[0,*]) + deltax*0.05, 0.895*deltay, dataID, charsize=1.4, color=fsc_color('Black')
	xyouts, min(hist_DR1[0,*]) + deltax*0.05, 0.840*deltay, textoidl('P_{oblate} = ') + strtrim(string(proba_DR1_oblate, format='(f5.1)'),2) + '%', charsize=1.4, color=fsc_color('Black')
	fimp

	print, 'Prolate Probability: ', proba_DR1_prolate
	print, 'Oblate Probability: ', proba_DR1_oblate
        openw, 3, dir_out + 'asphericity_proba.txt'
                printf, 3, 'Prolate Probability: ', proba_DR1_prolate
        printf ,3, 'Oblate Probability: ', proba_DR1_oblate
        close, 3	
save, eta0a1, beta0, hist_beta0, hist_eta0a1, hist_DR0, hist_DR1, cdf_DR1, cdf_DR0, filename=dir_out + 'beta0_etaa1_samples_pdfs_cdfs.sav'
	;stop
endif else begin
	posn0=where(stat_synthese_freq[3,0,*] gt 0)
	p=linfit(findgen(n_elements(stat_synthese_freq[3,0,posn0])), stat_synthese_freq[3,0,posn0], /double)
	Dnu=p[1]
	eta0=4.*!pi*(Dnu_sun/Dnu)^2 / (3. * G)
	eta0a1=eta0*(a1p*1d-6)^2
	beta0=eta0a1
endelse

name_out=identifier + '_'
;stop
if (stddev(a1p) ne 0 AND stddev(incp) ne 0) AND n_elements(a1p) gt 1 AND n_elements(incp) gt 1 then begin 
	;show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=a2p, a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=asymp, extra=extra, posextra=posextra, dir_out, name_out ; all
	;show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=a2p, a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=0, extra=extra, posextra=posextra, dir_out, name_out+'rotonly_' ; all except asymetry
	show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=1d5 *beta0*3./(8. * !pi), a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=0, extra=1d5 *eta0a1*3./(8. * !pi), dir_out, name_out+'rotonly_'  ; all except asymetry
	show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=1d5 *beta0*3./(8. * !pi), a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=asymp, extra=1d5 *eta0a1*3./(8. * !pi), dir_out, name_out ; all
;	show_hist_a1_a2_a3_mag_asym_inc_matrix, a1p, incp, a2=1d5 *beta0*3./(8. * !pi), a3=a3p, mag_b=mag_bp, mag_alfa=mag_alfap, asym=0, extra=1d5 *eta0a1*3./(8. * !pi), dir_out, name_out+'rotonly_' ; all except asymetry
endif else print, 'No correlation diagram to show: a1 and/or inclination are fixed'

end

function format_filename, index, idl_format
	if n_elements(idl_format) eq 0 then idl_format=0

	if idl_format eq 0 then pref='' else pref='param_'

	if index lt 10 then str=pref + '00' + strtrim(round(index),2)
	if index ge 10 AND index lt 100 then str=pref + '0' + strtrim(round(index),2)
	if index ge 100 then str=pref + strtrim(round(index),2)
	
return, str
end


function nice_hist1D_compare, hist1, hist2, xr=xr, y_max=y_max, $
	title=title, xtitle=xtitle, show_stats=show_stats, $
	legendsize=legendsize, legend_precision=legend_precision, $
	col_hists=col_hists, no_number=no_number

if n_elements(no_number) eq 0 then no_number=0 ; controls if we show the 1 sigma values
if n_elements(col_hists) eq 0 then col_hists=['Black', 'Red']
if n_elements(title) eq 0 then title=''
if n_elements(xtitle) eq 0 then xtitle=''
overplot=0
pmulti=0
;if n_elements(col) eq 0 then col=['Black', 'Red', 'Red']
fclose=0 ; a variable that forces to close the ps file after the plot... to use in case of overplot
if n_elements(show_stats) eq 0 then show_stats=1 ; defaut, we plot the stats
if n_elements(legendsize) eq 0 then size_char=0.85 else size_char=legendsize
if n_elements(legend_precision) eq 0 then txtformat='(f8.2)' else txtformat='(f8.' + strtrim(floor(legend_precision),2) +  ')'
smooth_coef=1

		xscale=(max(hist1[0,*]) - min(hist1[0,*]))*0.1
		if n_elements(y_max) eq 0 then y_max=max(hist1[1,*])*1.1
		if n_elements(xr) eq 0 then begin
			x_min=min(hist1[0,*]) - xscale
			x_max=max(hist1[0,*]) + xscale
		endif else begin
			x_min=xr[0]
			x_max=xr[1]
		endelse
		frac=0.02

		size_c=1.5
		
		plot, hist1[0,*], smooth(hist1[1,*], smooth_coef), psym=10,$
			xtitle=xtitle, ytitle=textoidl('\propto') + 'Probability density', title=title, charsize=size_c, $
			background=fsc_color('White'), color=fsc_color(col_hists[0]),yr=[0, y_max], xr=[x_min, x_max], /xst, /yst
		oplot, hist2[0,*], smooth(hist2[1,*], smooth_coef), psym=10, color=fsc_color(col_hists[1])

if show_stats eq 1 then begin
		tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

		stats1=estimate_1_sigma_error( hist1[0,*],hist1[1,*],68.3,2,tab_critere)

		b_max=interpol(smooth(hist1[1,*],smooth_coef), hist1[0,*], stats1[6])
		plots, [stats1[6],stats1[6]], [0, b_max] , color=fsc_color(col_hists[0]), thick=5. ; median
		if no_number eq 0 then xyouts, stats1[6]+(stats1[1]*0.01), b_max*(1-0.82), strtrim(string(stats1[6], format=txtformat),1), color=fsc_color(col_hists[0]),charsize=size_char

		b_max=interpol(smooth(hist1[1,*],smooth_coef), hist1[0,*], stats1[5])
		plots, [stats1[5],stats1[5]],[0, b_max],  color=fsc_color(col_hists[0]), thick=5.,linestyle=2 ; -1sigma
		if no_number eq 0 then xyouts, stats1[5]+(stats1[1]*0.03	), y_max*(1-0.94), strtrim(string(stats1[5], format=txtformat),1), color=fsc_color(col_hists[0]),charsize=size_char

		b_max=interpol(smooth(hist1[1,*],smooth_coef), hist1[0,*], stats1[7])
		plots,  [stats1[7],stats1[7]],[0, b_max], color=fsc_color(col_hists[0]), thick=5.,linestyle=2 ; +1sigma
		if no_number eq 0 then xyouts, stats1[7]+(stats1[1]*0.03), y_max*(1-0.94), strtrim(string(stats1[7], format=txtformat),1), color=fsc_color(col_hists[0]),charsize=size_char
		; ---------------------
		stats2=estimate_1_sigma_error( hist2[0,*],hist2[1,*],68.3,2,tab_critere)

		b_max=interpol(smooth(hist2[1,*],smooth_coef), hist2[0,*], stats2[6])
		plots, [stats2[6],stats2[6]], [0, b_max] , color=fsc_color(col_hists[1]), thick=5. ; median
		if no_number eq 0 then xyouts, stats2[6]+(stats2[1]*0.01), b_max*(1-0.82), strtrim(string(stats2[6], format=txtformat),1), color=fsc_color(col_hists[1]),charsize=size_char

		b_max=interpol(smooth(hist2[1,*],smooth_coef), hist2[0,*], stats2[5])
		plots, [stats2[5],stats2[5]],[0, b_max],  color=fsc_color(col_hists[1]), thick=5.,linestyle=2 ; -1sigma
		if no_number eq 0 then xyouts, stats2[5]+(stats2[1]*0.03	), y_max*(1-0.94), strtrim(string(stats2[5], format=txtformat),1), color=fsc_color(col_hists[1]),charsize=size_char

		b_max=interpol(smooth(hist2[1,*],smooth_coef), hist2[0,*], stats2[7])
		plots,  [stats2[7],stats2[7]],[0, b_max], color=fsc_color(col_hists[1]), thick=5.,linestyle=2 ; +1sigma
		if no_number eq 0 then xyouts, stats2[7]+(stats2[1]*0.03), y_max*(1-0.94), strtrim(string(stats2[7], format=txtformat),1), color=fsc_color(col_hists[1]),charsize=size_char
endif


	out={stats1:dblarr(n_elements(stats1)), stats2:dblarr(n_elements(stats2))}
	out.stats1=stats1
	out.stats2=stats2

	return, out
end
