function MS_Global_freq,file, ind0_nu, lmax, Nmax, Nf_ls, N_iter,Hist_length,output_ps

	tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

;*******************************************************************
;*************** PRELIMINAIRES AU CALCUL DE HAUTEUR ***************
;*******************************************************************
Nmax_tab=max(Nf_ls) ; Compatibility with mixed modes models

width=dblarr(lmax+1,Nmax_tab,N_iter)-1 ; initialize the table
freq=dblarr(lmax+1,Nmax_tab,N_iter)-1

Stat_Synthese_freq=dblarr(8,lmax+1,Nmax_tab)-1

; ---------- Determining the number of samples ------
skip_all=0

if lmax eq 1 then imax=ind0_nu+Nf_ls[0]+Nf_ls[1]-1
if lmax eq 2 then imax=ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2]-1
if lmax eq 3 then imax=ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2]+Nf_ls[3]-1

i=ind0_nu
N=0
print, 'Determning the number of samples...'
while (N lt 2) AND i lt imax do begin ; Finding the number of parameters using first free parameter
		file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
		restore,file_freq
		N=n_elements(param)
		i=i+1.
	endwhile	
	if i ge imax then begin
		print, 'No free frequencies found!' 
		N=10000. ; put a dummy value that will still allow the program to return stat_synthese_freq
	endif
; ----------------------------------------------------
	
;*********************************************************
;**************** CALCULS ET HISTOGRAMMES ****************
;*********************************************************
; classify all the samples for the frequencies... we need them to interpolate

en=0
for i=ind0_nu, ind0_nu+Nf_ls[0]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=0
	if n_elements(param) eq N then hist=build_histogram(param, Hist_length)
	if n_elements(param) eq 1 then  hist=build_histogram(replicate(param, N), Hist_length)
	if n_elements(param) ne 1 AND n_elements(param) ne N then begin
		print, 'Unexpected problem in the size of the param variables. Cannot proceed.'
		print, 'Debug Required. The program will stop now'
		stop
	endif
	distrib_stat=estimate_1_sigma_error(hist[0,*],hist[1,*],95.5,2,tab_critere) ; 
	Stat_synthese_freq[0:6,el,en]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
	Stat_synthese_freq[7,el,en]=distrib_stat[0] ; le max
	file_out=output_ps+'PDF_Freq_n' + addzeros(en+1) +  '_l' + strtrim(el, 2)  +'.eps'
	nice_hist1D, hist, ps=1, file_out=file_out, xr=xr, $
		;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
		title='!6n= n0+'+strtrim(en,1)+' l='+strtrim(el,1), xtitle='!6Frequency (microHz)', overplot=overplot, col=col, $
		fclose=fclose, show_stats=1, $
		legendsize=legendsize, legend_precision=legend_precision
	wait,0.01
	en=en+1
endfor
en=0
for i=ind0_nu+Nf_ls[0], ind0_nu+Nf_ls[0]+Nf_ls[1]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=1
	if n_elements(param) eq N then hist=build_histogram(param, Hist_length)
	if n_elements(param) eq 1 then  hist=build_histogram(replicate(param, N), Hist_length)
	if n_elements(param) ne 1 AND n_elements(param) ne N then begin
		print, 'Unexpected problem in the size of the param variables. Cannot proceed.'
		print, 'Debug Required. The program will stop now'
		stop
	endif
	distrib_stat=estimate_1_sigma_error(hist[0,*],hist[1,*],95.5,2,tab_critere) ; 
	Stat_synthese_freq[0:6,el,en]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
	Stat_synthese_freq[7,el,en]=distrib_stat[0] ; le max
	file_out=output_ps+'PDF_Freq_n' + addzeros(en+1) +  '_l' + strtrim(el, 2)  +'.eps'
	nice_hist1D, hist, ps=1, file_out=file_out, xr=xr, $
		;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
		title='!6n= n0+'+strtrim(en,1)+' l='+strtrim(el,1), xtitle='!6Frequency (microHz)', overplot=overplot, col=col, $
		fclose=fclose, show_stats=1, $
		legendsize=legendsize, legend_precision=legend_precision
	wait,0.01
	en=en+1
endfor
en=0
for i=Nmax+lmax+Nf_ls[0]+Nf_ls[1], ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=2
	if n_elements(param) eq N then hist=build_histogram(param, Hist_length)
	if n_elements(param) eq 1 then  hist=build_histogram(replicate(param, N), Hist_length)
	if n_elements(param) ne 1 AND n_elements(param) ne N then begin
		print, 'Unexpected problem in the size of the param variables. Cannot proceed.'
		print, 'Debug Required. The program will stop now'
		stop
	endif
	distrib_stat=estimate_1_sigma_error(hist[0,*],hist[1,*],95.5,2,tab_critere) ; 
	Stat_synthese_freq[0:6,el,en]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
	Stat_synthese_freq[7,el,en]=distrib_stat[0] ; le max
	file_out=output_ps+'PDF_Freq_n' + addzeros(en+1) + '_l' + strtrim(el, 2)  +'.eps'
	nice_hist1D, hist, ps=1, file_out=file_out, xr=xr, $
		;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
		title='!6n= n0+'+strtrim(en,1)+' l='+strtrim(el,1), xtitle='!6Frequency (microHz)', overplot=overplot, col=col, $
		fclose=fclose, show_stats=1, $
		legendsize=legendsize, legend_precision=legend_precision
	wait,0.01
	en=en+1
endfor
en=0
for i=ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2], ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2]+Nf_ls[3]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=3
	if n_elements(param) eq N then hist=build_histogram(param, Hist_length)
	if n_elements(param) eq 1 then  hist=build_histogram(replicate(param, N), Hist_length)
	if n_elements(param) ne 1 AND n_elements(param) ne N then begin
		print, 'Unexpected problem in the size of the param variables. Cannot proceed.'
		print, 'Debug Required. The program will stop now'
		stop
	endif
	distrib_stat=estimate_1_sigma_error(hist[0,*],hist[1,*],95.5,2,tab_critere) ; 
	Stat_synthese_freq[0:6,el,en]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
	Stat_synthese_freq[7,el,en]=distrib_stat[0] ; le max
	file_out=output_ps+'PDF_Freq_n' + addzeros(en+1) +  '_l' + strtrim(el, 2)  +'.eps'
	nice_hist1D, hist, ps=1, file_out=file_out, xr=xr, $
		;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
		title='!6n= n0+'+strtrim(en,1)+' l='+strtrim(el,1), xtitle='!6Frequency (microHz)', overplot=overplot, col=col, $
		fclose=fclose, show_stats=1, $
		legendsize=legendsize, legend_precision=legend_precision
	wait,0.01
	en=en+1
endfor


return, stat_synthese_freq
end
