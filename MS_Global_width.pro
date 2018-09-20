;fonction that calculate the width as in the model (ie, by interpolation)
; Recoit :
;		- file : noyau du nom de fichier contant les variables
;		- ind0_nu: First index pointing to frequencies
;       - ind0_W: First index pointing to Widths
;       - ind0_split: First index pointing to splittings. Used to compute the mode blending
;		- Nf_ls: Number of Frequency for each degree (array of size Max(lmax)=4)
;		- N_iter : nombre d'echantillons considÈrÈs
;		- Hist_length : nombre de classes pour l'histogramme
;		- output_ps : repertoire+noyau des fichers ps de sortie (si ps=1)
;		- Save_Samples : If you want to save all samples into IDL binary files (*.sav)
function MS_Global_width,file, ind0_nu, ind0_W, ind0_split, lmax, Nmax, Nf_ls, N_iter,Hist_length,output_ps, rules, rules_params, save_Samples=save_Samples

	tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

;*******************************************************************
;*************** PRELIMINAIRES AU CALCUL DE HAUTEUR ***************
;*******************************************************************
Nmax_tab=max(Nf_ls) ; Compatibility with mixed modes models

width=dblarr(lmax+1,Nmax_tab,N_iter) -1
freq=dblarr(lmax+1,Nmax_tab,N_iter) -1

hist_width=dblarr(lmax+1,Nmax_tab,N_iter/Hist_length) 
x_hist_width=dblarr(lmax+1,Nmax_tab,N_iter/Hist_length)
Stat_Synthese_width=dblarr(8,lmax+1,Nmax_tab) -1

;*********************************************************
;**************** CALCULS ET HISTOGRAMMES ****************
;*********************************************************
; classify all the samples for the frequencies... we need them to interpolate
en=0
eldone=dblarr(lmax+1) + 1000.
for i=ind0_nu, ind0_nu+Nf_ls[0]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=0
	freq[el,en, *]=param
	eldone[0]=1
	en=en+1
endfor
en=0
for i=ind0_nu+Nf_ls[0], ind0_nu+Nf_ls[0]+Nf_ls[1]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=1
	freq[el,en, *]=param
	eldone[1]=1
	en=en+1
endfor
en=0
for i=Nmax+lmax+Nf_ls[0]+Nf_ls[1], ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=2
	freq[el,en, *]=param
	eldone[2]=1
	en=en+1
endfor
en=0
for i=ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2], ind0_nu+Nf_ls[0]+Nf_ls[1]+Nf_ls[2]+Nf_ls[3]-1 do begin
	file_freq=file_syntax(i, file) ; detect automatically how many 0 I should put
	restore,file_freq
	el=3
	freq[el,en, *]=param
	eldone[3]=1
	en=en+1
endfor

for i=0,Nmax-1 do begin
	file_largeur=file_syntax(ind0_W+i, file) ; detect automatically how many 0 I should put
	restore,file_largeur
	width[0,i, *]=param
endfor

;for i=0,Nmax-1 do begin
;	for v_i=1, lmax do begin
;		for samp=long(0), N_iter-1 do $
;			width[v_i, i, samp]=abs(interpol(reform(width[0,*,samp]), reform(freq[0,*,samp]), reform(freq[v_i, i, samp])))
;	endfor	
;endfor
;for i=0, Nmax-1 do begin
;		for v_i=0, lmax do begin
;			hist_width[v_i,i,*]=HISTOGRAM(width[v_i,i, *],nbins=1.*N_iter/Hist_length)
;			x_hist_width[v_i,i,*]=findgen(n_elements(hist_width[v_i,i,*]))*(max(width[v_i,i, *])-min(width[v_i,i, *]))/n_elements(hist_width[v_i,i,*])+min(width[v_i,i, *])
;		endfor
;endfor

for v_i=1, lmax do begin
	if rules[v_i] eq 'p' then begin ; p modes with interpolated widths
		for i=0,Nf_ls[v_i]-1 do begin
			for samp=long(0), N_iter-1 do $
				width[v_i, i, samp]=abs(interpol(reform(width[0,*,samp]), reform(freq[0,*,samp]), reform(freq[v_i, i, samp])))
		endfor	
	endif 
	if rules[v_i] eq 'm' then begin ; mixed modes with independent widths
		for i=0,Nf_ls[v_i]-1 do begin
			file_largeur=file_syntax(rules_params[v_i,2]+i, file) ; rules_params for widths are on slot 2
			restore, file_largeur
			width[v_i, i, *]=param
			;print, 'debug stop in MS_Global_width'
			;stop
		endfor
	endif
	if rules[v_i] ne 'p' AND rules[v_i] ne 'm' then begin
		print, 'Unknown rule for this model. The program requires debuging'
		print, 'The program will stop now'
		stop
	endif
endfor

for v_i=0, lmax do begin
	for i=0, Nf_ls[v_i]-1 do begin
			hist_width[v_i,i,*]=HISTOGRAM(width[v_i,i, *],nbins=1.*N_iter/Hist_length)
			x_hist_width[v_i,i,*]=findgen(n_elements(hist_width[v_i,i,*]))*(max(width[v_i,i, *])-min(width[v_i,i, *]))/n_elements(hist_width[v_i,i,*])+min(width[v_i,i, *])
	endfor
endfor
hist_width=hist_width/N_iter

print, '    - Saving the processed widths for the fitted modes...'
if save_Samples eq 1 then save, width, hist_width, x_hist_width, filename=output_ps + 'Samples_Widths.sav'
if save_Samples ne 1 then save, hist_width, x_hist_width, filename=output_ps + 'Samples_Widths.sav'

print, '    - Computing/ploting/saving the ratio a1/Gamma sample-by-sample (modes blending ratio)...'
file_splitting=file_syntax(ind0_split, file)
restore, file_splitting
splitting_a1=param
loga1_ov_width=dblarr(lmax+1,Nmax_tab,N_iter)
hist_a1_ov_width=dblarr(lmax+1,Nmax_tab,N_iter/Hist_length)
x_hist_a1_ov_width=dblarr(lmax+1,Nmax_tab,N_iter/Hist_length)
Stat_Synthese_a1_ov_width=dblarr(8,lmax+1,Nmax_tab)-1
;for i=0,Nmax-1 do begin
;	for v_i=0, lmax do begin
;		log_a1w=alog(splitting_a1) - alog(width[v_i,i,*])
;		loga1_ov_width[v_i, i, *]=log_a1w
;		hist_a1_ov_width[v_i,i,*]=HISTOGRAM(log_a1w,nbins=1.*N_iter/Hist_length)
;		x_hist_a1_ov_width[v_i,i,*]=findgen(n_elements(hist_a1_ov_width[v_i,i,*]))*(max(log_a1w)-min(log_a1w))/n_elements(hist_a1_ov_width[v_i,i,*])+min(log_a1w)
;		hist0=dblarr(2, n_elements(x_hist_a1_ov_width[v_i,i,*]))
;		hist0[0,*]=x_hist_a1_ov_width[v_i,i,*]
;		hist0[1,*]=hist_a1_ov_width[v_i,i,*]
;		file_out=output_ps+'logPDF_a1_ov_Width_n'+addzeros(i) + '_l' + strtrim(v_i,2) +'.eps'
;		nice_hist1D, hist0, ps=1, file_out=file_out, xr=xr, $
;			;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6a1/Width (no unit)', overplot=overplot, col=col, $
;			title='!6n= n0+'+strtrim(i,1)+' l='+strtrim(v_i,1), xtitle='!6a1/Width (no unit)', overplot=overplot, col=col, $
;			fclose=fclose, show_stats=1, $
;			legendsize=legendsize, legend_precision=legend_precision
;
;		distrib_stat=estimate_1_sigma_error(x_hist_a1_ov_width[v_i,i,*],hist_a1_ov_width[v_i,i,*],95.5,2,tab_critere) 
;		Stat_synthese_a1_ov_width[0:6,v_i,i]=exp(distrib_stat[3:*]) ; les deciles, quartiles,mediane,...
;		Stat_synthese_a1_ov_width[7,v_i,i]=exp(distrib_stat[0]) ; le max
;	endfor	
;endfor
for v_i=0, lmax do begin
	for i=0,Nf_ls[v_i]-1 do begin
		log_a1w=alog(splitting_a1) - alog(width[v_i,i,*])
		loga1_ov_width[v_i, i, *]=log_a1w
		hist_a1_ov_width[v_i,i,*]=HISTOGRAM(log_a1w,nbins=1.*N_iter/Hist_length)
		x_hist_a1_ov_width[v_i,i,*]=findgen(n_elements(hist_a1_ov_width[v_i,i,*]))*(max(log_a1w)-min(log_a1w))/n_elements(hist_a1_ov_width[v_i,i,*])+min(log_a1w)
		hist0=dblarr(2, n_elements(x_hist_a1_ov_width[v_i,i,*]))
		hist0[0,*]=x_hist_a1_ov_width[v_i,i,*]
		hist0[1,*]=hist_a1_ov_width[v_i,i,*]
		file_out=output_ps+'logPDF_a1_ov_Width_n'+addzeros(i) + '_l' + strtrim(v_i,2) +'.eps'
		nice_hist1D, hist0, ps=1, file_out=file_out, xr=xr, $
			;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6a1/Width (no unit)', overplot=overplot, col=col, $
			title='!6n= n0+'+strtrim(i,1)+' l='+strtrim(v_i,1), xtitle='!6a1/Width (no unit)', overplot=overplot, col=col, $
			fclose=fclose, show_stats=1, $
			legendsize=legendsize, legend_precision=legend_precision

		distrib_stat=estimate_1_sigma_error(x_hist_a1_ov_width[v_i,i,*],hist_a1_ov_width[v_i,i,*],95.5,2,tab_critere) 
		Stat_synthese_a1_ov_width[0:6,v_i,i]=exp(distrib_stat[3:*]) ; les deciles, quartiles,mediane,...
		Stat_synthese_a1_ov_width[7,v_i,i]=exp(distrib_stat[0]) ; le max
	endfor	
endfor

info_samples='Values are calculated by alog(a1) - alog(width)'
save, loga1_ov_width, info_samples, filename=output_ps + 'Samples_a1_ov_Widths.sav'
info_synthese='all values of the histograms are shown in log. The Synthese array IS NOT in log'
save, Stat_Synthese_a1_ov_width, hist_a1_ov_width, x_hist_a1_ov_width, info_synthese, filename=output_ps + 'stat_synthese_a1_ov_Widths.sav'

width=0 & freq=0 & v=0
;*********************************************************
;**************** Affichage des Largeurs ***************
;*********************************************************
index=0 & j=0
!p.multi=0
;for i=0,Nmax-1 do begin
;	for v_i=0,lmax do begin
;		hist0=dblarr(2, n_elements(x_hist_width[v_i,i,*]))
;		hist0[0,*]=x_hist_width[v_i,i,*]
;		hist0[1,*]=hist_width[v_i,i,*]
;		file_out=output_ps+'PDF_Width_'+addzeros(index) +'.eps'
;		nice_hist1D, hist0, ps=1, file_out=file_out, xr=xr, $
;			;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
;			title='!6n= n0+'+strtrim(i,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
;			fclose=fclose, show_stats=1, $
;			legendsize=legendsize, legend_precision=legend_precision
;		index=index+1
;		wait,0.01
;	endfor
;endfor
for v_i=0,lmax do begin
	for i=0,Nf_ls[v_i]-1 do begin
		hist0=dblarr(2, n_elements(x_hist_width[v_i,i,*]))
		hist0[0,*]=x_hist_width[v_i,i,*]
		hist0[1,*]=hist_width[v_i,i,*]
		file_out=output_ps+'PDF_Width_'+addzeros(index) +'.eps'
		nice_hist1D, hist0, ps=1, file_out=file_out, xr=xr, $
			;title='!6n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
			title='!6n= n0+'+strtrim(i,1)+' l='+strtrim(v_i,1), xtitle='!6Width (microHz)', overplot=overplot, col=col, $
			fclose=fclose, show_stats=1, $
			legendsize=legendsize, legend_precision=legend_precision
		index=index+1
		wait,0.01
	endfor
endfor
;********************************************************
;******** CREATION DE LA SYNTHESE DES PARAMETRES ********
;********************************************************
;for i=0,Nmax-1 do begin
;	for j=0,lmax do begin
;		distrib_stat=estimate_1_sigma_error(x_hist_width[j,i,*],hist_width[j,i,*],95.5,2,tab_critere) ; dans l'ordre x,y puis tolerance (99.7% par def) et smooth valuewidth_proba[i]=distrib_stat[0]
;		Stat_synthese_width[0:6,j,i]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
;		Stat_synthese_width[7,j,i]=distrib_stat[0] ; le max
;	endfor
;endfor
for j=0,lmax do begin
	for i=0,Nf_ls[j]-1 do begin	
		distrib_stat=estimate_1_sigma_error(x_hist_width[j,i,*],hist_width[j,i,*],95.5,2,tab_critere) ; dans l'ordre x,y puis tolerance (99.7% par def) et smooth valuewidth_proba[i]=distrib_stat[0]
		Stat_synthese_width[0:6,j,i]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
		Stat_synthese_width[7,j,i]=distrib_stat[0] ; le max
	endfor
endfor
;------ Liberation de la memoire ------
hist_width=0 & x_hist_width=0

return, Stat_Synthese_width
end
