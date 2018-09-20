;fonction qui calcul la hauteur pour les l=1,l=2
; Recoit :
;		- file : noyau du nom de fichier contant les variables
;		- N_iter : nombre d'echantillons considérés
;		- Hist_length : nombre de classes pour l'histogramme
;		- output_ps : repertoire+noyau des fichers ps de sortie (si ps=1)
function MS_Global_height_l,file, lmax, Nmax, Nf_ls, N_iter,Hist_length,output_ps, rules, rules_params

	tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

;*******************************************************************
;*************** PRELIMINAIRES AU CALCUL DE HAUTEUR ***************
;*******************************************************************
Nmax_tab=max(Nf_ls) ; Compatibility with mixed modes models


h=dblarr(lmax+1,Nmax_tab,N_iter)
Stat_Synthese_Amplitude=dblarr(8,lmax+1,Nmax_tab) -1

if (N_iter ne 1) then begin

hist_amplitude=dblarr(lmax+1,Nmax_tab,N_iter/Hist_length)
x_hist_amplitude=dblarr(lmax+1,Nmax_tab,N_iter/Hist_length)

v=dblarr(lmax+1,N_iter)
v[0,*]=1.
for i=1,lmax do begin
	file_r=file_syntax(Nmax+i-1, file) ; detect automatically how many 0 I should put
	restore, file_r
	v[i,*]=param
endfor

;*********************************************************
;**************** CALCULS ET HISTOGRAMMES ****************
;*********************************************************
for i=0,Nmax-1 do begin
	file_hauteur=file_syntax(i, file) ; detect automatically how many 0 I should put

	restore,file_hauteur
	h=param

	hist_amplitude[0,i,*]=HISTOGRAM(h,nbins=1.*N_iter/Hist_length)
	x_hist_amplitude[0,i,*]=findgen(n_elements(hist_amplitude[0,i,*]))*(max(h)-min(h))/n_elements(hist_amplitude[0,i,*])+min(h)
	for v_i=1, lmax do begin
		if rules[v_i] eq 'p' then begin ; p modes with height defined by visibilities and H(l=0)
			h=h*v[v_i,*]
			hist_amplitude[v_i,i,*]=HISTOGRAM(h,nbins=1.*N_iter/Hist_length)
			x_hist_amplitude[v_i,i,*]=findgen(n_elements(hist_amplitude[v_i,i,*]))*(max(h)-min(h))/n_elements(hist_amplitude[v_i,i,*])+min(h)
		endif
	endfor
endfor
; --- rule for mixed modes ---
for v_i=1, lmax do begin
	if rules[v_i] eq 'm' then begin
		for i=0, Nf_ls[v_i]-1 do begin
			file_hauteur=file_syntax(rules_params[v_i,1]+i, file)
			restore,file_hauteur
			h=param
			hist_amplitude[v_i,i,*]=HISTOGRAM(h,nbins=1.*N_iter/Hist_length)
			x_hist_amplitude[v_i,i,*]=findgen(n_elements(hist_amplitude[0,i,*]))*(max(h)-min(h))/n_elements(hist_amplitude[0,i,*])+min(h)
			;stop
		endfor
	endif
endfor
; ----------------------------

hist_amplitude=hist_amplitude/N_iter
h=0 & h=0 & v=0

;*********************************************************
;**************** Affichage des Hauteurs ***************
;*********************************************************
!p.multi=0
index=0 & j=0
;for i=0,Nmax-1 do begin
;	for v_i=0,lmax do begin
;		file_out=output_ps+'PDF_Hl_'+ addzeros(index) +'.eps'
;		hist0=dblarr(2, n_elements(x_hist_amplitude[v_i,i,*]))
;		hist0[0,*]=x_hist_amplitude[v_i,i,*]
;		hist0[1,*]=hist_amplitude[v_i,i,*]
;		nice_hist1D, hist0, ps=1, file_out=file_out, xr=xr, $
;			;title='n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Height (ppm2/microHz)', overplot=overplot, col=col, $
;			title='n= n0 +'+strtrim(i,1)+' l='+strtrim(v_i,1), xtitle='!6Height (ppm2/microHz)', overplot=overplot, col=col, $
;			fclose=fclose, show_stats=1, $
;;			legendsize=legendsize, legend_precision=legend_precision
;		index=index+1
;		wait,0.01
;	endfor
;endfor
for v_i=0,lmax do begin
	for i=0,Nf_ls[v_i]-1 do begin
		file_out=output_ps+'PDF_Hl_'+ addzeros(index) +'.eps'
		hist0=dblarr(2, n_elements(x_hist_amplitude[v_i,i,*]))
		hist0[0,*]=x_hist_amplitude[v_i,i,*]
		hist0[1,*]=hist_amplitude[v_i,i,*]
		nice_hist1D, hist0, ps=1, file_out=file_out, xr=xr, $
			;title='n='+strtrim(i+n0-1,1)+' l='+strtrim(v_i,1), xtitle='!6Height (ppm2/microHz)', overplot=overplot, col=col, $
			title='n= n0 +'+strtrim(i,1)+' l='+strtrim(v_i,1), xtitle='!6Height (ppm2/microHz)', overplot=overplot, col=col, $
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
;distrib_stat=estimate_1_sigma_error(x_hist_amplitude[j,i,*],hist_amplitude[j,i,*],95.5,2,tab_critere) ; dans l'ordre x,y puis tolerance (99.7% par def) et smooth valueamplitude_proba[i]=distrib_stat[0]
;Stat_synthese_Amplitude[0:6,j,i]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
;Stat_synthese_Amplitude[7,j,i]=distrib_stat[0] ; le max
;	endfor
;endfor

for j=0,lmax do begin
	for i=0,Nf_ls[j]-1 do begin
		distrib_stat=estimate_1_sigma_error(x_hist_amplitude[j,i,*],hist_amplitude[j,i,*],95.5,2,tab_critere) ; dans l'ordre x,y puis tolerance (99.7% par def) et smooth valueamplitude_proba[i]=distrib_stat[0]
		Stat_synthese_Amplitude[0:6,j,i]=distrib_stat[3:*] ; les deciles, quartiles,mediane,...
		Stat_synthese_Amplitude[7,j,i]=distrib_stat[0] ; le max
	endfor
endfor


endif

;------ Liberation de la memoire ------
hist_amplitude=0 & x_hist_amplitude=0

return, Stat_Synthese_Amplitude
end
