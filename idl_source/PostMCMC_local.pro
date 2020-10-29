@build_histogram_v2
@histograms_bin2txt_params
@MS_Global_rotinc_correlations
@MS_Global_amplitudes
@MS_Global_height_l
@MS_Global_width
@MS_Global_freq
@MS_Global_profile_plots
@MS_Global_freqspacing
@MS_Global_localnoise
@MS_echellediagram
@synthese2stddev
@MS_Global_fitplot
;@MS_Global_rotinc_correlations
@local_rotinc_correlations
@ascii2sav
@estimate_1_sigma_error
;@/home/ob19/Programs/IDL_library/mes_routines/kepler/read_ascii_kepler.pro
;@legend.pro/
@/home/obenomar/Dropbox/IDL/IDL_library/astro/pro/legend.pro
@draw_boxandwiskers.pro
;Post Processing for TAMCMC-C
; Version 1.3.2b : Compatible with TAMCMC-C version 1.3.2 and below
pro iterative_PostMCMC_Local;, dir_outputs, dir_inputs

   
   dosavfilesonly=0 ; Default, you do everything from sav file until synthese files and plots 
   
    ;OS='Linux' ; choice between Mac or Linux
    OS=''
    ;if OS eq 'Mac' then dir_os='/Volumes/MCMC_RES/'
    ;if OS eq 'Linux' then dir_os='/home/obenomar/MCMC_RES/'
    dir_OS=''

    ;dir_outputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass2/TAMCMC-C-1.4.3-dev/data/outputs/'
    ;dir_inputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass2/TAMCMC-C-1.4.3-dev/data/inputs/'
    ;dir_out=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass2/postmcmc/'
    ;dir_outputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass3/TAMCMC-C-1.4.3-dev/data/outputs/'
    ;dir_inputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass3/TAMCMC-C-1.4.3-dev/data/inputs/'
    ;dir_out=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass3/postmcmc/'
	dir_outputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass4/TAMCMC-C-1.4.3-dev/data/outputs/'
    dir_inputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass4/TAMCMC-C-1.4.3-dev/data/inputs/'
    dir_out=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass4/postmcmc/'
;	dir_outputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass5/TAMCMC-C-1.4.31-dev-EXP/data/outputs/'
;    dir_inputs=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass5/TAMCMC-C-1.4.31-dev-EXP/data/inputs/'
;    dir_out=dir_os + '/Volumes/home/2020/RGB-depressed/Level2/pass5/postmcmc/'

    modelname='model_MS_local_basic'

    dir_filter='*' ; Used to choose which directory should be processed
    phase='A*'
    Nb_classes=100.;
    index0=000. ; index of the first entry that is kept
    keep_period=10. ; Keep 1 out of keep_period
 	
    dirs=file_search(dir_outputs + dir_filter) ; lists the directories
    Ndirs=n_elements(dirs)
    Ok=intarr(Ndirs)
    
    i0=0
    for i=long(0), Ndirs-1 do begin
    	print, '[' + strtrim(i,2) + '] ' + dirs[i]
    endfor
    print, "type '.cont' to proceed"
    stop
    for i=long(i0), Ndirs-1 do begin
    ;for i=14, Ndirs-1 do begin
    	print, '    ------- Processing ' + dirs[i] + ' --------'    
    	b=byte(dirs[i])
    	pos=max(where(b eq 47 OR b eq 92)) ; detect slashes
        if pos[0] eq -1 then starID=dirs[i] else starID=strtrim(b[pos+1:*], 2)
        Nslices=n_elements(file_search(dirs[i] + '/outputs/' + starID + '_*_' + phase + 'params.hdr' ))
        for slice=long(1), Nslices do begin
    		done=PostMCMC_Local_extract_data(dir_outputs, dir_inputs, dir_out, modelname, starID, phase, strtrim(slice,2), Nb_classes, index0, keep_period, dosavfilesonly)
    	endfor
    	print, '    -------------------------------------------'
	OK[i]=done
	endfor
	posNotOK=where(OK eq 0)
	if n_elements(posNotOK) gt 1 then begin
		print, "List of process that could not be processed:"
		for i=0, n_elements(posNotOK)-1 do print, dirs[posNotOK[i]]
	endif else begin
		if posNotOK[0] eq -1 then print, 'All done'
		if posNotOK[0] ne -1 then begin
			print, "List of process that could not be processed:"
			print, dirs[posNotOK[0]]
		endif
	endelse
end


; Used to interpret the results from all Local models
; starID: the id as define in the config_presets.cfg
; phase: phase (B, L, A) as defined in the config_presets.cfg
; Nb_classes: Number of classes for the histograms (default = 200)
; index0: first index for the kept samples (default = 0)
; keep_period: periodicity for the the kept samples (default = 1  ==> all samples are kept)
; delete_txt_outputs: if 1, deletes the ascii files for the samples output files. Only sav files will be kepts (default = 0)
function PostMCMC_Local_extract_data, root_outputs, root_inputs, dir_out, modelname, starID, phase, slice, Nb_classes, index0, keep_period, dosavfilesonly
	done=1

	subdir=''	
	dir_bin2txt='cpp_prg/' ; directory where the function that converts binaries into ascii is.
	dir_getmodel='cpp_prg/' ; directory where the function that generate the models is

	; --- Defining the directory/files using the strict rule for managing inputs/outputs ----
	dir_IDL_out=dir_out + starID + '_' + slice + '/'
	binresultdir=root_outputs + starID + '/outputs/'
	diagsdir=root_outputs  + starID + '/diags/'
	root_filename=binresultdir + starID + '_' + slice + '_' + phase + '_params' ; here phase may contain and asterix
	test=file_search(root_filename + '_chain-0.bin')
	if test eq '' then begin
		print, 'Could not find files compatible with requested phase'
		print, 'Change the phase name. The program will stop now'
		done=0
		goto, bypass
		stop
	endif else begin
		root_filename=detect_root_name(test) ; ensure that we use a valid root_filename
	endelse
	;stop
	data_file=file_search(root_inputs + starID + '.data' )
	if data_file eq '' then begin
			print, 'Warning: Input data file not found.' 
			print, 'Check that the given directory and file exist'
			print, 'The program will stop now'
			stop
	endif
	
	check_dir=file_search(binresultdir) ; check if the dir has no capital letter
	if check_dir[0] eq '' then begin
			print, 'Warning: Output data files not found in the specified directory.' 
			print, 'Check that you provided the correct output directory'
			print, 'The program will stop now'
			stop
	endif
	
	f=file_search(dir_IDL_out + '')
	if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + ''

	; ----- Convert binary files in text files, easier to read by IDL AND plot their pdf using nice1D_hist ----
	print, 'Convert binary files in text files and then into sav files, easier to read by IDL and plot their pdf...'
	f=file_search(dir_IDL_out + 'Files/')
	if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + 'Files/'	
	hists_info=histograms_bin2txt_params(dir_IDL_out+'Files/',  Nb_classes, root_filename, dir_bin2txt,  modelname,index0=index0, keep_period=keep_period)
	Nsamples=hists_info.Nsamples

	stat_synthese_unsorted=hists_info.stat_synthese_unsorted	

	; ---- Extract the value of parameters_length ----
	parameters_length=read_plength(dir_IDL_out +'Files/plength.txt')
	Nmax=parameters_length[0]
	lmax=3 ;parameters_length[1]
	Nf_ls=parameters_length[2:5] ; Nfl0, Nfl1, Nfl2, Nfl3
	ind0_nu=total(parameters_length[0:1])
	ind0_split=total(parameters_length[0:5])
	ind0_W=total(parameters_length[0:6])
	ind0_inc=total(parameters_length[0:8])
	ind0_trunc=total(parameters_length[0:9])
	

	; ---- Determine the Global Likelihood -----
	print, 'Evidence is assumed to be calculated by the CPP program: copy of the file in the IDL output directory...'
	if phase eq 'B*' then ph='B'
	if phase eq 'L*' then ph='L'
	if phase eq 'A*' then ph='A'

	spawn, 'cp ' + diagsdir + starID + '_' + slice + '_' + ph + '_evidence.txt ' + dir_IDL_out + starID + '_' + slice + '_' + ph + '_evidence.txt'

	; ---------- Plot of models ----
	; Skipping it because the diags are already showing model + data

	; --- NEED WORKING FROM HERE ----
	
	; ----- Correlation diagram for Rotation/Inclination -----
	print, 'Plots for correlations between a1,a2,a3,asymetry, magnetic effects parametrisation and inclination...'
	local_rotinc_correlations, parameters_length, dir_IDL_out +'Files/', dir_IDL_out, starID,0 ; can handle any number of correlation maps within the same plot
	;stop
	; ---- Save basic information on splitting and inclination in a text file ----
	show_split_inc_asym, dir_IDL_out, modelname, stat_synthese_unsorted, ind0_split, ind0_W, ind0_inc,Nmax
	bypass:
	return, done
end


; Used to interpret the results from all MS_Global models
; starID: the id as define in the config_presets.cfg
; phase: phase (B, L, A) as defined in the config_presets.cfg
; Nb_classes: Number of classes for the histograms (default = 200)
; index0: first index for the kept samples (default = 0)
; keep_period: periodicity for the the kept samples (default = 1  ==> all samples are kept)
; delete_txt_outputs: if 1, deletes the ascii files for the samples output files. Only sav files will be kepts (default = 0)
function PostMCMC_MS_Global_v2, root_outputs, root_inputs, dir_out, modelname, starID, phase, Nb_classes, index0, keep_period, dosavfilesonly
	done=1

	subdir=''	
	dir_bin2txt='cpp_prg/' ; directory where the function that converts binaries into ascii is.
	dir_getmodel='cpp_prg/' ; directory where the function that generate the models is

	; --- Defining the directory/files using the strict rule for managing inputs/outputs ----
	dir_IDL_out=dir_out + starID +'/'
	binresultdir=root_outputs + starID + '/outputs/'
	diagsdir=root_outputs  + starID + '/diags/'
	root_filename=binresultdir + starID + '_' + phase + '_params' ; here phase may contain and asterix
	test=file_search(root_filename + '_chain-0.bin')
	if test eq '' then begin
		print, 'Could not find files compatible with requested phase'
		print, 'Change the phase name. The program will stop now'
		done=0
		goto, bypass
		stop
	endif else begin
		root_filename=detect_root_name(test) ; ensure that we use a valid root_filename
	endelse
	;stop
	data_file=file_search(root_inputs + starID + '.data' )
	if data_file eq '' then begin
			print, 'Warning: Input data file not found.' 
			print, 'Check that the given directory and file exist'
			print, 'The program will stop now'
			stop
	endif
	
	check_dir=file_search(binresultdir) ; check if the dir has no capital letter
	if check_dir[0] eq '' then begin
			print, 'Warning: Output data files not found in the specified directory.' 
			print, 'Check that you provided the correct output directory'
			print, 'The program will stop now'
			stop
	endif
	
	f=file_search(dir_IDL_out + '')
	if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + ''

	; ----- Convert binary files in text files, easier to read by IDL AND plot their pdf using nice1D_hist ----
	print, 'Convert binary files in text files and then into sav files, easier to read by IDL and plot their pdf...'
	f=file_search(dir_IDL_out + 'Files/')
	if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + 'Files/'	
	hists_info=histograms_bin2txt_params(dir_IDL_out+'Files/',  Nb_classes, root_filename, dir_bin2txt,  modelname,index0=index0, keep_period=keep_period)

	if dosavfilesonly eq 0 then begin ; We may just want to convert files into save files and make just the basic plots (usually when limited free parameters)
		Nsamples=hists_info.Nsamples

		stat_synthese_unsorted=hists_info.stat_synthese_unsorted
	;	stop
		; ---- Extract the value of parameters_length ----
		parameters_length=read_plength(dir_IDL_out +'Files/plength.txt')
		Nmax=parameters_length[0]
		lmax=parameters_length[1]
		Nf_ls=parameters_length[2:5] ; Nfl0, Nfl1, Nfl2, Nfl3
		ind0_nu=total(parameters_length[0:1])
		ind0_split=total(parameters_length[0:5])
		ind0_W=total(parameters_length[0:6])
		ind0_inc=total(parameters_length[0:8])
		ind0_trunc=total(parameters_length[0:9])

		rules=['p','p','p','p'] ; defining default rules (p mode fitting for all modes)
		rules_params=dblarr(lmax+1, 4); no parameters for Frequencies, Heights, Width, splitting
		exclude_sep=0
		if modelname eq 'model_Evolved_Global_a1etaa3_l1mixed' then begin
			rules[1]='m' ; specify that l=1 are mixed in the case of this model
			rules_params[1,*]=[-1, total(parameters_length[0:10]), total(parameters_length[0:11]), total(parameters_length[0:12])] ; define the rules for l=1 mixed modes
			exclude_sep=1
		endif
		;print, 'debug stop for the rule definition'
		;stop


		; ---- Calculation of Amplitudes, Heights, Widths for each mode ----
		print, 'Determining Frequencies, Amplitudes, Heights and Widths for each mode...'
		f=file_search(dir_IDL_out + 'Frequencies/')
		if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + 'Frequencies/'
		print, 'Processing Frequencies...'				
		Stat_Synthese_Freq=MS_Global_freq(dir_IDL_out + 'Files/',ind0_nu, lmax, Nmax, Nf_ls, Nsamples,Nb_classes, dir_IDL_out + 'Frequencies/')
		f=file_search(dir_IDL_out + 'Amplitudes/')
		if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + 'Amplitudes/'
		print, 'Processing Amplitudes...'				
		Stat_Synthese_Amplitude=MS_Global_amplitudes(dir_IDL_out + 'Files/', ind0_W, lmax, Nmax, Nf_ls, Nsamples,Nb_classes,dir_IDL_out + 'Amplitudes/', rules, rules_params)
		f=file_search(dir_IDL_out + 'Heights/')
		if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + 'Heights/'
		print, 'Processing Heights...'
		Stat_Synthese_Height=MS_Global_height_l(dir_IDL_out + 'Files/', lmax, Nmax, Nf_ls, Nsamples,Nb_Classes, dir_IDL_out + 'Heights/', rules, rules_params)
		f=file_search(dir_IDL_out + 'Widths/')
		if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + 'Widths/'
		print, 'Processing Widths...'
		save_Samples=0; Switch to 1 if you need to keep the samples (Normaly just required for Inertia calculation)
		Stat_Synthese_width=MS_Global_width(dir_IDL_out + 'Files/', ind0_nu, ind0_W, ind0_split, lmax, Nmax, Nf_ls, Nsamples,Nb_classes, $
					    		dir_IDL_out + 'Widths/', rules, rules_params, save_Samples=save_Samples)
		print, 'Computing noise level below modes...'
		local_noise=MS_Global_localnoise(stat_synthese_freq, stat_synthese_unsorted, parameters_length)
		
		; ---- Plots of the profiles ----
		print, 'Ploting the profiles for the Heights, Widths and Amplitudes...'
		MS_Global_profile_plots, Stat_Synthese_freq,Stat_Synthese_Height,Stat_synthese_width,Stat_synthese_Amplitude,$
			local_noise,dir_IDL_out + ''
					
		; ---- Save the synthese results ----
		synthese_file=dir_IDL_out + 'synthese.sav' 
		bruit_local=local_noise ; for compatibility with my old programs
		print, 'Saving Stat_Synthese_freq, Stat_Synthese_Height, Stat_Synthese_Width, Stat_Synthese_Amplitude, Stat_Synthese_Unsorted, ...'	
		save, Stat_Synthese_freq, Stat_Synthese_Height, Stat_Synthese_Width, Stat_Synthese_Amplitude, Stat_Synthese_Unsorted, $
				parameters_length, bruit_local, local_noise, filename=synthese_file
	
		; ---- Determine the Global Likelihood -----
		print, 'Evidence is assumed to be calculated by the CPP program: copy of the file in the IDL output directory...'
		if phase eq 'B*' then ph='B'
		if phase eq 'L*' then ph='L'
		if phase eq 'A*' then ph='A'

		spawn, 'cp ' + diagsdir + starID + '_' + ph + '_evidence.txt ' + dir_IDL_out + ph + '_evidence.txt'
	
		; ---- Generate a short text output with relevant values for all the modes ----
		print, 'Computing frequency spacings...'
		if exclude_sep eq 0 then begin
			Dnus=MS_Global_freqspacing(synthese_file,dir_IDL_out, Nmax, lmax, show=0)
			ass_law_file=dir_IDL_out+'freq_spacings.sav'
		endif else begin
			l=0
			errors=sqrt( (stat_synthese_freq[3,l,0:Nmax-1]-stat_synthese_freq[2,l,0:Nmax-1])^2 + (stat_synthese_freq[4,l,0:Nmax-1]-stat_synthese_freq[3,l,0:Nmax-1])^2 ) / sqrt(2)
			if total(errors) ge 1d-3 then c=linfit(findgen(n_elements(stat_synthese_freq[3,l,0:Nmax-1])), reform(stat_synthese_freq[3,l,0:Nmax-1]), MEASURE_ERRORS=errors, sigma=err_asslaw)
			if total(errors) lt 1d-3 then c=linfit(findgen(n_elements(stat_synthese_freq[3,l,0:Nmax-1])), reform(stat_synthese_freq[3,l,0:Nmax-1]))
			Dnus=c[1]
			ass_law_file=''
		endelse
		print, 'Creating  ultra-simplistic outputs for modes parameters and generate global mode characteristic summary (e.g. \Delta\nu, ...)...'
		file_out=dir_IDL_out+'freq_std_dev'
		synthese2stddev,synthese_file,ass_law_file,file_out
	
		; ---- Save the model and the median parameters----
		do_model=1
		if do_model eq 1 then begin
			print, 'Modeled spectrum and median values...'	
			val_med=reform(stat_synthese_unsorted[3, *])
			val_med_m1s=reform(stat_synthese_unsorted[2, *]) ; median - 1sigma
			val_med_p1s=reform(stat_synthese_unsorted[4, *]) ; median + 1sigma	
			; Write a configuration file suitable for the compute_model.cpp function
			params_cfg= dir_IDL_out + 'best_models_params.txt'
			file_out=dir_IDL_out + 'best_models_fit.ascii'
			file_psfit=dir_IDL_out + 'best_models_fit.eps'
			openw, 3, params_cfg
				str='# This file contains in the first line, the parameters structure (plength). All following lines, correspond to a single vector of parameters'
				printf, 3, str
				str=''
				for i=long(0), n_elements(parameters_length)-1 do str=str+ '   ' + strtrim(parameters_length[i],2)
				printf, 3, str
				str=''
				for i=long(0), n_elements(val_med)-1 do str=str+ '   ' + strtrim(val_med[i],2)
				printf, 3, str
				str=''
				for i=long(0), n_elements(val_med_m1s)-1 do str=str+ '   ' + strtrim(val_med_m1s[i],2)
				printf, 3, str
				str=''
				for i=long(0), n_elements(val_med_p1s)-1 do str=str+ '   ' + strtrim(val_med_p1s[i],2)
				printf, 3, str
			close, 3
			; Use the in-built function of TAMCMC to get the median model
			spawn, dir_getmodel +'./getmodel ' + data_file + ' ' +  params_cfg + ' ' + modelname + ' ' + file_out
			;stop
			; read the file that was just created
			Ncols=detect_Ncolumns(file_out, skip=0)
			;Ncols=5 ; col[0]=freq, col[1]=spec_reg, col[2]=median_model, col[3]=median_minus1sigma_model, col[4]=median_plus1signa_model
			model_bestfit=read_Ncolumns(file_out, Ncols, 5d5, skip=0, ref_N=0)
			fmin=min(model_bestfit[0,*])
			fmax=max(model_bestfit[0,*])
			MS_Global_fitplot, model_bestfit[0, *], model_bestfit[1, *], model_bestfit[0, *], model_bestfit[2:*, *], Dnus[0], file_psfit, fmin, fmax
			;stop
	
			; ----- Echelle diagram ----
			params_cfg_noise= dir_IDL_out + 'noise_model_params.txt'
			file_noise_fit=dir_IDL_out + 'noise_models_fit.ascii'
			openw, 3, params_cfg_noise
				str='# This file contains in the first line, the parameters structure (plength). All following lines, correspond to a single vector of parameters'
				printf, 3, str
				str=''
				for i=long(0), n_elements(parameters_length)-1 do str=str+ '   ' + strtrim(parameters_length[i],2)
				printf, 3, str
				str=''
				val_med_noise=val_med
				val_med_noise[0:parameters_length[0]-1]=0 ; put heights to 0
				for i=long(0), n_elements(val_med_noise)-1 do str=str+ '   ' + strtrim(val_med_noise[i],2)
				printf, 3, str
			close, 3
			; Use the in-built function of TAMCMC to get the median model with noise only
			spawn, dir_getmodel +'./getmodel ' + data_file + ' ' +  params_cfg_noise + ' ' + modelname + ' ' + file_noise_fit
			; read the file that was just created
			;Ncols=5 ; col[0]=freq, col[1]=spec_reg, col[2]=median_model_noise
			Ncols=detect_Ncolumns(file_noise_fit, skip=0)
			model_noise=read_Ncolumns(file_noise_fit, Ncols, 5d5, skip=0, ref_N=0)
	
			file_ps_echelle=dir_IDL_out + 'Echelle_Diagram.eps' 
			nimp,name=file_ps_echelle,/paper,/eps
				show_ech_diag_v3_CPP, synthese_file, model_noise, ps=1, shifts=shifts, trunc_spec=trunc_spec
			fimp
			file_ps_echelle=dir_IDL_out + 'Echelle_Diagram_residuals.eps' 
			nimp,name=file_ps_echelle,/paper,/eps
				show_ech_diag_v3_CPP, synthese_file, model_bestfit[0:2,*], ps=1, shifts=shifts, trunc_spec=trunc_spec
			fimp
		endif
		; ----- Correlation diagram for Rotation/Inclination -----
		print, 'Plots for correlations between a1,a2,a3,asymetry, magnetic effects parametrisation and inclination...'
		MS_Global_rotinc_correlations, synthese_file, dir_IDL_out +'Files/', dir_IDL_out, starID,0 ; can handle any number of correlation maps within the same plot

		stop
		; ---- Save basic information on splitting and inclination in a text file ----
		show_split_inc_asym, dir_IDL_out, modelname, stat_synthese_unsorted, ind0_split, ind0_W, ind0_inc,Nmax
	endif ; endif for dosavfileonly
	bypass:
	return, done
end



; A procedure that deal with the various situations specific to models
pro show_split_inc_asym, dir_IDL_out, modelname, stat_synthese_unsorted, ind0_split, ind0_W, ind0_inc, Nmax

		split_vals=Stat_Synthese_Unsorted[0:6,ind0_split:ind0_W-1]
		inc_vals=Stat_Synthese_Unsorted[0:6,ind0_inc]

		if n_elements(split_vals[0,*]) eq 6 then split_names=['a1          ', 'eta         ', 'a3          ', 'sqrt(a1).cos(i)       ', 'sqrt(a1).sin(i)   ', 'asym        ']
		if n_elements(split_vals[0,*]) eq 7 then split_names=['a1(l=1)     ', 'eta         ', 'a3          ', 'empty slot       ', 'empty slot  ', 'asym        ', 'a1(l=2)']
		if n_elements(split_vals[0,*]) eq 6+Nmax then begin
			a1_names=strarr(Nmax)
			for i=long(0), Nmax-1 do a1_names[i]='a1(n0 +' + strtrim(i,2) + ')'
			split_names=['<a1>=0     ', 'eta         ', 'a3          ', 'empty slot       ', 'empty slot   ', 'asym        ', a1_names]
		endif
		if n_elements(split_vals[0,*]) eq 6+2*Nmax then begin
			a1_names=strarr(2*Nmax)
			for i=long(0), Nmax-1 do a1_names[i]='a1(n0 +' + strtrim(i,2) + ', l=1)'
			for i=long(0), Nmax-1 do a1_names[Nmax+i]='a1(n0 +' + strtrim(i,2) + ', l=2)'
			split_names=['<a1>=0     ', 'eta         ', 'a3          ', 'empty slot       ', 'empty slot   ', 'asym        ', a1_names]
		endif
		file_split_inc=dir_IDL_out + 'split_inc.txt'
		openw, 3, file_split_inc
			printf, 3, '# Short summary files with the splitting values and the stellar inclination'
			printf, 3, '# Format of the outputs'
			printf, 3, '# [Param name]   0%     2.25%   16%     50%    84%    97.75%   100%'
			for i=0, n_elements(split_vals[0,*])-1 do begin		
				str=split_names[i]
				for j=0, n_elements(split_vals[*,0])-1 do begin
					str=str + string(split_vals[j, i], format='(f20.9)')
				endfor
				printf, 3, str
			endfor
			str='inc         '
			for j=0, n_elements(inc_vals[*,0])-1 do str=str + string(inc_vals[j], format='(f20.9)')
			printf, 3, str
		close,3

	
end


; ------------------------------------------------------------------------------------
; A function that translate the variable names of the MCMC code into some nicer names
function interpret_varnames_Cpp, varname;, stopping

variable_name=strarr(4)
	
		;if stopping eq 1 then stop
		
		if varname eq 'Height_l' then begin
			variable_name[0]='!6Height (ppm!U2!N/!7l!m!6Hz)'
			variable_name[1]='!6l=0, n=?'
		endif
		if varname eq 'Height_l1' then begin
			variable_name[0]='!6Height (ppm!U2!N/!7l!m!6Hz)'
			variable_name[1]='!6l=1, n=?'
		endif
		if varname eq 'Frequency_l' then begin
			variable_name[0]='!6Frequency (!7l!X!6Hz)'
			variable_name[1]='!6l=? n=?'
		endif
		if varname eq 'Splitting_a1' then begin
			variable_name[0]='!6Rotational Splitting a1 (!7l!X!6Hz)'
			variable_name[1]=''
		endif
		if varname eq 'Splitting_a11' then begin
			variable_name[0]='!6Rotational Splitting a1(l=1) (!7l!X!6Hz)'
			variable_name[1]=''
		endif
		if varname eq 'Asphericity_eta' then begin
			variable_name[0]='!6Asphericity eta (no unit)'
			variable_name[1]=''
		endif
		if varname eq 'Splitting_a3' then begin
			variable_name[0]='!6Coefficient a3 (!7l!X!6Hz)'
			variable_name[1]=''
		endif		
		if varname eq '\sqrt(a1).cos(i)' then begin ; CHECK WHY WE DONT GET IN THIS
			variable_name[0]='!6Asphericity parameter Bmag (!7l!X!6Hz)'
			variable_name[1]=''
		endif		
		if varname eq '\sqrt(a1).sin(i)' then begin ; CHECK WHY WE DONT GET IN THIS
			variable_name[0]='!6Asphericity parameter alfa (no unit)'
			variable_name[1]=''
		endif		
		if varname eq 'Lorentzian_asymetry' then begin
			variable_name[0]='!6Lorentzian asymetry B (no unit)'
			variable_name[1]=''
		endif		
		if varname eq 'Inclination' then begin
			variable_name[0]='!6Stellar inclination (degree)'
			variable_name[1]=''
		endif		
		if varname eq 'Harvey-Noise_H' then begin
			variable_name[0]='!6Harvey Noise parameter H'
			variable_name[1]=''
		endif	
		if varname eq 'Harvey-Noise_tc' then begin
			variable_name[0]='!6Harvey Noise parameter tc (ksec)'
			variable_name[1]=''
		endif	
		if varname eq 'Harvey-Noise_p' then begin
			variable_name[0]='!6Harvey Noise parameter p (no unit)'
			variable_name[1]=''
		endif	
		if varname eq 'White_Noise_N0' then begin
			variable_name[0]='!6White noise parameter N0 (ppm!U2!N/!7l!X!6Hz)'
			variable_name[1]=''
		endif	
		if varname eq 'Width_l' then begin
			variable_name[0]='!6FWHM (!7l!X!6Hz)'
			variable_name[1]='!6l=0, n=?'
		endif
		if varname eq 'Width_l1' then begin
			variable_name[0]='!6FWHM(!7l!X!6Hz)'
			variable_name[1]='!6l=1, n=?'
		endif
return, variable_name
end

; A short function that read the plength.txt files generated by bin2txt
;function read_plength, file
;	a=''
;	openr, 3, file
;		readf, 3, a
;	close,3
;	uu=strsplit(a)
;    N_uu=N_elements(uu)-1
;    plength=intarr(N_uu)
;    for j=0,N_uu-1 do begin
;         plength(j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
;   endfor
;	plength(N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
;	
;	return, plength
;end

; A short function that read the plength.txt files generated by bin2txt
function read_plength, file
	a=''
	i=0.
	openr, 3, file
		while (eof(3) ne LOGICAL_TRUE(1)) do begin
			 readf, 3, a
			 if i eq 0 then plength=double(a) else plength=[plength, double(a)]
			 i=i+1
		endwhile
	close,3
	return, plength
end

; Determine how many columns exists in a file
function detect_Ncolumns, file, skip=skip

if n_elements(skip) eq 0 then skip=1 ; defaut we skip one line only

openr, 3, file

	
	a=''
	i=0d
    while i le skip +1 do begin
      	if i lt skip then readf,3,format='(q)'
        if i ge skip then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)
        endif
        i=i+1
    endwhile
close, 3

	;stop
return, N_uu
end

; N: number of columns
; K: maximum number of lines
function read_Ncolumns, file,N, K, skip=skip, ref_N=ref_N

if n_elements(skip) eq 0 then skip=1 ; defaut we skip one line only
if n_elements(ref_N) eq 0 then ref_N=1 ; defaut we identify the zero non-used tab elements with column 1
if n_elements(spectrum) eq 0 then spectrum=1
openr, 3, file

	param=dblarr(K,N)
	a=''
	i=0d
      while EOF(3) ne 1 do begin
      	if i lt skip then readf,3,format='(q)'
        if i ge skip then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)-1
          	for j=0,N_uu-1 do begin
          		param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          	endfor
			param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		endif
		i=i+1
      endwhile

close,3
param0=param
test=where(param[*,ref_N] ne 0)

param=param[test,*]
param=transpose(param)

return,param
end

function file_syntax, in, core

	if in lt 10 then file_out=core+'00'+strtrim(round(in),1)+'.sav'
	if in ge 10 AND in lt 100 then file_out=core+'0'+strtrim(round(in),1)+'.sav'
	if in ge 100 then file_out=core +strtrim(round(in),1)+'.sav'
	
return, file_out
end


function addzeros, in

	if in lt 10 then out='00'+strtrim(round(in),2)
	if in ge 10 AND in lt 100 then out='0'+strtrim(round(in),2)
	if in ge 100 then out=strtrim(round(in),2)
	
	return, out
end

; cut the file name just before _chain-*.bin
function detect_root_name, file_chain_bin

	b=byte(file_chain_bin)
	pos=max(where(b eq 95)) ;detect last '_'
	name=strtrim(b[0:pos-1],2)
	return, name
end
