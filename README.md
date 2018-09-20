# README #


### IDL PostProcessor for TAMCMC C++ ###

This is a suite of tools that interpret binary files generated by TAMCMC-C++ (https://github.com/OthmanB/TAMCMC-C). 
This program will generate all kind of outputs in an IDL sav format and/or in ASCII: 
	* Parameters list and their uncertainties (frequencies, width height of modes, noise background...).
	* Plots for the best fit, residuals, echelle diagram and correlation maps..
	* Results for the evidence, pdfs, etc.. 

* Current version is compatible with TAMCMC C++ version 1.3.2 or below.

### How do I get set up? ###

* All setup is made by editing the PostMCMC_MS_Global.pro file. Only the procedure 'iterative_PostMCMC_MS_Global' might require to be edited to suits your need.

* Dependencies: There is many IDL library dependencies. I am not sure all those required have been provided here. If any issues, contact me and I will see if I can help.

* No exhaustive documentation is yet available. 
Here below I give very quick overview of the inputs to be set:

1. You need to copy the getmodel and bin2txt program into the cpp_prg subdirectory.  Those files are generated when compiling TACMCMC-C++.
2. You need to setup  directories and options of the 'iterative_PostMCMC_MS_Global' procedure. Here are the key parameters to set:
    
	- Directory containing the binary output files for all the objects that you wish to be processed:
    	dir_outputs='/home/obenomar/Pro/PSM_WP128_Sept2018/Raw_Results/TAMCMC-C/Data/Outputs/'
    	
	- Directory containing the inputs files for the TAMCMC-C++ analyis (.model and .data files):
    	dir_inputs='/home/obenomar/Pro/PSM_WP128_Sept2018/Data/Finalized-setups/Inputs/'
    	
    - Directory that will contain the results from the post processing. 
    	dir_out='/home/obenomar/Pro/PSM_WP128_Sept2018/Level1/'
    	
    - Name of the model that was used to perform the fit of the data (must be same as the one used in TAMCMC-C++).
    
    - modelname='model_MS_Global_a1etaa3_HarveyLike'


### Contribution guidelines ###

No external contribution is expected. This project is constantly improved, so please contact me if you need to see some worthy functionnality implemented. 

### Who do I talk to? ###

* Owner: Othman Benomar (NYUAD research associate)
* Contact: ob19@nyu.edu
