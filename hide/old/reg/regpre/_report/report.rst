Node: regpre (fsl)
==================


 Hierarchy : reg.regpre
 Exec ID : regpre


Original Inputs
---------------


* angle_rep : <undefined>
* apply_isoxfm : <undefined>
* apply_xfm : <undefined>
* args : <undefined>
* bbrslope : <undefined>
* bbrtype : <undefined>
* bgvalue : <undefined>
* bins : <undefined>
* coarse_search : <undefined>
* cost : <undefined>
* cost_func : <undefined>
* datatype : <undefined>
* display_init : <undefined>
* dof : 6
* echospacing : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* fieldmap : <undefined>
* fieldmapmask : <undefined>
* fine_search : <undefined>
* force_scaling : <undefined>
* in_file : /Volumes/NewVolume/eggs/sub-094S6278_task-rest_bold_roi_mcf_mean_reg.nii.gz
* in_matrix_file : <undefined>
* in_weight : <undefined>
* interp : sinc
* min_sampling : <undefined>
* no_clamp : <undefined>
* no_resample : True
* no_resample_blur : <undefined>
* no_search : <undefined>
* out_file : <undefined>
* out_log : <undefined>
* out_matrix_file : <undefined>
* output_type : NIFTI_GZ
* padding_size : <undefined>
* pedir : <undefined>
* ref_weight : <undefined>
* reference : /Volumes/NewVolume/eggs/sub-094S6278_T1w_corrected_xform_masked.nii.gz
* rigid2D : <undefined>
* save_log : <undefined>
* schedule : <undefined>
* searchr_x : <undefined>
* searchr_y : <undefined>
* searchr_z : <undefined>
* sinc_width : <undefined>
* sinc_window : <undefined>
* uses_qform : <undefined>
* verbose : <undefined>
* wm_seg : <undefined>
* wmcoords : <undefined>
* wmnorms : <undefined>


Execution Inputs
----------------


* angle_rep : <undefined>
* apply_isoxfm : <undefined>
* apply_xfm : <undefined>
* args : <undefined>
* bbrslope : <undefined>
* bbrtype : <undefined>
* bgvalue : <undefined>
* bins : <undefined>
* coarse_search : <undefined>
* cost : <undefined>
* cost_func : <undefined>
* datatype : <undefined>
* display_init : <undefined>
* dof : 6
* echospacing : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* fieldmap : <undefined>
* fieldmapmask : <undefined>
* fine_search : <undefined>
* force_scaling : <undefined>
* in_file : /Volumes/NewVolume/eggs/sub-094S6278_task-rest_bold_roi_mcf_mean_reg.nii.gz
* in_matrix_file : <undefined>
* in_weight : <undefined>
* interp : sinc
* min_sampling : <undefined>
* no_clamp : <undefined>
* no_resample : True
* no_resample_blur : <undefined>
* no_search : <undefined>
* out_file : <undefined>
* out_log : <undefined>
* out_matrix_file : <undefined>
* output_type : NIFTI_GZ
* padding_size : <undefined>
* pedir : <undefined>
* ref_weight : <undefined>
* reference : /Volumes/NewVolume/eggs/sub-094S6278_T1w_corrected_xform_masked.nii.gz
* rigid2D : <undefined>
* save_log : <undefined>
* schedule : <undefined>
* searchr_x : <undefined>
* searchr_y : <undefined>
* searchr_z : <undefined>
* sinc_width : <undefined>
* sinc_window : <undefined>
* uses_qform : <undefined>
* verbose : <undefined>
* wm_seg : <undefined>
* wmcoords : <undefined>
* wmnorms : <undefined>


Execution Outputs
-----------------


* out_file : <undefined>
* out_log : <undefined>
* out_matrix_file : /Users/grahamseasons/fMRI/analysis_multiverse/hide/old/reg/regpre/sub-094S6278_task-rest_bold_roi_mcf_mean_reg_flirt.mat


Runtime info
------------


* cmdline : flirt -in /Volumes/NewVolume/eggs/sub-094S6278_task-rest_bold_roi_mcf_mean_reg.nii.gz -ref /Volumes/NewVolume/eggs/sub-094S6278_T1w_corrected_xform_masked.nii.gz -out sub-094S6278_task-rest_bold_roi_mcf_mean_reg_flirt.nii.gz -omat sub-094S6278_task-rest_bold_roi_mcf_mean_reg_flirt.mat -dof 6 -interp sinc -noresample
* duration : 677.339957
* hostname : Grahams-MBP.hitronhub.home
* prev_wd : /Users/grahamseasons/fMRI/analysis_multiverse/hide/old
* working_dir : /Users/grahamseasons/fMRI/analysis_multiverse/hide/old/reg/regpre


Terminal output
~~~~~~~~~~~~~~~


 


Terminal - standard output
~~~~~~~~~~~~~~~~~~~~~~~~~~


 


Terminal - standard error
~~~~~~~~~~~~~~~~~~~~~~~~~


 


Environment
~~~~~~~~~~~


* ANTSPATH : /Users/grahamseasons/antsbin
* Apple_PubSub_Socket_Render : /private/tmp/com.apple.launchd.2Sycx0lCyL/Render
* CLICOLOR : 1
* CONDA_DEFAULT_ENV : base
* CONDA_EXE : /opt/anaconda3/bin/conda
* CONDA_PREFIX : /opt/anaconda3
* CONDA_PREFIX_1 : /Users/grahamseasons/anaconda3
* CONDA_PROMPT_MODIFIER : (base) 
* CONDA_PYTHON_EXE : /opt/anaconda3/bin/python
* CONDA_SHLVL : 2
* DISPLAY : /private/tmp/com.apple.launchd.E1gxMFPioX/org.xquartz:0
* EVENT_NOKQUEUE : 1
* FSLDIR : /usr/local/fsl
* FSLGECUDAQ : cuda.q
* FSLLOCKDIR : 
* FSLMACHINELIST : 
* FSLMULTIFILEQUIT : TRUE
* FSLOUTPUTTYPE : NIFTI_GZ
* FSLREMOTECALL : 
* FSLTCLSH : /usr/local/fsl/bin/fsltclsh
* FSLWISH : /usr/local/fsl/bin/fslwish
* GIT_PAGER : cat
* HOME : /Users/grahamseasons
* JPY_PARENT_PID : 45754
* LANG : en_CA.UTF-8
* LANGUAGE : en
* LC_ALL : en_CA.UTF-8
* LOGNAME : grahamseasons
* MATLABCMD : /Applications/MATLAB_R2017b.app/bin/matlab
* MPLBACKEND : module://ipykernel.pylab.backend_inline
* PAGER : cat
* PATH : /Users/grahamseasons/antsbin/bin:/usr/local/fsl/bin:/Library/Frameworks/Python.framework/Versions/3.8/bin:/Library/Frameworks/Python.framework/Versions/3.8/bin:/opt/anaconda3/bin:/Users/grahamseasons/anaconda3/condabin:/Applications/MATLAB_R2017b.app/bin:/opt/local/bin:/opt/local/sbin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/Users/grahamseasons/abin:/Users/grahamseasons/abin
* PWD : /Users/grahamseasons
* PYTHONPATH : /Users/grahamseasons/fMRI/analysis_multiverse
* QT_SCALE_FACTOR : 
* QT_SCREEN_SCALE_FACTORS : 
* SECURITYSESSIONID : 186a7
* SHELL : /bin/bash
* SHLVL : 2
* SPYDER_ARGS : []
* SPY_AUTOCALL_O : 0
* SPY_AUTOLOAD_PYLAB_O : False
* SPY_BACKEND_O : 8
* SPY_BBOX_INCHES_O : True
* SPY_EXTERNAL_INTERPRETER : False
* SPY_FORMAT_O : 0
* SPY_GREEDY_O : False
* SPY_HEIGHT_O : 4
* SPY_HIDE_CMD : True
* SPY_JEDI_O : False
* SPY_PYLAB_O : True
* SPY_RESOLUTION_O : 72
* SPY_RUN_FILE_O : 
* SPY_RUN_LINES_O : 
* SPY_SYMPY_O : False
* SPY_TESTING : None
* SPY_UMR_ENABLED : True
* SPY_UMR_NAMELIST : 
* SPY_UMR_VERBOSE : True
* SPY_USE_FILE_O : False
* SPY_WIDTH_O : 6
* SSH_AUTH_SOCK : /private/tmp/com.apple.launchd.6Jajo3OBii/Listeners
* TERM : xterm-color
* TERM_PROGRAM : Apple_Terminal
* TERM_PROGRAM_VERSION : 421.2
* TERM_SESSION_ID : 5B213730-D500-4A4E-A3B8-0EBA31D46115
* TMPDIR : /var/folders/mx/mztbckq95hzc7px9341hsc480000gn/T/
* USER : grahamseasons
* XPC_FLAGS : 0x0
* XPC_SERVICE_NAME : 0
* _ : /opt/anaconda3/python.app/Contents/MacOS/python
* _CE_CONDA : 
* _CE_M : 

