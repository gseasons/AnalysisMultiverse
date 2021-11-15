Node: full_wm (utility)
=======================


 Hierarchy : brain_extraction_wf.full_wm
 Exec ID : full_wm


Original Inputs
---------------


* function_str : def _imsum(op1, op2, out_file=None):
    import nibabel as nb

    im1 = nb.load(op1)

    data = im1.get_fdata(dtype="float32") + nb.load(op2).get_fdata(dtype="float32")
    data /= data.max()
    nii = nb.Nifti1Image(data, im1.affine, im1.header)

    if out_file is None:
        from pathlib import Path

        out_file = str((Path() / "summap.nii.gz").absolute())

    nii.to_filename(out_file)
    return out_file

* op1 : /Users/grahamseasons/.cache/templateflow/tpl-OASIS30ANTs/tpl-OASIS30ANTs_res-01_label-WM_probseg.nii.gz
* op2 : /Users/grahamseasons/.cache/templateflow/tpl-OASIS30ANTs/tpl-OASIS30ANTs_res-01_label-BS_probseg.nii.gz
* out_file : <undefined>


Execution Inputs
----------------


* function_str : def _imsum(op1, op2, out_file=None):
    import nibabel as nb

    im1 = nb.load(op1)

    data = im1.get_fdata(dtype="float32") + nb.load(op2).get_fdata(dtype="float32")
    data /= data.max()
    nii = nb.Nifti1Image(data, im1.affine, im1.header)

    if out_file is None:
        from pathlib import Path

        out_file = str((Path() / "summap.nii.gz").absolute())

    nii.to_filename(out_file)
    return out_file

* op1 : /Users/grahamseasons/.cache/templateflow/tpl-OASIS30ANTs/tpl-OASIS30ANTs_res-01_label-WM_probseg.nii.gz
* op2 : /Users/grahamseasons/.cache/templateflow/tpl-OASIS30ANTs/tpl-OASIS30ANTs_res-01_label-BS_probseg.nii.gz
* out_file : <undefined>


Execution Outputs
-----------------


* out : /Users/grahamseasons/fMRI/analysis_multiverse/updated/preprocessing/brain_extraction_wf/full_wm/summap.nii.gz


Runtime info
------------


* duration : 0.856156
* hostname : Grahams-MBP.hitronhub.home
* prev_wd : /Users/grahamseasons/fMRI/analysis_multiverse/updated/preprocessing
* working_dir : /Users/grahamseasons/fMRI/analysis_multiverse/updated/preprocessing/brain_extraction_wf/full_wm


Environment
~~~~~~~~~~~


* ANTSPATH : /Users/grahamseasons/antsbin
* Apple_PubSub_Socket_Render : /private/tmp/com.apple.launchd.XbYuFiLAmb/Render
* CLICOLOR : 1
* CONDA_DEFAULT_ENV : base
* CONDA_EXE : /opt/anaconda3/bin/conda
* CONDA_PREFIX : /opt/anaconda3
* CONDA_PREFIX_1 : /Users/grahamseasons/anaconda3
* CONDA_PROMPT_MODIFIER : (base) 
* CONDA_PYTHON_EXE : /opt/anaconda3/bin/python
* CONDA_SHLVL : 2
* DISPLAY : /private/tmp/com.apple.launchd.FzhCKJ1Mzr/org.xquartz:0
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
* JPY_PARENT_PID : 1448
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
* SHELL : /bin/bash
* SHLVL : 2
* SPYDER_ARGS : ['--new-instance']
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
* SSH_AUTH_SOCK : /private/tmp/com.apple.launchd.0ucfl6iejx/Listeners
* TERM : xterm-color
* TERM_PROGRAM : Apple_Terminal
* TERM_PROGRAM_VERSION : 421.2
* TERM_SESSION_ID : 3849B059-386F-4A58-AD7A-B39848C972AE
* TMPDIR : /var/folders/mx/mztbckq95hzc7px9341hsc480000gn/T/
* USER : grahamseasons
* XPC_FLAGS : 0x0
* XPC_SERVICE_NAME : 0
* _ : /opt/anaconda3/python.app/Contents/MacOS/python
* _CE_CONDA : 
* _CE_M : 

