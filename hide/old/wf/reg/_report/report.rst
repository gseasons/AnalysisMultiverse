Node: reg (utility)
===================


 Hierarchy : wf.reg
 Exec ID : reg


Original Inputs
---------------


* T1w : /Volumes/NewVolume/a_copy/sub-031S0618_T1w.nii.gz
* applywarp_apply_isoxfm : 
* applywarp_apply_xfm : True
* applywarp_interp : spline
* applywarp_no_resample : True
* applywarp_uses_qform : 
* bbr : True
* bet : /Volumes/NewVolume/a_copy/sub-031S0618_T1w_corrected_xform_masked.nii.gz
* concatenate : True
* corrected_img : /Volumes/NewVolume/a_copy/sub-031S0618_task-rest_bold_roi_mcf_st.nii.gz
* function_str : def registration(T1w, mask, start_img, corrected_img, bet, wm_file, bbr, regpre_interp, regpre_no_resample, applywarp_apply_isoxfm, applywarp_apply_xfm, applywarp_uses_qform, concatenate, regbbr_interp, regbbr_no_resample, applywarp_interp, applywarp_no_resample, plugin_args):
    from nipype.interfaces.fsl import FLIRT
    from nipype.interfaces.base import Undefined
    from nipype import IdentityInterface
    from nipype import Workflow, Node
    from os.path import join as opj
    import os, glob, re

    reg = Workflow(name='reg')
    reg.base_dir = os.getcwd()

    bbr = vars().get('bbr', True)
    concatenate = vars().get('concatenate', True)
    #CHECK TO VERIFY THAT DOF IS CALLED DOF
    regpre = Node(FLIRT(in_file=start_img, dof=6, reference=bet, output_type='NIFTI_GZ'), name='regpre')

    applywarp = Node(FLIRT(in_file=corrected_img, reference=bet, output_type='NIFTI_GZ'), name='applywarp')

    outnode =  Node(IdentityInterface(fields=['warped','out_mat', 'mask']), name='outnode')

    if bbr:
        regbbr = Node(FLIRT(cost='bbr', reference=T1w, in_file=start_img, output_type='NIFTI_GZ',
                        schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch')), name='regbbr')
        regbbr.inputs.wm_seg = wm_file

        reg.connect([(regpre, regbbr, [('out_matrix_file', 'in_matrix_file')]),
                     (regbbr, outnode, [('out_matrix_file', 'out_mat')]),
                     ])

        if not concatenate:
            reg.connect([(regbbr, applywarp, [('out_matrix_file', 'in_matrix_file')])])

        node_reg = 'regbbr'
    else:
        reg.connect([(regpre, outnode, [('out_matrix_file', 'out_mat')]),
                     ])

        if not concatenate:
            reg.connect([(regpre, applywarp, [('out_matrix_file', 'in_matrix_file')])])

        node_reg = 'regpre'

    node_warp = 'applywarp'

    for param in ['regpre_interp', 'regpre_no_resample', 'applywarp_apply_isoxfm', 'applywarp_apply_xfm', 'applywarp_uses_qform', 'regbbr_interp', 'regbbr_no_resample', 'applywarp_interp', 'applywarp_no_resample']:
        search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)
        if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])
        else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)

    reg.run(plugin='IPython', plugin_args=plugin_args)

    out_mat = glob.glob(os.getcwd() + '/reg/' + node_reg + '/*.mat')[0]
    if concatenate:
        warped = corrected_img
    else:
        warped = glob.glob(os.getcwd() + '/reg/' + node_warp + '/*.nii*')[0]

    return out_mat, warped, glob.glob(os.getcwd() + '/reg/**', recursive=True)

* mask : /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
* plugin_args : {'profile': '', 'cluster_id': '1647026890-u1vi'}
* regbbr_interp : spline
* regbbr_no_resample : True
* regpre_interp : spline
* regpre_no_resample : True
* start_img : /Volumes/NewVolume/a_copy/sub-031S0618_task-rest_bold_roi_mcf_roi.nii.gz
* wm_file : /Volumes/NewVolume/a_copy/sub-031S0618_T1w_corrected_xform_masked_seg_2.nii.gz


Execution Inputs
----------------


* T1w : /Volumes/NewVolume/a_copy/sub-031S0618_T1w.nii.gz
* applywarp_apply_isoxfm : 
* applywarp_apply_xfm : True
* applywarp_interp : spline
* applywarp_no_resample : True
* applywarp_uses_qform : 
* bbr : True
* bet : /Volumes/NewVolume/a_copy/sub-031S0618_T1w_corrected_xform_masked.nii.gz
* concatenate : True
* corrected_img : /Volumes/NewVolume/a_copy/sub-031S0618_task-rest_bold_roi_mcf_st.nii.gz
* function_str : def registration(T1w, mask, start_img, corrected_img, bet, wm_file, bbr, regpre_interp, regpre_no_resample, applywarp_apply_isoxfm, applywarp_apply_xfm, applywarp_uses_qform, concatenate, regbbr_interp, regbbr_no_resample, applywarp_interp, applywarp_no_resample, plugin_args):
    from nipype.interfaces.fsl import FLIRT
    from nipype.interfaces.base import Undefined
    from nipype import IdentityInterface
    from nipype import Workflow, Node
    from os.path import join as opj
    import os, glob, re

    reg = Workflow(name='reg')
    reg.base_dir = os.getcwd()

    bbr = vars().get('bbr', True)
    concatenate = vars().get('concatenate', True)
    #CHECK TO VERIFY THAT DOF IS CALLED DOF
    regpre = Node(FLIRT(in_file=start_img, dof=6, reference=bet, output_type='NIFTI_GZ'), name='regpre')

    applywarp = Node(FLIRT(in_file=corrected_img, reference=bet, output_type='NIFTI_GZ'), name='applywarp')

    outnode =  Node(IdentityInterface(fields=['warped','out_mat', 'mask']), name='outnode')

    if bbr:
        regbbr = Node(FLIRT(cost='bbr', reference=T1w, in_file=start_img, output_type='NIFTI_GZ',
                        schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch')), name='regbbr')
        regbbr.inputs.wm_seg = wm_file

        reg.connect([(regpre, regbbr, [('out_matrix_file', 'in_matrix_file')]),
                     (regbbr, outnode, [('out_matrix_file', 'out_mat')]),
                     ])

        if not concatenate:
            reg.connect([(regbbr, applywarp, [('out_matrix_file', 'in_matrix_file')])])

        node_reg = 'regbbr'
    else:
        reg.connect([(regpre, outnode, [('out_matrix_file', 'out_mat')]),
                     ])

        if not concatenate:
            reg.connect([(regpre, applywarp, [('out_matrix_file', 'in_matrix_file')])])

        node_reg = 'regpre'

    node_warp = 'applywarp'

    for param in ['regpre_interp', 'regpre_no_resample', 'applywarp_apply_isoxfm', 'applywarp_apply_xfm', 'applywarp_uses_qform', 'regbbr_interp', 'regbbr_no_resample', 'applywarp_interp', 'applywarp_no_resample']:
        search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)
        if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])
        else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)

    reg.run(plugin='IPython', plugin_args=plugin_args)

    out_mat = glob.glob(os.getcwd() + '/reg/' + node_reg + '/*.mat')[0]
    if concatenate:
        warped = corrected_img
    else:
        warped = glob.glob(os.getcwd() + '/reg/' + node_warp + '/*.nii*')[0]

    return out_mat, warped, glob.glob(os.getcwd() + '/reg/**', recursive=True)

* mask : /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
* plugin_args : {'profile': '', 'cluster_id': '1647026890-u1vi'}
* regbbr_interp : spline
* regbbr_no_resample : True
* regpre_interp : spline
* regpre_no_resample : True
* start_img : /Volumes/NewVolume/a_copy/sub-031S0618_task-rest_bold_roi_mcf_roi.nii.gz
* wm_file : /Volumes/NewVolume/a_copy/sub-031S0618_T1w_corrected_xform_masked_seg_2.nii.gz


Execution Outputs
-----------------


* egg : ['/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/index.html', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/_report', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/_report/report.rst', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/_0xe03071bffdaec6c7dd8240cfd3507cc3.json', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/_inputs.pklz', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/result_regbbr.pklz', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/sub-031S0618_task-rest_bold_roi_mcf_roi_flirt.mat', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/_node.pklz', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/command.txt', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/_report', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/_report/report.rst', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/_inputs.pklz', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/_0x4b4ae41faeeb189fa5dfe05a7eeb4751.json', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/result_regpre.pklz', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/sub-031S0618_task-rest_bold_roi_mcf_roi_flirt.mat', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/_node.pklz', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regpre/command.txt', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/d3.js', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/graph1.json', '/Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/graph.json']
* out_mat : /Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg/reg/regbbr/sub-031S0618_task-rest_bold_roi_mcf_roi_flirt.mat
* warped : /Volumes/NewVolume/a_copy/sub-031S0618_task-rest_bold_roi_mcf_st.nii.gz


Runtime info
------------


* duration : 1131.514474
* hostname : Grahams-MBP.hitronhub.home
* prev_wd : /Users/grahamseasons/fMRI/analysis_multiverse/hide/old
* working_dir : /Users/grahamseasons/fMRI/analysis_multiverse/hide/old/wf/reg


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
* IPP_CLUSTER_ID : 1647026890-u1vi
* IPP_CONNECTION_INFO : {"ssh": "", "interface": "tcp://127.0.0.1", "registration": 49958, "control": 49959, "mux": 49960, "task": 49961, "iopub": 49962, "hb_ping": 49963, "hb_pong": 49964, "broadcast": [49965, 49966], "key": "9b0d3244-47713eecf78a8f7b8267f7fe", "curve_serverkey": null, "location": "Grahams-MBP.hitronhub.home", "pack": "json", "unpack": "json", "signature_scheme": "hmac-sha256"}
* IPP_PROFILE_DIR : /Users/grahamseasons/.ipython/profile_default
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

