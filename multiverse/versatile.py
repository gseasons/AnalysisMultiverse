#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 15:19:49 2021

@author: grahamseasons
"""

"""
MODIFIED NIPYPE SpecifyModel TO ACCOMODATE RESTING STATE fMRI
"""

"""
The modelgen module provides classes for specifying designs for individual
subject analysis of task-based fMRI experiments. In particular it also includes
algorithms for generating regressors for sparse and sparse-clustered acquisition
experiments.
"""
import csv, math, os

from nibabel import load
import numpy as np

from nipype.interfaces.base import (
    BaseInterface,
    TraitedSpec,
    InputMultiPath,
    traits,
    File,
    Bunch,
    BaseInterfaceInputSpec,
    isdefined,
)
from nipype.utils.filemanip import ensure_list
from nipype.utils.misc import normalize_mc_params
from nipype import logging

iflogger = logging.getLogger("nipype.interface")


def spm_hrf(RT, P=None, fMRI_T=16):
    """
    python implementation of spm_hrf

    See ``spm_hrf`` for implementation details::

      % RT   - scan repeat time
      % p    - parameters of the response function (two gamma
      % functions)
      % defaults  (seconds)
      % p(0) - delay of response (relative to onset)       6
      % p(1) - delay of undershoot (relative to onset)    16
      % p(2) - dispersion of response                      1
      % p(3) - dispersion of undershoot                    1
      % p(4) - ratio of response to undershoot             6
      % p(5) - onset (seconds)                             0
      % p(6) - length of kernel (seconds)                 32
      %
      % hrf  - hemodynamic response function
      % p    - parameters of the response function

    The following code using ``scipy.stats.distributions.gamma``
    doesn't return the same result as the ``spm_Gpdf`` function::

        hrf = gamma.pdf(u, p[0]/p[2], scale=dt/p[2]) -
              gamma.pdf(u, p[1]/p[3], scale=dt/p[3])/p[4]

    Example
    -------
    >>> print(spm_hrf(2))
    [  0.00000000e+00   8.65660810e-02   3.74888236e-01   3.84923382e-01
       2.16117316e-01   7.68695653e-02   1.62017720e-03  -3.06078117e-02
      -3.73060781e-02  -3.08373716e-02  -2.05161334e-02  -1.16441637e-02
      -5.82063147e-03  -2.61854250e-03  -1.07732374e-03  -4.10443522e-04
      -1.46257507e-04]

    """
    from scipy.special import gammaln

    p = np.array([6, 16, 1, 1, 6, 0, 32], dtype=float)
    if P is not None:
        p[0 : len(P)] = P

    _spm_Gpdf = lambda x, h, l: np.exp(
        h * np.log(l) + (h - 1) * np.log(x) - (l * x) - gammaln(h)
    )
    # modelled hemodynamic response function - {mixture of Gammas}
    dt = RT / float(fMRI_T)
    u = np.arange(0, int(p[6] / dt + 1)) - p[5] / dt
    with np.errstate(divide="ignore"):  # Known division-by-zero
        hrf = (
            _spm_Gpdf(u, p[0] / p[2], dt / p[2])
            - _spm_Gpdf(u, p[1] / p[3], dt / p[3]) / p[4]
        )
    idx = np.arange(0, int((p[6] / RT) + 1)) * fMRI_T
    hrf = hrf[idx]
    hrf = hrf / np.sum(hrf)
    return hrf


def orth(x_in, y_in):
    """Orthogonalize y_in with respect to x_in.

    >>> orth_expected = np.array([1.7142857142857144, 0.42857142857142883, \
                                  -0.85714285714285676])
    >>> err = np.abs(np.array(orth([1, 2, 3],[4, 5, 6]) - orth_expected))
    >>> all(err < np.finfo(float).eps)
    True

    """
    x = np.array(x_in)[:, None]
    y = np.array(y_in)[:, None]
    y = y - np.dot(x, np.dot(np.linalg.inv(np.dot(x.T, x)), np.dot(x.T, y)))
    if np.linalg.norm(y, 1) > np.exp(-32):
        y = y[:, 0].tolist()
    else:
        y = y_in
    return y


def scale_timings(timelist, input_units, output_units, time_repetition):
    """
    Scale timings given input and output units (scans/secs).

    Parameters
    ----------
    timelist: list of times to scale
    input_units: 'secs' or 'scans'
    output_units: Ibid.
    time_repetition: float in seconds

    """
    if input_units == output_units:
        _scalefactor = 1.0

    if (input_units == "scans") and (output_units == "secs"):
        _scalefactor = time_repetition

    if (input_units == "secs") and (output_units == "scans"):
        _scalefactor = 1.0 / time_repetition
    timelist = [np.max([0.0, _scalefactor * t]) for t in timelist]
    return timelist


def bids_gen_info(
    bids_event_files, condition_column="", amplitude_column=None, time_repetition=False
):
    """
    Generate a subject_info structure from a list of BIDS .tsv event files.

    Parameters
    ----------
    bids_event_files : list of str
        Filenames of BIDS .tsv event files containing columns including:
        'onset', 'duration', and 'trial_type' or the `condition_column` value.
    condition_column : str
        Column of files in `bids_event_files` based on the values of which
        events will be sorted into different regressors
    amplitude_column : str
        Column of files in `bids_event_files` based on the values of which
        to apply amplitudes to events. If unspecified, all events will be
        represented with an amplitude of 1.

    Returns
    -------
    subject_info: list of Bunch

    """
    info = []
    for bids_event_file in bids_event_files:
        with open(bids_event_file) as f:
            f_events = csv.DictReader(f, skipinitialspace=True, delimiter="\t")
            events = [{k: v for k, v in row.items()} for row in f_events]
        if not condition_column:
            condition_column = "_trial_type"
            for i in events:
                i.update({condition_column: "ev0"})
        conditions = sorted(set([i[condition_column] for i in events]))
        runinfo = Bunch(conditions=[], onsets=[], durations=[], amplitudes=[])
        for condition in conditions:
            selected_events = [i for i in events if i[condition_column] == condition]
            onsets = [float(i["onset"]) for i in selected_events]
            durations = [float(i["duration"]) for i in selected_events]
            if time_repetition:
                decimals = math.ceil(-math.log10(time_repetition))
                onsets = [np.round(i, decimals) for i in onsets]
                durations = [np.round(i, decimals) for i in durations]
            runinfo.conditions.append(condition)
            runinfo.onsets.append(onsets)
            runinfo.durations.append(durations)
            try:
                amplitudes = [float(i[amplitude_column]) for i in selected_events]
                runinfo.amplitudes.append(amplitudes)
            except KeyError:
                runinfo.amplitudes.append([1] * len(onsets))
        info.append(runinfo)
    return info


def gen_info(run_event_files):
    """Generate subject_info structure from a list of event files."""
    info = []
    for i, event_files in enumerate(run_event_files):
        runinfo = Bunch(conditions=[], onsets=[], durations=[], amplitudes=[])
        for event_file in event_files:
            _, name = os.path.split(event_file)
            if ".run" in name:
                name, _ = name.split(".run%03d" % (i + 1))
            elif ".txt" in name:
                name, _ = name.split(".txt")
                
            if "rest" in name:
                runinfo = Bunch(conditions=[], val=[])
                runinfo.conditions.append(name)
                event_info = np.atleast_2d(np.loadtxt(event_file))
                runinfo.val.append(event_info[0, :].tolist())
                continue

            runinfo.conditions.append(name)
            event_info = np.atleast_2d(np.loadtxt(event_file))
            runinfo.onsets.append(event_info[:, 0].tolist())
            if event_info.shape[1] > 1:
                runinfo.durations.append(event_info[:, 1].tolist())
            else:
                runinfo.durations.append([0])

            if event_info.shape[1] > 2:
                runinfo.amplitudes.append(event_info[:, 2].tolist())
            else:
                delattr(runinfo, "amplitudes")
        info.append(runinfo)
    return info


class SpecifyModelInputSpec(BaseInterfaceInputSpec):
    subject_info = InputMultiPath(
        Bunch,
        mandatory=True,
        xor=["subject_info", "event_files", "bids_event_file"],
        desc="Bunch or List(Bunch) subject-specific "
        "condition information. see "
        ":ref:`nipype.algorithms.modelgen.SpecifyModel` or for details",
    )
    event_files = InputMultiPath(
        traits.List(File(exists=True)),
        mandatory=True,
        xor=["subject_info", "event_files", "bids_event_file"],
        desc="List of event description files 1, 2 or 3 "
        "column format corresponding to onsets, "
        "durations and amplitudes",
    )
    bids_event_file = InputMultiPath(
        File(exists=True),
        mandatory=True,
        xor=["subject_info", "event_files", "bids_event_file"],
        desc="TSV event file containing common BIDS fields: `onset`,"
        "`duration`, and categorization and amplitude columns",
    )
    bids_condition_column = traits.Str(
        default_value="trial_type",
        usedefault=True,
        desc="Column of the file passed to ``bids_event_file`` to the "
        "unique values of which events will be assigned"
        "to regressors",
    )
    bids_amplitude_column = traits.Str(
        desc="Column of the file passed to ``bids_event_file`` "
        "according to which to assign amplitudes to events"
    )
    realignment_parameters = InputMultiPath(
        File(exists=True),
        desc="Realignment parameters returned by motion correction algorithm",
        copyfile=False,
    )
    parameter_source = traits.Enum(
        "SPM",
        "FSL",
        "AFNI",
        "FSFAST",
        "NIPY",
        usedefault=True,
        desc="Source of motion parameters",
    )
    outlier_files = InputMultiPath(
        File(exists=True),
        desc="Files containing scan outlier indices that should be tossed",
        copyfile=False,
    )
    functional_runs = InputMultiPath(
        traits.Either(traits.List(File(exists=True)), File(exists=True)),
        mandatory=True,
        desc="Data files for model. List of 4D "
        "files or list of list of 3D "
        "files per session",
        copyfile=False,
    )
    input_units = traits.Enum(
        "secs",
        "scans",
        mandatory=True,
        desc="Units of event onsets and durations (secs "
        "or scans). Output units are always in secs",
    )
    high_pass_filter_cutoff = traits.Float(
        mandatory=True, desc="High-pass filter cutoff in secs"
    )
    time_repetition = traits.Float(
        mandatory=True,
        desc="Time between the start of one volume "
        "to the start of  the next image volume.",
    )
    # Not implemented yet
    # polynomial_order = traits.Range(0, low=0,
    #        desc ='Number of polynomial functions to model high pass filter.')


class SpecifyModelOutputSpec(TraitedSpec):
    session_info = traits.Any(desc="Session info for level1designs")


class SpecifyModelVersatile(BaseInterface):
    """
    Makes a model specification compatible with spm/fsl designers.

    The subject_info field should contain paradigm information in the form of
    a Bunch or a list of Bunch. The Bunch should contain the following
    information::

        [Mandatory]
        conditions : list of names
        onsets : lists of onsets corresponding to each condition
        durations : lists of durations corresponding to each condition. Should be
            left to a single 0 if all events are being modelled as impulses.

        [Optional]
        regressor_names : list of str
            list of names corresponding to each column. Should be None if
            automatically assigned.
        regressors : list of lists
            values for each regressor - must correspond to the number of
            volumes in the functional run
        amplitudes : lists of amplitudes for each event. This will be ignored by
            SPM's Level1Design.

        The following two (tmod, pmod) will be ignored by any Level1Design class
        other than SPM:

        tmod : lists of conditions that should be temporally modulated. Should
            default to None if not being used.
        pmod : list of Bunch corresponding to conditions
          - name : name of parametric modulator
          - param : values of the modulator
          - poly : degree of modulation

    Alternatively, you can provide information through event files.

    The event files have to be in 1, 2 or 3 column format with the columns
    corresponding to Onsets, Durations and Amplitudes and they have to have the
    name event_name.runXXX... e.g.: Words.run001.txt. The event_name part will
    be used to create the condition names.

    Examples
    --------
    >>> from nipype.algorithms import modelgen
    >>> from nipype.interfaces.base import Bunch
    >>> s = modelgen.SpecifyModel()
    >>> s.inputs.input_units = 'secs'
    >>> s.inputs.functional_runs = ['functional2.nii', 'functional3.nii']
    >>> s.inputs.time_repetition = 6
    >>> s.inputs.high_pass_filter_cutoff = 128.
    >>> evs_run2 = Bunch(conditions=['cond1'], onsets=[[2, 50, 100, 180]], durations=[[1]])
    >>> evs_run3 = Bunch(conditions=['cond1'], onsets=[[30, 40, 100, 150]], durations=[[1]])
    >>> s.inputs.subject_info = [evs_run2, evs_run3]

    >>> # Using pmod
    >>> evs_run2 = Bunch(conditions=['cond1', 'cond2'], onsets=[[2, 50], [100, 180]], \
durations=[[0], [0]], pmod=[Bunch(name=['amp'], poly=[2], param=[[1, 2]]), \
None])
    >>> evs_run3 = Bunch(conditions=['cond1', 'cond2'], onsets=[[20, 120], [80, 160]], \
durations=[[0], [0]], pmod=[Bunch(name=['amp'], poly=[2], param=[[1, 2]]), \
None])
    >>> s.inputs.subject_info = [evs_run2, evs_run3]

    """

    input_spec = SpecifyModelInputSpec
    output_spec = SpecifyModelOutputSpec

    def _generate_standard_design(
        self, infolist, functional_runs=None, realignment_parameters=None, outliers=None
    ):
        """Generate a standard design matrix paradigm given information about each run."""
        sessinfo = []
        output_units = "secs"
        if "output_units" in self.inputs.traits():
            output_units = self.inputs.output_units

        for i, info in enumerate(infolist):
            sessinfo.insert(i, dict(cond=[]))
            if isdefined(self.inputs.high_pass_filter_cutoff):
                sessinfo[i]["hpf"] = float(self.inputs.high_pass_filter_cutoff)

            if hasattr(info, "conditions") and info.conditions is not None:
                for cid, cond in enumerate(info.conditions):
                    sessinfo[i]["cond"].insert(cid, dict())
                    sessinfo[i]["cond"][cid]["name"] = info.conditions[cid]
                    
                    if "rest" in info.conditions[cid]:
                        sessinfo[i]["cond"][cid]["val"] = info.val[cid]
                        continue
                    
                    scaled_onset = scale_timings(
                        info.onsets[cid],
                        self.inputs.input_units,
                        output_units,
                        self.inputs.time_repetition,
                    )
                    sessinfo[i]["cond"][cid]["onset"] = scaled_onset
                    scaled_duration = scale_timings(
                        info.durations[cid],
                        self.inputs.input_units,
                        output_units,
                        self.inputs.time_repetition,
                    )
                    sessinfo[i]["cond"][cid]["duration"] = scaled_duration
                    if hasattr(info, "amplitudes") and info.amplitudes:
                        sessinfo[i]["cond"][cid]["amplitudes"] = info.amplitudes[cid]

                    if hasattr(info, "tmod") and info.tmod and len(info.tmod) > cid:
                        sessinfo[i]["cond"][cid]["tmod"] = info.tmod[cid]

                    if hasattr(info, "pmod") and info.pmod and len(info.pmod) > cid:
                        if info.pmod[cid]:
                            sessinfo[i]["cond"][cid]["pmod"] = []
                            for j, name in enumerate(info.pmod[cid].name):
                                sessinfo[i]["cond"][cid]["pmod"].insert(j, {})
                                sessinfo[i]["cond"][cid]["pmod"][j]["name"] = name
                                sessinfo[i]["cond"][cid]["pmod"][j]["poly"] = info.pmod[
                                    cid
                                ].poly[j]
                                sessinfo[i]["cond"][cid]["pmod"][j][
                                    "param"
                                ] = info.pmod[cid].param[j]

            sessinfo[i]["regress"] = []
            if hasattr(info, "regressors") and info.regressors is not None:
                for j, r in enumerate(info.regressors):
                    sessinfo[i]["regress"].insert(j, dict(name="", val=[]))
                    if (
                        hasattr(info, "regressor_names")
                        and info.regressor_names is not None
                    ):
                        sessinfo[i]["regress"][j]["name"] = info.regressor_names[j]
                    else:
                        sessinfo[i]["regress"][j]["name"] = "UR%d" % (j + 1)
                    sessinfo[i]["regress"][j]["val"] = info.regressors[j]
            sessinfo[i]["scans"] = functional_runs[i]

        if realignment_parameters is not None:
            for i, rp in enumerate(realignment_parameters):
                mc = realignment_parameters[i]
                for col in range(mc.shape[1]):
                    colidx = len(sessinfo[i]["regress"])
                    sessinfo[i]["regress"].insert(colidx, dict(name="", val=[]))
                    sessinfo[i]["regress"][colidx]["name"] = "Realign%d" % (col + 1)
                    sessinfo[i]["regress"][colidx]["val"] = mc[:, col].tolist()

        if outliers is not None:
            for i, out in enumerate(outliers):
                numscans = 0
                for f in ensure_list(sessinfo[i]["scans"]):
                    shape = load(f).shape
                    if len(shape) == 3 or shape[3] == 1:
                        iflogger.warning(
                            "You are using 3D instead of 4D "
                            "files. Are you sure this was "
                            "intended?"
                        )
                        numscans += 1
                    else:
                        numscans += shape[3]

                for j, scanno in enumerate(out):
                    colidx = len(sessinfo[i]["regress"])
                    sessinfo[i]["regress"].insert(colidx, dict(name="", val=[]))
                    sessinfo[i]["regress"][colidx]["name"] = "Outlier%d" % (j + 1)
                    sessinfo[i]["regress"][colidx]["val"] = np.zeros((1, numscans))[
                        0
                    ].tolist()
                    sessinfo[i]["regress"][colidx]["val"][int(scanno)] = 1
        return sessinfo

    def _generate_design(self, infolist=None):
        """Generate design specification for a typical fmri paradigm"""
        realignment_parameters = []
        if isdefined(self.inputs.realignment_parameters):
            for parfile in self.inputs.realignment_parameters:
                realignment_parameters.append(
                    np.apply_along_axis(
                        func1d=normalize_mc_params,
                        axis=1,
                        arr=np.loadtxt(parfile),
                        source=self.inputs.parameter_source,
                    )
                )
        outliers = []
        if isdefined(self.inputs.outlier_files):
            for filename in self.inputs.outlier_files:
                try:
                    outindices = np.loadtxt(filename, dtype=int)
                except IOError:
                    outliers.append([])
                else:
                    if outindices.size == 1:
                        outliers.append([outindices.tolist()])
                    else:
                        outliers.append(outindices.tolist())

        if infolist is None:
            if isdefined(self.inputs.subject_info):
                infolist = self.inputs.subject_info
            elif isdefined(self.inputs.event_files):
                infolist = gen_info(self.inputs.event_files)
            elif isdefined(self.inputs.bids_event_file):
                infolist = bids_gen_info(
                    self.inputs.bids_event_file,
                    self.inputs.bids_condition_column,
                    self.inputs.bids_amplitude_column,
                    self.inputs.time_repetition,
                )
        self._sessinfo = self._generate_standard_design(
            infolist,
            functional_runs=self.inputs.functional_runs,
            realignment_parameters=realignment_parameters,
            outliers=outliers,
        )

    def _run_interface(self, runtime):
        """ """
        self._sessioninfo = None
        self._generate_design()
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        if not hasattr(self, "_sessinfo"):
            self._generate_design()
        outputs["session_info"] = self._sessinfo

        return outputs
    

from string import Template

from nipype.interfaces.base import (
    File,
    traits,
    TraitedSpec,
    BaseInterface,
    InputMultiPath,
    OutputMultiPath,
    BaseInterfaceInputSpec,
)


class Level1DesignInputSpec(BaseInterfaceInputSpec):
    interscan_interval = traits.Float(
        mandatory=True, desc="Interscan  interval (in secs)"
    )
    session_info = traits.Any(
        mandatory=True,
        desc=("Session specific information generated " "by ``modelgen.SpecifyModel``"),
    )
    bases = traits.Either(
        traits.Dict(
            traits.Enum("dgamma"), traits.Dict(traits.Enum("derivs"), traits.Bool)
        ),
        traits.Dict(
            traits.Enum("gamma"),
            traits.Dict(traits.Enum("derivs", "gammasigma", "gammadelay")),
        ),
        traits.Dict(
            traits.Enum("custom"), traits.Dict(traits.Enum("bfcustompath"), traits.Str)
        ),
        traits.Dict(traits.Enum("none"), traits.Dict()),
        traits.Dict(traits.Enum("none"), traits.Enum(None)),
        mandatory=True,
        desc=(
            "name of basis function and options e.g., " "{'dgamma': {'derivs': True}}"
        ),
    )
    orthogonalization = traits.Dict(
        traits.Int,
        traits.Dict(traits.Int, traits.Either(traits.Bool, traits.Int)),
        desc=(
            "which regressors to make orthogonal e.g., "
            "{1: {0:0,1:0,2:0}, 2: {0:1,1:1,2:0}} to make the second "
            "regressor in a 2-regressor model orthogonal to the first."
        ),
        usedefault=True,
    )
    model_serial_correlations = traits.Bool(
        desc="Option to model serial correlations using an \
autoregressive estimator (order 1). Setting this option is only \
useful in the context of the fsf file. If you set this to False, you need to \
repeat this option for FILMGLS by setting autocorr_noestimate to True",
        mandatory=True,
    )
    contrasts = traits.List(
        traits.Either(
            traits.Tuple(
                traits.Str,
                traits.Enum("T"),
                traits.List(traits.Str),
                traits.List(traits.Float),
            ),
            traits.Tuple(
                traits.Str,
                traits.Enum("T"),
                traits.List(traits.Str),
                traits.List(traits.Float),
                traits.List(traits.Float),
            ),
            traits.Tuple(
                traits.Str,
                traits.Enum("F"),
                traits.List(
                    traits.Either(
                        traits.Tuple(
                            traits.Str,
                            traits.Enum("T"),
                            traits.List(traits.Str),
                            traits.List(traits.Float),
                        ),
                        traits.Tuple(
                            traits.Str,
                            traits.Enum("T"),
                            traits.List(traits.Str),
                            traits.List(traits.Float),
                            traits.List(traits.Float),
                        ),
                    )
                ),
            ),
        ),
        desc="List of contrasts with each contrast being a list of the form - \
[('name', 'stat', [condition list], [weight list], [session list])]. if \
session list is None or not provided, all sessions are used. For F \
contrasts, the condition list should contain previously defined \
T-contrasts.",
    )


class Level1DesignOutputSpec(TraitedSpec):
    fsf_files = OutputMultiPath(File(exists=True), desc="FSL feat specification files")
    ev_files = OutputMultiPath(
        traits.List(File(exists=True)), desc="condition information files"
    )


class Level1DesignVersatile(BaseInterface):
    """Generate FEAT specific files

    Examples
    --------

    >>> level1design = Level1Design()
    >>> level1design.inputs.interscan_interval = 2.5
    >>> level1design.inputs.bases = {'dgamma':{'derivs': False}}
    >>> level1design.inputs.session_info = 'session_info.npz'
    >>> level1design.run() # doctest: +SKIP

    """

    input_spec = Level1DesignInputSpec
    output_spec = Level1DesignOutputSpec

    def _create_ev_file(self, evfname, evinfo):
        f = open(evfname, "wt")
        for i in evinfo:
            if len(i) == 3:
                f.write("%f %f %f\n" % (i[0], i[1], i[2]))
            else:
                f.write("%f\n" % i[0])
        f.close()

    def _create_ev_files(
        self,
        cwd,
        runinfo,
        runidx,
        ev_parameters,
        orthogonalization,
        contrasts,
        do_tempfilter,
        basis_key,
    ):
        """Creates EV files from condition and regressor information.

        Parameters:
        -----------

        runinfo : dict
            Generated by `SpecifyModel` and contains information
            about events and other regressors.
        runidx  : int
            Index to run number
        ev_parameters : dict
            A dictionary containing the model parameters for the
            given design type.
        orthogonalization : dict
            A dictionary of dictionaries specifying orthogonal EVs.
        contrasts : list of lists
            Information on contrasts to be evaluated
        """
        conds = {}
        evname = []
        if basis_key == "dgamma":
            basis_key = "hrf"
        elif basis_key == "gamma":
            try:
                _ = ev_parameters["gammasigma"]
            except KeyError:
                ev_parameters["gammasigma"] = 3
            try:
                _ = ev_parameters["gammadelay"]
            except KeyError:
                ev_parameters["gammadelay"] = 6
        ev_template = load_template("feat_ev_" + basis_key + ".tcl")
        ev_none = load_template("feat_ev_none.tcl")
        ev_ortho = load_template("feat_ev_ortho.tcl")
        ev_txt = ""
        # generate sections for conditions and other nuisance
        # regressors
        num_evs = [0, 0]
        for field in ["cond", "regress"]:
            for i, cond in enumerate(runinfo[field]):
                name = cond["name"]
                evname.append(name)
                evfname = os.path.join(
                    cwd, "ev_%s_%d_%d.txt" % (name, runidx, len(evname))
                )
                evinfo = []
                num_evs[0] += 1
                num_evs[1] += 1
                if field == "cond":
                    try:
                        for j, onset in enumerate(cond["onset"]):
                            try:
                                amplitudes = cond["amplitudes"]
                                if len(amplitudes) > 1:
                                    amp = amplitudes[j]
                                else:
                                    amp = amplitudes[0]
                            except KeyError:
                                amp = 1
                            if len(cond["duration"]) > 1:
                                evinfo.insert(j, [onset, cond["duration"][j], amp])
                            else:
                                evinfo.insert(j, [onset, cond["duration"][0], amp])
                        ev_parameters["cond_file"] = evfname
                        ev_parameters["ev_num"] = num_evs[0]
                        ev_parameters["ev_name"] = name
                        ev_parameters["tempfilt_yn"] = do_tempfilter
                        if "basisorth" not in ev_parameters:
                            ev_parameters["basisorth"] = 1
                        if "basisfnum" not in ev_parameters:
                            ev_parameters["basisfnum"] = 1
                        try:
                            ev_parameters["fsldir"] = os.environ["FSLDIR"]
                        except KeyError:
                            if basis_key == "flobs":
                                raise Exception("FSL environment variables not set")
                            else:
                                ev_parameters["fsldir"] = "/usr/share/fsl"
                        ev_parameters["temporalderiv"] = int(
                            bool(ev_parameters.get("derivs", False))
                        )
                        if ev_parameters["temporalderiv"]:
                            evname.append(name + "TD")
                            num_evs[1] += 1
                        ev_txt += ev_template.substitute(ev_parameters)
                    except:
                        evinfo = [[j] for j in cond["val"]]
                        ev_parameters["cond_file"] = evfname
                        ev_parameters["ev_num"] = num_evs[0]
                        ev_parameters["ev_name"] = name
                        ev_parameters["tempfilt_yn"] = do_tempfilter
                        if "basisorth" not in ev_parameters:
                            ev_parameters["basisorth"] = 1
                        if "basisfnum" not in ev_parameters:
                            ev_parameters["basisfnum"] = 1
                        try:
                            ev_parameters["fsldir"] = os.environ["FSLDIR"]
                        except KeyError:
                            if basis_key == "flobs":
                                raise Exception("FSL environment variables not set")
                            else:
                                ev_parameters["fsldir"] = "/usr/share/fsl"
                        ev_parameters["temporalderiv"] = int(
                            bool(ev_parameters.get("derivs", False))
                        )
                        #Resting state doesn't use temporal derivative
                        ev_parameters["temporalderiv"] = 0
                        #Don't convolve if resting state
                        ev_txt += ev_none.substitute(ev_parameters)
                        if ev_parameters["temporalderiv"]:
                            evname.append(name + "TD")
                            num_evs[1] += 1
                            ev_txt.replace("set fmri(deriv_yn1) 0", "set fmri(deriv_yn1) 1")
                        #ev_txt += ev_template.substitute(ev_parameters)
# =============================================================================
#                         ev_txt += ev_none.substitute(
#                         ev_num=num_evs[0],
#                         ev_name=name,
#                         tempfilt_yn=do_tempfilter,
#                         cond_file=evfname,
#                         )
# =============================================================================
                    
                elif field == "regress":
                    evinfo = [[j] for j in cond["val"]]
                    ev_txt += ev_none.substitute(
                        ev_num=num_evs[0],
                        ev_name=name,
                        tempfilt_yn=do_tempfilter,
                        cond_file=evfname,
                    )
                ev_txt += "\n"
                conds[name] = evfname
                self._create_ev_file(evfname, evinfo)
        # add ev orthogonalization
        for i in range(1, num_evs[0] + 1):
            initial = ev_ortho.substitute(c0=i, c1=0, orthogonal=1)
            for j in range(0, num_evs[0] + 1):
                try:
                    orthogonal = int(orthogonalization[i][j])
                except (KeyError, TypeError, ValueError, IndexError):
                    orthogonal = 0
                if orthogonal == 1 and initial not in ev_txt:
                    ev_txt += initial + "\n"
                ev_txt += ev_ortho.substitute(c0=i, c1=j, orthogonal=orthogonal)
                ev_txt += "\n"
        # add contrast info to fsf file
        if isdefined(contrasts):
            contrast_header = load_template("feat_contrast_header.tcl")
            contrast_prolog = load_template("feat_contrast_prolog.tcl")
            contrast_element = load_template("feat_contrast_element.tcl")
            contrast_ftest_element = load_template("feat_contrast_ftest_element.tcl")
            contrastmask_header = load_template("feat_contrastmask_header.tcl")
            contrastmask_footer = load_template("feat_contrastmask_footer.tcl")
            contrastmask_element = load_template("feat_contrastmask_element.tcl")
            # add t/f contrast info
            ev_txt += contrast_header.substitute()
            con_names = []
            for j, con in enumerate(contrasts):
                con_names.append(con[0])
            con_map = {}
            ftest_idx = []
            ttest_idx = []
            for j, con in enumerate(contrasts):
                if con[1] == "F":
                    ftest_idx.append(j)
                    for c in con[2]:
                        if c[0] not in list(con_map.keys()):
                            con_map[c[0]] = []
                        con_map[c[0]].append(j)
                else:
                    ttest_idx.append(j)

            for ctype in ["real", "orig"]:
                for j, con in enumerate(contrasts):
                    if con[1] == "F":
                        continue
                    tidx = ttest_idx.index(j) + 1
                    ev_txt += contrast_prolog.substitute(
                        cnum=tidx, ctype=ctype, cname=con[0]
                    )
                    count = 0
                    for c in range(1, len(evname) + 1):
                        if evname[c - 1].endswith("TD") and ctype == "orig":
                            continue
                        count = count + 1
                        if evname[c - 1] in con[2]:
                            val = con[3][con[2].index(evname[c - 1])]
                        else:
                            val = 0.0
                        ev_txt += contrast_element.substitute(
                            cnum=tidx, element=count, ctype=ctype, val=val
                        )
                        ev_txt += "\n"

                    for fconidx in ftest_idx:
                        fval = 0
                        if con[0] in con_map.keys() and fconidx in con_map[con[0]]:
                            fval = 1
                        ev_txt += contrast_ftest_element.substitute(
                            cnum=ftest_idx.index(fconidx) + 1,
                            element=tidx,
                            ctype=ctype,
                            val=fval,
                        )
                        ev_txt += "\n"

            # add contrast mask info
            ev_txt += contrastmask_header.substitute()
            for j, _ in enumerate(contrasts):
                for k, _ in enumerate(contrasts):
                    if j != k:
                        ev_txt += contrastmask_element.substitute(c1=j + 1, c2=k + 1)
            ev_txt += contrastmask_footer.substitute()
        return num_evs, ev_txt

    def _format_session_info(self, session_info):
        if isinstance(session_info, dict):
            session_info = [session_info]
        return session_info

    def _get_func_files(self, session_info):
        """Returns functional files in the order of runs"""
        func_files = []
        for i, info in enumerate(session_info):
            func_files.insert(i, info["scans"])
        return func_files

    def _run_interface(self, runtime):
        cwd = os.getcwd()
        fsf_header = load_template("feat_header_l1.tcl")
        fsf_postscript = load_template("feat_nongui.tcl")

        prewhiten = 0
        if isdefined(self.inputs.model_serial_correlations):
            prewhiten = int(self.inputs.model_serial_correlations)
        basis_key = list(self.inputs.bases.keys())[0]
        ev_parameters = dict(self.inputs.bases[basis_key])
        session_info = self._format_session_info(self.inputs.session_info)
        func_files = self._get_func_files(session_info)
        n_tcon = 0
        n_fcon = 0
        if isdefined(self.inputs.contrasts):
            for i, c in enumerate(self.inputs.contrasts):
                if c[1] == "T":
                    n_tcon += 1
                elif c[1] == "F":
                    n_fcon += 1

        for i, info in enumerate(session_info):
            do_tempfilter = 1
            if info["hpf"] == np.inf:
                do_tempfilter = 0
            num_evs, cond_txt = self._create_ev_files(
                cwd,
                info,
                i,
                ev_parameters,
                self.inputs.orthogonalization,
                self.inputs.contrasts,
                do_tempfilter,
                basis_key,
            )
            nim = load(func_files[i])
            (_, _, _, timepoints) = nim.shape
            fsf_txt = fsf_header.substitute(
                run_num=i,
                interscan_interval=self.inputs.interscan_interval,
                num_vols=timepoints,
                prewhiten=prewhiten,
                num_evs=num_evs[0],
                num_evs_real=num_evs[1],
                num_tcon=n_tcon,
                num_fcon=n_fcon,
                high_pass_filter_cutoff=info["hpf"],
                temphp_yn=do_tempfilter,
                func_file=func_files[i],
            )
            fsf_txt += cond_txt
            fsf_txt += fsf_postscript.substitute(overwrite=1)

            f = open(os.path.join(cwd, "run%d.fsf" % i), "w")
            f.write(fsf_txt)
            f.close()

        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        cwd = os.getcwd()
        outputs["fsf_files"] = []
        outputs["ev_files"] = []
        basis_key = list(self.inputs.bases.keys())[0]
        ev_parameters = dict(self.inputs.bases[basis_key])
        for runno, runinfo in enumerate(
            self._format_session_info(self.inputs.session_info)
        ):
            outputs["fsf_files"].append(os.path.join(cwd, "run%d.fsf" % runno))
            outputs["ev_files"].insert(runno, [])
            evname = []
            for field in ["cond", "regress"]:
                for i, cond in enumerate(runinfo[field]):
                    name = cond["name"]
                    evname.append(name)
                    evfname = os.path.join(
                        cwd, "ev_%s_%d_%d.txt" % (name, runno, len(evname))
                    )
                    if field == "cond":
                        ev_parameters["temporalderiv"] = int(
                            bool(ev_parameters.get("derivs", False))
                        )
                        if ev_parameters["temporalderiv"]:
                            evname.append(name + "TD")
                    outputs["ev_files"][runno].append(os.path.join(cwd, evfname))
        return outputs
    
def load_template(name):
    """Load a template from the model_templates directory

    Parameters
    ----------
    name : str
        The name of the file to load

    Returns
    -------
    template : string.Template

    """
    from pkg_resources import resource_filename as pkgrf

    full_fname = pkgrf(
        "nipype", os.path.join("interfaces", "fsl", "model_templates", name)
    )
    with open(full_fname) as template_file:
        template = Template(template_file.read())

    return template