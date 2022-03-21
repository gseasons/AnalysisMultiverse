#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 14:58:43 2021

@author: grahamseasons
"""
"""Common graph operations for execution."""
import sys
from copy import deepcopy
from glob import glob
import os
import shutil
from time import sleep, time
from traceback import format_exception

import numpy as np
from collections.abc import Iterable
from collections import Counter
import pickle, re

from ... import logging
from ...utils.misc import str2bool
from ..engine.utils import topological_sort, load_resultfile
from ..engine import MapNode
from .tools import report_crash, report_nodes_not_run, create_pyscript

logger = logging.getLogger("nipype.workflow")


class PluginBase(object):
    """Base class for plugins."""

    def __init__(self, plugin_args=None):
        if plugin_args is None:
            plugin_args = {}
        self.plugin_args = plugin_args
        self._config = None
        self._status_callback = plugin_args.get("status_callback")

    def run(self, graph, config, updatehash=False):
        """
        Instruct the plugin to execute the workflow graph.

        The core plugin member that should be implemented by
        all plugins.

        Parameters
        ----------
        graph :
            a networkx, flattened :abbr:`DAG (Directed Acyclic Graph)`
            to be executed
        config : :obj:`~nipype.config`
            a nipype.config object
        updatehash : :obj:`bool`
            whether cached nodes with stale hash should be just updated.

        """
        raise NotImplementedError


class DistributedPluginBase(PluginBase):
    """
    Execute workflow with a distribution engine

    Combinations of ``proc_done`` and ``proc_pending``:

    +------------+---------------+--------------------------------+
    | proc_done  | proc_pending  | outcome                        |
    +============+===============+================================+
    | True       | False         | Process is finished            |
    +------------+---------------+--------------------------------+
    | True       | True          | Process is currently being run |
    +------------+---------------+--------------------------------+
    | False      | False         | Process is queued              |
    +------------+---------------+--------------------------------+
    | False      | True          | INVALID COMBINATION            |
    +------------+---------------+--------------------------------+

    Attributes
    ----------
    procs : :obj:`list`
        list (N) of underlying interface elements to be processed
    proc_done : :obj:`numpy.ndarray`
        a boolean numpy array (N,) signifying whether a process has been
        submitted for execution
    proc_pending : :obj:`numpy.ndarray`
        a boolean numpy array (N,) signifying whether a
        process is currently running.
    depidx : :obj:`numpy.matrix`
        a boolean matrix (NxN) storing the dependency structure accross
        processes. Process dependencies are derived from each column.

    """

    def __init__(self, plugin_args=None):
        """
        Initialize runtime attributes to none

        """
        super(DistributedPluginBase, self).__init__(plugin_args=plugin_args)
        self.procs = None
        self.depidx = None
        self.refidx = None
        self.mapnodes = None
        self.mapnodesubids = None
        self.proc_done = None
        self.proc_pending = None
        self.pending_tasks = []
        self.max_jobs = self.plugin_args.get("max_jobs", np.inf)

    def _prerun_check(self, graph):
        """Stub method to validate/massage graph and nodes before running"""

    def _postrun_check(self):
        """Stub method to close any open resources"""
    
    def _load_state(self):
        files_missing = True
        
        task = self.__dict__['plugin_args']['task']
        batch = self.__dict__['plugin_args']['batch']
        
        checkpoints = glob('/scratch/processed/reproducibility/checkpoints_' + task + '_batch_' + str(batch) + '/checkpoint_*.pkl')
        if checkpoints:
            checkpoints = sorted(checkpoints, key=lambda val: int(re.search('.*_([0-9]+)', val).group(1)))
            count = 1
            _recent = checkpoints[-count]
            while files_missing:
                try:
                    with open(_recent, 'rb') as f:
                        saved_state = pickle.load(f)
                    files_missing = False
                except Exception as e:
                    try:
                        folder = re.search("a value of '(.*/).*.nii.*' <class", str(e)).group(1)
                        file = re.search("a value of '(.*)' <class", str(e)).group(1)
                        os.makedirs(folder)
                        open(file, 'a').close()
                    except:
                        count += 1
                        _recent = checkpoints[-count]
                        
            if 'iparallel' in self.__dict__:
                indices = range(saved_state['refidx'].shape[0])
                for idx in indices:
                    outdir = saved_state['procs'][idx].output_dir()
                    if saved_state['proc_done'][idx] and saved_state['proc_pending'][idx]:
                        saved_state['proc_done'][idx] = False
                        saved_state['proc_pending'][idx] = False
                        if os.path.exists(outdir + '/result_' + re.search('.*/(.*)$', outdir).group(1) + '.pklz'):
                            continue
                        
                        if os.path.exists(outdir):
                            shutil.rmtree(outdir)
                            
                    elif not saved_state['proc_done'][idx] and not saved_state['proc_pending'][idx]:
                        if os.path.exists(outdir + '/result_' + re.search('.*/(.*)$', outdir).group(1) + '.pklz'):
                            continue
                        
                        if os.path.exists(outdir):
                            shutil.rmtree(outdir)                        
                            
                saved_state['_taskid'] = saved_state['pending_tasks'][-1][0] - 1
                saved_state['_taskid']
                saved_state.pop('pending_tasks')
            else:
                indices = range(saved_state['refidx'].shape[0])
                for idx in indices:
                    outdir = saved_state['procs'][idx].output_dir()
                    if not saved_state['proc_done'][idx] and not saved_state['proc_pending'][idx]:
                        if os.path.exists(outdir + '/result_' + re.search('.*/(.*)$', outdir).group(1) + '.pklz'):
                            continue
                        
                        if os.path.exists(outdir):
                            shutil.rmtree(outdir)
                
            saved_state['_config']['local_hash_check'] = False
            self.__dict__.update(saved_state)
    
    def _save_state(self, stamp):
        temp = self.__dict__.copy()
        temp['proc_done'] = np.array(self.__dict__['proc_done'])
        temp['proc_pending'] = np.array(self.__dict__['proc_pending'])
        if 'pool' in temp:
            temp.pop('pool')
            temp.pop('processors')
            temp.pop('memory_gb')
            _task_obj = temp['_task_obj']
            temp['_task_obj'] = _task_obj.fromkeys(_task_obj, {})
        elif 'iparallel' in temp:
            temp.pop('client_args')
            temp.pop('iparallel')
            temp.pop('taskclient')
            temp.pop('taskmap')
        
        task = temp['plugin_args']['task']
        batch = temp['plugin_args']['batch']
        
        temp.pop('plugin_args')    
        temp.pop('timestamp')
        
        file_name = '/scratch/processed/reproducibility/checkpoints_' + task + '_batch_' + str(batch) + '/checkpoint_' + str(int(stamp)) + '.pkl'
        
        if not os.path.exists('/scratch/processed/reproducibility/checkpoints_' + task + '_batch_' + str(batch)):
            os.makedirs('/scratch/processed/reproducibility/checkpoints_' + task + '_batch_' + str(batch))
        
        with open(file_name, 'wb') as f:
            pickle.dump(temp, f)

    def run(self, graph, config, updatehash=False):
        """
        Executes a pre-defined pipeline using distributed approaches
        """
        self.timestamp = time()
        logger.info("Running in parallel.")
        self._config = config
        poll_sleep_secs = float(config["execution"]["poll_sleep_duration"])
        self.delete = []
        self._prerun_check(graph)
        # Generate appropriate structures for worker-manager model
        self._generate_dependency_list(graph)
        self.mapnodes = []
        self.mapnodesubids = {}
        # setup polling - TODO: change to threaded model
        notrun = []
        self.hours_elapsed = 0
        if str2bool(self._config["execution"]["remove_node_directories"]):
            self._load_state()
        
        old_progress_stats = None
        old_presub_stats = None
        while not np.all(self.proc_done) or np.any(self.proc_pending):
            loop_start = time()
            # Check if a job is available (jobs with all dependencies run)
            # https://github.com/nipy/nipype/pull/2200#discussion_r141605722
            jobs_ready = np.nonzero(~self.proc_done & (self.depidx.sum(0) == 0))[1]
            progress_stats = (
                len(self.proc_done),
                np.sum(self.proc_done ^ self.proc_pending),
                np.sum(self.proc_done & self.proc_pending),
                len(jobs_ready),
                len(self.pending_tasks),
                np.sum(~self.proc_done & ~self.proc_pending),
            )
            
            display_stats = progress_stats != old_progress_stats
            if display_stats:
                logger.debug(
                    "Progress: %d jobs, %d/%d/%d "
                    "(done/running/ready), %d/%d "
                    "(pending_tasks/waiting).",
                    *progress_stats
                )
                old_progress_stats = progress_stats
                
            toappend = []
            # trigger callbacks for any pending results
            while self.pending_tasks:
                taskid, jobid = self.pending_tasks.pop()
                try:
                    result = self._get_result(taskid)
                except Exception:
                    notrun.append(self._clean_queue(jobid, graph))
                else:
                    if result:
                        if result["traceback"]:
                            notrun.append(
                                self._clean_queue(jobid, graph, result=result)
                            )
                        else:
                            self._task_finished_cb(jobid)
                            self._remove_node_dirs()
                        self._clear_task(taskid)
                    else:
                        assert self.proc_done[jobid] and self.proc_pending[jobid]
                        toappend.insert(0, (taskid, jobid))
                
                num_jobs = len(self.pending_tasks)
                if num_jobs < self.max_jobs and (time() - loop_start) > poll_sleep_secs:
                    submitted_ = self._send_procs_to_workers(updatehash=updatehash, graph=graph)
                    if submitted_:
                        loop_start = time()
                    
            if toappend:
                self.pending_tasks.extend(toappend)

            num_jobs = len(self.pending_tasks)
            presub_stats = (num_jobs, np.sum(self.proc_done & self.proc_pending))
            display_stats = display_stats or presub_stats != old_presub_stats
            if display_stats:
                logger.debug("Tasks currently running: %d. Pending: %d.", *presub_stats)
                old_presub_stats = presub_stats
                
            if num_jobs < self.max_jobs and (time() - loop_start) > poll_sleep_secs:
                submitted_ = self._send_procs_to_workers(updatehash=updatehash, graph=graph)
            elif num_jobs > self.max_jobs:
                submitted_ = False
            else:
                submitted_ = False
                sleep_til = loop_start + poll_sleep_secs
                sleep(max(0, sleep_til - time()))
                
                
            if num_jobs < self.max_jobs and not submitted_:
                self._send_procs_to_workers(updatehash=updatehash, graph=graph)
            elif display_stats:
                logger.debug("Not submitting (max jobs reached)")
            
        self._remove_node_dirs()
        report_nodes_not_run(notrun)
        # close any open resources
        self._postrun_check()

    def _get_result(self, taskid):
        raise NotImplementedError

    def _submit_job(self, node, updatehash=False):
        raise NotImplementedError

    def _report_crash(self, node, result=None):
        tb = None
        if result is not None:
            node._result = result["result"]
            tb = result["traceback"]
            node._traceback = tb
        return report_crash(node, traceback=tb)

    def _clear_task(self, taskid):
        raise NotImplementedError

    def _clean_queue(self, jobid, graph, result=None):
        logger.debug("Clearing %d from queue", jobid)

        if self._status_callback:
            self._status_callback(self.procs[jobid], "exception")
        if result is None:
            result = {
                "result": None,
                "traceback": "\n".join(format_exception(*sys.exc_info())),
            }

        crashfile = self._report_crash(self.procs[jobid], result=result)
        if str2bool(self._config["execution"]["stop_on_first_crash"]):
            raise RuntimeError("".join(result["traceback"]))
        if jobid in self.mapnodesubids:
            # remove current jobid
            self.proc_pending[jobid] = False
            self.proc_done[jobid] = True
            # remove parent mapnode
            jobid = self.mapnodesubids[jobid]
            self.proc_pending[jobid] = False
            self.proc_done[jobid] = True
        # remove dependencies from queue
        return self._remove_node_deps(jobid, crashfile, graph)

    def _submit_mapnode(self, jobid):
        import scipy.sparse as ssp

        if jobid in self.mapnodes:
            return True
        self.mapnodes.append(jobid)
        mapnodesubids = self.procs[jobid].get_subnodes()
        numnodes = len(mapnodesubids)
        logger.debug("Adding %d jobs for mapnode %s", numnodes, self.procs[jobid])
        for i in range(numnodes):
            self.mapnodesubids[self.depidx.shape[0] + i] = jobid
        self.procs.extend(mapnodesubids)
        self.depidx = ssp.vstack(
            (self.depidx, ssp.lil_matrix(np.zeros((numnodes, self.depidx.shape[1])))),
            "lil",
        )
        self.depidx = ssp.hstack(
            (self.depidx, ssp.lil_matrix(np.zeros((self.depidx.shape[0], numnodes)))),
            "lil",
        )
        self.depidx[-numnodes:, jobid] = 1
        self.proc_done = np.concatenate(
            (self.proc_done, np.zeros(numnodes, dtype=bool))
        )
        return False

    def _send_procs_to_workers(self, updatehash=False, graph=None):
        """
        Sends jobs to workers
        """

        while not np.all(self.proc_done):
            num_jobs = len(self.pending_tasks)
            if np.isinf(self.max_jobs):
                slots = None
            else:
                slots = max(0, self.max_jobs - num_jobs)
            logger.debug("Slots available: %s", slots)
            if (num_jobs >= self.max_jobs) or (slots == 0):
                break

            # Check if a job is available (jobs with all dependencies run)
            # https://github.com/nipy/nipype/pull/2200#discussion_r141605722
            jobids = np.nonzero(~self.proc_done & (self.depidx.sum(0) == 0))[1]
            if len(jobids) > 0:
                # send all available jobs
                logger.info(
                    "Pending[%d] Submitting[%d] jobs Slots[%s]",
                    num_jobs,
                    len(jobids[:slots]),
                    slots or "inf",
                )

                for jobid in jobids[:slots]:
                    if isinstance(self.procs[jobid], MapNode):
                        try:
                            num_subnodes = self.procs[jobid].num_subnodes()
                        except Exception:
                            self._clean_queue(jobid, graph)
                            self.proc_pending[jobid] = False
                            continue
                        if num_subnodes > 1:
                            submit = self._submit_mapnode(jobid)
                            if not submit:
                                continue
                    # change job status in appropriate queues
                    self.proc_done[jobid] = True
                    self.proc_pending[jobid] = True
                    # Send job to task manager and add to pending tasks
                    logger.info("Submitting: %s ID: %d", self.procs[jobid], jobid)
                    if self._status_callback:
                        self._status_callback(self.procs[jobid], "start")

                    if not self._local_hash_check(jobid, graph):
                        if self.procs[jobid].run_without_submitting:
                            logger.debug(
                                "Running node %s on master thread", self.procs[jobid]
                            )
                            try:
                                self.procs[jobid].run()
                            except Exception:
                                self._clean_queue(jobid, graph)
                                
                            self._task_finished_cb(jobid)
                            self._remove_node_dirs()
                        else:
                            tid = self._submit_job(
                                deepcopy(self.procs[jobid]), updatehash=updatehash
                            )
                            if tid is None:
                                self.proc_done[jobid] = False
                                self.proc_pending[jobid] = False
                            else:
                                self.pending_tasks.insert(0, (tid, jobid))
                    logger.info(
                        "Finished submitting: %s ID: %d", self.procs[jobid], jobid
                    )
                    submitted = True
            else:
                submitted = False
                break
            
        return submitted

    def _local_hash_check(self, jobid, graph):
        if not str2bool(self.procs[jobid].config["execution"]["local_hash_check"]):
            return False

        try:
            cached, updated = self.procs[jobid].is_cached()
        except Exception:
            logger.warning(
                "Error while checking node hash, forcing re-run. "
                "Although this error may not prevent the workflow from running, "
                "it could indicate a major problem. Please report a new issue "
                "at https://github.com/nipy/nipype/issues adding the following "
                "information:\n\n\tNode: %s\n\tInterface: %s.%s\n\tTraceback:\n%s",
                self.procs[jobid],
                self.procs[jobid].interface.__module__,
                self.procs[jobid].interface.__class__.__name__,
                "\n".join(format_exception(*sys.exc_info())),
            )
            return False

        logger.debug(
            'Checking hash "%s" locally: cached=%s, updated=%s.',
            self.procs[jobid],
            cached,
            updated,
        )
        overwrite = self.procs[jobid].overwrite
        always_run = self.procs[jobid].interface.always_run

        if (
            cached
            and updated
            and (overwrite is False or overwrite is None and not always_run)
        ):
            logger.debug(
                "Skipping cached node %s with ID %s.", self.procs[jobid], jobid
            )
            try:
                self._task_finished_cb(jobid, cached=True)
                self._remove_node_dirs()
            except Exception:
                logger.debug(
                    "Error skipping cached node %s (%s).\n\n%s",
                    self.procs[jobid],
                    jobid,
                    "\n".join(format_exception(*sys.exc_info())),
                )
                self._clean_queue(jobid, graph)
                self.proc_pending[jobid] = False
            return True
        return False

    def _task_finished_cb(self, jobid, cached=False):
        """Extract outputs and assign to inputs of dependent tasks

        This is called when a job is completed.
        """
        logger.info(
            "[Job %d] %s (%s).",
            jobid,
            "Cached" if cached else "Completed",
            self.procs[jobid],
        )
        self._short_circuit_results(jobid)
        if self._status_callback:
            self._status_callback(self.procs[jobid], "end")
        # Update job and worker queues
        self.proc_pending[jobid] = False
        # update the job dependency structure
        rowview = self.depidx.getrowview(jobid)
        rowview[rowview.nonzero()] = 0
        if jobid not in self.mapnodesubids:
            self.refidx[self.refidx[:, jobid].nonzero()[0], jobid] = 0
        
        if str2bool(self._config["execution"]["remove_node_directories"]):
            if (time() - self.timestamp) / 3600 > 1:
                self.hours_elapsed += 1
                self._save_state(self.hours_elapsed)
                self.timestamp = time()

    def _generate_dependency_list(self, graph):
        """Generates a dependency list for a list of graphs."""
        import networkx as nx

        self.procs, _ = topological_sort(graph)
        self.depidx = nx.to_scipy_sparse_matrix(
            graph, nodelist=self.procs, format="lil"
        )
        self.refidx = self.depidx.astype(int)
        self.proc_done = np.zeros(len(self.procs), dtype=bool)
        self.proc_pending = np.zeros(len(self.procs), dtype=bool)

    def _remove_node_deps(self, jobid, crashfile, graph):
        import networkx as nx

        try:
            dfs_preorder = nx.dfs_preorder
        except AttributeError:
            dfs_preorder = nx.dfs_preorder_nodes
        subnodes = [s for s in dfs_preorder(graph, self.procs[jobid])]
        for node in subnodes:
            idx = self.procs.index(node)
            self.proc_done[idx] = True
            self.proc_pending[idx] = False

        return dict(node=self.procs[jobid], dependents=subnodes, crashfile=crashfile)
    
    def _short_circuit_results(self, jobid):
        node_dir = self.procs[jobid].output_dir()
        results_file = glob(os.path.join(node_dir, "result_*.pklz"))[0]
        result_data = load_resultfile(results_file)
        result_out = dict(result=None, traceback=None)
        if isinstance(result_data, dict):
            result_out["result"] = result_data["result"]
            result_out["traceback"] = result_data["traceback"]
            result_out["hostname"] = result_data["hostname"]
            if results_file:
                crash_file = os.path.join(node_dir, "crashstore.pklz")
                os.rename(results_file, crash_file)
        else:
            result_out["result"] = result_data
            
        self._check_results(result_out, jobid)

    def _flatten(self, l):
        if isinstance(l, dict):
            l_ = l.values()
        elif not isinstance(l, Iterable):
            l_ = [l]
        elif isinstance(l, (str, bytes)):
            l_ = [l]
        else:
            l_ = l
            
        for el in l_:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from self._flatten(el)
            else:
                yield el

    def _check_results(self, result, jobid):
        if str2bool(self._config["execution"]["remove_node_directories"]):
            if result:
                for item in list(result['result'].outputs.items()):
                    it = vars(result['result'].outputs).get(item[0])
                    if it == None:
                        it = result['result'].outputs.get()[item[0]]
                    it = list(self._flatten(it))
                    for key in self._flatten(result['result'].inputs):
                        if os.path.isfile(str(key)) or os.path.isdir(str(key)):
                            if isinstance(it, dict) and key in it.values():
                                vars(self)['_'+str(jobid)] = list(self.refidx[:, jobid].nonzero()[0])
                                self.delete += list(self.refidx[:, jobid].nonzero()[0]) #NODES THAT LEAD INTO CURRENT JOBID NODE
                                return
                            elif isinstance(it, (list, tuple, Iterable)) and (key in it or any(it_.startswith(key) for it_ in it if isinstance(it_, str))):
                                vars(self)['_'+str(jobid)] = list(self.refidx[:, jobid].nonzero()[0])
                                self.delete += list(self.refidx[:, jobid].nonzero()[0])
                                return
                            elif isinstance(it, str) and (key == it or it.startswith(key)):
                                vars(self)['_'+str(jobid)] = list(self.refidx[:, jobid].nonzero()[0])
                                self.delete += list(self.refidx[:, jobid].nonzero()[0])
                            else:
                                continue
                    

    def _graph_expansion(self, lst, searched):
        lst_out = lst.copy()
        for idx in lst:
            if idx in searched:
                continue
            lst_out += vars(self).get('_'+str(idx), [])
            searched.append(idx)
            
        return len(set(searched)) != len(set(lst_out)), lst_out, searched


    def _remove_node_dirs(self):
        """Removes directories whose outputs have already been used up"""
        if str2bool(self._config["execution"]["remove_node_directories"]):
            indices = np.nonzero((self.refidx.sum(axis=1) == 0).__array__())[0]
            indices_check = np.setdiff1d(list(set(self.delete)), indices)
            indices_count = Counter(self.delete)
            indices = np.setdiff1d(indices, list(set(self.delete)))
            for idx in indices:
                if idx in self.mapnodesubids:
                    continue
                if self.proc_done[idx] and (not self.proc_pending[idx]):
                    lst = [idx]
                    cont = True
                    searched = []
                    while cont:
                        cont, lst, searched = self._graph_expansion(lst, searched)
                        
                    for i in lst:
                        if self.proc_done[i] and (not self.proc_pending[i]):
                            if i in indices_check or i not in indices:
                                continue
                            
                            outdir = self.procs[i].output_dir()
                            if not os.path.exists(outdir):
                                continue
                            
                            if i in self.delete:
                                self.delete.remove(i)
                                if indices_count[i]:
                                    indices_count[i] -= 1
                                    continue
                            #MIGHT JUST REVERT TO DELETING EVERYTHING FROM SET
                            if vars(self).get('_'+str(i), '') != '':
                                #self.delete = list(np.setdiff1d(self.delete, vars(self)['_'+str(i)]))
                                for i_ in vars(self).get('_'+str(i)):
                                    if i_ in self.delete:
                                        self.delete.remove(i_)
                                        if indices_count[i_]:
                                            indices_count[i_] -= 1
                                del vars(self)['_'+str(i)]
                                
                            self.refidx[i, i] = -1
                            
                            logger.info(
                                (
                                    "[node dependencies finished] "
                                    "removing node: %s from directory %s"
                                )
                                % (self.procs[i]._id, outdir)
                            )
                            shutil.rmtree(outdir)
                            

class SGELikeBatchManagerBase(DistributedPluginBase):
    """Execute workflow with SGE/OGE/PBS like batch system"""

    def __init__(self, template, plugin_args=None):
        super(SGELikeBatchManagerBase, self).__init__(plugin_args=plugin_args)
        self._template = template
        self._qsub_args = None
        if plugin_args:
            if "template" in plugin_args:
                self._template = plugin_args["template"]
                if os.path.isfile(self._template):
                    with open(self._template) as tpl_file:
                        self._template = tpl_file.read()
            if "qsub_args" in plugin_args:
                self._qsub_args = plugin_args["qsub_args"]
        self._pending = {}

    def _is_pending(self, taskid):
        """Check if a task is pending in the batch system"""
        raise NotImplementedError

    def _submit_batchtask(self, scriptfile, node):
        """Submit a task to the batch system"""
        raise NotImplementedError

    def _get_result(self, taskid):
        if taskid not in self._pending:
            raise Exception("Task %d not found" % taskid)
        if self._is_pending(taskid):
            return None
        node_dir = self._pending[taskid]
        # MIT HACK
        # on the pbs system at mit the parent node directory needs to be
        # accessed before internal directories become available. there
        # is a disconnect when the queueing engine knows a job is
        # finished to when the directories become statable.
        t = time()
        timeout = float(self._config["execution"]["job_finished_timeout"])
        timed_out = True
        while (time() - t) < timeout:
            try:
                glob(os.path.join(node_dir, "result_*.pklz")).pop()
                timed_out = False
                break
            except Exception as e:
                logger.debug(e)
            sleep(2)
        if timed_out:
            result_data = {"hostname": "unknown", "result": None, "traceback": None}
            results_file = None
            try:
                error_message = (
                    "Job id ({0}) finished or terminated, but "
                    "results file does not exist after ({1}) "
                    "seconds. Batch dir contains crashdump file "
                    "if node raised an exception.\n"
                    "Node working directory: ({2}) ".format(taskid, timeout, node_dir)
                )
                raise IOError(error_message)
            except IOError as e:
                result_data["traceback"] = "\n".join(format_exception(*sys.exc_info()))
        else:
            results_file = glob(os.path.join(node_dir, "result_*.pklz"))[0]
            result_data = load_resultfile(results_file)
        result_out = dict(result=None, traceback=None)
        if isinstance(result_data, dict):
            result_out["result"] = result_data["result"]
            result_out["traceback"] = result_data["traceback"]
            result_out["hostname"] = result_data["hostname"]
            if results_file:
                crash_file = os.path.join(node_dir, "crashstore.pklz")
                os.rename(results_file, crash_file)
        else:
            result_out["result"] = result_data
        return result_out
    
    def _submit_job(self, node, updatehash=False):
        """submit job and return taskid"""
        pyscript = create_pyscript(node, updatehash=updatehash)
        batch_dir, name = os.path.split(pyscript)
        name = ".".join(name.split(".")[:-1])
        batchscript = "\n".join(
            (self._template.rstrip("\n"), "%s %s" % (sys.executable, pyscript))
        )
        batchscriptfile = os.path.join(batch_dir, "batchscript_%s.sh" % name)
        with open(batchscriptfile, "wt") as fp:
            fp.writelines(batchscript)
        return self._submit_batchtask(batchscriptfile, node)

    def _clear_task(self, taskid):
        del self._pending[taskid]


class GraphPluginBase(PluginBase):
    """Base class for plugins that distribute graphs to workflows"""

    def __init__(self, plugin_args=None):
        if plugin_args and plugin_args.get("status_callback"):
            logger.warning(
                "status_callback not supported for Graph submission" " plugins"
            )
        super(GraphPluginBase, self).__init__(plugin_args=plugin_args)

    def run(self, graph, config, updatehash=False):
        import networkx as nx

        pyfiles = []
        dependencies = {}
        self._config = config
        nodes = list(nx.topological_sort(graph))
        logger.debug("Creating executable python files for each node")
        for idx, node in enumerate(nodes):
            pyfiles.append(
                create_pyscript(node, updatehash=updatehash, store_exception=False)
            )
            dependencies[idx] = [
                nodes.index(prevnode) for prevnode in list(graph.predecessors(node))
            ]
        self._submit_graph(pyfiles, dependencies, nodes)

    def _get_args(self, node, keywords):
        values = ()
        for keyword in keywords:
            value = getattr(self, "_" + keyword)
            if keyword == "template" and os.path.isfile(value):
                with open(value) as f:
                    value = f.read()
            if (
                hasattr(node, "plugin_args")
                and isinstance(node.plugin_args, dict)
                and keyword in node.plugin_args
            ):
                if keyword == "template" and os.path.isfile(node.plugin_args[keyword]):
                    with open(node.plugin_args[keyword]) as f:
                        tmp_value = f.read()
                else:
                    tmp_value = node.plugin_args[keyword]

                if "overwrite" in node.plugin_args and node.plugin_args["overwrite"]:
                    value = tmp_value
                else:
                    value += tmp_value
            values += (value,)
        return values
    
    def _submit_graph(self, pyfiles, dependencies, nodes):
        """
        pyfiles: list of files corresponding to a topological sort
        dependencies: dictionary of dependencies based on the toplogical sort
        """
        raise NotImplementedError

    def _get_result(self, taskid):
        if taskid not in self._pending:
            raise Exception("Task %d not found" % taskid)
        if self._is_pending(taskid):
            return None
        node_dir = self._pending[taskid]

        glob(os.path.join(node_dir, "result_*.pklz")).pop()

        results_file = glob(os.path.join(node_dir, "result_*.pklz"))[0]
        result_data = load_resultfile(results_file)
        result_out = dict(result=None, traceback=None)

        if isinstance(result_data, dict):
            result_out["result"] = result_data["result"]
            result_out["traceback"] = result_data["traceback"]
            result_out["hostname"] = result_data["hostname"]
            if results_file:
                crash_file = os.path.join(node_dir, "crashstore.pklz")
                os.rename(results_file, crash_file)
        else:
            result_out["result"] = result_data

        return result_out

