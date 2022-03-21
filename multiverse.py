#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 14:45:59 2021

@author: grahamseasons
"""
import subprocess
import argparse, sys, os, json, pickle, re
from multiverse.gui.gui import MultiverseConfig
from os.path import join as opj
import numpy as np

dir = os.path.dirname(os.path.abspath(__file__))

def parse(start):
    args = start.parse_args()
    run_now = False
    process_mode = 'MultiProc'
    if args.config or args.run:
        if args.data == None and not args.config:
            print('Either -d must be set, or -c must be used')
            sys.exit()
        elif args.data != None and args.run and not os.path.isdir(args.data):
            print('The specified data path does not exist')
            sys.exit()
            
        if args.out == None and not args.config:
            print('Either -o must be set, or -c must be used')
            sys.exit()
        elif args.out != None and args.run and not os.path.isdir(args.out):
            print('The specified data path does not exist')
            print('Creating directory at specified path')
            os.makedirs(args.out)
            
        if args.run:
            if os.path.isdir(os.path.join(dir, 'multiverse', 'configuration')):
                if os.path.isfile(os.path.join(dir, 'multiverse', 'configuration', 'multiverse_configuration.pkl')) and os.path.isfile(os.path.join(dir, 'multiverse', 'configuration', 'general_configuration.pkl')):
                    run_now = True
                    
                    with open(opj(dir, 'multiverse', 'configuration', 'general_configuration.pkl'), 'rb') as f:
                        configure = pickle.load(f)
                    
                    process_mode = configure['processing']
                    print('Using previously defined configuration file')
            
            if not run_now:
                print('Configuration files missing, running configure')
                config = MultiverseConfig(args.rerun, args.data, args.out)
                run_now = config.run_now
                process_mode = config.configure['processing']
                configure = config.configure
        #MAKE SO THAT USING ARGS.RERUN WILL READ CONFIG FILE, AND SAVE WHATEVER TRUE OR FALSE WITH IT
        if args.config:
            config = MultiverseConfig(args.rerun, args.data, args.out)
            run_now = config.run_now
            process_mode = config.configure['processing']
            configure = config.configure
    else:
        print('Either the --run (-r) or --config (-c) flag must be specified\n')
        start.print_help()
        sys.exit()
    
    return run_now, args, process_mode, configure

def main():
    start = argparse.ArgumentParser()
    start.add_argument('-c', '--config', action='store_true', help='configure multiverse parameter file')
    start.add_argument('-r', '--run', action='store_true', help='run multiverse analysis')
    start.add_argument('-rr', '--rerun', action='store_true', help='re-run multiverse analysis using saved population files')
    start.add_argument('-d', '--data', type=str, metavar="DATA_DIR", action='store', help='path to BIDS formatted data directory')
    start.add_argument('-o', '--out', type=str, metavar="OUT_DIR", action='store', help='path to store outputs')

    run_now, args, process_mode, config = parse(start)
    
    if run_now:
        code_dir = os.path.join(dir, 'multiverse')
        print(dir)
        volumes = ['{code}:/multiverse/code'.format(code=code_dir), '{data}:/data'.format(data=args.data), 
                   '{work_dir}:/scratch'.format(work_dir=args.out), 'plugins_base.py:/opt/miniconda-latest/envs/multiverse/lib/python3.8/site-packages/nipype/pipeline/plugins/base.py']
        if process_mode != 'SLURM':
            try:
                import docker
            except:
                subprocess.call(['pip', 'install', 'docker'])
            client = docker.from_env()
            try:
                print('Running Container')
                container = client.containers.run('gseasons/multiverse:cluster', detach=True, tty=True, stdin_open=True, working_dir='/scratch', volumes=volumes, user='root')
                container.start()
                container.exec_run('sudo /bin/bash -c "source activate multiverse ; python /multiverse/code/run_multiverse.py {0}"'.format(args.rerun))
                container.stop()
                container.remove()
                client.volumes.prune()
                
            except docker.errors.ImageNotFound:
                print("Image not found, pulling from docker")
                for line in client.api.pull('gseasons/multiverse', stream=True, decode=True):
                    print(json.dumps(line, indent=4))
                    
                container = client.containers.run('gseasons/multiverse:cluster', detach=True, tty=True, stdin_open=True, working_dir='/scratch', volumes=volumes)
                container.start()
                container.exec_run('/bin/bash -c "source activate multiverse ; python /multiverse/code/run_multiverse.py {0}"'.format(args.rerun))
                container.stop()
                container.remove()
                client.volumes.prune()
        else:
            if config['account'] == 'def-':
                config['account'] = ''
            #config['ntasks'] = 24
            #config['account'] = 'def-emazerol'
            #config['time'] = '0-05:00'
            #config['mem'] = '40000'
            #REQUEST FULL NODES INSTEAD OF BUNCH OF CPUs
            #subprocess.call(['sbatch', '--ntasks={0}'.format(config['ntasks']), '--account={0}'.format(config['account']), '--time={0}'.format(config['time']), '--mem-per-cpu={0}'.format(config['mem']), 'multiverse/configuration/multiverse.sh', args.data, args.out, str(args.rerun)])
            if config['processing'] == 'SLURM' and 'num_generations' not in config:
                config['nodes'] = '1'
                config['ntasks'] = config['cpu_node']
                
                split_ = re.split('-|:', config['time'])
                print(split_)
                time_s = 0
                for k, sp in enumerate(split_):
                    if not k:
                        time_s += int(sp) * 24 * 3600
                    elif k == 1:
                        time_s += int(sp) * 3600
                    elif k == 2:
                        time_s += int(sp) * 60
                    else:
                        time_s += int(sp)
                time_s = time_s / config['pipelines'] * np.ceil(int(config['cpu_node']) / 4).astype(int)
                
                days = int(time_s / (24 * 3600))
                hours = max(int(time_s / 3600) - days * 24, 0)
                minutes = max(int(time_s / 60) - days * 24 * 60 - hours * 60, 0)
                seconds = max(int(time_s) - days * 24 * 3600 - hours * 3600 - minutes * 60, 0)
                
                config['time'] =  str(days) + '-' + str(hours) + ':' + str(minutes) + ':' + str(seconds)
                config['mem'] = np.ceil(float(config['mem']) * int(config['cpu_node'])).astype(int).astype(str)
                config['batch'] = np.ceil(int(config['cpu_node']) / 4).astype(int)
                
                subprocess.call(['{0}/configuration/intermediate.sh'.format(code_dir), args.data, args.out, str(args.rerun), config['nodes'], config['ntasks'], config['account'], config['time'], config['mem']])
            else:
                subprocess.call(['sbatch', '--nodes={0}'.format(config['nodes']), '--ntasks={0}'.format([config['ntasks']]), '--account={0}'.format(config['account']), '--time={0}'.format(config['time']), '--mem={0}'.format(config['mem']), 'multiverse/configuration/multiverse.sh', args.data, args.out, str(args.rerun)])
            

if __name__ == "__main__":
    main()