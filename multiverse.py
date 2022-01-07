#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 14:45:59 2021

@author: grahamseasons
"""
import subprocess
import argparse, sys, os, json
from code.gui.gui import MultiverseConfig

dir = os.path.dirname(os.path.abspath(__file__))

def parse(start):
    args = start.parse_args()
    run_now = False
    process_mode = 'MultiProc'
    if args.config or args.run:
        if args.data == None and not args.config:
            print('Either -d must be set, or -c must be used')
            sys.exit()
        elif not os.path.isdir(args.data):
            print('The specified data path does not exist')
            sys.exit()
            
        if args.out == None and not args.config:
            print('Either -o must be set, or -c must be used')
            sys.exit()
        elif not os.path.isdir(args.out):
            print('The specified data path does not exist')
            print('Creating directory at specified path')
            os.makedirs(args.out)
            
        if args.run:
            if os.path.isdir(os.path.join(dir, 'code', 'configuration')):
                if os.path.isfile(os.path.join(dir, 'code', 'configuration', 'multiverse_configuration.pkl')) and os.path.isfile(os.path.join(dir, 'code', 'configuration', 'general_configuration.pkl')):
                    run_now = True
                    print('Using previously defined configuration file')
            
            if not run_now:
                print('Configuration files missing, running configure')
                config = MultiverseConfig(args.rerun, args.data, args.out)
                run_now = config.run_now
                process_mode = config.configure['processing']
        #MAKE SO THAT USING ARGS.RERUN WILL READ CONFIG FILE, AND SAVE WHATEVER TRUE OR FALSE WITH IT
        if args.config:
            config = MultiverseConfig(args.rerun, args.data, args.out)
            run_now = config.run_now
            process_mode = config.configure['processing']
    else:
        print('Either the --run (-r) or --config (-c) flag must be specified\n')
        start.print_help()
        sys.exit()
    
    return run_now, args, process_mode

def main():
    start = argparse.ArgumentParser()
    start.add_argument('-c', '--config', action='store_true', help='configure multiverse parameter file')
    start.add_argument('-r', '--run', action='store_true', help='run multiverse analysis')
    start.add_argument('-rr', '--rerun', action='store_true', help='re-run multiverse analysis using saved population files')
    start.add_argument('-d', '--data', type=str, metavar="DATA_DIR", action='store', help='path to BIDS formatted data directory')
    start.add_argument('-o', '--out', type=str, metavar="OUT_DIR", action='store', help='path to store outputs')

    run_now, args, process_mode = parse(start)
    
    if run_now:
        code_dir = os.path.join(dir, 'code')
        print(dir)
        volumes = ['{code}:/multiverse/code'.format(code=code_dir), '{data}:/data'.format(data=args.data), '{work_dir}:/scratch'.format(work_dir=args.out)]
        if process_mode != 'SLURM':
            import docker
            client = docker.from_env()
            try:
                print('Running Container')
                container = client.containers.run('gseasons/multiverse', detach=True, tty=True, stdin_open=True, working_dir='/scratch', volumes=volumes, user='root')
                container.start()
                container.exec_run('sudo /bin/bash -c "source activate multiverse ; cp -f /multiverse/code/plugins_base.py /opt/miniconda-latest/envs/multiverse/lib/python3.8/site-packages/nipype/pipeline/plugins/base.py ; python /multiverse/code/run_multiverse.py"')
                container.stop()
                container.remove()
                client.volumes.prune()
                
            except docker.errors.ImageNotFound:
                print("Image not found, pulling from docker")
                for line in client.api.pull('gseasons/multiverse', stream=True, decode=True):
                    print(json.dumps(line, indent=4))
                    
                container = client.containers.run('gseasons/multiverse', detach=True, tty=True, stdin_open=True, working_dir='/scratch', volumes=volumes)
                container.start()
                container.exec_run('/bin/bash -c "source activate multiverse ; python /multiverse/code/run_multiverse.py"')
                container.stop()
                container.remove()
                client.volumes.prune()
        else:
            subprocess.call(['sbatch', 'code/configuration/multiverse.sh', args.data, args.out])
            

if __name__ == "__main__":
    main()