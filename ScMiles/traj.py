#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:20:21 2019

@author: Wei Wei

Running free trajectories.

"""

import os, time, re, sys
from run import *
from parameters import *
from colvar import *
from datetime import datetime
from log import log
from milestones import *
from shutil import copyfile
import inspect 
from itertools import combinations
import numpy as np
import math


class trajPool: 
    def __init__(self, dest) -> None:
        self.__dest__ = list(map(int, re.findall('\d+', dest)))
        
    def __enter__(self):
        return self
               
    def __exit__(self, exc_type, exc_value, traceback):
        return 
    
    def __repr__(self) -> str:
        return ('Current reached by {} milestones'
                .format(self.count_ms()))
        
    def count_ms(self):
        """
        return:
            1. the total number of milestones that have trajectories reached 
            current milestone
            2. a list of these milestones
        example;
            milestone 2_3 has trajectories from milestone 1_2 and 3_4
            return 2, ['1_2', '3_4']
        """
        attributes = inspect.getmembers(self, lambda a: not(inspect.isroutine(a)))
        att = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
        lst = []
        for i in range(len(att)):
            lst.append(att[i][0])
        return len(att), lst

    def count_total_traj(self):
        """
        count how many trajectories reached current milestone

        example:
            milestone 2_3 has trajectories from milestone 1_2 and 3_4, 120 trajectories totally.
            return 120
        """
        total = 0
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        att = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
        for i in range(len(att)):
            total += len(self.count_traj(att[i][0]))
        return total
    
    def count_traj(self, ms):
        """
        return the trajectories from milestone 'ms'
        """
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        att = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
        for i in range(len(att)):
            if att[i][0] == ms:
                return att[i][1] 
        return []
    
    def get_distribution(self, kij, index, q):
        """
        argv:
            kij: K matrix
              
            index: the index of K matrix. i.e. {0:[1,2]...} means 
                   the first row/column in K associates to milestone 1_2
                   
            q: flux
            
        return:
            the distribution of trajectories for next iteration.
            e.g. {'1_2':0.6, '3_4':0.4} indicates 60% of the trajectories 
            are from milestone 1_2, while 40% from milestone 3_4
        """
        total_orig, orig_list = self.count_ms() 
        distribution = {}
        distribution_prob = {}   
        # print(self.__dest__)
        # print(index)
        dest_idx = [i for i, a in enumerate(index) if a == self.__dest__]
        if not dest_idx:
            return 0
        else:
            dest_idx = dest_idx[0]
        for orig_ms in orig_list:
            orig = list(map(int, (re.findall('\d+', orig_ms))))
            orig_index = [i for i, a in enumerate(index) if a == orig][0]
            if len(orig_list) == 1:
                distribution[orig_ms] = 1.0
            else:
                if math.isnan(q[orig_index]):
                    distribution[orig_ms] = 0.0
                else:
                    distribution[orig_ms] = (kij[orig_index, dest_idx] * np.abs(q[orig_index]))
        if sum(distribution.values()) == 0:
            print("Wrong distribution on {}.".format(self.__dest__))
            # print("Exiting...")
            log("Wrong distribution on {}. Maybe due to no trajectory reached or wrong flux on the adjacent milestone"
                .format(self.__dest__))
        for orig_ms in orig_list:   
            distribution_prob[orig_ms] = distribution[orig_ms] / sum(distribution.values())
        return distribution_prob
        

class traj:
    def __init__(self, parameter, jobs):
        self.parameter = parameter
        self.jobs = jobs
        self.pool = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return
    
    def __repr__(self) -> str:
        return ('Free trajectories.')    

    def __snapshots(self, MSname):
        '''
        param MSname:
            milestone name
        return:
            next frame number
        example:
            250 free traj simulations for milestone MSname, return 251 as the next traj number for folder name
        '''
        lst = re.findall('\d+',MSname)
        name = lst[0] + '_' + lst[1]
        next_frame = 1
        while True:
            path = self.parameter.crdPath + '/' + name + '/' + str(self.parameter.iteration) + '/' + str(next_frame) 
            if os.path.exists(path):
                next_frame += 1
            else:
                return next_frame
    
    def check(self):
        '''check running jobs on cluster'''
        while True:
            name = 'MILES'
            if self.jobs.check(JobName=name) == False:
                pass
            else:
                return True
            print("Free trajectories are running... Next check in 200 seconds. {}".format(str(datetime.now())))
            time.sleep(200)  

    def __copy(self, orig, new, restart=1, enhanced=False):
        '''copy restart files for next iteration'''
        lst = ['.coor', '.vel', '.xsc']
        if not os.path.exists(new):
            os.makedirs(new)
        for ext in lst:
            newName = new + '/' + self.parameter.outputname + ext
            if restart == 1:
                oldName = orig + '/' + self.parameter.outputname + '.restart' + ext 
            else:
                oldName = orig + '/' + self.parameter.outputname + ext
            with open(new + '/original.txt', 'w+') as f:
                print('This trajectory is from', file=f)  
                print(orig, file=f) 
            if enhanced:
                with open(new + '/enhanced', 'w+') as f:
                    print('This trajectory is from enhancement', file=f)  
            copyfile(oldName, newName)

    def __iteration_prepare(self, path, final_ms, orig):
        '''add each trajectory to trajPool based on the ending milestone'''
        if final_ms == [0,0]:
            return
        if not os.path.isfile(path + '/' + self.parameter.outputname + '.restart.coor'):
            return
        ms = str(final_ms[0]) + '_' + str(final_ms[1])
        MSname = 'MS' + ms
        MSorig = 'MS' + orig
#        if self.parameter.ignorNewMS and MSname not in self.parameter.MS_list:
#            return
        self.add_to_trajPool(MSname, MSorig, path)

    def iteration_initialize(self):
        '''for each milestone, initialize trajectories based on distribution.'''
        import pandas as pd
        from compute import compute
#        print("Iteration # {}".format(self.parameter.iteration))
#        log("Iteration # {}".format(self.parameter.iteration))
        self.initialize_trajPool()
        
        for MSname in self.parameter.MS_list:
            name = re.findall('\d+', MSname)[0] + '_' + re.findall('\d+', MSname)[1]
            path = self.parameter.crdPath + '/' + name + '/' + str(self.parameter.iteration - 1) + '/'
            for j in range(self.parameter.trajPerLaunch):
                traj_path = path + str(j+1)
                if not os.path.exists(traj_path):
                    continue
                if self.parameter.restart == True:
                    if os.path.exists(self.parameter.crdPath + '/' + name + '/' + str(self.parameter.iteration) + '/distribution'):
                        continue
                end_file = path + str(j+1) + '/end.txt'
                if not os.path.isfile(end_file) or os.stat(end_file).st_size == 0:
                    continue
                final_ms = pd.read_csv(path + str(j+1) + '/end.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                final_ms = [int(x) for x in final_ms]
                ms = 'MS' + str(final_ms[0]) + '_' + str(final_ms[1])
                if ms not in self.parameter.MS_list:
                    continue   
#                print(final_ms)
                self.__iteration_prepare(traj_path, final_ms, name)

        kij, index, flux = compute(self.parameter)
        
#        print(MS10_11)
#        print(MS10_11.count_traj('MS11_12'))
        for MSname in self.parameter.MS_list:
            name = re.findall('\d+', MSname)[0] + '_' + re.findall('\d+', MSname)[1]
            path = self.parameter.crdPath + '/' + name + '/' + str(self.parameter.iteration) + '/'
            if self.parameter.restart == True:
                if os.path.exists(path + 'distribution'):
                    continue
            distribution = globals()[MSname].get_distribution(kij, index, flux)
            if distribution == 0:
                continue
            if not os.path.exists(path):
                os.makedirs(path)
            path = path + 'distribution'
            with open(path, 'w+') as f:
                print(distribution, file=f)  
#            print(MSname)
#            print(distribution)
            total_orig, orig_list = globals()[MSname].count_ms() 
            if total_orig == 0:
                print("No traj reached {}.".format(MSname))
                # print("Exitting...")
                log("No traj reached {}.".format(MSname))
                log("Skip {} for next iteration.".format(MSname))
                continue
                # sys.exit()
            total_traj = 0
#            print(MSname, globals()[MSname].count_ms())
            for ms in distribution.keys():
                traj_num = int(distribution[ms] * self.parameter.trajPerLaunch)
                if traj_num < 1:
                    traj_num = 1
                if ms == list(distribution.keys())[-1]:
                    traj_num = self.parameter.trajPerLaunch - total_traj
                if traj_num <= 0:
                    continue
                traj_lst = list(np.random.choice(globals()[MSname].count_traj(ms), traj_num))
                if len(traj_lst)  == 0:
                    with open('error', 'a+') as f:
                        print(MSname,ms, file=f)
                        print(globals()[MSname].count_traj(ms),file=f)
                        print(traj_num,file=f)
                self.copy_traj(MSname, traj_lst, total_traj)
                total_traj += traj_num
            del globals()[MSname]
        self.initialize_trajPool()
#        sys.exit()
    
    def copy_traj(self, MSname, traj_lst, total_traj):
        '''copy restart files'''
        [anchor1, anchor2] = list(map(int,(re.findall('\d+', MSname))))
        name = str(anchor1) + '_' + str(anchor2)
        filePath = os.path.dirname(os.path.abspath(__file__)) 
        pardir = os.path.abspath(os.path.join(filePath, os.pardir))
        path = pardir +  '/crd/' + name + '/' + str(self.parameter.iteration) + '/'  
        if not os.path.exists(path):
            os.makedirs(path)    
        continued = []
        for i, traj in enumerate(traj_lst):
            count = i + 1 + total_traj
            newPath = path + str(count)
            if not os.path.exists(newPath):
                os.makedirs(newPath)
            if traj not in continued:
                self.__copy(traj, newPath, restart=1)
                continued.append(traj)
            else:
                self.__copy(traj, newPath, restart=1, enhanced=True)
        with open(path + 'enhanced', 'a+') as f:
            lst = re.findall('\d+',traj_lst[0])
            ms = lst[-4] + '_' + lst[-3]
            print('Milestone {} has {} trajectories from enhancement.'\
                  .format(ms, len(traj_lst) - len(continued)), file=f)
    
    def edit_submit_scripts(self, path, traj_frame):
        #1 is done, 0 needs to restart
        import os
        from fileinput import FileInput
        if traj_frame == 1:
            beg_number = 1
            end_number = self.parameter.trajPerLaunch + 1
        else:
            beg_number = traj_frame - self.parameter.trajPerLaunch
            end_number = traj_frame
        for i in range(beg_number, end_number):
            next_line = False
            run = False
            current_path = path + '/' + str(i) + '/submit'
            if os.path.isfile(current_path):
                with FileInput(files=current_path, inplace=True) as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('##'):
                            print(line)
                            continue
                        if next_line == True:
                            line = '##' + line
                            next_line = False
                        elif 'cd' in line and 'crd' not in line:
                            info = line.split('/')
                            traj_num = int(info[-1])
                            if os.path.isfile(path + '/' + str(traj_num) + '/stop.colvars.state'):
                                line = '##' + line
                                next_line = True
                            else:
                                run = True
                        print(line)
                if run == False:
                    if not os.path.exists(path + '/' + str(i) + '/submit_done'):
                        os.rename(current_path, path + '/' + str(i) + '/submit_done')
                    
    
    def launch(self, joblist=None, lastLog=None):
        '''launch free trajectories'''
        print("Iteration # {}".format(self.parameter.iteration))
        log("Iteration # {}".format(self.parameter.iteration))
        
        if self.parameter.method == 1 and self.parameter.iteration > 1:
            if self.parameter.restart == False or (self.parameter.restart == True and "Iteration #" in lastLog):
                self.iteration_initialize()
                log('Iteration initialize complete')

        for ms in self.pool:
            self.pool[ms] = 0
        MS_list = self.parameter.MS_list.copy()
        colvar(self.parameter, free='yes').generate()
        
        for name in self.parameter.MS_list:
            [anchor1, anchor2] = list(map(int,(re.findall('\d+', name))))
            MSname = str(anchor1) + '_' + str(anchor2)
            if self.parameter.method == 1:
                next_frame = 1
            else:
                next_frame = self.__snapshots(name)
            next_script = next_frame
            if self.parameter.restart == True and os.path.exists(self.parameter.crdPath + '/' + MSname + '/' + str(self.parameter.iteration)):
                self.edit_submit_scripts(self.parameter.crdPath + '/' + MSname + '/' + str(self.parameter.iteration), next_frame)
            else:
                path_name = self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2) + '/' + str(self.parameter.iteration)
                create_folder(path_name)
                for j in range(self.parameter.trajPerLaunch):
                    folderPath = path_name + "/" + str(j+next_frame)
                    submit = False
                    frame = self.parameter.startTraj + next_frame // self.parameter.trajPerLaunch + j * self.parameter.interval
                    if frame > self.parameter.nframe:
                        break
                    if j + next_frame == next_script:
                        submit = True
                        next_script = next_script + self.parameter.traj_per_script[1]
                    create_folder(folderPath)
                    milestones(self.parameter).get_initial_ms(folderPath)
                    self.jobs.prepare_trajectory(a1=anchor1, a2=anchor2, snapshot=j+next_frame, frame=frame, script=submit)
            
        for name in self.parameter.MS_list:
            [anchor1, anchor2] = list(map(int,(re.findall('\d+', name))))
            MSname = str(anchor1) + '_' + str(anchor2)
            if self.parameter.method == 1:
                next_frame = 1
            else:
                next_frame = self.__snapshots(name)
            next_script = next_frame             
            launch = False
            if next_frame == 1:
                beg_number = 1
                end_number = self.parameter.trajPerLaunch + 1
            elif self.parameter.restart == True:
                beg_number = next_frame - self.parameter.trajPerLaunch
                end_number = next_frame
            else:
                beg_number = next_frame
                end_number = next_frame + self.parameter.trajPerLaunch
            for i in range(beg_number, end_number):
                traj_path = self.parameter.crdPath +'/' + str(anchor1) + '_' + str(anchor2) + '/' + str(self.parameter.iteration) + '/' + str(i)
                if os.path.isfile(traj_path + '/submit'):
                    self.jobs.submit(traj_path + '/submit')
                    launch = True
            if launch == True:
                print(str(datetime.now()))
                print("Short trajectories started from milestone {}...".format(name))
            else:
                print('No new trajectories started from {} this iteration'.format(name))
        
        log("{} short trajectories started from each milestone.".format(self.parameter.trajPerLaunch))
        for name in MS_list:
            [anchor1, anchor2] = list(map(int,(re.findall('\d+',name))))
            MSname = str(anchor1) + '_' + str(anchor2)
            MSpath = self.parameter.crdPath + '/' + MSname
            if self.parameter.restart == True:
                if next_frame == 1:
                    traj_frame = 1
                else:
                    traj_frame = self.__snapshots(MSname)
                    traj_frame = traj_frame - self.parameter.trajPerLaunch
                self.restore_scripts(MSpath + '/' + str(self.parameter.iteration), traj_frame)
        self.check()
        count = 0
        new_milestones = []
        for name in MS_list:
            [anchor1, anchor2] = list(map(int,(re.findall('\d+', name))))
            MSname = str(anchor1) + '_' + str(anchor2)
            MSpath = self.parameter.crdPath +  '/' + MSname
            last_frame = self.__snapshots(name)
            #if self.parameter.restart == True:
                #if self.parameter.method == 1:
                    #traj_frame = 1
                #else:
                    #traj_frame = self.__snapshots(MSname)
                    #traj_frame = traj_frame - self.parameter.trajPerLaunch
                #self.restore_scripts(MSpath + '/' + str(self.parameter.iteration), traj_frame)
            for j in range(self.parameter.nframe):
                path = MSpath + "/" + str(self.parameter.iteration) + "/" + str(j+1)
                if not os.path.exists(path):
                    continue
                time, final_ms = milestones(self.parameter).get_final_ms(path)
                count += 1
            if self.parameter.ignorNewMS == False:
                skip_compute, new_milestones = self.find_new_milestones(last_frame, MSpath + '/' + str(self.parameter.iteration), new_milestones)            
            else:
                skip_compute = False
        #if self.parameter.method == 1:
            for item in new_milestones:
                if item[1] >= self.parameter.new_ms_trajs or self.parameter.iteration <= self.parameter.new_ms_iterations:
                    if self.paramter.method == 1:
                        self.parameter.finished_constain.add(item[0])
                    self.parameter.MS_list.add(item[0])
                    [anchor1, anchor2] = list(map(int,(re.findall('\d+', item[0]))))
                    create_folder(self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2))
                    skip_compute = True
        return count / len(MS_list) + self.parameter.startTraj, skip_compute, new_milestones
    
    def restore_scripts(self, path, traj_frame):
        import os
        from fileinput import FileInput
        for i in range(self.parameter.trajPerLaunch):
            current_path = path + '/' + str(i + traj_frame)
            print(current_path)
            if os.path.isfile(current_path + '/submit_done'):
                os.rename(current_path + '/submit_done', current_path + '/submit')
            if os.path.isfile(current_path + '/submit'):
                with FileInput(files=current_path + '/submit', inplace=True) as r:
                    for line in r:
                        line = line.strip()
                        if line.startswith('##') and 'source' not in line:
                            line = line.replace('##', '')
                        print(line)
        
    def find_new_milestones(self, last_frame, MSpath, new_milestones):                
        import pandas as pd
        skip_compute = False
        if self.parameter.new_ms_iterations < self.parameter.iteration and self.parameter.new_ms_iterations != 0:
            self.parameter.ignorNewMS = True
            return False, new_milestones
        for i in range(1, last_frame):
            end_info = MSpath + '/' + str(i) + '/end.txt'
            if os.path.isfile(end_info):
                end = pd.read_csv(end_info, header=None, delimiter=r'\s+').values.tolist()[0]
            else:
                continue
            end_ms = 'MS' + str(end[0]) + '_' + str(end[1])
            if end_ms not in self.parameter.MS_list:
                new = True
                for item in new_milestones:
                    if end_ms in item[0]:
                        new = False
                        item[1] += 1
                if new == True:
                    new_milestones.append([end_ms, 1])
            #checking against user input
        '''
        if self.parameter.new_ms_trajs != 0 and len(new_milestones) > 0:
            for ms in new_milestones:
                if ms[1] >= self.parameter.new_ms_trajs:
                    skip_compute = True
                #maybe add an else delete here
            if skip_compute == False:
                self.parameter.ignorNewMS = True
        if self.parameter.new_ms_iterations != 0 and len(new_milestones) > 0:
            skip_compute = True
        if skip_compute == False:
            self.parameter.ignorNewMS = True
        '''
        return skip_compute, new_milestones
    
    
    def add_to_trajPool(self, dest, orig, path):
        '''add trajectory info to pool'''
        try:
            getattr(globals()[dest], orig).append(path)
        except:
            setattr(globals()[dest], orig, [path])

    def initialize_trajPool(self):
        '''initialize pool'''
        for ms in self.parameter.MS_list:
            globals()[ms] = trajPool(ms)

#    def fail_traj_restart(self, path):
#        import subprocess
#        state_path = path + '/' + self.parameter.outputname + '.colvars.state'
#        time, final_ms = milestones(self.parameter).read_state(state_path)
#        if time > 0.8 * self.parameter.freeTraj_walltime:
#            return
#        subprocess.run([self.parameter.jobsubmit, path + '/submit'])
        

if __name__ == '__main__':
#    MS2_3 = trajPool('2_3')
#    traj.add_to_trajPool(traj,'MS2_3', '1_2', 'path')
    from parameters import *
    from run import *
    from colvar import *


    new = parameters()
    new.initialize()
    jobs = run(new)
    new.iteration = 1
    new.traj_per_script = [2,2]
    new.trajPerLaunch = 6
    new.AnchorNum = 2
    new.MS_list = ['MS1_2']
    #traj(new, jobs).launch()
    traj = traj(new, jobs)
    traj.restore_scripts(new.crdPath + '/1_2/1', 1)
