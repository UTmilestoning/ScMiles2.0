#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:20:21 2019

@author: Wei Wei

Running free trajectories.

"""

from math import isnan
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
            #self.parameter.MS_list.remove(self.__dest__)
            log("Wrong distribution on {}. Maybe due to no trajectory reached or wrong flux on the adjacent milestone"
                .format(self.__dest__))
        for orig_ms in orig_list:   
            distribution_prob[orig_ms] = distribution[orig_ms] / sum(distribution.values())
        return distribution_prob
        

class traj:
    def __init__(self, parameter):
        self.parameter = parameter
        #self.jobs = jobs
        self.pool = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return
    
    def __repr__(self) -> str:
        return ('Free trajectories.')    

    def __snapshots(self, MSname, ignore_method=False):
        '''
        param MSname:
            milestone name
        return:
            next frame number
        example:
            250 free traj simulations for milestone MSname, return 251 as the next traj number for folder name
        '''
        next_frame = 1
        if self.parameter.method == 1 and self.parameter.restart == False and ignore_method == False:
            return next_frame
        lst = re.findall('\d+',MSname)
        name = lst[0] + '_' + lst[1]
        while True:
            path = self.parameter.crdPath + '/' + name + '/' + str(self.parameter.iteration) + '/' + str(next_frame) 
            if os.path.exists(path):
                next_frame += 1
            else:
                return next_frame
    
    def __copy(self, orig, new, restart=1, enhanced=False):
        '''copy restart files for next iteration'''
        create_folder(new)
        if self.parameter.software == 'namd': #different since namd has 3 files
            lst = ['.coor', '.vel', '.xsc']
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
                try:  
                    copyfile(oldName, newName)
                except:
                    pass
        else:
            if self.parameter.software == 'gromacs':
                newName = new + '/' + self.parameter.outputname + '.cpt'
                oldName = orig + '/' + 'state.cpt'
            else: #lammps
                newName = new + '/' + self.parameter.outputname + '.restart'
                oldName = orig + '/' + self.parameter.outputname + '.restart'
            with open(new + '/original.txt', 'w+') as f:
                print('This trajectory is from',file=f)
                print(orig, file=f)            
            if enhanced:
                with open(new + '/enhanced','w+') as f:
                    print('This trajectory is from enhancement', file=f)
            try:
                copyfile(oldName, newName)
            except:
                pass 


    def __iteration_prepare(self, path, final_ms, orig):
        '''add each trajectory to trajPool based on the ending milestone'''
        if final_ms == [0,0]:
            return
        if self.parameter.software == 'namd' and not os.path.isfile(path + '/' + self.parameter.outputname + '.restart.coor'):
            return
        if self.parameter.software in ('lammps','gromacs') and not os.path.isfile(path + '/' + self.parameter.outputname + '.colvars.traj'):
            return
        ms = str(final_ms[0]) + '_' + str(final_ms[1])
        MSname = 'MS' + ms
        if MSname not in self.parameter.MS_list:
            return
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
                if self.parameter.ignorNewMS == True:
                    final = 'MS' + str(final_ms[0]) + '_' + str(final_ms[1])
                    if final not in self.parameter.MS_list:
                        continue   
#                print(final_ms)
                self.__iteration_prepare(traj_path, final_ms, name)

        kij, index, flux = compute(self.parameter, partial_compute = True)
        
#        print(MS10_11)
#        print(MS10_11.count_traj('MS11_12'))
        ms_list = self.parameter.MS_list.copy()
        for MSname in ms_list:
            name = re.findall('\d+', MSname)[0] + '_' + re.findall('\d+', MSname)[1]
            path = self.parameter.crdPath + '/' + name + '/' + str(self.parameter.iteration) + '/'
            if self.parameter.restart == True:
                if os.path.exists(path + 'distribution'):
                    continue
            distribution = globals()[MSname].get_distribution(kij, index, flux)
            if distribution == 0:
                self.parameter.skip_MS.append(MSname)
                continue
            create_folder(path)
            path = path + 'distribution'
            with open(path, 'w+') as f:
                print(distribution, file=f)  
#            print(MSname)
#            print(distribution)
            total_orig, orig_list = globals()[MSname].count_ms() 
            if total_orig == 0:
                print("No traj reached {}.".format(MSname))
                # print("Exiting...")
                log("No traj reached {}.".format(MSname))
                log("Skip {} for next iteration.".format(MSname))
                self.parameter.skip_MS.append(MSname)
                continue
                # sys.exit()
            total_traj = 0
#            print(MSname, globals()[MSname].count_ms())
            for ms in distribution.keys():
                if isnan(distribution[ms]):
                    continue 
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
        path = self.parameter.crdPath +'/' + name + '/' + str(self.parameter.iteration) + '/'  
        create_folder(path)    
        continued = []
        for i, traj in enumerate(traj_lst):
            count = i + 1 + total_traj
            newPath = path + str(count)
            create_folder(newPath)
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
        if (len(traj_lst) - len(continued)) > 100:
            log("Warning: for iteration {}, trajectories from {} have been enhanced more than 100 times for ms {}".format(self.parameter.iteration, name, ms))
    
    def launch(self, joblist=None, lastLog=None):
        '''launch free trajectories'''
        print("Iteration # {}".format(self.parameter.iteration))
        log("Iteration # {}".format(self.parameter.iteration))
        
        if self.parameter.method == 1 and self.parameter.iteration > 1:
            if self.parameter.restart == False or (self.parameter.restart == True and "Iteration #" in lastLog):
                self.iteration_initialize()
                log('Iteration initialize complete')
        elif self.parameter.method == 0 and self.parameter.skip_compute == True or self.parameter.restart == True:
            launch(self.parameter, 'sample').launch_trajectories()

        self.parameter.skip_compute = False

        for ms in self.pool:
            self.pool[ms] = 0
        
        count = launch(self.parameter, 'free').launch_trajectories(lastLog = lastLog)
        return count / len(self.parameter.MS_list) + self.parameter.startTraj
    
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
    new.MS_list = ['MS1_2', 'MS3_9']
    traj(new, jobs).launch()
    #traj = traj(new, jobs)
    #traj.restore_scripts(new.crdPath + '/1_2/1', 1)
