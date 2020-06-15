#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:35:54 2019

@author: Wei Wei

This subroutine stores the milestone information.
It initializes the milestone list, also contains the function that provides initial and final milestone.

"""

import os, re, time
import pandas as pd
import numpy as np
from run import *
from colvar import *
from log import log
from parameters import *
from network_check import *


class milestones: 

    def __init__(self, parameter) -> None:
        self.parameter = parameter
        
    def __enter__(self):
        return self
               
    def __exit__(self, exc_type, exc_value, traceback):
        return 
            
    def __repr__(self) -> str:
        return ('miletstones details.'
                .format(self.anchor_orig, self.anchor_dest, self.lifetime))

    def __get_next_frame_num(self, struPath):
        next_frame = 1
        while True:
            pdbPath = struPath + '/' + str(next_frame) 
            if os.path.exists(pdbPath):
                next_frame += 1
            else:
                return next_frame
    
    def initialize(self, status=0):
        self.parameter.iteration = 0
        MS_list = set()
        if status == 1:
            MS_list = self.read_milestone_folder()
            self.parameter.finished_constain = MS_list.copy()
            return MS_list
        elif self.parameter.milestone_search == 0:
            for i in range(1, self.parameter.AnchorNum):
                MS_list.add('MS' + str(i) + '_' + str(i + 1))
            if self.parameter.pbc:
                MS_list.add('MS' + str(self.parameter.pbc[0]) + '_' + str(self.parameter.pbc[1]))
            return MS_list
        else:
            while True:    
            #   free runs from each anchors, markdown once it reaches another cell (i.e. closer to another anchor ).
                MS_list = self.read_milestone_folder()
                if MS_list: 
                    if network_check(self.parameter, MS_list=MS_list) == True:
                        break
                self.seek_milestones()
            # check if reactant and product are connected
            # read folders to get the milestones list 
            return MS_list

    def get_initial_ms(self, path):
        path_split = path.split("/")
        initial_ms = list(map(int,(re.findall('\d+', path_split[-3]))))
        with open(path + '/start.txt', 'w+') as f1:
            print(initial_ms[0], initial_ms[1], file=f1)    
        return initial_ms
    
    def get_final_ms(self, path):
        state = path + "/stop.colvars.state"
        if not os.path.isfile(state):
            print("HERE")
            return -1, [0, 0]
        final_ms = [0, 0]
        # read state file generated at termination point 
        # smallest rmsd indicates new cell #
        RMSDs, lifetime = self.read_state(state)
        final_ms[0] = RMSDs.index(sorted(RMSDs)[0]) + 1
        
        if self.parameter.pbc:
            if final_ms[0] == self.parameter.AnchorNum or final_ms[0] == 1:
                # open traj file to read the very last output
                # smallest rmsd indicates previous cell #
                traj = path + "/" + self.parameter.outputname + ".colvars.traj"
                firstRMSD = self.parameter.colvarsNum + 1
                try:
                    RMSDs_prev = pd.read_fwf(traj, widths=self.parameter.trajWidths).values[-1,firstRMSD:].astype(np.float64).tolist()
                except:
                    print('HERE2')
                    return -1, [0, 0]
                final_ms[1] = RMSDs_prev.index(sorted(RMSDs_prev)[0]) + 1 
            else:
                final_ms[1] = RMSDs.index(sorted(RMSDs)[1]) + 1 
        else:
            # use the second min value for previous cell #
            final_ms[1] = RMSDs.index(sorted(RMSDs)[1]) + 1 
        
        if self.parameter.iteration >= 1:
            try:
                start_ms = pd.read_csv(path + '/start.txt', header=None, delimiter=r'\s+').values.tolist()[0]
            except:
                print("here 3")
                return -1, [0, 0]
            if start_ms[0] not in final_ms and start_ms[1] not in final_ms:
                print("Here 4")
                return -1, [0, 0]
        final_ms.sort()
        final_info = path + '/end.txt'
        if not os.path.isfile(final_info):
            with open(final_info, 'w+') as f1:
                print(final_ms[0], final_ms[1], file=f1)    
        time_info = path + '/lifetime.txt'
        lifetime *= self.parameter.timeFactor
        if not os.path.isfile(time_info):
            with open(time_info, 'w+') as f1:
                print(lifetime, file=f1)
        print("FINISHED")    
        return lifetime, final_ms            
        
    def seek_milestones(self):
        from shutil import copy
        milestones = set()
        launch = []
        for i in range(self.parameter.AnchorNum):
            launch.append(False)
        #prepare scripts
        colvar(self.parameter, free='yes', initial='yes').generate()  
        for an in range(1, self.parameter.AnchorNum + 1):
            initialNum = self.__get_next_frame_num(self.parameter.seekPath + '/structure' + str(an))
            if self.parameter.restart == True and initialNum > 1:
                self.__restart_seek(an, initialNum)
            else:
                next_script = initialNum
                for i in range(self.parameter.initial_traj):
                    submit = False
                    if i + initialNum == next_script:
                        next_script = next_script + self.parameter.traj_per_script[0]
                        submit = True
                    create_folder(self.parameter.seekPath + '/structure' + str(an) + '/' + str(i + initialNum))
                    run(self.parameter).prepare_trajectory(a1=an, a2=999, initial='yes', initialNum=i+initialNum, script=submit)
        #submit jobs
        for an in range(1,self.parameter.AnchorNum + 1):
            last_frame = self.__get_next_frame_num(self.parameter.seekPath + '/structure' + str(an))
            for i in range(last_frame - self.parameter.initial_traj, last_frame):
                if os.path.exists(self.parameter.seekPath + '/structure' + str(an) + '/' + str(i) + '/submit'):
                    run(self.parameter).submit(self.parameter.seekPath + '/structure' + str(an) + '/' + str(i) + '/submit')
                    launch[an - 1] = True
        if self.parameter.restart == False:        
            log("{} trajectories started from each anchor, run for {} ps.".format(self.parameter.initial_traj, self.parameter.initialTime))
        elif True in launch:
            log("Remaining trajectories have been launched")
        else:
            print('No new trajectories launched for this iteration. Starting next iteration')        

        if True in launch:
            finished = []
            while True:
                for an in range(1, self.parameter.AnchorNum + 1):
                    MSname = 'a' + str(an)
                   #if launch[an-1] == False:
                        #finished.append(MSname)
                    if not run(self.parameter).check(MSname=MSname):
                        continue
                    elif MSname in finished:
                        continue
                    else:
                        finished.append(MSname)
                    print(finished)
                if len(finished) == self.parameter.AnchorNum:
                    break
                time.sleep(60)
        
        for i in range(1, self.parameter.AnchorNum + 1):
            curt_frame = self.__get_next_frame_num(self.parameter.seekPath + '/structure' + str(i))
            for traj in range(1, curt_frame):
                path = self.parameter.seekPath + '/structure' + str(i) + '/' + str(traj)
                if not os.path.exists(path):
                    continue
                if os.path.isfile(path + '/end.txt'):
                    continue
                timetmp, final_ms = self.get_final_ms(path)
                print(final_ms)
                if final_ms == [0, 0]:
                    continue
                name = 'MS' + str(final_ms[0]) + '_' + str(final_ms[1])
                if name in milestones:
                    continue
                ms_path = self.parameter.crdPath + '/' + str(final_ms[0]) + '_' + str(final_ms[1])
                if os.path.exists(ms_path):
                    continue          
                elif not os.path.isfile(path + '/' + self.parameter.outputname + '.restart.coor'):
                    continue
                else:
                    os.makedirs(ms_path)
                    copy(path + '/' + self.parameter.outputname + '.restart.coor', 
                         ms_path + '/seek.ms.pdb')
                    if self.parameter.namd_conf == True:
                        copy(path + '/' + self.parameter.outputname + '.xst', ms_path + '/sample.xsc')
                    milestones.add(name)
        for i in range(1,self.parameter.AnchorNum + 1):
            current_frame = self.__get_next_frame_num(self.parameter.seekPath + '/structure' + str(i))     
            traj_frame = current_frame - self.parameter.initial_traj
            self.__restore_scripts(self.parameter.seekPath + '/structure' + str(i), traj_frame)

        log("{} milestones have been identified.".format(len(milestones)))  
        self.parameter.restart = False
        return milestones  

    def __restore_scripts(self, path, traj_frame):
        import os
        from fileinput import FileInput
        for i in range(self.parameter.initial_traj):
            current_path = path + '/' + str(i + traj_frame)
            if os.path.isfile(current_path + '/submit_done'):
                os.rename(current_path + '/submit_done', current_path + '/submit')
            if os.path.isfile(current_path + '/submit'):
                with FileInput(files=current_path + '/submit', inplace=True) as r:
                    for line in r:
                        line = line.strip()
                        if line.startswith('##') and 'source' not in line:
                            line = line.replace('##', '')
                        print(line)
                        
                        
    def read_milestone_folder(self):
        MS_list = set()
        for i in range(1, self.parameter.AnchorNum):
            for j in range(i, self.parameter.AnchorNum + 1):
                name = str(i) + '_' + str(j)
                if os.path.exists(self.parameter.crdPath + '/' + name):
                    MS_list.add('MS' + name)
        return MS_list

    def read_state(self, path):
        file = open(path, 'r').read()
        time = int(re.findall(r"[-+]?\d*\.\d+|\d+", file)[0])
        numbers = re.findall(r"[-+]?\ *[0-9]+\.?[0-9]*(?:[eE]\ *[-+]?\ *[0-9]+)", file)
        numbers.pop(0)
        numbers = numbers[self.parameter.colvarsNum:]
        rmsd = [abs(float(x)) for x in numbers]
        return rmsd, time    
 
    def __restart_seek(self, anchor, frame):
        submit_status = []
        anchor_path = self.parameter.seekPath + '/structure' + str(anchor)
        for i in range(frame - self.parameter.initial_traj, frame):
            if os.path.exists(anchor_path + '/' + str(i) + '/stop.colvars.state'):
                submit_status.append(1)
            else:
                submit_status.append(0)
        self.edit_submit_scripts(anchor_path, submit_status)
    
    def edit_submit_scripts(self, path, restart_scripts):
        #1 is done, 0 needs to restart
        import os
        from fileinput import FileInput

        for i in range(1, self.parameter.trajPerLaunch + 1):
            next_line = False
            run = False
            current_path = path + '/' + str(i) + '/submit'
            if os.path.exists(current_path):
                with FileInput(files=current_path, inplace=True) as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('##'):
                            print(line)
                            continue
                        if next_line == True:
                            line = '##' + line
                            next_line = False
                        elif 'cd' in line and ('crd' not in line or ('crd' in line and self.parameter.traj_per_script[0] == 1)):
                            info = line.split('/')
                            traj_num = int(info[-1]) - 1
                            if os.path.isfile(path + '/' + str(i) + '/stop.colvars.state'):
                                line = '##' + line
                                next_line = True
                            else:
                                run = True
                        print(line)
                if run == False:
                    if not os.path.exists(path + '/' + str(i) + '/submit_done'):
                        os.rename(current_path, path + '/' + str(i) + '/submit_done')
                        
                        
if __name__ == '__main__':
    from parameters import *
    from run import *
    from colvar import *
    
    new = parameters()
    new.initialize()
    jobs = run(new)
    new.iteration = 1
    new.traj_per_script = [2,2]
    new.initial_traj = 6
    new.AnchorNum = 2
    colvar(new, free='yes', initial='yes').generate()
    milestones(new).seek_milestones()
