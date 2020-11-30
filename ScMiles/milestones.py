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
        '''
        Description: This function is used to find the next non-existent folder
            So if we have the folders 
            crd/1_2/1, crd/1_2/2, crd/1_2/3
            This would return 4 since it is the first one to not exist.
        Arguments: struPath, the path that we are looking at. So for the example in the 
            description, this would be 'crd/1_2'
        Returns: next_frame, which is the next non-existent folder (4 in the example)
        '''
        next_frame = 1
        while True:
            pdbPath = struPath + '/' + str(next_frame) 
            if os.path.exists(pdbPath):
                next_frame += 1
            else:
                return next_frame
    
    def initialize(self, status=0):
        '''
        Description: This is used to initialize our MS_list
            We have the options of using: 
            - status=1 to just read in the milestones from crd and adds them to both MS_list and finished_constain (so we assume they have already been sampled on)
              this is mostly used when restarting the simulation
            - milestone_search = 0: traverse, Voronoi Cell
            - milestone_search = 1: seek, Voronoi Cell
            - milestone_search = 2: traverse, Grid
            - milestone_search = 3: seek, Grid
        Arguments: status: if status = 1 we just use the read_milestone_folder to initialize MS_list
           if status=0, we use the value for milestone_search to dictate how we inialize MS_list
        Returns: MS_list, which is stored in our main parameter instance
        '''
        MS_list = set()
        if status == 1:
            if self.parameter.milestone_search == 3:
                self.parameter.milestone_search = 2
            MS_list = self.read_milestone_folder()
            self.parameter.finished_constain = MS_list.copy()
            return MS_list
        elif self.parameter.milestone_search == 0:
            for i in range(1, self.parameter.AnchorNum):
                MS_list.add('MS' + str(i) + '_' + str(i + 1))
            if self.parameter.pbc:
                MS_list.add('MS' + str(self.parameter.pbc[0]) + '_' + str(self.parameter.pbc[1]))
            return MS_list
        elif self.parameter.milestone_search == 1 or self.parameter.milestone_search == 3:
            self.parameter.iteration = 0
            possible_ms = None
            if self.parameter.milestone_search == 3:
                self.parameter.milestone_search = 2
                possible_ms = self.grid_ms()
            while True:    
            #   free runs from each anchors, markdown once it reaches another cell (i.e. closer to another anchor ).
                MS_list = self.read_milestone_folder()
                if MS_list: 
                    if network_check(self.parameter, MS_list=MS_list) == True:
                        break
                self.seek_milestones(possible_ms)
            # check if reactant and product are connected
            # read folders to get the milestones list 
            return MS_list
        else: #grid
            MS_list = self.grid_ms(ms_list_form = 'yes')
            return MS_list

    def grid_ms(self, ms_list_form = None):
        '''
        Description: If we are using grid, this finds all of the possible miletones based on the anchors given
            This is any two anchors that exist that are adjacent to each other. 
        Arguments: ms_list_form
            If ms_list_form = 'yes', we store the values as a set with the form of our MS_list in parameters: set('MS1_2','MS2_3',...'MS11_12')
            If ms_list_form is omitted, we store in a  list of lists that just holds the values [[1,2],[2,3],...[11,12]]
        Returns: MS_list in either form
        '''
        '''
        #print(self.parameter.anchors)
        max_and_min = []
        if self.parameter.pbc:
            for i in range(len(self.parameter.anchors[0])):
                min_value = self.parameter.anchors[0][0]
                max_value = self.parameter.anchors[0][0]
                for j in range(len(self.parameter.anchors)):
                    if self.parameter.anchors[j][i] < min_value:
                        min_value = self.parameter.anchors[j][i]
                    elif self.parameter.anchors[j][i] > max_value:
                        max_value = self.parameter.anchors[j][i]
                if abs(min_value) + abs(max_value) == 360:
                    max_and_min.append([min_value, max_value])
                else:
                    max_and_min.append(['NA','NA'])
        '''
        if ms_list_form == 'yes':
            MS_list = set()
        else:
            MS_list = []   
        for i in range(len(self.parameter.anchors)):
            for j in range(i, len(self.parameter.anchors)):
                length_count = 0
                neighbor_count = 0
                for k in range(len(self.parameter.anchors[i])):
                    if self.parameter.anchors[i][k] == self.parameter.anchors[j][k]:
                        length_count += 1
                    elif self.parameter.anchors[i][k] - self.parameter.anchors[j][k] in [self.parameter.deltas[k],-self.parameter.deltas[k]]:
                        neighbor_count += 1
                    elif self.parameter.pbc and abs(self.parameter.anchors[i][k]) + abs(self.parameter.anchors[j][k]) == 360:
                        length_count += 1
                    elif self.parameter.pbc and self.parameter.anchors[i][k] + self.parameter.anchors[j][k] == self.parameter.deltas[k]:
                        neighbor_count += 1
                if neighbor_count == 1 and length_count == len(self.parameter.deltas) - 1:
                    if list_form:
                        MS_list.add('MS' + str(i+1) + '_' + str(j+1))
                    else:
                        MS_list.append([i+1,j+1])
        return MS_list

    def get_initial_ms(self, path):
        '''
        Description: This finds the initial MS (the folder that we are in) and prints it to the file start.txt
        Arguments: path: this is the path we are looking at. So if we are in crd/1_2/1/1, this will find that
            our initial_ms is 1_2
        Returns: initial_ms: this is the ms that the trajectory started from, whihc is just the folder we are in
        '''
        path_split = path.split("/")
        initial_ms = list(map(int,(re.findall('\d+', path_split[-3]))))
        with open(path + '/start.txt', 'w+') as f1:
            print(initial_ms[0], initial_ms[1], file=f1)    
        return initial_ms
        
    def get_final_ms(self, path, anchor=None):
        '''
        Description: After launching trajectories, this looks for a stop.colvars.state file
            If the file exists, it is then read in and we find the ending milestone
            If we return -1, [0,0], this means that for some reason there is no ending milestone and
            end.txt and lifetime.txt are not printed. Then at this point the trajectory is not considered
            in the final computations.
        Arguments:
            path: This is the path that we are looking at, and where we will extract the stop.colvars.state
                information from
            anchor: This is
        '''
        
        state = path + "/stop.colvars.state"
        if not os.path.isfile(state):
            print('1')
            return -1, [0, 0]
        
        final_ms = [0, 0]
        
        # read state file generated at termination point 
        # smallest rmsd indicates new cell #
        RMSDs, lifetime = self.read_state(state)
        final_ms[0] = RMSDs.index(sorted(RMSDs)[0]) + 1
        
        if self.parameter.milestone_search == 2:
            if anchor:
                colvar_number = len(self.parameter.anchors[0])
                anchor_values = self.parameter.anchors[anchor-1]
                count = 0
                values = []
                for i in range(colvar_number):
                    if i == 0:
                        values.append(1)
                    else:
                        values.append(2)
                end_value = None
                for i in range(len(RMSDs)-len(values) + 1):
                    if end_value:
                        break
                    count = 0
                    for j in range(len(values)):
                        if int(RMSDs[i+j]) == values[j]:
                            count += 1
                        else:
                            break
                        if count == colvar_number:
                            end_value = i + colvar_number #this is the equation number 
                total_equations = 3*colvar_number
                #end_value = 18
                #final = total_equations / end_value
                new_end = int(end_value / colvar_number)
                possible = []
                new_anchor = []
                for i in range(colvar_number):
                    a = [anchor_values[i] - self.parameter.deltas[i], anchor_values[i] + self.parameter.deltas[i]]
                    possible.append([min(a), i])
                    possible.append([max(a), i])
                    new_anchor.append(0)
                #print(possible)
                neighbor = possible[new_end-1][0]
                neighbor_name = possible[new_end-1][1]
                for i in range(colvar_number):
                    if neighbor_name == i:
                        new_anchor[i] = neighbor
                    else:
                        new_anchor[i] = anchor_values[i]
                final_ms = [anchor, 0]
                #print(anchor_values)
                #print(new_anchor)
                for i in range(self.parameter.AnchorNum):
                    for j in range(colvar_number):
                        if self.parameter.anchors[i][j] == new_anchor[j]:
                            if j == len(new_anchor) - 1:
                                final_ms[1]= i+1
                        else:
                            break
                #print(final_ms)
                final_ms.sort()
                
            else:
                for i in range(len(self.parameter.anchors)):
                    for j in range(len(self.parameter.anchors[0])):
                        delta_values = self.parameter.deltas[j]/2
                        values = [self.parameter.anchors]
                
                colvar_number = len(self.parameter.anchors[0])
                start_ms = pd.read_csv(path + '/start.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                start_anchors = [self.parameter.anchors[start_ms[0] - 1], self.parameter.anchors[start_ms[1] - 1]]
                count = 0
                values = []
                for i in range(colvar_number):
                    if i == 0:
                        values.append(1)
                    else:
                        values.append(2)
                end_value = None
                for i in range(len(RMSDs)-len(values) + 1):
                    if end_value:
                        break
                    count = 0
                    for j in range(len(values)):
                        if int(RMSDs[i+j]) == values[j]:
                            count += 1
                        else:
                            break
                        if count == colvar_number:
                            end_value = i + colvar_number #this is the equation number            
                length = []
                #print(end_value)
                for i in range(len(start_anchors[0])):
                    if start_anchors[0][i] != start_anchors[1][i]:
                        if start_anchors[0][i] > start_anchors[1][i]:
                            upper_ms = 0
                            lower_ms = 1
                        else:
                            upper_ms = 1
                            lower_ms = 0
                        neighbor_values = [start_anchors[0][i], start_anchors[1][i]]
                        neighbor_values.sort()
                        neighbor_index = i
                        length.append("n")
                    else:
                        length.append(start_anchors[0][i])
                new_neighbor = None  
                # first are our 2 parallel planes, where our lengths are the same and our neighbor changes
                if end_value == colvar_number*1:
                    new_neighbor = neighbor_values[1] + self.parameter.deltas[neighbor_index]
                    final_ms[0] = start_ms[upper_ms]
                elif end_value == colvar_number*2:
                    new_neighbor = neighbor_values[0] - self.parameter.deltas[neighbor_index]
                    final_ms[0] = start_ms[lower_ms]
                if new_neighbor is not None:
                    length[neighbor_index] = new_neighbor
                #now we have our perpendicular planes in each direction.
                else:
                    total_statements = (2 + (len(self.parameter.anchors[0])-1)*4)*colvar_number
                    remaining = total_statements - 2*colvar_number #24
                    length_cut = remaining / (colvar_number-1) #so 24/2 = 12
                    mod_end = (end_value - colvar_number*2) % length_cut
                    end_var_index = (end_value - 2*colvar_number - 1) // length_cut
                    count = -1
                    for i in range(len(length)):
                        if length[i] == 'n':
                            continue
                        else:
                            count += 1
                            if count == end_var_index:
                                changing_length = i
                    #print(path)
                    if mod_end == colvar_number*1: #first option after parallel
                        new_neighbor = length[changing_length] + self.parameter.deltas[changing_length]
                        new_length = neighbor_values[1] 
                        final_ms[0] = start_ms[upper_ms]
                    if mod_end == colvar_number*2:
                        new_neighbor = length[changing_length] + self.parameter.deltas[changing_length]
                        new_length = neighbor_values[0] 
                        final_ms[0] = start_ms[lower_ms]
                    if mod_end == colvar_number*3:
                        new_neighbor = length[changing_length] - self.parameter.deltas[changing_length]
                        new_length = neighbor_values[1]
                        final_ms[0] = start_ms[upper_ms]
                    if mod_end == 0:
                        new_neighbor = length[changing_length] - self.parameter.deltas[changing_length]
                        new_length = neighbor_values[0]
                        final_ms[0] = start_ms[lower_ms]
                    
                    for i in range(len(length)):
                        if length[i] == 'n':
                            length[i] = new_length
                        if i == changing_length:
                            length[i] = new_neighbor
                    
                for i in range(self.parameter.AnchorNum):
                    for j in range(len(length)):
                        if self.parameter.anchors[i][j] == length[j]:
                            if j == len(length) - 1:
                                final_ms[1]= i+1
                        else:
                            break
                final_ms.sort()

            
        else:
            if self.parameter.pbc:
                if final_ms[0] == self.parameter.AnchorNum or final_ms[0] == 1:
                    # open traj file to read the very last output
                    # smallest rmsd indicates previous cell #
                    traj = path + "/" + self.parameter.outputname + ".colvars.traj"
                    firstRMSD = self.parameter.colvarsNum + 1
                    try:
                        RMSDs_prev = pd.read_fwf(traj, widths=self.parameter.trajWidths).values[-1,firstRMSD:].astype(np.float64).tolist()
                    except:
                        return -1, [0, 0]
                    final_ms[1] = RMSDs_prev.index(sorted(RMSDs_prev)[0]) + 1 
                else:
                    final_ms[1] = RMSDs.index(sorted(RMSDs)[1]) + 1 
            else:
                # use the second min value for previous cell #
                final_ms[1] = RMSDs.index(sorted(RMSDs)[1]) + 1 
            
            if 'seek' not in path:
                try:
                    start_ms = pd.read_csv(path + '/start.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                except:
                    return -1, [0, 0]
                if start_ms[0] not in final_ms and start_ms[1] not in final_ms:
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
        return lifetime, final_ms            
        
    def seek_milestones(self, possible_ms = None):
        '''
        This is the main workflow for the 'seek' option
        This function is used when milestone_search = 1 for Voronoi Cell
        or milestone_search = 3 for Grid
        '''
        from shutil import copy
        from run import run
        milestones = set()
        launch = []
        for i in range(self.parameter.AnchorNum):
            launch.append(False)
        #prepare scripts
        if self.parameter.milestone_search  == 1:
            colvar(self.parameter, free='yes', initial='yes').generate()  
        for an in range(1, self.parameter.AnchorNum + 1):
            initialNum = self.__get_next_frame_num(self.parameter.seekPath + '/structure' + str(an))
            if self.parameter.restart == True and initialNum > 1:
                self.edit_submit_scripts(self.parameter.seekPath + '/structure' + str(an),None,initialNum)
            else:
                next_script = initialNum
                for i in range(self.parameter.initial_traj):
                    submit = False
                    if i + initialNum == next_script:
                        next_script = next_script + self.parameter.traj_per_script[0]
                        submit = True
                    if self.parameter.milestone_search != 1:
                        colvar(self.parameter, free='yes',initial='yes', anchor1=an).generate()
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
                if len(finished) == self.parameter.AnchorNum:
                    break
                time.sleep(60)
        
        for i in range(1, self.parameter.AnchorNum + 1):
            curt_frame = self.__get_next_frame_num(self.parameter.seekPath + '/structure' + str(i))
            for traj in range(1, curt_frame):
                path = self.parameter.seekPath + '/structure' + str(i) + '/' + str(traj)
                #print(path)
                if not os.path.exists(path):
                    continue
                if os.path.isfile(path + '/end.txt'):
                    continue
                if self.parameter.milestone_search == 1:
                    timetmp, final_ms = self.get_final_ms(path)
                else:
                    timetmp, final_ms = self.get_final_ms(path, anchor=i)
                #print(final_ms)
                if final_ms == [0, 0]:
                    continue
                name = 'MS' + str(final_ms[0]) + '_' + str(final_ms[1])
                if name in milestones:
                    continue
                ms_path = self.parameter.crdPath + '/' + str(final_ms[0]) + '_' + str(final_ms[1])
                if self.parameter.milestone_search != 1:
                    if final_ms[0] == 0 or final_ms[1] == 0 or final_ms not in possible_ms:
                        continue
                if os.path.exists(ms_path):
                    continue          
                elif not os.path.isfile(path + '/' + self.parameter.outputname + '.restart.coor'):
                    continue
                if self.parameter.dist_cut != 0:
                    keep_ms = self.__check_distance(name)
                    if keep_ms == False:
                        continue
                os.makedirs(ms_path)
                copy(path + '/' + self.parameter.outputname + '.restart.coor', 
                    ms_path + '/seek.ms.pdb')
                #print(self.parameter.namd_conf)
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
                        
    def __check_distance(self, milestone):
        #This only works for Voronoi cell for now.
        from math import sqrt
        anchors = []
        total = 0
        dist = 0
        if self.parameter.milestone_search == 1: #Just in case we are adding this to grid as well
            for i in milestone:
                anchors.append(self.parameter.anchors[i-1])
            #calculate distance
             #this is for each anchor in the milestone
            total = 0
            for i in range(len(anchors[0])): #this is for each colvar
                total += (anchors[0][i] - anchors[1][i])**2
            dist = sqrt(total)
            if dist > self.parameter.dist_cut:
                return False
            else:
                return True
    
    def read_milestone_folder(self):
        '''
        This function reads in the names of the folders in crd
        Returns a set() with every foldername in crd 
        '''
        MS_list = set()
        for i in range(1, self.parameter.AnchorNum):
            for j in range(i, self.parameter.AnchorNum + 1):
                name = str(i) + '_' + str(j)
                if os.path.exists(self.parameter.crdPath + '/' + name):
                    MS_list.add('MS' + name)
        return MS_list

    def read_state(self, path):
        '''
        Opens the stop.colvars.state file and reads in the values.
        Saves the step value as lifetime and the other values as the variable values
        We get rid of the first values that are the ending values of our colvars,
        and just keep the values for rmsd (if we are using Voronoi) or
        step functions (if we are using grid).
        '''
        file = open(path, 'r').read()
        time = int(re.findall(r"[-+]?\d*\.\d+|\d+", file)[0])
        numbers = re.findall(r"[-+]?\ *[0-9]+\.?[0-9]*(?:[eE]\ *[-+]?\ *[0-9]+)", file)
        numbers.pop(0)
        numbers = numbers[self.parameter.colvarsNum:]
        values = [abs(float(x)) for x in numbers]
        return values, time    
    
    def edit_submit_scripts(self, path, restart_scripts,frame_number):
        #1 is done, 0 needs to restart
        import os
        from fileinput import FileInput

        for i in range(1, frame_number + 1):
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
    #from run import *
    from colvar import *
    
    new = parameters()
    new.initialize()
    #jobs = run(new)
    new.iteration = 1
    new.pbc = True
    new.traj_per_script = [2,2]
    new.initial_traj = 6
    new.milestone_search = 2
    new.deltas = [30,30,30]
    new.AnchorNum = 42*3
    new.MS_list = milestones(new).initialize(status=0)
    #milestones(new).get_final_ms(new.crdPath + '15_16/1/14')
    print(new.MS_list)
    #milestones(new).get_final_ms_grid(new.ScMilesPath)
