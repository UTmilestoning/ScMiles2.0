#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:35:54 2019

@author: Wei Wei

This subroutine stores the milestone information.
It initializes the milestone list, also contains the function that provides initial and final milestone.

"""

import os, re, time
import itertools
import pandas as pd
import numpy as np
from run import *
from colvar import *
from log import log
from parameters import *
from network_check import *
from plumed import *
#from test_new_launch import *

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
        
    
    def initialize(self, status=0, lastLog=""):
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
        #seek = launch(self.parameter, 'seek')
        
        MS_list = set()
        if status == 1: #status set to 1 just has us skip this entire step. It just populates MS_list with the names of our folders in crd
            if self.parameter.milestone_search == 3:
                self.parameter.milestone_search = 2
            MS_list = self.read_milestone_folder()
            self.parameter.finished_constain = MS_list.copy()
        elif self.parameter.milestone_search == 0: #This is traverse, so each adjacent anchor makes a milestone, like 1_2,2_3,3_4...,11_12
            for i in range(1, self.parameter.AnchorNum):
                MS_list.add('MS' + str(i) + '_' + str(i + 1))
            if self.parameter.pbc: #if pbc, we add the milestone of the values in the user input pbc. In Alanine Dipeptide this is the milestone 1_12
                MS_list.add('MS' + str(self.parameter.pbc[0]) + '_' + str(self.parameter.pbc[1]))
        elif self.parameter.milestone_search == 1 or self.parameter.milestone_search == 3: 
            #This is for seek, 1 is for Voronoi Cell, 3 is for Grid
            self.parameter.iteration = 0
            possible_ms = None
            if self.parameter.milestone_search == 3: #After this point 2 and 3 are the same thing, so I decided to just make them the same 
                self.parameter.milestone_search = 2 
                self.parameter.all_grid = self.grid_ms()
            while True: #this keeps going until we have a network where the reactant and product are connected 
            #   free runs from each anchors, markdown once it reaches another cell (i.e. closer to another anchor ).
                MS_list = self.read_milestone_folder() #We read milestone folder at the beginning of a loop because if we are restarting our script, seek will be skipped over because a network already exists
                if MS_list:
                    if network_check(self.parameter, MS_list=MS_list) == True:
                        break
                launch(self.parameter, 'seek').launch_trajectories(possible_ms, lastLog) #This is where we launch seek trajectories
                self.parameter.restart = False
                lastLog=None
            # check if reactant and product are connected
            # read folders to get the milestones list 
            return MS_list
        else: #grid
            MS_list = self.grid_ms(ms_list_form = 'yes')
        for ms in MS_list:
            milestones = get_anchors(ms)
            milestones = sorted(milestones)
            create_folder(self.parameter.crdPath + '/' + str(milestones[0]) + '_' + str(milestones[1]))  
        if self.parameter.milestone_search == 2 and not self.parameter.all_grid:
            self.parameter.all_grid = self.grid_ms()
        return MS_list

    def grid_ms(self, ms_list_form = None):
        '''
        Description: If we are using grid, this finds all of the possible miletones based on the anchors given
            This is any two anchors that exist that are adjacent to each other. 
        Arguments: ms_list_form
            If ms_list_form = 'yes', we store the values as a set with the form of our MS_list in parameters: set('MS1_2','MS2_3',...'MS11_12')
            If ms_list_form is omitted, we store in a  list of lists that just holds the values [[1,2],[2,3],...[11,12]]
        Returns: MS_list in either form. Defaulys to a list (I realize that calling it ms_list_form is kind of confusing, I just named it that because
            our main set we use is called MS_list. So ms_list_form is referring to this, not the python data type of list.
        '''
        if ms_list_form == 'yes':
            MS_list = set()
        else:
            MS_list = dict() 
        values = []
        #colvars, names = colvar(self.parameter).get_colvars()
        for i in range(len(self.parameter.anchors)):
            for j in range(len(self.parameter.anchors)):
                if i == j:
                    continue
                length_count = 0
                neighbor_count = 0
                for k in range(len(self.parameter.anchors[i])): #This is the amount of colvars
                    if round(self.parameter.anchors[i][k],4) == round(self.parameter.anchors[j][k],4): #The rounding is more for gromacs, it may need to be updated
                        #If they are equal, this is (in my notation) considered a 'length', so there are two step functions, so it has to be between a certain value
                        '''         30
                                     '  
                                 l-------l So for these anchors, our length would be 30 degrees, and spans from 15 to 45 (and is our x variable)   
                           90 -  l   o   l Our neighbor is the middle point between the two differing values (60 + 90)/2 = 75
                                 l-------l There should be ONE neighbor, and the rest should be lengths.
                           60 -  l   o   l 
                                 l-------l
                                           Also, for this iteration of the script, I have updated the way I am having grid work. I am doing it more like a 
                                           network rather than having it build the step functions from scratch every time. The main thing that is kind of weird
                                           is that for the purpose of the step functions, I have the milestones 5_6 and 6_5 being different things. So like for
                                           the picture, lets say 5 is at 30,60 and 6 is at 90,60. So for the step functions 5_6, Crossing this milestone would be
                                           when our y > 75. But from the other side, if we are starting at 6 and going to 5, we cross the milestone when y < 75. So
                                           it basically just has to do with the directionality. For all other purposes, they are sorted and are both treated as 5_6,
                                           this is just for pulling the functions for colvar and plumed files.
                        '''
                        if self.parameter.names[k] in self.parameter.grid_ignore: #if the user wants to ignore it, we don't add it to our grid dictionary
                            continue
                        length_count += 1
                        values.append(('step(' + self.parameter.names[k] + '-' + str(self.parameter.anchors[i][k] - self.parameter.deltas[k]/2) + ')' +
                                           ' + step(' +  str(self.parameter.anchors[i][k] + self.parameter.deltas[k]/2) + '-' + self.parameter.names[k]) + ')')
                        #values.append([self.parameter.anchors[i][k] + self.parameter.deltas[k]/2, self.parameter.anchors[i][k] - self.parameter.deltas[k]/2])
                    elif round(self.parameter.anchors[i][k] - self.parameter.anchors[j][k],4) in [self.parameter.deltas[k],-self.parameter.deltas[k]]:
                        #These are multiplied by 2 just so that they come out to the right value (2 x number of colvars) in any situation
                        neighbor_count += 1
                        neighbor_value = (self.parameter.anchors[i][k] + self.parameter.anchors[j][k])/2
                        if self.parameter.max_grid_value:
                            if self.parameter.max_grid_value[k]:
                                if neighbor_value > self.parameter.max_grid_value[k]:
                                    neighbor_value = self.parameter.max_grid_value[k]
                        if self.parameter.min_grid_value:
                            if self.parameter.min_grid_value[k]:
                                if neighbor_value < self.parameter.min_grid_value[k]:
                                    neighbor_value = self.parameter.min_grid_value[k] 
                                        
                        if not self.parameter.milestone_delta:
                            if self.parameter.anchors[i][k] > self.parameter.anchors[j][k]: #less than
                                equation = 'step(' + str(neighbor_value) + '-' + str(self.parameter.names[k]) + ')*2'
                            else:
                                equation = 'step(' + str(self.parameter.names[k]) + '-' + str(neighbor_value) + ')*2'
                        else:
                            neighbor_array = [neighbor_value - self.parameter.milestone_delta, neighbor_value + self.parameter.milestone_delta]
                            equation = 'step(' + str(neighbor_array[1]) + '-' + str(self.parameter.names[k]) + ') + step(' + str(self.parameter.names[k]) + '-' + str(neighbor_array[0]) + ')'   
                        values.append(equation)
                        #values.append((self.parameter.anchors[i][k] + self.parameter.anchors[j][k])/2)
                    elif self.parameter.grid_pbc != False and abs(self.parameter.anchors[i][k]) + abs(self.parameter.anchors[j][k]) + self.parameter.deltas[k] == 360:
                        neighbor_count += 1
                    #elif self.parameter.pbc and self.parameter.anchors[i][k] + self.parameter.anchors[j][k] == self.parameter.deltas[k]:
                    #    neighbor_count += 1                    
                if neighbor_count == 1 and length_count == len(self.parameter.deltas) - 1:
                    if ms_list_form:
                        milestone = sorted([i+1, j+1])
                        milestone = 'MS' + str(i+1) + '_' + str(j+1)
                        if milestone in self.parameter.grid_ignore:
                            pass
                        else:
                            MS_list.add(milestone)
                    else:
                        MS_list[str(i+1) + '_' + str(j+1)] = (neighbor_value, values, 'g')
                values = []
        if not ms_list_form and self.parameter.grid_caps == True:
            MS_list = self.grid_caps(MS_list)
        if self.parameter.corners and not ms_list_form:
                MS_list = self.find_corners(MS_list)
        return MS_list
    
    
    def grid_caps(self, MS_list):
        '''
        Grid caps are what I am calling a situation like this: 
        l-----l-----l So we have three anchors, but there is sort of a gap at the edges. Since there is not a fourth anchor, there is not a set of 
        l  o  l  o  l step equations in our dictionary to terminate these. These "caps" (which are illustrated with the x's), just kind of make it so 
        l-----xxxxxxx trajectories don't escape and run wild over to the right. When these are not hit, the trajectories are terminated, but since there 
        l  o  x       is not a "real" anchor it is not taken into consideration for the results. It is just ignored in our computations.
        l-----x
        '''

        for i in range(len(self.parameter.anchors)):
            for j in range(len(self.parameter.anchors[0])):
                anchor_values_high = self.parameter.anchors[i].copy()
                anchor_values_low = self.parameter.anchors[i].copy()
                anchor_values_low[j] = self.parameter.anchors[i][j] - self.parameter.deltas[j]
                anchor_values_high[j] = self.parameter.anchors[i][j] + self.parameter.deltas[j]
                low_neighbor = round(anchor_values_low[j] + self.parameter.deltas[j]/2,4)
                high_neighbor = round(anchor_values_high[j] - self.parameter.deltas[j]/2,4)
                if self.parameter.max_grid_value and not self.parameter.pbc_names:
                    if self.parameter.max_grid_value[j]:
                        if high_neighbor > self.parameter.max_grid_value[j]:
                            high_neighbor = self.parameter.max_grid_value[j]
                        if anchor_values_high[j] > self.parameter.max_grid_value[j]:
                            anchor_values_high[j] = self.parameter.max_grid_value[j]
                if self.parameter.min_grid_value and not self.parameter.pbc_names:
                    if self.parameter.min_grid_value[j]:
                        if low_neighbor < self.parameter.min_grid_value[j]:
                            low_neighbor = self.parameter.min_grid_value[j] 
                        if anchor_values_low[j] > self.parameter.min_grid_value[j]:
                            anchor_values_low[j] = self.parameter.min_grid_value[j]
                neighbor_name = self.parameter.names[j]
                keep_high = True
                keep_low = True
                for k in range(len(self.parameter.anchors)):
                    low_count = 0
                    high_count = 0
                    count = 0
                    values = []
                    for l in range(len(anchor_values_high)):
                        count += 1
                        if round(anchor_values_high[l],3) == round(self.parameter.anchors[k][l],3):
                            high_count += 1
                        if round(anchor_values_low[l],3) == round(self.parameter.anchors[k][l],3):
                            low_count += 1                                              
                        if count == len(anchor_values_high):
                            if high_count == count:
                                keep_high = False
                            if low_count == count:
                                keep_low = False
                length = []
                for m in range(len(self.parameter.names)):
                    if self.parameter.names[m] == neighbor_name:
                        continue
                    length.append('step(' + self.parameter.names[m] + '-' + str(self.parameter.anchors[i][m] - self.parameter.deltas[m]/2) + ')' + ' + step(' +  str(self.parameter.anchors[i][m] + self.parameter.deltas[m]/2) + '-' + self.parameter.names[m] + ')')
                if len(length) == 1:
                    length = length[0]
                else:
                    '+'.join(length)
                if keep_high == True:
                    if not self.parameter.milestone_delta:
                        equation = 'step(' + str(self.parameter.names[j]) + '-' + str(high_neighbor) + ')*2'
                    else:
                        neighbor_array = [high_neighbor - self.parameter.milestone_delta,high_neighbor + self.parameter.milestone_delta]
                        equation = 'step(' + str(self.parameter.names[j]) + '-' + str(neighbor_array[0]) +') + step(' + str(neighbor_array[1]) + '-' + str(self.parameter.names[j]) + ')'
                    MS_list[str(i+1) + '_' + self.parameter.names[j] + '-high'] = [high_neighbor,[equation, length],'g']
                if keep_low == True:
                    if not self.parameter.milestone_delta:
                        equation = 'step(' + str(low_neighbor) + '-' + str(self.parameter.names[j]) + ')*2'
                    else:
                        neighbor_array = sorted([low_neighbor - self.parameter.milestone_delta,low_neighbor + self.parameter.milestone_delta])
                        equation = 'step(' + str(self.parameter.names[j]) + '-' + str(neighbor_array[0]) +') + step(' + str(neighbor_array[1]) + '-' + str(self.parameter.names[j]) + ')'
                    MS_list[str(i+1) + '_' + self.parameter.names[j] + '-low'] = [low_neighbor, [equation, length],'g']
        return MS_list

        
    def get_rmsd(self, anchors):
        equation = []
        for i in range(len(anchors)):
            equation.append('(' + self.parameter.names[i] + '-' + str(anchors[i]) + ')^2')
        equation = " + ".join(equation)
        equation = 'sqrt(' + equation + ')'
        return equation
    
        
    def read_milestone_folder(self):
        '''
        This function reads in the names of the folders in crd
        Arguments: None
        Returns a set() with every folder name in crd  (with MS in front so it is the same
                                                        format as MS_list)
        '''
        MS_list = set()
        for i in range(1, self.parameter.AnchorNum):
            for j in range(i, self.parameter.AnchorNum + 1):
                name = str(i) + '_' + str(j)
                if os.path.exists(self.parameter.crdPath + '/' + name):
                    MS_list.add('MS' + name)
        return MS_list
       
    def find_corners(self, MS_list):
        combinations = list(itertools.product([0,1], repeat=len(self.parameter.anchors[0])))
        corners = dict()
        anchors_list = []
        for i in self.parameter.anchors:
            anchors_list.append(list(i))
        for anchor in range(len(self.parameter.anchors)):
            values = list(self.parameter.anchors[anchor])
            total = dict()
            for i in range(len(combinations)):
                values = list(self.parameter.anchors[anchor])
                for j in range(len(combinations[0])):
                    if combinations[i][j] == 0:
                        continue
                    else:
                        values[j] += self.parameter.deltas[j]
                for k in range(len(anchors_list)):
                    if values == anchors_list[k]:
                        total[k+1] = values
            if len(total) != len(combinations):
                continue
            else:
                delta_sum = 0
                for i in range(len(self.parameter.deltas)):
                    delta_sum += self.parameter.deltas[i]
                delta_sum = delta_sum
                for i in total.keys():
                    for j in total.keys():
                        anchor_difference = 0
                        if i == j:
                            continue
                        for anchor_values in range(len(total[i])):
                            anchor_difference += abs(total[i][anchor_values] - total[j][anchor_values])
                        if anchor_difference != delta_sum:
                            continue #not diagonal
                        milestone = str(i) + '_' + str(j)
                        #going from i to j
                        step_functions = []
                        center = np.zeros(self.parameter.colvarsNum)
                        for anchor_values in range(len(total[i])):
                            center[anchor_values] = (total[i][anchor_values] + total[j][anchor_values])/2
                            if total[i][anchor_values] > center[anchor_values]:
                                step_functions.append('2*step(' + str(center[anchor_values]) + '-' + self.parameter.names[anchor_values] + ')')
                            else:
                                step_functions.append('2*step(' + self.parameter.names[anchor_values] + '-' + str(center[anchor_values]) + ')')
                        #step_functions.append('-1')
                        MS_list[milestone] = [list(center), step_functions,'c']
        return MS_list
        '''
        for milestone in MS_list.keys():
            anchors = get_anchors(milestone)
            values = []
            for i in anchors:
                values.append(list(self.parameter.anchors[i-1]))
            center = values[0]
            for i in range(self.parameter.colvarsNum):
                if values[0][i] != values[1][i]:
                    neighbor = i   
                    center[i] = (values[0][i] + values[1][i])/2
            for i in range(self.parameter.colvarsNum):
                for j in range(0,2):
                    if i == neighbor:
                        continue
                    tmp = center.copy()
                    if j == 0:
                        tmp[i] = values[0][i] + self.parameter.deltas[i]/2
                    else:
                        tmp[i] = values[0][i] - self.parameter.deltas[i]/2
                    tmp = list(tmp)
                    print(combinations)
                    name = [0,0]
                    for l in range(len(combinations)):
                        naming_corner = tmp.copy()
                        for k in range(len(combinations[0])):
                            if combinations[l][k] == 0:
                                naming_corner[k] -= self.parameter.deltas[k]/2
                            else:
                                naming_corner[k] += self.parameter.deltas[k]/2
                        total.append(naming_corner)
                    options = list(itertools)
                    for k in total:
                        if k ==
                        for k in range(len(name)):
                            for anchor_values in range(len(self.parameter.anchors)):
                                if naming_corner == list(self.parameter.anchors[anchor_values]):
                                    name.append(anchor_values + 1)
                    if 0 in name:
                        continue
                    step_function = []
                    for step in range(self.parameter.colvarsNum):
                        step_function.append('step(' + self.parameter.names[step] + '-' + str(tmp[step]) + ')')
                    corner_name = ('_').join([str(m) for m in name])
                    step_function = '+'.join(step_function)
                    corners[corner_name] = (tmp,step_function,'c')
                    name = name.reverse()
                    
        MS_list.update(corners)
        print(corners.keys())
        '''
        return MS_list
        
            
    def get_rmsd(self, anchors):
        equation = []
        for i in range(len(anchors)):
            equation.append('(' + self.parameter.names[i] + '-' + str(anchors[i]) + ')^2')
        equation = " + ".join(equation)
        equation = 'sqrt(' + equation + ')'
        return equation
               
if __name__ == '__main__':
    from parameters import *
    #from run import *
    from colvar import *
    
    new = parameters()
    new.initialize()
    #jobs = run(new)
    new.iteration = 1
    new.pbc = True
    new.traj_per_script = [2,2,2]
    new.deltas = [30,30]
    new.initial_traj = 6
    new.milestone_search = 2
    new.max_grid_value = [None,180]
    new.min_grid_value = [None,-180]
    new.pbc_names = ['psi']
    #new.corners = True
    new.MS_list = milestones(new).initialize(status=0)
    print(new.all_grid['1_2'])
    #milestones(new).find_corners()
    #milestones(new).get_final_ms(new.crdPath + '15_16/1/14')
    #milestones(new).get_final_ms_grid(new.ScMilesPath)
