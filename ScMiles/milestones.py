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
                launch(self.parameter, 'seek').launch_trajectories(possible_ms, lastLog)
                self.parameter.restart = False
                lastLog=None
            # check if reactant and product are connected
            # read folders to get the milestones list 
            return MS_list
        else: #grid
            MS_list = self.grid_ms(ms_list_form = 'yes')
        print(MS_list)
        for ms in MS_list:
            [anchor1, anchor2] = get_anchors(ms)
            create_folder(self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2))            
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
                    elif self.parameter.grid_pbc != False and abs(self.parameter.anchors[i][k]) + abs(self.parameter.anchors[j][k]) + self.parameter.deltas[k] == 360:
                        neighbor_count += 1
                    #elif self.parameter.pbc and self.parameter.anchors[i][k] + self.parameter.anchors[j][k] == self.parameter.deltas[k]:
                    #    neighbor_count += 1
                if neighbor_count == 1 and length_count == len(self.parameter.deltas) - 1:
                    if ms_list_form:
                        MS_list.add('MS' + str(i+1) + '_' + str(j+1))
                    else:
                        MS_list.append([i+1,j+1])
        return MS_list
        
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
    new.initial_traj = 6
    new.milestone_search = 1
    new.AnchorNum = 12
    new.MS_list = milestones(new).initialize(status=0)
    #milestones(new).get_final_ms(new.crdPath + '15_16/1/14')
    #milestones(new).get_final_ms_grid(new.ScMilesPath)
