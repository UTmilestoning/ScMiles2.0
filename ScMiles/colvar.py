# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 09:58:01 2020

@author: allis
"""

# -*- coding: utf-8 -*-
'''
This code generates the colvar configuration file that required by NAMD.
Two constraints will be considered:
    1. RMSD(x, anchor_a) = RMSD(x, anchor_b).
    2. RMSD(x, any_anchors_besides_a_or_b) > RMSD(x, anchor_a) &&
       RMSD(x, any_anchors_besides_a_or_b) > RMSD(x, anchor_b).
       
Note:        
    RMSD(x, anchor_a): the root mean square displacement from anchor_a to x
'''

#from log import log
import os
import re
import pandas as pd
import itertools
from fileinput import FileInput
import numpy as np
from additional_functions import *


class colvar:
    def __init__(self, parameter, anchor1=None, anchor2=None, 
                 free=None, initial=None, step=None, 
                 config_path=None, colvar_text=None, pdb_sampling=None):
        self.parameter = parameter
        self.anchor1 = anchor1
        self.anchor2 = anchor2
        self.colvars_number = len(self.parameter.anchors[0])
        self.step = step
        self.anchors_list = []
        self.sub_colvar = None
        self.config_path = self.parameter.ScMilesPath + "/colvar_free.conf" if self.step != 'sample' else self.parameter.ScMilesPath + "/colvar.conf"
        if self.parameter.pdb_sampling == True:
            self.pdb_sampling = True
        else:
            self.pdb_sampling = False        

    def __exit__(self, exc_type, exc_value, traceback):
        return 

    def __repr__(self) -> str:
        return ('Colvar generator')             
                     
    def __collective_vari(self, name=None, coeff=None, space=0, anchor1=None, anchor2=None):
        '''Saves all text from colvar.txt for each name (so not rmsd)'''
        tmp = []
        count = 0
        section = 1
        section = []
        with open(file=self.parameter.inputPath+'/colvar.txt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                if line.split() == []:
                    continue
                if 'colvar' in line.split() or 'rmsd' in line.split() or 'dihedral' in line.split(): #may need to add other things here
                    tmp.append(section)
                    section = []
                section.append(line)
        tmp.append(section)
        repeat = []
        fconf = open(self.config_path, 'a')
        for i in range(len(tmp)):
            repeat_colvar = False
            for line in tmp[i]:
                if 'anchornumber' in line:
                    repeat_colvar = True
                    line = line.replace('anchornumber',str(anchor1))
                print("  " * space + "  " + line, file=fconf)
            if repeat_colvar == True:
                repeat.append(tmp[i])
        for i in range(len(repeat)):
            for line in repeat[i]:
                if 'anchornumber' in line:
                    line = line.replace('anchornumber',str(anchor2))
                print("  " * space + "  " + line, file=fconf)
        fconf.close()

        '''
        for i in range(1,self.colvars_number + 1):
            if self.variables[i-1] == '':
                log("Colvar Error. Please name your colvars")
        '''

    def __substitution(self):
        sub = None
        name = None
        count = 2
        with open(self.parameter.inputPath + '/custom.colvar') as f:
            for line in f:
                line = line.split()
                if 'name' in line:
                    count += 1
                    name = line[1]
                if 'customFunction' in line:
                    sub = line[1]
                    break
        if name and sub:
            with FileInput(files=self.config_path, inplace=True) as f:
                for line in f:
                    line = line.strip()
                    if name in line and 'customFunction' in line:
                        line = line.replace(name, sub)
                    print(line)
        if self.parameter.milestone_search in (2,3) and self.step == 'sample':
            with FileInput(files=self.config_path, inplace=True) as f:
                for line in f:
                    line = line.strip()
                    if '#' in line:
                        line = 'customFunction ' + sub
                    print(line)

    '''    
    def __rmsd_to_anchor(self, anchor, coeff=None, space=0):
        #Used in "free" case, replaces "anchor" with the corresponding number in anchors.txt
        # scriptPath = os.path.dirname(os.path.abspath(__file__))
        # inputdir = os.path.abspath(os.path.join(scriptPath, os.pardir)) + '/my_project_input'
        tmp = []
        count = 0
        first = True
        section = 1
        name_get = False
        with open(file=self.parameter.inputPath + '/colvar.txt') as f:
            for line in f:
                if '{' in line:
                    first = False
                    count += 1
                if '}' in line:
                    count -= 1
                    if count == 0 and first == False:
                        if section == self.colvars_number + 1: 
                            tmp.append(line + '\n')
                        section += 1
                        continue
                if section == self.colvars_number + 1: 
                    if 'name' in line and name_get == False:
                        line = "  name rmsd" + str(anchor)
                        name_get = True
                    if 'anchor' in line:
                        if self.parameter.pbc_names:
                            line = self.__pbc_distance(anchor)
                        else:
                            for i in range(0, self.colvars_number): 
                                line = line.replace("anchor", '('+str(self.parameter.anchors[anchor-1][i])+')', 1)
                    tmp.append(line + '\n')
                
        fconf = open(self.config_path, 'a')
        for line in tmp:
            print("  " * space + "  " + line, file=fconf)
        fconf.close()
    '''
    def __rmsd_to_anchor(self, anchor, coeff=None, space=0):
        '''Used in "free" case, replaces "anchor" with the corresponding number in anchors.txt'''
        tmp = []
        keep = False
        rmsd = []
        tmp.append('colvar {')
        tmp.append('name rmsd' + str(anchor))
        function = self.__pbc_distance(anchor)
        tmp.append('customFunction ' + str(function))
        
        with open(file=self.parameter.inputPath + '/colvar.txt') as f:
            for line in f:
                line = line.strip()
                tmp.append(line)
        tmp.append('}')    
        fconf = open(self.config_path, 'a')
        for line in tmp:
            if 'anchornumber' in line:
                line = line.replace('anchornumber',str(anchor))
            print("  " * space + "  " + line, file=fconf)
        fconf.close()
     
    def __pbc_distance(self,anchor):
        line = ''
        line += 'sqrt('
        for colvar_index in range(self.colvars_number):
            name = self.parameter.names[colvar_index]
            anchor_value = self.parameter.anchors[anchor-1][colvar_index]
            absolute_value = '(' + name + '-' + str(anchor_value) + ')'
            if self.parameter.scale_rmsd:
                scale = str(self.parameter.scale_rmsd[colvar_index]) + '*'
            else:
                scale = ''
            if name in self.parameter.pbc_names:
                line += '('
                name_index = self.parameter.pbc_names.index(name)
                l_value = self.parameter.l[name_index]  
                l_half = l_value/2
                absolute_value = 'abs' + absolute_value
                line += scale + '(' + str(l_value) + '-' + absolute_value + ')*step(' + absolute_value + '-' + str(l_half) + ') + ' + scale + absolute_value + '*step(' + str(l_half) + '-' + absolute_value + '))^2'

                #e = 'e^'
                #reg_lambda = '(-' + str(self.parameter.distance_lambda) + '*' + absolute_value + ')'
                #pbc_lambda = '(-' + str(self.parameter.distance_lambda) + '*(360 -' + absolute_value + '))'
                #z = 'e^' + reg_lambda + '+e^' + pbc_lambda + ''
                #line += '(' + absolute_value + '*' + e + reg_lambda + '/(' + z + ')+(360-' + absolute_value + ')*' + e + pbc_lambda + '/(' + z + '))^2'
            else:
                if self.parameter.scale_rmsd:
                    absolute_value = scale+ absolute_value
                absolute_value = absolute_value + '^2'
                line += absolute_value
            line += '+'
        line = line.replace('--','+')
        line = line[0:-1]
        line += ')'
        return line   
     
    def __min_pbc(self, anchor):
        ms_values = []
        rmsd = []
        #Finding the names of our colvars, so phi and psi for Alanine Dipeptide in a vacuum
        colvar, names = self.get_colvars()
        #Looping through the number of colvars and appending each midpoint.
        #So like if we are doing milestone 1_2, we use the midpoint between these two as our reference point
        for i in range(len(self.parameter.anchors[0])):
            ms_values.append((self.parameter.anchors[self.anchor1-1][i] + self.parameter.anchors[self.anchor2-1][i])/2)
        #Now we check for l values (entered by the user as a list)
        #So if our first variable is periodic and our second is not, it would look something like this:
        #l_values = [360, 0] Zero just means it is not periodic
        for i in range(len(self.parameter.l_values)):
            #if our l value is zero, our variable is not periodic so we just do the regular rmsd
            if self.parameter.l_values[i] == 0:
                rmsd.append('(' + str(names[i]) + ' - ' + str(self.parameter.anchors[anchor-1][i]) + ')^2')
            #If we do have an l value at the index, we calculate:
            #regular: which is the straight distance between the two
            #pbc: which is our L value - regular, so our value if we were to go around the boundary
            #Then we see which one is smaller. If regular is smaller, we just do the regular rmsd
            #If pbc is smaller, we use the version with the L value
            else:
                regular = abs(ms_values[i] - self.parameter.anchors[anchor-1][i])
                pbc = (self.parameter.l_values[i]) - regular
                print(regular,pbc)
                if regular <= pbc:
                    rmsd.append('(' + str(names[i]) + ' - ' + str(self.parameter.anchors[anchor-1][i]) + ')^2')
                else:
                    rmsd.append('(' + str(self.parameter.l_values[i]) + ' - abs('+ str(names[i]) + ' - ' + str(self.parameter.anchors[anchor-1][i]) + '))^2')
        
        #Then this just puts our pieces of rmsd together (one term for each colvar)
        line = "  customFunction sqrt("
        for i in rmsd:
            i = i.replace('- -', '+ ')
            i = i.replace('--', '+ ')
            line += i + ' + '
        line = line[:-3]
        line += ')'
        return line
        
                  
    def __find_rmsd(self, milestone):
        #This only works for Voronoi cell for now.
        from math import sqrt
        anchors = []
        total = 0
        dist = 0
        for i in milestone:
            anchors.append(self.parameter.anchors[i-1])
        total = 0
        for i in range(len(anchors[0])): #this is for each colvar
            total += (anchors[0][i] - anchors[1][i])**2
        dist = sqrt(total)
        if dist > self.parameter.dist_cut:
            return False
        else:
            return True

    def generate(self):   
        '''This is the main function that generates colvars '''
        if self.step == 'seek': 
            outputFrequency = 1
        else:
            outputFrequency = self.parameter.colvarsTrajFrequency
            
        fconf = open(self.config_path, 'w+')
        print("colvarsTrajFrequency      {}".format(outputFrequency), file=fconf)
        print("colvarsRestartFrequency	 {}".format(self.parameter.colvarsRestartFrequency), file=fconf)
        if self.step == 'free' or self.step == 'seek':
            print("scriptedColvarForces on", file=fconf)
        if self.parameter.customColvars == True:
            print("", file=fconf)
            with open(file=self.parameter.inputPath + '/custom.colvar') as f_custom:
                for line in     f_custom:
                    print(line, file=fconf)
        fconf.close()
        
        if self.step == 'sample':
            crd_path = self.parameter.crdPath + '/' + str(self.anchor1) + '_' + str(self.anchor2)
            if not os.path.isfile(crd_path + '/seek.ms.pdb') and not os.path.isfile(crd_path + '/' + str(self.parameter.outputname) + '.coor'):
                self.pdb_sampling = True

        if self.step != 'sample':
            if self.parameter.milestone_search == 2 or self.parameter.milestone_search == 3:
                self.__grid_colvars()
            else:
                for i in range(self.parameter.AnchorNum):
                    self.__rmsd_to_anchor(i+1)
        else:
            if self.colvars_number == 1:
                self.__constraint1D1()
                self.__harmonic1D()
            elif self.parameter.milestone_search == 2 or self.parameter.milestone_search == 3:
                self.__grid_sampling()
            else:
                self.__constraint2D1()
                colvarList, centers = self.__constraint2D2()
                self.__harmonic2D()
                self.__harmonicWalls(colvarList, centers)
        if self.parameter.substitution == True:
            self.__substitution()

    def __constraint1D1(self):
        fconf = open(self.config_path, 'a')
        print("\ncolvar {", file=fconf)
        print("  name colv", file=fconf)
        fconf.close()
        self.__collective_vari(anchor1=self.anchor1, anchor2=self.anchor2)
        fconf = open(self.config_path, 'a')
        print("}\n\n", file=fconf)
        fconf.close()

    def __harmonic1D(self):
        fconf = open(self.config_path, 'a')
        print("\nharmonic {", file=fconf)
        print("  colvars colv", file=fconf)
        center = (self.parameter.anchors[self.anchor1-1][0] + self.parameter.anchors[self.anchor2-1][0]) / 2
        if self.parameter.pbc != [] and abs(self.anchor1 - self.anchor2) > 1:
            center = 180
        print("  centers {}".format(center), file=fconf)
        print("  forceConstant {}".format(self.parameter.forceConst[0]), file=fconf)
        print("}", file=fconf)
        fconf.close()

    def __constraint2D1(self):
        fconf = open(self.config_path, 'a')
        print("\ncolvar {", file=fconf)
        print("  name neighbor", file=fconf)
        customFunc = self.__custom_function(self.anchor1-1, self.anchor2-1)
        print(customFunc, file=fconf)
        fconf.close()   
        self.__collective_vari(space=1, anchor1=self.anchor1, anchor2=self.anchor2)    
        fconf = open(self.config_path, 'a')
        print("}\n\n", file=fconf)
        fconf.close()
                
    def __custom_function(self, anchor1, anchor2):
        '''Creates the customFunction for cases with more than one colvar'''
        customFunc = "  customFunction "
        for section in (1,2):
            if section == 1:
                anchor = anchor1
            else:
                anchor = anchor2
            if self.parameter.pbc_names:
                function = self.__pbc_distance(anchor + 1)
                customFunc = customFunc + function
                if section == 1:
                    customFunc += '-'
            else:
                customFunc = customFunc + 'sqrt('
                for i in range(1, self.colvars_number + 1):
                    customFunc = customFunc + '(' + self.parameter.names[i-1] + '-(' + str(self.parameter.anchors[anchor][i-1]) + '))^2'
                    if i != self.colvars_number:
                        customFunc = customFunc + ' + '
                if section == 1:
                    customFunc = customFunc + ') - '
                else:
                    customFunc = customFunc + ')'
            if 'anchornumber' in customFunc:
                customFunc = customFunc.replace('anchornumber',str(anchor+1))
        return customFunc          

    def __constraint2D2(self):
        colvarList = ""
        centers = ""
        for i in range(self.parameter.AnchorNum):
            if i + 1 != self.anchor1 and i + 1 != self.anchor2:
                fconf = open(self.config_path, 'a')
                print("colvar {", file=fconf)
                print("  name {}_{}".format(i + 1, self.anchor1), file=fconf)
                customFunc = self.__custom_function(i, self.anchor1-1)
                print(customFunc, file=fconf)
                colvarList += str(i + 1) + "_" + str(self.anchor1) + " "
                centers += "0 "
                fconf.close()
                self.__collective_vari(space=2, anchor1=i+1, anchor2 =self.anchor1)

                fconf = open(self.config_path, 'a')
                print("}\n", file=fconf)       
    
                print("colvar {", file=fconf)
                print("  name {}_{}".format(i + 1, self.anchor2), file=fconf)
                customFunc = self.__custom_function(i,self.anchor2-1)
                print(customFunc, file=fconf)
    
                colvarList += str(i + 1) + "_" + str(self.anchor2) + " "
                centers += "0 "
                fconf.close()
                self.__collective_vari(space=2, anchor1=i+1, anchor2=self.anchor2)

                fconf = open(self.config_path, 'a')
                print("}\n", file=fconf)
                fconf.close()
        return colvarList, centers
    
    def __harmonic2D(self, center=None, colvar=0):
        fconf = open(self.config_path, 'a')
        print("harmonic {", file=fconf)
        print("  colvars neighbor", file=fconf)
        if not center:
            center = 0
        print("  centers {}".format(str(center)), file=fconf)
        if self.pdb_sampling:
            print("  forceConstant {}".format((self.parameter.forceConst[colvar])/10), file=fconf)
            print("  targetForceConstant {}".format(self.parameter.forceConst[colvar]), file=fconf)
            print("  targetNumSteps {}".format(self.parameter.targetNumSteps), file=fconf)   
        else:
            print("  forceConstant {}".format(self.parameter.forceConst[colvar]), file=fconf)
        print("}", file=fconf)
        fconf.close()

    def __harmonicWalls(self, colvarList, centers, colvar=0):
        fconf = open(self.config_path, 'a')
        print("\n", file=fconf)
        print("harmonicWalls {", file=fconf)
        print("  colvars {}".format(colvarList), file=fconf)
        print("  lowerWalls {}".format(centers), file=fconf)
        if self.pdb_sampling:
            print("  lowerWallConstant {}".format((self.parameter.forceConst[colvar])/10), file=fconf)
            print("  targetForceConstant {}".format(self.parameter.forceConst[colvar]), file=fconf)
            print("  targetNumSteps {}".format(self.parameter.targetNumSteps), file=fconf)   
        else:
            print("  lowerWallConstant {}".format(self.parameter.forceConst[colvar]), file=fconf)
        print("}", file=fconf)
        fconf.close()

    def __grid_colvars(self):
        colvar_list, names = self.get_colvars()
        equations,grid = self.get_equations()
        #if self.parameter.corners:
            #equations = self.__grid_corners(equations)
        fconf = open(self.config_path, 'a')
        for eq in equations:
            print('colvar {', file = fconf)
            print('  name ' + eq[1], file = fconf)
            print('  customFunction ' + eq[0], file=fconf)
            for colvar in colvar_list:
                print("    " + colvar, file=fconf)
            print("}", file=fconf)
        fconf.close() 
        
        
    def get_equations(self):
        done = []
        ms_network = set()
        if self.step == 'free':
            ms_dne = [self.anchor1, self.anchor2]
        else:
            ms_dne = [self.anchor1]
        while ms_dne:
            current_anchor = ms_dne.pop(0)
            done.append(current_anchor)
            for i in self.parameter.all_grid.keys():
                if self.parameter.all_grid[i][2] == 'c':
                    continue
                i_list = i.split('_')
                for item in range(len(i_list)):
                    try:
                        i_list[item] = int(i_list[item])
                    except:
                        continue
                if current_anchor == i_list[0]:
                    i_reverse = str(i_list[1]) + '_' + str(i_list[0])
                    if i in ms_network or i_reverse in ms_network:
                        continue
                    elif i in self.parameter.grid_ignore or i_reverse in self.parameter.grid_ignore:
                        for j in i_list:
                            if j == current_anchor:
                                continue
                            elif j in done:
                                continue
                            else:
                                ms_dne.append(j)
                    else:
                        for item in range(len(i_list)):
                            i_list[item] = str(i_list[item])
                        if i_list == [str(self.anchor1), str(self.anchor2)] or i_list == [str(self.anchor2), str(self.anchor1)]:
                            continue
                          
                        ms_network.add('_'.join(i_list))
                        if self.parameter.corners:
                            for corner in self.parameter.all_grid.keys():
                                if self.parameter.all_grid[corner][2] == 'c':
                                    inside_anchor = corner.split('_')[0]
                                    if int(inside_anchor) == int(i_list[0]):
                                        ms_network.add(corner)
                                    
        equations = []
        grid = []
        count = 0
        for i in ms_network:
            equations.append([' + '.join(self.parameter.all_grid[i][1]),i])
            grid.append([count, i])
            count += 1
        return equations, grid
        
    def __grid_sampling(self):

        anchor1_coor = self.parameter.anchors[self.anchor1-1]
        anchor2_coor = self.parameter.anchors[self.anchor2-1]
        
        colvars, names = self.get_colvars()
        for i in range(len(colvars)):
            name_index = []
            for j in range(len(colvars[i])):
                if 'name' in colvars[i][j]:
                    name_index.append(j)
            if len(name_index) > 1:
                continue
            else:
                colvars[i].pop(name_index[0])
                
        names = self.parameter.names
        for i in range(len(anchor1_coor)):
            if anchor1_coor[i] != anchor2_coor[i]:
                neighbor = i
        delta_value = self.parameter.deltas.copy() #change this later, or have user input value
        center = (anchor1_coor[neighbor] + anchor2_coor[neighbor])/2
        if self.parameter.grid_pbc != False:
            if names[neighbor] == self.parameter.grid_pbc:
                if abs(anchor1_coor[neighbor]) + delta_value[neighbor]/2 == 180 and abs(anchor1_coor[neighbor]) == abs(anchor2_coor[neighbor]):
                    center = 180
        fconf = open(self.config_path, 'a')
        print("\ncolvar {", file=fconf)
        print("  name neighbor", file=fconf)
        for item in colvars[neighbor]:
            print("   {}".format(item), file=fconf)
        print("  }", file=fconf)
        count = 0
        for i in range(len(names)):
            if i == neighbor:
                continue
            count += 1
            print("\ncolvar {", file=fconf)
            print("  name length{}".format(count), file=fconf)
            for item in colvars[i]:
                print("   {}".format(item), file=fconf)
            print("  }", file=fconf)
        fconf.close()
        count = 0
        print(neighbor)
        self.__harmonic2D(center, neighbor)
        fconf = open(self.config_path, 'a')    
        for i in range(len(names)):
            if i == neighbor:
                continue
            count += 1
            walls = [anchor1_coor[i] + delta_value[i]/2, anchor1_coor[i] - delta_value[i]/2]
            walls.sort()
            lower_wall = walls[0]
            upper_wall = walls[1]
            walls_epsilon = abs(upper_wall-lower_wall)*0.1
            for j in range(len(walls)):
                if self.parameter.min_grid_value:
                    if walls[j] <= self.parameter.min_grid_value[i]:
                        walls[j] = int(walls[j] + walls_epsilon)
                if self.parameter.max_grid_value:
                    if walls[j] >= self.parameter.max_grid_value[i]:
                        walls[j] = int(walls[j] - walls_epsilon)
            print("\n", file=fconf)
            print("harmonicWalls {", file=fconf)
            print("  colvars length{}".format(count), file=fconf)
            print("  lowerWalls {}".format(walls[0]), file=fconf)
            print("  upperWalls {}".format(walls[1]), file=fconf)        
            print("  lowerWallConstant {}".format(self.parameter.forceConst[i]), file=fconf)
            print("  UpperWallConstant {}".format(self.parameter.forceConst[i]), file=fconf)
            print("} \n", file=fconf)           
        fconf.close()  
        
    def get_colvars(self):
        count = 0
        tmp = []
        colvar_list = []
        names = []
        keep = False
        with open(file = self.parameter.inputPath + '/colvar.txt') as f:
            for line in f:
                if line == '' or line == '\n':
                    continue
                if self.step != 'sample' and line.startswith('#'):
                    continue
                if self.step == 'sample':
                    if '#' in line and 'customFunction' not in line:                            
                        colvar_list.append(tmp)
                        tmp = []
                        continue
                tmp.append(line)
        print(colvar_list)
        if self.step == 'sample':
            return colvar_list, names
        return tmp, names
    
    
    def __grid_corners(self, equations):
        #find corner values
        a1_values = set()
        a2_values = set()
        for i in range(self.colvars_number):
            values = []
            for j in (self.anchor1,self.anchor2):
                a1_values.add(self.parameter.anchors[self.anchor1-1][i] + self.parameter.deltas[i]/2)
                a1_values.add(self.parameter.anchors[self.anchor1-1][i] - self.parameter.deltas[i]/2)
                a2_values.add(self.parameter.anchors[self.anchor2-1][i] + self.parameter.deltas[i]/2)
                a2_values.add(self.parameter.anchors[self.anchor2-1][i] - self.parameter.deltas[i]/2)
        all_corners = list(itertools.product(a1_values,a2_values))
        for i in range(len(all_corners)):
            all_corners[i] = list(all_corners[i])
            for j in self.parameter.all_grid.keys():
                if j[0] == all_corners[i]:
                    print(j)
        print(self.parameter.all_grid)
        print(all_corners)
        return equations
                
            
        
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
        '''
        if self.parameter.software == 'namd':
            state = path + "/stop.colvars.state"
        else:
            state = path + '/' + self.parameter.outputname + '.colvars.traj'
        if not os.path.isfile(state):
            return -1, [0, 0]  
        final_ms = [0, 0]
        equation = 0
        # read state file generated at termination point 
        # smallest rmsd indicates new cell #
        if self.parameter.software == 'namd':
            RMSDs, lifetime = self.read_state(state)
        else:
            RMSDs, lifetime, colvar_values = self.read_traj(state)
            if RMSDs == None:
                return -1, [0,0]
        final_ms[0] = RMSDs.index(sorted(RMSDs)[0]) + 1
        firstRMSD = self.parameter.colvarsNum + 1
        colvar_number = len(self.parameter.anchors[0])
        if self.parameter.milestone_search == 2:
            if not self.parameter.milestone_delta:
                value = 2*(colvar_number)
            else:
                value = 2*colvar_number
            for i in range(len(RMSDs)):
                if int(RMSDs[i]) == value:
                    equation = i
            count = 0
            if self.parameter.software == 'namd':
                equation += colvar_number
                if RMSDs.count(value) == 0:
                    return -1,0
                if RMSDs.count(value) == 1:
                    with open(path + '/stop.colvars.state') as f:
                        for line in f:
                            if 'name' in line:
                                if count == equation:
                                    line = line.split()
                                    print(line)
                                    end = line[1]
                                    if 'low' in end or 'high' in end:
                                        return -1, [0,0]
                                    final_ms = get_anchors(end)
                                    break
                                count += 1
                else:
                    #print(path)
                    possible = []
                    final_step = []
                    prev_step = []
                    count = 0
                    next_line = False
                    with open(path + '/stop.colvars.state') as f:
                        for line in f:
                            if next_line == True:
                                line = line.split()
                                final_step.append(float(line[1]))
                                next_line = False
                            if 'name' in line:
                                if count < colvar_number:
                                    next_line = True                                
                                elif RMSDs[count - colvar_number] == value:
                                    line = line.split()
                                    possible.append(line[1])
                                count += 1
                            if count >= colvar_number + len(RMSDs):
                                break

                    if not os.path.isfile(path + '/' + self.parameter.outputname + '.colvars.traj') or os.stat(path + '/' + self.parameter.outputname + '.colvars.traj').st_size == 0:
                        return -1, [0,0]
                    with open(path + '/' + self.parameter.outputname + '.colvars.traj') as f:
                        for line in f:
                            last_line = line
                    last_line = last_line.split()
                    last_line.pop(0)
                    end = None
                    for i in range(colvar_number):
                        prev_step.append(float(last_line[i]))
                    for i in range(colvar_number):
                        for j in range(len(possible)):
                            range_values = sorted([final_step[i],prev_step[i]])
                            if not self.parameter.milestone_delta:
                                if range_values[0] < self.parameter.all_grid[possible[j]][0] < range_values[1]:
                                    end = possible[j]
                                    break
                            else:
                                if range_values[0] + self.parameter.milestone_delta < self.parameter.all_grid[possible[j]][0] < range_values[1] + self.parameter.milestone_delta:
                                    end = possible[j]
                                    break
                    if end == None:
                        return -1, [0,0]
                    if 'low' in end or 'high' in end:
                        return -1, [0,0]
                    
                    final_ms = get_anchors(end)                        

            else:
                equation += 1
                if RMSDs.count(value) == 1:
                    with open(path + '/plumed_free.dat') as f:
                        for line in f:
                            if 'LABEL=stepfunction' in line:
                                count += 1
                            if count == equation:
                                end = line.replace('LABEL=stepfunction','')
                                end = end.replace(' ', '')
                                #print(end)
                                if 'low' in end or 'high' in end:
                                    return -1, [0,0]
                                #print(end)
                                final_ms = get_anchors(end)
                                #print(end)
                                break
                                #print(end)
                else:
                    possible = []
                    final_step = colvar_values[0]
                    prev_step = colvar_values[1]
                    with open(path + '/plumed_free.dat') as f:
                        for line in f:
                            if 'LABEL=stepfunction' in line:
                                count += 1
                            if count == equation:
                                end = line.replace('LABEL=stepfunction','')
                                end = end.replace(' ', '') 
                                end = end.strip()
                                final_ms = get_anchors(end)
                                possible.append(end)
                                break
                    for i in range(colvar_number):
                        range_values = sorted([final_step[i],prev_step[i]])
                        #print(range_values)
                        for j in range(len(possible)):
                            if not self.parameter.milestone_delta:
                                if range_values[0] < self.parameter.all_grid[possible[j]][0] < range_values[1]:
                                    end = possible[j]
                                    break
                            else:
                                if range_values[0] + self.parameter.milestone_delta < self.parameter.all_grid[possible[j]][0] < range_values[1] + self.parameter.milestone_delta:
                                    end = possible[j]
                                    break
                    
                    if end == None:
                        return -1, [0,0]
                    if 'low' in end or 'high' in end:
                        return -1, [0,0]
                    final_ms = get_anchors(end)
        else:
            if self.parameter.pbc:
                if final_ms[0] == self.parameter.AnchorNum or final_ms[0] == 1:
                    # open traj file to read the very last output
                    # smallest rmsd indicates previous cell #
                    traj = path + "/" + self.parameter.outputname + ".colvars.traj"
                    #print(traj)
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
                    #print(path + '/start.txt')
                    start_ms = pd.read_csv(path + '/start.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                except:
                    return -1, [0, 0]
                if start_ms[0] not in final_ms and start_ms[1] not in final_ms:
                    return -1, [0, 0]
        final_ms.sort()
        final_info = path + '/end.txt'
        if str(final_ms[0]) == str(final_ms[1]):
            return -1, [0,0]
        if not os.path.isfile(final_info):
            with open(final_info, 'w+') as f1:
                print(final_ms[0], final_ms[1], file=f1)    
        time_info = path + '/lifetime.txt'
        lifetime *= self.parameter.timeFactor
        if not os.path.isfile(time_info):
            with open(time_info, 'w+') as f1:
                print(lifetime, file=f1)    
        return lifetime, final_ms
    
    def read_state(self, path):
        file = open(path, 'r').read()
        time = int(re.findall(r"[-+]?\d*\.\d+|\d+", file)[0])
        numbers = re.findall(r"[-+]?\ *[0-9]+\.?[0-9]*(?:[eE]\ *[-+]?\ *[0-9]+)", file)
        numbers.pop(0)
        numbers = numbers[self.parameter.colvarsNum:]
        values = [abs(float(x)) for x in numbers]
        return values, time

    def read_traj(self, path):
        last_line = None
        colvar_values = [[],[]]
        second_to_last = None
        third_to_last = None
        with open(path) as r:
            for line in r:
                 if line.startswith('#'):
                     continue
                 line = line.split()
                 if line != []:
                    third_to_last = second_to_last
                    second_to_last = last_line
                    last_line = line
        #print(path)
        try:
            time = last_line.pop(0)
        except:
            return None, None, None
        for i in range(len(self.parameter.anchors[0])):
            last_line.pop(-1)
        values = [abs(float(x)) for x in last_line]
#        print(second_to_last)
#        print(third_to_last)
        if not third_to_last:
            return None, None, None
        if self.parameter.milestone_search == 2:
            time = second_to_last.pop(0)
            for i in range(len(self.parameter.anchors[0])):
                colvar_values[1].append(float(second_to_last.pop(-1)))
                colvar_values[0].append(float(third_to_last.pop(-1)))
            values = [abs(float(x)) for x in second_to_last]
        colvar_values[1].reverse()
        colvar_values[0].reverse()
        time = float(time) * 1000
        return values, time, colvar_values

    
if __name__ == '__main__':
    from parameters import *
    from milestones import *
    new = parameters()
    new.restart = False
    new.initialize()
    new.AnchorNum = 42
    new.software = 'gromacs'
    new.deltas = [30,30]
    new.substitution = True
    #new.deltas = [30,30]
    #new.grid_ignore= ['13_20','20_21','20_27','19_20','18_19','17_18']
    new.milestone_search = 2
    #new.scale_rmsd = [0.4,1.0]
    #new.all_grid = milestones(new).grid_ms()
    #new.distance_lambda = 0.1
    #new.pbc_names = ['chi']
    #new.l = [360]
    #print(new.all_grid)
    #new.names = ['chi','covrmsdanchornumber']
    #new.MS_list = milestones(new).initialize()
    #print(new.forceConst)
    #new.deltas = [30,30,30]
    #colvar(new, anchor1=1, anchor2=2).generate()
    #colvar(new, step='sample', anchor1=19, anchor2=20).generate()
    #colvar(new, step='seek', anchor1=5, anchor2=6).generate()
    new.milestone_search = 2
    print(colvar(new, step='sample', anchor1=1,anchor2=2).generate())
