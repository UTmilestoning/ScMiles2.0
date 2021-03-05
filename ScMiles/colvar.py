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
import numpy as np

class colvar:
    def __init__(self, parameter, anchor1=None, anchor2=None, 
                 free=None, initial=None, step=None, 
                 config_path=None, variables=None):
        self.parameter = parameter
        self.anchor1 = anchor1
        self.anchor2 = anchor2
        self.variables = []
        self.free = free
        self.colvars_number = len(self.parameter.anchors[0])
        for i in range(1, self.colvars_number + 1):
            self.variables.append("")
        self.initial = initial
        self.step = step
        self.config_path = self.parameter.ScMilesPath + "/colvar_free.conf" if self.step != 'sample' else self.parameter.ScMilesPath + "/colvar.conf"
        

    def __exit__(self, exc_type, exc_value, traceback):
        return 

    def __repr__(self) -> str:
        return ('Colvar generator')             

#    def __collective_vari_1(self, name=None, coeff=None, space=0):
#        '''
#        Change this function for different 1D case.
#        Follow the format in Colvars to define a collective variable.
#        For the commands below, it generates the following.
#        
#          dihedral {
#            name psi
#            group1 atomNumbers 7
#            group2 atomNumbers 9
#            group3 atomNumbers 15
#            group4 atomNumbers 17
#          }
#          
#        '''
#        fconf = open(self.config_path, 'a')
#        print("  " * space + "  dihedral {", file=fconf)
#        if name:
#            print("  " * space + "    name {}".format(name), file=fconf) 
#        if coeff:
#            print("  " * space + "    componentCoeff {}".format(coeff), file=fconf) 
#        print("  " * space + "    group1 atomNumbers 7", file=fconf)
#        print("  " * space + "    group2 atomNumbers 9", file=fconf)
#        print("  " * space + "    group3 atomNumbers 15", file=fconf)
#        print("  " * space + "    group4 atomNumbers 17", file=fconf)
#        print("  " * space + "  }", file=fconf)
#        fconf.close()

    def __get_colvar_names(self):
        '''Stores colvar names in array "variables"'''
        count = 0
        section = 1
        with open(file=self.parameter.inputPath + '/colvar.txt') as f:
            for line in f:
                if '{' in line:
                    count += 1
                if '}' in line:
                    count -= 1
                    if count == 0:
                        section += 1
                if "name" in line:
                    info = line.split("#")[0].split()
                    if len(info) >= 2 and info[0] == "name":
                        self.variables[section-1] = str(info[1])
                        if self.colvars_number == section:
                            break
                     
    def __collective_vari(self, name=None, coeff=None, space=0):
        '''Saves all text from colvar.txt for each name (so not rmsd)'''
        tmp = []
        count = 0
        section = 1
        with open(file=self.parameter.inputPath+'/colvar.txt') as f:
            for line in f:
                if '{' in line:
                    count += 1
                if '}' in line:
                    count -= 1
                    if count == 0:
                        section += 1
                        if section > self.colvars_number:
                            tmp.append(line + '\n')
                            break
                tmp.append(line + '\n')
        fconf = open(self.config_path, 'a')
        for line in tmp:
            print("  " * space + "  " + line, file=fconf)
        fconf.close()
        for i in range(1,self.colvars_number + 1):
            if self.variables[i-1] == '':
                log("Colvar Error. Please name your colvars")

        
        
        
    def __rmsd_to_anchor(self, anchor, coeff=None, space=0):
        '''Used in "free" case, replaces "anchor" with the corresponding number in anchors.txt'''
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
                        for i in range(0, self.colvars_number):
                            line = line.replace("anchor", '('+str(self.parameter.anchors[anchor-1][i])+')', 1)
                    tmp.append(line + '\n')
                
        fconf = open(self.config_path, 'a')
        for line in tmp:
            print("  " * space + "  " + line, file=fconf)
        fconf.close()
        

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
                for line in f_custom:
                    print(line, file=fconf)
        fconf.close()
        
        if self.step != 'sample':
            if self.parameter.milestone_search == 2:
                if self.step == 'seek':
                    self.__grid_seek()
                else:
                    self.__grid_colvars()
            else:
                for i in range(self.parameter.AnchorNum):
                    self.__rmsd_to_anchor(i+1)       
        else:
            self.__get_colvar_names()
            if self.colvars_number == 1:
                self.__constraint1D1()
                self.__harmonic1D()
            elif self.parameter.milestone_search == 2:
                self.__grid_sampling()
            else:
                self.__constraint2D1()
                colvarList, centers = self.__constraint2D2()
                self.__harmonic2D()
                self.__harmonicWalls(colvarList, centers)

    def __grid_seek(self):
        equations = []
        values = []
        count = 0
        colvar_list, names = self.__get_colvars()
        anchor_coordinates = self.parameter.anchors[self.anchor1 - 1]
        for i in range(len(anchor_coordinates)):
            values.append([anchor_coordinates[i] - self.parameter.deltas[i]/2, anchor_coordinates[i] + self.parameter.deltas[i]/2])
        for i in values:
            i.sort()
        for i in range(len(values)):
            equations.append(self.__create_step(values[i][0], names[i]))
            for k in range(len(values)):
                if  k == i:
                    continue
                equations.append(self.__create_step(values[k][1], names[k], names[k], values[k][0]))
            equations.append(self.__create_step(names[i], values[i][1]))
            for k in range(len(values)):
                if  k == i:
                    continue
                equations.append(self.__create_step(values[k][1], names[k], names[k], values[k][0]))
                
        fconf = open(self.config_path, 'a')
        for eq in equations:
            count += 1
            eq = eq.replace('- -', '+ ')
            eq = eq.replace('--', '+')
            print('colvar {', file = fconf)
            print('  name ' + str(count), file = fconf)
            print('  customFunction ' + eq, file=fconf)
            for i in range(len(names)):
                if names[i] in eq:
                    for line in colvar_list[i]:
                        print("    " + line, file=fconf)
            print("}", file=fconf)   
        fconf.close()

    def __constraint1D1(self):
        fconf = open(self.config_path, 'a')
        print("\ncolvar {", file=fconf)
        print("  name colv", file=fconf)
        fconf.close()
        self.__collective_vari()
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
        self.__collective_vari(space=1)    
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
            customFunc = customFunc + 'sqrt('
            for i in range(1, self.colvars_number + 1):
                customFunc = customFunc + '(' + self.variables[i-1] + '-(' + \
                    str(self.parameter.anchors[anchor][i-1]) + '))^2'
                if i != self.colvars_number:
                    customFunc = customFunc + ' + '
            if section == 1:
                customFunc = customFunc + ') - '
            else:
                customFunc = customFunc + ')'
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
                self.__collective_vari(space=2)

                fconf = open(self.config_path, 'a')
                print("}\n", file=fconf)       
    
                print("colvar {", file=fconf)
                print("  name {}_{}".format(i + 1, self.anchor2), file=fconf)
                customFunc = self.__custom_function(i,self.anchor2-1)
                print(customFunc, file=fconf)
    
                colvarList += str(i + 1) + "_" + str(self.anchor2) + " "
                centers += "0 "
                fconf.close()
                self.__collective_vari(space=2)

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
        print("  forceConstant {}".format(self.parameter.forceConst[colvar]), file=fconf)
        print("}", file=fconf)
        fconf.close()

    def __harmonicWalls(self, colvarList, centers, colvar=0):
        fconf = open(self.config_path, 'a')
        print("\n", file=fconf)
        print("harmonicWalls {", file=fconf)
        print("  colvars {}".format(colvarList), file=fconf)
        print("  lowerWalls {}".format(centers), file=fconf)
        print("  lowerWallConstant {}".format(self.parameter.forceConst[colvar]), file=fconf)
        print("}", file=fconf)
        fconf.close()

    def __grid_colvars(self):
        anchor1_coors = self.parameter.anchors[self.anchor1-1]
        anchor2_coors = self.parameter.anchors[self.anchor2-1]
         
        colvar_list, names = self.__get_colvars()
        count = 0
        #IF WE HAVE 3 anchors
        
        deltas = []
        length = []
        length_names = []

         #finding our delta values 
        '''
        for j in range(len(names)):
            values = []
            for i in range(len(self.parameter.anchors)):
                if self.parameter.anchors[i][j] not in values:
                    values.append(self.parameter.anchors[i][j])
            values.sort()
            deltas.append(abs(values[0] - values[1]))
        '''
        deltas = self.parameter.deltas.copy()
                
        for i in range(len(names)):
            if anchor1_coors[i] != anchor2_coors[i]:
                value = (anchor1_coors[i] + anchor2_coors[i])/2
                neighbor = [value - deltas[i], value, value + deltas[i]]
                neighbor.sort()
                if value == 0.0:
                    if 180 - deltas[i]/2 in [anchor1_coors[i],anchor2_coors[i]]:
                        value = 180
                        neighbor = [180 - deltas[i]/2, 180, -180 + deltas[i]/2]
                print(value, self.anchor1, self.anchor2)
                neighbor_name = names[i]
            else:
                length_values = [anchor1_coors[i], anchor1_coors[i] + deltas[i]/2, anchor1_coors[i] - deltas[i]/2]
                length_values.sort()
                length.append(length_values)
                length_names.append(names[i])
                
        #print(neighbor_name)
        equations = []
        #greater than, name first -- less than, name last
        #always parallel first
        equations.append(self.__create_step(neighbor_name, neighbor[2])) #parallel max
        for i in range(len(length)):
            equations.append(self.__create_step(length[i][2], length_names[i], length_names[i], length[i][0]))
        equations.append(self.__create_step(neighbor[0], neighbor_name))
        for i in range(len(length)):
            equations.append(self.__create_step(length[i][2], length_names[i], length_names[i], length[i][0]))
        for i in range(len(length)):
            equations.append(self.__create_step(length_names[i], length[i][2])) #max length (right upper)
            equations.append(self.__create_step(neighbor_name, neighbor[1], neighbor[2], neighbor_name))
            for k in range(len(length)):
                if i != k:
                    equations.append(self.__create_step(length_names[k], length[k][0], length[k][2], length_names[k]))
                
            equations.append(self.__create_step(length_names[i], length[i][2])) #right lower
            equations.append(self.__create_step(neighbor[1], neighbor_name, neighbor_name, neighbor[0]))            
            for k in range(len(length)):
                if i != k:
                    equations.append(self.__create_step(length_names[k], length[k][0], length[k][2], length_names[k]))
            
            equations.append(self.__create_step(length[i][0], length_names[i]))#min length (left upper)
            equations.append(self.__create_step(neighbor_name, neighbor[1], neighbor[2], neighbor_name))
            for k in range(len(length)):
                if i != k:
                    equations.append(self.__create_step(length_names[k], length[k][0], length[k][2], length_names[k]))
                
            equations.append(self.__create_step(length[i][0], length_names[i])) #left lower
            equations.append(self.__create_step(neighbor[1], neighbor_name, neighbor_name, neighbor[0]))
            for k in range(len(length)):
                if i != k:
                    equations.append(self.__create_step(length_names[k], length[k][0], length[k][2], length_names[k]))

        fconf = open(self.config_path, 'a')
        for eq in equations:
            count += 1
            eq = eq.replace('- -', '+ ')
            eq = eq.replace('--', '+')
            if self.parameter.grid_pbc in eq:
                if eq.count(self.parameter.grid_pbc) == 1:
                    if '180.0' in eq:
                        if 'step(-180.0-' + self.parameter.grid_pbc + ')' in eq:
                            eq = 'step(' + self.parameter.grid_pbc + ')*step(180.0-' + self.parameter.grid_pbc + ')'
                        elif 'step(' + self.parameter.grid_pbc + '-180.0' + ')' in eq:
                            eq = 'step(-' + self.parameter.grid_pbc + ')*step(' + self.parameter.grid_pbc + '+180.0)'  
                        elif 'step(180.0-' + self.parameter.grid_pbc + ')' in eq:
                            eq = 'step(' + self.parameter.grid_pbc + ')*step(180.0-' + self.parameter.grid_pbc + ')'
                        elif 'step(' + self.parameter.grid_pbc + '+180.0)' in eq:
                            eq = 'step(-' + self.parameter.grid_pbc + ')*step(' + self.parameter.grid_pbc + '+180.0)'  
                    elif str(180 - self.parameter.deltas[names.index(self.parameter.grid_pbc)]) in eq:
                        second_value = str(180 - self.parameter.deltas[names.index(self.parameter.grid_pbc)])
                        if '.' not in second_value:
                            second_value += '.0'
                        if 'step(-' + second_value + '-' + self.parameter.grid_pbc + ')' in eq:
                            eq = 'step(-' + self.parameter.grid_pbc + ')*' + eq
                        elif 'step(' + self.parameter.grid_pbc + '-' + second_value + ')' in eq:
                            eq = 'step(' + self.parameter.grid_pbc + ')*' + eq  
                        elif 'step(' + second_value + '-' + self.parameter.grid_pbc + ')' in eq:
                            eq = 'step(' + self.parameter.grid_pbc + ')*' + eq
                        elif 'step(' + self.parameter.grid_pbc + '+' + second_value + ')' in eq:
                            eq = 'step(-' + self.parameter.grid_pbc + ')*' + eq  
                    elif str(180 - self.parameter.deltas[names.index(self.parameter.grid_pbc)]/2) in eq:
                        second_value = str(180 - self.parameter.deltas[names.index(self.parameter.grid_pbc)]/2)
                        if '.' not in second_value:
                            second_value += '.0'
                        if 'step(-' + second_value + '-' + self.parameter.grid_pbc + ')' in eq:
                            eq = 'step(-' + self.parameter.grid_pbc + ')*' + eq
                        elif 'step(' + self.parameter.grid_pbc + '-' + second_value + ')' in eq:
                            eq = 'step(' + self.parameter.grid_pbc + ')*' + eq  
                        elif 'step(' + second_value + '-' + self.parameter.grid_pbc + ')' in eq:
                            eq = 'step(' + self.parameter.grid_pbc + ')*' + eq
                        elif 'step(' + self.parameter.grid_pbc + '+' + second_value + ')' in eq:
                            eq = 'step(-' + self.parameter.grid_pbc + ')*' + eq  

                        
            print('colvar {', file = fconf)
            print('  name ' + str(count), file = fconf)
            print('  customFunction ' + eq, file=fconf)
            for i in range(len(names)):
                if names[i] in eq:
                    for line in colvar_list[i]:
                        print("    " + line, file=fconf)
            print("}", file=fconf)   
        fconf.close()
            
    def __create_step(self, value1, value2, value3=None, value4=None):
        if value3 == None:
            return 'step(' + str(value1) + '-' + str(value2) + ')'
        else:
            return 'step(' + str(value1) + '-' + str(value2) + ')' + ' + ' + 'step(' + str(value3) + '-' + str(value4) + ')'
        
        
    def __grid_sampling(self):

        anchor1_coor = self.parameter.anchors[self.anchor1-1]
        anchor2_coor = self.parameter.anchors[self.anchor2-1]
        
        colvars, names = self.__get_colvars()
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
            print("\n", file=fconf)
            print("harmonicWalls {", file=fconf)
            print("  colvars length{}".format(count), file=fconf)
            print("  lowerWalls {}".format(walls[0]), file=fconf)
            print("  upperWalls {}".format(walls[1]), file=fconf)        
            print("  lowerWallConstant {}".format(self.parameter.forceConst[i]), file=fconf)
            print("  UpperWallConstant {}".format(self.parameter.forceConst[i]), file=fconf)
            print("} \n", file=fconf)
            
        fconf.close()  
        
    def __get_colvars(self):
        count = 0
        tmp = []
        colvar_list = []
        names = []
        with open(file = self.parameter.inputPath + '/colvar.txt') as f:
            for line in f:
                if line == '' or line == '\n':
                    continue
                tmp.append(line)
                if 'name' in line:
                    names.append(line.split()[1])
                if '{' in line:
                    count += 1
                if '}' in line:
                    count -= 1
                if count == 0:
                    colvar_list.append(tmp)
                    tmp = []
        return colvar_list, names
    
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
        if self.parameter.software == 'namd':
            state = path + "/stop.colvars.state"
        else:
            state = path + '/' + self.parameter.outputname + '.colvars.traj'
        if not os.path.isfile(state):
            return -1, [0, 0]  
        final_ms = [0, 0]
        
        # read state file generated at termination point 
        # smallest rmsd indicates new cell #
        if self.parameter.software == 'namd':
            RMSDs, lifetime = self.read_state(state)
        else:
            RMSDs, lifetime = self.read_traj(state)
            if RMSDs == None:
                return -1, [0,0]
        final_ms[0] = RMSDs.index(sorted(RMSDs)[0]) + 1
        
        if self.parameter.milestone_search == 2 or self.parameter.milestone_search == 3:
            if anchor: #this just means if it is seek
                if self.parameter.grid_pbc:
                    colvars, names = self.__get_colvars()
                colvar_number = len(self.parameter.anchors[0])
                anchor_values = self.parameter.anchors[anchor-1]
                count = 0
                values = []
                if self.parameter.software == 'namd':
                    for i in range(colvar_number):
                        if i == 0:
                            values.append(1)
                        else:
                            values.append(2)
                else:
                    values.append(colvar_number*2-1)
                #print(path)
                end_value = None
                if self.parameter.software == 'namd':
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
                else:
                    for i in range(len(RMSDs)):
                        if RMSDs[i] == values[0]:
                            end_value = i
                            break
                #print(end_value) 
                #total_equations = 3*colvar_number
                #end_value = 18
                #final = total_equations / end_value
                if self.parameter.software == 'namd':
                    new_end = int(end_value / colvar_number)
                else:
                    new_end = end_value + 1
                possible = []
                new_anchor = []
                for i in range(colvar_number):
                    a = [round(anchor_values[i] - self.parameter.deltas[i],5), round(anchor_values[i] + self.parameter.deltas[i],5)]
                    possible.append([min(a), i])
                    possible.append([max(a), i])
                    new_anchor.append(0)
                #print(possible)
                #print
                #print(new_anchor)
                #print(type(self.parameter.anchors))
                neighbor = possible[new_end-1][0]
                neighbor_name = possible[new_end-1][1]
                #print(neighbor_name)
                for i in range(colvar_number):
                    if neighbor_name == i:
                        new_anchor[i] = neighbor
                    else:
                        new_anchor[i] = anchor_values[i]
                final_ms = [anchor, 0]
                if self.parameter.software == 'gromacs':
                    for i in range(len(new_anchor)):
                       new_anchor[i] = round(new_anchor[i],5)
                #print(anchor_values)
                #print(new_anchor)
                if self.parameter.grid_pbc != False:
                        if names[neighbor_name] == self.parameter.grid_pbc:
                            if new_anchor[neighbor_name] > 180:
                                new_anchor[neighbor_name] = -180 + self.parameter.deltas[neighbor_name]/2
                            elif new_anchor[neighbor_name] < -180:
                                new_anchor[neighbor_name] = 180 - self.parameter.deltas[neighbor_name]/2
                print(new_anchor)
                for i in range(self.parameter.AnchorNum):
                    for j in range(colvar_number):
                        if round(self.parameter.anchors[i][j],5) == round(new_anchor[j],5):
                            if j == len(new_anchor) - 1:
                                final_ms[1]= i+1
                                break
                        else:
                            break
                final_ms.sort()
                if anchor not in final_ms:
                    return -1, [0,0]
                
            else:
                for i in range(len(self.parameter.anchors)):
                    for j in range(len(self.parameter.anchors[0])):
                        delta_values = self.parameter.deltas[j]/2
                        values = [self.parameter.anchors]
                if self.parameter.grid_pbc:
                   colvars, names = self.__get_colvars()                
                colvar_number = len(self.parameter.anchors[0])
                start_ms = pd.read_csv(path + '/start.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                start_anchors = [self.parameter.anchors[start_ms[0] - 1], self.parameter.anchors[start_ms[1] - 1]]
                count = 0
                values = []
                end_value = None
                if self.parameter.software == 'namd':
                    for i in range(colvar_number):
                        if i == 0:
                            values.append(1)
                        else:
                            values.append(2)
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
                else:
                    for i in range(len(RMSDs)):
                        if RMSDs[i] == (colvar_number*2) - 1:
                            end_value = (i+1)*colvar_number
                
                length = []
                #print(end_value)
                #print(start_anchors)
                for i in range(len(start_anchors[0])):
                    if round(start_anchors[0][i],5) != round(start_anchors[1][i],5):
                        if start_anchors[0][i] > start_anchors[1][i]:
                            upper_ms = 0
                            lower_ms = 1
                        else:
                            upper_ms = 1
                            lower_ms = 0
                        neighbor_values = [round(start_anchors[0][i],5), round(start_anchors[1][i],5)]
                        neighbor_values.sort()
                        neighbor_index = i
                        length.append("n")
                    else:
                        length.append(round(start_anchors[0][i],5))
                if end_value == None:
                    return -1, [0,0]
                #print(neighbor_index)
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
                    #print(path)
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
                

                if self.parameter.grid_pbc != False:                    
                    for i in range(len(names)):
                        if names[i] == self.parameter.grid_pbc:
                            if length[i] > 180:
                                length[i] = -180 + self.parameter.deltas[i]/2
                            if length[i] < -180:
                                length[i] = 180 - deltas[i]/2
                #print(final_ms, length)
                for i in range(self.parameter.AnchorNum):
                    for j in range(len(length)):
                        if round(self.parameter.anchors[i][j],5) == round(length[j],5):
                            if j == len(length) - 1:
                                final_ms[1]= i+1
                        else:
                            break
                final_ms.sort()
                if 0 in final_ms or final_ms[0] == final_ms[1]:
                    return -1, [0,0]            
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
        if float(lifetime) > 1.0:
            print(final_ms, lifetime, path)
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
        second_to_last = None
        with open(path) as r:
            for l in r:
                 line = l.split()
                 if line != []:
                    second_to_last = last_line
                    last_line = line
        #print(last_line)
        #last_line = last_line.split()
        #time = second_to_last.pop(0)
        time = last_line.pop(0)
        for i in range(len(self.parameter.anchors[0])):
            last_line.pop(-1)
        values = [abs(float(x)) for x in last_line]
        if self.parameter.milestone_search == 2:
            time = second_to_last.pop(0)
            for i in range(len(self.parameter.anchors[0])):
                second_to_last.pop(-1)
            values = [abs(float(x)) for x in second_to_last]
            
        return values, time
        

    
if __name__ == '__main__':
    from parameters import *
    new = parameters()
    new.restart = False
    new.initialize()
    new.milestone_search = 2
    print(new.forceConst)
    #new.deltas = [30,30,30]
    #colvar(new, anchor1=1, anchor2=2).generate()
    #colvar(new, anchor1=19, anchor2=20).generate()
    print(colvar(new, step='sample', anchor1=1, anchor2=2).generate())
    
