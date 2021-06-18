#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 13:39:31 2020

@author: allis
"""

from math import sqrt
import os
from colvar import *


class plumed:
    def __init__(self, parameter, anchor1=None, anchor2=None,  
                 config_path=None, step=None):
        self.parameter = parameter
        self.anchor1 = anchor1
        self.anchor2 = anchor2
        self.step = step
        #self.free = free
        self.colvars_number = len(self.parameter.anchors[0])
        self.config_path = self.parameter.ScMilesPath + "/plumed_free.dat"  if self.step != 'sample' else self.parameter.ScMilesPath + "/plumed.dat"


    def __exit__(self, exc_type, exc_value, traceback):
        return 

    def __repr__(self) -> str:
        return ('Plumed generator')  
    
    def generate_plumed(self):
        '''
        Description: This dictates how the plumed file is created based on where we are in the simulation 
            (seek, sampling, or free trajectories) and user input (Grid vs Voronoi Cell)
        Arguments: None
        Returns: None
        '''
        plumed = []
        if self.step == 'seek': 
            outputFrequency = 1
        else:
            outputFrequency = self.parameter.colvarsTrajFrequency
        #if self.step == 'sample':
        filename = '/plumed.dat'
        #else:
        #    filename = '/plumed_free.dat'
        with open(self.parameter.inputPath + filename) as f:
            for line in f:
                plumed.append(line)
        fconf = open(self.config_path, 'w+')
        for line in plumed:
            if line == "":
                continue
            line = line.replace('\n','')
            print(line, file=fconf)
        fconf.close()
        if self.step == 'sample':
            if self.parameter.milestone_search in (0,1):
                self.__sampling_voronoi()
            else:
                self.__sampling_grid()
        else:
            if self.parameter.milestone_search in (2,3):
                self.__plumed_grid()
            else:
                self.__rmsd()

        

    def __sampling_grid(self):
        '''
        Description: Used to generate sampling files when using grid
        '''  
        #this is just finding the coordinates of our two anchors using the anchor numbers, and the corresponding 
        #values in anchor.txt
        anchor1_coor = self.parameter.anchors[self.anchor1-1]
        anchor2_coor = self.parameter.anchors[self.anchor2-1]
       
        names = self.__get_names()
        for i in range(len(anchor1_coor)):
            if anchor1_coor[i] != anchor2_coor[i]:
                neighbor = i
        delta_value = self.parameter.deltas.copy()
        center = (anchor1_coor[neighbor] + anchor2_coor[neighbor])/2
        
        fconf = open(self.config_path, 'a')
        print('nei: RESTRAINT ...', file=fconf)
        print('ARG=' + names[neighbor], file=fconf)
        print('AT=' + str(center), file=fconf)
        print('KAPPA=' + str(self.parameter.neighbor_kappa), file=fconf)
        print('... nei:', file=fconf)

        for i in range(len(names)):
            if i == neighbor:
                continue
            print('uwall: UPPER_WALLS ARG=' + str(names[i]) + ' AT=' + str(round(anchor1_coor[i] + self.parameter.deltas[i]/2,5)) + ' KAPPA=' + str(self.parameter.neighbor_kappa), file=fconf)
            print('lwall: LOWER_WALLS ARG=' + str(names[i]) + ' AT=' + str(round(anchor1_coor[i] - self.parameter.deltas[i]/2,5)) + ' KAPPA=' + str(self.parameter.neighbor_kappa), file=fconf)
        print('\n',file=fconf)
        print('PRINT STRIDE=1 ARG=' + ','.join(names) + ' FILE=' + self.parameter.outputname + '.colvars.traj', file=fconf)
        fconf.close()

    def __plumed_grid(self):
        '''
        This creates the plumed file for the free trajectory step when using Grid rather than Voronoi
        Arguments: None
        Returns: None
        '''
        count = 0  
        equations, ms_network = colvar(self.parameter, self.anchor1, self.anchor2, step=self.step).get_equations()
        variables = []
        print(ms_network)
        names = self.__get_names()
        for i in range(len(names)):
            variables.append('x' + str(i+1))
        fconf = open(self.config_path, 'a')
        args = []
        for eq in equations:
            count += 1
            eq[0] = eq[0].replace(' ', '')
            eq[0] = eq[0].replace('--', '+')
            eq[0] = eq[0].replace('CV','x')
            if eq[0][-1] == '+':
                eq[0] = eq[0][slice(0,-1)]
            print("CUSTOM ...", file=fconf)
            print("LABEL=stepfunction{}".format(eq[1]),file=fconf)
            print("PERIODIC=NO",file=fconf)
            print("ARG={}".format(",".join(names)), file=fconf)
            print("VAR={}".format(",".join(variables)), file=fconf)
            print("FUNC={}".format(eq[0]), file=fconf)
            print("... CUSTOM\n", file=fconf)
            args.append(eq[1])
        print('\n', file=fconf)
        step_amount = 2*(len(self.parameter.anchors[0])) -0.1
        for i in range(len(args)):
            committor = self.__committor('stepfunction'+str(args[i]),str(step_amount),
                                         str(step_amount+0.2),'-10','-9',str(i+1))
            for i in committor:
                print(i,file=fconf)
        for i in range(len(args)):
           args[i] = 'stepfunction' + args[i]
        print("PRINT STRIDE=1 ARG=" + ",".join(args) + ',' + ",".join(names) + ' FILE=' + self.parameter.outputname + '.colvars.traj', file=fconf)
        fconf.close()
        
        
    def __committor(self, name, basin_ll1, basin_ul1, basin_ll2, basin_ul2, file):
        '''
        This creates the PLUMED COMMITTOR function
        Arguments:
            name: what the function is labeled
            basin_ll1: Lower limit for first basin
            basin_ul1: Upper limit for first basin
            basin_ll2: Lower limit for second basin (may be an arbitrary value/placeholder, the function requires
                                                     two basins)
            basin_ul2: Upper limit for second basin (may be an arbitrary value/placeholder, the function requires
                                                     two basins)
        Returns: 
            Entire committor function as a list where each element is a string
        '''
        full_committor = []
        full_committor.append("COMMITTOR ...")
        full_committor.append("STRIDE=1")
        full_committor.append("ARG={}".format(name))
        full_committor.append("BASIN_LL1={}".format(basin_ll1))
        full_committor.append("BASIN_UL1={}".format(basin_ul1))
        full_committor.append("BASIN_LL2={}".format(basin_ll2))
        full_committor.append("BASIN_UL2={}".format(basin_ul2))
        full_committor.append("FILE=committor" + str(file) + ".dat")
        full_committor.append("... COMMITTOR\n")  
        return full_committor

    
    def __create_step(self, value1, value2, value3=None, value4=None):
        '''
        Description: This creates our step functions used in our plumed_free files
        Arguments:
            value1:This is a necessary value. It is the first value in one of our step functions
            value2: Second alue in the step function
            value3: This is used if we have a compound step function
            value4: This is used if we have a compound step function
        '''
        if value3 == None:
            return 'step(' + str(value1) + '-' + str(value2) + ')'
        else:
            return 'step(' + str(value1) + '-' + str(value2) + ')' + \
                ' + ' + 'step(' + str(value3) + '-' + str(value4) + ')'
                        
    
    def __rmsd(self):
        '''
        Description: Creates the plumed files used in seek and free trajectories when
            using Voronoi Cell
        Arguments: None
        Returns: None, but prints the text to self.config_path (located in the same folder as the python files,
            which is later copied over to the appropriate file) The config_path when using
            rmsd is always plumed_free.dat
        '''
        names = self.__get_names()
        full_function = []
        walls=[]
        '''
        for i in range(len(names)):
            full_function.append("CUSTOM...")
            full_function.append("LABEL=cos" + str(names[i]))
            full_function.append("ARG=" + names[i])
            full_function.append("VAR=x1")
            full_function.append("FUNC=cos(x1)")
            full_function.append("PERIODIC=-pi,pi")
            full_function.append("...CUSTOM\n")
            walls.append("cos"+str(names[i]))
        for i in walls:
            full_function.append("LOWER_WALLS ARG=" + i + " AT=" + " KAPPA= ")
        '''
        full_args = []
        arguments = ",".join(names)
        for an in range(len(self.parameter.anchors)):
            full_args.append('rmsd' + str(an+1))
            full_function.append("CUSTOM ...")
            full_function.append("LABEL=" + 'rmsd' + str(an+1))
            arg = ".".join(names)
            arg = 'ARG=' + arguments
            var = 'VAR='            
            for i in range(len(names)):
                var += 'x' + str(i+1)
                if i != len(names)-1:
                    var += ','   
            full_function.append(arg)
            full_function.append(var)
            if len(names) == 1:
                custom_function = 'abs(x1-' + str(self.parameter.anchors[an][0]) + ')'
            else:
                custom_function=""
                for i in range(len(names)):
                    custom_function += '(x' + str(i+1) + '-' + str(self.parameter.anchors[an][i]) + ')^2'
                    if i != len(names) - 1:
                        custom_function += '+'   
                custom_function = 'sqrt(' + custom_function + ')'
            full_function.append("FUNC=" + custom_function)
            if self.parameter.pbc:
                full_function.append("PERIODIC=0,6.2832")
            else:
                full_function.append("PERIODIC=NO")
            full_function.append("... CUSTOM\n")
        for an in range(len(self.parameter.anchors)):
            if self.anchor1 == an + 1 or (self.anchor2 == an + 1 and self.step != 'seek'):
                continue
            full_function.append("CUSTOM ...")
            full_function.append("LABEL=stepvalue" + str(an+1))
            if self.step == 'seek':
                full_function.append('ARG='+ 'rmsd' + str(self.anchor1) + ',rmsd' + str(an+1))
                full_function.append('VAR=a,b')
                full_function.append('FUNC=step(a-b)')
                low = '0.9'
                high = '1.1'
            else:
                full_function.append('VAR=a,b,c')
                full_function.append('ARG='+ 'rmsd' + str(self.anchor1) +  \
                                     ',rmsd' + str(self.anchor2) + ',rmsd' + str(an+1))
                full_function.append('FUNC=step(a-c)+step(b-c)')
                #low = '0.9'
                low = '1.9'
                high = '2.1'                

            full_function.append('PERIODIC=NO') #maybe change this
            full_function.append("... CUSTOM\n")
        for an in range(len(self.parameter.anchors)):
            if self.anchor1 == an + 1 or (self.anchor2 == an + 1 and self.step == 'free'):
                continue 
            committor = self.__committor('stepvalue' + str(an+1), str(low), 
                                         str(high), str(-10), str(-9), str(an+1))
            for i in committor:
                full_function.append(i)
        for i in names:
            full_args.append(i)
        full_args = ",".join(full_args)
        full_function.append("PRINT STRIDE={} ARG={} FILE={}.colvars.traj"
                             .format(1, full_args, self.parameter.outputname))
        fconf = open(self.config_path, 'a')
        for i in full_function:
            i=i.replace('--','+')
            print(i,file=fconf)
        fconf.close()
                    
    def __get_names(self):
        '''
        Description: Retrieves the names of our variables from the input file plumed.txt
        Arguments: None
        Returns: Our variables names as strings in a list, in order they are presented in the file
        '''
        names = []
        with open(self.parameter.inputPath + '/plumed.dat') as f:
            for line in f:
                if 'ATOMS' in line and 'CV' in line:
                    line = line.split()
                    line[0] = line[0].replace(':','')
                    names.append(line[0])
        if self.parameter.CV_suffixes == []:
            return names
        for i in range(len(self.parameter.CV_suffixes)):
            if self.parameter.CV_suffixes[i] == 'None':
                continue
            names[i] += '.' + self.parameter.CV_suffixes[i]
        return names        
    
    def __sampling_voronoi(self):
        names = self.__get_names()
        centers = 0
        arg = ",".join(names)
        var = []
        label=[]
        centers = []
        kappa = []
        for i in range(len(names)):
            var.append('x' + str(i+1))
        var = ",".join(var)
        full_function = []
        full_function.append("CUSTOM ...")
        full_function.append("LABEL=neighbor")
        full_function.append("ARG=" + arg)
        full_function.append("VAR=" + var)
        full_function.append("FUNC=" + self.__custom_function(self.anchor1-1,
                                                              self.anchor2-1))
        full_function.append("PERIODIC=NO")
        full_function.append("... CUSTOM\n")
        for i in range(self.parameter.AnchorNum):
            if i+1 == self.anchor1 or i+1 == self.anchor2:
                continue
            for j in [self.anchor1, self.anchor2]:
                label.append(str(i+1) + '_' + str(j))
                full_function.append("CUSTOM ...")
                full_function.append("LABEL=" + str(i+1) + '_' + str(j))
                full_function.append("ARG=" + arg)
                full_function.append("VAR=" + var)
                full_function.append("FUNC=" + self.__custom_function(i, j-1))
                full_function.append("PERIODIC=NO")
                full_function.append("... CUSTOM\n")
                centers.append('0.0')
                kappa.append(str(self.parameter.walls_kappa))
        full_function.append("nei: RESTRAINT ...")
        full_function.append("ARG=neighbor")
        full_function.append("AT=0.0")
        full_function.append("KAPPA=" + str(self.parameter.neighbor_kappa))
        full_function.append("... nei:")
        full_function.append("walls: LOWER_WALLS ...")
        full_function.append("ARG=" + ",".join(label))
        full_function.append("AT=" + ",".join(centers))
        full_function.append("KAPPA=" + ",".join(kappa))
        full_function.append("... walls:\n")
        full_function.append("PRINT STRIDE=1 ARG=" + ",".join(names) + ",neighbor," + ",".join(label) + " FILE=" + self.parameter.outputname + '.colvars.traj')
        fconf = open(self.config_path, 'a')
        for i in full_function:
            print(i,file=fconf)
        fconf.close()

    def __custom_function(self, anchor1, anchor2):
        '''
        Description: Creates the customFunction for cases with more than one collective variable
        Arguments:
            anchor1: The index of the first anchor we use to compute
            anchor2: The index of the second anchor we use to compute
            Note: These indexes are are one less than the anchor we are going for, because
            we are indexing through self.parameter.anchor, where we start at index 0
            so, for example, the anchor 1 is at index 0, or the anchor 45 is at index 44.
            So this is the index, not the actual anchor
        Returns: The custom function as a string to be printed in the file
        
        '''
        customFunc = ""
        names = self.__get_names()
        for i in range(0,2):
            if i == 0:
                anchor = anchor1
            else:
                anchor = anchor2
            customFunc = customFunc + 'sqrt('
            for j in range(self.colvars_number):
                customFunc = customFunc + '(x' + str(j+1) + '-' + \
                    str(self.parameter.anchors[anchor][j]) + ')^2'
                if j != self.colvars_number - 1:
                    customFunc = customFunc + '+'
            if i == 0:
                customFunc = customFunc + ')-'
            else:
                customFunc = customFunc + ')'
        customFunc=customFunc.replace('--','+')
        return customFunc
        

    def __find_rmsd(self, anchor1, anchor2, newanchor):
        '''
        Description: This calculates the rmsd of our two original anchors and sees how
            far each is from our new anchor
        Arguments:
            anchor1: This is the first anchor in our starting milestone
            anchor2: This is the second anchor in our starting milestone
            newanchor:This is the new anchor that was hit when running trajectories
        Returns: The anchor from our two starting anchors that is closer to the new anchor
        Example: If we started at milestone 1_3 and ended up hitting milestone 2,
            we would need to see if our ending milestone would be 1_2 or 2_3. We compute the
            rmsd from anchor 1 to 2 annd from anchors 2 to 3, and see which one is smaller. 
            So if the rmsd from anchor 1 to 2 was smaller, this would return 1, which would be
            the one that is part of our final_ms
        '''
        rmsd1 = 0
        rmsd2 = 0
        anchor1_values = self.parameter.anchors[anchor1 - 1]
        anchor2_values = self.parameter.anchors[anchor2 - 1]
        newanchor_values = self.parameter.anchors[newanchor - 1]
        if len(self.parameter.anchors[0]) == 1:
            rmsd1 = abs(anchor1_values[0] - newanchor_values[0])
            rmsd2 = abs(anchor2_values[0] - newanchor_values[0])
        else:
            for i in range(len(anchor1_values)):
                rmsd1 += (anchor1_values[i] - newanchor_values[i])**2
                rmsd2 += (anchor2_values[i] - newanchor_values[i])**2
            rmsd1 = sqrt(rmsd1)
            rmsd2 = sqrt(rmsd2)
        if rmsd1 < rmsd2:
            return anchor1
        else:
            return anchor2
                    
                    
            
if __name__ == '__main__':
    from parameters import *
    from milestones import *
    new = parameters()
    new.restart = False
    new.initialize()
    new.deltas = [30,30]
    #new.grid_ignore= ['13_20','20_21','20_27','19_20','18_19','17_18']
    new.milestone_search = 2
    new.all_grid = milestones(new).grid_ms()
    #print(new.all_grid)
    new.MS_list = milestones(new).initialize()
    #plumed(new, anchor1=1,anchor2=2,step='sample').generate_plumed()
    plumed(new, anchor1=5, anchor2=13, step='free').generate_plumed()
    

        
        
        
        
        
