#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 13:39:31 2020

@author: allis
"""

from math import sqrt
import os


class plumed:
    def __init__(self, parameter, anchor1=None, anchor2=None,  
                 config_path=None, variables=None, step=None):
        self.parameter = parameter
        self.anchor1 = anchor1
        self.anchor2 = anchor2
        self.step = step
        self.variables = []
        #self.free = free
        self.colvars_number = len(self.parameter.anchors[0])
        for i in range(1, self.colvars_number + 1):
            self.variables.append("")
        #self.initial = initial
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
        elif self.step == 'seek':
            if self.parameter.milestone_search in (2,3):
                self.__plumed_seek_grid()
            else:
                self.__rmsd()
        else:
            if self.parameter.milestone_search == 2:
                self.__plumed_grid()
            else:
                self.__rmsd()
        
        

    def __sampling_grid(self):
        '''
        Description: Used to generate sampling files when using grid
        '''    
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


    def __plumed_seek_grid(self):
        '''
        Description: Used to generate the plumed file used when using milestone_search = 3
        Arguments: None
        Returns:None
        '''
        equations = []
        borders = []
        count = 0
        names = self.__get_names()

        #print(self.parameter.deltas)
        anchor_coordinates = self.parameter.anchors[self.anchor1 - 1]
        for i in range(len(anchor_coordinates)):
            borders.append([anchor_coordinates[i] - (self.parameter.deltas[i]/2), anchor_coordinates[i] + 
                           (self.parameter.deltas[i]/2)])
        for i in borders:
            i.sort()
        for i in range(len(borders)):
            equations.append((self.__create_step(borders[i][0], names[i])))
            for k in range(len(borders)):
                if  k == i:
                    continue
                equations.append((self.__create_step(borders[k][1], names[k], names[k], borders[k][0])))
            equations.append((self.__create_step(names[i], borders[i][1])))
            for k in range(len(borders)):
                if  k == i:
                    continue
                equations.append((self.__create_step(borders[k][1], names[k], names[k], borders[k][0])) + '+')
                
        combined_equations = []
        for eq in range(1,len(equations)+1,2):
            combined_equations.append(equations[eq-1] + '+' + equations[eq])
        fconf = open(self.config_path, 'a')
        count = 0
        vars = []
        args='PRINT STRIDE=1 ARG='
        args_list=[]
        for i in range(len(names)):
            vars.append('x' + str(i+1))
        for eq in combined_equations:
            count += 1
            args_list.append('stepfunction' + str(count))
            eq = eq.replace(' ', '')
            eq = eq.replace('--', '+')
            eq = eq.replace('CV', 'x')
            if eq[-1] == '+':
                eq = eq[slice(0,-1)]
            print("CUSTOM ...", file=fconf)
            print("LABEL=stepfunction{}".format(count),file=fconf)
            print("PERIODIC=NO", file=fconf)
            print("ARG={}".format(",".join(names)), file=fconf)
            print("VAR={}".format(",".join(vars)), file=fconf)
            print("FUNC={}".format(eq), file=fconf)
            print("... CUSTOM\n", file=fconf)  
        count = 0
        lower_basin = 2*(len(names)-1) + 0.9
        for i in range(2*len(names)):
            committor = self.__committor('stepfunction'+str(i+1),str(lower_basin),
                                         str(lower_basin+0.2),'-10','-9',str(i+1))
            for i in committor:
                print(i,file=fconf)
        print(args + ",".join(args_list) + ',' + ','.join(names) + " FILE=" + str(self.parameter.outputname) + '.colvars.traj', file=fconf)
        fconf.close()
    
    def __plumed_grid(self):
        '''
        This creates the plumed file for the free trajectory step when using Grid rather than Voronoi
        Arguments: None
        Returns: None
        '''
        anchor1_coors = self.parameter.anchors[self.anchor1-1]
        anchor2_coors = self.parameter.anchors[self.anchor2-1]
         
        count = 0
        length = []
        length_names = []
        names = self.__get_names()
        
        full_function = []
        deltas = self.parameter.deltas.copy()
        
        for i in range(len(names)):
            if anchor1_coors[i] != anchor2_coors[i]:
                value = (anchor1_coors[i] + anchor2_coors[i])/2
                neighbor = [round(value - deltas[i],5), round(value,5), round(value + deltas[i],5)]
                neighbor.sort()
                neighbor_name = names[i]
            else:
                length_values = [round(anchor1_coors[i],5), round(anchor1_coors[i] + deltas[i]/2,5), 
                                 round(anchor1_coors[i] - deltas[i]/2,5)]
                length_values.sort()
                length.append(length_values)
                length_names.append(names[i])
        equations = []
        #greater than, name first -- less than, name last
        
        '''
        The first two before the loop are for the two parallel milestones
        No matter what dimension we are working with, there will always be two parallel
        Based on my notation, these ones use the same "neighbor", and are just shifted up and down,
        changing the length
        '''
        tmp = ''
        tmp += self.__create_step(neighbor_name, neighbor[2]) + '+'
        for i in range(len(length)):
            tmp += self.__create_step(length[i][2], length_names[i], 
                                      length_names[i], length[i][0])
        equations.append(tmp)
        tmp= ''
        tmp += self.__create_step(neighbor[0], neighbor_name) + '+'

        for i in range(len(length)):
            tmp += self.__create_step(length[i][2], length_names[i], 
                                      length_names[i], length[i][0]) + '+'
        equations.append(tmp)
        '''
        Now we enter the long for loop
        This is for every one that is perpendicular to our original milestone
        The amount of times this loops through is equal to our number of colvars minus one
        So if we have two colvars, this loops through one time.
        It basically assigns a main "length" that becomes our new "neighbor".
        Our neighbor then becomes like a length and we have another loop that is repeated within this big loop
        These "k" loops cycle through all of our original "lengths", skipping over our "main length" (when i=k)
        '''
        for i in range(len(length)):
            tmp = ''
            tmp += self.__create_step(length_names[i], length[i][2]) + '+' #max length (right upper)
            tmp += (self.__create_step(neighbor_name, neighbor[1], 
                                       neighbor[2], neighbor_name)) + '+'
            for k in range(len(length)):
                if i != k:
                    tmp += (self.__create_step(length_names[k], length[k][0], 
                                               length[k][2], length_names[k])) + '+'
            equations.append(tmp)  
            tmp = ''
            tmp += (self.__create_step(length_names[i], length[i][2])) + '+' #right lower
            tmp += (self.__create_step(neighbor[1], neighbor_name, 
                                       neighbor_name, neighbor[0])) + '+'    
            for k in range(len(length)):
                if i != k:
                    tmp += (self.__create_step(length_names[k], length[k][0],
                                               length[k][2], length_names[k])) + '+'
            equations.append(tmp)
            tmp=''
            tmp += (self.__create_step(length[i][0], length_names[i])) + '+' #min length (left upper)
            tmp += (self.__create_step(neighbor_name, neighbor[1], 
                                       neighbor[2], neighbor_name)) + '+'
            for k in range(len(length)):
                if i != k:
                    tmp(self.__create_step(length_names[k], length[k][0],
                                           length[k][2], length_names[k])) + '+'
            equations.append(tmp)
            tmp = ''
            tmp += self.__create_step(length[i][0], length_names[i]) + '+' #left lower
            tmp += (self.__create_step(neighbor[1], neighbor_name, 
                                       neighbor_name, neighbor[0])) + '+'
            for k in range(len(length)):
                if i != k:
                    tmp += (self.__create_step(length_names[k], length[k][0], 
                                               length[k][2], length_names[k])) + '+'
            equations.append(tmp)
        
        '''
        Now we loop through our total list of equations and write everything to our 
        configuration file, which is plumed_free.dat
        '''
        vars = []
        for i in range(len(names)):
            vars.append('x' + str(i+1))
        fconf = open(self.config_path, 'a')
        for eq in equations:
            count += 1
            eq = eq.replace(' ', '')
            eq = eq.replace('--', '+')
            eq = eq.replace('CV','x')
            if eq[-1] == '+':
                eq = eq[slice(0,-1)]
            print("CUSTOM ...", file=fconf)
            print("LABEL=stepfunction{}".format(count),file=fconf)
            print("PERIODIC=NO",file=fconf)
            print("ARG={}".format(",".join(names)), file=fconf)
            print("VAR={}".format(",".join(vars)), file=fconf)
            print("FUNC={}".format(eq), file=fconf)
            print("... CUSTOM\n", file=fconf)
        print('\n', file=fconf)
        args = []
        step_amount = 2*(len(names)-1) + 0.9
        for i in range(len(names)*2 + 2):
            args.append('stepfunction' + str(i+1))
            committor = self.__committor('stepfunction'+str(i+1),str(step_amount),
                                         str(step_amount+0.2),'-10','-9',str(i+1))
            for i in committor:
                print(i,file=fconf)
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
    new.initialize()
    new.restart = False
    #milestones(new).initialize(status=1)
    new.iteration=1
    new.milestone_search = 2
    #plumed(new, anchor1=1,anchor2=2,step='sample').generate_plumed()
    plumed(new, anchor1=19, anchor2=20, step='sample').final_ms_plumed(self.parameter.crdPath + '/19_20')
    
    

        
        
        
        
        
