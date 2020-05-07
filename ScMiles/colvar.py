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


class colvar:
    def __init__(self, parameter, anchor1=None, anchor2=None, 
                 free=None, initial=None, 
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
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.input_dir = self.parent_path + '/my_project_input'
        self.config_path = self.path + "/colvar_free.conf" if self.free == 'yes' else self.path + "/colvar.conf"


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
        with open(file=self.input_dir + '/colvar.txt') as f:
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
        with open(file=self.input_dir+'/colvar.txt') as f:
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
        with open(file=self.input_dir+'/colvar.txt') as f:
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
#                    if 'anchor' in line:
#                        line = line.replace("anchor", '('+str(self.parameter.anchors[anchor-1][0])+')')
                    if 'anchor' in line:
                        for i in range(0, self.colvars_number):
                            line = line.replace("anchor", '('+str(self.parameter.anchors[anchor-1][i])+')', 1)
                    tmp.append(line + '\n')
                
        fconf = open(self.config_path, 'a')
        for line in tmp:
            print("  " * space + "  " + line, file=fconf)
        fconf.close()
        
#    def __rmsd_to_anchor(self, anchor, coeff=None, space=0):
#        '''
#        Change this function for different 1D case.
#        Follow the format in Colvars to define the distance measurement.
#        For the commands below, it generates the following.
#        
#        colvar {
#          name rmsd1
#          customFunction abs(psi - (-165.0))
#          dihedral {
#            name psi
#            group1 atomNumbers 7
#            group2 atomNumbers 9
#            group3 atomNumbers 15
#            group4 atomNumbers 17
#          }
#        }
#          
#        '''
#        fconf = open(self.config_path, 'a')
#        name = "rmsd" + str(anchor)
#        print("\n" + "  " * space + "colvar {", file=fconf)
#        print("  " * space + "  name {:5}".format(name), file=fconf)
#        func = "abs(psi - (" + str(self.parameter.anchors[anchor-1][0]) + "))"
#        print("  " * space + "  customFunction {}".format(func), file=fconf)
#        fconf.close()
#        self.__collective_vari_1()
#        fconf = open(self.config_path, 'a')
#        print("  " * space + "}", file=fconf)
#        fconf.close()

    def generate(self):    
        '''This is the main function that generates colvars '''
        # scriptPath = os.path.dirname(os.path.abspath(__file__))
        if self.initial == 'yes': 
            outputFrequency = 1
        else:
            outputFrequency = self.parameter.colvarsTrajFrequency
            
        fconf = open(self.config_path, 'w+')
        print("colvarsTrajFrequency      {}".format(outputFrequency), file=fconf)
        print("colvarsRestartFrequency	 {}".format(self.parameter.colvarsRestartFrequency), file=fconf)
        if self.free == 'yes':
            print("scriptedColvarForces on", file=fconf)
        if self.parameter.customColvars == 1:
            print("", file=fconf)
            with open(file=self.input_dir + '/custom.colvar') as f_custom:
                for line in f_custom:
                    print(line, file=fconf)
        fconf.close()
        
        if self.free == 'yes':
            for i in range(self.parameter.AnchorNum):
                self.__rmsd_to_anchor(i+1)       
        else:
            self.__get_colvar_names()
            if self.colvars_number == 1:
                self.__constraint1D1()
                self.__harmonic1D()
            else:
                self.__constraint2D1()
                colvarList, centers = self.__constraint2D2()
                self.__harmonic2D()
                self.__harmonicWalls(colvarList, centers)



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
        print("  forceConstant {}".format(self.parameter.forceConst), file=fconf)
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
    
    def __harmonic2D(self):
        fconf = open(self.config_path, 'a')
        print("harmonic {", file=fconf)
        print("  colvars neighbor", file=fconf)
        center = 0
        print("  centers {}".format(str(center)), file=fconf)
        print("  forceConstant {}".format(self.parameter.forceConst), file=fconf)
        print("}", file=fconf)
        fconf.close()

    def __harmonicWalls(self, colvarList, centers):
        fconf = open(self.config_path, 'a')
        print("\n", file=fconf)
        print("harmonicWalls {", file=fconf)
        print("  colvars {}".format(colvarList), file=fconf)
        print("  lowerWalls {}".format(centers), file=fconf)
        print("  lowerWallConstant {}".format(self.parameter.forceConst), file=fconf)
        print("}", file=fconf)
        fconf.close()


if __name__ == '__main__':
    from parameters import *
    new = parameters()
    new.initialize()
    print(new.anchors)
    colvar(new, anchor1=1, anchor2=2).generate()
    colvar(new, anchor1=1, anchor2=2, free='yes').generate()
    
