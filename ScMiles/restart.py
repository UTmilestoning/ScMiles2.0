# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:41:42 2020

@author: allis
"""
import glob
from milestones import *

def read_log(parameter, status):
    import os
    fullLog = []
    if parameter.milestone_search == 1 or parameter.milestone_search == 3:
        seek = False
    else:
        seek = True
    if status == 1:
        seek = True
        sample = True
    else:
        sample = False
    #uses log to find latest iteration
    with open(file = os.path.join(parameter.currentPath, 'log')) as r:
        for line in r:
            line = line.rstrip()
            info = line.split()
            info = info[2:]
            if 'Iteration' in line and 'complete' not in line and 'created' not in line:
                parameter.iteration = int(info[2])
            elif 'Reactant and product are connected' in line:
                #We are done seeking
                seek = True
            elif 'Ready to launch free trajectories.' in line:
                #We are done sampling
                sample = True
            lastLog = (" ".join(str(x) for x in info))
    if parameter.method == 0 or parameter.iteration == 0:
        parameter.iteration = 1
    if (seek == True and sample == False) or (sample == False and parameter.method == 0) or ("Sampling on each" in lastLog):
        for folder in glob.glob(parameter.crdPath + '/*'):
            #print(folder)
            if os.path.exists(folder + '/' + parameter.outputname + '.colvars.state') or os.path.exists(folder + '/restarts'):
                parameter.finished_constain.add('MS' + folder.split('/')[-1])
    return seek, sample, lastLog

            

if __name__ == '__main__':
    from parameters import *
    from namd_conf_custom import *
    new = parameters()
    new.initialize()
    restart_simulation(new, 0)
    print(new.iteration)
    
      
