# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:41:42 2020

@author: allis
"""

def read_log(parameter, status):
    import os
    fullLog = []
    if parameter.milestone_search == 1:
        seek = False
    else:
        seek = True
    if status == 1:
        sample = True
    else:
        sample = False
    #uses log to find latest iteration
    with open(file = os.path.join(parameter.currentPath, 'log')) as r:
        for line in r:
            line = line.rstrip()
            info = line.split()
            info = info[2:]
            if 'Iteration' in line and 'complete' not in line:
                parameter.iteration = int(info[2])
            elif 'Reactant and product are connected' in line:
                #We are done seeking
                seek = True
            elif 'Ready to launch free trajectories.' in line:
                #We are done sampling
                sample = True
            fullLog.append(" ".join(str(x) for x in info))
    try:
        lastLog = fullLog[-1]
    except:
        print('log file is empty. Cannot restart simulation')
    if parameter.method == 0 or parameter.iteration == 0:
        parameter.iteration = 1     
    return seek, sample, lastLog


if __name__ == '__main__':
    from parameters import *
    from namd_conf_custom import *
    new = parameters()
    new.initialize()
    restart_simulation(new, 0)
    print(new.iteration)
    
      
