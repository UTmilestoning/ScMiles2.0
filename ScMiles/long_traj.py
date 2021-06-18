# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 10:55:52 2021

@author: allis
"""
from additional_functions import *
import pandas as pd
import numpy as np
import os
from parameters import *
from compute import *
from math import sqrt

parentDirectory = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
inputPath = os.path.join(parentDirectory, 'my_project_input')
outputPath = os.path.join(parentDirectory, 'my_project_output')
currentPath = os.path.join(outputPath, 'current')

create_folder(outputPath)
create_folder(currentPath)

raw_data = pd.read_csv(inputPath + "/path_history.dat", header=None, delimiter=r"\s+").values.tolist()
milestones = []
anchors = [0,0]
times = [0,0]
transitions = []
lifetimes = dict()
ms_index = dict()

previous_milestone = str(int(raw_data[0][1]) + 1) + '_' + str(int(raw_data[1][1] + 1))
times[0] = (raw_data[1][0] * 1000000)
for i in range(len(raw_data)):
    anchors[0] = anchors[1]
    anchors[1] = int(raw_data[i][1] + 1)
    times[1] = (raw_data[i][0] * 1000000)
    ms_int_form = sorted(anchors)
    ms = str(ms_int_form[0]) + '_' + str(ms_int_form[1])
    if 0 in ms_int_form:
        continue
    if ms_int_form not in milestones:
        milestones.append(ms_int_form)
    if ms != previous_milestone:
        transitions.append([previous_milestone, ms, round(times[1]-times[0],5)])
        times[0] = times[1]
        previous_milestone = ms
        
milestones = sorted(milestones)
for i in range(len(milestones)):
    milestones[i] = str(milestones[i][0]) + '_' + str(milestones[i][1])
    lifetimes[milestones[i]] = []
matrix = np.zeros((len(milestones),len(milestones)))
for i in range(len(transitions)):
    matrix[milestones.index(transitions[i][0]),milestones.index(transitions[i][1])] += int(1)
    lifetimes[transitions[i][0]].append(transitions[i][2])
printing_lifetimes = [[],[],[]]
for i in range(len(milestones)):
    printing_lifetimes[0].append(milestones[i])
    printing_lifetimes[1].append(np.mean(lifetimes[milestones[i]]))
    if len(lifetimes[milestones[i]]) > 1:
        printing_lifetimes[2].append((np.std(lifetimes[milestones[i]], ddof=1))/sqrt(len(lifetimes[milestones[i]])))
    else:
        printing_lifetimes[2].append(0)

remove_index = []
if k_min_sum:
    for i in range(len(milestones)):
        total = 0
        for j in range(len(milestones)):
            total += matrix[j][i]
        if total <= k_min_sum:
            remove_index.append(i)

if remove_index:
    remove_index.reverse()
    for i in remove_index:
        for j in range(len(printing_lifetimes)):
            printing_lifetimes[j].pop(i)
        np.delete(matrix,i,0)
        np.delete(matrix,i,1)
        milestones.pop(i)

with open(currentPath + '/life_time.txt', 'w+') as f:  
    print(''.join(['{:>10}'.format(item) for item in printing_lifetimes[0]]),file=f)
    print('',file=f)
    print('\n'.join([''.join(['{:10.2f}'.format(item) for item in np.squeeze(printing_lifetimes[1])])]),file=f)
    print('\n'.join([''.join(['{:10.2f}'.format(item) for item in np.squeeze(printing_lifetimes[2])])]),file=f)
    
with open(currentPath + '/k.txt', 'w+') as f:
    print(''.join(['{:>10}'.format(item) for item in printing_lifetimes[0]]),file=f)
    print('',file=f)
    print('\n'.join([''.join(['{:10d}'.format(int(item)) for item in row])for row in matrix]),file=f)
  
k_ave = k_average(matrix)
with open(outputPath + '/k_norm.txt','w+') as f:
    f.write('\n'.join([''.join(['{:10.5f}'.format(item) for item in row])for row in k_ave]))

for i in range(len(milestones)):
    [anchor1, anchor2] = get_anchors(milestones[i])
    ms_index[i] = sorted([anchor1,anchor2])
np.save(currentPath + '/ms_index.npy', ms_index)


new = parameters()
new.reactant = [15,20]
new.product = [6,11]
compute(new)
    
