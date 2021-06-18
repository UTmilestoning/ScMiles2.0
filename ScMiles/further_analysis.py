# -*- coding: utf-8 -*-
"""
Created on Tue May  4 10:33:35 2021

@author: allis
"""

#further analysis
import os
import networkx as nx
from compute import *
from milestoning_mp import *
from parameters import *
from network_check import *
from math import sqrt



class analysis:
    def __init__(self, parameter, k_cutoff=None, max_lifetime=None, custom_k=None, max_flux=None):
        self.parameter = parameter
        self.k_min_sum = None
        self.max_lifetime = None
        self.data_file = False
        self.source = False
        self.iteration = 1
        self.reactant = None
        self.product = None
        self.max_flux = False

    def read_analysis(self):
        if not os.path.isfile(self.parameter.inputPath + '/analysis.txt'): 
            print('No file called analysis.txt in my_project_input. Unable to do analysis')
            return
        self.parameter.initialize()

        with open(self.parameter.inputPath + '/analysis.txt') as f:
            for line in f:
                line = line.split()
                if line == []:
                    continue
                if 'source' in line[0]:
                    if line[1].lower() == 'scmiles':
                        self.source = 'scmiles'
                    elif line[1].lower() == 'path_history.dat': 
                        self.source = 'path_history'
                    elif line[1].lower == 'custom':
                        self.source = 'custom'
                    elif '.colvars.traj' in line[1].lower():
                        self.source = 'colvar_file'
                if 'k_min_sum' in line[0]:
                    self.k_min_sum = int(line[1])
                    self.parameter.k_min_sum = int(line[1])
                if 'ignore_milestones' in line:
                    self.parameter.analysis_ignore_milestones = line[1].split(',')
                if 'ignore_transitions' in line: 
                    transitions = line[1].split(',')
                    for i in range(len(transitions)):
                       transitions[i] = transitions[i].split('-')
                    self.parameter.ignore_transitions = transitions 
                if 'max_lifetime' in line[0]:
                    self.parameter.max_lifetime = float(line[1])
                if 'all_iterations' in line[0]:
                    self.parameter.running_average = True
                if 'data_file' in line[0]:
                    self.data_file=True
                if 'iteration' in line[0]:
                    self.iteration = int(line[1])
                    self.parameter.iteration = int(line[1])
                if 'reactant' in line[0]:
                    r = line[1].split(',')
                    for i in range(len(r)):
                        r[i] = int(r[i])
                    self.parameter.reactant = r
                if 'product' in line[0]:
                    r = line[1].split(',')
                    for i in range(len(r)):
                        r[i] = int(r[i])
                    self.parameter.product = r  
                if 'max_flux' in line:
                    self.max_flux = True
        if self.data_file == True:
            self.make_data_file()
        if self.source == 'path_history':
            self.long_traj()
        elif self.source == 'scmiles':
            self.parameter.MS_list = milestones(self.parameter).initialize(status=1)
            print(self.parameter.MS_list)
            if self.parameter.analysis_ignore_milestones:
                for i in self.parameter.MS_list:
                    if i in self.parameter.analysis_ignore_milestones:
                        self.parameter.MS_list.remove(i)
            #network_check(self.parameter.MS_list)
            print(self.parameter.MS_list)
            print(self.parameter.k_min_sum)
            milestoning(self.parameter, False)
        elif self.source == 'colvar_file':
            self.colvar_file()
        #if self.max_flux == True:
        #    self.max_flux_path()
                   
    def find_last_iteration(self, path):
        iteration = 1
        while True:
            if os.path.exists(path + '/' + str(iteration)):
                iteration += 1
            else:
                return iteration
            
    def make_data_file(self):
        import pandas as pd
        import os
        import re
        data = []
        for ms in self.parameter.MS_list:
            lst = re.findall('\d+',ms)
            name = lst[0] + '_' + lst[1]
            ms_path = self.parameter.crdPath + '/' + name
            iteration = self.find_last_iteration(ms_path)
            next_frame = get_next_frame_num(ms_path)
            for i in range(1, iteration):
                for traj in range(1, next_frame):
                    tmp = []
                    traj_path = ms_path + str(i) + '/' + str(traj) 
                    if os.path.isfile(traj_path + '/start.txt'):
                        start = pd.read_csv(traj_path + '/start.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                    else:
                        start = ['N','N']
                    if os.path.isfile(traj_path + '/end.txt'):
                        end = pd.read_csv(traj_path + '/end.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                    else:
                        end = ['N','N']
                    if os.path.isfile(traj_path + '/lifetime.txt'):
                        lifetime = pd.read_csv(traj_path + '/lifetime.txt', header=None, delimiter=r'\s+').values.tolist()[0]
                    else:
                        lifetime = 'N'
                    if os.path.isfile(traj_path + '/enhanced'):
                        enhanced = 'Enhanced'
                    else:
                        enhanced = 'NotEnhanced'
                    data.append([parameter.iteration, start[0], start[1], end[0], end[1], lifetime[0], enhanced, traj_path])
        with open(parameter.currentPath + '/all_data.txt', 'w+') as f1:
            for item in data:
                f1.write(" ".join(map(str,item)) + '\n')
                
    def colvar_file(self):
        raw_data = pd.read_csv(self.parameter.inputPath + "/analysis.colvars.traj", header=1, delimiter=r"\s+").values.tolist()
        anchors = pd.read_csv(self.parameter.inputPath + "/anchors.txt", header=None, delimiter=r"\s+").values.tolist()       
        anchor_number = len(anchors)
        current_anchors = [0,0]
        transitions = []
        milestones = []
        lifetime = 0
        if self.parameter.software == 'gromacs':
            for i in range(len(raw_data)):
                raw_data[i][0] = raw_data[i][0]*1000
                
        for line in raw_data:
            all_rmsd = [0] * anchor_number
            for an in range(anchor_number):
                rmsd = 0
                for colvar in range(1,len(line)):
                    rmsd += (line[colvar] - anchors[an][colvar-1])**2
                rmsd = sqrt(rmsd)
                all_rmsd[an] = rmsd
            anchor = all_rmsd.index(min(all_rmsd)) + 1
            all_rmsd[anchor-1] = max(all_rmsd)
            second_anchor = all_rmsd.index(min(all_rmsd)) + 1
            if anchor in current_anchors:
                continue
            else:
                current_time = line[0] - lifetime
                if 0 in current_anchors:
                    current_anchors = [anchor,second_anchor]
                    continue
                previous_ms = sorted(current_anchors)
                milestones.append(str(previous_ms[0]) + '_' + str(previous_ms[1]))
                new_ms = sorted([anchor, second_anchor]) 
                transitions.append([str(previous_ms[0]) + '_' + str(previous_ms[1]), str(new_ms[0]) + '_' + str(new_ms[1]),current_time])
                current_anchors = new_ms
                lifetime = current_time
        self.create_files(milestones, transitions)
        
    def long_traj(self):
        create_folder(self.parameter.outputPath)
        create_folder(self.parameter.currentPath)
        raw_data = pd.read_csv(self.parameter.inputPath + "/path_history.dat", header=None, delimiter=r"\s+").values.tolist()
        milestones = []
        anchors = [0,0]
        times = [0,0]
        transitions = []
        
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
        self.create_files(milestones, transitions)
        
    def create_files(self, milestones, transitions):
        lifetimes = dict()
        ms_index = dict()        
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
        milestones_to_delete = []
        remove_index = []
        print(len(matrix[0]))
        if self.k_min_sum:
            for i in range(len(milestones)):
                total = 0
                for j in range(len(milestones)):
                    total += matrix[j][i]
                if total <= self.k_min_sum:
                    remove_index.append(i)
                    milestones_to_delete.append(milestones[i])  
        if remove_index:
            remove_index.reverse()
            for i in remove_index:
                for j in range(len(printing_lifetimes)):
                    printing_lifetimes[j].pop(i)
                matrix = np.delete(matrix,i,0)
                matrix = np.delete(matrix,i,1)
                milestones.pop(i)
        with open(self.parameter.currentPath + '/life_time.txt', 'w+') as f:  
            print(''.join(['{:>10}'.format(item) for item in printing_lifetimes[0]]),file=f)
            print('',file=f)
            print('\n'.join([''.join(['{:10.2f}'.format(item) for item in np.squeeze(printing_lifetimes[1])])]),file=f)
            print('\n'.join([''.join(['{:10.2f}'.format(item) for item in np.squeeze(printing_lifetimes[2])])]),file=f)
            
        with open(self.parameter.currentPath + '/k.txt', 'w+') as f:
            print(''.join(['{:>10}'.format(item) for item in printing_lifetimes[0]]),file=f)
            print('',file=f)
            print('\n'.join([''.join(['{:10d}'.format(int(item)) for item in row])for row in matrix]),file=f)
          
        k_ave = k_average(matrix)
        with open(self.parameter.currentPath + '/k_norm.txt','w+') as f:
            f.write('\n'.join([''.join(['{:10.5f}'.format(item) for item in row])for row in k_ave]))
        
        for i in range(len(milestones)):
            [anchor1, anchor2] = get_anchors(milestones[i])
            ms_index[i] = sorted([anchor1,anchor2])
            
        np.save(self.parameter.currentPath + '/ms_index.npy', ms_index)
        k, index, q = compute(parameter)
        G = nx.Graph()
        for i in transitions:
            if i[0] in milestones_to_delete or i[1] in milestones_to_delete:
                continue
            G.add_node(i[0])
            G.add_edge(i[0],i[1])
        nx.draw(G, with_labels=True)

         

if __name__ == '__main__':
    from milestoning_mp import *
    from milestones import *

    parameter = parameters()
    analyze = analysis(parameter)
    analyze.read_analysis()
    
    
