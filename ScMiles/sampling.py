#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:06:01 2019

@author: Wei Wei

sampling class:
This is a class that is used to generate and submit sampling files, as well as 
    move sampling files to the /restart folder.
When initialized, we accept two arguments, parameter and jobs
    parameter is the instance of our parameter class that has all of the input data stored
    jobs is an instance of the 'run' class that we will use to run the submit scripts

"""

import time
import re
import os
import glob
from shutil import move
from colvar import *
from parameters import *
from log import *
from network_check import *
from log import log
from run import *


class sampling:
    def __init__(self, parameter, jobs):
        self.parameter = parameter
        self.jobs = jobs

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        return
    
    def __repr__(self) -> str:
        return ('Sampling on milestones.')    

    def constrain_to_ms(self):
        '''
        Description: This generates and our sampling files, and submits the jobs
        Parameters: None, this class is initialized with parameter and jobs (which is an instance of the run class)
        Returns: None
        '''
        import os
        MS_list = self.parameter.MS_list.copy()
        finished = self.parameter.Finished.copy()
        sleep = False
        script_milestones = []
        count = 0
        for name in MS_list:
            count += 1
            submit = False
            if name in finished or name in self.parameter.finished_constain:
                continue
            #if not self.jobs.check(MSname=name):
            #    continue   
            lst = re.findall('\d+', name)
            anchor1 = int(lst[0])
            anchor2 = int(lst[1])
            script_milestones.append(str(anchor1) + '_' + str(anchor2))
            if len(script_milestones) == self.parameter.traj_per_script[1]:
                submit = script_milestones
            elif len(MS_list) == count:
                submit = script_milestones
            restartsPath = self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2) + '/restarts'
            create_folder(self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2))
            if self.parameter.restart == True:
                restart_path = self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2)
                if self.parameter.software == 'namd' and os.path.exists(restart_path + '/' + self.parameter.outputname + '.colvars.state'): 
                    continue
                elif self.parameter.software == 'gromacs' and os.path.exists(restart_path + '/' + self.parameter.outputname + '.colvars.traj'):
                    self.parameter.finished_constain.add(name)
                    #move(restart_path + '/submit',restart_path + '/submit_finished')
                    continue
            if not os.path.exists(restartsPath) or self.parameter.additional_sampling == True:
                if self.parameter.software == 'namd':
                    colvar(self.parameter, anchor1, anchor2).generate()
                else:
                    plumed(self.parameter, anchor1, anchor2).generate_plumed()
                if self.parameter.additional_sampling == True:
                    new_sampling = 0
                    for folder in glob.glob(restartsPath + '/*'):
                        try:
                            folder_num = int(folder.split('.')[-2])
                        except:
                            continue
                        if folder_num > new_sampling:
                            new_sampling = folder_num
                    self.jobs.prepare_trajectory(anchor1,anchor2,script=submit,new_sampling=new_sampling)
                else:
                    self.jobs.prepare_trajectory(anchor1,anchor2,script=submit)
            if submit is not False:
                script_milestones = []
        for name in MS_list:
            lst = re.findall('\d+', name)
            anchor1 = int(lst[0])
            anchor2 = int(lst[1])
            submit_path = self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2) + '/submit'
            if os.path.exists(submit_path):
                self.jobs.submit(submit_path)
                #move(submit_path, self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2) + '/submit_submitted')
                sleep = True
        log("{} milestones identified.".format(str(len(MS_list))))             
        if sleep == True:
            log("Sampling on each milestone...")  
            time.sleep(60)
        #self.parameter.restart = False

    def check_sampling(self):
        '''
        Description: check the progress of our sampling
        Parameters: None
        Returns: self.parameter.Finished, which is a list of all of the milestones that have finished sampling
        '''
        import re
        from datetime import datetime
        finished = set()
        MS_list = self.parameter.MS_list.copy()
        while True: #we stay in this this loop until all of our sampling jobs are done running
            for name in self.parameter.MS_list:
                MSname = 'MS'
                if not self.jobs.check(SampleName=MSname):
                    continue
                elif name in self.parameter.finished_constain:
                    continue
                elif name in finished:
                    continue
                else:
                    [anchor1, anchor2] = list(map(int,(re.findall('\d+', name))))
                    finished.add(name)
                    if name not in self.parameter.finished_constain:
                        log("Finished sampling on milestone ({},{}).".format(anchor1, anchor2))  
                        self.parameter.finished_constain.add(name)
                            
            if finished | self.parameter.finished_constain == MS_list:
                self.parameter.Finished = finished.copy()
                log("Finished sampling on all milestones.")    
                self.move_restart(self.parameter.MS_list)
                log("Ready to launch free trajectories.") 
                return self.parameter.Finished
            print("Next check in 600 seconds. {}".format(str(datetime.now())))
            time.sleep(600)   # 600 seconds

    def move_restart(self, names):
        '''
        Description: moves all files with the extension .coor, .vel, or .ext to a folder, restarts and
            renames our submit file to 'submit_finished' to show that sampling on that script has finished
        Parameters: names, which is just a copy of our MS_list
        Returns: None
        '''

        for name in names:
            [anchor1, anchor2] = list(map(int,(re.findall('\d+', name))))
            ms = str(anchor1) + '_' + str(anchor2)
            filePath = self.parameter.crdPath + '/' + ms
            restartFolder = filePath + '/restarts'
            create_folder(restartFolder)
            #failed_sampling = glob.glob(filePath + '/core.*')
            #if len(failed_sampling) != 0:
            #    print('A NAMD error possibly occured in the folder {}'.format(filePath))
            if os.path.exists(filePath + '/submit'):
                move(filePath + '/submit',filePath + '/submit_finished')
            for ext in ["coor", "vel", "xsc"]:
                for file in glob.glob(filePath + '/*.' + ext):
                    try:
                        move(file, restartFolder)
                    except:
                        continue             
            
if __name__ == '__main__':
    from parameters import *
    from run import *
    new = parameters()
    new.initialize()
    new.restart = False
    new.MS_list = milestones(new).read_milestone_folder()
    jobs = run(new)
    test = sampling(new, jobs)
    test.constrain_to_ms()
#    print(new.anchors)
#    test.generate()
