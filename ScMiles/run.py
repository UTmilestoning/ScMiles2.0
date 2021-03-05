
#!/usr/bin/env python3 -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 16:33:01 2020

@author: allis
"""
import os, time
from shutil import copy
from datetime import datetime
import re
import subprocess
import glob
from shutil import move
from shutil import copyfile
from fileinput import FileInput
from additional_functions import *
from colvar import *
from plumed import *
from log import *
#from milestones import *
from namd_conf_custom import *
        
class launch:
    def __init__(self, parameter, step, loop1=None, loop2=None, colvar_file=None,
                 traj_per_script=None, start_frame=None):
        self.parameter = parameter
        self.step = step
        self.resubmit = False
        self.finished = set()
        
        if self.parameter.software == 'namd':
            if self.step == 'seek' or self.step == 'free':
                self.colvar_file = 'colvar_free.conf'
                self.config_file = 'free.namd'
            elif step == 'sample':
                self.colvar_file = 'colvar.conf'
                self.config_file = 'sample.namd'
        elif self.parameter.software == 'gromacs':
            if self.step == 'seek' or self.step == 'free':
                self.colvar_file = 'plumed_free.dat'
                self.config_file = 'free.mdp'
            elif self.step == 'sample':
                self.colvar_file = 'plumed.dat'
                self.config_file = 'sample.mdp'      
        elif self.parameter.software == 'lammps':
            if self.step == 'seek' or self.step == 'free':
                self.colvar_file = 'plumed_free.dat'
                self.config_file = 'free.in'
            elif self.step == 'sample':
                self.colvar_file = 'plumed.dat'
                self.config_file = 'sample.in'
        if self.step == 'free':
            self.time = 200
            self.script_name = 'MILES'
            self.loop1 = self.parameter.MS_list.copy()
            self.loop2 = range(1, self.parameter.trajPerLaunch + 1)
            self.traj_per_script = self.parameter.traj_per_script[2]
            self.total_trajectories = self.parameter.trajPerLaunch
        elif self.step == 'sample':
            self.time = 600
            self.script_name = 'MS'
            self.loop1 = self.parameter.MS_list.copy()
            self.loop2 = range(0,1)
            self.traj_per_script = self.parameter.traj_per_script[1]
            self.total_trajectories = None
        elif self.step == 'seek':
            self.time = 60
            self.script_name = 'a'
            self.loop1 = range(1, self.parameter.AnchorNum + 1)
            self.loop2 = range(1, self.parameter.initial_traj + 1)
            self.total_trajectories = self.parameter.initial_traj
            self.traj_per_script = self.parameter.traj_per_script[0]
            if self.parameter.software == 'namd':
                self.colvar_file = 'colvar_free.conf'
            elif self.parameter.software == 'gromacs':
                self.colvar_file = 'plumed_free.dat'
        self.start_frame = {}

    def launch_trajectories(self, possible_ms=None, lastLog=None):
        '''
        Here is the main workflow of submitting jobs
        After each step a statement is logged
        Arguments:
            possible_ms: This is used with grid. It is used if new milestones are found.
               Only new milestones that are possible based on the original anchors that were provided
               by the user in anchor.txt
            lastLog: This is used if restart is used. It tells the last thing the script did
               before stoppinng, so we know where we need to start
        '''
        log("Beginning next step: {}".format(self.step))
        if lastLog != None:
            if self.step == 'sample':
                pass    
            else:
                self.resubmit = True
                lastLog = None
        print_log = self.__prepare_trajectory_files()
        if print_log == True:
            log("All files prepared for {} step".format(self.step))      
        print_log = self.__submit_trajectories()
        if print_log == True:
            log("All trajectories launched")  
            self.__check_progress()
            log("Trajectories finished running")
        if self.step != 'sample':
            count = self.__post_launch(possible_ms)
            log("Post launch analysis of trajectories completed")
            return count
        elif self.step == 'sample': # and self.parameter.software == 'namd':
            log('Restarts moved to restart folders')
            self.move_files()
        log("Finished launching {} trajectories. Moving to next step".format(self.step))
        self.parameter.restart = False

    def __start_frames(self):
        '''
        This finds the starting frame depending on the method, whether or not we are restarting, and the step 
        '''
        for i in self.loop1:
            path = self.find_parent_path(i)
            next_frame = get_next_frame_num(path)
            if self.step == 'sample':
                self.start_frame[i] = 0
            elif self.parameter.method == 1 and self.step == 'free':
                self.start_frame[i] = 1
            elif self.resubmit == False: # or self.step == 'seek':
                self.start_frame[i] = next_frame
            else:
                self.start_frame[i] = next_frame - self.total_trajectories
            if self.start_frame[i] <= 1:
                self.start_frame[i] = 1
            self.parameter.restart = False
            #print(i, self.start_frame[i])


            
    def __prepare_trajectory_files(self):
        '''
        This prepares all of our trajectory folders and files
        '''
        print_log = False
        self.__start_frames()
        script_milestones = []
        for i in self.loop1:
            if self.step == 'seek' and self.parameter.active_anchors != None:
                if i not in self.parameter.active_anchors:
                    continue
            #if self.step == 'free' and i in self.parameter.skip_MS:
            #    continue
            count = 0
            for j in self.loop2:
                trajectory_number = 0
                [anchor1, anchor2] = get_anchors(i)
                count += 1
                submit = False
                if self.step != 'sample':
                    trajectory_number = self.start_frame[i] + j - 1
                #print(self.start_frame)
                path = self.find_parent_path(i, trajectory_number) 
                if os.path.isfile(path + '/submit'):
                    path_num = 0
                    found = False
                    while found == False:
                        path_num += 1
                        if not os.path.isfile(path + '/submit_done' + str(path_num)):
                            move(path + '/submit', path + '/submit_done' + str(path_num))
                            found = True        
                finished = self.__check_finished(path)
                if finished == True:
                    continue
                if self.step == 'sample':
                    script_milestones.append(str(anchor1) + '_' + str(anchor2))
                else:
                    script_milestones.append(trajectory_number)
                if len(script_milestones) == self.traj_per_script:
                    submit = script_milestones
                elif self.step == 'sample' and len(self.parameter.MS_list) == count:
                    submit = script_milestones
                elif count == self.total_trajectories and self.step != 'sample':
                    submit = script_milestones
                if submit:
                    script_milestones = []
                if self.step == 'free':
                    frame = self.parameter.startTraj + self.start_frame[i] // self.parameter.trajPerLaunch + j * self.parameter.interval
                    if frame > self.parameter.nframe:
                        break
                else:
                    frame = None
                create_folder(path)
                if not os.path.isfile(path + '/start.txt') and self.step == 'free':
                    with open(path + '/start.txt', 'w+') as f1:
                        print(anchor1, anchor2, file = f1)
                #if os.path.isfile(os.path.join(path, self.colvar_file)): #maybe add a more sophisticated check finished here
                #    continue
                if self.step == 'sample' and self.parameter.additional_sampling == True and self.parameter.software == 'namd':
                    new_sampling = self.__new_sampling(path)
                else:
                    new_sampling = None
                if self.parameter.software == 'namd':
                    colvar(self.parameter, step=self.step, anchor1=anchor1, anchor2=anchor2).generate()
                elif self.parameter.software == 'gromacs' or self.parameter.software == 'lampps':
                    plumed(self.parameter, step=self.step, anchor1=anchor1, anchor2=anchor2).generate_plumed()    
                if not os.path.isfile(path + '/submit') and submit != False:
                    self.__prepare_script(anchor1, anchor2, trajectory_number, submit, frame)
                if not os.path.isfile(path + '/' + str(self.config_file)) or new_sampling is not None:
                    if self.parameter.software == 'namd':
                        print_log = True
                        self.__prepare_namd(anchor1, anchor2, trajectory_number, submit, frame, new_sampling)
                    elif self.parameter.software == 'gromacs':
                        print_log = True
                        copy(self.parameter.inputPath + '/' + self.config_file, 
                             path + '/' + str(self.config_file))
                        if self.parameter.iteration > 1:
                            self.__gromacs_exact(path)
                copy(self.parameter.ScMilesPath + '/' + str(self.colvar_file), 
                     path + '/' + str(self.colvar_file))
        return print_log

    def __check_finished(self, path):
        '''
        When using restart, this checks a trajectory to see if it has already been done or not.
        Arguments: path: This is the path we are looking at
        Returns: True if our trajectory has finished running, False if it needs to be run again
        '''
        if self.parameter.software == 'namd':
            if os.path.isfile(path + '/stop.colvars.state'):
                return True
            elif os.path.isfile(path + '/1.log'):
                with open(path + '/1.log') as f:
                    for line in f:
                        if 'End of program' in line:
                            return True
            if self.step == 'sample' and self.parameter.additional_sampling == False and os.path.exists(path + '/restarts'):
                return True
            return False
        if self.parameter.software == 'gromacs':
            if self.step == 'sample' and os.path.isfile(path + '/state.cpt'):
                return True
            else:
                for i in range(self.parameter.AnchorNum):
                    if os.path.isfile(path + '/committor' + str(i+1) + '.dat'):
                        return True
                return False
                            
    def __new_sampling(self, path):
        '''
        If the user wants to sample on milestones to get even more restart files
        than they got originally, they can use new_sampling. This finds our largest
        restart file so we can use that as the starting point.
        '''
        new_sampling = 0
        for folder in glob.glob(path + '/restarts/*'):
            try:
                folder_num = int(folder.split('.')[-2])
            except:
                continue
            if folder_num > new_sampling:
                new_sampling = folder_num
        return new_sampling
            
    def __prepare_script(self, a1, a2, trajectory_number, script=None, frame=None):
        '''modify job submission file'''
        path = self.find_parent_path([a1,a2], trajectory_number)
        if self.step != 'sample':
            script.sort()
            submit_path = self.find_parent_path([a1,a2], str(script[0]))
            start_script = str(script[0])
            end_script = str(script[-1])
            if end_script != start_script:
                script_number = str(start_script) + '_' + str(end_script)
            else:
                script_number = str(start_script)
        else:
            submit_path = path
        if self.step != 'seek':
            name = self.script_name + '_' + str(a1) + '_' + str(a2)
            if self.step == 'free':
                name += '_' + script_number       
        else:
            name = self.script_name + '_' + str(a1) + '_' + script_number
            
        create_folder(path)
        copyfile(self.parameter.inputPath + '/submit', submit_path +'/submit')
        
        if self.parameter.software == 'namd':
            self.__script_namd(submit_path +'/submit', name, script, path)
        else:
            self.__script_gromacs(submit_path, name, script, a1, frame)
            
    def __script_gromacs(self, newScriptPath, name, script, a1, snapshot):
        count = 0
        ndx = ""
        gromacs_lines = []
        frame_value = []
        if self.step == 'seek':
            frame_value = [0]
            relative_input = '../../../../my_project_input'
            gro_option = ' -c '
            config = ''
            if self.parameter.ndx_file is not None:
                ndx = ' -n ' + relative_input + '/' + self.parameter.ndx_file
        elif self.step == 'free':
            if self.parameter.iteration == 1:
                for i in script:
                    frame_value.append(round(self.parameter.gromacs_timestep * (i + self.parameter.startTraj - 1),2))
            else:
                frame_value = [0]
            #print(frame_value)
            relative_input = '../../../../my_project_input'
            gro_option = ' -t'
            if self.parameter.iteration == 1:
                config = '-time frameval -c ' + relative_input + '/conf.gro'     
            else:
                config = '-c ' + relative_input + '/conf.gro'     
            if self.parameter.ndx_file is not None:
                ndx = ' -n ' + relative_input + '/' + self.parameter.ndx_file        
        else:
            frame_value = [0]
            relative_input = '../../my_project_input'
            #if os.path.exists(newScriptPath + '/statefromseek.cpt'):
            gro_option = ' -t'
            #else:
            #    gro_option = ' -c'
            config = ' -c ' + relative_input + '/conf.gro'
            if self.parameter.ndx_file is not None:
                ndx = ' -n' + relative_input + '/conf.gro'
            
        if self.step != 'free':
            if os.path.exists(self.parameter.seekPath) and self.step == 'sample':
                new_path = 'statefromseek.cpt'
            else:
                new_path = relative_input + '/pdb/' + str(a1) + '.pdb'
        else:
            if self.parameter.iteration == 1:
                new_path = '../../traj.trr'
            else:
                new_path = self.parameter.outputname + '.cpt'
                
        with FileInput(files=newScriptPath + '/submit', inplace=True) as f:
            for line in f:
                line = line.strip()
                info = line.split()
                if not info:
                    continue
                if 'name' in line:
                    line = line.replace('name',name)
                if 'path' in line:
                    line = line.replace('path', newScriptPath)
                if '/bin/' in line and 'bash' not in line:
                    count += 1
                    if count == 1:    
                        line += ' grompp -f ' + self.config_file + gro_option + ' ' + str(new_path) + ' ' + config + ' -p ' + relative_input + '/topol.top' + ndx + ' -o topol.tpr'
                        gromacs_lines.append(line)
                        continue
                    else:
                        if 'mdrun' not in line:
                            line += ' mdrun'
                        line += ' -s topol.tpr -plumed ' + self.colvar_file + ' '
                        gromacs_lines.append(line)
                        continue
                print(line)
        newScript = newScriptPath + '/submit'
        fconf = open(newScript, 'a')
        for i in range(len(script)):
            print('cd ../' + str(script[i]), file=fconf)
            for j in gromacs_lines:  
                if 'frameval' in j:
                    j = j.replace('frameval', str(round(frame_value[i],2)))
                fconf.write(j + '\n')
        fconf.close()

    def __script_namd(self, newScript, name, script, path):
        with FileInput(files=newScript, inplace=True) as f:
            for line in f:
                line = line.strip()
                info = line.split()
                if not info:
                    continue
                keyword = None
                for item in ['source', 'name', 'path', 'namd']:
                    if item in info:
                        keyword = item
                        place = info.index(item)
                        break
                if keyword == None:
                    print(" ".join(str(x) for x in info))
                    continue
                if keyword == 'source':
                    if self.parameter.nodes:
                        import numpy as np
                        rand = np.random.uniform(0,len(self.parameter.nodes),1)
                        info[2] = 'hostname="' + self.parameter.nodes[int(rand)] + '"'
                    else:
                        info.insert(0, '#')  
                elif keyword == 'name':
                    info[place] = name
                elif keyword == 'path':
                    info[place] = path
                elif keyword == 'namd':
                    info[place] = self.config_file
                    namd_line = (" ".join(str(x) for x in info))
                    continue
                line = (" ".join(str(x) for x in info))
                print(line)
        fconf = open(newScript, 'a')
        for i in script:
            print('cd ../' + str(i), file=fconf)
            print(namd_line, file=fconf)
        fconf.close()

    def __gromacs_exact(self, path):
        if os.path.isfile(path + '/enhanced'):
            return
        with FileInput(files=path + '/' + str(self.config_file), inplace=True) as f:
            for line in f:
                if 'gen_vel' in line:
                    info = line.split()
                    info[-1] = 'no'
                    line = " ".join(info)
                print(line)
        
        
    def __prepare_namd(self, a1=None, a2=None, trajectory_number=None, script=None, frame=None, new_sampling=None):
        '''modify namd configuration file'''
        from fileinput import FileInput
        from random import randrange as rand 
        import re
        
        enhanced = 0
        
        template = self.parameter.inputPath + '/' + self.config_file
        MSpath = self.find_parent_path([a1, a2], trajectory_number)
        newNamd = MSpath + '/' + self.config_file
        copyfile(template, newNamd)
        if os.path.exists(MSpath + '/enhanced') and self.step == 'free':
            enhanced = 1
            
        tmp = []
        colvar_commands = False
        with open(newNamd, 'r') as f:
            for line in f:
#                line = line.lower()
#                 info = line.split("#")[0].split()
                info = line.split()
                if info == []:
                    continue                
                if "colvars" in info and "on" in info:
                    colvar_commands = True
                if "colvarsConfig" in info and colvar_commands:
                    info[1] = self.colvar_file
                    l = " ".join(str(x) for x in info)+"\n"
                    tmp.append(l)
                    continue
                
                if "run" in info or 'minimize' in info:
#                    if info[0] == '#':
#                        continue
                    if not colvar_commands:
                        tmp.append('colvars on\n')
                        info[0] = 'colvarsConfig'
                        info[1] = self.colvar_file + '\n\n'
                        l = " ".join(str(x) for x in info)
                        tmp.append(l)
                        colvar_commands = True
                    if self.step == 'sample':
                        filename = None
                    elif self.parameter.milestone_search == 2:
                        filename = '/tclScript_grid.txt'
                    elif self.step == 'seek':
                        filename = '/tclScript_seek.txt'
                    elif self.step == 'free':
                        filename = '/tclScript_step2.txt'
                    if filename is not None:      
                        with open(file=self.parameter.ScMilesPath + '/' + filename) as f_tcl:
                            for l in f_tcl:
                                if "qsub" in l:
                                    kill = l.strip()
                                    killswitch = kill.split()
                                    if self.parameter.jobsubmit == "qsub":
                                        killswitch.pop(0)
                                    else:
                                        killswitch[0] = '#'
                                    a = " ".join(str(x) for x in killswitch)
                                    tmp.append(a +'\n')   
                                elif "sbatch" in l:
                                    kill = l.strip()
                                    killswitch = kill.split()
                                    if self.parameter.jobsubmit == "sbatch":
                                        killswitch.pop(0)
                                    else:
                                        killswitch[0] = '#'
                                    a = " ".join(str(x) for x in killswitch)
                                    tmp.append(a +'\n')      
                                else:
                                    tmp.append(l)
                        tmp.append('\n')
                tmp.append(line)     
                
        with open(newNamd, 'w') as f:
            for i in range(len(tmp)):
                f.write(tmp[i])
                
        if self.parameter.namd_conf == True:
            if self.step != 'free' and (self.step == 'seek' or self.parameter.milestone_search == 0 or self.parameter.milestone_search == 2):
                namd_conf_mod(self.parameter.inputPath, newNamd, a1)
            
        with FileInput(files=newNamd, inplace=True) as f:
            for line in f:
                line = line.strip()
#                line = line.lower()
                info = line.split()
                
                if "coordinates" in line and "bincoordinates" not in line.lower():
                    info[1] = self.parameter.inputPath + '/pdb/' + str(a1) + ".pdb"
                    if self.step == 'sample' and (os.path.exists(self.parameter.seekPath)) and self.parameter.iteration == 0:
                        info[1] = "./seek.ms.pdb"
                        
                if "outputname" in line:
                    info[1] = self.parameter.outputname
                    
                if "seed" in line:
                    info[1] = rand(10000000, 99999999)
                if 'restartfreq' in line:
                    if self.step == 'seek':
                        info[1] = 2
                if "bincoordinates" in line or "binCoordinates" in line:
                    if self.step == 'free':
                        info[0] = 'bincoordinates'
                        if self.parameter.iteration == 1:
                            info[1] = '../../restarts/' + self.parameter.outputname + '.' + \
                                      str(frame*self.parameter.sampling_interval) + '.coor'
                        else:
                            info[1] = self.parameter.outputname + '.coor'
                    elif new_sampling is not None:
                        info[0] = 'bincoordinates'
                        info[1] = './restarts/' +str(self.parameter.outputname) + '.' + str(new_sampling) + '.coor'
                
                if "binvelocities" in line or "binVelocities" in line:
                    if self.step == 'free':
                        info[0] = 'binvelocities'
                        if self.parameter.iteration == 1:
                            info[0] = '#binvelocities'
                            info[1] = '../../restarts/' + self.parameter.outputname + '.' + \
                                      str(frame*self.parameter.sampling_interval) + '.vel'
                        else:
                            if not self.parameter.NVT:
                                if enhanced == 0:
                                    info[0] = 'binvelocities'
                                else:
                                    info[0] = '#binvelocities'
                            if enhanced == 1:
                                info[0] = '#binvelocities'
                            info[1] = self.parameter.outputname + '.vel'
                    elif new_sampling is not None:
                        #info[0] = 'binvelocities'
                        info[1] = './restarts/' + str(self.parameter.outputname) + '.' + str(new_sampling) + '.vel'
            
                if "extendedSystem" in line or "extendedsystem" in line:
                    if self.step == 'free':
                        info[0] = 'extendedSystem'
                        if self.parameter.iteration == 1:
                            info[1] = '../../restarts/' + self.parameter.outputname + '.' + \
                                      str(frame * self.parameter.sampling_interval) + '.xsc'
                        else:
                            info[1] = self.parameter.outputname + '.xsc'
                    elif new_sampling is not None:
                        info[0] = 'extendedSystem'
                        info[1] = './restarts/' + str(self.parameter.outputname) + '.' + str(new_sampling) + '.xsc'
                    elif self.parameter.namd_conf == True and self.step != 'seek' and os.path.exists(self.parameter.seekPath) and self.parameter.iteration == 0:
                        info[0] = 'extendedSystem'
                        info[1] = './sample.xsc'
                        
                if "restartsave" in line:
                    if self.step != 'sample':
                        info[1] = "off"
    
                if "binaryrestart" in line:
                    if self.step == 'seek':
                        info[1] = "no"    
                if "set" in line and "valTotal" in line:
                    info[2] = 2*(len(self.parameter.anchors[0]))-1                         
                if "temperature" in line and "pressure" not in line:
                    if self.parameter.iteration > 1:
                        info[0] = '#temperature'
                    else:
                        info[0] = 'temperature'
                    if self.step == 'seek':
                        info[0] = 'temperature'
                    if enhanced == 1:
                        info[0] = 'temperature'

                # if "langevin" in line and self.parameter.nve and snapshot is not None::
                #     info[0] = '#'
                if "lreplace" in line:
#                    line = re.sub(r'[^\w]', ' ', line)
                    if self.parameter.colvarsNum == 0:
                        info[0] = '#set'
                    else:
                        if ']' in info[-1]:
                            info[-1] = str(self.parameter.colvarsNum - 1) 
                            info.append(']')
                        else:
                            info[-2] = str(self.parameter.colvarsNum - 1) 
                            
                if 'firsttimestep' in line and new_sampling is not None:
                    info[0] = 'firsttimestep'
                    info[1] = str(new_sampling)

                if "a111" in line:
                #As is this part is kind of messy, but leaving it like this just in case I am missing something
                #by commenting out the 'get_initial_ms' part for a111 and a222
                    if self.step != 'free':
                        info[2] = str(a1)
                    else:# snapshot != None:
                        info[2] = str(a1)
                        #if self.parameter.iteration == 1:
                        #    info[2] = str(a1) 
                        #else:
                        #    info[2] = str(get_initial_ms(MSpath)[0])
                        
                if "a222" in line:
                    if self.step != 'free':
                        info[2] = str(a2) 
                    else: # snapshot != None:
                        info[2] = str(a2)
                        #if self.parameter.iteration == 1:
                        #    info[2] = str(a2) 
                        #else:
                        #    info[2] = str(get_initial_ms(MSpath)[1]) 
                
                if self.step == 'seek' and "run" in info:
                    info[1] = str(self.parameter.initialTime * 1000)

                line = " ".join(str(x) for x in info)
                print(line)
    
    def move_files(self):
        '''
        When using namd, this moves all of our restart files to a folder called restart.
        When using Gromacs, there are no restart files so moving them is not neccessary,
        but the existance of the the restart folder is used to indicate that sampling has
        been done on that milestone.
        '''
        for name in self.loop1:
            [anchor1, anchor2] = get_anchors(name)
            filePath = self.parameter.crdPath + '/' + str(anchor1) + '_' + str(anchor2)
            restartFolder = filePath + '/restarts'
            create_folder(restartFolder)
            if os.path.exists(filePath + '/submit'):
                move(filePath + '/submit',filePath + '/submit_finished')
            for ext in ["coor", "vel", "xsc"]:
                for file in glob.glob(filePath + '/*.' + ext):
                    try:
                        move(file, restartFolder)
                    except:
                        continue
        
    def __submit_trajectories(self):
        '''
        This is where our jobs themselves are submitted. 
        '''
        print_log = False
        for i in self.loop1:
            print_individual = False
            for j in self.loop2:
                if self.step != 'sample':
                    trajectory_number = self.start_frame[i] + j - 1
                else:
                    trajectory_number = False
                if self.step == 'seek':
                    anchor1 = i
                    anchor2 = None
                else:
                    [anchor1, anchor2] = get_anchors(i)
                path = self.find_parent_path([anchor1, anchor2], trajectory_number)
                path += '/submit'
                if os.path.isfile(path):
                    print_log = True
                    print_individual = True
                    while True:
                        out = subprocess.check_output([self.parameter.jobcheck, '-u', self.parameter.username]).decode("utf-8").split('\n')
                        if len(out) - 2 < self.parameter.max_jobs:
                            subprocess.run([self.parameter.jobsubmit,path])
                            break
                        else:
                            time.sleep(60)
                    move(path, path + '_done')
            if self.step == 'free':
                if print_individual:
                   print(str(datetime.now()))
                   print("Short trajectories started from milestone {}...".format(i))
                else:
                   print("No new trajectories started from milestone {}".format(i))   
       
        return print_log

    def __check_progress(self):
        '''check running jobs on cluster'''
        finished = False
        while finished == False:
            finished = True
            out = subprocess.check_output([self.parameter.jobcheck,'-u', self.parameter.username]).decode("utf-8").split('\n')
            job = []
            if self.parameter.jobcheck == 'squeue' and len(out) > 2:
                job.append(list(filter(None, out[1].split(' '))))
            for i in range(2, len(out)-1):
                job.append(list(filter(None, out[i].split(' '))))
            for i in range(len(job)):
                if job[i][2].split('_')[0] == self.script_name:
                    finished = False
            if finished == False:
                print("Trajectories are running... Next check in {} seconds. {}".format(str(self.time), str(datetime.now())))
                time.sleep(self.time)
    
    def __post_launch(self, possible_ms = None):
        ''' After launching trajectories, we find our ending milestones here '''
        count = 0
        new_milestones = []
        for i in self.loop1:
            for j in self.loop2:
                count += 1
                [anchor1, anchor2] = get_anchors(i)
                trajectory_number = self.start_frame[i] + j - 1
                path = self.find_parent_path(i, trajectory_number)
                if not os.path.exists(path):
                    continue
                if self.step == 'seek':
                    anchor = i
                else:
                    anchor = None
                #if self.parameter.software == 'namd' or self.parameter.milestone_search in (2,3):
                time, final_ms = colvar(self.parameter).get_final_ms(path, anchor)
                #elif self.parameter.software == 'gromacs':
                #    time, final_ms = plumed(self.parameter, anchor1=anchor1, anchor2=anchor2, step=self.step).get_final_ms(path)
                final_ms.sort()
                if 0 in final_ms:
                    continue
                #print(final_ms)
                name = 'MS' + str(final_ms[0]) + '_' + str(final_ms[1])
                #if possible_ms is not None:
                #    if final_ms not in possible_ms:
                #        continue
                if self.step == 'seek':
                    ms = self.__post_seek(path,final_ms)
                    if ms:
                        self.parameter.MS_list.add(ms)
                else:
                    if self.parameter.new_ms_iterations < self.parameter.iteration and self.parameter.new_ms_iterations != 0:
                        self.parameter.ignorNewMS = True
                        continue
                    if self.step == 'free' and name not in self.parameter.MS_list:
                        new_milestones = self.__find_new_milestones(final_ms, new_milestones, path)
        
        #if self.step == 'free' and self.parameter.ignorNewMS == False and self.parameter.skip_compute == False:
        #    self.parameter.ignorNewMS = True
        return count
    
    def __find_new_milestones(self, final_ms, new_milestones, path):  
        new = True              
        final_ms.sort()
        if final_ms[0] == final_ms[1]:
            return new_milestones
        if new_milestones != []:
            for item in new_milestones:
                if item[0] == final_ms:
                    item[1] += 1
                    new = False
        if new == True:
            new_milestones.append([final_ms, 1, path])
        if self.parameter.new_ms_trajs != 0 and len(new_milestones) > 0:
            for ms in new_milestones:
                if ms[1] >= self.parameter.new_ms_trajs:
                    name = str(ms[0][0]) + '_' + str(ms[0][1])
                    if self.parameter.ignorNewMS == True:
                        log("A new milestone, {} was found and hit by {} trajectories but was not kept based on user parameters.".format(name, str(ms[1])))
                    else:
                        self.parameter.MS_list.add('MS' + name)
                        self.parameter.skip_compute = True
                        create_folder(self.find_parent_path(final_ms, no_iteration=True))
                        if self.parameter.method == 0 and not os.path.isfile(self.find_parent_path(final_ms, no_iteration = True) + '/statefromseek.cpt') and self.parameter.software == 'gromacs':
                            copy(path + '/state.cpt', self.find_parent_path(final_ms, no_iteration=True) + '/statefromseek.cpt')
        #print(new_milestones)
        return new_milestones
    
    def __post_seek(self, path, final_ms):
        sampling_path = self.parameter.crdPath + '/' + str(final_ms[0]) + '_' + str(final_ms[1])
        if self.parameter.software == 'namd':
            copy_path = path + '/' + self.parameter.outputname + '.restart.coor'
            dest_path = sampling_path + '/seek.ms.pdb'
        elif self.parameter.software == 'gromacs':
            copy_path = path + '/state.cpt'
            dest_path = sampling_path + '/statefromseek.cpt'
        if os.path.exists(sampling_path) or not os.path.exists(copy_path):
            return            
        if self.parameter.dist_cut != 0:
            if self.__check_distance(final_ms) == False:
                return
        create_folder(sampling_path)
        if self.parameter.software == 'namd' and self.parameter.namd_conf == True:
            copy(path + '/' + self.parameter.outputname + '.restart.coor', 
                 sampling_path + '/seek.ms.pdb')
            copy(path + '/' + self.parameter.outputname + '.xst', sampling_path + '/sample.xsc')
        copy(copy_path, dest_path)
        return 'MS' + str(final_ms[0]) + '_' + str(final_ms[1])
    
    def __check_distance(self, milestone):
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
                
    def find_parent_path(self, anchors, trajectory_number=False, no_iteration=False):
        import os
        if type(anchors) == list and self.step == 'seek':
            for i in range(len(anchors)):
                if anchors[i] == 0 or anchors[i] == None:
                    anchors.pop(i)
            anchors = anchors[0]
        if self.step == 'seek':
            path = os.path.join(self.parameter.seekPath, 'structure' + str(anchors))
        else:
            if type(anchors) == str:
                [anchor1, anchor2] = get_anchors(anchors)
            else:
                anchor1 = anchors[0]
                anchor2 = anchors[1]
            path = os.path.join(self.parameter.crdPath, str(anchor1) + '_' + str(anchor2))
            if self.step == 'free' and no_iteration == False:
                path = os.path.join(path, str(self.parameter.iteration))
        if trajectory_number != False:
            path = os.path.join(path, str(trajectory_number))
        return path
            
        
  
if __name__ == '__main__':
    from parameters import *
    from run import *
    from milestones import *
    
    new = parameters()
    new.initialize()
    new.method = 1
    new.milestone_search = 1
    #new.trajPerLaunch = 5
    #new.traj_per_script = [1,1,4]
    #new.initial_traj = 5
    new.MS_list = milestones(new).read_milestone_folder()
    #new.iteration = 1
    launch(new, step='sample').move_files()
    #free = launch(new, step='free')
    #free.post_launch()        
            
