# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 08:27:01 2020

@author: allis
"""
from additional_functions import *
import os

class parameters:
       

        
    def __init__(self, MS_list=None, Finished=None, MS_new=None, 
                 ignorNewMS=None, coor=None, NVT=None,
                 nodes=None, timeFactor=None, current_iteration_time=None,
                 sampling_interval=None, maxIteration=None, network=None,
                 boundary=None, iteration=None, method=None, pbs=None,
                 reactant_milestone=None, product_milestone=None,
                 reactant=None,product=None, bincoordinates=None, binvelocities=None, 
                 extendedSystem=None, customColvars=None, colvarsTrajFrequency=None, 
                 colvarsRestartFrequency=None, colvarsNum=None,
                 forceConst=None, trajWidths=None, username=None, tolerance=None, 
                 err_sampling=None, initial_traj=None, initialTime=None,
                 startTraj=None, trajPerLaunch=None, interval=None, freeTraj_walltime=None,
                 nframe=None, structure=None, coordinates=None, finished_constain=None,
                 outputname=None, namd_conf=None, AnchorPath=None, AnchorNum=None,  
                 new_anchor=None, anchor_dist=None, jobsubmit=None, jobcheck=None,
                 anchors=None, atomNumbers=None, error=None, MFPT=None, kij=None, 
                 index=None, flux=None, sing=None, seek_restartfreq=None, max_jobs=None,
                 colvarNumber=None, deltas=None, pause=None, traj_per_script=None, 
                 new_ms_iterations=None, new_ms_trajs=None, dist_cut=None, not_finish_trajs=None,
                 data_file=False, customMS_list = False) -> None:

        self.iteration = 0    # current iteration number
        self.method = 0       # 0 for classic milestoning; 1 for exact milestoning: iteration
        self.maxIteration = 100    # max iteration
        self.network = {}  # network for free trajectories
        self.milestone_search = 0   # 0 for traverse; 1 for seek
        self.pbc = []       # periodic boundary for 1D at the moment
        self.current_iteration_time = {}   # current iteration time
        self.timeFactor = 1    # fs per step
        self.sampling_interval = 1000   # output restarts every *** timesteps
        self.jobsubmit = jobsubmit    # command for job submission
        self.jobcheck = jobcheck    # command for job check 
        self.nodes = []  # available node list
        self.initial_traj = initial_traj    # number of trajs start from each anchor
        self.initialTime = initialTime     # time step in ps for initial trajs
        self.MS_list = set()    # milestone list
        self.ignorNewMS = False    # ignore new milestones found by free traj
        self.Finished = set()   # milestones that finished free trajs
        self.finished_constain = set()  # milestones that finished sampling
        self.MS_new = set()    # new milestone reached
        self.error = error     #
        self.NVT = False            # NVT for free trajectories
        self.boundary = [-1, -1]     # milestone number for reactant and product
        self.reactant_milestone = []   # milestone list for reactant
        self.product_milestone = []    # milestone list for product
        self.reactant = reactant     # voronoi cell for reactant
        self.product = product       # voronoi cell for product
        self.MFPT = MFPT       
        self.namd_conf = False       # additional modification required for NAMD configuration
        self.colvarsNum = colvarsNum  # number of colvars, used for read command
        self.forceConst = 1  # force constant for harmonic constrain
        self.trajWidths = trajWidths  # width for each traj output, default = 13
        self.customColvars = 0       # customized colvar: 0 - no; 1- yes
        self.colvarsTrajFrequency = colvarsTrajFrequency   # how often to output colvars
        self.colvarsRestartFrequency = colvarsRestartFrequency    # how often to write to file
        self.AnchorNum = AnchorNum     # total number of anchor
        self.new_anchor = False         # find new anchor
        self.anchor_dist = 100.0          # new anchor seperation distance
        self.bincoordinates = bincoordinates  # coordinates file name
        self.binvelocities = binvelocities # velocity file name
        self.extendedSystem = extendedSystem    # extended system file name
        self.coor = coor     # coordinates
        self.nframe = nframe    # number of sampling frame
        self.startTraj = startTraj    # initial sampling frame 
        self.trajPerLaunch = trajPerLaunch    # number of free trajs to launch each time
        self.freeTraj_walltime = freeTraj_walltime  # walltime for free trajectories
        self.interval = interval    # interval between each frame taken
        self.structure = structure   # structure file name; psf
        self.coordinates = coordinates   # structure file name; pdb
        self.outputname = outputname    # output name
        self.username = username    # user name on cluster
        self.tolerance = tolerance   # tolerance for MFPT convergency
        self.err_sampling = 1000    # number of error sampling
        self.sing = True      # k matrix singularity
        self.kij = []       #kij matrix
        self.index = []     # milestone index
        self.flux = []      #flux
        self.colvarNumber = colvarNumber
        self.seek_restartfreq = 1000
        self.anchors = anchors
        self.max_jobs = 9999999
        self.split_jobs = False
        self.ScMilesPath = os.path.dirname(os.path.abspath(__file__))
        self.parentDirectory = os.path.abspath(os.path.join(self.ScMilesPath, os.pardir))
        self.crdPath = os.path.join(self.parentDirectory, 'crd')
        self.seekPath = os.path.join(self.crdPath, 'seek')
        self.inputPath = os.path.join(self.parentDirectory, 'my_project_input')
        self.outputPath = os.path.join(self.parentDirectory, 'my_project_output')
        self.currentPath = os.path.join(self.outputPath, 'current')
        self.AnchorPath = os.path.join(self.inputPath, 'anchors.txt') # file path for anchor
        self.restart = False
        self.correctParameters = True
        self.traj_per_script = [1,1] #seek and free
        self.new_ms_iterations = 0
        self.new_ms_trajs = 0
        self.deltas = None
        self.dist_cut = 0
        self.not_finish_trajs = None
        self.additional_sampling = False
        self.data_file = False
        self.customMS_list = None

    def initialize(self):
        '''
        In this function, we read in all of the user input from input.txt (as well as a few things from free.namd
        '''
        import os
        import pandas as pd
        import re
        from log import log
        

        parameter_list= (('method', 'method', 'integer'),
                        ('max_iteration', 'maxIteration', 'integer'),
                        ('milestoneSearch','milestone_search', 'integer'),
                        ('pbc', 'pbc', 'replace_comma'),
                        ('structure', 'structure', 'string'),
                        ('coordinates', 'coordinates', 'string'),
                        ('outputname', 'outputname', 'string'),
                        ('NVT', 'NVT', 'yes_or_no'),
                        ('time_step','timeFactor','float'),
                        ('initial_traj','initial_traj','integer'),
                        ('initial_time', 'initialTime', 'integer'),
                        ('ignore_new_ms','ignorNewMS', 'yes_or_no'),
                        ('colvarNumber', 'colvarNumber', 'integer'),
                        ('colvarType', 'colvarType', 'string'),
                        ('custom_colvar', 'colvarsNum', 'integer'),
                        ('colvarsTrajFrequency', 'colvarsTrajFrequency', 'string'),
                        ('colvarsRestartFrequency', 'colvarsRestartFrequency', 'string'),
                        ('customColvars','customColvars', 'yes_or_no'),
                        ('force_const', 'forceConst', 'integer'),
                        ('anchorsNum', 'AnchorNum', 'integer'),
                        ('find_new_anchor','new_anchor', 'yes_or_no'),
                        ('new_anchor_dist', 'anchor_dist', 'float'),
                        ('reactant', 'reactant', 'replace_comma'),
                        ('product', 'product', 'replace_comma'),
                        ('total_trajs', 'nframe', 'integer'),
                        ('start_traj', 'startTraj', 'integer'),
                        ('traj_per_launch', 'trajPerLaunch', 'integer'),
                        ('interval', 'interval', 'integer'),
                        ('tolerance', 'tolerance', 'float'),
                        ('error_sampling','err_sampling','integer'),
                        ('jobsubmission','jobsubmit', 'string'),
                        ('jobcheck','jobcheck','string'),
                        ('username','username','string'),
                        ('namd_conf_custom', 'namd_conf', 'yes_or_no'),
                        ('restart', 'restart', 'yes_or_no'),
                        ('seek_restartfreq', 'seek_restartfreq', 'integer'),
                        ('max_jobs', 'max_jobs', 'integer'),
                        ('split_jobs', 'split_jobs', 'integer'),
                        ('traj_per_script', 'traj_per_script', 'replace_comma'),
                        ('new_ms_trajs','new_ms_trajs', 'integer'),
                        ('new_ms_iterations','new_ms_iterations','integer'),
                        ('deltas','deltas','replace_comma'),
                        ('MS_list','MS_list','replace_comma'),
                        ('additional_sampling','additional_sampling', 'yes_or_no'),
                        ('data_file','data_file','yes_or_no'),
                        ('dist_cut','dist_cut','float'))
                              
        with open(file = self.inputPath +'/input.txt') as r:
            for line in r:
                line = line.rstrip()
                info = line.split(" ")
                if info == [] or line.startswith('#'):
                    continue
                for item in parameter_list:
                    if item[0] in info:
                        if item[2] == 'integer':
                            setattr(self, item[1], int(info[1]))
                        elif item[2] == 'float':
                            setattr(self, item[1], float(info[1]))
                        elif item[2] == 'yes_or_no':
                            if str(info[1]).lower() in ('true','yes','on','1'):
                                setattr(self, item[1],True)
                            else:
                                setattr(self, item[1],False)
                        elif item[2] == 'replace_comma':
                            rm = line.replace(","," ").replace("  "," ").split(" ")
                            rm.pop(0)
                            if item[0] == 'deltas':
                                for i in range(len(rm)):
                                    rm[i] == float(rm[i])
                                setattr(self, item[1], rm)
                            try:
                                setattr(self, item[1], list(map(int, rm)))
                            except:
                                MS_list = set()
                                for i in rm:
                                    MS_list.add(i)
                                self.ignorNewMS = True
                                self.customMS_list = MS_list
                                setattr(self, item[1], MS_list) 
                        elif item[2] == 'string':
                            setattr(self, item[1], str(info[1]))
            
        self.trajWidths = [13]
        for i in range(self.colvarsNum + self.AnchorNum):
            self.trajWidths.append(23)
        
        if os.path.isfile(os.path.join(self.ScMilesPath, 'nodelist')):
            with open(file=nodelist) as f:
                for line in f:
                    if "#" in line:
                        continue
                    line = line.split("\n")
                    self.nodes.append(str(line[0]))
                            
        self.anchors = pd.read_table(self.AnchorPath, delimiter='\s+', header=None).values
        '''
        if self.colvarNumber != len(self.anchors[0]):
            self.correctParameters = False
            print('Number of columns in anchors.txt ({}) does not match the colvarNumber specified ({}) in input.txt. Please make sure these numbers match'
                  .format(len(self.anchors[0]), self.colvarNumber))
        '''
        create_folder(self.crdPath)
        create_folder(self.outputPath)
        create_folder(self.currentPath)
        # read initial run time for seek and time step setup
        with open(self.inputPath + '/free.namd', 'r') as f:   
            for line in f:
                info = line.split("#")[0].split()
            # info = line.split()
                if len(info) < 1 or line.startswith('#'):
                    continue
                if 'run' in info[0].lower():
                    self.freeTraj_walltime = int(re.findall(r"[-+]?\d*\.\d+|\d+", info[1])[0])
                elif 'timestep' in info[0].lower():
                    try: 
                        self.timeFactor = float(re.findall(r"[-+]?\d*\.\d+|\d+", info[1])[0])
                    except:
                        continue             
        # read restart frequency to get the name of restart files
        # such restart files will be used as the initial position for free traj
        with open(os.path.join(self.inputPath,'sample.namd'), 'r') as f:
            for line in f:
                info = line.split("#")[0].split()
                if len(info) < 1 or line.startswith('#'):
                    continue
                if "restartfreq" in info[0].lower():
                    self.sampling_interval = int(re.findall(r"[-+]?\d*\.\d+|\d+", info[1])[0])
        # initial log file
        if self.restart == False:
            if os.path.exists(os.path.join(self.currentPath, 'log')):
                os.remove(os.path.join(self.currentPath, 'log'))           
            log("Initialized with {} anchors.".format(self.AnchorNum))
            
    def ms_to_path(ms):
        import re
        [anchor1, anchor2] = list(map(int, (re.findall('\d+', ms))))
        return self.crdPath + '/' + str(anchor1) + '_' + str(anchor2)
                        
if __name__ == '__main__':
    new = parameters()
    new.initialize()
    print(new.customColvars)
    print(new.anchors)

    
