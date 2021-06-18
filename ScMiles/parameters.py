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
                 data_file=False, customMS_list = False, software=None, MS_discarded=None,
                 CV_suffixes=None, gromacs_timestep=None, ndx_file=None, skip_compute=None, active_anchors=None,
                 plots=None, grid_pbc=None, skip_MS=None, l_values=None, all_grid=None, grid_ignore=None,
                 grid_dict=None, grid_caps=None, milestone_delta=None, corners=None, names=None, 
                 distance_lambda=None, pbc_names=None, l=None, substitution=None, k_cutoff=None,
                 max_lifetime=None, compute_only=None, skip_sampling=None, targetNumSteps=None, pdb_sampling=False, restartfreq=None,
                 ignore_transitions = None, analysis_ignore_milestones=None, k_min_sum=None, max_grid_value=None,
                 min_grid_value=None) -> None:

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
        self.initial_traj = 10    # number of trajs start from each anchor
        self.initialTime = 50     # time step in ps for initial trajs
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
        #These filepaths are just to simplify so we don't have to find them over and over again
        self.ScMilesPath = os.path.dirname(os.path.abspath(__file__))
        self.parentDirectory = os.path.abspath(os.path.join(self.ScMilesPath, os.pardir))
        self.crdPath = os.path.join(self.parentDirectory, 'crd')
        self.seekPath = os.path.join(self.crdPath, 'seek')
        self.inputPath = os.path.join(self.parentDirectory, 'my_project_input')
        self.outputPath = os.path.join(self.parentDirectory, 'my_project_output')
        self.currentPath = os.path.join(self.outputPath, 'current')
        self.AnchorPath = os.path.join(self.inputPath, 'anchors.txt') # file path for anchor
        self.restart = False
        self.correctParameters = True #Returns False and terminates simulation if parameters are not correct.
        self.traj_per_script = [1,1,1] #seek and free
        self.new_ms_iterations = 0 
        self.new_ms_trajs = 1 #Number of hits needed to keep a new milestone
        self.deltas = None #Used in grid, difference between two anchors. Entered as a list
        self.dist_cut = 0 #Used with voronoi, only keeps mielstones containing anchors of a certain distance
        self.not_finish_trajs = None
        self.additional_sampling = False #If user wants to do additional sampling after the original sampling
        self.customMS_list = None #User can enter a custom milestone list to use
        self.software = 'namd' #Namd, gromacs or lammps (lammps is untested as of now)
        self.neighbor_kappa = 0 #Used with gromacs in plumed file
        self.walls_kappa = 0 #Used with gromacs in plumed file
        self.MS_discarded = []
        self.CV_suffixes = [] 
        self.gromacs_timestep = 0.1
        self.ndx_file = None
        self.skip_compute = False #If we find new milestones we don't compute because it breaks the code
        self.active_anchors = None 
        self.plots = False
        self.grid_pbc = False #I don't think this is working as of now. It was kind of a temporary solution for something
        self.skip_MS = []
        self.l_values = [0] #Used with periodicity
        self.all_grid = dict() #The script fills this with the grid definitions (done in milestones.py)
        self.grid_ignore = set() #The user specifies if certain milestones should be ignored (will not be added to all_grid)
        self.grid_dict = None 
        self.grid_caps = False #Used if there are ends that do not have milestones.
        self.milestone_delta = None #Creates 'thick' milestones that are more of a zone than just a line
        self.corners = None #Used for grid, allows corner milestones to be used
        self.names = [] #This is our variable names, inputted by the user
        self.distance_lambda = None 
        self.pbc_names = [] #
        self.l = None
        self.substitution = None
        self.origin = 1
        self.k_cutoff = 0
        self.max_lifetime = None
        self.velocity_weighting = False
        self.scale_rmsd = False
        self.compute_only = False
        self.targetNumSteps = 50000
        self.pdb_sampling = False
        self.restartfreq = None
        self.ignore_transitions = None
        self.analysis_ignore_milestones = None    
        self.k_min_sum = None
        self.min_grid_value = None
        self.max_grid_value = None

    def initialize(self):
        '''
        In this function, we read in all of the user input from input.txt (as well as a few things from free.namd)
        I changed this to this matrix style just to make adding things easier.
        The first column is what it looks for in input, second column is the parameter name, third is how to save the variable
        '''
        import os
        import pandas as pd
        import re
        from log import log 
        create_folder(self.outputPath)
        create_folder(self.currentPath)

        incorrect_input = []
        parameter_list = (('method', 'method', 'integer'),
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
                        ('force_const', 'forceConst', 'replace_comma'),
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
                        ('dist_cut','dist_cut','float'),
                        ('software','software','string'),
                        ('walls_kappa','walls_kappa','float'),
                        ('neighbor_kappa','neighbor_kappa','float'),
                        ('CV_suffixes','CV_suffixes','replace_comma'),
                        ('gromacs_timestep','gromacs_timestep','float'),
                        ('ndx_file','ndx_file','string'),
                        ('active_anchors', 'active_anchors', 'replace_comma'),
                        ('plots','plots','yes_or_no'),
                        ('grid_pbc','grid_pbc','string'),
                        ('l_values','l_values','replace_comma'),
                        ('grid_caps','grid_caps','yes_or_no'),
                        ('grid_ignore','grid_ignore','replace_comma'),
                        ('milestone_delta','milestone_delta','float'),
                        ('corners','corners','yes_or_no'),
                        ('distance_lambda','distance_lambda','float'),
                        ('pbc_names', 'pbc_names','replace_comma'),
                        ('L','l','replace_comma'),
                        ('substitution', 'substitution', 'yes_or_no'),
                        ('origin','origin','integer'),
                        ('k_cutoff','k_cutoff','integer'),
                        ('velocity_weighting','velocity_weighting','yes_or_no'),
                        ('max_liftime','max_lifetime','float'),
                        ('scale_rmsd','scale_rmsd','replace_comma'),
                        ('colvar_names','names','replace_comma'),
                        ('compute_only','compute_only','yes_or_no'),
                        ('targetNumSteps','targetNumSteps','integer'),
                        ('pbc_sampling','pbc_sampling','yes_or_no'),
                        ('k_min_sum','k_min_sum','integer'),
                        ('max_grid_value','max_grid_value','replace_comma'),
                        ('min_grid_value','min_grid_value','replace_comma'))
                              
        with open(file = self.inputPath +'/input.txt') as r:
            for line in r:
                found = False
                line = line.rstrip()
                info = line.split(" ")
                if info == [] or line.startswith('#') or info == ['']:
                    continue
                print(info)
                for item in parameter_list:
                    if item[0] == info[0]:
                        print(item[0])
                        found = True
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
                            if item[0] == 'deltas' or item[0] == 'L' or item[0] == 'scale_rmsd':
                                for i in range(len(rm)):
                                    rm[i] = float(rm[i])
                                setattr(self, item[1], rm)
                                continue
                            elif item[0] == 'CV_suffixes' or item[0] == 'grid_ignore' or item[0] == 'pbc_names' or item[0] == 'colvar_names' or '_grid_value' in item[0]:
                                setattr(self, item[1], rm)
                                continue
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
                if found == False:
                    incorrect_input.append(info[0])
        
        if self.restart == False:
            if os.path.exists(os.path.join(self.currentPath, 'log')):
                os.remove(os.path.join(self.currentPath, 'log'))         
        if incorrect_input:
            log('The keywords ' + ', '.join(incorrect_input) + ' do not exist and will be ignored')
        #if self.active_anchors != None:
        #    self.active_anchors = range(self.active_anchors[0], self.active_anchors[1] + 1)     
        self.trajWidths = [13]
        for i in range(self.colvarsNum + self.AnchorNum):
            self.trajWidths.append(23)
        
        if self.min_grid_value is not None:
            for i in range(len(self.min_grid_value)):
                try:
                    self.min_grid_value[i] = float(self.min_grid_value[i])
                except:
                    self.min_grid_values[i] = None
        if self.max_grid_value is not None:
            for i in range(len(self.max_grid_value)):
                try:
                    self.max_grid_value[i] = float(self.max_grid_value[i])
                except:
                    self.max_grid_values[i] = None
                    
        if os.path.isfile(os.path.join(self.ScMilesPath, 'nodelist')):
            with open(file=nodelist) as f:
                for line in f:
                    if "#" in line:
                        continue
                    line = line.split("\n")
                    self.nodes.append(str(line[0]))
                            
        self.anchors = pd.read_table(self.AnchorPath, delimiter='\s+', header=None).values
        if self.software == 'gromacs':
            for i in range(len(self.anchors)):
                for j in range(len(self.anchors[0])):
                    self.anchors[i][j] = round(self.anchors[i][j],5)
        #print(self.anchors)
        
        self.check_parameters()
        
        create_folder(self.crdPath)
        # read initial run time for seek and time step setup
        if self.milestone_search == 0:
            self.pbc_sampling = True          
        if len(self.forceConst) == 1 and len(self.anchors[0]) != 1:
            for i in range(len(self.anchors[0])):
                self.forceConst.append(self.forceConst[0])

        if self.software == 'namd':
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
                    elif 'restartfreq' in info[0].lower():
                        self.restartfreq = int(info[1])             

        # read restart frequency to get the name of restart files
        # such restart files will be used as the initial position for free traj
        if self.software == 'namd':
            ext = '.namd'
        else:
            ext = '.mdp'
        with open(os.path.join(self.inputPath,'sample' + ext), 'r') as f:
                for line in f:
                    info = line.split("#")[0].split()
                    if len(info) < 1 or line.startswith('#'):
                        continue
                    if "restartfreq" in info[0].lower() and self.software == 'namd':
                        self.sampling_interval = int(re.findall(r"[-+]?\d*\.\d+|\d+", info[1])[0])
                    elif "nstxout" in info[0].lower() and self.software == 'gromacs':
                        self.sampling_interval = int(re.findall(r"[-+]?\d*\.\d+|\d+", info[2])[0])
        # initial log file         
        log("Initialized with {} anchors.".format(self.AnchorNum))
        
    def check_parameters(self):
        #If any of these return False, the script will terminate and tell the user to fix their inputs.
        if len(self.names) != len(self.anchors[0]):
            self.correctParameters = False
            print('Number of columns in anchors.txt ({}) does not match the number of colvars specified ({}) in input.txt (colvar_names). Please make sure these numbers match'
                  .format(len(self.anchors[0]), self.colvarNumber))
        if self.software != 'namd' and self.velocity_weighting == True:
            self.correctParameters = False
            print('Velocity Weighting is currently only available with namd. Please either turn off velocity weighting or change your software to namd')
        if self.pbc_names:
            if len(self.pbc_names) != len(self.l):
                self.correctParameters = False
                print('There must be an L value for each coarse variable in pbc_names.')
        if self.milestone_search in (0,1) and self.deltas != None:
            log('Using milestone_search = 0 or milestone_search = 1 means that you are using Voronoi cell, so delta values are not used and will be ignored')
            print('Using milestone_search = 0 or milestone_search = 1 means that you are using Voronoi cell, so delta values are not used and will be ignored')
        if self.software not in ('namd','gromacs','lammps'):
            self.correctParameters = False
            print('The options for software are namd, lammps, and gromacs. Please update your input.txt to use one of these options')
        if self.scale_rmsd != False and len(self.scale_rmsd) != len(self.names):
            self.correctParameters = False
            print('If you want to use the input scale_rmsd, please make sure you use the same value as the number of colvars. If you do not want one to be scaled, please just use the value 1.0 for that colvar')
        #if self.restart == True and not os.path.isfile(self.crdPath):
        #    self.restart = False
        #    print('ScMiles is not far enough for the restart option to be used. It will start from the beginning')
        if self.AnchorNum != len(self.anchors):
            self.correctParameters = False
            print('The number of rows in anchor.txt does not match the anchorsNum in input.txt. Please update your files so you have one row for each anchor')

if __name__ == '__main__':
    new = parameters()
    new.initialize()

    
