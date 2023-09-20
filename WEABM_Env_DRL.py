# WEABM DRL Environment-V1
# Created by Dale Larie 7/8/2021

# This version was updated/cleaned on 6/13/2023

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import wrapper_setup as wrapper_setup

import gym
import mpy4py as MPI

"""Setting up all the global/settings variables"""

# labels of the cytokines in the simulation, in order
cytokine_labels = ["IL8", "MIP1b", "TNF", "HGF", "PAF", "IL1", "PDGF", "IL13", "TGFb", "VEGF", "GCSF", "IL10", "IL4", "IL17", "Myf5", "IL6", "IL1ra", "sIL1r", "sTNFr", "endotoxin", "IFNg", "IL12", "MRF4", "cytotox", "IGF", "DAMP", "MCP1", "FNE"]

MAX_STEPS = 1400                        #Max number of steps the simulation will play
TOTAL_OBSERVABLE_CYTOKINES = 28         #Number of cytokines in the WEABM
TOTAL_CONTROLLABLE_CYTOKINES = 28       #Number of cytokines that can be controlled in the WEABM

# OBSERVATION_CYTOKINES = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]       #Looking at all cytokines
OBSERVATION_CYTOKINES = [0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0]       # Looking at M1 Genes (TNF, IL1, IL6, IL1ra, IFNg, IL12)
# print("OBSERVED CYTOKINES " + str([observed_cytokine for observed_cytokine, index in zip(cytokine_labels, OBSERVATION_CYTOKINES) if index == 1]))


SUM_M1_GENES = False
if SUM_M1_GENES:
    NUM_OBSERVABLE_CYTOKINES=1      #summing M1 genes to represent NO
else:
    NUM_OBSERVABLE_CYTOKINES = sum(OBSERVATION_CYTOKINES) #if all cytokines are being observed individually

OBSERVE_TGFb = True    #index =8, observing seperately, dont include in the observation vector because that will add it to the NO sum

# ACTION_INDECES = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]  # using all available actions
ACTION_INDECES = [0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]  # using {TGFb, TNF, IL10, HGF, cytotox}
# ACTION_INDECES = [0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  # using only TNF, IL10 and IL4
# ACTION_INDECES = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  # totally uncontrolled
# print([observed_cytokine for observed_cytokine, index in zip(cytokine_labels, OBSERVATION_CYTOKINES) if index == 1])
# print("ACTION CYTOKINES " + str([observed_cytokine for observed_cytokine, index in zip(cytokine_labels, ACTION_INDECES) if index == 1]))

NUM_CYTOKINES_CONTROLLED = sum(ACTION_INDECES)   #Number of cytokines controlled by the agent

# empirically observed maximum values for each cytokine with random control. These values are used for normalization/ input scaling
observed_signal_max = [6507.83349609,   51.18336487, 1526.96228027 ,  56.05470657 ,1716.56884766,
   65.22439575,  105.38170624, 2129.74145508, 2632.51220703,  100.20832825,
 1759.00390625, 1949.2902832 , 2390.04394531,  640.81066895 , 113.35812378,
 1698.16662598,   53.98086929,   52.61504745,   97.77040863  , 97.39041138,
 2115.72216797, 2548.54370117,   95.35305786,  283.16220093,   95.18818665,
 2248.67578125, 2238.12207031,  939.54016113]
signal_max = observed_signal_max
TGFb_max_obs = signal_max[8]    # tracking TGFb for observing TGFb seperately

signal_max = [signal_max[i] for i in range(len(OBSERVATION_CYTOKINES)) if OBSERVATION_CYTOKINES[i]==1] #getting just observed cytos max values

# Min values for all cytokines
signal_min = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
TGFb_min_obs = signal_min[8]   # tracking TGFb for observing TGFb seperately

signal_min = [signal_min[i] for i in range(len(OBSERVATION_CYTOKINES)) if OBSERVATION_CYTOKINES[i]==1] #getting just observed cytos min values

if SUM_M1_GENES:
    signal_max = sum(signal_max)    #when summing M1 genes
    signal_min = sum(signal_min)    #when summing M1 genes

# initializing the python ctypes wrapper. Simulation dimensions are specified in wrapper_setup.py to get correct return array shapes
SIM = wrapper_setup.setUpWrapper()


def createWEABM(seed):
    """ function that will create an instance of the WEABM simulation and set it up to run.
        seed (int): Takes an int that specifies the RNG seed in the C++ code
    """

    # getting rule matrix
    IP=np.genfromtxt('RuleMat3.csv',delimiter=',')
    IP=IP.flatten()
    numMatrixElements=IP.shape[0]
    array_type = ctypes.c_float*numMatrixElements

    # creating and setting up the WEABM instance
    instance = SIM.CreateInstance(array_type(*IP))
    SIM.setSeed(instance, seed)
    SIM.InitializeSimulation(instance)
    SIM.DoInjury(instance)
    return instance

class WEABM_Environment(gym.Env):
    """ WEABM_Environment class
        The actual gym object that the DRL agents will interact with, extends the OpenAI gym class
        The agent will observe cytokine levels present on the surface of the simulated wound, then take an action based on that obervation. The action will either add an amount of a particular cytokine evenly to every cell on the surface, or take away up to all of a particular cytokine present in the surface cells.

    """
    def __init__(self, action_scaling_factor=1, SAVE_LOCATION=None, rank=None, seed=0):
        """ action_scaling_factor (float): The amount that the agents actions will be scaled when doing the action.
                Agent action choices are between -1 and 1, so this factor scales the chosen number
            SAVE_LOCATION (string): path like object to determine where the simulation data will be saved
        """
        super(WEABM_Environment, self).__init__()
        self.save_path = SAVE_LOCATION
        self.rank = rank                                            # rank for saving data if running envs in parallel
        self.MAX_STEPS = MAX_STEPS                                  # Max steps for simulation episode. Set as a global
        self.seed = seed                                               # Seed for the simulation
        self.current_step = 0                                       # The current step of the simulation environment
        self.current_cytos = np.zeros(TOTAL_OBSERVABLE_CYTOKINES)   # Placeholder array to store current oberved surface cytokines
        self.current_viable_muscle = 0                              # counters to track current totals of muslce, collagen, etc
        self.current_collagen1 = 0
        self.current_collagen3 = 0
        self.current_necrosis = 0
        self.current_life = 0
        self.previous_step_life = 0
        self.action_scaling = action_scaling_factor

        self.ptrToEnv = createWEABM(self.seed)        #creating a dummy instance of the simulation to get x,y dims
        self.xDim = SIM.get_xDim(self.ptrToEnv)       #<--- set in reset
        self.yDim = SIM.get_yDim(self.ptrToEnv)       #<---

        """Setting up arrays of histories for the simulation episode. These will be filled out with value(s) on each step"""
        self.observation_history = np.zeros(shape=(MAX_STEPS,2+NUM_OBSERVABLE_CYTOKINES+ int(OBSERVE_TGFb))) #normal env; (observing current muscle, current collagen, current agent observed cytokines, TGFb seperately)
        self.reward_history = np.zeros(MAX_STEPS)
        self.viable_muscle_history = np.zeros(MAX_STEPS)
        self.collagen1_history = np.zeros(MAX_STEPS)
        self.collagen3_history = np.zeros(MAX_STEPS)
        self.necrosis_history = np.zeros(MAX_STEPS)
        self.life_history = np.zeros(MAX_STEPS)
        self.cytokine_history = np.zeros(shape=(MAX_STEPS, TOTAL_OBSERVABLE_CYTOKINES)) # saving all cytokines in sim
        self.norm_hist = np.zeros(shape=(MAX_STEPS, TOTAL_OBSERVABLE_CYTOKINES))        # saving normalized cytokines
        self.action_history = np.zeros(shape=(MAX_STEPS, TOTAL_CONTROLLABLE_CYTOKINES)) # saving actions taken by agent

        self.allSignalsReturn_output = np.zeros(43*MAX_STEPS)   #placeholder for output of "end_simulation()", returns the standard C++ output of the simulation
        self.current_action = np.zeros(TOTAL_CONTROLLABLE_CYTOKINES) #placeholder for agent chosen action

        # defining agent action space
        self.action_space = gym.spaces.Box(
            low=-1,
            high=1,
            shape=(NUM_CYTOKINES_CONTROLLED,),
            dtype=np.float32)

        # observation space is TotalViableMuscle, TotalCollagen3 and 28 cytokine values plus tgfb if observing seperately
        obs_space_high = np.zeros(2+NUM_OBSERVABLE_CYTOKINES + int(OBSERVE_TGFb))
        obs_space_high[0] = self.xDim*self.yDim*20
        obs_space_high[1] = self.xDim*self.yDim*2
        if OBSERVE_TGFb:
            obs_space_high[-1] = TGFb_max_obs
            obs_space_high[2:-1] = np.array(signal_max)
        else:
            obs_space_high[2:] = np.array(signal_max)
        obs_max = obs_space_high

        obs_space_low = np.zeros(2+NUM_OBSERVABLE_CYTOKINES + int(OBSERVE_TGFb))
        obs_space_low[0] = 0
        obs_space_low[1] = 0
        if OBSERVE_TGFb:
            obs_space_low[-1] = TGFb_min_obs
            obs_space_low[2:-1] = np.array(signal_min)
        else:
            obs_space_low[2:] = np.array(signal_min)

        obs_min = obs_space_low

        self.observation_space = gym.spaces.Box(
            low=obs_space_low,
            high=obs_space_high,
            shape=obs_space_high.shape,
            dtype=np.float32)

        # call reset to initialze all the variables so it is ready to simulatie
        self.reset()

    def setSeed(self, new_seed):
        """sets the seed in the C++ simulation
           new_seed (int): new seed for RNG
        """
        SIM.setSeed(self.ptrToEnv, new_seed)

    def step(self, action):
        """steps the simulation forward one step, based on the action chosen by the agent
           action (array[floats]): an array of floats between -1 and 1 that is NUM_CYTOKINES_CONTROLLED in length

           returns:
           observation (array[floats]): the observation returned to the agent
           reward (float): The reward for this step based on the state of the model
           done (bool): Whether the episode has finished or not
           info (dict): extra info about the environment at the current step

        """
        self.take_action(action)
        observation = self.next_observation()
        done = self.calculate_done()
        self.done = done
        reward = self.calculate_reward()
        self.reward = reward

        done = bool(done) #done starts as an int for what type of finish. convert to a bool now for return
        if done:
            if self.save_path is not None:
                self.save_training_data()

        # If we havent reached the end of the time limit, add current step and reward to histories
        if self.current_step < MAX_STEPS:
            self.observation_history[self.current_step,:] = observation
            self.reward_history[self.current_step] = reward

        # return basically all the info we can about the simulation, including full histories this can probably be cut down if the user isnt interested in them on current step.
        info = {"muscle":self.current_viable_muscle,
                "step":self.current_step,
                "necrosis":np.sum(np.array(SIM.getNecrosis(self.ptrToEnv))),
                "collagen1":np.sum(np.array(SIM.getCollagen1(self.ptrToEnv))),
                "collagen3":np.sum(np.array(SIM.getCollagen3(self.ptrToEnv))),
                "life":np.sum(np.array(SIM.getLife(self.ptrToEnv))),
                "muscle_hist":self.viable_muscle_history,
                "collagen1_hist":self.collagen1_history,
                "collagen3_hist":self.collagen3_history,
                "necrosis_hist":self.necrosis_history,
                "life_hist":self.life_history,
                "observation_hist": self.observation_history,
                "action_hist": self.action_history,
                "cytokine_hist":self.cytokine_history,
                "norm_hist":self.norm_hist,
                "allSignalsReturn_output": self.allSignalsReturn_output,
                "muscle_spatial": self.muscle_spatial,
                "collagen1_spatial": self.collagen1_spatial,
                "collagen3_spatial": self.collagen3_spatial,
                "necrosis_spatial": self.necrosis_spatial
                }
        return observation, reward, done, info

    def save_training_data(self):
        """ saves the relevant data from the simulation from training
        """
        if rank is None:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
        else:
            rank = self.rank
        with open(self.save_path+"/" + str(rank) + "_training_data.csv", "ab") as f:
            saveArray = self.viable_muscle_history
            saveArray = np.vstack((saveArray, self.collagen3_history))
            saveArray = np.vstack((saveArray, self.necrosis_history))
            saveArray = np.vstack((saveArray, self.collagen1_history))
            saveArray = np.vstack((saveArray, self.norm_hist.T))
            saveArray = np.vstack((saveArray, self.action_history.T))
            np.savetxt(f, saveArray, delimiter=",")

    def take_action(self, action):
        """ implements the action that the DRL agent chose to take
            action (array[floats]) of size NUM_CYTOKINES_CONTROLLED
        """
        action_num = 0
        action *= self.action_scaling   # scaling chosen actions

        # implementing the actions. going through all potential actions, and only acting on the ones indicated in the ACTION_INDECES vector
        # this is adding the amount of cytokine to the surface of the wound
        if ACTION_INDECES[0] != 0:
            SIM.apply_IL8(self.ptrToEnv, action[action_num])
            self.current_action[0] = action[action_num]
            action_num += 1

        if ACTION_INDECES[1] != 0:
            SIM.apply_MIP1b(self.ptrToEnv, action[action_num])
            self.current_action[1] = action[action_num]
            action_num += 1

        if ACTION_INDECES[2] != 0:
            SIM.apply_TNF(self.ptrToEnv, action[action_num])
            self.current_action[2] = action[action_num]
            action_num += 1

        if ACTION_INDECES[3] != 0:
            SIM.apply_HGF(self.ptrToEnv, action[action_num])
            self.current_action[3] = action[action_num]
            action_num += 1

        if ACTION_INDECES[4] != 0:
            SIM.apply_PAF(self.ptrToEnv, action[action_num])
            self.current_action[4] = action[action_num]
            action_num += 1

        if ACTION_INDECES[5] != 0:
            SIM.apply_IL1(self.ptrToEnv, action[action_num])
            self.current_action[5] = action[action_num]
            action_num += 1

        if ACTION_INDECES[6] != 0:
            SIM.apply_PDGF(self.ptrToEnv, action[action_num])
            self.current_action[6] = action[action_num]
            action_num += 1

        if ACTION_INDECES[7] != 0:
            SIM.apply_IL13(self.ptrToEnv, action[action_num])
            self.current_action[7] = action[action_num]
            action_num += 1

        if ACTION_INDECES[8] != 0:
            SIM.apply_TGFb(self.ptrToEnv, action[action_num])
            self.current_action[8] = action[action_num]
            action_num += 1

        if ACTION_INDECES[9] != 0:
            SIM.apply_VEGF(self.ptrToEnv, action[action_num])
            self.current_action[9] = action[action_num]
            action_num += 1

        if ACTION_INDECES[10] != 0:
            SIM.apply_GCSF(self.ptrToEnv, action[action_num])
            self.current_action[10] = action[action_num]
            action_num += 1

        if ACTION_INDECES[11] != 0:
            SIM.apply_IL10(self.ptrToEnv, action[action_num])
            self.current_action[11] = action[action_num]
            action_num += 1

        if ACTION_INDECES[12] != 0:
            SIM.apply_IL4(self.ptrToEnv, action[action_num])
            self.current_action[12] = action[action_num]
            action_num += 1

        if ACTION_INDECES[13] != 0:
            SIM.apply_IL17(self.ptrToEnv, action[action_num])
            self.current_action[13] = action[action_num]
            action_num += 1

        if ACTION_INDECES[14] !=0:
            SIM.apply_Myf5(self.ptrToEnv, action[action_num])
            self.current_action[14] = action[action_num]
            action_num += 1

        if ACTION_INDECES[15] !=0:
            SIM.apply_IL6(self.ptrToEnv, action[action_num])
            self.current_action[15] = action[action_num]
            action_num += 1

        if ACTION_INDECES[16] !=0:
            SIM.apply_IL1ra(self.ptrToEnv, action[action_num])
            self.current_action[16] = action[action_num]
            action_num += 1

        if ACTION_INDECES[17] !=0:
            SIM.apply_sIL1r(self.ptrToEnv, action[action_num])
            self.current_action[17] = action[action_num]
            action_num += 1

        if ACTION_INDECES[18] !=0:
            SIM.apply_sTNFr(self.ptrToEnv, action[action_num])
            self.current_action[18] = action[action_num]
            action_num += 1

        if ACTION_INDECES[19] !=0:
            SIM.apply_endotoxin(self.ptrToEnv, action[action_num])
            self.current_action[19] = action[action_num]
            action_num += 1

        if ACTION_INDECES[20] !=0:
            SIM.apply_IFNg(self.ptrToEnv, action[action_num])
            self.current_action[20] = action[action_num]
            action_num += 1

        if ACTION_INDECES[21] !=0:
            SIM.apply_IL12(self.ptrToEnv, action[action_num])
            self.current_action[21] = action[action_num]
            action_num += 1

        if ACTION_INDECES[22] !=0:
            SIM.apply_MRF4(self.ptrToEnv, action[action_num])
            self.current_action[22] = action[action_num]
            action_num += 1

        if ACTION_INDECES[23] !=0:
            SIM.apply_cytotox(self.ptrToEnv, action[action_num])
            self.current_action[23] = action[action_num]
            action_num += 1

        if ACTION_INDECES[24] !=0:
            SIM.apply_IGF(self.ptrToEnv, action[action_num])
            self.current_action[24] = action[action_num]
            action_num += 1

        if ACTION_INDECES[25] !=0:
            SIM.apply_DAMP(self.ptrToEnv, action[action_num])
            self.current_action[25] = action[action_num]
            action_num += 1

        if ACTION_INDECES[26] !=0:
            SIM.apply_MCP1(self.ptrToEnv, action[action_num])
            self.current_action[26] = action[action_num]
            action_num += 1

        if ACTION_INDECES[27] !=0:
            SIM.apply_FNE(self.ptrToEnv, action[action_num])
            self.current_action[27] = action[action_num]
            action_num += 1

        self.action_history[self.current_step,:] = self.current_action  #adding actions to action history

        SIM.StepSimulation(self.ptrToEnv)   # actually stepping simulation forward
        self.current_step+=1

    def next_observation(self):
        """ Function that will return the observation of the env at the current step, scaled and ordered according to the observation space

        returns:
        obs(array[floats]): returns an array of floats that is the observation
        (observing current muscle, current collagen, current agent observed cytokines, TGFb seperately)
        """

        # spatial data for muscle, collagen and necrosis
        self.muscle_spatial = SIM.getViableMuscle(self.ptrToEnv)
        self.collagen1_spatial = SIM.getCollagen1(self.ptrToEnv)
        self.collagen3_spatial = SIM.getCollagen3(self.ptrToEnv)
        self.necrosis_spatial = SIM.getNecrosis(self.ptrToEnv)

        # single value totals for mmuscle, collagen, necrosis
        self.previous_step_viable_muscle = self.current_viable_muscle
        self.current_viable_muscle = np.sum(self.muscle_spatial)
        self.current_collagen1 = np.sum(self.collagen1_spatial)
        self.current_collagen3 = np.sum(self.collagen3_spatial)
        self.current_necrosis = np.sum(self.necrosis_spatial)

        # for every step after the start, get current life, starting with 10000 life
        if self.current_step !=0:
            self.previous_step_life = self.current_life
            self.current_life = np.sum(np.array(SIM.getLife(self.ptrToEnv)))
        else:
            self.current_life = 10000

        # getting current surface cytokine values
        self.current_cytos = SIM.getSurfaceCytokineLevels(self.ptrToEnv)
        if self.current_step < MAX_STEPS:
            self.viable_muscle_history[self.current_step] = self.current_viable_muscle
            self.collagen1_history[self.current_step] = self.current_collagen1
            self.collagen3_history[self.current_step] = self.current_collagen3
            self.necrosis_history[self.current_step] = self.current_necrosis
            self.life_history[self.current_step] = self.current_life
            self.cytokine_history[self.current_step,:] = self.current_cytos
            self.norm_hist[self.current_step,:] = self.current_cytos / observed_signal_max


        # get and scaling observations
        obs = np.zeros(shape=(2+NUM_OBSERVABLE_CYTOKINES + int(OBSERVE_TGFb)))
        obs[0] = (self.current_viable_muscle - self.xDim*self.yDim*10) /(self.xDim*self.yDim*20)    #current viable muscle
        obs[1] = self.current_collagen3/(self.xDim*self.yDim*2)                                     #current collagen3
        # some logic to set up the observation correctly, and scale apropriately
        if SUM_M1_GENES:
            # if summing, only one number to worry about, so set it explicitly
            total = np.array(self.current_cytos)[np.array(OBSERVATION_CYTOKINES, dtype=bool)]
            obs[2] = np.sum(total)/signal_max
        else:
            # if not summing, fill up the rest of the array until the last one if observing TGFb, otherwise fill them all
            if OBSERVE_TGFb:
                obs[2:-1] = np.array(self.current_cytos)[np.array(OBSERVATION_CYTOKINES, dtype=bool)] / signal_max
            else:
                obs[2:] = np.array(self.current_cytos)[np.array(OBSERVATION_CYTOKINES, dtype=bool)] / signal_max

        if OBSERVE_TGFb:
            # if observing TGFb, add it to the end
            obs[-1] = self.current_cytos[8]/TGFb_max_obs
        return obs

    def calculate_reward(self):
        """ Calculating the reward from this current step
            returns:
            reward(float): reward for the current step based on the current state of the env vs last step
        """
        reward = (self.current_viable_muscle - self.previous_step_viable_muscle)
        if self.done ==2 or self.done ==3 or self.done==5:    #If pmn array hits 25k or Macro array hits 150k
            reward -=1000
        if self.done ==4 or self.done ==1:                    # normal finish or completely filled with muscle
            reward += (self.current_viable_muscle - self.xDim*self.yDim*10)
        return float(reward)

    def calculate_done(self):
        """ determine if the simulation should be done or not based on timesteps, and cell counts, and total collagen+muscle counts

        returns:
        (int): if the simulation should finish, and with what finish code. if not finished, returns False, othewise return int for finish code
        """
        if SIM.getNumPMN(self.ptrToEnv) > 25000: #10x the max in an uncontrolled run, kill it
            return 2
        if SIM.getNumMacrophage(self.ptrToEnv) > 150000: #10x the max in an uncontrolled run, kill it
            return 3
        if self.current_life >= (self.xDim*self.yDim) * 200 - 1000: #filled up the simulation with either muscle of collagen or both
            return 4
        if SIM.getNumFibroblast(self.ptrToEnv) > 150000: #10x the max in an uncontrolled run, kill it
            return 5
        if self.current_step >= MAX_STEPS-1:    # run out of time, kill it
            return 1
        return False


    def reset(self):
        """ Resets the simulation so it is ready to run again from the start (eg start of a new episode)

        returns:
        observation (array[floats]): initial observation of the simulation state after reset
        """
        self.seed+=1
        self.ptrToEnv = createWEABM(self.seed)
        self.xDim = SIM.get_xDim(self.ptrToEnv)
        self.yDim = SIM.get_yDim(self.ptrToEnv)
        done = False
        self.current_step = 0                                       # The current step of the simulation environmen
        self.current_cytos = np.zeros(NUM_OBSERVABLE_CYTOKINES)
        self.current_viable_muscle = np.sum(SIM.getViableMuscle(self.ptrToEnv))
        self.previous_step_life = 9900
        self.previous_step_viable_muscle = self.current_viable_muscle

        self.reset_history()

        return self.next_observation()

    def reset_history(self):
        """ resets the history trackers for the new episode"""

        self.observation_history = np.zeros(shape=(MAX_STEPS,2+NUM_OBSERVABLE_CYTOKINES+ int(OBSERVE_TGFb)))
        # self.observation_history = np.zeros(shape=(MAX_STEPS,3))        # M1 sum experiment
        self.reward_history = np.zeros(MAX_STEPS)
        self.viable_muscle_history = np.zeros(MAX_STEPS)
        self.collagen1_history = np.zeros(MAX_STEPS)
        self.collagen3_history = np.zeros(MAX_STEPS)
        self.necrosis_history = np.zeros(MAX_STEPS)
        self.life_history = np.zeros(MAX_STEPS)
        self.cytokine_history = np.zeros(shape=(MAX_STEPS, TOTAL_OBSERVABLE_CYTOKINES))
        self.norm_hist = np.zeros(shape=(MAX_STEPS, TOTAL_OBSERVABLE_CYTOKINES))
        self.action_history = np.zeros(shape=(MAX_STEPS, TOTAL_CONTROLLABLE_CYTOKINES))


    def render(self,kward=None):
        """console rendering for every step. renders current cytokines, actions taken, cell counts, and tissue counts"""
        if self.current_step == 0:
            return
        print("Stepping Simulation, Step: " + str(self.current_step) + "              ")
        print()
        print("Surface Cytokines: ")
        print("IL8:{: 1.2f}, MIP1b:{: 1.2f}, TNF:{: 1.2f}, HGF:{: 1.2f}, PAF:{: 1.2f}, IL1:{: 1.2f}, PDGF:{: 1.2f}, IL13:{: 1.2f}, TGFb:{: 1.2f} ".format(self.current_cytos[0],self.current_cytos[1],self.current_cytos[2],self.current_cytos[3],self.current_cytos[4],self.current_cytos[5],self.current_cytos[6],self.current_cytos[7],self.current_cytos[8]))
        print("VEGF:{: 1.2f}, GCSF:{: 1.2f}, IL10:{: 1.2f}, IL4:{: 1.2f}, IL17:{: 1.2f}, Myf5:{: 1.2f}, IL6:{: 1.2f}, IL1ra:{: 1.2f}, sIL1r:{: 1.2f}, sTNFr:{: 1.2f} ".format(self.current_cytos[9],self.current_cytos[10],self.current_cytos[11],self.current_cytos[12],self.current_cytos[13],self.current_cytos[14],self.current_cytos[15],self.current_cytos[16],self.current_cytos[17],self.current_cytos[18]))
        print("endotoxin:{: 1.2f}, IFNg:{: 1.2f}, IL12:{: 1.2f}, MRF4:{: 1.2f}, cytotox:{: 1.2f}, IGF:{: 1.2f}, DAMP:{: 1.2f}, MCP1:{: 1.2f}, FNE:{: 1.2f} ".format(self.current_cytos[19],self.current_cytos[20],self.current_cytos[21],self.current_cytos[22],self.current_cytos[23],self.current_cytos[24],self.current_cytos[25],self.current_cytos[26], self.current_cytos[27]))

        print()
        print("Actions: ")
        print("IL8:{: 1.2f}, MIP1b:{: 1.2f}, TNF:{: 1.2f}, HGF:{: 1.2f}, PAF:{: 1.2f}, IL1:{: 1.2f}, PDGF:{: 1.2f}, IL13:{: 1.2f}, TGFb:{: 1.2f} ".format(self.current_action[0],self.current_action[1],self.current_action[2],self.current_action[3],self.current_action[4],self.current_action[5],self.current_action[6],self.current_action[7],self.current_action[8]))
        print("VEGF:{: 1.2f}, GCSF:{: 1.2f}, IL10:{: 1.2f}, IL4:{: 1.2f}, IL17:{: 1.2f}, Myf5:{: 1.2f}, IL6:{: 1.2f}, IL1ra:{: 1.2f}, sIL1r:{: 1.2f}, sTNFr:{: 1.2f} ".format(self.current_action[9],self.current_action[10],self.current_action[11],self.current_action[12],self.current_action[13],self.current_action[14],self.current_action[15],self.current_action[16],self.current_action[17],self.current_action[18]))
        print("endotoxin:{: 1.2f}, IFNg:{: 1.2f}, IL12:{: 1.2f}, MRF4:{: 1.2f}, cytotox:{: 1.2f}, IGF:{: 1.2f}, DAMP:{: 1.2f}, MCP1:{: 1.2f}, FNE:{: 1.2f}   ".format(self.current_action[19],self.current_action[20],self.current_action[21],self.current_action[22],self.current_action[23],self.current_action[24],self.current_action[25],self.current_action[26], self.current_action[27]))
        print()
        print("Viable Muscle: {:d}       ".format(int(self.current_viable_muscle)))
        print("Total Collagen1: {:d}       ".format(int(self.current_collagen1)))
        print("Total Collagen3: {:d}       ".format(int(self.current_collagen3)))
        print("Total Necrosis: {:d}       ".format(int(self.current_necrosis)))
        print("Total Life: {:d}          ".format(int(self.current_life)))
        print("Num PMNs: {:d}       ".format(SIM.getNumPMN(self.ptrToEnv)))
        print("Num Macrophages: {:d}      ".format(SIM.getNumMacrophage(self.ptrToEnv)))
        print("Num Fibroblasts: {:d}       ".format(SIM.getNumFibroblast(self.ptrToEnv)))
        print("Num Satellites: {:d}       ".format(SIM.getNumSatellite(self.ptrToEnv)))
        print("Num Myoblasts: {:d}       ".format(SIM.getNumMyoblast(self.ptrToEnv)))
        print("Num T0: {:d}, T1: {:d}, T2: {:d}, T17: {:d}".format(SIM.getNumTH0(self.ptrToEnv), SIM.getNumTH1(self.ptrToEnv), SIM.getNumTH2(self.ptrToEnv), SIM.getNumTH17(self.ptrToEnv)))
        np.set_printoptions(precision=2, suppress=True)
        print(self.next_observation(),end="\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F")
