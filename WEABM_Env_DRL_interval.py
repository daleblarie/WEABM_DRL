# WEABM DRL Environment-V2, interval
# Created by Dale Larie 9/9/2022

# this version will allow the agent full control, but only at limited intervals. It can vary its actions while it is turned on, but when it is turned off it inputs no control
# change control hours and uncontrol hours to determine the control intervals

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import wrapper_setup
from mpi4py import MPI

import gym

SENSE_NO = True


cytokine_labels = ["IL8", "MIP1b", "TNF", "HGF", "PAF", "IL1", "PDGF", "IL13", "TGFb", "VEGF", "GCSF", "IL10", "IL4", "IL17", "Myf5", "IL6", "IL1ra", "sIL1r", "sTNFr", "endotoxin", "IFNg", "IL12", "MRF4", "cytotox", "IGF", "DAMP", "MCP1", "FNE"]
MAX_STEPS = 1400               #Max number of steps the simulation will play
NUM_OBSERVTAIONS = 1             #Number of steps in the past the agent is allowed to see
MINUTES_PER_STEP = 15           #15 minutes of real time per simulation step
# observation cytokines are only individual cytokines. this does not affect whether aggregated cytokines (NO) are observed
# OBSERVATION_CYTOKINES = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
OBSERVATION_CYTOKINES = [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]       # Looking at only TGFb
# OBSERVATION_CYTOKINES = [0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0]       # Looking at only M1 Genes

NUM_OBSERVABLE_CYTOKINES = sum(OBSERVATION_CYTOKINES) + int(SENSE_NO)
TOTAL_OBSERVABLE_CYTOKINES = 28         #Number of cytokines in the WEABM

TOTAL_CONTROLLABLE_CYTOKINES = 28       #Number of cytokines that can be controlled
# ACTION_INDECES = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
ACTION_INDECES = [0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]  # using {TGFb, TNF, IL10, HGF, cytotox}
# ACTION_INDECES = [0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  #using only TNF, IL10 and IL4
# ACTION_INDECES = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  # totally uncontrolled
NUM_CYTOKINES_CONTROLLED = sum(ACTION_INDECES)   #Number of cytokines controlled by the agent

# empirically observed maximum values for each cytokine with random control
signal_max = [6507.83349609,   51.18336487, 1526.96228027 ,  56.05470657 ,1716.56884766,
   65.22439575,  105.38170624, 2129.74145508, 2632.51220703,  100.20832825,
 1759.00390625, 1949.2902832 , 2390.04394531,  640.81066895 , 113.35812378,
 1698.16662598,   53.98086929,   52.61504745,   97.77040863  , 97.39041138,
 2115.72216797, 2548.54370117,   95.35305786,  283.16220093,   95.18818665,
 2248.67578125, 2238.12207031,  939.54016113]
signal_max = [signal_max[i] for i in range(len(OBSERVATION_CYTOKINES)) if OBSERVATION_CYTOKINES[i]==1]
signal_max = sum(signal_max)    #when summing M1 genes

signal_min = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
signal_min = [signal_min[i] for i in range(len(OBSERVATION_CYTOKINES)) if OBSERVATION_CYTOKINES[i]==1]
signal_min = sum(signal_min)    #when summing M1 genes

if not (isinstance(signal_min, int) and isinstance(signal_max, float)):
    signal_range = [a_i - b_i for a_i, b_i in zip(signal_max, signal_min)]

SIM = wrapper_setup.setUpWrapper()



def createWEABM(seed):

    IP=np.genfromtxt('RuleMat3.csv',delimiter=',')
    #np.save('baseParameterization1.npy',x)
    #IP=np.load('baseParameterization1.npy')
    IP=IP.flatten()
    numMatrixElements=IP.shape[0]
    array_type = ctypes.c_float*numMatrixElements


    instance = SIM.CreateInstance(array_type(*IP))
    SIM.setSeed(instance, seed)
    SIM.InitializeSimulation(instance)
    SIM.DoInjury(instance)
    return instance

class WEABM_Environment(gym.Env):

    def __init__(self, SAVE_LOCATION="./", action_scaling_factor=1, interval_training=True):
        super(WEABM_Environment, self).__init__()
        self.action_scaling_factor = action_scaling_factor
        self.interval_training = interval_training
        self.save_path = SAVE_LOCATION
        self.MAX_STEPS = MAX_STEPS
        self.seed = 0
        self.current_step = 0                                       # The current step of the simulation environment
        self.RL_step = 0                                            # The current step of the agent
        self.control_hours = 6
        self.control_steps = self.control_hours*60/MINUTES_PER_STEP
        self.uncontrol_hours = 42
        self.uncontrol_steps = self.uncontrol_hours*60/MINUTES_PER_STEP
        self.full_interval = self.control_steps + self.uncontrol_steps
        self.current_cytos = np.zeros(TOTAL_OBSERVABLE_CYTOKINES)
        self.current_viable_muscle = 0
        self.current_collagen1 = 0
        self.current_collagen3 = 0
        self.current_necrosis = 0
        self.current_life = 0
        self.previous_step_life = 0
        self.xDim = 0       #<--- set in reset
        self.ydim = 0       #<---

        self.observation_history = np.zeros(shape=(MAX_STEPS,2+NUM_OBSERVABLE_CYTOKINES)) #normal env
        # self.observation_history = np.zeros(shape=(MAX_STEPS,3))        # m1 sum experiment
        self.reward_history = np.zeros(MAX_STEPS)
        self.viable_muscle_history = np.zeros(MAX_STEPS)
        self.collagen1_history = np.zeros(MAX_STEPS)
        self.collagen3_history = np.zeros(MAX_STEPS)
        self.necrosis_history = np.zeros(MAX_STEPS)
        self.life_history = np.zeros(MAX_STEPS)
        self.cytokine_history = np.zeros(shape=(MAX_STEPS, TOTAL_OBSERVABLE_CYTOKINES))
        self.norm_hist = np.zeros(shape=(MAX_STEPS, TOTAL_OBSERVABLE_CYTOKINES))
        self.action_history = np.zeros(shape=(MAX_STEPS, TOTAL_CONTROLLABLE_CYTOKINES))

        self.allSignalsReturn_output = np.zeros(43*MAX_STEPS)

        self.current_action = np.zeros(TOTAL_CONTROLLABLE_CYTOKINES)

        self.added_cytokines = np.zeros(TOTAL_CONTROLLABLE_CYTOKINES)

        self.action_space = gym.spaces.Box(
            low=-1,
            high=1,
            shape=(NUM_CYTOKINES_CONTROLLED,),
            dtype=np.float32)

        # observation space is TotalViableMuscle, TotalCollagen3 and 27 cytokine values
        obs_space_high = np.zeros(2+NUM_OBSERVABLE_CYTOKINES)
        obs_space_high[0] = 30*30*20
        obs_space_high[1] = 30*30*2
        obs_space_high[2:] = np.array(signal_max)
        obs_max = obs_space_high

        obs_space_low = np.zeros(2+NUM_OBSERVABLE_CYTOKINES)
        obs_space_low[0] = 0
        obs_space_low[1] = 0
        obs_space_low[2:] = np.array(signal_min)
        obs_min = obs_space_low

        # stack the observation space with the number of frames to be observed
        for i in range(NUM_OBSERVTAIONS-1):
            obs_space_high = np.vstack((obs_space_high,obs_max))
        obs_space_high = obs_space_high.T.flatten()


        self.observation_space = gym.spaces.Box(
            low=obs_space_low,
            high=obs_space_high,
            shape=obs_space_high.shape,
            dtype=np.float32)

        # # test obs space for summing M1 genes
        # self.observation_space = gym.spaces.Box(
        #     low=np.array([0,0,-1]),
        #     high=np.array([30*30*30, 30*30*2, 100000]),
        #     shape=(3,),
        #     dtype=np.float32)


        # call reset to initialze all the variables so it is ready to simulatie
        self.reset()


    def setSeed(self, new_seed):
        self.seed = new_seed
        SIM.setSeed(self.ptrToEnv, new_seed)
        # pass

    def interval_training_step(self):
        if self.interval_training:
            step_in_interval = self.current_step % self.full_interval
            if self.current_step > self.control_steps:
                blank_action = np.zeros(shape=NUM_CYTOKINES_CONTROLLED)
                self.step(blank_action)


    def step(self, action):
        if self.current_step == 0:
            self.reset_history()
        self.RL_step += 1
        self.take_action(action)

        observation = self.next_observation()

        done = self.calculate_done()
        self.done = done
        reward = self.calculate_reward()
        self.reward = reward

        done = bool(done)
        if done:
            self.allSignalsReturn_output = SIM.endSimulation(self.ptrToEnv)
            self.save_training_data()

        if self.current_step < MAX_STEPS:
            self.observation_history[self.current_step,:] = observation
            self.reward_history[self.current_step] = reward


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

        if self.interval_training and not self.done:
            step_in_interval = self.current_step % self.full_interval
            if step_in_interval > self.control_steps:
                blank_action = np.zeros(shape=NUM_CYTOKINES_CONTROLLED)
                _, self.reward, self.done, _ = self.step(blank_action)

        return observation, reward, done, info

    def save_training_data(self):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        with open(self.save_path+"/" + str(rank) + "_training_data.csv", "ab") as f:
            saveArray = self.viable_muscle_history
            saveArray = np.vstack((saveArray, self.collagen3_history))
            saveArray = np.vstack((saveArray, self.necrosis_history))
            saveArray = np.vstack((saveArray, self.collagen1_history))
            saveArray = np.vstack((saveArray, self.norm_hist.T))
            saveArray = np.vstack((saveArray, self.action_history.T))
            np.savetxt(f, saveArray, delimiter=",")


    def take_action(self, action):
        action_num = 0
        action *= self.action_scaling_factor
        # action +=1
        # action/=2
        if ACTION_INDECES[0] != 0:
            SIM.apply_IL8(self.ptrToEnv, action[action_num])
            # print("Applied IL8")
            self.current_action[0] = action[action_num]
            self.added_cytokines[0] = action[action_num]
            action_num += 1
        if ACTION_INDECES[1] != 0:
            SIM.apply_MIP1b(self.ptrToEnv, action[action_num])
            # print("Applied mip1b")
            self.current_action[1] = action[action_num]
            self.added_cytokines[1] = action[action_num]
            action_num += 1
        if ACTION_INDECES[2] != 0:
            SIM.apply_TNF(self.ptrToEnv, action[action_num])
            # print("Applied tnf")
            self.current_action[2] = action[action_num]
            self.added_cytokines[2] = action[action_num]
            action_num += 1
        if ACTION_INDECES[3] != 0:
            SIM.apply_HGF(self.ptrToEnv, action[action_num])
            # print("Applied hgf")
            self.current_action[3] = action[action_num]
            self.added_cytokines[3] = action[action_num]
            action_num += 1
        if ACTION_INDECES[4] != 0:
            SIM.apply_PAF(self.ptrToEnv, action[action_num])
            # print("Applied paf")
            self.current_action[4] = action[action_num]
            self.added_cytokines[4] = action[action_num]
            action_num += 1
        if ACTION_INDECES[5] != 0:
            SIM.apply_IL1(self.ptrToEnv, action[action_num])
            # print("Applied IL1")
            self.current_action[5] = action[action_num]
            self.added_cytokines[5] = action[action_num]
            action_num += 1
        if ACTION_INDECES[6] != 0:
            SIM.apply_PDGF(self.ptrToEnv, action[action_num])
            # print("Applied pdgf")
            self.current_action[6] = action[action_num]
            self.added_cytokines[6] = action[action_num]
            action_num += 1
        if ACTION_INDECES[7] != 0:
            SIM.apply_IL13(self.ptrToEnv, action[action_num])
            # print("Applied IL13")
            self.current_action[7] = action[action_num]
            self.added_cytokines[7] = action[action_num]
            action_num += 1
        if ACTION_INDECES[8] != 0:
            SIM.apply_TGFb(self.ptrToEnv, action[action_num])
            # print("Applied tgfb")
            self.current_action[8] = action[action_num]
            self.added_cytokines[8] = action[action_num]
            action_num += 1
        if ACTION_INDECES[9] != 0:
            SIM.apply_VEGF(self.ptrToEnv, action[action_num])
            # print("Applied vegf")
            self.current_action[9] = action[action_num]
            self.added_cytokines[9] = action[action_num]
            action_num += 1
        if ACTION_INDECES[10] != 0:
            SIM.apply_GCSF(self.ptrToEnv, action[action_num])
            # print("Applied gcsf")
            self.current_action[10] = action[action_num]
            self.added_cytokines[10] = action[action_num]
            action_num += 1
        if ACTION_INDECES[11] != 0:
            SIM.apply_IL10(self.ptrToEnv, action[action_num])
            # print("Applied IL10")
            self.current_action[11] = action[action_num]
            self.added_cytokines[11] = action[action_num]
            action_num += 1
        if ACTION_INDECES[12] != 0:
            SIM.apply_IL4(self.ptrToEnv, action[action_num])
            # print("Applied IL4")
            self.current_action[12] = action[action_num]
            self.added_cytokines[12] = action[action_num]
            action_num += 1
        if ACTION_INDECES[13] != 0:
            SIM.apply_IL17(self.ptrToEnv, action[action_num])
            # print("Applied IL17")
            self.current_action[13] = action[action_num]
            self.added_cytokines[13] = action[action_num]
            action_num += 1
        if ACTION_INDECES[14] !=0:
            SIM.apply_Myf5(self.ptrToEnv, action[action_num])
            # print("Applied myf5")
            self.current_action[14] = action[action_num]
            self.added_cytokines[14] = action[action_num]
            action_num += 1
        if ACTION_INDECES[15] !=0:
            SIM.apply_IL6(self.ptrToEnv, action[action_num])
            # print("Applied IL6")
            self.current_action[15] = action[action_num]
            self.added_cytokines[15] = action[action_num]
            action_num += 1
        if ACTION_INDECES[16] !=0:
            SIM.apply_IL1ra(self.ptrToEnv, action[action_num])
            # print("Applied IL1ra")
            self.current_action[16] = action[action_num]
            self.added_cytokines[16] = action[action_num]
            action_num += 1
        if ACTION_INDECES[17] !=0:
            SIM.apply_sIL1r(self.ptrToEnv, action[action_num])
            # print("Applied sil1r")
            self.current_action[17] = action[action_num]
            self.added_cytokines[17] = action[action_num]
            action_num += 1
        if ACTION_INDECES[18] !=0:
            SIM.apply_sTNFr(self.ptrToEnv, action[action_num])
            # print("Applied stnfr")
            self.current_action[18] = action[action_num]
            self.added_cytokines[18] = action[action_num]
            action_num += 1
        if ACTION_INDECES[19] !=0:
            SIM.apply_endotoxin(self.ptrToEnv, action[action_num])
            # print("Applied endo")
            self.current_action[19] = action[action_num]
            self.added_cytokines[19] = action[action_num]
            action_num += 1
        if ACTION_INDECES[20] !=0:
            SIM.apply_IFNg(self.ptrToEnv, action[action_num])
            # print("Applied ifng")
            self.current_action[20] = action[action_num]
            self.added_cytokines[20] = action[action_num]
            action_num += 1
        if ACTION_INDECES[21] !=0:
            SIM.apply_IL12(self.ptrToEnv, action[action_num])
            # print("Applied il12")
            self.current_action[21] = action[action_num]
            self.added_cytokines[21] = action[action_num]
            action_num += 1
        if ACTION_INDECES[22] !=0:
            SIM.apply_MRF4(self.ptrToEnv, action[action_num])
            # print("Applied mrf4")
            self.current_action[22] = action[action_num]
            self.added_cytokines[22] = action[action_num]
            action_num += 1
        if ACTION_INDECES[23] !=0:
            SIM.apply_cytotox(self.ptrToEnv, action[action_num])
            # print("Applied cytotox")
            self.current_action[23] = action[action_num]
            self.added_cytokines[23] = action[action_num]
            action_num += 1
        if ACTION_INDECES[24] !=0:
            SIM.apply_IGF(self.ptrToEnv, action[action_num])
            # print("Applied igf")
            self.current_action[24] = action[action_num]
            self.added_cytokines[24] = action[action_num]
            action_num += 1
        if ACTION_INDECES[25] !=0:
            SIM.apply_DAMP(self.ptrToEnv, action[action_num])
            # print("Applied damp")
            self.current_action[25] = action[action_num]
            self.added_cytokines[25] = action[action_num]
            action_num += 1
        if ACTION_INDECES[26] !=0:
            SIM.apply_MCP1(self.ptrToEnv, action[action_num])
            # print("Applied mcp1")
            self.current_action[26] = action[action_num]
            self.added_cytokines[26] = action[action_num]
            action_num += 1
        if ACTION_INDECES[27] !=0:
            SIM.apply_FNE(self.ptrToEnv, action[action_num])
            # print("Applied mcp1")
            self.current_action[27] = action[action_num]
            self.added_cytokines[27] = action[action_num]
            action_num += 1

        self.action_history[self.current_step,:] = self.current_action

        SIM.StepSimulation(self.ptrToEnv)
        # self.adjust_tracked_cytos()
        self.current_step+=1

    def adjust_tracked_cytos(self):
        # This implementation ignores cytokine that diffuses back up from the 2nd layer.
        # DAMP is the cytokine where the largest amount would come back up,
        # with 1 unit added, the amount that diffuses back up the next turn is .008, 0.8% - small enough to ignore for now

        # IL8
        # diffusing out of central cell, and tracking the amount that stays on the top layer.
        # If it goes up it leaves the sim, if it goes down it is not read in the observation
        self.added_cytokines[0] = self.added_cytokines[0]*0.4 + (8/26 * self.added_cytokines[0]*0.6)
        # now evaporate
        self.added_cytokines[0] = self.added_cytokines[0]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))

        # TNF
        self.added_cytokines[2] = self.added_cytokines[2]*0.4 + (8/26 * self.added_cytokines[2]*0.6)
        self.added_cytokines[2] = self.added_cytokines[2]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # PAF
        self.added_cytokines[4] = self.added_cytokines[4]*0.4 + (8/26 * self.added_cytokines[4]*0.6)
        self.added_cytokines[4] = self.added_cytokines[4]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # IL1
        self.added_cytokines[5] = self.added_cytokines[5]*0.4 + (8/26 * self.added_cytokines[5]*0.6)
        self.added_cytokines[5] = self.added_cytokines[5]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # TGFb
        self.added_cytokines[8] = self.added_cytokines[8]*0.4 + (8/26 * self.added_cytokines[8]*0.6)
        self.added_cytokines[8] = self.added_cytokines[8]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # IL10
        self.added_cytokines[11] = self.added_cytokines[11]*0.4 + (8/26 * self.added_cytokines[11]*0.6)
        self.added_cytokines[11] = self.added_cytokines[11]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # IL4
        self.added_cytokines[12] = self.added_cytokines[12]*0.6 + (8/26 * self.added_cytokines[12]*0.4)
        self.added_cytokines[12] = self.added_cytokines[12]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # IL1ra
        self.added_cytokines[16] = self.added_cytokines[16]*0.4 + (8/26 * self.added_cytokines[16]*0.6)
        self.added_cytokines[16] = self.added_cytokines[16]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # sIL1r
        self.added_cytokines[17] = self.added_cytokines[17]*0.4 + (8/26 * self.added_cytokines[17]*0.6)
        self.added_cytokines[17] = self.added_cytokines[17]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # IFNg
        self.added_cytokines[20] = self.added_cytokines[20]*0.4 + (8/26 * self.added_cytokines[20]*0.6)
        self.added_cytokines[20] = self.added_cytokines[20]*0.1* (1-self.current_necrosis/(self.xDim*self.yDim))
        # cytotox
        self.added_cytokines[23] = self.added_cytokines[23]*0.25 + (8/26 * self.added_cytokines[23]*0.75)
        self.added_cytokines[23] = self.added_cytokines[23]*0.01* (1-self.current_necrosis/(self.xDim*self.yDim))
        # DAMP
        self.added_cytokines[25] = self.added_cytokines[25]*0.5 + (8/26 * self.added_cytokines[25]*0.5)
        self.added_cytokines[25] = self.added_cytokines[25]*0.3* (1-self.current_necrosis/(self.xDim*self.yDim))
        # MCP1
        self.added_cytokines[26] = self.added_cytokines[26]*0.4 + (8/26 * self.added_cytokines[26]*0.6)
        self.added_cytokines[26] = self.added_cytokines[26]*0.1 * (1-self.current_necrosis/(self.xDim*self.yDim))

    def sense_NO(self):
        surface_vals = SIM.getSurfaceCytokineLevels(self.ptrToEnv)
        NO_cytokines = [2,5,15,16,20,21]
        NO = np.sum(surface_vals[NO_cytokines])
        return NO

    def next_observation(self):
        self.muscle_spatial = SIM.getViableMuscle(self.ptrToEnv)
        self.collagen1_spatial = SIM.getCollagen1(self.ptrToEnv)
        self.collagen3_spatial = SIM.getCollagen3(self.ptrToEnv)
        self.necrosis_spatial = SIM.getNecrosis(self.ptrToEnv)

        self.current_collagen1 = np.sum(self.collagen1_spatial)
        self.current_collagen3 = np.sum(self.collagen3_spatial)
        self.current_necrosis = np.sum(self.necrosis_spatial)

        self.previous_step_viable_muscle = self.current_viable_muscle
        self.current_viable_muscle = np.sum(self.muscle_spatial)
        if self.current_step !=0:
            self.previous_step_life = self.current_life
            self.current_life = np.sum(np.array(SIM.getLife(self.ptrToEnv)))
        else:
            self.current_life = 10000
        self.current_cytos = SIM.getSurfaceCytokineLevels(self.ptrToEnv)
        if self.current_step < MAX_STEPS:
            self.viable_muscle_history[self.current_step] = self.current_viable_muscle
            self.collagen1_history[self.current_step] = self.current_collagen1
            self.collagen3_history[self.current_step] = self.current_collagen3
            self.necrosis_history[self.current_step] = self.current_necrosis
            self.life_history[self.current_step] = self.current_life
            self.cytokine_history[self.current_step,:] = self.current_cytos
            self.norm_hist[self.current_step,:] = self.current_cytos# - self.added_cytokines*(self.xDim*self.yDim)


        if SENSE_NO:
            # M1 summed experiment obs
            obs = np.zeros(shape=(2+NUM_OBSERVABLE_CYTOKINES))
            obs[0] = (self.current_viable_muscle - self.xDim*self.yDim*10) /(self.xDim*self.yDim*20)
            obs[1] = self.current_collagen3/(self.xDim*self.yDim*2)
            obs[2] = self.sense_NO()
            obs[3:] = np.array(self.current_cytos)[np.array(OBSERVATION_CYTOKINES, dtype=bool)]
        else:
            obs = np.zeros(shape=(2+NUM_OBSERVABLE_CYTOKINES))
            # obs[0] = self.current_viable_muscle/25000
            # obs[1:] = np.divide(self.current_cytos, signal_max)
            obs[0] = (self.current_viable_muscle - self.xDim*self.yDim*10) /(self.xDim*self.yDim*20)
            obs[1] = self.current_collagen3/(self.xDim*self.yDim*2)
            obs[2:] = np.array(self.current_cytos)[np.array(OBSERVATION_CYTOKINES, dtype=bool)]# - self.added_cytokines*(self.xDim*self.yDim)
            # obs[2:] = obs[2:] - np.array(signal_min)[np.array(OBSERVATION_CYTOKINES, dtype=bool)]
            # obs[2:] = 2*np.divide(obs[2:], np.array(signal_range)[np.array(OBSERVATION_CYTOKINES, dtype=bool)])-1
            # # obs = np.clip(obs, -1, 1)
            # # print(obs)
        return obs

    def calculate_reward(self):
        reward = (self.current_viable_muscle - self.previous_step_viable_muscle)
        # reward=  0
        if self.done ==2 or self.done ==3 or self.done==5:    #If pmn array hits 25k or Macro array hits 150k
            reward -=1000
        if self.done ==4 or self.done ==1:                    # normal finish or completely filled with muscle
            reward += (self.current_viable_muscle - self.xDim*self.yDim*10)
        return float(reward)

    def calculate_done(self):
        if SIM.getNumPMN(self.ptrToEnv) > 25000: #10x the max in an uncontrolled run
            return 2
        if SIM.getNumMacrophage(self.ptrToEnv) > 150000: #10x the max in an uncontrolled run
            return 3
        if self.current_life >= (self.xDim*self.yDim) * 200 - 1000:
            return 4
        if SIM.getNumFibroblast(self.ptrToEnv) > 150000: #10x the max in an uncontrolled run
            return 5
        if self.current_step >= MAX_STEPS-1:
            return 1
        return False


    def reset(self):
        # self.seed+=1
        # print("SEED " + str(self.seed))
        self.ptrToEnv = createWEABM(self.seed)
        self.xDim = SIM.get_xDim(self.ptrToEnv)
        self.yDim = SIM.get_yDim(self.ptrToEnv)
        done = False
        self.current_step = 0                                       # The current step of the simulation environmen
        self.RL_step = 0                                            # The current step of the agent
        self.current_cytos = np.zeros(NUM_OBSERVABLE_CYTOKINES)
        self.current_viable_muscle = np.sum(SIM.getViableMuscle(self.ptrToEnv))
        self.previous_step_life = 9900
        self.previous_step_viable_muscle = self.current_viable_muscle
        self.added_cytokines = np.zeros(TOTAL_CONTROLLABLE_CYTOKINES)

        return self.next_observation()

    def reset_history(self):
        self.observation_history = np.zeros(shape=(MAX_STEPS,2+NUM_OBSERVABLE_CYTOKINES))
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
        print(self.next_observation(),end="\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F\033[F")
