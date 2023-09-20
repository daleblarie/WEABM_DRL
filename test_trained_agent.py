import gym
import sys
from stable_baselines3 import DDPG
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
import numpy as np
from mpi4py import MPI
from callbacks import CustomCallback
from WEABM_Env_DRL import WEABM_Environment

LOAD_LOCATION = "WEABM_DRL_Experiments/Experiment1"
SAVE_LOCATION="WEABM_DRL_Tests/Experiment1"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# specifying save location, should save at the end of every episode when done==True
env = WEABM_Environment(action_scaling_factor=1, SAVE_LOCATION=SAVE_LOCATION)
n_actions = env.action_space.shape[-1]
action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))

print("env created")
try:
    model = DDPG.load(LOAD_LOCATION + "/" + str(rank) + "_Best_Agent", env=env)
except:
    try:
        model = DDPG.load(LOAD_LOCATION + "/" + str(rank) + "_Agent", env=env)
    except:
        print("No pretrained agents were found at " + LOAD_LOCATION +". Exiting now")
        sys.exit(1)


for episode in range(5):
    print(rank, episode)
    env.setSeed(episode)
    obs = env.reset()
    done = False
    while not done:
        action, _states = model.predict(obs)
        obs, rewards, done, info = env.step(action)
        # env.render()
