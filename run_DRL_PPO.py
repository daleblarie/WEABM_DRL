import gym
import sys
from stable_baselines3 import DDPG, PPO
from stable_baselines3.common.vec_env import SubprocVecEnv, DummyVecEnv

from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
import numpy as np
from mpi4py import MPI

from callbacks import CustomCallback
from WEABM_Env_DRL import WEABM_Environment



SAVE_LOCATION="WEABM_DRL_Experiments/PPO1"

# Create vectorized environment for parallel processing
num_cpu = 4  # Number of processes to use

def make_env(rank, seed=0):
    def _init():
        return WEABM_Environment(action_scaling_factor=1, rank=rank, seed=seed, SAVE_LOCATION=SAVE_LOCATION)
    return _init


if __name__ =="__main__":
    env = SubprocVecEnv([make_env(i, seed=i*1000) for i in range(num_cpu)])  

    print("env created")

    try:
        model = PPO.load(SAVE_LOCATION + "/" + str(rank) + "_Best_Agent", env=env)
    except:
        try:
            model = PPO.load(SAVE_LOCATION + "/" + str(rank) + "_Agent", env=env)
        except:
            model = PPO('MlpPolicy', env,  device="cpu")

    print("model created")

    # callback will enable rendering and do some data saving for spatial data and reward histories and such
    cb = CustomCallback(SAVE_FILE_NAME=SAVE_LOCATION, rendering=False)
    print("CB created")

    model.learn(total_timesteps=1401100, callback=cb)
    model.save(SAVE_LOCATION + "/Agent")