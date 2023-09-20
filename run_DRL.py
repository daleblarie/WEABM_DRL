import gym
import sys
from stable_baselines3 import DDPG
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
import numpy as np
from mpi4py import MPI

from callbacks import CustomCallback
from WEABM_Env_DRL import WEABM_Environment



SAVE_LOCATION="WEABM_DRL_Experiments/Experiment1"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# if you specify SAVE_LOCATION in WEABM_Environment then then environment will save only the training data after each episode
# otherwise, create and specify the CustomCallback and that will take care of saving for you
env = WEABM_Environment(action_scaling_factor=1)
n_actions = env.action_space.shape[-1]
action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))

print("env created")
try:
    model = DDPG.load(SAVE_LOCATION + "/" + str(rank) + "_Best_Agent", env=env)
except:
    try:
        model = DDPG.load(SAVE_LOCATION + "/" + str(rank) + "_Agent", env=env)
    except:
        model = DDPG('MlpPolicy', env, train_freq=(1,"episode"), learning_starts=10, action_noise = action_noise, device="cpu")

print("model created")
# callback will enable rendering and do some data saving for spatial data and reward histories and such
cb = CustomCallback(SAVE_FILE_NAME=SAVE_LOCATION, rendering=False)
print("CB created")
model.learn(total_timesteps=1401100, callback=cb)
model.save(SAVE_LOCATION + "/Agent")
