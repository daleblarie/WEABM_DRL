import gym
import sys
from stable_baselines3 import DDPG
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
import numpy as np
from mpi4py import MPI

from callbacks_interval import CustomCallback
from WEABM_Env_DRL_interval import WEABM_Environment


SAVE_LOCATION="WEABM_DRL_Experiments/Interval_1"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# specifying save location in the environment because of interval saving- it gets weirder than just handing everyting to the callback handler
env = WEABM_Environment(action_scaling_factor=10, SAVE_LOCATION=SAVE_LOCATION)
n_actions = env.action_space.shape[-1]
action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))

print("env created")
try:
    model = DDPG.load(SAVE_LOCATION + "/" + str(rank) + "_Best_Agent", env=env)
except:
    try:
        model = DDPG.load(SAVE_LOCATION + "/" + str(rank) + "_Agent", env=env)
    except:
        model = DDPG('MlpPolicy', env, train_freq=(1,"episode"), learning_starts=1000, action_noise = action_noise, device="cpu")

print("model created")
cb = CustomCallback(SAVE_FILE_NAME=SAVE_LOCATION, rendering=False)
print("CB created")
model.learn(total_timesteps=1401100, callback=cb)
model.save(SAVE_LOCATION + "/Agent")
