import gym
import sys
from stable_baselines3 import DDPG
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
import numpy as np
from mpi4py import MPI

from callbacks_interval import CustomCallback
from WEABM_Env_DRL_interval import WEABM_Environment

# sys.path.insert(1,"WEABM_DRL_Experiments")
# from WEABM_Env_Summed_m1 import WEABM_Environment as WEABM_sum_m1


Model_SAVE_LOCATION="WEABM_DRL_Experiments/Interval_1"
Test_data_SAVE_LOCATION="WEABM_DRL_Tests/Interval_1"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

env = WEABM_Environment(action_scaling_factor=10, SAVE_LOCATION=Test_data_SAVE_LOCATION)
n_actions = env.action_space.shape[-1]
action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))

print("env created")
try:
    model = DDPG.load(Model_SAVE_LOCATION + "/" + str(rank) + "_Best_Agent.zip", env=env, device="cpu")
    print("loaded " + str(rank) + "_Best_Agent")
except Exception as e:
    try:
        print(e)
        model = DDPG.load(Model_SAVE_LOCATION + "/" + str(rank) + "_Agent", env=env)
    except:
        print("CREATING NEW MODEL")
        model = DDPG('MlpPolicy', env, train_freq=(1,"episode"), learning_starts=1000, action_noise = action_noise, device="cpu")



for episode in range(5):
    print(rank, episode)
    env.setSeed(episode)
    obs = env.reset()
    done = False
    while not done:
        action, _states = model.predict(obs)
        obs, rewards, done, info = env.step(action)
        # env.render()

    env.save_training_data()
