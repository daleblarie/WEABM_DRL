from stable_baselines3.common.callbacks import BaseCallback
import h5py
import os
import shutil
import numpy as np
from mpi4py import MPI

class CustomCallback(BaseCallback):
    """
    A custom callback that derives from ``BaseCallback``.

    :param verbose: (int) Verbosity level 0: not output 1: info 2: debug
    """
    def __init__(self, verbose=0, SAVE_FILE_NAME="training_data", rendering = False):
        super(CustomCallback, self).__init__(verbose)
        self.created_spatial_dataset = False
        self.save_path = SAVE_FILE_NAME
        self.render = rendering
        self.episode_num = 0
        self.cum_reward = 0
        self.best_reward = -10000
        print("callback Created")
        # Those variables will be accessible in the callback
        # (they are defined in the base class)
        # The RL model
        # self.model = None  # type: BaseRLModel
        # An alias for self.model.get_env(), the environment used for training
        # self.training_env = None  # type: Union[gym.Env, VecEnv, None]
        # Number of time the callback was called
        # self.n_calls = 0  # type: int
        # self.num_timesteps = 0  # type: int
        # local and global variables
        # self.locals = None  # type: Dict[str, Any]
        # self.globals = None  # type: Dict[str, Any]
        # The logger object, used to report things in the terminal
        # self.logger = None  # type: logger.Logger
        # # Sometimes, for event callback, it is useful
        # # to have access to the parent object
        # self.parent = None  # type: Optional[BaseCallback]

    def _on_training_start(self) -> None:
        """
        This method is called before the first rollout starts.
        """
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank==0:
            # try:
            #     shutil.rmtree(self.save_path)
            #     os.remove(self.save_path+"/training_data.csv")
            #     os.remove(self.save_path+"/spatial_data.npy")
            # except Exception as e:
            #     print(e)
            try:
                os.mkdir(self.save_path)
            except Exception as e:
                print(e)
        print("starting training")

        pass

    def _on_rollout_start(self) -> None:
        """
        A rollout is the collection of environment interaction
        using the current policy.
        This event is triggered before collecting new samples.
        """
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        self.model.save(self.save_path + "/" + str(rank) + "_Agent")

    def _on_step(self) -> bool:
        """
        This method will be called by the model after each call to `env.step()`.

        For child callback (of an `EventCallback`), this will be called
        when the event is triggered.

        :return: (bool) If the callback returns False, training is aborted early.
        """
        self.cum_reward += self.training_env.get_attr("reward")[0]
        if self.training_env.get_attr("done")[0]:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            # with open(self.save_path+"/" + str(rank) + "_training_data.csv", "ab") as f:
            #     saveArray = self.training_env.get_attr("viable_muscle_history")[0]
            #     saveArray = np.vstack((saveArray, self.training_env.get_attr("collagen3_history")[0]))
            #     saveArray = np.vstack((saveArray, self.training_env.get_attr("necrosis_history")[0]))
            #     saveArray = np.vstack((saveArray, self.training_env.get_attr("collagen1_history")[0]))
            #     saveArray = np.vstack((saveArray, self.training_env.get_attr("norm_hist")[0].T))
            #     saveArray = np.vstack((saveArray, self.training_env.get_attr("action_history")[0].T))
            #     np.savetxt(f, saveArray, delimiter=",")

            with open(self.save_path+"/" + str(rank)+"_simulation_output.csv", "ab") as f:
                np.savetxt(f, self.training_env.get_attr("allSignalsReturn_output")[0], delimiter=",")

            try:
                save_data = np.load(self.save_path+"/" + str(rank) + "_muscle_spatial_data.npy", allow_pickle=True)
                save_data = np.concatenate((save_data, np.expand_dims(self.training_env.get_attr("muscle_spatial")[0], axis=0)), axis=0)
                np.save(self.save_path+"/" + str(rank) + "_muscle_spatial_data.npy", save_data)
            except FileNotFoundError:
                np.save(self.save_path+"/" + str(rank) + "_muscle_spatial_data.npy", np.expand_dims(self.training_env.get_attr("muscle_spatial")[0], axis=0))

            try:
                save_data = np.load(self.save_path+"/" + str(rank) + "_collagen1_spatial_data.npy", allow_pickle=True)
                save_data = np.concatenate((save_data, np.expand_dims(self.training_env.get_attr("collagen1_spatial")[0], axis=0)), axis=0)
                np.save(self.save_path+"/" + str(rank) + "_collagen1_spatial_data.npy", save_data)
            except FileNotFoundError:
                np.save(self.save_path+"/" + str(rank) + "_collagen1_spatial_data.npy", np.expand_dims(self.training_env.get_attr("collagen1_spatial")[0], axis=0))

            try:
                save_data = np.load(self.save_path+"/" + str(rank) + "_collagen3_spatial_data.npy", allow_pickle=True)
                save_data = np.concatenate((save_data, np.expand_dims(self.training_env.get_attr("collagen3_spatial")[0], axis=0)), axis=0)
                np.save(self.save_path+"/" + str(rank) + "_collagen3_spatial_data.npy", save_data)
            except FileNotFoundError:
                np.save(self.save_path+"/" + str(rank) + "_collagen3_spatial_data.npy", np.expand_dims(self.training_env.get_attr("collagen3_spatial")[0], axis=0))

            try:
                save_data = np.load(self.save_path+"/" + str(rank) + "_necrosis_spatial_data.npy", allow_pickle=True)
                save_data = np.concatenate((save_data, np.expand_dims(self.training_env.get_attr("necrosis_spatial")[0], axis=0)), axis=0)
                np.save(self.save_path+"/" + str(rank) + "_necrosis_spatial_data.npy", save_data)
            except FileNotFoundError:
                np.save(self.save_path+"/" + str(rank) + "_necrosis_spatial_data.npy", np.expand_dims(self.training_env.get_attr("necrosis_spatial")[0], axis=0))

            if self.render:
                print("\033[FCompleted Episode "+ str(self.episode_num))
                print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
            self.episode_num+=1

            print(str(rank) +", " + str(self.episode_num - 1) + ", " + str(np.max(self.training_env.get_attr("viable_muscle_history")[0])))

        if self.render:
            self.training_env.render()
        return True

    def _on_rollout_end(self) -> None:
        """
        This event is triggered before updating the policy.
        """
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if os.path.isfile(self.save_path+"/" + str(rank) + "_reward.csv"):
            with open(self.save_path+"/" + str(rank) + "_reward.csv", "ab") as f:
                np.savetxt(f, [self.cum_reward], delimiter=",")
        else:
            np.savetxt(self.save_path+"/" + str(rank) + "_reward.csv", [self.cum_reward], delimiter=",")
        if self.cum_reward >= self.best_reward:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            self.model.save(self.save_path + "/" + str(rank) + "_Best_Agent")
            self.best_reward = self.cum_reward
        self.cum_reward = 0

    def _on_training_end(self) -> None:
        """
        This event is triggered before exiting the `learn()` method.
        """
        pass
