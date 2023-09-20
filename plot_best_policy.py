import numpy as np
import matplotlib.pyplot as plt
import os
import math
plt.rc('font', weight='bold')

cytokine_labels = ["IL8", "MIP1b", "TNF", "HGF", "PAF", "IL1", "PDGF", "IL13", "TGFb", "VEGF", "GCSF", "IL10", "IL4", "IL17", "Myf5", "IL6", "IL1ra", "sIL1r", "sTNFr", "endotoxin", "IFNg", "IL12", "MRF4", "Reactive Oxygen Species", "IGF", "DAMP", "MCP1", "FNE"]

actions_scaling_factor = 10
top_n_runs = 50 #saves top number of runs here for faster loading in the future for replotting
top_n_to_plot = 5
# actions_scaling_factor = 10
action_mags = ["Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5"]
# data_dirs = ["SenseM1_95/", "SenseM1_96/", "SenseM1_97/", "SenseM1_98/", "SenseM1_99/", "SenseM1_100/", "SenseM1_101/"]
data_dirs = ["everything_else_local/Interval_4/"]

# actions_to_keep =  [*range(len(cytokine_labels))]
# actions_to_keep = [3] #{HGF}
actions_to_keep = [2,8,11,23] #{TNF, TGFb, IL10, cytotox}
# actions_to_keep = [2,3,8,11,23] #{TNF, HGF, TGFb, IL10, cytotox}
# actions_to_keep = [2,3,8,11] #{TNF, HGF, TGFb, IL10}
# actions_to_keep = [2, 11, 12]   #{TNF IL10, IL4}
cytos_to_keep = [cytokine_labels[i] for i in actions_to_keep]
print(cytos_to_keep)

best_run_ind = np.zeros(shape=(top_n_to_plot))
best_muscle= np.zeros(shape=(top_n_to_plot))
best_actions = np.zeros(shape=(28,1398,top_n_to_plot))
best_cytos = np.zeros(shape=(28,1398,top_n_to_plot))
best_muscle_data = np.zeros(shape=(1398,top_n_to_plot))
best_col1_data = np.zeros(shape=(1398,top_n_to_plot))
best_col3_data = np.zeros(shape=(1398,top_n_to_plot))
best_necrosis_data = np.zeros(shape=(1398,top_n_to_plot))
best_data_index = 0

for dir in data_dirs:
    data_directory = "/media/dale/T7/WEABM_DRL_Experiments/" + dir
    print(data_directory)
    try:
        best_cytos_temp = np.load(data_directory + "top_" + str(top_n_runs) + "_cyto_data.npy")
        best_actions_temp = np.load(data_directory + "top_" + str(top_n_runs) + "_action_data.npy")
        best_muscle_data_temp = np.load(data_directory + "top_" + str(top_n_runs) + "_muscle_data.npy")
        best_col3_data_temp = np.load(data_directory + "top_" + str(top_n_runs) + "_col3_data.npy")
        best_necrosis_data_temp = np.load(data_directory + "top_" + str(top_n_runs) + "_necrosis_data.npy")

    except:
        # data = np.genfromtxt("/media/dale/T7/test/training_data.csv", delimiter=",")
        data = None
        file_names = sorted(os.listdir(data_directory))

        for file in file_names:
            if file.endswith("training_data.csv"):
                if data is not None:
                    temp_data = np.genfromtxt(data_directory + file, delimiter=",")
                    data = np.concatenate((data,temp_data), axis=0)
                else:
                    data = np.genfromtxt(data_directory+file, delimiter=",")

        print(data.shape)

        best_run_ind_temp = np.zeros(shape=(top_n_runs))
        best_muscle_temp= np.zeros(shape=(top_n_runs))
        best_actions_temp = np.zeros(shape=(28,1398,top_n_runs))
        best_cytos_temp = np.zeros(shape=(28,1398,top_n_runs))
        best_muscle_data_temp = np.zeros(shape=(1398,top_n_runs))
        best_col1_data_temp = np.zeros(shape=(1398,top_n_runs))
        best_col3_data_temp = np.zeros(shape=(1398,top_n_runs))
        best_necrosis_data_temp = np.zeros(shape=(1398,top_n_runs))
        for i in range(int(data.shape[0]/60)):
            current_run = data[i*60:(i+1)*60,1:]
            muscle = current_run[0,-1]
            # if current_run[-27,10] != 0:
            #     continue
            if np.any(current_run[:,:1400] > 10000000) or np.any(current_run[:,:1400] < -10000000):
                continue
            duplicate_run_flag = False
            for j in range(top_n_runs):
                if np.array_equal(current_run[0,:250], best_muscle_data_temp[:250,j]):
                    duplicate_run_flag = True
            if duplicate_run_flag:
                continue
            for j in range(top_n_runs):
                if muscle >= best_muscle_temp[j]:
                    # shift all the best ones down after the current best
                    if j >0:
                        for k in reversed(range(j,top_n_runs)):
                            best_muscle_temp[k] = best_muscle_temp[k-1]
                            best_run_ind_temp[k] = best_run_ind_temp[k-1]
                            best_cytos_temp[:,:,k] = best_cytos_temp[:,:,k-1]
                            best_actions_temp[:,:,k] = best_actions_temp[:,:,k-1]
                            best_muscle_data_temp[:,k] = best_muscle_data_temp[:,k-1]
                            best_col3_data_temp[:,k] = best_col3_data_temp[:,k-1]
                            best_necrosis_data_temp[:,k] = best_necrosis_data_temp[:,k-1]
                    else:   # if the best is in the first slot, shift everything down by one
                        for k in reversed(range(1,top_n_runs)):
                            best_muscle_temp[k] = best_muscle_temp[k-1]
                            best_run_ind_temp[k] = best_run_ind_temp[k-1]
                            best_cytos_temp[:,:,k] = best_cytos_temp[:,:,k-1]
                            best_actions_temp[:,:,k] = best_actions_temp[:,:,k-1]
                            best_muscle_data_temp[:,k] = best_muscle_data_temp[:,k-1]
                            best_necrosis_data_temp[:,k] = best_necrosis_data_temp[:,k-1]


                    best_muscle_temp[j] = muscle
                    best_run_ind_temp[j] = i
                    best_cytos_temp[:,:,j] = current_run[4:-28,:-1]
                    # print(best_cytos.shape)
                    best_actions_temp[:,:,j] = current_run[-28:,:-1]
                    best_muscle_data_temp[:,j] = current_run[0,:-1]
                    best_col3_data_temp[:,j] = current_run[1,:-1]
                    best_necrosis_data_temp[:,j] = current_run[2,:-1]
                    break
        print(best_muscle_temp, best_run_ind_temp)

        np.save(data_directory + "top_" + str(top_n_runs) + "_cyto_data", best_cytos_temp)
        np.save(data_directory + "top_" + str(top_n_runs) + "_action_data", best_actions_temp)
        np.save(data_directory + "top_" + str(top_n_runs) + "_muscle_data", best_muscle_data_temp)
        np.save(data_directory + "top_" + str(top_n_runs) + "_col3_data", best_col3_data_temp)
        np.save(data_directory + "top_" + str(top_n_runs) + "_necrosis_data", best_necrosis_data_temp)

    best_cytos[:,:,:] = best_cytos_temp[:,:,:top_n_to_plot]
    best_actions[:,:,:] = best_actions_temp[:,:,:top_n_to_plot]
    best_muscle_data[:,:] = best_muscle_data_temp[:,:top_n_to_plot]
    best_col3_data[:,:] = best_col3_data_temp[:,:top_n_to_plot]
    best_necrosis_data[:,:] = best_necrosis_data_temp[:,:top_n_to_plot]
    best_data_index += 1

best_actions = best_actions

groups = np.any(best_muscle_data<0, axis=0)

best_muscle_data = best_muscle_data-9900

plt.figure(figsize=(8,10))
plt.subplot(3,1,1)
# plt.plot(best_muscle_data[:,groups], c="r")
plt.plot(best_muscle_data[:,~groups])
# plt.plot(baseline_muscle, c="k", linewidth=2.5)
plt.title("Muscle over time", fontweight="bold")
plt.subplot(3,1,2)
# plt.plot(best_col3_data[:,groups], c="r")
plt.plot(best_col3_data[:,~groups])
# plt.plot(baseline_col3, c="k", linewidth=2.5)
plt.legend(labels=action_mags, loc="center left", bbox_to_anchor=(1,0.5))
plt.title("Col3 over time", fontweight="bold")
plt.subplot(3,1,3)
plt.title("Total Closure %", fontweight="bold")
closure = 100*(best_muscle_data + best_col3_data)/(30*30*200)
# plt.plot(baseline_closure, c="k", linewidth=2.5)
plt.plot(closure)
# # plt.plot(best_necrosis_data[:,groups], c="r")
# plt.plot(best_necrosis_data[:,~groups])
# plt.title("Necrosis over time")

plt.tight_layout()

plt.figure(figsize=(10,20))
plt.suptitle("Cytokine Values at Surface (30x30 grid)")
plt.subplot(math.ceil(len(cytos_to_keep)/2),2,1)
ind = 0
for cyto in cytos_to_keep:
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    # set the spacing between subplots
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.4)

    # plt.plot(best_cytos[ind,:,:])
    # print(best_cytos.shape)

    # r = best_cytos[ind,:,groups]
    # print(r.shape)
    # for i in range(r.shape[0]):
        # plt.plot(r[i,:], c="r")
    r = best_cytos[actions_to_keep[ind],:,~groups]
    for i in range(r.shape[0]):
        plt.plot(r[i,:])

    plt.title(cyto, fontweight="bold")
    ind+=1
    if ind < len(cytos_to_keep):
        plt.subplot(math.ceil(len(cytos_to_keep)/2),2,ind+1)
plt.legend(labels=action_mags, loc="center left", bbox_to_anchor=(1,0.5))



ind = 0

cytokines_augment = []
cytokines_diminish = []
for cyto in cytokine_labels:
    if np.mean(best_actions[ind,:,0]) > 0:
        cytokines_augment.append(cyto)
    else:
        cytokines_diminish.append(cyto)
    ind +=1

print("Augmented cytokines: ")
print(cytokines_augment)

print("Diminished cytokines: ")
print(cytokines_diminish)

plt.figure(figsize=(8,6))
plt.title("Policy Ranking Heatmap")
mean_actions = np.mean(best_actions, axis=1).T
mean_actions = mean_actions[:,actions_to_keep]
print(mean_actions.shape)
plt.imshow(mean_actions, cmap="seismic")
plt.ylabel("Policy Ranking")
plt.xticks(np.arange(len(cytos_to_keep)), labels=cytos_to_keep, rotation=90, fontweight="bold")
yticks_pos = np.concatenate((np.array([0]), np.arange(9,mean_actions.shape[0], 10)))
yticks_label = np.concatenate((np.array([1]), np.arange(10,mean_actions.shape[0]+1, 10)))
plt.yticks(yticks_pos, labels = yticks_label, fontweight="bold")
# plt.legend(loc="right")


plt.figure(figsize=(10,20))
plt.suptitle("Cytokine Actions", fontweight="bold")
if len(cytos_to_keep) < 2:
    num_cols = 1
else:
    num_cols = 2
num_cols
plt.subplot(math.ceil(len(cytos_to_keep)/2),num_cols,1)
ind = 0
for cyto in cytos_to_keep:
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    # set the spacing between subplots
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.4)

    plt.ylim(-1*actions_scaling_factor - .5, actions_scaling_factor+.5)

    # plt.ylim(-1*actions_scaling_factor - .5, actions_scaling_factor+.5)
    # plt.plot(best_actions[ind,:,groups], c="r")
    # plt.plot(best_actions[ind,:,~groups], c="b")
    # r = best_actions[ind,:,groups]
    # for i in range(r.shape[0]):
        # plt.plot(r[i,:], c="r")
    r = best_actions[actions_to_keep[ind],:,~groups]
    for i in range(r.shape[0]):
        plt.plot(r[i,:])

    plt.title(cyto, fontweight="bold")
    ind+=1
    if ind < len(cytos_to_keep):
        plt.subplot(math.ceil(len(cytos_to_keep)/2),num_cols,ind+1)
plt.legend(labels=action_mags, loc = "center left", bbox_to_anchor=(1,0.5))

# 100% of wound closed = 30*30*200
plt.figure(figsize=(8,10))
plt.subplot(3,1,1)
plt.suptitle("Wound Closure", fontweight="bold")
plt.title("Total Closure %", fontweight="bold")
closure = 100*(best_muscle_data + best_col3_data)/(30*30*200)
plt.plot(closure)

plt.subplot(3,1,2)
plt.title("Muscle %", fontweight="bold")
plt.plot(100*(best_muscle_data)/(30*30*200))
plt.legend(labels=action_mags, loc="center left", bbox_to_anchor=(1,0.5))

plt.subplot(3,1,3)
plt.title("Collagen3 %", fontweight="bold")
plt.plot(100*(best_col3_data)/(30*30*200))
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.show()
