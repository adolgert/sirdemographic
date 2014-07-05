import numpy as np
import matplotlib.pyplot as plt


def plot_single(traj_dset):
    cnt=len(traj_dset)
    # The dataset is an array of tuples, not a 2d array.
    s=np.zeros(cnt, np.int64)
    i=np.zeros(cnt, np.int64)
    r=np.zeros(cnt, np.int64)
    t=np.zeros(cnt, np.double)
    for idx, val in enumerate(traj_dset):
        s[idx], i[idx], r[idx], t[idx]=val

    # Does t sometimes go backwards? By how much?
    tmax=t[0] # largest seen
    timax=0   # largest difference between now and max.
    for idx in range(1,len(t)):
        if t[idx]<tmax:
            timax=max(tmax-t[idx], timax)
        tmax=max(t[idx], tmax)
    print("A value was {0} less than a previous value".format(timax))

    fix, ax=plt.subplots()
    ax.plot(t, s, 'k')
    ax.plot(t, i, 'r')
    plt.show()



def foreach_trajectory(filename, func):
    import h5py
    f=h5py.File(filename, "r")
    trajectories=f['/trajectory']

    for trajectory_name in trajectories:
        dset=trajectories[trajectory_name]
        func(dset)



if __name__ == "__main__":
    foreach_trajectory("sirexp.h5", print)
    foreach_trajectory("sirexp.h5", plot_single)
