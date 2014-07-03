import h5py
f=h5py.File("sirexp.h5", "r")
ds=f['/trajectory/dset1']
for i in range(50):
    print(ds[i])
