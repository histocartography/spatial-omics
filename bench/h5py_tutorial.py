# %%
import h5py
import numpy as np
import pandas as pd
import networkx as nx

# %% creat file
f = h5py.File("mytestfile.hdf5", "w")
# f = h5py.File("mytestfile.hdf5", "r")

# %% create data set
dset = f.create_dataset("mydataset", (100,), dtype='i')
dset2 = f.create_dataset("myotherdataset", (200,), dtype='i')
dset.name
f.name

# %% create groups
f = h5py.File('mydataset.hdf5', 'a')  # open file in append mode
grp = f.create_group("subgroup")

dset2 = grp.create_dataset("another_dataset", (50,), dtype='f')
dset2.name

dset3 = f.create_dataset('subgroup2/dataset_three', (10,), dtype='i')  # specify full path
dset3.name

# %% retrieve data
dataset_three = f['subgroup2/dataset_three']

for name in f: print(name)
f.visit(lambda x: print(x))

d = f['mydataset'][...]

# %% pandas
df = pd.DataFrame({'a': [1,2,3], 'b': [5,6,7]})

df.to_hdf('mytestfile.hdf5', 'spl1/df')

# %% setup h5py data set
import pickle as pk
import os

user = os.path.expanduser('~')
with open(os.path.join(user,'X.pkl'), 'rb') as f:
    X = pk.load(f)
with open(os.path.join(user,'obs.pkl'), 'rb') as f:
    obs = pk.load(f)
with open(os.path.join(user,'spl.pkl'), 'rb') as f:
    spl = pk.load(f)
with open(os.path.join(user,'meta.pkl'), 'rb') as f:
    meta = pk.load(f)
with open(os.path.join(user,'G.pkl'), 'rb') as f:
    G = pk.load(f)
with open(os.path.join(user,'var.pkl'), 'rb') as f:
    var = pk.load(f)

with open(os.path.join(user,'cellmasks.pkl'), 'rb') as f:
    cellmasks = pk.load(f)
with open(os.path.join(user,'stromamask.pkl'), 'rb') as f:
    stromamask = pk.load(f)
with open(os.path.join(user,'images.pkl'), 'rb') as f:
    images = pk.load(f)

# %%
fname = "spatialOmics.hdf5"
h5pyf = h5py.File(fname, "a")

# %%
h5pyf['graph_engine'] = 'networkx'
h5pyf['random_seed'] = 42
h5pyf['pickle_file'] = ""
h5pyf['h5py_file'] = "spatialOmics.hdf5"

# %%
grp = 'X'
for sample in X.keys():
    path = os.path.join(grp, sample)
    X[sample].to_hdf("spatialOmics.hdf5", path)

# %%
grp = 'obs'
for sample in X.keys():
    path = os.path.join(grp, sample)
    obs[sample].to_hdf("spatialOmics.hdf5", path, format="table")

# %%
grp = 'G'
for sample in X.keys():
    path = os.path.join(grp, sample)
    tmp = h5pyf.create_group(path)
    g = nx.to_numpy_array(G[sample])
    tmp.create_dataset(sample, data=g)

# %%
grp = 'var'
for sample in X.keys():
    path = os.path.join(grp, sample)
    var[sample].to_hdf("spatialOmics.hdf5", path, format="table")

# %%
grp = 'spl'
path = os.path.join(grp)
spl.to_hdf("spatialOmics.hdf5", path, format="table")

# %%
grp = 'images'
for sample in X.keys():
    path = os.path.join(grp, sample)
    tmp = h5pyf.create_group(path)
    tmp.create_dataset(sample, data = images[sample])

# %%
grp = 'masks'
for sample in X.keys():
    path = os.path.join(grp, sample, 'cellmasks')
    tmp = h5pyf.create_group(path)
    tmp.create_dataset(sample, data = cellmasks[sample])

# %%
grp = 'masks'
for sample in X.keys():
    path = os.path.join(grp, sample, 'stromamask')
    tmp = h5pyf.create_group(path)
    tmp.create_dataset(sample, data = stromamask[sample])

# %%
h5pyf.close()
f.close()

# %%
import h5py
import os
f = h5py.File(fname, 'r')
