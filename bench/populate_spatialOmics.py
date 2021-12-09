
import os
from skimage import io
from smmhm.IMCData.imc_data import IMCData
ad = IMCData.from_pickle('ad_basel_knn_out.pkl','basel')
ad.set_root(os.path.expanduser('~/Documents/thesis/data/SingleCellPathologyLandscapeOfBreastCancer/basel/'))

from spatialOmics.spatialOmics import SpatialOmics
so = SpatialOmics()

# %%
cores = ad.cores[:5]

so.spl = ad.meta.loc[cores]
so.uns['cmaps'] = ad.uns['cmaps']
so.uns['cmap_labels'] = ad.uns['cmap_labels']
for c in cores:
    so.X[c] = ad.X[c]
    so.obs[c] = ad.obs[c]
    so.var[c] = ad.var[c]
    so.G[c] = ad.G[c]
    so.images[c] = io.imread(ad.meta.loc[c].file_fullstack + '.tiff')
    so.masks[c] = {'cellmask':io.imread(ad.meta.loc[c].file_cellmask),
                   'stromamask':io.imread(ad.meta.loc[c].file_tumor_stroma_mask)}

import pickle
with open('spatialOmics.pkl', 'wb') as f:
    pickle.dump(so,f)

# %%
import h5py
# import numpy as np
import pandas as pd

fname = "/Users/art/Documents/spatial-omics/spatialOmics.hdf5"
so = SpatialOmics.from_h5py(fname)
