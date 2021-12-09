import numpy as np
from spatialOmics import SpatialOmics
from tifffile import imsave
import os
import shutil

path = '/Users/art/Box/thesis/data/SingleCellPathologyLandscapeOfBreastCancer/basel/tiff_stacks'
file = 'BaselTMA_SP41_15.475kx12.665ky_10000x8500_5_20170905_90_88_X11Y5_242_a0_full.tiff'

# %%
def image_file(tmpdir):
    img = np.zeros((5,20,20), 'uint16')
    fn = os.path.join(tmpdir, 'image.tiff')
    imsave(fn, img)
    return fn

def mask_file(tmpdir):
    img = np.zeros((20,20))
    fn = os.path.join(tmpdir, 'mask.tiff')
    imsave(fn, img)
    return fn

# %%

tmpdir=os.path.expanduser('~/tmp')
os.mkdir(tmpdir)
img = image_file(tmpdir)
msk = mask_file(tmpdir)

so = SpatialOmics()
so.add_image('spl1', img)
so.add_image('spl2', img)
so.add_mask('spl1', 'msk1', msk)
so.add_mask('spl1', 'msk2', msk)
so.add_mask('spl3', 'msk4', msk)
so.to_h5py()

shutil.rmtree(tmpdir)
tmp = SpatialOmics.from_h5py(so.h5py_file)
# tmp = so.from_h5py(so.h5py_file)
#
# import h5py
# f = h5py.File(so.h5py_file, 'r')

