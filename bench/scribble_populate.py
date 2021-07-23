from spatialOmics import SpatialOmics
import os

path = '/Users/art/Box/thesis/data/SingleCellPathologyLandscapeOfBreastCancer/basel/tiff_stacks'
file = 'BaselTMA_SP41_15.475kx12.665ky_10000x8500_5_20170905_90_88_X11Y5_242_a0_full.tiff'

# %%

so = SpatialOmics()
so.add_image('spl1', os.path.join(path,file))