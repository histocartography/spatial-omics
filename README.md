# `SpatialOmics`
The `SpatialOmics` class is designed to accommodate storing and processing spatial omics datasets in a technology-agnostic and memory-efficient way. A `SpatialOmics` instance incorporates multiple attributes that bundle together the multiplexed raw images with the segmentation masks, cell-cell graphs, single-cell values, and sample-, feature- and cell-level annotations, as outlined in the figure below. Since ATHENA works with multiplexed images, memory complexity is a problem. `SpatialOmics` stores data in a HDF5 file and lazily loads the required images on the fly to keep the memory consumption low. The `SpatialOmics` structure is sample-centric, i.e., all samples from a spatial omics experiment are stored separately by heavily using Python dictionaries. 

![overview](img/spatialOmics.png)

Specifically, each `SpatialOmics` instance contains the following attributes:
1. `.images`: A Python dictionary (length: `#samples`) of raw multiplexed images, where each sample is mapped to a [numpy](https://numpy.org/) array of shape: `#features x image_width x image_height`.
2. `.masks`: A nested Python dictionary (length: `#samples`) supporting different types of segmentation masks (e.g., cell and tissue masks), where each sample is mapped to an inner dictionary (length: `#mask_types`), and each value of the inner dictionary is a binary [numpy](https://numpy.org/) array of shape: `#image_width x image_height`.
3. `.G`: A nested Python dictionary (length: `#samples`) supporting different topologies of graphs (e.g., knn, contact or radius graph), where each sample is mapped to an inner dictionary (length: `#graph_types`), and each value of the inner dictionary is a [networkx](https://networkx.org/) graph. 
4. `.X`: A Python dictionary of single-cell measurements (length: `#samples`), where each sample is mapped to a [pandas](https://pandas.pydata.org/) dataframe of shape: `#single_cells x #features`. The values in `.X` can either be uploaded or directly computed from `.images` and `.masks`.
5. `.spl`: A [pandas](https://pandas.pydata.org/) dataframe containing sample-level annotations (e.g., patient clinical data) of shape: `#samples x #annotations`.
6. `.obs`: A Python dictionary (length: `#samples`) containing single-cell-level annotations (e.g., cluster id, cell type, morphological fatures), where each sample is mapped to a [pandas](https://pandas.pydata.org/) dataframe of shape: `#single_cells x #annotations`. 
7. `.var`: A Python dictionary (length: `#samples`) containing feature-level annotations (e.g., name of protein/transcript), where each sample is mapped to a [pandas](https://pandas.pydata.org/) dataframe of shape: `#features x #annotations`. 
8. `.uns`: A Python dictionary containing unstructed data, e.g. various colormaps, experiment properties etc.

## Usage
```python
import tarfile
import tempfile
from skimage import io
import os
import pandas as pd
```

```python
from spatialOmics import SpatialOmics

# create empty instance
so = SpatialOmics()
```

```python
import urllib.request
import tarfile

# url from which we download example images
url = 'https://ndownloader.figshare.com/files/29006556'
filehandle, _ = urllib.request.urlretrieve(url)
```

```python
# extract images from tar archive
fimg = 'BaselTMA_SP41_15.475kx12.665ky_10000x8500_5_20170905_122_166_X15Y4_231_a0_full.tiff'
fmask = 'BaselTMA_SP41_15.475kx12.665ky_10000x8500_5_20170905_122_166_X15Y4_231_a0_full_maks.tiff'
fmeta = 'meta_data.csv'
root = 'spatialOmics-tutorial'

with tempfile.TemporaryDirectory() as tmpdir:
    with tarfile.open(filehandle, 'r:gz') as tar:
        tar.extractall(tmpdir)
        
        img = io.imread(os.path.join(tmpdir, root, fimg))
        mask = io.imread(os.path.join(tmpdir, root, fmask))
        meta = pd.read_csv(os.path.join(tmpdir, root, fmeta)).set_index('core')
        
        # set sample data of spatialOmics
        so.spl = meta[[fimg in i for i in meta.filename_fullstack]]
        
        # add high-dimensional tiff image
        so.add_image(so.spl.index[0], os.path.join(tmpdir, root, fimg), to_store=False)
        
        # add segmentation mask
        so.add_mask(so.spl.index[0], 'cellmasks', os.path.join(tmpdir, root, fmask), to_store=False)
```

## Installation
```{bash}
pip install "git+https://github.com/AI4SCR/spatial-omics.git@master"
```