import pandas as pd
import numpy as np
import networkx as nx
from skimage import io
import time

import os
import copy

# from tqdm import tqdm
import pickle
import h5py

# %%
class SpatialOmics():

    def __init__(self):

        self.graph_engine = 'networkx'
        self.random_seed = 42  # for reproducibility
        self.pickle_file = ''  # backend store
        self.h5py_file = ''  # backend store

        self.obs = {}  # container for observation level features
        self.obsm = {}  # container for multidimensional observation level features
        self.spl = pd.DataFrame()  # container for sample level features
        self.splm = {}  # container for multidimensional sample level features
        self.var = {}  # container with variable descriptions of X
        self.X = {}  # container for cell level expression data of each spl
        self.uns = {}  # unstructured container
        self.G = {}  # graphs

        # self.obs_keys = None  # list of observations in self.obs.index
        self.spl_keys = None  # list of samples in self.spl.index

        self.images = {}
        self.masks = {}

    def add_image(self, spl, file):
        """Add the image for a given sample"""
        if spl not in self.images:
            self.masks[spl] = {}
        self.images[spl] = io.imread(file)

    def get_image(self, spl):
        """Get the image of a given sample"""
        if spl in self.images:
            return self.images[spl]
        elif True:
            f = self.spl.loc[spl].file_fullstack
            self.add_image(spl, f)
            return self.images[spl]
        else:
            with h5py.File(self.h5py_file) as f:
                path = f'masks/{spl}/'
                if path in f:
                    return f[path][spl][...]
                else:
                    raise KeyError(f'no images exists for {spl}.')

    def add_mask(self, spl, mask, file):
        """Add a mask for a given sample"""
        if spl not in self.masks:
            self.masks[spl] = {}
        self.masks[spl].update({mask:io.imread(file)})

    def get_mask(self, spl, mask):
        """Get a particular mask of a given sample"""
        if spl in self.masks and mask in self.masks[spl]:
            return self.masks[spl][mask]
        elif True:
            f = self.spl.loc[spl].file_cellmask
            self.add_mask(spl, mask, f)
            return self.masks[spl][mask]
            # if spl not in self.mask:
            #     self.masks[spl] = {}
            # self.masks[spl].update({'cellmasks': io.imread(f)})
        else:
            with h5py.File(self.h5py_file) as f:
                path = f'masks/{spl}/{mask}'
                if path in f:
                    return f[path][spl][...]
                else:
                    raise KeyError(f'no {mask} mask exists for {spl}.')
                return self.masks[spl][mask]

    def __str__(self):
        s = """SpatialOmics object
        With this info.
        """
        return s

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.spl)

    def to_h5py(self, file=None):
        if file is None:
            file = self.h5py_file

        print(f'File `{os.path.basename(file)}` save at {os.path.dirname(file)}')
        print(f'File size: {os.path.getsize(file) / (1024 * 1024):.2f} MB')

    @classmethod
    def form_h5py(cls, file=None, include_images=False, include_mask=False):
        """

        Args:
            file: h5py file from which to reconstruct SpatialOmics instance
            include_images: Whether to load images into memory
            include_mask: Whether to load masks into memory

        Returns:
            SpatialOmics instance

        """

        so = SpatialOmics()
        # all class attributes
        attrs = ['graph_engine', 'random_seed', 'pickle_file', 'h5py_file', 'obs', 'spl', 'var', 'X', 'G', 'uns', 'spl_keys', 'images', 'masks']

        # sample-level attributes
        spl_attrs = ['G', 'X', 'images', 'masks', 'obs', 'var']
        pd_attrs = ['X', 'obs', 'var']  # panda dataframes
        nx_attrs = ['G']  # networkx graphs
        uns_attrs = ['masks', 'uns']
        with h5py.File(file, 'r') as f:
            for attr in attrs:
                if attr not in f:
                    continue
                if attr == 'masks' and include_mask is False:
                    continue
                if attr == 'images' and include_images is False:
                    continue

                if attr in spl_attrs:
                    spls = f[attr].keys()
                    tmp = {}
                    for spl in spls:
                        if attr in pd_attrs:
                            tmp.update({spl:pd.read_hdf(file, f'{attr}/{spl}/')})
                        elif attr in nx_attrs:
                            # TODO: For the reconsturction of the graphs we need to save the observation ids as well otherwise we start at 0
                            m = f[attr][spl][spl][...]
                            tmp.update({spl:nx.from_numpy_array(m)})
                        elif attr in uns_attrs:
                            key_item = {}
                            for key in f[attr][spl].keys():
                                m = f[attr][spl][key][spl][...]
                                key_item.update({key: m})
                            tmp.update({spl: key_item})
                        else:
                            m = f[attr][spl][spl][...]
                            tmp.update({spl:m})
                    setattr(so, attr, tmp)
                elif attr == 'spl':
                    tmp = pd.read_hdf(file, f'{attr}/')
                    setattr(so, attr, tmp)
                else:
                    m = f[attr][...]
                    if attr == 'random_seed':
                        setattr(so, attr, int(m))
                    else:
                        setattr(so, attr, m)
        return so

    def to_pickle(self, file=None):
        if file is None:
            file = self.pickle_file

        print(f'File `{os.path.basename(file)}` save at {os.path.dirname(file)}')
        print(f'File size: {os.path.getsize(file) / (1024 * 1024):.2f} MB')

    @classmethod
    def from_pickle(cls, file=None):
        with open(file, 'rb') as f:
            so = pickle.load(f)
        return so

    def copy(self):
        """copy IMCData without copying graphs, masks and tiffstacks"""
        c = copy.copy(self)
        c.obs = copy.deepcopy(self.obs)
        c.var = copy.deepcopy(self.var)
        c.X = copy.deepcopy(self.X)
        c.uns = copy.deepcopy(self.uns)
        return c

    def deepcopy(self):
        return copy.deepcopy(self)