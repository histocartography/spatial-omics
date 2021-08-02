import pandas as pd
import networkx as nx
from skimage import io

import os
import copy

import h5py


# %%
class SpatialOmics:

    def __init__(self):

        self.graph_engine = 'networkx'
        self.random_seed = 42  # for reproducibility
        self.pickle_file = ''  # backend store
        self.h5py_file = 'spatialOmics.h5py'  # backend store

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

    def add_image(self, spl, file, in_memory=True, to_store=True):
        """Add the image for a given sample"""
        im = io.imread(file)

        if to_store:
            path = f'images/{spl}'
            with h5py.File(self.h5py_file, 'a') as f:
                if path in f:
                    del f[path]
                f.create_dataset(path, data=im)

        if in_memory:
            if spl not in self.images:
                self.images[spl] = {}
            self.images[spl] = im

    def get_image(self, spl):
        """Get the image of a given sample"""
        if spl in self.images:
            return self.images[spl]
        else:
            with h5py.File(self.h5py_file, 'r') as f:
                path = f'images/{spl}'
                if path in f:
                    return f[path][:]
                else:
                    raise KeyError(f'no images exists for {spl}.')

    def add_mask(self, spl, mask, file, in_memory=True, to_store=True):
        """Add a mask for a given sample"""
        im = io.imread(file)

        if to_store:
            path = f'masks/{spl}/{mask}'
            with h5py.File(self.h5py_file, 'a') as f:
                if path in f:
                    del f[path]
                f.create_dataset(path, data=im)

        if in_memory:
            if spl not in self.masks:
                self.masks[spl] = {}
            self.masks[spl][mask] = im

    def get_mask(self, spl, mask):
        """Get a particular mask of a given sample"""
        if spl in self.masks and mask in self.masks[spl]:
            return self.masks[spl][mask]
        else:
            with h5py.File(self.h5py_file, 'r') as f:
                path = f'masks/{spl}/{mask}'
                if path in f:
                    return f[path][...]
                else:
                    raise KeyError(f'no {mask} mask exists for {spl}.')

    def __str__(self):
        s = f"""
        SpatialOmics object with
        {len(self.spl)} samples
        """
        return s

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.spl)

    def to_h5py(self, file: str = None) -> None:
        """

        Args:
            file: file to write to, defaults to self.h5py_file

        Returns:

        """
        if file is None:
            file = self.h5py_file

        with h5py.File(file, mode='w') as f:
            f.create_dataset('h5py_file', data=self.h5py_file)

            # obsm
            for spl in self.obsm:
                f.create_dataset(f'obsm/{spl}', data=self.obsm[spl])

            # images
            for spl in self.images:
                img = self.images[spl]
                f.create_dataset(f'images/{spl}', data=img)

            # masks
            for spl in self.masks:
                for key in self.masks[spl].keys():
                    msk = self.masks[spl][key]
                    f.create_dataset(f'masks/{spl}/{key}', data=msk)

            # uns
            # TODO: currently we do not support storing uns to h5py due to datatype restrictions

        # we need to write the dataframes outside the context manager because the file is locked
        # spl
        self.spl.to_hdf(file, 'spl', format="table")

        # var
        for spl in self.var:
            # use pandas function
            self.var[spl].to_hdf(file, f'var/{spl}', format="table")

        # X
        for spl in self.X:
            # use pandas function
            self.X[spl].to_hdf(file, f'X/{spl}', format="table")

        # G
        for spl in self.G:
            for key in self.G[spl]:
                g = self.G[spl][key]
                df = nx.to_pandas_edgelist(g)

                # use pandas function
                df.to_hdf(file, f'G/{spl}/{key}', format="table")

        # obs
        for spl in self.obs:
            # use pandas function
            self.obs[spl].to_hdf(file, f'obs/{spl}', format="table")

        print(f'File `{os.path.basename(file)}` saved to {os.path.abspath(file)}')
        print(f'File size: {os.path.getsize(file) / (1024 * 1024):.2f} MB')

    @classmethod
    def from_h5py(cls, file=None):
        """

        Args:
            file: h5py file from which to reconstruct SpatialOmics instance
            include_images: Whether to load images into memory
            include_mask: Whether to load masks into memory

        Returns:
            SpatialOmics instance

        """

        so = SpatialOmics()

        with h5py.File(file, 'r') as f:
            so.h5py_file = str(f['h5py_file'][...])

            # obs
            if 'obs' in f:
                for spl in f['obs'].keys():
                    so.obs[spl] = pd.read_hdf(file, f'obs/{spl}')

            # obsm
            if 'obsm' in f:
                for spl in f['obsm'].keys():
                    so.obsm[spl] = f[f'obsm/{spl}'][...]

            # spl
            so.spl = pd.read_hdf(file, 'spl')

            # var
            if 'var' in f:
                for spl in f['var'].keys():
                    so.var[spl] = pd.read_hdf(file, f'var/{spl}')

            # X
            if 'X' in f:
                for spl in f['X'].keys():
                    so.X[spl] = pd.read_hdf(file, f'X/{spl}')

            # G
            if 'G' in f:
                for spl in f['G'].keys():
                    for key in f[f'G/{spl}'].keys():
                        if key not in so.G:
                            so.G[spl] = {}
                        so.G[spl][key] = nx.from_pandas_edgelist(pd.read_hdf(file, f'G/{spl}/{key}'))

            # images
            if 'images' in f:
                for spl in f['images'].keys():
                    if spl not in so.images:
                        so.images[spl] = {}
                    so.images[spl] = f[f'images/{spl}'][...]

            # masks
            if 'masks' in f:
                for spl in f['masks'].keys():
                    for key in f[f'masks/{spl}']:
                        if spl not in so.masks:
                            so.masks[spl] = {}
                        so.masks[spl][key] = f[f'masks/{spl}/{key}'][...]

        return so

    def to_pickle(self, file: str = None) -> None:
        """Save spatialOmics instance to pickle.

        Args:
            file: file to which instance is saved.

        Returns:

        """
        raise NotImplementedError('This version does not yet support saving pickled instances')

    @classmethod
    def from_pickle(cls, file: str = None) -> None:
        """Load spatialOmics instance from pickled file.

        Args:
            file: file to un-pickle

        Returns:
            spatialOmics instance

        """
        raise NotImplementedError('This version does not yet support reading pickled instances')

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
