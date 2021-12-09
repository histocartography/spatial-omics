import pandas as pd
import networkx as nx
from skimage import io
import seaborn as sns

import os
import copy

import h5py


# %%
class SpatialOmics:
    cmaps = {
        'default': sns.color_palette('Reds', as_cmap=True),
        'category': sns.color_palette('Set3', as_cmap=True)
    }

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
        self.G = {}  # graphs
        self.uns = {}  # unstructured container
        self.uns.update({'cmaps':self.cmaps,
                         'cmap_labels': {}})

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

    @staticmethod
    def from_annData(ad: AnnData, img_container: ImageContainer = None,
          sample_id: str = 'sample_id',
          img_layer: str = 'image',
          segmentation_layers: list = ['segmented_watershed']) -> SpatialOmics:
        """Converts a AnnData instance to a SpatialOmics instance.

        Args:
            ad: AnnData object
            img_container: Squidpy Image Container
            sample_id: column name that identifies different libraries in ad.obs

        Returns:
            SpatialOmics
        """

        if sample_id not in ad.obs:
            sample_name = 'spl_0'
            ad.obs[sample_id] = sample_name
            ad.obs[sample_id] = ad.obs[sample_id].astype('category')

        if len(ad.obs[sample_id].unique()) > 1:
            raise ValueError("""more than 1 sample_id present in ad.obs[sample_id].
            Please process each each sample individually.""")
        else:
            sample_name = ad.obs[sample_id][0]

        so = SpatialOmics()
        x = pd.DataFrame(ad.X.A, columns=ad.var.index)

        so.X = {sample_name: x}
        so.obs = {sample_name: ad.obs}
        so.var = {sample_name: ad.var}
        so.spl = pd.DataFrame(index=[sample_name])

        if 'spatial' in ad.obsm:
            coord = ad.obsm['spatial']
            coord = pd.DataFrame(coord, index=so.obs[sample_name].index, columns=['x','y'])
            so.obs[sample_name] = pd.concat((so.obs[sample_name], coord), 1)

        if img_container is not None:
            img = img_container[img_layer]
            so.images = {sample_name: img}

            segmentations = {}
            for i in segmentation_layers:
                if i in img_container:
                    segmentations.update({i:img_container[i]})
            so.masks.update({sample_name:segmentations})

        return so

    def to_annData(self,
        one_adata=True,
        spatial_keys_so=['x', 'y'],
        spatial_key_ad='spatial'):
        """Converts the current SpatialOmics instance into a AnnData instance.
        Does only takes .X, .obs and .var attributes into account.

        Args:
            one_adata: bool whether for each sample a individual AnnData should be created
            spatial_keys_so: tuple column names of spatial coordinates of observations in so.obs[spl]
            spatial_key_ad: str key added to ad.obsm to store the spatial coordinates

        Returns:

        """

        so = self

        if one_adata:
            keys = list(so.obs.keys())
            # we iterate through the keys to ensure that we have the order of the different dicts aligned
            # we could apply pd.concat directly on the dicts

            X = pd.concat([so.X[i] for i in keys])
            obs = pd.concat([so.obs[i] for i in keys])
            obs.index = range(len(obs))
            var = so.var[keys[0]]

            # create AnnData
            ad = AnnData(X=X.values, obs=obs, var=var)

            # import spatial coordinates
            if all([i in obs for i in spatial_keys_so]):
                spatial_coord = ad.obs[spatial_keys_so]
                ad.obs = ad.obs.drop(columns=spatial_keys_so)
                ad.obsm.update({spatial_key_ad: spatial_coord.values})

            return ad

        else:
            ads = []
            for spl in so.obs.keys():

                # create AnnData
                ad = AnnData(X=so.X[spl].values, obs=so.obs[spl], var=so.var[spl])

                # import spatial coordinates
                if all([i in ad.obs for i in spatial_keys_so]):
                    spatial_coord = ad.obs[spatial_keys_so]
                    ad.obs = ad.obs.drop(columns=spatial_keys_so)
                    ad.obsm.update({spatial_key_ad: spatial_coord.values})

                ads.append(ad)
            return ads

