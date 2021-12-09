# %%
import logging

import matplotlib.colors
import scanpy as sc
import squidpy as sq
from squidpy.im import ImageContainer
from spatialOmics import SpatialOmics
import spatialHeterogeneity as sh
from anndata import AnnData

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, hex2color
import seaborn as sns
import numpy as np
import pandas as pd

# %%

def ad2so(ad: AnnData, img_container: ImageContainer = None,
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

def so2ad(so: SpatialOmics,
          one_adata=True,
          spatial_keys_so=['x', 'y'],
          spatial_key_ad='spatial'):
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

# %% load the pre-processed dataset

# imc = sq.datasets.imc()
# fish = sq.datasets.seqfish()
# slide = sq.datasets.slideseqv2()
# mibi = sq.datasets.mibitof()
# ad4i = sq.datasets.four_i()

# # visium with segmentation masks
# visium = sq.datasets.visium_fluo_adata_crop()
# visium_img = sq.datasets.visium_fluo_image_crop()
# sq.im.process(img=visium_img,layer="image",method="smooth",)
# sq.im.segment(img=visium_img, layer="image_smooth", method="watershed", channel=0, chunks=1000)

# # hne with segmentation masks
# hne_img = sq.datasets.visium_hne_image()
# hne = sq.datasets.visium_hne_adata()
# sq.im.process(hne_img, layer="image", method="smooth", sigma=4)
# sq.im.segment(img=hne_img, layer="image_smooth", method="watershed", thresh=90, geq=False)

# so = sh.dataset.imc()
# spl = list(so.X.keys())[0]
# for i in so.obs.keys():
#     sh.pp.extract_centroids(so, i)

# %% plot images and segmentation masks

fig,axs = plt.subplots(2,3)
for img,idx in zip([hne_img, visium_img], range(2)):
    img.show(layer='image', interpolation='none', ax=axs[idx,0])
    img.show(layer='segmented_watershed', interpolation='none', cmap='gray', ax=axs[idx,1])
    axs[idx,2].imshow(img['segmented_watershed'].squeeze()>0, interpolation='none', cmap='gray')
    axs[idx,2].set_axis_off()
fig.tight_layout()
fig.show()

cmap = sns.color_palette('Set3')
cmap = cmap * np.ceil(len(np.unique(crop['segmented_watershed'])) / len(cmap)).astype(int)
cmap.insert(0, (0,0,0))
cmap = ListedColormap(cmap)

# %% convert imc
imc = sq.datasets.imc()
so = ad2so(imc)
sc.pl.spatial(imc, color="cell type", spot_size=10)

# we have some overhead here as we need to convert to numeric types for our framework
spl = list(so.obs.keys())[0]
so.obs[spl]['cell_type_id'] = so.obs[spl].groupby('cell type').ngroup().astype('category')

# generate colormap
labs = so.obs[spl].groupby(['cell_type_id']).head(1)[['cell_type_id', 'cell type']].set_index('cell_type_id').to_dict()
cmap = ListedColormap([hex2color(i) for i in imc.uns['cell type_colors']])
so.uns['cmaps'].update({'cell_type_id':cmap})
so.uns['cmap_labels'].update({'cell_type_id': labs['cell type']})

# graph building, metrics, plot
sh.graph.build_graph(so, spl, mask_key=None)
sh.pl.spatial(so, spl, attr='cell_type_id', edges=True)
sh.metrics.shannon(so, spl, attr='cell_type_id')
sh.pl.spatial(so, spl, attr='shannon_cell_type_id_knn')

# %% convert fish
fish = sq.datasets.seqfish()
so = ad2so(fish)
spl = list(so.obs.keys())[0]

# plot AnnData
sc.pl.spatial(fish, color="celltype_mapped_refined", spot_size=0.03)

# we have some overhead here as we need to convert to numeric types for our framework
col_uns_name = 'celltype_mapped_refined'
so.obs[spl]['cell_type_id'] = so.obs[spl].groupby(col_uns_name).ngroup().astype('category')

# generate colormap
from matplotlib import cm
labs = so.obs[spl].groupby(['cell_type_id']).head(1)[['cell_type_id', col_uns_name]].set_index('cell_type_id').to_dict()
cmap = ListedColormap([hex2color(i) for i in fish.uns[col_uns_name+'_colors']])
so.uns['cmaps'].update({'cell_type_id':cmap})
so.uns['cmap_labels'].update({'cell_type_id': labs[col_uns_name]})
so.uns['cmaps'].update({'default': cm.plasma})

# graph building, metrics, plot
sh.graph.build_graph(so, spl, mask_key=None)
sh.pl.spatial(so, spl, attr='cell_type_id', edges=True)
sh.metrics.shannon(so, spl, attr='cell_type_id')
sh.pl.spatial(so, spl, attr='shannon_cell_type_id_knn', background_color='black', node_size= 1)

# neighborhood enrichment squidpy
sq.gr.spatial_neighbors(fish, coord_type="generic")
sq.gr.nhood_enrichment(fish, cluster_key="celltype_mapped_refined")
sq.pl.nhood_enrichment(fish, cluster_key="celltype_mapped_refined", method="ward");plt.tight_layout();plt.show()

# neighborhood enrichment ATHENA, observation
sh.neigh.interactions(so, spl, 'cell_type_id', mode='proportion')
sh.pl.interactions(so, spl, 'cell_type_id', mode='proportion', prediction_type='observation')
# neighborhood enrichment ATHENA, diff
sh.neigh.interactions(so, spl, 'cell_type_id', mode='proportion', prediction_type='diff')
fig, ax = plt.subplots(figsize=(8,6))
sh.pl.interactions(so, spl, 'cell_type_id', mode='proportion', prediction_type='diff', ax=ax)
ylab = ax.get_ymajorticklabels()
newlab = []
for lab in ylab:
    n = so.uns['cmap_labels']['cell_type_id'][int(lab.get_text())] + f',{lab.get_text()}'
    newlab.append(n)
ax.set_yticklabels(newlab)
fig.tight_layout()
fig.show()

# %% convert visium with segmentation masks
visium = sq.datasets.visium_fluo_adata_crop()
visium_img = sq.datasets.visium_fluo_image_crop()
sq.im.process(img=visium_img,layer="image",method="smooth")
sq.im.segment(img=visium_img, layer="image_smooth", method="watershed", channel=0, chunks=1000)

# convert
so = ad2so(visium, img_container=visium_img)

# plot AnnData
sc.pl.spatial(visium, color="cluster")

# we have some overhead here as we need to convert to numeric types for our framework
col_uns_name = 'cluster'
so.obs[spl]['cluster_id'] = so.obs[spl].groupby(col_uns_name).ngroup().astype('category')

# generate colormap
from matplotlib import cm
labs = so.obs[spl].groupby(['cluster_id']).head(1)[['cluster_id', col_uns_name]].set_index('cluster_id').to_dict()
cmap = ListedColormap([hex2color(i) for i in visium.uns[col_uns_name+'_colors']])
so.uns['cmaps'].update({'cluster_id':cmap})
so.uns['cmap_labels'].update({'cluster_id': labs[col_uns_name]})
so.uns['cmaps'].update({'default': cm.plasma})

# graph building, metrics, plot
sh.graph.build_graph(so, spl, mask_key=None)
sh.pl.spatial(so, spl, attr='cluster_id', edges=True, node_size=50, background_color='black')
sh.metrics.shannon(so, spl, attr='cluster_id')
sh.pl.spatial(so, spl, attr='shannon_cluster_id_knn', background_color='black', node_size= 50)

# %% alternative processing of Visium Fluorescence on cell-level
so2 = SpatialOmics()
spl = 'spl_0'
so2 = so.deepcopy()

# remove segmentation mask artefacts
mask = so2.masks[spl]['segmented_watershed'].squeeze().to_numpy().copy()
fig, ax = plt.subplots(1,2)
ax[0].imshow(mask > 0); fig.show()
mask[:1000, :1000] = 0
ax[1].imshow(mask > 0); fig.show()
fig.show()

# area of segmentations
area = pd.Series(mask.flatten()).value_counts()
area = counts[counts.index != 0]
plt.hist(area, bins=100);plt.show()

# color pixle according to the size of the cell they belong to
# segmentation is not good in Dentate region
mapping = area.to_dict()
mapping.update({0:0})
func = np.vectorize(lambda x: mapping[x], otypes=[int])
im = func(mask)
plt.imshow(im, cmap='jet');plt.show()

# determine quantile values
q995 = np.quantile(area, q=.995)
plt.hist(area[area < q995], bins=100);plt.show()

# highlight small objects
small_objs = area.index[(area < 200) & (area > 100)]
tmp = mask.copy().astype(int)
tmp[np.isin(mask, small_objs)] = -1 
plt.imshow(tmp < 0, cmap='gray');plt.show()

# highlight the environment of such a object
x, y = np.where(mask == 28662)
xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
pad = 200
plt.imshow(mask[xmin-pad:xmax+pad, ymin-pad:ymax+pad] > 0, cmap='gray');plt.show()

# based on previous inspection discard every segmentation < 100 and >q995
excl = area.index[(area > q995) | (area < 100)]
mask[np.isin(mask, excl)] = 0
print(len(np.unique(mask)))
so2.masks[spl]['cellmasks'] = mask
fig, axs = plt.subplots(1,2)
axs[0].imshow(so2.masks[spl]['segmented_watershed'].squeeze() > 0, cmap='gray')
axs[1].imshow(so2.masks[spl]['cellmasks'] > 0, cmap='gray')
fig.tight_layout();fig.show()
pd.Series(so2.masks[spl]['cellmasks'].reshape(-1)).value_counts()

# %%  extract image featuers for each object
from tqdm import tqdm

expr = so2.images[spl].squeeze().to_numpy()  # get rid of the z-dimension
mask = so2.masks[spl]['cellmasks']

ids = np.unique(mask)
ids = ids[ids != 0]

# extract single-cell expression values
res = []
for i in tqdm(ids):
    res.append(expr[mask == i].mean(0))

import pickle as pk
with open('single_cell_values_visium_cellmasks.pkl', 'wb') as f:
    pk.dump(res, f)

so2.X[spl] = pd.DataFrame(np.stack(res, axis=0), index=ids)
so2.X[spl].columns = ['channel'+str(i) for i in so2.X[spl].columns]

# %% compute cell-graph based on segmentation mask

sh.graph.build_graph(so2, spl)

# construct obs
g = so2.G[spl]['knn']
so2.obs[spl] = pd.DataFrame(index=g.nodes)
sh.pp.extract_centroids(so2, spl)

# scale to rgb range
img = so2.images[spl].squeeze().to_numpy().astype(float)
img = np.log1p(img)
IMG = img.copy()
for i in range(3):
    tmp = img[:,:,i]
    IMG[:,:,i] = tmp / tmp.max()

so2.uns['cmaps'].update({'default': cm.plasma})

fig, axs = plt.subplots(1, 3, figsize=(12,4), dpi=700)
visium_img.show('image', ax=axs[0])
axs[0].invert_yaxis()
sh.pl.spatial(so2, spl, attr='channel0', ax=axs[1], node_size=1)
sh.pl.spatial(so2, spl, attr='channel0', mode='mask', ax=axs[2], background_color='black')
fig.tight_layout()
fig.show()

# %% so to ad
so = sh.dataset.imc()
for spl in so.obs:
    sh.pp.extract_centroids(so,spl)
ad = so2ad(so)
ads = so2ad(so, one_adata=False)

# %%