# %%
import pandas as pd
import warnings
import numpy as np
from matplotlib import cm

from .imc_data import IMCData
from smmhm.utils.general import make_iterable, is_categorical


# TODO: rewrite this robustly. Need to change _compute_metrics, neighbors to fetch data differently
def _get_obsX_data(ad: IMCData, core: str, attr, squeeze: bool = True, not_all_attr_found='raise'):
    obs = None
    attr = make_iterable(attr)
    attr = set(attr)
    if 'X' in attr:  # NOTE: this only works because we make the atts iterable (!), otherwise not save since 'prefix' in 'prefix_attr' is True and only 'prefix' in ['prefix_attr'] is False
        obs = ad.X[core]
        attr.remove('X')
    if 'obs' in attr:
        obs = ad.obs[core]
        attr.remove('obs')

    # accumulate single attrs from obs and X
    cols = ad.obs[core].columns
    if cols.isin(attr).any():
        obs = pd.concat((obs, ad.obs[core].loc[:, cols.isin(attr)]), axis=1)
        for i in cols[cols.isin(attr)]: attr.remove(i)

    cols = ad.X[core].columns
    if cols.isin(attr).any():
        obs = pd.concat((obs, ad.X[core].loc[:, cols.isin(attr)]), axis=1)
        for i in cols[cols.isin(attr)]: attr.remove(i)

    # check if we fond all attr
    if len(attr) > 0:
        if not_all_attr_found == 'suppress':
            pass
        elif not_all_attr_found == 'warn':
            warnings.warn(f'not all attributes ({attr}) found.')
        else:
            raise KeyError(f'not all attributes ({attr}) found.')

    return obs.squeeze() if squeeze else obs


def _get_core_data(ad: IMCData, core: str, attr, squeeze: bool = False, not_all_attr_found='raise'):
    # add core level data
    obs = None
    attr = make_iterable(attr)
    attr = set(attr)

    cols = ad.obsc.columns
    if cols.isin(attr).any():
        obs = pd.concat((obs, ad.obsc.loc[core, cols.isin(attr)]))
        for i in cols[cols.isin(attr)]: attr.remove(i)

    cols = ad.meta.columns
    if cols.isin(attr).any():
        obs = pd.concat((obs, ad.meta.loc[core, cols.isin(attr)]))
        for i in cols[cols.isin(attr)]: attr.remove(i)

    # check if we fond all attr
    if len(attr) > 0:
        if not_all_attr_found == 'suppress':
            pass
        elif not_all_attr_found == 'warn':
            warnings.warn(f'not all attributes ({attr}) found.')
        else:
            raise KeyError(f'not all attributes ({attr}) found.')

    return obs.squeeze() if squeeze else obs


def _get_data(ad: IMCData, core: str, attr, squeeze: bool = True, source='both', not_all_attr_found='raise'):
    VALID_SOURCES = ['local', 'global', 'both', 'local_first', 'global_first']

    # #TODO: should we us
    # if attr == 'X': #NOTE: for backwards compatibility with quadratic entropy
    #     # container = [ad.X[c] for c in ad.X.keys()]
    #     # return pd.concat(container)
    #     return ad.X[core]

    if source == 'local':
        data = _get_obsX_data(ad, core, attr, squeeze=False, not_all_attr_found='raise')
    elif source == 'global':
        data = _get_core_data(ad, core, attr, squeeze=False, not_all_attr_found='raise')
    elif source == 'both':
        # NOTE: if we have the same attribute name in local and global data structures they are combined!
        local_dat = _get_obsX_data(ad, core, attr, squeeze=False, not_all_attr_found='suppress')
        core_dat = _get_core_data(ad, core, attr, squeeze=False, not_all_attr_found='suppress')

        if local_dat is not None and core_dat is not None:
            # reshape to concat
            tmp = np.atleast_2d(core_dat.values).repeat(len(local_dat), axis=0)
            df = pd.DataFrame(tmp, index=local_dat.index, columns=core_dat.index)

            # concat
            data = pd.concat((local_dat, df), axis=1)

            # check if some attributes were in local and core data
            if (data.columns.value_counts() > 1).any():
                # duplicated_feat = data.columns[(data.columns.value_counts() > 1)]
                warnings.warn(f'features {(data.columns)} have been found in local and core data containers.')
        elif local_dat is not None:
            data = local_dat
        elif core_dat is not None:
            data = core_dat
        else:
            raise KeyError(f'not all attributes ({attr}) found.')
    else:
        raise ValueError(f'{source} is not a valid source. Available are {VALID_SOURCES}.')

    # check if all features found
    # cols = data.columns if hasattr(data, 'columns') else data.index
    # attr = set(make_iterable(attr))
    # for i in cols: attr.remove(i)
    # if len(attr) > 0:
    #     if not_all_attr_found == 'suppress':
    #         pass
    #     elif not_all_attr_found == 'warn':
    #         warnings.warn(f'not all attributes ({attr}) found.')
    #     else:
    #         raise KeyError(f'not all attributes ({attr}) found.')

    return data.squeeze() if squeeze else data


def _get_cmap(ad: IMCData, attr: str, data):
    '''
    Return the cmap and cmap labels for a given attribute if available, else a default

    Parameters
    ----------
    ad: IMCData
        ad object form which to fetch the data
    core: str
        core for which to get data
    attr: str
        attribute for which to get the cmap and cmap labels if available

    Returns
    -------
    cmap and cmap labels for attribute

    '''

    # TODO: recycle cmap if more observations than colors
    cmap, cmap_labels = None, None
    if attr in ad.uns['cmaps'].keys():
        cmap = ad.uns['cmaps'][attr]
    elif is_categorical(attr):  # check if data is cateogrical
        cmap = ad.uns['cmaps']['category']
    else:
        cmap = ad.uns['cmaps']['default']

    # get callable cmap if string
    if isinstance(cmap, str):
        cmap = cm.get_cmap(cmap)

    if attr in ad.uns['cmap_labels'].keys():
        cmap_labels = ad.uns['cmap_labels'][attr]

    return cmap, cmap_labels


def _recycle_colormap(cmap, n):
    pass


def _get_loc(ad: IMCData, core: str, cell_ids=None):
    if cell_ids:
        return ad.obs[core].loc[cell_ids, ['x', 'y']]
    else:
        return ad.obs[core][['x', 'y']]


def get_uns(self, metric=None, key=None, uns_path=None, cores=None):
    if cores is None:
        cores = self.cores

    cores = make_iterable(cores)

    if uns_path:
        path = uns_path.split('/')
    else:
        path = [metric, key]

    data = []
    for core in cores:
        cur = self.uns[core]
        for i in range(len(path) - 1):
            cur = cur[path[i]]

        data.append(cur[path[-1]].squeeze())

    df = pd.concat(data, axis=1)
    df.columns = cores
    return df


# %%