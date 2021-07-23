"""Unit testing for module."""
import pytest
import numpy as np
from spatialOmics import SpatialOmics
from tifffile import imsave

@pytest.fixture(scope='module')
def so():
    return SpatialOmics()

@pytest.fixture(scope="session")
def image_file(tmpdir_factory):
    img = np.zeros((5,20,20), 'uint16')
    fn = tmpdir_factory.mktemp("data").join('image.tiff')
    imsave(fn, img)
    return fn

@pytest.fixture(scope="session")
def mask_file(tmpdir_factory):
    img = np.zeros((20,20))
    fn = tmpdir_factory.mktemp("data").join("mask.tiff")
    imsave(fn, img)
    return fn

@pytest.mark.parametrize('in_memory, to_store', [(True, True), (True, False), (False, True)])
def test_add_image(so: SpatialOmics, image_file, in_memory, to_store):
    so.add_image('spl1', str(image_file), in_memory=in_memory, to_store=to_store)
    so.get_image('spl1')

def test_add_mask(so: SpatialOmics, mask_file):
    so.add_mask('spl1', 'mask1', str(mask_file))
    so.get_mask('spl1', 'mask1')