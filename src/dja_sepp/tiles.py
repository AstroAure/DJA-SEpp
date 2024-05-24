import os
import re
import fnmatch
import glob
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, pixel_to_skycoord
from astropy.visualization import ImageNormalize, MinMaxInterval, ZScaleInterval, LogStretch
from astropy.coordinates import SkyCoord
from astropy.table import vstack
import astropy.units as u
from reproject import reproject_interp, reproject_exact
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from . import utils

def create_tiles(data, wcs, centers, sizes):
    """
    Create the tiles out of a big image and given center position and sizes of tiles

    data : data of the main image
    wcs : WCS of the main image
    centers : list of SkyCoord for the center of each tile
    sizes : list of tuple of angular sizes for each tile

    Returns : list of Cutout2D representing the tiles
    """
    tiles = []
    for i in range(len(centers)):
        cutout = Cutout2D(data, position=centers[i], size=sizes[i], wcs=wcs, mode='partial', fill_value=0.0)
        tiles.append(cutout)
    return tiles

def tile_positions(wcs, nx, ny, overlap):
    """
    Calculates the centers and sizes of tiles based on a grid

    wcs : WCS of the main image
    nx/y : number of tiles on the x/y axis
    overlap : angular size of the overlap between tiles

    Returns : lists of centers and sizes for the tiles
    """
    width, height = wcs.pixel_shape
    pixscale = abs(proj_plane_pixel_scales(wcs)[0])*u.degree
    overlap_px = int(overlap/pixscale)
    width_tile = (width-overlap_px)/nx + overlap_px
    height_tile = (height-overlap_px)/ny + overlap_px
    centers = []
    sizes = []
    for j in range(ny):
        for i in range(nx):
            center = pixel_to_skycoord(0.5*width_tile + i*(width_tile-overlap_px), 0.5*height_tile + j*(height_tile-overlap_px), wcs)
            centers.append(center)
            sizes.append((width_tile*pixscale, height_tile*pixscale))
    return centers, sizes

def tile_grid(wcs, max_size_x, max_size_y, overlap):
    """
    Calculates the grid size for a given maximum angular tile size

    wcs : WCS of the main image
    max_size_x/y : maximum angular size on the x/y axis for the tiles (without overlap)
    overlap : angular size of the overlap between tiles

    Returns : number of tiles on the x and y axis
    """
    width_px, height_px = wcs.pixel_shape
    pixscale = abs(proj_plane_pixel_scales(wcs)[0])*u.degree
    width, height = width_px*pixscale, height_px*pixscale
    nx, ny = int(width/max_size_x)+1, int(height/max_size_y)+1
    return nx, ny

def save_tiles(tiles, main_filename, save_folder=None, verbose=False):
    """
    Saves the tiles as FITS images

    tiles : list of Cutout2D representing the tiles
    main_filename : filename of the main image
    save_folder : folder to save the cutouts in.
        By default, tiles are saved in a `cutout` sub-folder where the main image is.
    """
    hdu = fits.open(main_filename, memmap=True)[0]
    for i, tile in enumerate(tiles):
        if verbose: print(f"Tile {i+1} / {len(tiles)}")
        hdu.data = tile.data
        hdu.header.update(tile.wcs.to_header())
        parsed_name = main_filename.split("/")
        folder = "/".join(parsed_name[:-1])
        init_name = parsed_name[-1]
        name = f"{'.'.join(init_name.split('.')[:-1])}_tile-{i}.fits"
        if save_folder is None: save_folder = f"{folder}/tiles"
        if not os.path.exists(save_folder): os.mkdir(save_folder) 
        hdu.writeto(f"{save_folder}/{name}", overwrite=True)

def plot_tiles(nx, ny, tiles, plot_main=False, data=None, wcs=None):
    """
    Plot the tiles for visual inspection

    nx/y : number of tiles on the x/y axis
    tiles : list of Cutout2D representing the tiles
    plot_main : whether to plot the main image or not
    data : data array of the main image (required if plot_main==True)
    wcs : WCS of the main image (required if plot_main==True)
    """
    h, w = tiles[0].data.shape
    plot_x_size = min(12, 2*nx)
    plot_y_size = min(8, plot_x_size*ny/nx*h/w)
    plot_x_size = min(12, plot_y_size*nx/ny*w/h)
    if plot_main:
        fig = plt.figure(layout='constrained', figsize=(plot_x_size,plot_y_size/2))
        gs = GridSpec(ny, nx*2, figure=fig)
        ax_main = fig.add_subplot(gs[:,nx:], projection=wcs)
        norm = ImageNormalize(data, interval=ZScaleInterval())
        ax_main.imshow(data, cmap='gray', origin='lower', norm=norm)
        utils.add_good_scalebar(ax_main, wcs)
        ax_main.set_xlabel('Right Ascension (J2000)')
        ax_main.set_ylabel('Declination (J2000)')
    else:
        fig = plt.figure(layout='constrained', figsize=(plot_x_size,plot_y_size))
        gs = GridSpec(ny, nx, figure=fig)
        norm = ImageNormalize(tiles[0].data, interval=ZScaleInterval())
    for j in range(ny):
        for i in range(nx):
            tile = tiles[i+nx*j]
            ax = fig.add_subplot(gs[ny-j-1,i], projection=tile.wcs)
            ax.imshow(tile.data, cmap='gray', origin='lower', norm=norm)
            ax.set_axis_off()
            if (j==ny-1) & (i==nx-1):
                utils.add_good_scalebar(ax, tile.wcs)
    return fig

def batch_tiling(generic_filename   : str,
                 tile_max_size      : u.Quantity,
                 overlap            : u.Quantity,
                 save_folder        : str = None, 
                 plot               : bool = False,
                 plot_str           : str = None,
                 verbose            : bool = False) -> None :
    """
    Create tiles for all selected images, 
    centered on the same points and with the same angular size.

    * `generic_filename` (`str` with wildcard) : 
        Generic filename for the images to process.
        Can use wildcards (*,?) readable by `glob.glob`.
    * `tile_max_size` (`u.Quantity`) :
        Maximum angular size of the tiles (without overlap)
    * `overlap` (`u.Quantity`) :
        Angular size of overlap between tiles
    * `save_folder` (`str`, optional):
        Folder to save the cutouts in.
        By default, the images will be saved in a `cutout` sub-folder of the initial folder.
    * `plot` (`bool`, optional):
        Plots the tiles for the first band with MinMax and Log stretching.
    * `plot_str` (`str`, optional):
        Regex to choose the image to plot. Can use wildcards (*,?) readable by `glob.glob`.
    * `verbose` (`bool`, optional)
    """
    image_list = glob.glob(generic_filename)
    n = len(image_list)
    if verbose: print(f"Number of images found : {n}")
    if verbose: print(f"Images found :")
    if verbose: 
        for img in image_list: print(img)
    for i,img in enumerate(image_list):
        print(f"Image {i+1} / {n}")
        # Load image
        hdu = fits.open(img, memmap=True)[0]
        wcs = WCS(hdu.header)
        data = hdu.data
        # Calculate tiling on first image
        if i==0:
            nx, ny = tile_grid(wcs, tile_max_size, tile_max_size, overlap)
            centers, sizes = tile_positions(wcs, nx, ny, overlap)
        # Tiles the image
        tiles = create_tiles(data, wcs, centers, sizes)
        # Save tiles
        save_tiles(tiles, img, save_folder, verbose=verbose)
        if plot & fnmatch.fnmatch(img, plot_str):
            fig = plot_tiles(nx, ny, tiles, plot_main=False)
    if plot: plt.show()
    if plot: return fig

def merge_catalogs(cat1, cat2, max_sep=0.25*u.arcsec, filter_merge='f200w', return_matches=False):
    """
    Merge two tile catalog by finding matches and selecting between these matches by lowest mag error

    cat1/cat2 : catalogs (astropy.table.Table) to merge
    max_sep : maximum separation between sources for cross match
    filter_merge : filter to use for the selection using lowest mag error
    return_matches : whether to return the indices of the sources in cat1/cat2 matched

    Returns : merged catalog (astropy.table.Table)
    """
    # Creating columns that are in one catalog but not the other, and in list-columns (missing filters)
    keys1 = cat1.keys()
    filter_list1 = utils.get_filter_list(keys1)
    keys2 = cat2.keys()
    filter_list2 = utils.get_filter_list(keys2)
    filter_list = list(set(filter_list1 + filter_list2))
    filter_list.sort()
    list_cols = []
    for key in cat1.keys():
        if len(cat1[key].shape)>1:
            if cat1[key].shape[1]==len(filter_list1):
                list_cols.append(key)
    for i, filter in enumerate(filter_list):
        if filter not in filter_list1:
            for col in list_cols:
                cat1[col] = np.insert(ma.array(cat1[col]), i, ma.masked, axis=1)
        if filter not in filter_list2:
            for col in list_cols:
                cat2[col] = np.insert(ma.array(cat2[col]), i, ma.masked, axis=1)
    # Matching sources
    coord1 = SkyCoord(cat1['world_centroid_alpha']*u.degree, cat1['world_centroid_delta']*u.degree)
    coord2 = SkyCoord(cat2['world_centroid_alpha']*u.degree, cat2['world_centroid_delta']*u.degree)
    idx, d2d, _ = coord1.match_to_catalog_sky(coord2)
    match1 = np.where(d2d<max_sep)[0]
    match2 = idx[match1]
    # Selecting matched sources based on mag uncertainty
    merge_choice = cat1[match1][f'MAG_MODEL_{filter_merge.upper()}_err']<cat2[match2][f'MAG_MODEL_{filter_merge.upper()}_err']
    merged_1 = cat1[match1][merge_choice]
    merged_2 = cat2[match2][~merge_choice]
    # Selecting unmatched sources
    cat1_unmatched = cat1[[i for i in range(len(cat1)) if i not in match1]]
    cat2_unmatched = cat2[[i for i in range(len(cat2)) if i not in match2]]
    # Merging catalogs
    merged = vstack([cat1_unmatched, cat2_unmatched, merged_1, merged_2], join_type='outer')
    if return_matches: return merged, match1, match2
    return merged
    # return cat1

def merge_tiles(catalogs, max_sep=0.25*u.arcsec, filter_merge='f200w'):
    """
    Merge tile catalogs by finding matches and selecting between these matches by lowest mag error

    catalogs : catalogs (astropy.table.Table) to merge
    max_sep : maximum separation between sources for cross match
    filter_merge : filter to use for the selection using lowest mag error

    Returns : merged catalog (astropy.table.Table)
    """
    merged = catalogs[0]
    for i in range(1, len(catalogs)):
        merged = merge_catalogs(catalogs[i], merged, max_sep, filter_merge)
    return merged

def merge_images(folder, filter_list, type, wcs=None, shape=None, exact=False, out_folder=None, verbose=False):
    """
    Mosaic tile images (data, model, residual) by reprojection and co-addition

    folder : folder of the images to mosaic
    filter_list : list of filters (in image names) to use for mosaicing process
    type : type of images to mosaic (fromerly, prefix of image name)
    wcs : WCS for final image
    shape : shape for final image
    exact : use `reproject.reproject_exact` for reprojection (by default, `reproject.reproject_interp`)
    out_folder : folder to save the mosaiced images. By default, same folder as `folder`
    verbose : verbose
    """
    if (wcs is None) or (shape is None):
        if verbose: print("Finding optimal WCS")
        images = []
        for filter in filter_list:
            images.extend(glob.glob(f"{folder}/*{type}*{filter}*sci*[!tile-full]*"))
        wcs, shape = find_optimal_celestial_wcs(images, auto_rotate=True)
    for filter in filter_list:
        if verbose: print(f"---- {filter.upper()} ----")
        images = glob.glob(f"{folder}/*{type}*{filter}*sci*tile-[!full]*")
        if verbose: print(f"{filter.upper()} : {images}")
        if verbose: print("Mosaicing")
        reproject_function = reproject_exact if exact else reproject_interp
        mosaic, footprint = reproject_and_coadd(images, wcs, shape_out=shape, 
                                                reproject_function=reproject_function,
                                                match_background=False)
        if verbose: print("Saving")
        full = fits.PrimaryHDU(mosaic, wcs.to_header())
        out_folder = folder if out_folder is None else out_folder
        name = re.sub('tile-\d+', 'tile-full', images[0].split('/')[-1])
        full.writeto(f"{out_folder}/{name}", overwrite=True)
    return wcs, shape