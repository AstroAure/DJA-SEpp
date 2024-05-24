import importlib.resources
import os
import fnmatch
import glob
from typing import Union
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from astropy.io import fits, ascii
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, MinMaxInterval, ZScaleInterval, LogStretch
from astropy.visualization.wcsaxes import add_scalebar
import astropy.units as u

# Color dictionaries to make the plot nicer with colored filter names
color_dict = {'F090W':'#0000E3', 
              'F115W':'#1020FF',
              'F150W':'#0162FF',
              'F182M':'#00B1FF',
              'F200W':'#04D1FF',
              'F210M':'#02E8F4',
              'F250M':'#55FFAE',
              'F277W':'#87FF7F',
              'F300M':'#AAFF57',
              'F335M':'#EBFF0C',
              'F356W':'#FFD807',
              'F410M':'#FF6D03',
              'F430M':'#FF3F00',
              'F444W':'#FF330C',
              'F460M':'#DF0606',
              'F480M':'#B20000',}
channel_dict = {'F090W':'SW', 'F115W':'SW', 'F150W':'SW', 'F182M':'SW', 'F200W':'SW', 'F210M':'SW',
                'F250M':'LW', 'F277W':'LW', 'F300M':'LW', 'F335M':'LW', 'F356W':'LW', 'F410M':'LW', 'F430M':'LW', 'F444W':'LW', 'F460M':'LW', 'F480M':'LW'}
channel_color_dict = {'SW' : '#0162FF', 'LW' : '#FF330C'}

filters_waveband = {'F090W': {'pivot': 0.901, 'band': 0.194},
                    'F115W': {'pivot': 1.154, 'band': 0.225},
                    'F150W': {'pivot': 1.501, 'band': 0.318},
                    'F182M': {'pivot': 1.845, 'band': 0.238},
                    'F200W': {'pivot': 1.990, 'band': 0.461},
                    'F210M': {'pivot': 2.093, 'band': 0.205},
                    'F250M': {'pivot': 2.503, 'band': 0.181},
                    'F277W': {'pivot': 2.786, 'band': 0.672},
                    'F300M': {'pivot': 2.996, 'band': 0.318},
                    'F335M': {'pivot': 3.365, 'band': 0.347},
                    'F356W': {'pivot': 3.563, 'band': 0.787},
                    'F410M': {'pivot': 4.092, 'band': 0.436},
                    'F430M': {'pivot': 4.280, 'band': 0.228},
                    'F444W': {'pivot': 4.421, 'band': 1.024},
                    'F460M': {'pivot': 4.624, 'band': 0.228},
                    'F480M': {'pivot': 4.834, 'band': 0.303},}
# Source : https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters

def add_good_scalebar(ax, wcs, color='white', fraction=0.25, corner='bottom right', pad=0.1, fontsize='x-large'):
    """
    Utility function for astropy.visualization.wcsaxes.add_scalebar
    to automaticaly resize the scalebar to the size of the image.

    ax : matplotlib.pyplot.Axes where to draw the scale bar
    wcs : WCS of the image
    color : color of the scale bar
    fraction : fraction of the image to aim for the scale bar (will always be smaller)
    corner : where to place the scale bar
    pad : padding around the scale bar
    fontsize : fontsize of the scale text
    """
    width = abs(proj_plane_pixel_scales(wcs)[0]*wcs.pixel_shape[0])*u.degree
    good_values = [0.5*u.arcsec, 1*u.arcsec, 2*u.arcsec, 3*u.arcsec, 5*u.arcsec, 10*u.arcsec, 20*u.arcsec,
                0.5*u.arcmin, 1*u.arcmin, 2*u.arcmin, 3*u.arcmin, 5*u.arcmin, 10*u.arcmin, 20*u.arcmin,
                0.5*u.degree, 1*u.degree, 2*u.degree, 5*u.degree, 10*u.degree]
    dist = np.array([(fraction*width-val).value for val in good_values])
    dist[dist<0] = np.inf
    size = good_values[np.argmin(dist)]
    add_scalebar(ax, size, label=f"{size:latex}", color=color, fontproperties=FontProperties(size=fontsize), label_top=True, pad=pad, corner=corner)

def save_cutouts(generic_filename   : str, 
                 center             : Union[tuple[str],SkyCoord] , 
                 size               : Union[int,u.Quantity],
                 save_folder        : str = None, 
                 plot               : bool = False,
                 plot_str           : str = None,
                 verbose            : bool = False) -> None :
    """
    Create cutouts for all selected images, 
    centered on the same point and with the same size.

    * `generic_filename` (`str` with wildcard) : 
        Generic filename for the images to process.
        Can use wildcards (*,?) readable by `glob.glob`.
    * `center` (`tuple[str]` or `SkyCoord`) :
        Center of the cutouts.
        Can be in pixel units, or an `astropy.coordinates.SkyCoord` object (e.g. RA-DEC).
    * `size` (`int` or `Quantity`):
        Size of the cutouts.
        Can be in pixel units, or an `astropy.units.Quantity` object (e.g. angular size).
        WARNING : If pixel scale is different between images, giving a size in pixel units will result in non-aligned images.
    * `save_folder` (`str`, optional):
        Folder to save the cutouts in.
        By default, the images will be saved in a `cutout` sub_folder of the initial folder.
    * `plot` (`bool`, optional):
        Plots the cutouts with MinMax and Log stretching.
    * `plot_str` (`str`, optional):
        String to be found in image name to plot it. 
        Only active if `plot=True`. By default, all images are plotted.
    * `verbose` (`bool`, optional)
    """
    image_list = glob.glob(generic_filename)
    n = len(image_list)
    if verbose: print(f"Number of images found : {n}")
    if verbose: print(f"Images found :")
    if verbose: 
        for img in image_list: print(img)
    if plot: 
        len_plot = np.sum([fnmatch.fnmatch(img, plot_str) for img in image_list])
        size_plot = 1
        while size_plot*size_plot < len_plot: size_plot += 1
        fig, j = plt.figure(figsize=(3*size_plot,3*size_plot)), 1
    for i,img in enumerate(image_list):
        print(f"Image {i+1} / {n}")
        # Load image
        hdu = fits.open(img, memmap=True)[0]
        wcs = WCS(hdu.header)
        # Cutout
        cutout = Cutout2D(hdu.data, 
                    position=center, size=size, 
                    wcs=wcs,
                    mode='partial', fill_value=0.0)
        # Save cutout
        hdu.data = cutout.data
        hdu.header.update(cutout.wcs.to_header())
        parsed_name = img.split("/")
        folder = "/".join(parsed_name[:-1])
        init_name = parsed_name[-1]
        name = f"{'.'.join(init_name.split('.')[:-1])}_cutout.fits"
        if save_folder is None: save_folder = f"{folder}/cutout"
        if not os.path.exists(save_folder): os.mkdir(save_folder) 
        hdu.writeto(f"{save_folder}/{name}", overwrite=True)
        # Plot cutout
        if plot:
            if plot_str is None or fnmatch.fnmatch(img, plot_str):
                ax = fig.add_subplot(size_plot,size_plot,j, projection=cutout.wcs)
                norm = ImageNormalize(cutout.data, interval=ZScaleInterval())
                ax.imshow(cutout.data, cmap='gray', origin='lower', norm=norm)
                ax.set_title(init_name, size='xx-small')
                ax.set_axis_off()
                add_good_scalebar(ax, cutout.wcs)
                j+=1
    if plot: fig.tight_layout()
    if plot: plt.show()
    if plot: return fig

def show_source(id, cat, filter_list,
                data_folder,
                model_folder,
                resid_folder,
                data_suffix="sci.fits",
                model_suffix="sci_1.fits",
                resid_suffix="sci_1.fits",
                fov=3*u.arcsec, 
                ra_name='world_centroid_alpha', de_name='world_centroid_delta',):
    """
    Displays data, model and residual images for each band,
    for the source at a given index in a given catalog.

    id : index of the source in cat
    cat : catalog from SE++
    filter_list : list of the filters (in lowercase) to use for the display
    fov : field of view of the cutouts around the source.
          Must be a astropy.units.Quantity
    ra/de_name : column name in cat for right ascension and declination
    data/model/resid_folder : folder to the data/model/resid images
    data/model/resid_suffix : suffix for image name selection (can contain glob.glob wildcards)
    """
    target = SkyCoord(cat[ra_name][id], cat[de_name][id], frame='icrs', unit='deg')
    size = (1.1*fov,1.1*fov)

    fig = plt.figure(figsize=(2*len(filter_list)+1,2*3+0.5))
    for i, filter in enumerate(filter_list):
        # Loads images
        hdu       = fits.open(glob.glob(f"{data_folder}/*{filter}*{data_suffix}")[0], memmap=True)[0]
        model_hdu = fits.open(glob.glob(f"{model_folder}/model*{filter}*{model_suffix}")[0], memmap=True)[0]
        resid_hdu = fits.open(glob.glob(f"{resid_folder}/resid*{filter}*{resid_suffix}")[0], memmap=True)[0]
        wcs = WCS(hdu.header)
        # Cutouts in images
        cutout       = Cutout2D(hdu.data,       position=target, size=size, wcs=wcs, mode='partial', fill_value=0.0)
        model_cutout = Cutout2D(model_hdu.data, position=target, size=size, wcs=wcs, mode='partial', fill_value=0.0)
        resid_cutout = Cutout2D(resid_hdu.data, position=target, size=size, wcs=wcs, mode='partial', fill_value=0.0)
        # Plots data
        ax = fig.add_subplot(3,len(filter_list),i+1, projection=cutout.wcs)
        norm = ImageNormalize(cutout.data, interval=ZScaleInterval())
        # norm = ImageNormalize(cutout.data, interval=MinMaxInterval(), stretch=LogStretch())
        ax.imshow(cutout.data, cmap='gray', origin='lower', norm=norm)
        add_good_scalebar(ax, cutout.wcs, fraction=1.0, corner='bottom')
        ax.set_title(filter.upper(), fontsize=24)
        ax.set_axis_off()
        # Plots model
        ax_model = fig.add_subplot(3,len(filter_list),i+len(filter_list)+1, projection=cutout.wcs)
        ax_model.imshow(model_cutout.data, cmap='gray', origin='lower', norm=norm)
        ax_model.set_axis_off()
        # Plots residual
        ax_resid = fig.add_subplot(3,len(filter_list),i+2*len(filter_list)+1, projection=cutout.wcs)
        ax_resid.imshow(resid_cutout.data, cmap='gray', origin='lower', norm=norm)
        ax_resid.set_axis_off()
    plt.figtext(1.0, 0.79, "Data", figure=fig, va='center', ha='left', rotation='vertical', fontsize='xx-large')
    plt.figtext(1.0, 0.48, "Model", figure=fig, va='center', ha='left',rotation='vertical', fontsize='xx-large')
    plt.figtext(1.0, 0.17, "Residual", figure=fig, va='center', ha='left',rotation='vertical', fontsize='xx-large')
    plt.figtext(0.0, 0.48, f"Source #{id}", figure=fig, va='center', ha='right', rotation='vertical', fontsize=24)
    fig.tight_layout()
    plt.show()
    return fig

def plot_filters(ax, filter_list=None, scale=1, names=True, 
                 throughput_folder="data/NIRCam_Filters-Throughput"):
    """
    Plots normalized filter throughputs on ax
    Source : https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters

    ax : matplotlib.pyplot.Axes to plot on
    filter_list : list of the filters to plot
    scale : maximum value of the plotted throughput curves
    names : whether to display the names of the filters
    throughput_folder : folder with the throughput data downloaded from the JWST documentation
    """
    throughput_folder = importlib.resources.files('dja_sepp').joinpath(throughput_folder)
    filter_list = [filter.lower() for filter in color_dict.keys()] if filter_list is None else filter_list
    for filter in filter_list:
        table = ascii.read(f"{throughput_folder}/{filter.upper()}_mean_system_throughput.txt")
        color = color_dict[filter.upper()] if filter.upper() in color_dict else '#000000'
        max_val = scale if 'w' in filter else scale/2
        ax.plot(table['Microns'], max_val*table['Throughput']/table['Throughput'].max(), color=color)
        ax.fill_between(table['Microns'], max_val*table['Throughput']/table['Throughput'].max(), alpha=0.4, color=color)
        if names: ax.text(np.average(table['Microns'], weights=table['Throughput']), 0.55*max_val, filter.upper(), ha='center', va='center', c='k', size='x-small')

def plot_photometric_spectrum(id, cat, filter_list, mag=False, custom_ax=None, title=None):
    """
    Plots photometric spectrum of cat[id] in flux or mag depending of `mag` (default, flux).
    The value of the wavelength and bandwidth associated to each band comes from the JWST documentation.
    Source : https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters

    id : index of the source in cat
    cat : catalog from SE++
    filter_list : list of the filters to plot
    mag : whether to plot the spectrum in AB mag or flux (default)
    custom_ax : matplotlib.pyplot.Axes to plot on.
                By default, the function creates a new figure
    title : title to set for the plot
    """

    fig, ax = plt.subplots(figsize=(12,3)) if custom_ax is None else (None, custom_ax)
    # DJA images pixel value are in 10 nJy
    flux_target = [10*cat[id-1][f"FLUX_MODEL_{filter.upper()}"] for filter in filter_list]
    flux_err_target = [10*cat[id-1][f"FLUX_MODEL_{filter.upper()}_err"] for filter in filter_list]
    mag_target = [cat[id-1][f"MAG_MODEL_{filter.upper()}"] for filter in filter_list]
    mag_err_target = [cat[id-1][f"MAG_MODEL_{filter.upper()}_err"] for filter in filter_list]

    plot_filters(ax, filter_list, (max(mag_target) if mag else max(flux_target))/5)
    if mag:
        ax.errorbar([filters_waveband[filter.upper()]['pivot'] for filter in filter_list], mag_target, 
                    xerr=np.array([filters_waveband[filter.upper()]['band'] for filter in filter_list])/2, yerr=mag_err_target, 
                    fmt='o', ms=3, c='k', capsize=3, capthick=1)
    else:
        ax.errorbar([filters_waveband[filter.upper()]['pivot'] for filter in filter_list], flux_target, 
                    xerr=np.array([filters_waveband[filter.upper()]['band'] for filter in filter_list])/2, yerr=flux_err_target, 
                    fmt='o', ms=3, c='k', capsize=3, capthick=1)
    ax.set_xlabel("Wavelength (µm)")
    if mag:
        ax.set_ylabel("AB Mag")
    else:
        ax.set_ylabel(r"$F_{\nu}$ (nJy)")
    ax.set_title(title)
    if custom_ax is None: plt.show()
    if custom_ax is None: return fig

def plot_group_filter(filter_list, plot_func, **kwargs):
    """
    Make plots for different filters

    filter_list : list of filter names (in lowercase)
    plot_func : function defining what to plot for one filter.
            Its first two arguments must be a plt.Axes and the filter name
    **kwargs : additionnal arguments to pass to plot_func

    Returns : fig, axs of the figure
    """
    channel, count = np.unique([channel_dict[filter.upper()] for filter in filter_list], return_counts=True)
    channel_count = dict(zip(channel, count))
    w = min(5,channel_count[max(channel_count, key=channel_count.get)])
    h = sum([-(channel_count[channel]//-w) for channel in channel_count])
    fig, axs = plt.subplots(h,w,figsize=(3*w+1,3*h+1), sharex=True, sharey=True, gridspec_kw = {'wspace':0, 'hspace':0})
    indices = [i for i in range(len(axs.flatten()))]
    for i, filter in enumerate(filter_list):
        if (len(channel_count)==2) & (channel_dict[filter.upper()]=="LW"):
            mod = channel_count["SW"]%w
            mod = w if mod==0 else mod
            i += max(0, w-mod)
        indices.remove(i)
        ax = axs.flatten()[i] if type(axs)==np.ndarray else axs
        plot_func(ax, filter, **kwargs)
    for i in indices:
        axs.flatten()[i].set_axis_off()
    return fig, axs

def get_filter_list(keys):
    """
    Find filter list from column names of a SE++ catalog

    keys : list of column names (obtained throught cat.keys())

    Returns : list of filter names (in lowercase)
    """
    filter_list = [key.split("_")[-1].lower() for key in keys if (('MAG_MODEL' in key) and ('err' not in key) and ('BULGE' not in key) and ('DISK' not in key))]
    filter_list.sort()
    return filter_list


def main():
    fig, ax = plt.subplots(figsize=(12,5))
    plot_filters(ax, ['f090w', 'f115w', 'f200w', 'f277w', 'f356w', 'f444w'])
    plt.show()

if __name__=="__main__":
    main()