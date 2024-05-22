import sys
from astropy.coordinates import SkyCoord
import astropy.units as u
import dja_sepp

# python3 tile.py ceers-full-grizli-v7.2 5 0.5 /FlashStorage/DJA-SEpp true

def main():
    field = sys.argv[1]
    tile_size_arcmin = sys.argv[2]
    overlap_arcmin = sys.argv[3]
    home = sys.argv[4] if len(sys.argv)>4 else "/tmp"
    plot = (sys.argv[5].lower()=='true') if len(sys.argv)>5 else False
    img_subpath = sys.argv[6] if len(sys.argv)>6 else ""

    fig = dja_sepp.tiles.batch_tiling(generic_filename=f"{home}/fields/{field}/image/{img_subpath}/*.fits",
                                      tile_max_size=tile_size_arcmin*u.arcmin, overlap=overlap_arcmin*u.arcmin,
                                      plot=plot, plot_str='*ir*sci*', verbose=True)
    if plot: fig.savefig(f"{home}/fields/{field}/image/{img_subpath}/tiles/{field}_tile.png", bbox_inches=0, pad_inches=0, dpi=200)

if __name__=='__main__':
    main()