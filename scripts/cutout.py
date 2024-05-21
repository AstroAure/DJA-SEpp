import sys
from astropy.coordinates import SkyCoord
import astropy.units as u
import dja_sepp

# python3 cutout.py ceers-full-grizli-v7.2 14h19m40s 52d51m30s 2 /FlashStorage/DJA-SEpp

def main():
    field = sys.argv[1]
    ra = sys.argv[2]
    dec = sys.argv[3]
    size_arcmin = float(sys.argv[4])
    home = sys.argv[5] if len(sys.argv)>5 else "/tmp"

    center = SkyCoord(f"{ra} {dec}", frame='icrs', unit='deg')
    size=u.Quantity((size_arcmin,size_arcmin), u.arcmin)
    dja_sepp.utils.save_cutouts(generic_filename=f"{home}/fields/{field}/image/*.fits", 
                                center=center, size=size,
                                plot=False, plot_str="sci", verbose=True)

if __name__=='__main__':
    main()