# DJA x SourceXtractor++

## Introduction

This project looks into implementing SourceXtractor++ on the [DAWN JWST Archive](https://dawn-cph.github.io/dja/) (DJA). The DJA is a repository of public JWST galaxy data reduced and ready for science. This project aims at expanding the DJA with model fitting and more precise measurements on the different sets of photometric images.

For this, we use different well-known and documented programs :
* [SExtractor](https://www.astromatic.net/software/sextractor/)
* [PSFEx](https://www.astromatic.net/software/psfex/)
* [SourceXtractor++](https://github.com/astrorama/SourceXtractorPlusPlus)
* [Astropy](https://www.astropy.org/index.html)

## Usage

0. *(optional)* Cutouts in the image : [`00.2_Cutout.ipynb`](00.2_Cutout.ipynb)
1. Detect point-like sources with the F200W band + create the PSF for every band : [`03.1_SingleBand-PSF.ipynb`](03.1_SingleBand-PSF.ipynb)
2. Create association catalog and segmentation map based on the DJA catalog : [`04.1_AssociationCatalog.ipynb`](04.1_AssociationCatalog.ipynb)
3. Select the images (sci, wht, psf) for SourceXtractor++ : [`sepp-config.py`](config/sepp-config.py)
4. Run SourceXtractor++ in association mode : [`04.2_SE++Assoc.ipynb`](04.2_SE++Assoc.ipynb)
5. View and analyse the results : [`05_Analysis.ipynb`](05_Analysis.ipynb)
6. Compare the results to the DJA : [`06_Validation.ipynb`](06_Validation.ipynb)