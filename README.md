# DJA x SourceXtractor++

## Introduction

This project looks into implementing SourceXtractor++ on the [DAWN JWST Archive](https://dawn-cph.github.io/dja/) (DJA). The DJA is a repository of public JWST galaxy data reduced and ready for science. This project aims at expanding the DJA with model fitting and more precise measurements on the different sets of photometric images.

## Usage

0. *(optional)* Cutouts in the image : [`00.2_Cutout.ipynb`](00.2_Cutout.ipynb)
1. Detect point-like sources with the F200W band + create the PSF for every band : [`03.1_SingleBand-PSF.ipynb`](03.1_SingleBand-PSF.ipynb)
2. Run SourceXtractor++ in detection mode : [`04_SE++.ipynb`](04_SE++.ipynb)
3. Compare the results to the DJA : [`06_Validation.ipynb`](06_Validation.ipynb)
4. View and analyse the results : [`05_Analysis.ipynb`](05_Analysis.ipynb)

## Dependencies

* [SExtractor](https://www.astromatic.net/software/sextractor/)
* [PSFEx](https://www.astromatic.net/software/psfex/)
* [SourceXtractor++](https://github.com/astrorama/SourceXtractorPlusPlus)
* [Astropy](https://www.astropy.org/index.html)

## Utilization with AWS EC2

One of the goal of this project is also to run SourceXtractor++ on [AWS EC2](https://aws.amazon.com/ec2/) to make it run in the background and faster than on your personal computer. 

A benchmark is currently being run to find how to set the parameters in SE++ and AWS EC2 to minimize runtime depending on the size of your images and the number of bands wanted.

To know how to use easily AWS EC2 in combination with Jupyter, please refer to [this series of tutorial](https://github.com/AstroAure/VSJupytEC2).