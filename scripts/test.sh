#!/usr/bin/env bash

FILTER=$1

aws s3 mv s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/B+D/checkimages/model_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/B+D/checkimages/model_B+D_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz
aws s3 mv s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/B+D/checkimages/resid_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/B+D/checkimages/resid_B+D_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz
aws s3 mv s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/sersic_rg4/checkimages/model_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/sersic_rg4/checkimages/model_sersic_rg4_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz
aws s3 mv s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/sersic_rg4/checkimages/resid_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz s3://aurelien-sepp/ceers-full-grizli-v7.2/sepp/sersic_rg4/checkimages/resid_sersic_rg4_ceers-full-grizli-v7.2-$FILTER-clear_drc_sci_tile-full_1.fits.gz