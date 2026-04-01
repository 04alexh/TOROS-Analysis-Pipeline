# TOROS-Analysis-Pipeline
Full variability analysis pipeline for the TOROS telescope in Argentina. These are updated versions of the previous programs I posted before; assume any TOROS programs I make that are not within this repo are depreciated.

## Structure
Programs have been divided into two categories: image reduction programs and image analysis programs. The image reduction programs are made in Python (and older versions have some components also made in C++) and are meant to be run on images first. These programs will callibrate images, detect stars, and collate star data per image into organized CSV files. The image analysis programs are made in R and can do a variety of different analysis on individual stars and on the entire field as a whole. All the programs are broken down more through the documentation.

## Image Analysis Programs

### TOROSphotoCalibration
This program will calibrate raw TOROS images. It performs a scaled bias+dark subtraction, clips the overscan regions, performs flat-fielding, and does a background subtraction using the _SExtractor_ function from photutils. \
\
**Inputs** \
\
*raw_file:* String, path to raw image. \
*bias_file:* String, path to bias image. \
*dark_file:* String, path to dark image. \
*flat_file:* String, path to flat image. \
*science_file:* String, path to calibrated image.
*write_bkg:* Boolean, if true the program will write the generated background image to disk. Set to False by default. \
*bkg_file:* String, path to background image. Set to "" by default. \
*write_nobkg:* Boolean, if true the program will write the science image before doing background subtraction. Set to False by default. \
*nobkg_file:* String, path to pre-background subtracted image. Set to "" by default. \
*use_mask:* Boolean, if true the program will put a circular mask on image before doing background subtraction. Set to False by default. \
*cx:* Float, the x-center of the mask. Set to 0 by default. \
*cy:* Float, the y-center of the mask. Set to 0 by default. \
*r:* Float, the radius of the mask. Set to 0 by default. \


