# TOROS-Analysis-Pipeline
Full variability analysis pipeline for the TOROS telescope in Argentina. These are updated versions of the previous programs I posted before; assume any TOROS programs I make that are not within this repo are depreciated.

## Structure
Programs have been divided into two categories: image reduction programs and image analysis programs. The image reduction programs are made in Python (and older versions have some components also made in C++) and are meant to be run on images first. These programs will callibrate images, detect stars, and collate star data per image into organized CSV files. The image analysis programs are made in R and can do a variety of different analysis on individual stars and on the entire field as a whole. All the programs are broken down more through the documentation.

## Image Reduction Programs

### TOROSphotoCalibration
This program will calibrate raw TOROS images. It performs a scaled bias+dark subtraction, clips the overscan regions, performs flat-fielding, and does a background subtraction using the _SExtractor_ function from photutils. \
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
*r:* Float, the radius of the mask. Set to 0 by default.


### TOROSaperturePhotometry
This program does many things at once. First, the science image is loaded in and a preliminary star detection is done with a high sigma. Using these first preliminary detections, a FWHM for the image is fit and the true star detection is then run (both star detections use the DAOStarFinder algorithm). Next, to avoid multiple centroids being assigned to the same star due to photometric defects, a DBSCAN is run which mergers very nearby centroids together. In order to avoid effects from vignetting, any data too near the edges of the images are cut off (controlled by the edge_buffer parameter). Because TOROS does not natively save its WCS data, it must be created manually. This is done by assuming the user has astrometry.net installed natively on their device. For this program to run correctly, you must change the WSL paths so that they reflect the locations of you "solve-field" file and the matching field files. Once the WCS object is made, the actual photometry is performed on the detected centroids. The photometry is done with both an aperture and annuli to do local background subtraction. Once all statistics are made, a table is build and returned. This table can also be written to disk. \
\
*science_file:* String, path to science image. \
*temp_dir:* String, path to folder where you would like temporary files to be stored. \
*starList_file:* String, path to folder where you want tables to be saved (if that is desired). Set to "" by default. \
*write_table:* Boolean, if true the program will write the generated statistics table to disk. Set to False by default. \
*edge_buffer:* Float, this is the amount of pixels that will be cut off from all edges of the image to avoid vignetting. Set to 200 by default. \
*use_mask:* Boolean, if true the program will put a circular mask on the image before doing star detection. Set to False by default. \
*cx:* Float, the x-center of the mask. Set to 0 by default. \
*cy:* Float, the y-center of the mask. Set to 0 by default. \
*r:* Float, the radius of the mask. Set to 0 by default. 


### TOROSphotometryAlign
This program will ensure that StarIDs for all images reference the same star/centroid. This is done by matching the RA and DECs of stars from a "master" table to stars from a "comparator" table and changing the matched "comparator" StarIDs to the "master" StarID. \
\
*master:* Astropy Table, contains StarIDs that will be retained. \
*comparator:* Astropy Table, contains StarIDs that will be changed. \
*aligned_file:* String, path to the final photometry table that will have modified StarIDs. \
*snr_threshold:* Integer, the SNR cutoff when filtering stars before alignment. Set to 10 by default. 


## Image Analysis Programs




