# TOROS-Analysis-Pipeline
Full variability analysis pipeline for the TOROS telescope in Argentina. These are updated versions of the previous programs I posted before; assume any TOROS programs I make that are not within this repo are depreciated.

## Structure
Programs have been divided into two categories: image reduction programs and image analysis programs. The image reduction programs are made in Python (and older versions have some components also made in C++) and are meant to be run on images first. These programs will callibrate images, detect stars, and collate star data per image into organized CSV files. The image analysis programs are made in R and can do a variety of different analysis on individual stars and on the entire field as a whole. All the programs are broken down more through the documentation.

## Image Reduction Programs

### TOROSphotoCalibration
This program will calibrate raw TOROS images. It performs a scaled bias+dark subtraction, clips the overscan regions, performs flat-fielding, and does a background subtraction using the _SExtractor_ function from photutils. \
\
**INPUTS** \
*raw_file*: String, path to raw image. \
*bias_file*: String, path to bias image. \
*dark_file*: String, path to dark image. \
*flat_file*: String, path to flat image. \
*science_file*: String, path to calibrated image.
*write_bkg*: Boolean, if true the program will write the generated background image to disk. Set to False by default. \
*bkg_file*: String, path to background image. Set to "" by default. \
*write_nobkg*: Boolean, if true the program will write the science image before doing background subtraction. Set to False by default. \
*nobkg_file*: String, path to pre-background subtracted image. Set to "" by default. \
*use_mask*: Boolean, if true the program will put a circular mask on image before doing background subtraction. Set to False by default. \
*cx*: Float, the x-center of the mask. Set to 0 by default. \
*cy*: Float, the y-center of the mask. Set to 0 by default. \
*r*: Float, the radius of the mask. Set to 0 by default.
\
\
**OUTPUTS** \
Output is a tuple. \
*[0]*: Astropy table, contains the .fits data of the science image. \
*[1]*: fits header, the .fits header file for the science image.

---

### TOROSaperturePhotometry
This program does many things at once. First, the science image is loaded in and a preliminary star detection is done with a high sigma. Using these first preliminary detections, a FWHM for the image is fit and the true star detection is then run (both star detections use the DAOStarFinder algorithm). Next, to avoid multiple centroids being assigned to the same star due to photometric defects, a DBSCAN is run which mergers very nearby centroids together. In order to avoid effects from vignetting, any data too near the edges of the images are cut off (controlled by the edge_buffer parameter). Because TOROS does not natively save its WCS data, it must be created manually. This is done by assuming the user has astrometry.net installed natively on their device. For this program to run correctly, you must change the WSL paths so that they reflect the locations of you "solve-field" file and the matching field files. Once the WCS object is made, the actual photometry is performed on the detected centroids. The photometry is done with both an aperture and annuli to do local background subtraction. Once all statistics are made, a table is build and returned. This table can also be written to disk. \
\
**INPUTS** \
*science_file*: String, path to science image. \
*temp_dir*: String, path to folder where you would like temporary files to be stored. \
*starList_file*: String, path to folder where you want tables to be saved (if that is desired). Set to "" by default. \
*write_table*: Boolean, if true the program will write the generated statistics table to disk. Set to False by default. \
*edge_buffer*: Float, this is the amount of pixels that will be cut off from all edges of the image to avoid vignetting. Set to 200 by default. \
*use_mask*: Boolean, if true the program will put a circular mask on the image before doing star detection. Set to False by default. \
*cx*: Float, the x-center of the mask. Set to 0 by default. \
*cy*: Float, the y-center of the mask. Set to 0 by default. \
*r*: Float, the radius of the mask. Set to 0 by default. 
\
\
**OUTPUTS** \
*phot_table*: Astropy Table, contains the photometry data for all stars detected in the field.

---


### TOROSphotometryAlign
This program will ensure that StarIDs for all images reference the same star/centroid. This is done by matching the RA and DECs of stars from a "master" table to stars from a "comparator" table and changing the matched "comparator" StarIDs to the "master" StarID. \
\
**INPUTS** \
*master*: Astropy Table, contains StarIDs that will be retained. \
*comparator*: Astropy Table, contains StarIDs that will be changed. \
*aligned_file*: String, path to the final photometry table that will have modified StarIDs. \
*snr_threshold*: Integer, the SNR cutoff when filtering stars before alignment. Set to 10 by default. 
\
\
**OUTPUTS** \
*aligned_tbl*: Astropy Table, contains the photometry data for the comparator field with the new master ids appended.


## Image Analysis Programs

### TOROSloadInNights
This program will load in individual folders containing nights worth of images into an R list. The R list will contain the photometric data for each night as a dataframe, a list of all the caught StarIDs, and a list of all StarIDs present in every image. \
\
**INPUTS** \
*night_folders*: Vector, contains the paths to folders. These folders should contain the tables created from **TOROSphotometryAlign** that you want to analyze.
\
\
**OUTPUTS** \
Output is a list. \
*$CommonIDs*: Vector, contains StarIDs detected among all fields.
*$AllIDs*: Vector, contains every detected StarID.
*$Photometry*: List, contains the photometry data of all stars detected.

___

### TOROSmakeStarStatistics
This program will seperate each night's dataframe into individual dataframe of information for every individual detected star. This program will also correct for any night offsets in magnitude that may be present and will perform a median-built systematic removal. The systematic removal runs in parallel to prevent the program from taking an absurd time to run for large data. Finally, this program returns a list containing the night offset data and a list of dataframes with each star's lightcurve. \
\
**INPUTS** \
*data_list:* List, contains the photometry data from the **$Photometry** component of the **TOROSloadInNights** output.
\
\
**OUTPUTS** \
Output is a list. \
*$LCs*: List, contains all of the adjusted photometric data for each star separated into individual dataframes.
*$NightOffsets*: Vector, contains the magnitude offsets for each night in chronological order.

---


### TOROSshowLightCurve
This program will create an unphase-folded light curve of a specified star. It will also plot systematic data if available.
\
**INPUTS** \
*lcs*: List, contains the photometry data from the **$LCs** component of the **TOROSmakeStarStatistics** output. \
*lc_targetID*: Integer, the StarID of the star you want the light curve of. \
*observed_field*: String, the name of the field you are observing (purely for cosmetic purposes). Set to "Observed Field" by default. \
*use_nameCustom*: Boolean, if true the program will use a custom name for the star. Set to F by default. \
*nameCustom*: String, the name of the star being plotted (purely for cosmetic purposes). Set to "" by default. 
\
\
**OUTPUTS** \
None

---

### TOROSshowPhaseFoldCurve
This program will create a phase-folded light curve of a specified star.
\
**INPUTS** \
*lcs*: List, contains the photometry data from the **$LCs** component of the **TOROSmakeStarStatistics** output. \
*lc_targetID*: Integer, the StarID of the star you want the light curve of. \
*period*: Float, the period of the star being plotted; used to fold the curve. \
*observed_field*: String, the name of the field you are observing (purely for cosmetic purposes). Set to "Observed Field" by default. \
*use_nameCustom*: Boolean, if true the program will use a custom name for the star. Set to F by default. \
*nameCustom*: String, the name of the star being plotted (purely for cosmetic purposes). Set to "" by default. \
*bin*: Boolean, if true the program will bin the data points by taking the median of each night. Set to F by default.
\
\
**OUTPUTS** \
None

---

### TOROSlomgScarg
This program will create an LSP of a specified star using the lomb package in R.
\
**INPUTS** \
*lcs*: List, contains the photometry data from the **$LCs** component of the **TOROSmakeStarStatistics** output. \
*lsp_targetID*: Integer, the StarID of the star you want the Lomb-Scargle Periodogram of. \
*observed_field*: String, the name of the field you are observing (purely for cosmetic purposes). Set to "Observed Field" by default. \
*use_nameCustom*: Boolean, if true the program will use a custom name for the star. Set to F by default. \
*nameCustom*: String, the name of the star being plotted (purely for cosmetic purposes). Set to "" by default. \
*phaseFold*: Boolean, if true the program will use the periods from the LSP to create phase folded light curves of the target star. Set to F by default. \
*phaseFold_num*: Integer, the number of LSP periods you want to be used in phase folding (ie if this number =3, then phase folded light curves will be made using the 1st, 2nd, and 3rd most significant periods). Set to 0 by default. \
\
**OUTPUTS** \
*lomb*: List, contains all data from the lsp function of the *lomb* package. \

---

### TOROSvariabilityAndRMS
This program will do RMS and variability analysis of the stars found in the field.
\
**INPUTS** \
*lcs*: List, contains the photometry data from the **$LCs** component of the **TOROSmakeStarStatistics** output. \
*observed_field*: String, the name of the field you are observing (purely for cosmetic purposes). Set to "Observed Field" by default. \
\
**OUTPUTS** \
Output is a list. \
*$Variability*: Dataframe, contains the variability analysis data of the stars in the field. \
*$RMSData*: Dataframe, contains the RMS analysis data of the stars in the field. \

---

### TOROSobsTimeLSSTpredict
This program will perform a sampling simulation to determine the optimal number of observations and time intervals of observations to recover the variability of a star. \
**INPUTS** \
*lcs*: List, contains the photometry data from the **$LCs** component of the **TOROSmakeStarStatistics** output. \
*test_targetID*: Integer, the StarID of the star you want to test. \
*sigma_level*: Integer, the sigma level above which a star is considered variable. Set to 2 by default.\
*N*: Integer, the number of interations for the observation test. Set to 100 by default. \
*M*: Integer, the number of interations for the interval test. Set to 10000 by default. \
*use_nameCustom*: Boolean, if true the program will use a custom name for the star. Set to F by default. \
*nameCustom*: String, the name of the star being plotted (purely for cosmetic purposes). Set to "" by default. \
*observed_field*: String, the name of the field you are observing (purely for cosmetic purposes). Set to "Observed Field" by default. \
\
**OUTPUTS** \
Output is a list. \
*$LSSTPred1*: List, contains the observation number test data. \
*$LSSTPred2*: List, contains the time interval test data. \

## Misc

### TOROSsimbadPositionList
This program will create a .txt file containing the RA and DEC values of the flagged stars in your field, ordered by brightness. \
**INPUTS** \
*lcs*: List, contains the photometry data from the **$LCs** component of the **TOROSmakeStarStatistics** output. \
*var_df*: Dataframe, contains the variability data from the **$Variability** component of the **TOROSvariabilityAndRMS** output. \
\
**OUTPUTS** \
None










