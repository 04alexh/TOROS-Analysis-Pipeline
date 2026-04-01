# TOROS-Analysis-Pipeline
Full variability analysis pipeline for the TOROS telescope in Argentina. These are updated versions of the previous programs I posted before; assume any TOROS programs I make that are not within this repo are depreciated.

## Structure
Programs have been divided into two categories: image reduction programs and image analysis programs. The image reduction programs are made in Python (and older versions have some components also made in C++) and are meant to be run on images first. These programs will callibrate images, detect stars, and collate star data per image into organized CSV files. The image analysis programs are made in R and can do a variety of different analysis on individual stars and on the entire field as a whole. All the programs are broken down more through the documentation.




