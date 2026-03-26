def TOROSphotoCalibration(raw_file , bias_file , flat_file , dark_file , science_file ,
                          write_bkg = False , bkg_file = "" ,
                          write_nobkg = False , nobkg_file = "" ,
                          use_mask = False , cx = 0.0 , cy = 0.0 , r = 0.0):

    """
    This program will reduce raw TOROS images by performing flat fielding, bias subtraction, dark subtraction, and a background subtraction.
    :param raw_file: String, file path to raw TOROS image
    :param bias_file: String, file path to bias image
    :param flat_file: String, file path to flat image
    :param dark_file: String, file path to dark image
    :param science_file: String, file path to science image
    :param write_bkg: Boolean, whether to write background image
    :param bkg_file: String, file path to background image
    :param write_nobkg: Boolean, whether to write image pre-background subtraction
    :param nobkg_file: String, file path to pre-background subtraction image
    :param use_mask: Boolean, whether to mask image during background subtraction
    :param cx: Float, x center of mask
    :param cy: Float, y center of mask
    :param r: Float, radius of mask
    :return: (science , science_header)
    """


    ###Required packages
    import os
    from astropy.io import fits
    from photutils.background import Background2D , SExtractorBackground
    from astropy.table import Table
    import numpy as np
    from astropy.stats import sigma_clip , SigmaClip


    ###Subfunctions (written by Dr. Ryan Oelkers)
    def clip_image(img, header):
        """ This function will clip the image and remove any overscan regions. This function is written for TOROS
        specifically, and you will need to update it for any given CCD.

        :parameter - image - The image to clip
        :parameter - header - the header of the image

        :return image_clip, header - The clipped image and the new header
        """

        # make the clipped image
        image_clip = np.zeros((10560, 10560))

        # what are the sizes of the over scan?
        ovs_x_size = 180
        ovs_y_size = 40

        # what are the sizes of each chip (not including the over scan)?
        chip_x_size = 1320
        chip_y_size = 5280

        # what is the full size of the chip (including over scan)
        full_chip_x = chip_x_size + ovs_x_size
        full_chip_y = chip_y_size + ovs_y_size

        # move through x and y
        idx = 0
        for x in range(0, 12000, full_chip_x):

            idy = 0
            for y in range(0, 10600, full_chip_y):
                # put the clipped image into the holder image
                image_clip[idy:idy + chip_y_size, idx:idx + chip_x_size] = img[y:y + chip_y_size, x:x + chip_x_size]

                # increase the size of the yclip
                idy = idy + chip_y_size

            # increase the size of the xclip
            idx = idx + chip_x_size

        # update the header
        header['OVERSCAN'] = 'removed'
        header['X_CLIP'] = ovs_x_size
        header['Y_CLIP'] = ovs_y_size

        return image_clip, header

    def subtract_scaled_bias_dark(img, header, bias_file, dark_file):
        """ This function will scale the bias frame to match the overscan and then remove the dark level.

        :parameter - img - The image to remove the bias and dark frame from
        :parameter - header - The header of the image

        :return img_bias_dark, header
        """

        # read in the bias frame and dark frame
        bias = fits.getdata(bias_file)
        dark = fits.getdata(dark_file)

        # what are the sizes of the over scan?
        ovs_x_size = 180
        ovs_y_size = 20

        # what are the sizes of each chip (not including the over scan)?
        chip_x_size = 1320
        chip_y_size = 5280

        # what is the full size of the chip (including over scan)
        full_chip_x = chip_x_size + ovs_x_size
        full_chip_y = chip_y_size + ovs_y_size

        # make a copy of the bias frame
        bias_scl = bias.copy()

        # move through x and y to mask the "image" parts of hte image
        for x in range(0, img.shape[1], full_chip_x):  # img.shape[1] gives x size
            for y in range(0, img.shape[0], full_chip_y):

                if y == 0:
                    # pull out the overscan from the raw image and set the "image" to 0
                    img_slice = img[y:y + full_chip_y, x:x + full_chip_x].copy()
                    img_slice[0:chip_y_size, 0:chip_x_size] = 0
                else:
                    # pull out the overscan from the raw image and set the "image" to 0
                    img_slice = img[y:y + full_chip_y, x:x + full_chip_x].copy()
                    img_slice[ovs_y_size:ovs_y_size + chip_y_size, 0:chip_x_size] = 0

                # put the clipped image into the holder image
                bias_scl[y:y + full_chip_y, x:x + full_chip_x] = (bias_scl[y:y + full_chip_y, x:x + full_chip_x] +
                                                                  np.median(img_slice[img_slice > 0]))

        # now remove the bias from the image
        img_bias = img - bias_scl

        # now remove the dark current from the image
        img_bias_dark = img_bias - dark

        # now update the image header
        # update the header
        header['BIAS_SUBT'] = 'Y'
        header['BIAS_TYPE'] = 'SCALE'
        header['DARK_SUBT'] = 'Y'

        return img_bias_dark, header


    ###Check if files exist
    if not(os.path.exists(raw_file) and os.path.exists(bias_file) and os.path.exists(dark_file) and os.path.exists(flat_file)):
        print("Missing file(s)!")
        return -1


    ###Load in raw and flat
    (raw , raw_header) = fits.getdata(filename = raw_file , header = True)
    (flat , flat_header) = fits.getdata(filename = flat_file , header = True)


    ###Perform dark+bias subtraction using subtract_scaled_bias_dark()
    (bdsub , bdsub_header) = subtract_scaled_bias_dark(img = raw , header = raw_header , bias_file = bias_file , dark_file = dark_file)


    ###Peform flat fielding
    #First, clip the overscans from the bdsub and the flat
    (bdsub_clipped , bdsub_clipped_header) = clip_image(img = bdsub , header = bdsub_header)
    (flat_clipped , flat_clipped_header) = clip_image(img = flat , header = flat_header )

    #Normalize flat
    flat_clipped_norm = flat_clipped / np.median(flat_clipped)

    #Now do the flat fielding using the clipped images
    bdf = bdsub_clipped / flat_clipped_norm
    bdf_header = bdsub_clipped_header.copy()
    bdf_header['CAL'] = "bias_subbed , dark_subbed , flat_fielded"

    #Remove bad values if they occur
    bdf = np.nan_to_num(bdf , nan = 0.0 , posinf = 0.0 , neginf = 0.0)


    ###If desired, the image before background subtraction can be saved
    if write_nobkg:

        hdu = fits.PrimaryHDU(data = bdf , header = bdf_header)
        hdu.writeto(nobkg_file, overwrite = True)


    ###Now perform background subtraction on image
    #Create background removal objects and set parameters
    sigma_clipper = SigmaClip(sigma = 3.0)
    bkg_estimator = SExtractorBackground(sigma_clip = sigma_clipper)

    #If mask is enabled, mask image before making background
    if use_mask:

        mask = np.zeros(bdf.shape , dtype = bool)
        (yy , xx) = np.indices(bdf.shape)
        mask_arg = (xx - cx)**2 + (yy - cy)**2 < r**2
        mask[mask_arg] = True

        bkg = Background2D(bdf , box_size = (40 , 40) ,
                           filter_size = (3 , 3) ,
                           sigma_clip = sigma_clipper ,
                           bkg_estimator = bkg_estimator ,
                           mask = mask)

    else:

        bkg = Background2D(bdf , box_size = (40 , 40) ,
                           filter_size = (3 , 3) ,
                           sigma_clip = sigma_clipper ,
                           bkg_estimator = bkg_estimator)

    #Subtract background
    science = bdf - bkg.background
    science_header = bdf_header.copy()

    #Update header accordingly
    if use_mask:

        science_header['CAL'] = "bias_subbed , dark_subbed , flat_fielded , bkg_sub with mask"

    else:

        science_header['CAL'] = "bias_subbed , dark_subbed , flat_fielded , bkg_sub without mask"


    ###Write science file to disk
    hdu = fits.PrimaryHDU(data = science, header = science_header)
    hdu.writeto(science_file , overwrite = True)


    ###Save the background if desired
    if write_bkg:
        hdu = fits.PrimaryHDU(data = bkg.background)
        hdu.writeto(bkg_file, overwrite = True)


    ###End function
    return (science , science_header)
