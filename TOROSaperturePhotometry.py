def TOROSaperturePhotometry(science_file , temp_dir , starList_file = "" ,
                            write_table = False ,
                            edge_buffer = 200 ,
                            use_mask = False ,
                            cx = 0.0 , cy = 0.0 , r = 0.0):
    """
    This function will perform aperture photometry on a TOROS science image and save the information as a csv file.
    :param science_file: String, file path to science image
    :param starList_file: String, file path to star list (.csv file)
    :param write_table: Boolean, if True will write the aperture photometry information to a table
    :param edge_buffer: Int, buffer filter so that stars near edge are not analyzed
    :param use_mask: Bool, whether to mask image during star finding
    :param cx: Float, x center of mask
    :param cy: Float, y center of mask
    :param r: Float, radius of mask
    :return: phot_table
    """
    ###Required packages
    import os
    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats , SigmaClip
    from photutils.detection import DAOStarFinder
    import numpy as np
    from photutils.aperture import CircularAperture , aperture_photometry , CircularAnnulus , ApertureStats
    from astropy.table import Table , Column
    from astropy.wcs import WCS
    from sklearn.cluster import DBSCAN
    import subprocess
    from photutils.psf import fit_fwhm
    import glob
    import gc

    ###Check if files exist
    if not(os.path.exists(science_file)):

        print("Missing science file!")
        return -1


    ###Load in science image
    (science , science_header) = fits.getdata(filename = science_file , header = True , memmap = True)


    ###Obtain stddev for star_finder
    (mean , median , std) = sigma_clipped_stats(data = science , sigma = 3.0)
    starFind = DAOStarFinder(fwhm = 3.0 , threshold = 8 * std)


    ###We now will do a preliminary source detection to fit a fwhm
    if use_mask:

        mask = np.zeros(science.shape , dtype = bool)
        (yy , xx) = np.indices(science.shape)
        mask_arg = (xx - cx)**2 + (yy - cy)**2 < r**2
        mask[mask_arg] = True

        print("Finding stars [prelim]...")
        sources_rough = starFind(science , mask = mask)

    else:

        print("Finding stars [prelim]...")
        sources_rough = starFind(science)

    #Make sure sources_rough has values
    if sources_rough is None or len(sources_rough) == 0:
        print("No stars found in initial rough detection...")
        return -1

    #Pick good stars for the fwhm fitting
    idx = np.argsort(sources_rough['flux'])[::-1]
    candidates = sources_rough[idx][20:min(200 , len(idx))]

    #Fit fwhm
    fwhm_list = []

    for candidate in candidates:

        x = int(candidate['xcentroid'])
        y = int(candidate['ycentroid'])


        if x < 16 or y < 16 or x > science.shape[1]-16 or y > science.shape[0]-16:
            continue #Reject candiates too close to the edge

        try:
            cutout = science[y - 15 : y + 16, x - 15 : x + 16]
            fwhm_current = fit_fwhm(data = cutout)
            if np.isfinite(fwhm_current):
                fwhm_list.append(fwhm_current)

        except:
            continue

    #Check if fwhm_list has actual values, if so use median
    if len(fwhm_list) == 0:
        print("FWHM fitting failed, using fallback of FWHM = 3.0")
        fwhm = 3.0

    else:
        fwhm = np.median(fwhm_list)


    ###Use global_fwhm for starfinding now
    starFind = DAOStarFinder(fwhm = fwhm , threshold = 5 * std)
    if use_mask:

        mask = np.zeros(science.shape , dtype = bool)
        (yy , xx) = np.indices(science.shape)
        mask_arg = (xx - cx)**2 + (yy - cy)**2 < r**2
        mask[mask_arg] = True

        print("Finding stars...")
        sources = starFind(science , mask = mask)

    else:

        print("Finding stars...")
        sources = starFind(science)

    #Make sure sources has values
    if sources is None or len(sources) == 0:
        print("No stars found in detection...")
        return -1


    ###To avoid multiple centroids associated to the same star, nearby centroids will be merged using DBSCAN
    print("Merging...")
    source_centroids = np.vstack([sources['xcentroid'] , sources['ycentroid']]).T

    clustering = DBSCAN(eps = .7 * fwhm , min_samples = 1).fit(source_centroids)
    labels = clustering.labels_

    #Collate merged sources
    merged_ids = []

    for label in np.unique(labels):

        group = sources[labels == label]

        #Get mean of centroids labeling same star
        x = np.mean(group['xcentroid'])
        y = np.mean(group['ycentroid'])

        #Choose maximum flux from group
        flux = np.max(group['flux'])

        #Append to merged_ids
        merged_ids.append((x , y , flux))

    #Create lists of merged x, y, and flux
    if len(merged_ids) == 0:
        print("DBSCAN fail.")
        return -1

    (merged_x , merged_y , merged_flux) = zip(*merged_ids)
    merged_sources = Table()
    merged_sources['xcentroid'] = merged_x
    merged_sources['ycentroid'] = merged_y
    merged_sources['flux'] = merged_flux

    #Ignore stars below specified SNR
    f = merged_sources['flux']
    snr_est = f / np.sqrt(np.abs(f)) #Assume poisson
    mask = snr_est > 65
    merged_sources = merged_sources[mask]
    print(merged_sources)

    ###Ignore sources on fringes of image to avoid vignetting effects
    (h , w) = science.shape
    merged_sources = merged_sources[
        (merged_sources['xcentroid'] > edge_buffer ) &
        (merged_sources['xcentroid'] < w - edge_buffer ) &
        (merged_sources['ycentroid'] > edge_buffer ) &
        (merged_sources['ycentroid'] < h - edge_buffer )
    ]


    ###TOROS does not save WCS data so it will have to be made natively. This program assumes the user has astrometry.net natively installed.
    #Use the brightest 500 stars to construct a list of sources
    N = 500
    bright_idx = np.argsort(merged_sources['flux'])[::-1][:N]
    astrometry_sources = merged_sources[bright_idx]

    #Convert astrometric sources into a table of just x and y
    astrometry_table = Table()
    astrometry_table['x'] = astrometry_sources['xcentroid']
    astrometry_table['y'] = astrometry_sources['ycentroid']

    #Save this file to disk
    base = os.path.splitext(os.path.basename(science_file))[0]
    astrometry_file = os.path.join(temp_dir , base + "_astrometry.fits")
    astrometry_table.write(astrometry_file , overwrite = True)

    #Convert Windows path to WSL path (since astrometry.net runs natively on Linux)
    wsl_astrometry_file = astrometry_file.replace("E:\\" , "/mnt/e/").replace("\\" , "/") #NOTE: one must change this to reflect where they are saving their science files and where their WSL is mounted, it will change per computer

    #Run the astrometric matching
    (h , w) = science.shape

    print("Creating WCS...")
    try:
        subprocess.run(
            ['wsl' , '/usr/bin/solve-field' , wsl_astrometry_file ,
             '--x-column' , 'x' , '--y-column' , 'y' ,
             '--width' , str(w) , '--height' , str(h) ,
            '--overwrite' , '--no-plots' ] ,
            check = True)

    except subprocess.CalledProcessError:
        print("Astrometry failed. Womp Womp.")
        return -1


    ###Also save date and time of event
    date = str(science_header["DATE"].split("T")[0])
    time = science_header["DATE"].split("T")[1]

    #Make lists of date and time with same length as sources
    size = len(merged_sources)

    date_list = []
    time_list = []

    for i in range(size):

        date_list.append(str(date))
        time_list.append(str(time))

    #Put columns with date and time for each star
    date_column = Column(date_list , name = "Date")
    time_column = Column(time_list , name = "Time")


    ###Perfom aperture photometry
    #Define positions, aperture, and annuli
    print("Performing photometry...")
    positions = np.transpose((merged_sources['xcentroid'] , merged_sources['ycentroid']))

    r_ap = 2 * fwhm
    r_in = 4 * fwhm
    r_out = 6 * fwhm

    apertures = CircularAperture(positions , r = r_ap)
    aperture_annuli = CircularAnnulus(positions , r_in = r_in , r_out = r_out)

    #Create statistics of local backgrounds
    sigclip = SigmaClip(sigma = 3.0 , maxiters = 10)
    locBkg_stats = ApertureStats(science , aperture_annuli , sigma_clip = sigclip)
    locBkg_med = locBkg_stats.median
    locBkg_std = locBkg_stats.std

    #Create statistics for apertures
    aperture_stats = ApertureStats(science , apertures , sum_method = "center" , local_bkg = locBkg_med)

    #Calculate a more robust SNR assuming Poisson noise and other errors
    ap_flux = aperture_stats.sum
    ap_area = aperture_stats.sum_aper_area.value
    flux_err = np.sqrt(np.abs(ap_flux) + ap_area * (locBkg_std**2))
    snr = ap_flux / flux_err


    ###Load in WCS data created previously
    wcs_file = astrometry_file.replace(".fits" , ".wcs")
    if not os.path.exists(wcs_file):
        print("WCS file not found!")
        return -1

    wcs_header = fits.getheader(wcs_file)
    telescope = WCS(wcs_header)

    sky_apertures = apertures.to_sky(telescope)
    sky_positions = sky_apertures.positions #RA and DEC

    RAcolumn = Column(sky_positions.ra.deg , name = "RA")
    DEcolumn = Column(sky_positions.dec.deg , name = "DEC")


    ###Clean up temporary astrometry files
    print("Cleaning...")
    base = astrometry_file.replace(".fits" , "")
    for extension in [".axy" , ".corr" , ".match" , ".rdls" , ".solved" , ".wcs" , ".xyls"]:
        f = base + extension
        if os.path.exists(f):
            os.remove(f)


    ###Put all data into an astropy table and garbage collect
    print("Making table...")

    phot_table = Table({
        "id" : aperture_stats.id ,
        "xcenter" : aperture_stats.xcentroid ,
        "ycenter" : aperture_stats.ycentroid ,
        "flux" : aperture_stats.sum ,
        "fluxerr" : flux_err ,
        "SNR" : snr ,
        "area" : aperture_stats.sum_aper_area ,
        "bkg_median" : locBkg_med ,
        "bkg_std" : locBkg_std ,
    })

    phot_table.add_column(date_column , name = "Date")
    phot_table.add_column(time_column , name = "Time")
    phot_table.add_column(RAcolumn , name = "RA")
    phot_table.add_column(DEcolumn , name = "DEC")

    if write_table:

        phot_table.write(starList_file , delimiter = ',' , format = "csv" , overwrite = True)

    del science
    del science_header
    del sources
    del merged_sources
    del apertures
    del aperture_annuli
    del aperture_stats
    del locBkg_stats

    gc.collect()
    return phot_table
