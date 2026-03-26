def TOROSphotometryAlign(master , comparator , aligned_file ,
                         snr_threshold = 100):
    """

    :param master_file: String, file path to star list you want as basis for star ids
    :param comparator_file: String, file path to star list you want to match to master_file star ids
    :param aligned_file: String, file path to aligned star list (.csv)
    :param snr_threshold: Integer, signal-to-noise filter for stars
    :return: aligned_tbl
    """

    ###Required packages
    import os
    import numpy as np
    from astropy.table import Table
    from astropy.coordinates import SkyCoord
    import astropy.units as u


    ###Read the .csv star lists
    master_tbl = master
    comparator_tbl = comparator


    ###Filter out data to remove any remaining bad stars
    #Remove negative flux caused by background subtraction
    master_tbl = master_tbl[master_tbl['flux'] >= 0]
    comparator_tbl = comparator_tbl[comparator_tbl['flux'] >= 0]
    print(master_tbl , "MASTER FIXED")
    print(comparator_tbl , "COMPARATOR FIXED")

    #Filter out fluxes that do not meet the SNR requirement
    master_tbl = master_tbl[master_tbl['SNR'] >= snr_threshold]
    comparator_tbl = comparator_tbl[comparator_tbl['SNR'] >= snr_threshold]
    print(master_tbl , "MASTER FILTERED")
    print(comparator_tbl , "COMPARATOR FILTERED")

    #Sort tables from highest flux to lowest
    master_tbl = master_tbl[np.argsort(master_tbl['flux'])[::-1]]
    comparator_tbl = comparator_tbl[np.argsort(comparator_tbl['flux'])[::-1]]
    print(master_tbl , "MASTER SORTED")
    print(comparator_tbl , "COMPARATOR SORTED")


    ###Obtain RA and DEC from master and comparator
    master_RA = master_tbl['RA']
    master_DEC = master_tbl['DEC']

    comparator_RA = comparator_tbl['RA']
    comparator_DEC = comparator_tbl['DEC']

    eqc_master = SkyCoord(master_RA * u.deg , master_DEC * u.deg)
    eqc_comparator = SkyCoord(comparator_RA * u.deg , comparator_DEC * u.deg)


    ###Match starIDs between comparator to master using RA and DEC
    #Find indexes in comparator that correspond to master stars
    (idx , d2d , _) = eqc_master.match_to_catalog_sky(eqc_comparator)
    good_match = d2d < 1.0 * u.arcsec

    #Prevent duplicate assignments and make matching more robust
    (unique_idx , unique_mask) = np.unique(idx[good_match] , return_index = True)

    #Create matched lists
    master_matched = master_tbl[good_match][unique_mask]
    comparator_matched = comparator_tbl[idx[good_match]][unique_mask] #Now master_matched[i] = comparator_matched[i]

    if len(master_matched) == 0:
        print("No matches found.")
        return -1


    ###Build new astropy table
    aligned_tbl = Table()

    aligned_tbl['Aligned_ID'] = master_matched['id']
    aligned_tbl['Old_ID'] = comparator_matched['id']
    aligned_tbl['x'] = comparator_matched['xcenter']
    aligned_tbl['y'] = comparator_matched['ycenter']
    aligned_tbl['RA'] = comparator_matched['RA']
    aligned_tbl['DEC'] = comparator_matched['DEC']
    aligned_tbl['Flux'] = comparator_matched['flux']
    aligned_tbl['FluxError'] = comparator_matched['fluxerr']
    aligned_tbl['SNR'] = comparator_matched['SNR']
    aligned_tbl['PixArea'] = comparator_matched['area']
    aligned_tbl['LocBkgMed'] = comparator_matched['bkg_median']
    aligned_tbl['LocBkgStd'] = comparator_matched['bkg_std']
    aligned_tbl['Date'] = comparator_matched['Date']
    aligned_tbl['Time'] = comparator_matched['Time']

    #Sort by id
    aligned_tbl.sort('Aligned_ID')

    #Write table
    aligned_tbl.write(aligned_file , delimiter = ',' , format = "csv" , overwrite = True)
    return aligned_tbl
