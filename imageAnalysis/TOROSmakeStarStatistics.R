TOROSmakeStarStatistics <- function(data_list , id_list)
{

  ###Load in required packages
  library(dplyr)


  ###Combine all dataframe from list into one
  all_data <- do.call(rbind , data_list)
  all_data$Aligned_ID <- as.character(all_data$Aligned_ID)


  ###Clean and compute required columns
  #Remove bad fluxes
  all_data <- all_data[all_data$Flux > 0 , ]

  #Create time stamps
  all_data$Date <- as.Date(all_data$Date)
  all_data$TimeStamp <- as.POSIXct(paste(all_data$Date , all_data$Time) ,
               format = "%Y-%m-%d %H:%M:%OS" ,
               tz = "UTC"
               )

  #Create instrumental magnitude column
  all_data$Mag <- 25 - 2.5*log10(all_data$Flux)


  ###Each night may have a median magnitude offset that must be corrected
  #Obtain median magnitude per star per night
  night_medians <- all_data %>%
    group_by(Date , Aligned_ID) %>% #group all_data by Aligned_ID and Date
    summarize(star_med = median(Mag) , .groups = "drop") %>% #Calculate the median magnitudes of all the groups
    group_by(Date) %>% #now group by just date
    summarize(night_med = median(star_med , na.rm = T) , .groups = "drop") #take the median of the median magnitudes per star for all dates

  #Calculate offsets relative to first night
  baseline <- night_medians$night_med[1]
  night_medians$offset <- baseline - night_medians$night_med

  #Map the offsets back to the dates
  offset_map <- setNames(night_medians$offset , night_medians$Date)
  all_data$Mag <- all_data$Mag + offset_map[as.character(all_data$Date)]


  #Next systematics will be removed using a median built approach
  print("systematics")

  process_star <<- function(i , star_summary , lcs)
  {

    print("getting target")
    target <- star_summary[i , ]
    target_id <- target$StarID

    #Get deltaRA and deltaDEC of target from all stars
    print("calcing deltas")
    deltaRA <- (star_summary$ra - target$ra) * cos(target$dec * (pi / 180))
    deltaDEC <- star_summary$dec - target$dec

    #Get the angular distance (assuming SAA) and the magnitude difference from all stars
    delta_pos <- sqrt(deltaRA^2 + deltaDEC^2)
    delta_mag <- abs(star_summary$med_mag - target$med_mag)

    #Create a mask that only allows stars that are nearby and similar (and not too variable)
    print("make mask")
    mask <- delta_pos > 0 & delta_pos < .13775 & delta_mag < 1.2 #& star_summary$rms < .05
    similar_ids <- star_summary$StarID[mask]

    #Get light curves of similar_ids
    print("make similarlcs")
    target_lc <- lcs[[as.character(target_id)]]
    similar_lcs <- lcs[similar_ids]

    #If too few neighbors, dont systematic remove
    print("small check")
    if (length(similar_lcs) < 2 || all(sapply(similar_lcs , is.null)))
    {

      target_lc$SMag <- 0
      target_lc$AMag <- target_lc$Mag
      target_lc$SysRemoved <- F

      return(list(id = target_id , lc = target_lc))

    }

    #Make matrix where columns are stars and rows are timesteps
    print("mag matrix make")
    all_times <- sort(unique(unlist(lapply(similar_lcs , function(df) df$TimeStamp))))
    mag_mat <- matrix(NA , nrow = length(all_times) , ncol = length(similar_lcs))

    print("small loop")
    for (j in seq_along(similar_lcs))
    {

      df <- similar_lcs[[j]]
      idx <- match(df$TimeStamp , all_times)
      mag_mat[idx , j] <- df$Mag

    }

    #Normalize the magnitude matrix columns by the median mag of star lcs
    print("normalize")
    norm_mag_mat <- sweep(mag_mat , 2 , colMedians(mag_mat , na.rm = T) , FUN = "-")

    #Systematic is the median normalized changes in magnitude among all light curves
    print("systematic made")
    systematic <- rowMedians(norm_mag_mat , na.rm = T)

    #Align systematic with target
    sys_df <- data.frame(
      TimeStamp = all_times ,
      SMag = systematic
    )
    target_lc <- merge(target_lc , sys_df , by = "TimeStamp" , all.x = T)
    target_lc <- target_lc[order(target_lc$TimeStamp) , ]

    #Apply systematic correction
    target_lc$SMag[is.na(target_lc$SMag)] <- 0
    target_lc$AMag <- target_lc$Mag - target_lc$SMag
    target_lc$SysRemoved <- T

    return(list(id = target_id , lc = target_lc))

  }

  process_star_safe <<- function(i , star_summary , lcs)
  {

  tryCatch(
    process_star(i , star_summary , lcs) ,
    error = function(e) {
      warning(paste0("Star " , i , " failed: " , e$message))
      return(NULL)
    }
  )

}

  removeSystematics <- function(all_data)
  {

    #Load in packages
    library(matrixStats)
    library(parallel)

    #Locally split all_data into lcs per star
    lcs <<- split(all_data , all_data$Aligned_ID)

    #Calculate star properties needed to determine nearby and similar neighbors (pixel scale of TOROS -> .4959 arcsec/pix)
    star_summary <<- do.call(rbind , lapply(lcs , function(df){
      data.frame(
        StarID = df$Aligned_ID[1] ,
        med_mag = median(df$Mag) ,
        ra = mean(df$RA) ,
        dec = mean(df$DEC) ,
        rms = sd(df$Mag)
      )
    }))

    #Run process star in parallel for speed
    num_cores <- detectCores() - 1
    print("make clust")
    cl <- makeCluster(num_cores)

    print("cluster exp")
    clusterExport(cl , varlist = c("star_summary" , "lcs" , "colMedians" , "rowMedians" , "process_star" , "process_star_safe"))
    print("clustereval")
    clusterEvalQ(cl , library(matrixStats))

    print("starreults start")
    star_results <- parLapply(
      cl , 1:nrow(star_summary) , function(i) {
        process_star_safe(i , star_summary , lcs)
      })
    print("stop clust")
    stopCluster(cl)

    #Rebuild lcs list
    for (result in star_results)
    {

      lcs[[as.character(result$id)]] <- result$lc

    }

    #Return all data
    all_data_new <- do.call(rbind , lcs)
    return(all_data_new)

  }
  all_data <- removeSystematics(all_data = all_data)


  ###Split all data into individual light curves
  starlc <- split(all_data , all_data$Aligned_ID)

  #Filter and finalize the stars
  avg_numNights <- nrow(all_data) / length(lcs)
  print(avg_numNights)
  starlc <- lapply(starlc , function(df)
                   {

                     if (nrow(df) < ceiling(avg_numNights)) {return(NULL)}

                     df <- df[order(df$TimeStamp) , ] #Orders by date and time
                     df$RunningDay <- as.numeric(difftime(df$TimeStamp , min(df$TimeStamp) , units = "days"))

                     #Reformat columns for cleanliness (yaya)
                     return(data.frame(
                       StarID = df$Aligned_ID ,
                       Date = df$Date ,
                       Time = df$Time ,
                       TimeStamp = df$TimeStamp ,
                       RunningDay = df$RunningDay ,
                       X_Position = df$x ,
                       Y_Position = df$y ,
                       RA = df$RA ,
                       DEC = df$DEC ,
                       Flux = df$Flux ,
                       FluxEr_Pois = sqrt(df$Flux) ,
                       FluxEr_Real = df$FluxError ,
                       SNR = df$SNR ,
                       PixelArea = df$PixArea ,
                       LocalBkgMedian = df$LocBkgMed ,
                       LocalBkgStdev = df$LocBkgStd ,
                       RMag = df$Mag ,
                       AMag = df$AMag ,
                       SMag = df$SMag ,
                       SysRem = df$SysRemoved
                     ))
                   })

  #Remove NULL entries (stars with <10 observations)
  starlc <- starlc[!sapply(starlc , is.null)]


  ###Return data
  return(list(LCs = starlc , NightOffsets = night_medians$offset))

}
