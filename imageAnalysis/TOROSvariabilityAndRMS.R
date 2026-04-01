TOROSVariabilityAndRMS <- function(lcs , observed_field = "Observed Field")
{

  ###Load in packages
  library(magicaxis)
  library(ggplot2)
  library(ie2misc)


  ###Load in function to calculate Stetson J and D90
  jstetcalc <- function(df)
  {

    df <- df[order(df$TimeStamp) , ]
    mags <- df$AMag
    #mags <- sort(mags)
    counts <- 10^((-(mags - 25))/2.5)
    er <- 1 / sqrt(counts)
    rms <- sd(mags)
    nms <- length(mags)

    cut <- as.integer(.1 * nms)
    d90 <- abs(sort(mags)[cut + 1] - sort(mags)[length(mags) - cut])

    Jt <- rep(0.0 , nms)
    Jb <- rep(0.0 , nms)

    wk <- 1
    MeanMag <- mean(mags)

    for (i in seq(from = 1 ,to = nms - 2 , by = 2))
      {

        Sigi <- (mags[i] - MeanMag) / (er[i])*(sqrt(nms / (nms - 1)))
        Sigj <- (mags[i+1] - MeanMag) / (er[i+1])*(sqrt(nms / (nms - 1)))

        Pk <- Sigi*Sigj

        if (Pk > 0.0)
        {
          sgnPk <- 1.0
        }
        if (Pk == 0.0)
        {
          sgnPk <- 0.0
        }
        if (Pk < 0.0)
        {
          sgnPk <- -1.0
        }

        Jt[i] <- wk*sgnPk*(sqrt(abs(Pk)))
        Jb[i] <- (wk)

      }

    jstetson <- sum(Jt) / sum(Jb)
    result <- list(Jstetson = jstetson , D90 = d90)

    # print("Returning JStet result")
    # print(result)
    return(result)



  }

  jstetcalc_by_night <- function(df)
  {

    #Split df by dates
    night_split <- split(df , df$Date)
    night_split <- night_split[sapply(night_split , nrow) >= 3]

    #Calculate the J value for the nights
    J_vals <- sapply(night_split , function(df) {jstetcalc(df)[[1]]})
    pair_counts <- sapply(night_split , function(x) {floor(nrow(x) / 2)})

    #Get weighted average
    J_total <- sum(J_vals * pair_counts) / sum(pair_counts)
    return(J_total)

  }


  ###Iterate through light curves, determine if they are variable with J and D90
  id_list <- NULL
  j_list <- NULL
  d90_list <- NULL

  for (star in lcs)
  {

    #Obtain variability of this star
    var_results <- jstetcalc(star)

    #Save needed data
    id_list <- c(id_list , star$StarID[1])
    j_list <- c(j_list , jstetcalc_by_night(star))
    #j_list <- c(j_list , var_results$Jstetson)
    d90_list <- c(d90_list , var_results$D90)

  }


  ###Calculate thresholds for variability using J
  #Sigma clip data
  sigmaClip_j_list <- magclip(j_list , sigma = 3)$x

  #Find median and std dev for Js
  median_j <- median(sigmaClip_j_list)
  sd_j <- sd(sigmaClip_j_list)

  #Calculate threshold J as 2sigma above median
  threshJ <- median_j + 2 * sd_j


  ###Constructing variability dataframe
  var_df <- data.frame(StarID = id_list , Stetson_J = j_list , D90 = d90_list)

  #Flag stars with J higher than threshold
  var_df$J_Flag <- var_df$Stetson_J > threshJ

  #Initialize D90 flag as booleans
  var_df$D90_Flag <- F


  ###D90 flagging is more complicated since it requires binning stars by magnitude
  #Get median magnitudes of stars
  mean_mags <- sapply(lcs , function(df) {mean(df$AMag)})

  #Create bins of stars by magnitude
  segmentations <- quantile(mean_mags , probs = seq(0 , 1 , .1))
  bins <- cut(mean_mags , breaks = segmentations , include.lowest = T , labels = F)

  #Go through each 1% of magnitude and flag that way
  for (bin_idx in sort(unique(bins)))
  {

    #Get this current percentile
    lc_thisBin <- lcs[bins == bin_idx]

    #Calculate threshold D90 for this grouping of LCs
    ids_thisBin <- sapply(lc_thisBin , function(df) return(df$StarID[1]))
    var_df_thisBin <- var_df[var_df$StarID %in% ids_thisBin , ]
    medianD90_thisBin <- median(magclip(var_df_thisBin$D90 , sigma = 3)$x)
    sdD90_thisBin <- sd(magclip(var_df_thisBin$D90 , sigma = 3)$x)
    threshd90 <- medianD90_thisBin + 2*sdD90_thisBin

    #Update the D90 flag status of stars in this bin
    idx <- var_df$StarID %in% ids_thisBin
    var_df$D90_Flag[idx] <- var_df$D90[idx] > threshd90

  }


  ###Create variability histograms
  #J index
  Jplot <- ggplot(data = var_df , aes(x = Stetson_J)) +
    geom_histogram(binwidth = .5) +
    geom_vline(xintercept = threshJ , color = 'red' , linetype = 'solid' , size = .7) +
    labs(
      title = paste0("Stetson J Variability Index Values for " , observed_field) ,
      subtitle = paste0("Threshold: " , threshJ) ,
      x = "J Value" ,
      y = "Count"
    ) +
    theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
    )

  d90plot <- ggplot(data = var_df , aes(x = D90)) +
    geom_histogram(binwidth = .005) +
    labs(title = paste0("D90 Values for " , observed_field) ,
         x = "D90 Value" ,
         y = "Count") +
    theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
    )


  ###Now it is also possible to see how RMS of stars changes with magnitude
  #Initialize lists
  id_list <- NULL
  rms_list <- NULL
  medMag_list <- NULL
  status_list <- NULL

  for (star in lcs)
  {

    #Find median and RMS
    med_mag <- median(star$AMag)
    rms_mag <- 1.4826 * mad(star$AMag)

    #Append to lists
    id_list <- c(id_list , star$StarID[1])
    medMag_list <- c(medMag_list , med_mag)
    rms_list <- c(rms_list , rms_mag)

    #Create variability status based on the flag from var_df
    status <- "No Flag"

    if (var_df[var_df$StarID == star$StarID[1] , ]$J_Flag == T && var_df[var_df$StarID == star$StarID[1] , ]$D90_Flag == T)
    {
      status <- "J + D90 Flag"
    }

    else if (var_df[var_df$StarID == star$StarID[1] , ]$J_Flag == T)
    {
      status <- "J Flag"
    }

    else if (var_df[var_df$StarID == star$StarID[1] , ]$D90_Flag == T)
    {
      status <- "D90 Flag"
    }

    status_list <- c(status_list , status)

  }


  ###Create a data frame to hold RMS info
  rms_df <- data.frame(
    StarID = id_list ,
    Median_Mag = medMag_list ,
    RMS = rms_list ,
    Variability_Status = status_list ,
    stringsAsFactors = F
  )

  #Plot RMS vs. Mag
  RMSplot <- ggplot(data = rms_df , aes(x = Median_Mag , y = RMS , color = Variability_Status)) +
    geom_point() +
    scale_color_manual(values = c('No Flag' = 'black' , 'J + D90 Flag' = 'purple' , 'J Flag' = 'red' , 'D90 Flag' = 'blue' )) +
    scale_y_log10() +
    labs(
      title = paste0("RMS vs. Instrumental Magnitude") ,
      subtitle = observed_field ,
      x = "Instrumental Magnitude" ,
      y = "RMS" ,
      color = "Variability Flag"
    ) +
    theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
    )


  ###Print all plots and return data
  print(Jplot)
  print(d90plot)
  print(RMSplot)
  return(list(Variability = var_df , RMSData = rms_df))

}
