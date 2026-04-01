TOROSobsTimeLSSTpredict <- function(lcs , test_targetID , sigma_level = 2 , N = 1000 , M = 10000 ,
                                    use_nameCustom = F , nameCustom = "" , observed_field = "Observed Field")
{

  ###Load in packages
  library(magicaxis)
  library(ggplot2)
  library(dplyr)


  ###Load in required functions
  ###Load in function to calculate Stetson J and D90
  jstetcalc <- function(df)
  {

    mags <- df$AMag
    mags <- sort(mags)
    counts <- 10^((-(mags - 25))/2.5)
    er <- 1 / sqrt(counts)
    rms <- sd(mags)
    nms <- length(mags)

    cut <- as.integer(.1 * nms)
    d90 <- abs(mags[cut + 1] - mags[length(mags) - cut])

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


  #####THIS SIM WILL FIND THE MINIMUM NUMBER OF OBSERVATIONS REQUIRED TO RECOVER VARIABILITY AND THE TIME INTERVALS OF OBSERVATION


  ###Load in target data
  target_star <- lcs[[as.character(test_targetID)]]
  lcs[[as.character(test_targetID)]] <- NULL #Remove target to avoid contamination

  #Initialize vectors
  thresh_list <- NULL
  J_Stet_Comp_Table <- data.frame()


  ###Iterate through sample sizes 3 to the number of observations of the target
  for (i in 3:nrow(target_star))
  {

    print(paste0("Sample size = " , i , " observations"))

    ###Loop through light curves of other stars to compute the threshold J for this sample size
    print("Computing threshold J")
    compJ_list <- NULL

    for (lc in lcs)
    {

      print(paste0("Analyzing star " , lc$StarID[1]))

      #Sample lc with sample size i
      if (nrow(lc) < i) {next}
      sample_lc_idx <- sample(seq_len(nrow(lc)) , size = i)
      cut_lc <- lc[sample_lc_idx , ]

      #Calculate J value of this sample
      J <- jstetcalc(cut_lc)[[1]]
      compJ_list <- c(compJ_list , J)

    }

    ###Compute threshold J for this sample of size i
    cur_medJ <- median(magclip(compJ_list , sigma = 2)$x)
    cur_sdJ <- sd(magclip(compJ_list , sigma = 2)$x)

    cur_threshJ <- cur_medJ + (cur_sdJ * sigma_level)
    thresh_list <- c(thresh_list , cur_threshJ)

    ###Now sample from the target star N times seeing if it passes through this threshold
    print("Computing target Js")
    tarJ_List <- numeric(N)

    for (j in 1:N)
    {

      #Same as previous process
      sample_tar_idx <- sample(seq_len(nrow(target_star)) , size = i)
      cut_tarlc <- target_star[sample_tar_idx , ]

      tarJ <- jstetcalc(cut_tarlc)[[1]]
      tarJ_List[j] <- tarJ

    }

    J_Stet_Comp_Table <- rbind(J_Stet_Comp_Table , data.frame(t(as.matrix(tarJ_List))))

  }


  ###Bind threshold values to last column of table
  J_Stet_Comp_Table <- cbind(J_Stet_Comp_Table , thresh_list)

  #Name columns and rows
  colnames(J_Stet_Comp_Table) <- paste0("Sample Number: " , 1:N)
  colnames(J_Stet_Comp_Table)[ncol(J_Stet_Comp_Table)] <- "Threshold J"
  rownames(J_Stet_Comp_Table) <- paste0("Sample Size: " , 3:nrow(target_star))


  ###Create a plot of sample-pass% vs sample size
  #Create a vector containg values of %pass data for each sample size
  counts_passed <- NULL
  sample_df <- J_Stet_Comp_Table

  for (ii in 1:nrow(sample_df))
  {

    sample_df[ii , ] <- as.numeric(sample_df[ii , 1:ncol(sample_df)]) > sample_df[ii , ncol(sample_df)]
    counts_passed <- c(counts_passed , sum(sample_df[ii , ]))

  }

  #Create a df to make plotting easier
  sample_sizes <- seq(from = 3 , to = nrow(sample_df) + 2)
  new_df <- data.frame(Sample_Sizes = sample_sizes , Passed = counts_passed , Percent_Passed = (counts_passed) / N)

   #Name the star for labels
  if (use_nameCustom == T)
  {

    star_name <- nameCustom

  }

  else
  {

    star_name <- paste0(observed_field , " " , tar_star$StarID[1])

  }

  #Find where samples first hit 100%
  perf_idx <- which(new_df$Percent_Passed == 1)[1]
  xmax <- new_df$Sample_Sizes[perf_idx + 20]

  #Create plot
  obsplot <- ggplot(data = new_df , aes(x = Sample_Sizes)) +
    geom_line(aes(y = Percent_Passed) , size = .8) +
    geom_hline(aes(yintercept = 1) , color = 'red' , linetype = 'dotted' , size = .9) +
    xlim(0 , xmax) +
    labs(
      title = "% of Passing Samples vs. Sample Sizes" ,
      subtitle = star_name ,
      x = "Sample Size (Number of Observations)" ,
      y = "Percentage Passed"
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

  ###Now iterate through the target using the optimal observation number and record J vs. time interval
  I <- new_df$Sample_Sizes[perf_idx]
  J_list <- NULL
  tint_list <- NULL

  print("Beginning interval stage...")
  for (j in 1:M)
  {

    print(paste0("Iteration " , j))

    #Take sample
    sample_tar_idx <- sample(seq_len(nrow(target_star)) , size = I)
    cut_tarlc <- target_star[sample_tar_idx , ]
    cut_tarlc <- cut_tarlc %>%
      arrange(TimeStamp)

    #Calculate J
    J <- jstetcalc(cut_tarlc)[[1]]
    J_list <- c(J_list , J)

    #Find time interval between earliest and latest selected lc point
    time_int <- cut_tarlc$TimeStamp[nrow(cut_tarlc)] - cut_tarlc$TimeStamp[1]
    time_int <- as.numeric(time_int , units = "days")
    tint_list <- c(tint_list , time_int)

  }

  ###Build another dataframe
  J_timeint_df <- data.frame(Time_Interval = tint_list , J = J_list)

  #Truncate time since it only matters how many nights of observations, not how many hours,mins,sec,etc
  J_timeint_df$Time_Interval <- trunc(J_timeint_df$Time_Interval , 2)


  ###Now obtain variability measures for each time interval
  uniq_tint <- unique(J_timeint_df$Time_Interval)
  Jm_list <- NULL
  Jsd_list <- NULL
  Jtime_list <- NULL

  for (time in uniq_tint)
  {

    #Obtain rows with this time interval
    cut_J_df <- J_timeint_df[J_timeint_df$Time_Interval == time , ]
    if (nrow(cut_J_df) < 2) {next} #Ignore cuts with not enough data

    J_median <- median(magclip(cut_J_df$J , sigma = 2)$x)
    J_sd <- sd(magclip(cut_J_df$J , sigma = 2)$x)
    J_time <- cut_J_df$Time_Interval[1]

    Jm_list <- c(Jm_list , J_median)
    Jsd_list <- c(Jsd_list , J_sd)
    Jtime_list <- c(Jtime_list , J_time)

  }

  #Construct new data frame for easy plotting
  new_Jtime_df <- data.frame(Time_Interval = Jtime_list , J_median = Jm_list , J_sd = Jsd_list)

  ###Plot
  intPlot <- ggplot(data = new_Jtime_df , aes(x = Time_Interval , y = J_median)) +
    geom_point(aes(color = "Median") , size = 1) +
    geom_errorbar(aes(ymin = J_median - J_sd , ymax = J_median + J_sd , color = "1 StdDev") , width = .01) +
    labs(
      title = paste0("J vs. Time Interval of Observations") ,
      subtitle = paste0(star_name , " | Sample Size: " , I , ", Iterations: " , M),
      x = "Time Interval [days]" ,
      y = "Median J" ,
    ) +
    scale_color_manual(name = "Legend" , values = c(
      "Median" = "black" , "1 StdDev" = "red"
    )) +
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

  ###Print everything and return results
  print(obsplot)
  print(intPlot)

  LSSTPred1 <- list(J_Comp = J_Stet_Comp_Table , Sample_Data = new_df)
  LSSTPred2 <- list(Time_Int_Data = new_Jtime_df)

  return(list(LSSTPred1 , LSSTPred2))
}
