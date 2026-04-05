TOROSshowPhaseFoldCurve <- function(lcs , lc_targetID , period , observed_field = "Observed Field" ,
                                    use_nameCustom = F , nameCustom = "" , bin = F)
{

  ###Load in packages
  library(ggplot2)


  ###Load in target star data
  tar_star <- lcs[[as.character(lc_targetID)]]

  #Create a phase column
  t0 <- tar_star$RunningDay[which.min(tar_star$AMag)]
  tar_star$Phase <- ((tar_star$RunningDay - t0) / period) %% 1


  ###Create ggplot object
  #Name the star for labels
  if (use_nameCustom == T)
  {

    star_name <- nameCustom

  }

  else
  {

    star_name <- paste0(observed_field , " " , tar_star$StarID[1])

  }


  #If bin is true, the magnitudes for each unique date will be averaged and plotted at the first phase value of that date
  if (bin == T)
  {

    #Iterate through unique dates
    meds <- NULL
    phases <- NULL
    for (date in unique(tar_star$Date))
    {

      #Get data for this date
      cur_df <- tar_star[tar_star$Date == date , ]

      #Get median mag and phase
      med <- median(cur_df$AMag)
      phase <- median(cur_df$Phase)

      #Append
      meds <- c(meds , med)
      phases <- c(phases , phase)

    }

    #Create dataframe
    binned_data <- data.frame(Phase = phases , AMag = meds , Date = unique(tar_star$Date))

    #Double the phase plot to see pattern easier
    binned_data2 <- binned_data
    binned_data2$Phase <- binned_data2$Phase + 1
    binned_data <- rbind(binned_data , binned_data2)

    #Make plot
    binned_phaseFold <- ggplot(data = binned_data , aes(x = Phase , y = AMag , color = as.character(Date))) +
      geom_point() +
      scale_y_reverse() +
      labs(
        title = "Phase Folded Light Curve" ,
        subtitle = star_name ,
        x = "Phase" ,
        y = "Instrumental Magnitude" ,
        color = "Observation Date"
      ) +
      scale_color_discrete() +
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

    print(binned_phaseFold)
    return()

  }

  #Repeat plot once so cycle easier to see
  tar_star2 <- tar_star
  tar_star2$Phase <- tar_star$Phase + 1
  tar_star <- rbind(tar_star , tar_star2)

  #Make plot
  phaseFold <- ggplot(data = tar_star , aes(x = Phase , y = AMag , color = as.character(Date))) +
    geom_point() +
    scale_y_reverse() +
    labs(
      title = "Phase Folded Light Curve" ,
      subtitle = star_name ,
      x = "Phase" ,
      y = "Instrumental Magnitude" ,
      color = "Observation Date"
    ) +
    scale_color_discrete() +
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

  print(phaseFold)

}
