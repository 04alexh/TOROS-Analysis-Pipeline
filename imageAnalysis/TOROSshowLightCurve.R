TOROSshowLightCurve <- function(lcs , lc_targetID , observed_field = "Observed Field" ,
                                use_nameCustom = F  , nameCustom = "")
{

  ###Load in packages
  library(ggplot2)


  ###Pull target data
  tar_star <- lcs[[as.character(lc_targetID)]]


  ###Create light curves using ggplot
  #Name the star for labels
  if (use_nameCustom == T)
  {

    star_name <- nameCustom

  }

  else
  {

    star_name <- paste0(observed_field , " " , tar_star$StarID[1])

  }

  #Make light curve using raw magnitude and adjusted magnitude
  raw <- ggplot(data = tar_star , aes(color = as.character(Date))) +
    geom_point(aes(x = RunningDay , y = RMag , shape = "Raw Magnitude")) +
    geom_point(aes(x = RunningDay , y = AMag , shape = "Adjusted Magnitude")) +
    scale_y_reverse(limits = c(min(tar_star$AMag) - 1 , max(tar_star$AMag) + 1)) +
    scale_color_discrete() +
    scale_shape_manual(
      values = c(
        "Raw Magnitude" = 4 ,
        "Adjusted Magnitude" = 16
      )
    ) +
    labs(
      title = paste0("Light Curve") ,
      subtitle = star_name ,
      x = "Running Time [days]" ,
      y = "Instrumental Magnitude" ,
      color = "Observation Date" ,
      shape = "Magnitude Type"
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
  print(raw)

  #If systematic data exists, also plot it
  if (tar_star$SysRem[1] == T)
  {

    observation_number <- seq(from = 1 , to = length(tar_star$SMag))

    sys <- ggplot(data = tar_star , aes(color = as.character(Date))) +
      geom_line(aes(x = observation_number , y = SMag)) +
      scale_y_reverse() +
      scale_color_discrete() +
      geom_hline(aes(yintercept = 0)) +
      labs(
        title = "Systematic Light Curve" ,
        subtitle = star_name ,
        x = "Observation Number" ,
        y = "Instrumental Magnitude" ,
        color = "Observation Date"
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

    print(sys)

  }

}
