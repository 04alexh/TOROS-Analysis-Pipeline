TOROSlombScarg <- function(lcs , lsp_targetID , observed_field = "Observed Field" ,
                           use_nameCustom = F , nameCustom = "" ,
                           phaseFold = F , phaseFold_num = 0)
{

  ###Load in packages
  library(ggplot2)
  library(lomb)


  ###Create lomb scargle object
  #Load into target data
  tar_star <- lcs[[as.character(lsp_targetID)]]

  #Run LSP on tar_star
  lomb <- lsp(x = tar_star$AMag , times = tar_star$RunningDay , type = "period" , ofac = 3 , plot = F)

  #Find top ranked periods
  ord <- order(lomb$power , decreasing = T)
  top_periods <- data.frame(Period = lomb$scanned[ord] , Power = lomb$power[ord])

  #Clean up harmonics
  clean_periods <- top_periods

  for (i in 1:(nrow(clean_periods) - 1))
  {

    ratios <- clean_periods$Period[i] / clean_periods$Period[(i + 1):nrow(clean_periods)] #Get ratios with all other periods
    bad_ratios <- which(abs(ratios - round(ratios)) < .01) #Get indices of periods with near integer ratios

    #Remove bad ratios
    if (length(bad_ratios) != 0)
    {

      clean_periods <- clean_periods[-(i + bad_ratios)]

    }

  }

  ###Plot LSP
  #Name the star for labels
  if (use_nameCustom == T)
  {

    star_name <- nameCustom

  }

  else
  {

    star_name <- paste0(observed_field , " " , tar_star$StarID[1])

  }


  ###Use LSP to phase fold if desired
  if (phaseFold == T)
  {

    #Load phase fold program
    TOROSshowPhaseFoldCurve <- function(lcs , lc_targetID , period , observed_field = "Observed Field" ,
                                    use_nameCustom = F , nameCustom = "")
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

    #Use top N periods to phase fold
    for (i in 1:phaseFold_num)
    {

      period <- top_periods$Period[i]
      TOROSshowPhaseFoldCurve(
        lcs = lcs , lc_targetID = lsp_targetID , period = period ,
        use_nameCustom = T , observed_field = observed_field , nameCustom = paste0(star_name , " LSP Period " , i , ":" , period)
      )

    }




  }

  lsp_plot <- ggplot() +
    geom_line(aes(x = clean_periods$Period , y = clean_periods$Power , color = "LSP") , size = 1) +
    geom_hline(aes(yintercept = lomb$sig.level , color = "P < .01") , linetype = "dotted" , size = 1.5) +
    geom_vline(aes(xintercept = clean_periods$Period[1] , color = "1st Period") , linetype = "dashed" , size = .5) +
    geom_vline(aes(xintercept = clean_periods$Period[2] , color = "2nd Period") , linetype = "dashed" , size = .5) +
    geom_vline(aes(xintercept = clean_periods$Period[3] , color = "3rd Period") , linetype = "dashed" , size = .5) +
    labs(
      title = "Lomb-Scargle Periodogram" ,
      subtitle = star_name ,
      x = "Period [days]" ,
      y = "Power"
    ) +
    scale_color_manual(
      name = "Legend" ,
      values = c(
        "LSP" = "black" ,
        "P < .01" = "red" ,
        "1st Period" = "purple" ,
        "2nd Period" = "blue" ,
        "3rd Period" = "cyan"
      )
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

  print(lsp_plot)
  return(lomb)

}
