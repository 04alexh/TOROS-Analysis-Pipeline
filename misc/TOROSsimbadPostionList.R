TOROSsimbadPositionList <- function(lcs , var_df)
{

  ### This program will make an ascii file suitable to be loaded into a Simbad Position Query using the variability dataframe
  ### made from TOROSvariabilityAndRMS

  ###Sort var_df by J
  a1 <- var_df[order(var_df$Stetson_J , decreasing = T) , ]


  ###Cut stars with false J flags
  a2 <- a1[a1$J_Flag == T , ]


  ###Pull only the StarIDs
  a3 <- a2$StarID


  ###Iterate through a3 and pull RA and DEC from each star in LCs
  RAs <- NULL
  DECs <- NULL
  for (StarID in a3)
  {

    current <- lcs[[as.character(StarID)]]
    RA <- current$RA[1]
    DEC <- current$DEC[1]

    #Append
    RAs <- c(RAs , RA)
    DECs <- c(DECs , DEC)

  }


  ###Make dataframe
  simbad_df <- data.frame(RA = RAs , DEC = DECs , ID = as.integer(a3))


  ###Write simbad_df as a ascii
  write.csv(simbad_df , "simbadpos.csv" , row.names = F , col.names = T)

}
