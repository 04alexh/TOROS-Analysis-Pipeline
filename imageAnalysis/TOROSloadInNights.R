TOROSloadInNights <- function(night_folders)
{

  ###Iterate through folders, then for each folder read the csvs inside
  image_photometry_list <- list()
  caught_ids <- list()
  idx <- 1

  for (folder in night_folders)
  {

    #Load files in folder into a list
    night_file_list <- list.files(path = folder , full.names = T)

    for (file in night_file_list)
    {

      #Read in csv as a dataframe
      df <- tryCatch(
               read.csv(file = file) ,
               error = function(e) {
                 warning(paste0("Failed to read file: " , file))
                 return(NULL)
               })

      if (is.null(df)) {next}

      #Append lists
      image_photometry_list[[idx]] <- df
      caught_ids[[idx]] <- df$Aligned_ID
      idx <- idx + 1

    }

  }

  #Create lists of star ids in all fields and all unique ids
  common_ids <- Reduce(f = intersect , caught_ids)
  all_ids <- unique(unlist(caught_ids))


  ###Return list containing the ID lists and the photometry
  return(list(CommonIDs = common_ids , AllIDs = all_ids , Photometry = image_photometry_list))

}
