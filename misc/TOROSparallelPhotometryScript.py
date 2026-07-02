import os
import gc
from astropy.table import Table
from multiprocessing import Pool , cpu_count

### REPLACE r"" STRINGS WITH APPROPRIATE PATHS!!!

science_folder = r""
starList_folder = r""
master = Table.read(
  r"" , format = "csv" , delimiter = ","
)

def process_file(file):

  # Will get photometry and align files in parallel yaya

  science_directory = os.path.join(science_folder , file)
  aligned_directory = os.path.join(starList_folder , file.replace(".fits" , ".csv"))

  tab = TOROSaperturePhotometry(
    science_file = science_directory ,
    temp_dir = r"" , #temporary directory for intermediate files, you can make this whatevs I just have it so my folders are organized :)
    write_table = False  ,
    use_mask = True ,
    cx = 0 , cy , 0 , r = 0 # define your masks accordingly if u need them
  )

  TOROSphotometryAlign(
    master = master ,
    comparator = tab ,
    aligned_file = aligned_directory ,
    snr_threshold = 10
  )

# Run process_file in parallel for speed :)
if __name__ = "__main__":
  with Pool(4) as p: #I have a 16core laptop with 32 gigs of ram and I literally couldnt run with more than 4 threads change with caution!!!!
    p.map(process_file , file_list)
  
