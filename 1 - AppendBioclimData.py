#To run this script on Terminal, give the taxon key and name as arguments like:
#   python3 AppendBioclimData.py 3007383 'Rosa alba'

import datetime, sys, gc, pandas as pd, glob
from osgeo import gdal

taxon_key = sys.argv[1]
taxon_name = sys.argv[2]

#######################################################################################################
#######################################SET SCRIPT VARIABLES############################################
#######################################################################################################

#Specific to running environment 
thinned_occ_folder_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/2-Thinned/'
bioclim_occ_folder_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/3-Bioclim/'
bioclim_folder_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Bio1980-2010/'

#######################################################################################################
############################################GET RAW DATA FROM GBIF#####################################
#######################################################################################################



def get_all_bioclim(taxon_key,taxon_name):
    print(str(datetime.datetime.now()) + ": Starting getting Bioclim features for " + taxon_name + "!")
    print(str(datetime.datetime.now()) + ": Loading thinned occurrences of " + taxon_name + "!")
    df_int=pd.read_csv(glob.glob(thinned_occ_folder_path + str(taxon_key) + '*.txt')[0], sep='@', encoding='utf-8',dtype=float)
    print(str(datetime.datetime.now()) + ": " + taxon_name + " has " + str(len(df_int)) + " thinned occurrences")
    occ_list = [(x,y) for x,y in zip(df_int['decimalLongitude'].tolist(),df_int['decimalLatitude'].tolist())]
    df = pd.DataFrame()
    for n in list(range(1,20)):
        print(str(datetime.datetime.now()) +': Start with Bioclim' + f'{n:02}')
        ##https://chelsa-climate.org/downloads/
        ##https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
        bioclim = bioclim_folder_path + 'CHELSA_bio' + str(n) + '_1981-2010_V.2.1.tif'
        dataset = gdal.Open(bioclim)
        band = dataset.GetRasterBand(1)
        cols = dataset.RasterXSize
        rows = dataset.RasterYSize
        transform = dataset.GetGeoTransform()
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = -transform[5]
        data = band.ReadAsArray(0, 0, cols, rows)
        col_values =[]
        for point in occ_list:
            row = int((yOrigin - point[1] ) / pixelHeight)
            col = int((point[0] - xOrigin) / pixelWidth)
            col_values.append(data[row][col])
        df['Bioclim' + f'{n:02}'] = col_values
        print(str(datetime.datetime.now()) +': Finished with Bioclim' + f'{n:02}')
    df['Taxon'] = taxon_name
    df = df.reset_index(drop=True)
    df.to_csv(bioclim_occ_folder_path + str(taxon_key) + '_' + taxon_name + '_'+ datetime.date.today().strftime("%d%m%Y") + '.txt', sep='@', index=False, encoding='utf-8')
    print(str(datetime.datetime.now()) + ": Saved "  + str(df.shape[0]) + " Bioclim features  for " + taxon_name + "!")
    del df_int,occ_list,df,data,col_values,bioclim,dataset,band,cols,rows,transform
    gc.collect()
    return


get_all_bioclim(taxon_key,taxon_name)
