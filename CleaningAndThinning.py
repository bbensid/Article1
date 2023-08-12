#To run this script on Terminal, give the taxon name as argument like:
#   python3 CleaningAndThinning.py 'Populus alba'


import datetime, sys, os, time, gc, pandas as pd, geopandas, random, numpy as np
from pygbif import species
from pygbif import occurrences as occ
import shapely.wkt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from zipfile import ZipFile
import shapefile as shp

taxon_name = sys.argv[1]

#######################################################################################################
#######################################SET SCRIPT VARIABLES############################################
#######################################################################################################

#Specific to running environment and user
gbif_user = 'badisbens'
gbif_password ='cocobanane'
raw_occ_folder_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/0-Raw/'
cleaned_occ_folder_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/1-Cleaned/'
thinned_occ_folder_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/2-Thinned/'
institutions_path = '/Users/home/Documents/PhD/GBIF_source/CoordinateCleanerPy/institutions_df.txt'
country_centroids_path = '/Users/home/Documents/PhD/GBIF_source/CoordinateCleanerPy/country_centroid_df.txt'
capitals_path = '/Users/home/Documents/PhD/GBIF_source/CoordinateCleanerPy/capital_df.txt'
land_path = '/Users/home/Documents/PhD/GBIF_source/biomes/official/wwf_terr_ecos.shp'
bioclim_1_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Bio1980-2010/CHELSA_bio1_1981-2010_V.2.1.tif'

#Always true
crs = '+proj=longlat +datum=WGS84' 

#######################################################################################################
############################################GET RAW DATA FROM GBIF#####################################
#######################################################################################################

##Get the Taxon GBIF key

def get_key(taxon_name):
    i=0
    while i<3:
        try:
            taxon_key = species.name_backbone(name=taxon_name, kingdom='plants', genus = taxon_name.split(' ')[0])['speciesKey']
            return taxon_key
            break
        except:
            print(str(datetime.datetime.now()) + ": Error connecting to GBIF, waiting...")
            time.sleep(20)
            i+=1
    print(str(datetime.datetime.now()) + ": " + taxon_name + " not found...")
    sys.exit()


taxon_key = get_key(taxon_name)

##Get the raw occurrences fom GBIF

def get_raw(taxon_key,taxon_name):
    print(str(datetime.datetime.now()) + ": Fetching occurrences for " + taxon_name + "!")
    occ.count(taxonKey = str(taxon_key),isGeoreferenced = True)
    a = occ.download(queries=['taxonKey = ' + str(taxon_key), 'hasCoordinate = TRUE','hasGeospatialIssue = False','decimalLatitude < 90','decimalLongitude < 180','decimalLatitude > -90','decimalLongitude > -180'],user=gbif_user, pwd=gbif_password, email='')
    while True:
        try:
            b=occ.download_meta(key=a[0])
        except:
            print(str(datetime.datetime.now()) + ": Error connecting to GBIF, waiting...")
            time.sleep(20)
        if b['status'] == 'SUCCEEDED':
            print(str(datetime.datetime.now()) + ": Your download is ready!")
            x=occ.download_get(key=a[0],path=raw_occ_folder_path)
            raw_occ_path = raw_occ_folder_path + str(taxon_key) + '_' + taxon_name + '_doi_' + b['doi'].replace('.','point').replace('/','fslash') + '_'+ datetime.datetime.strptime(b['created'], '%Y-%m-%dT%H:%M:%S.%f%z').strftime("%d%m%Y") + '.zip'
            os.rename(raw_occ_folder_path + b['downloadLink'].split('/')[-1], raw_occ_path)
            break
        else:
            print(str(datetime.datetime.now()) + ": Status: '" + str(b['status']) + "'")
            time.sleep(20)
    print(str(datetime.datetime.now()) + ": Occurrences fetched for " + str(taxon_key))
    return raw_occ_path


raw_occ_path = get_raw(taxon_key,taxon_name)

#######################################################################################################
############################################CLEANING RAW DATA##########################################
#######################################################################################################

def read_shapefile(sf):
    """
    Read a shapefile into a Pandas dataframe with a 'coords' 
    column holding the geometry information. This uses the pyshp
    package
    reference: https://towardsdatascience.com/mapping-geograph-data-in-python-610a963d2d7f
    """
    fields = [x[0] for x in sf.fields][1:]
    records = sf.records()
    shps = [s.points for s in sf.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    return df


def get_outmost_polygons(coord_list):
    i=0
    index_list=[i]
    while i < len(coord_list):
        index_list.append(coord_list[i+1:].index(coord_list[i]) + i + 2)
        i=index_list[-1] 
    coord_list_new=[]
    while len(index_list) > 1:
        coord_list_new.append(coord_list[index_list[0]:index_list[1]])
        del index_list[0]
    return coord_list_new[0]


def clean_occ(taxon_key,taxon_name):
    print(str(datetime.datetime.now()) + ": Start cleaning raw occurences for " + taxon_name)
    print(str(datetime.datetime.now()) + ": Starting preparing reference dataframes for cleaning...")
    #Get instituions, country centroids and world capitals data
    institution_df = pd.read_csv(institutions_path, sep='@', encoding='utf-8')[['polygon','polygon_min_lat','polygon_max_lat','polygon_min_lon','polygon_max_lon']]
    country_centroid_df = pd.read_csv(country_centroids_path, sep='@', encoding='utf-8')[['polygon','polygon_min_lat','polygon_max_lat','polygon_min_lon','polygon_max_lon']]
    capital_df = pd.read_csv(capitals_path, sep='@', encoding='utf-8')[['polygon','polygon_min_lat','polygon_max_lat','polygon_min_lon','polygon_max_lon']]
    institution_df = institution_df.drop_duplicates()
    #Drop duplicates
    country_centroid_df = country_centroid_df.drop_duplicates()
    capital_df = capital_df.drop_duplicates()
    #Turn the polygon column into actual Polygon object
    institution_df['polygon'] = institution_df.apply(lambda x: shapely.wkt.loads(x['polygon']), axis=1)
    country_centroid_df['polygon'] = country_centroid_df.apply(lambda x: shapely.wkt.loads(x['polygon']), axis=1)
    capital_df['polygon'] = capital_df.apply(lambda x: shapely.wkt.loads(x['polygon']), axis=1)
    #Transform dataframes into GeoDataframes
    institution_gdf = geopandas.GeoDataFrame(institution_df, geometry=institution_df['polygon'],crs=crs)
    country_centroid_gdf = geopandas.GeoDataFrame(country_centroid_df, geometry=country_centroid_df['polygon'],crs=crs)
    capital_gdf = geopandas.GeoDataFrame(capital_df, geometry=capital_df['polygon'],crs=crs)
    #Get the WWF land data
    sf = shp.Reader(land_path, encoding = 'ISO8859-1')
    df_biomes = read_shapefile(sf)
    df_biomes['new_coords'] = df_biomes.apply(lambda p: get_outmost_polygons(p.coords), axis=1)
    df_biomes['polygon'] = df_biomes.apply(lambda p: Polygon(p.new_coords), axis=1)
    ###the fact is the biome is at 30' precision, this means that island less than 1km2 won't have a polygon, or 1km from any shoreline could be out/in a polygon, this is the limitation, but let's live with this...
    df_biomes = df_biomes.drop(df_biomes[df_biomes.ECO_NAME == 'Lake'].index)
    gdf_biomes = geopandas.GeoDataFrame(df_biomes, geometry=df_biomes['polygon'],crs=crs)
    print(str(datetime.datetime.now()) + ": Finished preparing reference dataframes for cleaning")
    del institution_df,country_centroid_df,capital_df,sf,df_biomes
    gc.collect()
    zip_file = ZipFile(raw_occ_path)
    print(str(datetime.datetime.now()) + ": Loading raw data...")
    df = pd.read_csv(zip_file.open('occurrence.txt'),delimiter='\t',on_bad_lines='skip', low_memory=False)[['decimalLatitude','decimalLongitude', 'datasetKey','basisOfRecord']]
    print(str(datetime.datetime.now()) + ": " + taxon_name + " has " + str(len(df)) + " raw occurrences")
    print(str(datetime.datetime.now()) + ": Starting cleaning...")
    ##remove any NaN
    if not df.empty:
        df = df.drop(df[np.isnan(df.decimalLatitude)|np.isnan(df.decimalLongitude)].index)
    ##remove fossils
    if not df.empty:
        df = df.drop(df[df.basisOfRecord == 'FOSSIL_SPECIMEN'].index)
    ##remove datasets from iNaturalist and Pl@ntnet
    if not df.empty:
        df = df.drop(df[df.datasetKey.isin(['50c9509d-22c7-4a22-a47d-8c48425ef4a7','7a3679ef-5582-4aaa-81f0-8c2545cafc81','14d5676a-2c54-4f94-9023-1e8dcd822aa0'])].index)
    ##remove lat and long equal, or equal to zero
    if not df.empty:
        df = df.drop(df[(df.decimalLongitude == df.decimalLatitude) | (df.decimalLongitude == 0) | (df.decimalLatitude == 0)].index)
    ##remove duplicates
    df = df[['decimalLatitude','decimalLongitude']]
    df = df.drop_duplicates(ignore_index=True)
    if df.empty:
        print(str(datetime.datetime.now()) + ": Finished cleaning and saving...")
        cleaned_occ_path = cleaned_occ_folder_path + str(taxon_key) + '_' + taxon_name + '_'+ datetime.date.today().strftime("%d%m%Y") + '.txt'
        df.to_csv(cleaned_occ_path, sep='@', index=False, encoding='utf-8')
        print(str(datetime.datetime.now()) + ": Saved "  + str(len(df)) + " cleaned occurences for " + taxon_name)
        return cleaned_occ_path
    gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df.decimalLongitude, df.decimalLatitude),crs=crs)
    ##get occurrences too close to institutions
    merged_institution = geopandas.sjoin(gdf, institution_gdf, how='left', predicate='within')
    gdf = merged_institution.loc[np.isnan(merged_institution.index_right)][['decimalLongitude','decimalLatitude','geometry']]
    ##get occurrences too close to country centroids
    merged_country_centroid = geopandas.sjoin(gdf, country_centroid_gdf, how='left', predicate='within')
    gdf = merged_country_centroid.loc[np.isnan(merged_country_centroid.index_right)][['decimalLongitude','decimalLatitude','geometry']]
    ##get occurrences too close to world capitals
    merged_capital = geopandas.sjoin(gdf, capital_gdf, how='left', predicate='within')
    gdf = merged_capital.loc[np.isnan(merged_capital.index_right)][['decimalLongitude','decimalLatitude','geometry']]
    ##get occurrences in sea
    merged_sea = geopandas.sjoin(gdf, gdf_biomes, how='left', predicate='within')
    df = merged_sea[~np.isnan(merged_sea['index_right'])][['decimalLatitude', 'decimalLongitude']].drop_duplicates(ignore_index=True)
    ##save
    print(str(datetime.datetime.now()) + ": Finished cleaning and saving...")
    cleaned_occ_path = cleaned_occ_folder_path + str(taxon_key) + '_' + taxon_name + '_'+ datetime.date.today().strftime("%d%m%Y") + '.txt'
    df.to_csv(cleaned_occ_path, sep='@', index=False, encoding='utf-8')
    print(str(datetime.datetime.now()) + ": Saved "  + str(df.shape[0]) + " cleaned occurences for " + taxon_name)
    return cleaned_occ_path



cleaned_occ_path = clean_occ(taxon_key,taxon_name)

#######################################################################################################
############################################THINNING CLEANED DATA######################################
#######################################################################################################

###PREPARE ALL FUNCTIONS    

def get_first_grid(geopoint):
    xmin, ymin, xmax, ymax= geopoint.total_bounds
    xdiff=xmax-xmin
    ydiff=ymax-ymin
    original_cell_size=min(xdiff,ydiff)/2
    n=np.log(original_cell_size/min_cell_size)/np.log(2)
    decimal_part = n%1
    n_round_for_next=n//1
    cell_size =original_cell_size/(2**decimal_part)
    grid_cells = []
    for x0 in np.linspace(xmin, xmin + (int(xdiff/cell_size))*cell_size,int(xdiff/cell_size)+1):
        for y0 in np.linspace(ymin, ymin + (int(ydiff/cell_size))*cell_size,int(ydiff/cell_size)+1):
            x1 = x0+cell_size
            y1 = y0+cell_size
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
    cell = geopandas.GeoDataFrame(grid_cells, columns=['geometry'],crs=crs)
    return cell,cell_size,n_round_for_next

def get_selected_cells(cell,geopoint):
    merged = geopandas.sjoin(geopoint, cell, how='left', predicate='within')
    df_on_grid = merged[np.isnan(merged['index_right'])][['decimalLongitude','decimalLatitude','geometry']]
    merged = merged.drop(merged[np.isnan(merged['index_right'])].index)
    if not df_on_grid.empty:
        merged2 = geopandas.sjoin(df_on_grid, cell, how='left', predicate='touches')
        merged2=merged2.reset_index(drop=True)
        merged2 = merged2.loc[merged2.groupby(['decimalLongitude','decimalLatitude']).index_right.idxmin()].reset_index(drop=True)
        if not merged2.empty:
            merged=pd.concat([merged,merged2])
    return cell.loc[cell.index.isin(merged['index_right'].tolist())]

def get_next_grid(cell):
    cell2=cell.copy()
    cell2['x_list'] = cell2.apply(lambda p: list(p['geometry'].exterior.coords.xy[0]), axis=1)
    cell2['y_list'] = cell2.apply(lambda p: list(p['geometry'].exterior.coords.xy[1]), axis=1)
    cell2['center_x'] = cell2.apply(lambda x: min(x['x_list']) + (max(x['x_list'])-min(x['x_list']))/2, axis=1)
    cell2['center_y'] = cell2.apply(lambda x: min(x['y_list']) + (max(x['y_list'])-min(x['y_list']))/2, axis=1)
    cell2['xy_list'] = cell2.apply(lambda p: [(a,b) for a,b in zip(p['x_list'],p['y_list'])], axis=1)
    new_cells = cell2.explode(column='xy_list',ignore_index=True)[['xy_list','center_x','center_y']]
    new_cells=new_cells.drop_duplicates()
    new_cells['side']=new_cells.apply(lambda p: min(max(p['xy_list'][0],p['center_x'])-min(p['xy_list'][0],p['center_x']),max(p['xy_list'][1],p['center_y'])-min(p['xy_list'][1],p['center_y'])), axis=1)
    new_cells['geometry']=new_cells.apply(lambda p: shapely.geometry.box(min(p['xy_list'][0],p['center_x']),min(p['xy_list'][1],p['center_y']),max(p['xy_list'][0],p['center_x']),max(p['xy_list'][1],p['center_y'])), axis=1)
    cell_size=min(new_cells['side'].tolist())
    new_cells = geopandas.GeoDataFrame(new_cells.drop(['xy_list', 'center_x', 'center_y', 'side'],axis=1).reset_index(drop=True),crs=crs)
    return new_cells, cell_size

def thin_occ(taxon_key,taxon_name):
    print(str(datetime.datetime.now()) + ": Start thinning cleaned occurences for " + taxon_name)
    ###GET MIN CELL SIZE BASED ON ANY OF THE BIOCLIM RASTERS
    dataset = gdal.Open(bioclim_1_path)
    transform = dataset.GetGeoTransform()
    min_cell_size = transform[1]
    del dataset
    gc.collect()
    print(str(datetime.datetime.now()) + ": Minimum cell size for thinning based on the bioclim rasters precision is: " + str(min_cell_size))
    df = pd.read_csv(cleaned_occ_path, sep='@', encoding='utf-8')
    print(str(datetime.datetime.now()) + ": Loaded "  + str(df.shape[0]) + " cleaned occurences for " + taxon_name)
    gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df.decimalLongitude, df.decimalLatitude),crs=crs)
    print(str(datetime.datetime.now()) + ": " + taxon_name + " has "+ str(gdf.shape[0]) +" cleaned occurrences")
    del df
    gc.collect()
    if len(gdf) > 19:
        cell,cell_size,n_round_for_next = get_first_grid(gdf)
        print(str(datetime.datetime.now()) + ": original cell df has " + str(cell.shape[0]) + " rows")
        ##Recursively reduce grid cell size and select only cells that encapsulate data points
        while n_round_for_next>0:
            next_cell,cell_size=get_next_grid(cell)
            print(str(datetime.datetime.now()) + ": next_cell df has " + str(next_cell.shape[0]) + " rows")
            cell = get_selected_cells(next_cell,gdf)
            print(str(datetime.datetime.now()) + ": cell_size="+str(cell_size) + ": updated cell df has " + str(cell.shape[0]) + " rows")
            n_round_for_next -= 1
        ##Now directly transform these cells into a list of their centroids as thinned occurrences:
        cell['decimalLongitude'] = cell.apply(lambda p: p['geometry'].centroid.coords[0][0], axis=1)
        cell['decimalLatitude'] = cell.apply(lambda p: p['geometry'].centroid.coords[0][1], axis=1)
        final_df = cell[['decimalLatitude','decimalLongitude']]
        print(str(datetime.datetime.now()) + ":final_df df has "+ str(final_df.shape[0]) +" rows")
        thinned_occ_path = thinned_occ_folder_path + str(taxon_key) + '_' + taxon_name + '_'+ datetime.date.today().strftime("%d%m%Y") + '.txt'
        if len(final_df) > 19:
            final_df.to_csv(thinned_occ_path, sep='@', index=False, encoding='utf-8')
            print(str(datetime.datetime.now()) + ": Saved "  + str(final_df.shape[0]) + " thinned occurences for " + taxon_name)
            del final_df
            gc.collect()
        else:
            print(str(datetime.datetime.now()) + ": Did not save any thinned occurences for " + taxon_name + " because too few are left (only " + str(len(final_df)) + ")")
    else:
        print(str(datetime.datetime.now()) + ": Did not start thinning for " + taxon_name + " because it has too few cleaned occurrences (only " + str(len(gdf)) + ")")
    return


thin_occ(taxon_key,taxon_name)

