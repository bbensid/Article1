###make sure to use the final definition of M with positive definite covariant matrix !
###global PCA -> Get M19 and K
###specific PCA -> Get M19, K, PCrel (with the permutations), MPCrel
###

import glob, pandas as pd, pickle as pk, datetime, gc, numpy as np, scipy as sp, os
import matplotlib.pyplot as plt, seaborn as sns
from shapely.geometry.polygon import Polygon
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
from scipy.stats import chi2


#######################################################################################################
#######################################SET SCRIPT VARIABLES############################################
#######################################################################################################

#Specific to running environment 
bioclim_occ_folder_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/3-Bioclim/'
final_file_path = '/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/FINAL.txt'

#######################################################################################################
############################################GET ALL DATA#####################################
#######################################################################################################

bioclim_files = glob.glob(bioclim_occ_folder_path +'*.txt')

df_list = []

for x in bioclim_files:
	if os.path.getsize(x) <= 3000000: #only files less than 3MB for computation capacity reasons 
		df_list.append(pd.read_csv(x, sep='@', encoding='utf-8'))

df_all=pd.concat(df_list).reset_index(drop=True) #  ==> 276 (78 for exceeding file size limit of 3MB, 58 for being under occurrences limit of 19)

all_taxons = list(set(df_all.Taxon))

#####to delete after manual run

a = pd.read_csv(final_file_path).values.tolist()
already_run = [(x[0],x[2]) for x in a]

# Get all possible pairs of taxa =>75'900
pairs = []
for taxon1 in all_taxons:
    for taxon2 in all_taxons:
        if taxon1 != taxon2 and (taxon1,taxon2) not in already_run:
           pairs.append((taxon1,taxon2)) 

list_to_iterate = pairs if len(pairs) < 100 else pairs[:100]

############################################
################M AND K FUNCTIONS###########
############################################

##ref: https://gis.stackexchange.com/questions/99917/converting-matplotlib-contour-objects-to-shapely-objects
def get_poly(df,thresh):
    fig,ax = plt.subplots()
    kde = sns.kdeplot(
        data = df
        ,x='principal component 1'
        ,y='principal component 2'
        ,levels = 1
        ,thresh=thresh
        )
    col = kde.collections[0]
    paths = []
    # Loop through all polygons that have the same intensity level
    for contour in col.get_paths(): 
        if len(contour)==1:
        ## Contour is just one point, skipping##
            continue
        # Create a polygon for the countour
        # First polygon is the main countour, the rest are holes if any (should not happen as those polygons are from kernel-density based contours)
        for ncp,cp in enumerate(contour.to_polygons()):
            if type(cp) is not np.ndarray:
                cp = np.array(cp)
            x = cp[:,0]
            y = cp[:,1]
            new_shape = Polygon([(i[0], i[1]) for i in zip(x,y)])
            if ncp == 0:
                poly = new_shape
            else:
                # Remove holes, if any (should not happen as those polygons from kernel-density based contours)
                poly = poly.difference(new_shape)
        # Append polygon to list
        paths.append(poly)
    plt.close()
    del kde
    gc.collect()
    return paths

def removing_holes(p_list): 
    ## This function gets a list of polyon pairs, and gives a new list of polygons with actual holes in them ##
    ## when one polygon is inside another in the original list ##
    if len(p_list)==0:
        ##list empty,returning None##
        return None
    elif len(p_list)==1:
        ##return the single polygon in case there is only one polygon##
        return p_list
    else:
        area_list=[]
        for p in p_list:
            area_list.append(p.area)
        ##Create an ordered list of the polygons by decreasing area##
        order_p=[]
        p_list_copy=p_list.copy()
        area_list_copy=area_list.copy()
        for i in range(len(p_list)):
            order_p.append(p_list_copy[area_list_copy.index(max(area_list_copy))])
            p_list_copy.pop(area_list_copy.index(max(area_list_copy)))
            area_list_copy.pop(area_list_copy.index(max(area_list_copy)))
        ##checking if any of the smaller polygons are contained in the bigger polygon,in that case it is properly a hole and it is removed##
        ##this way takes also care of the (rare) cases with "russian dolls" polygons, where a holes have actual areas in them that should be counted as the niche area##
        ##by recursively removing the holes, hence the subsequent "plain" contour inside the hole is not anymore contained in the new polygon with the hole##
        final_list=[]
        to_remove=[]
        for p in order_p:
            if order_p.index(p)<len(order_p)-1 and p not in to_remove:
                result_polygon=p
                for smaller_p in order_p[order_p.index(p)+1:]:
                    if result_polygon.contains(smaller_p):
                        to_remove.append(smaller_p)
                        result_polygon=result_polygon.difference(smaller_p)
                final_list.append(result_polygon)
            elif order_p.index(p)==len(order_p)-1 and p not in to_remove:
                final_list.append(p)
        return final_list


def get_K(p1_list,p2_list):
    if p1_list==[] or p2_list==[]:
            ## Exiting with None, as both polygon list should be non-empty ##
            return None
    else:
        p1_list_new = removing_holes(p1_list)
        p2_list_new = removing_holes(p2_list)
        p1_area = 0
        intersection_area = 0
        for p1 in p1_list_new:
            p1_area += p1.area
            for p2 in p2_list_new:
                intersection_area += p1.intersection(p2).area
        ## the rounding to 10 decimals is because, the intersection aera is a sum over calculation that have a certain error,this adds uo andcan get over 100%##
        ##usually at the 16th decimal, so this is really just not to have that###
        return round(intersection_area/p1_area,10)

def chunking_dot(big_matrix, small_matrix, chunk_size=5000):
    # Make a copy if the array is not already contiguous
    small_matrix = np.ascontiguousarray(small_matrix)
    diag = np.empty((big_matrix.shape[0],))
    for j in range(0, big_matrix.shape[0], chunk_size):
        end = j + chunk_size
        f = np.dot(big_matrix[j:end], small_matrix)
        diag[j:end] = f.diagonal()
        del f
        gc.collect()
    return diag


def mahalanobis(x=None, data=None, cov=None):
    """Compute the Mahalanobis Distance between each row of x and the data  
    x    : vector or matrix of data with, say, p columns.
    data : ndarray of the distribution from which Mahalanobis distance of each observation of x is to be computed.
    cov  : covariance matrix (p x p) of the distribution. If None, will be computed from data.
    """
    x_minus_mu = x - np.mean(data)
    if not cov:
        cov = np.cov(data.values.T)
    # Making sure the matric is Definite Positive, else return None
    if np.all(np.linalg.eigvals(cov) > 0):
        inv_covmat = np.linalg.pinv(cov)
        del cov
        gc.collect()
        if np.all(np.linalg.eigvals(inv_covmat) > 0):
            left_term = np.dot(x_minus_mu, inv_covmat)
            trans = x_minus_mu.T.to_numpy()
            del inv_covmat,x_minus_mu
            gc.collect()
            mahal = chunking_dot(left_term, trans)
            del left_term
            gc.collect()
            return mahal
    else:
        return None

def get_M(df1,df2,n):
    maha_df = df1.copy()
    if mahalanobis(x=df1,data=df2) is None:
        return None
    else:
        maha_df['mahala']=mahalanobis(x=df1,data=df2)
        maha_df['p_value'] = 1 - chi2.cdf(maha_df['mahala'], n-1)
        maha_df['outlier'] = maha_df['p_value'] <0.05
        overlap=1-sum(maha_df['outlier'])/len(maha_df['p_value'])
        del maha_df
        gc.collect()
        return overlap

####################GLOBAL PCA########################


# Run the global PCA on all the taxon list

df_features = df_all.drop(columns=['Taxon'])
# Standardizing the features
scaler = MinMaxScaler(feature_range=(0, 1))
x = scaler.fit_transform(df_features)

# Run the PCA
pca = PCA()
principalComponents = pca.fit_transform(x)

# Save global Standard Scaler and PCA

pk.dump(scaler, open("/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/4-StandardScaler/global_minmaxscaler_"+ datetime.date.today().strftime("%d%m%Y") + ".pkl" ,"wb"))

pk.dump(pca, open("/Users/home/Documents/PhD/GBIF_source/Bioclim/Occurrences/5-PCA/global_pca_"+ datetime.date.today().strftime("%d%m%Y") + ".pkl" ,"wb"))



all_list = ['principal component 1', 'principal component 2', 'principal component 3', 'principal component 4', 'principal component 5', 'principal component 6', 'principal component 7', 'principal component 8', 'principal component 9', 'principal component 10', 'principal component 11', 'principal component 12', 'principal component 13', 'principal component 14', 'principal component 15', 'principal component 16', 'principal component 17', 'principal component 18', 'principal component 19']
principalDf = pd.DataFrame(data = principalComponents, columns = all_list)


principalDf['tag'] = df_all.Taxon

del principalComponents,df_features
gc.collect()

##########SPECIFIC PCA######

def de_correlate_df(df):
	X_aux = df.copy()
	for col in df.columns:
		X_aux[col] = df[col].sample(len(df)).values
	return X_aux

def get_pair_pca(taxon1,taxon2):
	df_pair  = df_all[(df_all.Taxon==taxon1)|(df_all.Taxon==taxon2)].reset_index(drop=True)
	df_features = df_pair.drop(columns=['Taxon'])
	# Standardizing the features
	scaler = MinMaxScaler(feature_range=(0, 1))
	x = scaler.fit_transform(df_features)
	pca = PCA()
	principalComponents = pca.fit_transform(x)
	pair_principalDf = pd.DataFrame(data = principalComponents, columns = all_list)
	pair_principalDf['tag'] = df_pair.Taxon
	d1=pair_principalDf[pair_principalDf.tag==taxon1].reset_index(drop=True).drop(columns=['tag'])
	d2=pair_principalDf[pair_principalDf.tag==taxon2].reset_index(drop=True).drop(columns=['tag'])
	del principalComponents,df_pair
	gc.collect()
	# Do the permutations so we can get the number of relevant Principal Components
	original_variance = pca.explained_variance_ratio_
	N_permutations = 100
	variance = np.zeros((N_permutations, len(df_features.columns)))
	for i in range(N_permutations):
		X_aux = de_correlate_df(df_features)
		x = scaler.fit_transform(X_aux)
		pca = PCA()
		principalComponents = pca.fit_transform(x)
		variance[i, :] = pca.explained_variance_ratio_
	p_val = np.sum(variance > original_variance, axis=0) / N_permutations
	PC = sum(map(lambda x : x == 0, p_val.tolist()))
	PC = PC if PC > 1 else 2 # Taking the first 2 even if only the first is the only relevant one
	del principalComponents,df_features
	gc.collect()
	return d1,d2, PC

all_rows = []

i=0

for pair in list_to_iterate:
	print(str(datetime.datetime.now()) + ": " + str(i) + ": Starting with pair: " + str(pair))
	df1=principalDf[principalDf.tag==pair[0]].reset_index(drop=True).drop(columns=['tag'])
	df2=principalDf[principalDf.tag==pair[1]].reset_index(drop=True).drop(columns=['tag'])
	p1_list = get_poly(df1[all_list[:2]],thresh=0.05)
	p2_list = get_poly(df2[all_list[:2]],thresh=0.05)
	d1, d2, PC =get_pair_pca(pair[0],pair[1])
	p1_list_pair = get_poly(d1[all_list[:2]],thresh=0.05)
	p2_list_pair = get_poly(d2[all_list[:2]],thresh=0.05)
	all_rows.append([pair[0],len(df1),pair[1],len(df2),get_M(df1,df2,19),get_K(p1_list,p2_list),get_M(d1,d2,19),get_K(p1_list_pair,p2_list_pair),PC,get_M(d1[all_list[:PC]],d2[all_list[:PC]],PC)])
	print(str(datetime.datetime.now()) + ": " + str(i) + ": Finished with pair: " + pair[0] + " (" + str(len(df1)) +") and " + pair[1] + " (" + str(len(df2)) +")")
	i+=1
	# Save for every 10 done (pre-created csv file with one row containing comma separated column names)
	if i%10 == 0 :
		all_rows_push = pd.read_csv(final_file_path).values.tolist() + all_rows
		new_rows = pd.DataFrame(all_rows_push,columns=['taxon1','size1','taxon2','size2','M19_global','K_global','M19_specific','K_specific','PC_relevant','M_PC_relevant'])
		new_rows.to_csv(final_file_path,index=False)
		del all_rows_push,new_rows,all_rows
		gc.collect()
		all_rows = []



