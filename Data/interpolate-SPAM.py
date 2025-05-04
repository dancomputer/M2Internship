'''
Create interpolated harvested area dataset from SPAM2020.
'''

def process_SPAM_data(folder_path,crop):
    '''
    Upscales and processes SPAM data. From a folder, extracts data, and creates a Pandas DataFrame.'''

    all_files = os.listdir(folder_path)
    all_files = [file for file in all_files if file[-3:] == "tif"] #filter only files with the specified crop
    all_files = [file for file in all_files if file[:-4].split("_")[-2] == crop] #filter only files with the specified crop
    data = []
    for file in all_files:
        technology = file[:-4].split("_")[-1]
        if technology == "R" or technology == "I":
            ds = xr.open_dataset(folder_path+"/"+file)
            ds = ds['band_data'].drop("band").sel(band=0)
            ds = ds.coarsen(x=6, y=6, boundary='trim').sum()
            ds['x'] = np.round(ds['x'].values,2)
            ds['y'] = np.round(ds['y'].values,2)
            data.append({
                'technology': technology,
                'data': ds,
                'year': file[:-4].split("_")[0][4:8]
            })
    df = pd.DataFrame(data)
    return df

#Load the harvested area and production data
#for the 2000 dataset, it needs to be processed sliced to [0:308,:]
#for the rest of the datsets they should first all be loaded together, for rice: 
#the file paths are as follwos
#1. r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2000V3r107_global_harvested_area-geotiff"
#2. r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2005V3r2_global_harvested_area-geotiff"
#3. r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2020V1r0_global_harvested_area-geotiff" 
#
#Next should come an interpolation to a yearly timescale, i.e. bridign the gaps between the years of downlaoded data linearly.
# 

import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def interpolate_SPAM_data(df):
    """
    Interpolates SPAM data linearly across missing years.
    """
    years = sorted(df['year'].unique())
    technologies = df['technology'].unique()
    
    interpolated_data = []
    
    for tech in technologies:
        df_tech = df[df['technology'] == tech]
        ds_list = []
        for y in years:
            ds = df_tech[df_tech['year'] == y]['data'].values[0]
            ds = ds.expand_dims({'year':[int(y)]})
            ds_list.append(ds)#[:,:309,:])
        interpolated = xr.concat(ds_list, dim="year", coords="minimal")
        interpolated = interpolated.interp(year=np.arange(int(years[0]), int(years[-1]) + 1))
        interpolated = interpolated.sortby("y", ascending=False)
        interpolated_data.append({
                'technology': tech,
                'data': interpolated
            })
    
    return pd.DataFrame(interpolated_data)

# Load and interpolate the data
folders = {
    2000: r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2000V3r107_global_harvested_area-geotiff",
    2005: r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2005V3r2_global_harvested_area-geotiff",
    #2010: r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2010V2_global_harvested_area-geotiff",
    2020: r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2020V1r0_global_harvested_area-geotiff"
}

dataframes = [process_SPAM_data(folder, crop="MAIZ") for folder in folders.values()]
combined_df = pd.concat(dataframes)
interpolated_df = interpolate_SPAM_data(combined_df)

#
#Next, check the data to see if its reasonable.
#

#1.
lats,lons = np.where((interpolated_df.iloc[0].data.x >= -20) & (interpolated_df.iloc[0].data.x <= 55) & (interpolated_df.iloc[0].data.y >= -35) & (interpolated_df.iloc[0].data.y <= 15))
print("change in rainfed area 2000-2020, SSA: " + str(round(np.nansum(interpolated_df[interpolated_df.technology=='R'].iloc[0].data.values[-1,lons[0]:lons[-1],lats[0]:lats[-1]])/np.nansum(interpolated_df[interpolated_df.technology=='R'].iloc[0].data.values[0,lons[0]:lons[-1],lats[0]:lats[-1]]),2)))
for i in len(interpolated_df):
    ds = interpolated_df.iloc[i].data
    ds.to_netcdf("spam_2000-2020_linearinterpolation_global_harvestedarea_"+"MAIZE_"+interpolated_df.iloc[i].technology+".nc")
