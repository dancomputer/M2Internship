import os
import pandas as pd
import xarray as xr
import numpy as np

def process_netcdf_data(folder_path):
    """
    Processes NetCDF files in a folder, extracts data, and creates a Pandas DataFrame.

    Args:
        folder_path: The path to the folder containing the .nc files.

    Returns:
        A Pandas DataFrame containing the processed data.  Returns an empty dataframe if no .nc files are found.
    """
    all_files = os.listdir(folder_path)
    all_files = [file for file in all_files if file[-3:] == ".nc"] #filter only .nc files
    print(all_files)
    if not all_files:
        print("No .nc files found in the specified folder.")
        return pd.DataFrame() # Return an empty dataframe

    data = []

    for file_path in all_files:
        try:
            # 1. Load the NetCDF file using xarray
            ds = xr.open_dataset(folder_path+"/"+file_path,decode_times=False)
            #print(ds)
            # 2. Get the yield variable (first variable)
            variable_names = list(ds.data_vars)
            if not variable_names:
                print(f"Warning: No data variables found in {file_path}. Skipping.")
                #ds.close()
                #continue  # Skip this file
            
            yield_variable_name = variable_names[0]
            yield_data = ds[yield_variable_name]
            #print(yield_variable_name)

            # 3. & 4. Parse the filename
            filename = os.path.basename(file_path)
            parts = filename[:-3].split("_")  # Remove '.nc' and split

            # Handle potential errors in filename parsing robustly.
            if len(parts) != 10:
                print(f"Warning: Unexpected filename format for {filename}. Skipping.")
                #ds.close()
                #continue

            model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, variable_crop_irrigation, region, time_step, start_year, end_year = parts
            variable, crop_irrigation = variable_crop_irrigation.split("-", 1) #split just in two on first occurrence
            crop, irrigation = crop_irrigation.rsplit("-",1) # Split only the last occurrence

            #print(crop)
            # 5. Append to the data list
            data.append({
                'model': model,
                'climate_forcing': climate_forcing,
                'climate_scenario': climate_scenario,
                'soc_scenario': soc_scenario,
                'sens_scenario': sens_scenario,
                'variable': variable,
                'crop': crop,
                'irrigation': irrigation,
                'region': region,
                'time_step': time_step,
                'start_year': start_year,
                'end_year': 2016, #end_year,
                'yield_data': yield_data,  # Store as xarray object
            })
            #ds.close() #close after reading

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            # Consider if you want to continue or raise the exception.
            # Here, we continue to the next file.

    # 0. Create the Pandas DataFrame
    df = pd.DataFrame(data)
    #Create a median model
    median_df = df.groupby(['climate_forcing', 'climate_scenario', 'soc_scenario',
                    'sens_scenario', 'variable', 'crop', 'irrigation', 'region',
                    'time_step','start_year','end_year']).agg(
    yield_data=('yield_data', lambda x: xr.concat(x, dim="model").median(dim="model"))).reset_index()
    
    # Add the 'model' column with the value 'median'
    median_df['model'] = 'median'
    # Reorder columns to match the original DataFrame's structure
    median_df = median_df[['model', 'climate_forcing', 'climate_scenario', 'soc_scenario',
                           'sens_scenario', 'variable', 'crop', 'irrigation', 'region',
                           'time_step', 'start_year', 'end_year', 'yield_data']]
    final_df = pd.concat([df, median_df], ignore_index=True)

    return final_df

#Load and process data from the ISIMIP3b-YieldModels folder
folder_path = r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\ISIMIP3a-Yield"  
ISIMIP_Yields_df = process_netcdf_data(folder_path)


###
#Load SPAM2020, observational data on harvested areas and production
###

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
            data.append({
                'technology': technology,
                'data': ds,
                'year': file[:-4].split("_")[0][4:8]
            })
    df = pd.DataFrame(data)
    return df

#Load the harvested area and production data

SPAM_MAIZ_production = process_SPAM_data(r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2000V3r107_global_production-geotiff",crop="MAIZ")
SPAM_MAIZ_harvested_area = process_SPAM_data(r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2000V3r107_global_harvested_area-geotiff",crop="MAIZ")
SPAM_MAIZ_production_2005 = process_SPAM_data(r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2005V3r2_global_production-geotiff",crop="MAIZ")
SPAM_MAIZ_harvested_area_2005 = process_SPAM_data(r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\Global_Geotiff\spam2005V3r2_global_harvested_area-geotiff",crop="MAIZ")

##Compare the difference in production of SPAM vs ISIMIP for each model for the year 2000

ISIMIP_prod_R = {model:np.zeros(SPAM_MAIZ_production.iloc[0][1].shape) for model in ISIMIP_Yields_df.model.unique()}
ISIMIP_prod_I = {model:np.zeros(SPAM_MAIZ_production.iloc[0][1].shape) for model in ISIMIP_Yields_df.model.unique()}
ISIMIP_diff_R = {model:np.zeros(SPAM_MAIZ_production.iloc[0][1].shape) for model in ISIMIP_Yields_df.model.unique()}
ISIMIP_diff_I = {model:np.zeros(SPAM_MAIZ_production.iloc[0][1].shape) for model in ISIMIP_Yields_df.model.unique()}

for i in range(len(ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")])):
    #[]
    #if ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")][i]['irrigation'] == "firr":

    if ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")]['irrigation'].loc[i] == "noirr":
        year2000_yield =ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")]['yield_data'].iloc[i].sel(time= 2000-int(ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")]['start_year'].loc[i])).values
        print(len(year2000_yield))    
        year2000_area = SPAM_MAIZ_harvested_area['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "R")  )[0][0]].values
        year2000_SPAM_production = SPAM_MAIZ_production['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "R")  )[0][0]].values
        year2000_ISIMIP_production = year2000_yield * year2000_area
       
        ISIMIP_prod_R [ISIMIP_Yields_df.model.iloc[i]] =  year2000_ISIMIP_production
        ISIMIP_diff_R [ISIMIP_Yields_df.model.iloc[i]] = (year2000_ISIMIP_production - year2000_SPAM_production)/year2000_SPAM_production
   
    if ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")]['irrigation'].loc[i] == "firr":
        year2000_yield =ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")]['yield_data'].iloc[i].sel(time= 2000-int(ISIMIP_Yields_df[(ISIMIP_Yields_df['crop'] == "mai")]['start_year'].loc[i])).values
        print(len(year2000_yield))    
        year2000_area = SPAM_MAIZ_harvested_area['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "I")  )[0][0]].values
        year2000_SPAM_production = SPAM_MAIZ_production['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "I")  )[0][0]].values
        year2000_ISIMIP_production = year2000_yield * year2000_area
       
        ISIMIP_prod_I [ISIMIP_Yields_df.model.iloc[i]] = ISIMIP_prod_I [ISIMIP_Yields_df.model.iloc[i]] + year2000_ISIMIP_production
        ISIMIP_diff_I [ISIMIP_Yields_df.model.iloc[i]] = (year2000_ISIMIP_production - year2000_SPAM_production)/year2000_SPAM_production
#  
#Compare results
#
import matplotlib.pyplot as plt
plt.plot([ np.nansum(ISIMIP_diff_I[key]) for key in list(ISIMIP_diff_I.keys())] )
plt.plot([ np.nanstd(ISIMIP_diff_R[key]) for key in list(ISIMIP_diff_R.keys())] )
for model in list(ISIMIP_diff_R.keys()):
    plt.imshow(ISIMIP_diff_R[model].clip(-1,1),cmap='RdYlGn',vmin=-1,vmax=1 )
    plt.colorbar()
    plt.title(model+" Rainfed, % Diff with SPAM")
    plt.show()

#for model in list(ISIMIP_diff_I.keys()):
#    plt.imshow(ISIMIP_diff_I[model].clip(-1,1),cmap='RdYlGn',vmin=-1,vmax=1 )
#    plt.colorbar()
#    plt.title(model+" Irrigated, % Diff with SPAM")
#    plt.show()

#
#Trying option D
#
# a) compute difference in yields 2000 vs 2005 from SPAM data
# b) compute difference in yields 2000 vs 2005 from ISIMIP data
# c) compute difference in differences

# a) compute difference in yields 2000 vs 2005 from SPAM data
year2000_SPAM_yield = SPAM_MAIZ_production['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "R")  )[0][0]].values[0:308,:]/SPAM_MAIZ_harvested_area['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "R")  )[0][0]].values[0:308,:]
year2005_SPAM_yield = SPAM_MAIZ_production_2005['data'].iloc[np.where( (SPAM_MAIZ_harvested_area_2005['year'].astype(int)==2005) & (SPAM_MAIZ_harvested_area_2005['technology'] == "R")  )[0][0]].values/SPAM_MAIZ_harvested_area_2005['data'].iloc[np.where( (SPAM_MAIZ_harvested_area_2005['year'].astype(int)==2005) & (SPAM_MAIZ_harvested_area_2005['technology'] == "R")  )[0][0]].values
SPAM_diff_R = (year2005_SPAM_yield - year2000_SPAM_yield)/year2000_SPAM_yield

year2000_SPAM_yield = SPAM_MAIZ_production['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "I")  )[0][0]].values[0:308,:]/SPAM_MAIZ_harvested_area['data'].iloc[np.where( (SPAM_MAIZ_harvested_area['year'].astype(int)==2000) & (SPAM_MAIZ_harvested_area['technology'] == "I")  )[0][0]].values[0:308,:]
year2005_SPAM_yield = SPAM_MAIZ_production_2005['data'].iloc[np.where( (SPAM_MAIZ_harvested_area_2005['year'].astype(int)==2005) & (SPAM_MAIZ_harvested_area_2005['technology'] == "I")  )[0][0]].values/SPAM_MAIZ_harvested_area_2005['data'].iloc[np.where( (SPAM_MAIZ_harvested_area_2005['year'].astype(int)==2005) & (SPAM_MAIZ_harvested_area_2005['technology'] == "I")  )[0][0]].values
SPAM_diff_I = (year2005_SPAM_yield - year2000_SPAM_yield)/year2000_SPAM_yield

# b) compute difference in yields 2000 vs 2005 from ISIMIP data
ISIMIP_diff_R = {model:  (ISIMIP_Yields_df.iloc[np.where( 
    (ISIMIP_Yields_df.model==model) & 
    (ISIMIP_Yields_df.irrigation=="noirr") & 
    (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
    ]["yield_data"].sel(time= 
            2005 - int(ISIMIP_Yields_df.iloc[np.where( 
            (ISIMIP_Yields_df.model==model) & 
            (ISIMIP_Yields_df.irrigation=="noirr") & 
            (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
            ]["start_year"])) 
        - 
        ISIMIP_Yields_df.iloc[np.where( 
        (ISIMIP_Yields_df.model==model) & 
        (ISIMIP_Yields_df.irrigation=="noirr") & 
        (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
        ]["yield_data"].sel(time= 
            2000 - int(ISIMIP_Yields_df.iloc[np.where( 
            (ISIMIP_Yields_df.model==model) & 
            (ISIMIP_Yields_df.irrigation=="noirr") & 
            (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
            ]["start_year"])))
        /
        ISIMIP_Yields_df.iloc[np.where( 
        (ISIMIP_Yields_df.model==model) & 
        (ISIMIP_Yields_df.irrigation=="noirr") & 
        (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
            ]["yield_data"].sel(time= 
                2000 - int(ISIMIP_Yields_df.iloc[np.where( 
                (ISIMIP_Yields_df.model==model) & 
                (ISIMIP_Yields_df.irrigation=="noirr") & 
                (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
                ]["start_year"]))
 for model in ISIMIP_Yields_df.model.unique()}

ISIMIP_diff_I = {model:  (ISIMIP_Yields_df.iloc[np.where( 
    (ISIMIP_Yields_df.model==model) & 
    (ISIMIP_Yields_df.irrigation=="firr") & 
    (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
    ]["yield_data"].sel(time= 
        2005 - int(ISIMIP_Yields_df.iloc[np.where( 
        (ISIMIP_Yields_df.model==model) & 
        (ISIMIP_Yields_df.irrigation=="firr") & 
        (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
        ]["start_year"])) - ISIMIP_Yields_df.iloc[np.where( 
    (ISIMIP_Yields_df.model==model) & 
    (ISIMIP_Yields_df.irrigation=="firr") & 
    (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
    ]["yield_data"].sel(time= 
        2000 - int(ISIMIP_Yields_df.iloc[np.where( 
        (ISIMIP_Yields_df.model==model) & 
        (ISIMIP_Yields_df.irrigation=="firr") & 
        (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
        ]["start_year"])))/ISIMIP_Yields_df.iloc[np.where( 
    (ISIMIP_Yields_df.model==model) & 
    (ISIMIP_Yields_df.irrigation=="firr") & 
    (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
    ]["yield_data"].sel(time= 
        2000 - int(ISIMIP_Yields_df.iloc[np.where( 
        (ISIMIP_Yields_df.model==model) & 
        (ISIMIP_Yields_df.irrigation=="firr") & 
        (ISIMIP_Yields_df['crop'] == "mai"))[0][0]
        ]["start_year"]))
 for model in ISIMIP_Yields_df.model.unique()}

# c) compute difference in differences
for model in list(ISIMIP_diff_R.keys()):
    plt.imshow(ISIMIP_diff_R[model][:308,:] - SPAM_diff_R ,cmap='RdYlGn',vmin=-1,vmax=1 )
    plt.colorbar()
    plt.title(model+" Rainfed, % Diff year 2000->2005 yields, vs SPAM")
    plt.show()

for model in list(ISIMIP_diff_I.keys()):
    plt.imshow(ISIMIP_diff_I[model].clip(-1,1),cmap='RdYlGn',vmin=-1,vmax=1 )
    plt.colorbar()
    plt.title(model+" Irrigated, % Diff year 2000->2005 yields, vs SPAM")
    plt.show()
    
