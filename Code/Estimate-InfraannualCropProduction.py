import os
import pandas as pd
import xarray as xr


def process_netcdf_data(yield_path, maturity_path):
    """
    Processes NetCDF files in a folder, extracts yield and maturity data, and creates a Pandas DataFrame.

    Args:
        yield_path: The path to the folder containing the yield .nc files.
        maturity_path: The path to the folder containing the maturity .nc files.

    Returns:
        A Pandas DataFrame containing the processed data. Returns an empty DataFrame if no .nc files are found.
    """
    yield_files = [f for f in os.listdir(yield_path) if f.endswith(".nc")]
    maturity_files = [f for f in os.listdir(maturity_path) if f.endswith(".nc")]

    if not yield_files:
        print("No .nc files found in the specified yield folder.")
        return pd.DataFrame()  # Return an empty DataFrame

    # Create a dictionary to store parsed metadata of maturity files
    maturity_dict = {}

    for file_name in maturity_files:
        parts = file_name[:-3].split("_")  # Remove .nc and split
        if len(parts) != 10:
            continue  # Skip files with incorrect formatting

        model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, variable_crop_irrigation, region, time_step, start_year, end_year = parts
        variable, crop_irrigation = variable_crop_irrigation.split("-", 1)
        crop, irrigation = crop_irrigation.rsplit("-", 1)

        # Store the file path using parsed metadata as key
        key = (model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, crop, irrigation, region, time_step, start_year, end_year)
        maturity_dict[key] = os.path.join(maturity_path, file_name)

    data = []

    for file_name in yield_files:
        try:
            yield_file_path = os.path.join(yield_path, file_name)
            ds_yield = xr.open_dataset(yield_file_path, decode_times=False)

            # Extract yield variable
            variable_names = list(ds_yield.data_vars)
            if not variable_names:
                print(f"Warning: No data variables found in {file_name}. Skipping.")
                continue

            yield_variable_name = variable_names[0]
            yield_data = ds_yield[yield_variable_name]

            # Parse filename
            parts = file_name[:-3].split("_")
            if len(parts) != 10:
                print(f"Warning: Unexpected filename format for {file_name}. Skipping.")
                continue

            model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, variable_crop_irrigation, region, time_step, start_year, end_year = parts
            variable, crop_irrigation = variable_crop_irrigation.split("-", 1)
            crop, irrigation = crop_irrigation.rsplit("-", 1)

            # Try to find the matching maturity file
            key = (model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, crop, irrigation, region, time_step, start_year, end_year)
            maturity_data = None

            if key in maturity_dict:
                ds_maturity = xr.open_dataset(maturity_dict[key], decode_times=False)
                maturity_variable_name = list(ds_maturity.data_vars)[0]  # Assuming the first variable is correct
                maturity_data = ds_maturity[maturity_variable_name]
                ds_maturity.close()

            # Append to data list
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
                'end_year': 2016,
                'yield_data': yield_data,
                'maturity_data': maturity_data  # Include maturity data
            })

            ds_yield.close()

        except Exception as e:
            print(f"Error processing {file_name}: {e}")

    # Create DataFrame
    df = pd.DataFrame(data)

    median_df = df.groupby([
        'climate_forcing', 'climate_scenario', 'soc_scenario', 'sens_scenario',
        'variable', 'crop', 'irrigation', 'region', 'time_step', 'start_year', 'end_year'
    ]).agg(
        yield_data=('yield_data', lambda x: xr.concat(x, dim="model").median(dim="model")),
        maturity_data=('maturity_data', lambda x: xr.concat([d for d in x if d is not None], dim="model").median(dim="model") 
                    if any(d is not None for d in x) else None)
    ).reset_index()

    # Add median model entry
    median_df['model'] = 'median'

    final_df = pd.concat([df, median_df], ignore_index=True)

    return final_df

#Load and process data from the ISIMIP3b-YieldModels folder
yield_path = r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\ISIMIP3a-Yield"  
maturity_path = r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\ISIMIP3a-MatyDay"  
ISIMIP_Yields_df = process_netcdf_data(yield_path, maturity_path)

    yield_files = [f for f in os.listdir(yield_path) if f.endswith(".nc")]
    maturity_files = [f for f in os.listdir(maturity_path) if f.endswith(".nc")]

    # Create a dictionary to store parsed metadata of maturity files
    maturity_dict = {}

    for file_name in maturity_files:
        parts = file_name[:-3].split("_")  # Remove .nc and split
        if len(parts) != 10:
            continue  # Skip files with incorrect formatting

        model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, variable_crop_irrigation, region, time_step, start_year, end_year = parts
        variable, crop_irrigation = variable_crop_irrigation.split("-", 1)
        crop, irrigation = crop_irrigation.rsplit("-", 1)

        # Store the file path using parsed metadata as key
        key = (model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, crop, irrigation, region, time_step, start_year, end_year)
        maturity_dict[key] = os.path.join(maturity_path, file_name)
    data = []
    for file_name in yield_files:

            yield_file_path = os.path.join(yield_path, file_name)
            ds_yield = xr.open_dataset(yield_file_path, decode_times=False)

            # Extract yield variable
            variable_names = list(ds_yield.data_vars)
            if not variable_names:
                print(f"Warning: No data variables found in {file_name}. Skipping.")

            yield_variable_name = variable_names[0]
            yield_data = ds_yield[yield_variable_name]

            # Parse filename
            parts = file_name[:-3].split("_")
            if len(parts) != 10:
                print(f"Warning: Unexpected filename format for {file_name}. Skipping.")

            model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, variable_crop_irrigation, region, time_step, start_year, end_year = parts
            variable, crop_irrigation = variable_crop_irrigation.split("-", 1)
            crop, irrigation = crop_irrigation.rsplit("-", 1)

            # Try to find the matching maturity file
            key = (model, climate_forcing, climate_scenario, soc_scenario, sens_scenario, crop, irrigation, region, time_step, start_year, end_year)
            maturity_data = None

            if key in maturity_dict:
                ds_maturity = xr.open_dataset(maturity_dict[key], decode_times=False)
                maturity_variable_name = list(ds_maturity.data_vars)[0]  # Assuming the first variable is correct
                maturity_data = ds_maturity[maturity_variable_name]
                ds_maturity.close()

import numpy as np

from dask.distributed import Client
import dask.delayed
import time
import pandas as pd 
import dask
from dask.distributed import LocalCluster
#os.chdir("/tank/crop_modelling/Daniel/NecessaryM1InternshipCode/ProjectRice")

#start cluster
cluster = LocalCluster(n_workers=os.cpu_count() - 5 ,threads_per_worker=1)          # Fully-featured local Dask cluster: forward port 8787 to via vscode to see dashboard
client = cluster.get_client()


pixel_list = np.arange(0,len(lons))    
lat_index = xr.DataArray(lats, dims="z")
lon_index = xr.DataArray(lons, dims="z")
tsmax_selected=tsmax_ds.isel(lat=lat_index,lon=lon_index)

#
#start_index= [yearstart+plantingdate_selected.values[:,i].astype(int)- 1 + 0 for i  in range(len(lons)) ] #from day 0 of the season
#end_index = [yearstart+plantingdate_selected.values[:,i].astype(int)- 1 + 0 + n_days for i in range(len(lons))] #compute up to 400 days after. 


results = []

#
#First calculate daily cumulative GDDs, then the average GDDs to estimate cultivar PHUs, and divide to get the "% of growing season" variable.
#

for i in pixel_list:
    #print(i)
    res = calculate_GDD(tsmaxi = tsmax_selected[:,i], tsmini = tsmin_selected[:,i], start_index=start_index[i], end_index=end_index[i], t_max=t_max, t_min=t_min)
    results.append(res)
#compute
#cum_gdd_result = [dask.compute(computation) for computation in results]
cum_gdd_results = client.compute(results)
cum_GDD = np.array([ress.result() for ress in cum_gdd_results])


"""
Generate a 365-day time series based on a normal distribution centered on maturity day.

Args:
    yield_data (xarray.DataArray): The yield for the year.
    maturity_day (xarray.DataArray): The day of year when the crop matures.

Returns:
    xarray.DataArray: A 365-day time series for the year.
"""
days = np.arange(0, 365)  # Days of the year

# Standard deviation to achieve ~20-day spread
std_dev = 20 / 2  # 95% confidence within ~20 days
lats,lons = np.where(np.isnan(maturity_day).sum(axis=0)<maturity_day.shape[0])
maturitys = np.array([maturity_day[:,i,j] for i,j in zip(lats,lons)])
yields = np.array([yield_data.values[:,i,j] for i,j in zip(lats,lons)])
results = []

from scipy.stats import norm

#@dask.delayed

def custom_normpdf(x,loc,scale):
    result = norm.pdf(x, loc=loc, scale=scale)
    if ~np.isnan(loc):
        min = np.minimum(10,loc).astype(int)
        max = np.minimum(10,365-loc).astype(int)
        result[0:loc-min] = 0
        result[loc+max:] = 0
        result /= np.sum(result)
    return result
    
@dask.delayed
def temporallydownscale_cropyield(maturity_data,yield_data):
    results = np.zeros((365,len(yield_data)),dtype=int)
    days = np.arange(1, 366)  # Days of the year
    results = [ norm.pdf(days, loc=maturity_day, scale=10)*yields for maturity_day,yields in zip(maturity_data,yield_data)]

    return np.array(results)

start =time.time()
for pixel in pixel_list[1:100]:

    results.append(temporallydownscale_cropyield(maturity_data=maturitys[pixel,:],yield_data=yields[pixel,:]))
end =time.time()
output = client.compute(results)
final_outputs = np.array([ress.result() for ress in output])
print(start-end)

weights /= weights.sum()  # Normalize to ensure area = 1

# Multiply by yield to distribute it over 365 days
time_series = yield_data * weights