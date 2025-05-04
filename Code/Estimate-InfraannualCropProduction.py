import numpy as np
import xarray as xr
from dask.distributed import Client
import dask.delayed
import time
import pandas as pd 
import dask
import os
from dask.distributed import LocalCluster
#os.chdir("/tank/crop_modelling/Daniel/NecessaryM1InternshipCode/ProjectRice")

#start cluster
#cluster = LocalCluster(n_workers=os.cpu_count() - 5 ,threads_per_worker=1)          # Fully-featured local Dask cluster: forward port 8787 to via vscode to see dashboard
#client = cluster.get_client()

'''

This code's purpose is to estimate the infra-annual crop production based on the ISIMIP3a model which was selected in ISIMIP3a_SpatialYieldIntercomparison.

The final temporal resolution will be either daily or monthly.

'''

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

def Temporally_Downscale(maturity_data,yield_data):
    """
    Generate a 365-day time series based on a normal distribution centered on maturity day.

    Args:
        yield_data (xarray.DataArray): The yield for the year.
        maturity_day (xarray.DataArray): The day of year when the crop matures.

    Returns:
        xarray.DataArray: A 365-day time series for the year.
    """
 
    # Standard deviation to achieve ~20-day spread
    std_dev = 20 / 2  # 95% confidence within ~20 days
    lats,lons = np.where(np.isnan(maturity_data.values).sum(axis=0)<maturity_data.values.shape[0]) #lats, lons are the indices of the pixels that have no missing values   
    latitudes,longitudes = maturity_data.lat.values[lats],maturity_data.lon.values[lons] #latitudes and longitudes of the pixels that have no missing values
    pixel_list = np.arange(len(lons))
    maturitys = np.array([maturity_data.values[:,i,j] for i,j in zip(lats,lons)]) #use the indices to get the maturity data for the pixels that have no missing values
    yields = np.array([yield_data.values[:,i,j] for i,j in zip(lats,lons)]) #use the indices to get the yield data for the pixels that have no missing values
    results = []

    from scipy.stats import norm

    def custom_normpdf(x,loc,scale):
        '''
        Custom normal pdf function that sets the values of the pdf to 0 before and after the maturity day, then normalizes the pdf.   
        '''
        result = norm.pdf(x, loc=loc, scale=scale)
        if ~np.isnan(loc):
            loc = int(loc)
            min = np.minimum(10,loc).astype(int)
            max = np.minimum(10,365-loc).astype(int)
            print(loc-min)
            result[0:loc-min] = 0
            result[loc+max:] = 0
            #now upscale the daily to monthly pdf
            days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            # Create an index to slice the data based on month lengths
            start_idx = 0
            sum_months = []
            # Loop through each month and sum the data for that month
            for days in days_in_months:
                end_idx = start_idx + days
                sum_months.append(np.sum(result[start_idx:end_idx]))
                start_idx = end_idx
            # Convert the list of sums into a numpy array for convenience
            sum_months = np.array(sum_months)
            sum_months /= np.sum(sum_months)
            return sum_months

    #@dask.delayed
    def temporallydownscale_cropyield(maturity_data,yield_data):
        '''
        Temporally downscales the crop yield data to a daily time series based on the maturity day.
        '''
        #results = np.zeros((12,len(yield_data)),dtype=int)
        days = np.arange(1, 366)  # Days of the year
        results = [ custom_normpdf(days, loc=maturity_day, scale=10)*yields for maturity_day,yields in zip(maturity_data,yield_data)]
        return np.array(results)

    #
    #Compute, either with dask or without dask, the temporally downscaled crop yield data for the pixels in the subset_pixel_list
    #
    BurkinaFaso_pixel_list = pixel_list[ np.where( (latitudes >= 9) & (latitudes <= 15) & (longitudes >= -6) & (longitudes <= 3) ) ]
    #africa_pixel_list = pixel_list[ np.where( (latitudes >= -35) & (latitudes <= 37) & (longitudes >= -25) & (longitudes <= 55) ) ]

    start =time.time()
    for pixel in BurkinaFaso_pixel_list:
        results.append(temporallydownscale_cropyield(maturity_data=maturitys[pixel,-30:-1],yield_data=yields[pixel,-30:-1]))

    #output = client.compute(results)
    #final_outputs = np.array([ress.result() for ress in output])
    final_outputs = np.array(results)
    end = time.time()
    print(start-end)


    ##
    ##Save the temporally downscaled crop yield data as a netcdf file
    ##

    import numpy as np

    # Generate the latitude and longitude grid
    # Define years and days
    years = np.arange(1901, 2016)[-30:-1].astype(int)  # Reduced above earlier. no need for such a long time series.
    origin_year = years[0]  # Or any other year you are using.
    origin_date = pd.to_datetime(f'{origin_year}-01-01')
    months = np.arange(1,12+1,1).astype(int)
    #
    #finally, save as netcdf4 filae
    #
    subset_latitudes = latitudes[BurkinaFaso_pixel_list]
    subset_longitudes = longitudes[BurkinaFaso_pixel_list]
    #create gridded data
    yieldts_gridded = np.empty((len(np.unique(subset_latitudes)),len(np.unique(subset_longitudes)),len(years),len(months)))*np.nan
    for i in range(len(BurkinaFaso_pixel_list)): #length pixel list because that is the spatial dimension collapsed, and we fill time dimensions all at once
        yieldts_gridded[np.where(np.unique(subset_latitudes)==subset_latitudes[i])[0][0],np.where(np.unique(subset_longitudes)==subset_longitudes[i])[0][0],:,:] = final_outputs[i,:,:]
    print('part')
    # Setting up coordinates and dimensions
    coords = {
        'lat': np.unique(subset_latitudes),
        'lon': np.unique(subset_longitudes),
        'year': years,
        'months': months,
    }
    dims = ('lat','lon', 'year', 'months')
    # Creating DataArray
    data_array = xr.DataArray(yieldts_gridded, coords=coords, dims=dims) #multiply by 100 to convert from fraction to percent.
    # Creating Dataset
    dataset = xr.Dataset({
        'yield': data_array
    })
    return dataset
results = Temporally_Downscale(ISIMIP_Yields_df)