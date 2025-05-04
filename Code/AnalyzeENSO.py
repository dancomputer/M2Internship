import xarray as xr
import glob
import re

# Step 1: Find all yield.nc4 files
file_list = sorted(glob.glob(r"C:\Users\danie\OneDrive\Desktop\M2 Internship\Data\AgriculturalProduction\GDHYv12-13\soybean\yield*.nc4"))  # Adjust path

# Step 2: Extract years from filenames
datasets = []
for file in file_list:
    match = re.search(r"yield_(\d{4})\.nc4", file)  # Extract year
    if match:
        year = int(match.group(1))  # Convert to integer
        ds = xr.open_dataset(file)  # Load dataset
        ds = ds.expand_dims(year=[year])  # Add 'year' as a new dimension
        datasets.append(ds)

# Step 3: Concatenate along 'year' dimension
merged_ds = xr.concat(datasets, dim="year")
# left off at making the shape suitable for an M-SSA analysis
#also commented off the nan sum filter because only about 2200 pixels had fully non na data.
merged_ds['var'].values.reshape(36,-1)
import numpy as np
all_nonna_signals = merged_ds['var'].values.reshape(36,-1)[:,np.sum(np.isnan(merged_ds['var'].values.reshape(36,-1)),axis=0)<36]
normed_all_nonna_signals=(all_nonna_signals-np.nanmean(all_nonna_signals,axis=0))/np.nanstd(all_nonna_signals,axis=0)
normed_all_nonna_signals[np.isnan(normed_all_nonna_signals)] = 0
mean_signal = np.nanmean(all_nonna_signals,axis=1)
np.savetxt('GDHY_soyyields.txt',normed_all_nonna_signals)