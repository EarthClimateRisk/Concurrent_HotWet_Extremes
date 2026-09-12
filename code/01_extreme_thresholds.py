# -*- coding: utf-8 -*-
"""
Extreme-event threshold calculation.

Calculates grid-cell-specific percentile thresholds from daily NetCDF files
of maximum wet-bulb temperature and precipitation. The baseline period can
be specified, allowing the same functions to be used for ERA5-Land and
CMIP6/ISIMIP3b data.

Input:
    Monthly NetCDF files containing daily data of maximum wet-bulb temperature
    and precipitation.
Main thresholds:
    Heatwave: 95th percentile of daily maximum wet-bulb temperature.
    Extreme precipitation: 90th percentile of wet-day precipitation (> 0).

Historical baseline period:
    1981–2010.
"""

import os
import xarray as xr


def calculate_wetbulb_threshold(input_root, output_root, start_year=1981, end_year=2010, variable_name="wetbulb", time_dim="time", file_prefix="wetbulb", percentiles=(95,)):
    """Calculate grid-cell-specific wet-bulb temperature thresholds."""

    os.makedirs(output_root, exist_ok=True)
    all_data = []

    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            month_str = f"{month:02d}"
            input_file = f"{file_prefix}_{year}_{month_str}.nc"
            input_path = os.path.join(input_root, input_file)

            if os.path.exists(input_path):
                try:
                    ds = xr.open_dataset(input_path)
                    all_data.append(ds[variable_name])
                except Exception as e:
                    print(f"Error processing: {input_path}")
                    print(f"Error message: {e}")
            else:
                print(f"File not found: {input_path}")

    if len(all_data) == 0:
        raise ValueError("No wet-bulb temperature files were found.")

    all_data = xr.concat(all_data, dim=time_dim)
    latitude = all_data.latitude
    longitude = all_data.longitude

    for p in percentiles:
        threshold = all_data.quantile(p / 100, dim=time_dim, skipna=True)
        threshold_array = xr.DataArray(threshold.values, dims=["latitude", "longitude"], coords={"latitude": latitude, "longitude": longitude}, name=variable_name)
        output_path = os.path.join(output_root, f"MaxWetbulb_q{p}.nc")
        threshold_array.to_netcdf(output_path)
        print(f"Saved: {output_path}")


def calculate_precipitation_threshold(input_root, output_root, start_year=1981, end_year=2010, variable_name="tp", time_dim="valid_time", file_prefix="Total_Precipitation", percentiles=(90,)):
    """Calculate grid-cell-specific wet-day precipitation thresholds."""

    os.makedirs(output_root, exist_ok=True)
    all_data = []

    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            month_str = f"{month:02d}"
            input_file = f"{file_prefix}_{year}_{month_str}.nc"
            input_path = os.path.join(input_root, input_file)

            if os.path.exists(input_path):
                try:
                    ds = xr.open_dataset(input_path)
                    all_data.append(ds[variable_name])
                except Exception as e:
                    print(f"Error processing: {input_path}")
                    print(f"Error message: {e}")
            else:
                print(f"File not found: {input_path}")

    if len(all_data) == 0:
        raise ValueError("No precipitation files were found.")

    all_data = xr.concat(all_data, dim=time_dim)
    data_positive = all_data.where(all_data > 0)

    latitude = all_data.latitude
    longitude = all_data.longitude

    for p in percentiles:
        threshold = data_positive.quantile(p / 100, dim=time_dim, skipna=True)
        threshold_array = xr.DataArray(threshold.values, dims=["latitude", "longitude"], coords={"latitude": latitude, "longitude": longitude}, name=variable_name)
        output_path = os.path.join(output_root, f"Precip_q{p}.nc")
        threshold_array.to_netcdf(output_path)
        print(f"Saved: {output_path}")
