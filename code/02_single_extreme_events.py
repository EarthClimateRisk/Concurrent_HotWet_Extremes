# -*- coding: utf-8 -*-
"""
Single extreme-event identification.

Identifies extreme precipitation events and heatwaves from monthly NetCDF
files containing daily precipitation and daily maximum wet-bulb temperature.

Input:
    Precipitation: monthly NetCDF files containing daily precipitation data.
    Wet-bulb temperature: monthly NetCDF files containing daily maximum
    wet-bulb temperature data.
    Thresholds: grid-cell-specific thresholds calculated from the historical
    baseline period (1981–2010), including the P90 wet-day precipitation
    threshold and P95 daily maximum wet-bulb temperature threshold.
    Mask: 1° land mask defining the study area.

Output:
    Event: daily occurrence of extreme precipitation events or heatwaves.
    Duration: event duration recorded on the first day of each event.

Future periods:
    SSP2-4.5: 2015–2044 and 2036–2065.
    SSP5-8.5: 2014–2043 and 2028–2057.
"""

import os
import xarray as xr
import numpy as np
import pandas as pd
import rasterio
import logging


def process_highPrecip(events_array):
    """Identify extreme precipitation events and their durations."""

    events_array1 = np.copy(events_array)
    events_array_for_endss = np.copy(events_array)

    for i in range(events_array_for_endss.shape[0] - 2, -1, -1):
        events_array_for_endss[i, events_array_for_endss[i, ...] > 0] = events_array_for_endss[i + 1, events_array_for_endss[i, ...] > 0] + 1

    diff = np.zeros(events_array_for_endss.shape, dtype=np.int32)
    diff[0, ...] = events_array_for_endss[0, ...]
    diff[1:, ...] = np.diff(events_array_for_endss, axis=0)

    endss = np.zeros(events_array_for_endss.shape, dtype=int)
    endss[diff >= 1] = events_array_for_endss[diff >= 1]

    return events_array1, endss


def process_heatwaves(events_array):
    """Identify heatwave events lasting at least three consecutive days."""

    for i in range(events_array.shape[0] - 2, -1, -1):
        events_array[i, events_array[i, ...] > 0] = events_array[i + 1, events_array[i, ...] > 0] + 1

    diff = np.zeros(events_array.shape)
    diff[0, ...] = events_array[0, ...]
    diff[1:, ...] = np.diff(events_array, axis=0)

    events_array[diff == 2] = 0
    events_array[np.roll(diff == 2, 1, axis=0)] = 0
    events_array[diff == 1] = 0
    events_array[events_array > 0] = 1

    events_array_actual = np.copy(events_array)
    events_array_idea = np.copy(events_array_actual)

    for i in range(events_array_idea.shape[0] - 2, -1, -1):
        events_array_idea[i, events_array_idea[i, ...] > 0] = events_array_idea[i + 1, events_array_idea[i, ...] > 0] + 1

    diff = np.zeros(events_array_idea.shape)
    diff[0, ...] = events_array_idea[0, ...]
    diff[1:, ...] = np.diff(events_array_idea, axis=0)

    duration_idea = np.zeros(events_array_idea.shape, dtype=int)
    duration_idea[diff > 2] = events_array_idea[diff > 2]
    duration_idea[duration_idea < 3] = 0

    return events_array_actual, duration_idea


def identify_precipitation_events(input_root, threshold_root, output_root_event, output_root_duration, mask_path, scenarios):
    """Identify extreme precipitation events for historical and future periods."""

    with rasterio.open(mask_path) as src:
        land_mask = src.read(1)

    valid_mask = land_mask == 1

    threshold_path = os.path.join(threshold_root, "historical_Precip_q90.nc")
    ds_thr = xr.open_dataset(threshold_path)
    threshold = ds_thr["pr"].values
    lat = ds_thr["lat"].values
    lon = ds_thr["lon"].values

    for scenario in scenarios:
        if scenario == "historical":
            scenario_periods = {"1981-2010": list(range(1981, 2011))}
        elif scenario == "ssp245":
            scenario_periods = {"2015-2044": list(range(2015, 2045)), "2036-2065": list(range(2036, 2066))}
        elif scenario == "ssp585":
            scenario_periods = {"2014-2043": list(range(2014, 2044)), "2028-2057": list(range(2028, 2058))}

        for period_name, years in scenario_periods.items():
            print(f"Processing precipitation: {scenario} - {period_name}")

            pr_list = []
            time_list = []

            for year in years:
                for month in range(1, 13):
                    folder = scenario if year >= 2015 else "historical"
                    file = os.path.join(input_root, folder, "precipitation", f"{folder}_pr_{year}_{month:02d}.nc")

                    if os.path.exists(file):
                        ds = xr.open_dataset(file)
                        pr_list.append(ds["pr"].values)
                        time_list.append(pd.date_range(f"{year}-{month:02d}-01", periods=ds.sizes["time"], freq="D"))
                    else:
                        print(f"File not found: {file}")

            if not pr_list:
                logging.info(f"Skip {scenario}-{period_name}: no precipitation data")
                continue

            pr_all = np.concatenate(pr_list, axis=0)
            valid_mask_3d = np.broadcast_to(valid_mask, pr_all.shape)
            pr_all[~valid_mask_3d] = np.nan
            time_all = pd.DatetimeIndex(np.concatenate([t.values for t in time_list]))

            event_actual, duration = process_highPrecip((pr_all > threshold).astype(int))

            event_actual = event_actual.astype(np.float32)
            duration = duration.astype(np.float32)
            event_actual[~valid_mask_3d] = np.nan
            duration[~valid_mask_3d] = np.nan

            for year in years:
                year_mask = time_all.year == year

                if not np.any(year_mask):
                    continue

                coords = {"time": time_all[year_mask], "latitude": lat, "longitude": lon}

                event_path = os.path.join(output_root_event, scenario, period_name)
                duration_path = os.path.join(output_root_duration, scenario, period_name)
                os.makedirs(event_path, exist_ok=True)
                os.makedirs(duration_path, exist_ok=True)

                event_da = xr.DataArray(event_actual[year_mask, :, :], coords=coords, dims=["time", "latitude", "longitude"], name="event")
                duration_da = xr.DataArray(duration[year_mask, :, :], coords=coords, dims=["time", "latitude", "longitude"], name="duration")

                event_da.to_netcdf(os.path.join(event_path, f"{scenario}_event_{year}.nc"))
                duration_da.to_netcdf(os.path.join(duration_path, f"{scenario}_duration_{year}.nc"))


def identify_heatwave_events(input_root, threshold_root, output_root_event, output_root_duration, mask_path, scenarios):
    """Identify heatwaves for historical and future periods."""

    with rasterio.open(mask_path) as src:
        land_mask = src.read(1)

    valid_mask = land_mask == 1

    threshold_path = os.path.join(threshold_root, "historical_MaxWetbulb_q95.nc")
    ds_thr = xr.open_dataset(threshold_path)
    threshold = ds_thr["wetbulb"].values
    lat = ds_thr["latitude"].values
    lon = ds_thr["longitude"].values

    for scenario in scenarios:
        if scenario == "historical":
            scenario_periods = {"1981-2010": list(range(1981, 2011))}
        elif scenario == "ssp245":
            scenario_periods = {"2015-2044": list(range(2015, 2045)), "2036-2065": list(range(2036, 2066))}
        elif scenario == "ssp585":
            scenario_periods = {"2014-2043": list(range(2014, 2044)), "2028-2057": list(range(2028, 2058))}

        for period_name, years in scenario_periods.items():
            print(f"Processing heatwave: {scenario} - {period_name}")

            wetbulb_list = []
            time_list = []

            for year in years:
                for month in range(1, 13):
                    folder = scenario if year >= 2015 else "historical"
                    file = os.path.join(input_root, folder, "wetbulb", f"wetbulb_{folder}_{year}_{month:02d}.nc")

                    if os.path.exists(file):
                        ds = xr.open_dataset(file)
                        wetbulb_list.append(ds["wetbulb"].values)
                        time_list.append(pd.date_range(f"{year}-{month:02d}-01", periods=ds.sizes["time"], freq="D"))
                    else:
                        print(f"File not found: {file}")

            if not wetbulb_list:
                print(f"Skip {scenario}-{period_name}: no wet-bulb temperature data")
                continue

            wetbulb_all = np.concatenate(wetbulb_list, axis=0)
            valid_mask_3d = np.broadcast_to(valid_mask, wetbulb_all.shape)
            wetbulb_all[~valid_mask_3d] = np.nan
            time_all = pd.DatetimeIndex(np.concatenate([t.values for t in time_list]))

            event_actual, duration = process_heatwaves((wetbulb_all > threshold).astype(int))

            event_actual = event_actual.astype(np.float32)
            duration = duration.astype(np.float32)
            event_actual[~valid_mask_3d] = np.nan
            duration[~valid_mask_3d] = np.nan

            for year in years:
                year_mask = time_all.year == year

                if not np.any(year_mask):
                    continue

                coords = {"time": time_all[year_mask], "latitude": lat, "longitude": lon}

                event_path = os.path.join(output_root_event, scenario, period_name)
                duration_path = os.path.join(output_root_duration, scenario, period_name)
                os.makedirs(event_path, exist_ok=True)
                os.makedirs(duration_path, exist_ok=True)

                event_da = xr.DataArray(event_actual[year_mask, :, :], coords=coords, dims=["time", "latitude", "longitude"], name="event")
                duration_da = xr.DataArray(duration[year_mask, :, :], coords=coords, dims=["time", "latitude", "longitude"], name="duration")

                event_da.to_netcdf(os.path.join(event_path, f"{scenario}_event_{year}.nc"))
                duration_da.to_netcdf(os.path.join(duration_path, f"{scenario}_duration_{year}.nc"))
