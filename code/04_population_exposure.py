# -*- coding: utf-8 -*-
"""
Population exposure calculation.

Calculates population exposure to SEPH, SHEP, and CHEP as event duration
multiplied by population, with exposure expressed in person-days.

Input:
    Event duration: GeoTIFF files containing compound-event duration (days).

    ERA5 population: annual GeoTIFF files containing population counts for
    each year from 1981 to 2020.

    CMIP6 population: GeoTIFF files containing population counts for the
    representative year of each warming period.

    For ERA5, annual event duration is combined with population in the
    corresponding year. For CMIP6, mean annual event duration over each
    30-year warming period is combined with population for the representative
    year of that period.

Output:
    GeoTIFF files of population exposure in person-days.
"""

import os
import rasterio
import numpy as np


def calculate_population_exposure(duration_path, population_path, output_path, invalid_below=None, output_nodata=np.nan):
    """Calculate population exposure from event duration and population."""

    with rasterio.open(duration_path) as duration_src, rasterio.open(population_path) as pop_src:
        duration_data = duration_src.read(1).astype(np.float32)
        population_data = pop_src.read(1).astype(np.float32)

        if duration_src.nodata is not None:
            duration_data = np.where(duration_data == duration_src.nodata, np.nan, duration_data)

        if pop_src.nodata is not None:
            population_data = np.where(population_data == pop_src.nodata, np.nan, population_data)

        if invalid_below is not None:
            duration_data[duration_data < invalid_below] = np.nan

        exposure = duration_data * population_data

        profile = duration_src.profile
        profile.update(dtype=rasterio.float32, nodata=output_nodata, compress="lzw")

        if not np.isnan(output_nodata):
            exposure[np.isnan(exposure)] = output_nodata

        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(exposure.astype(np.float32), 1)


def calculate_era5_exposure(duration_dir, population_dir, output_dir, event_type, years=range(1981, 2021)):
    """Calculate annual population exposure for ERA5-based compound events."""

    os.makedirs(output_dir, exist_ok=True)

    for year in years:
        duration_path = os.path.join(duration_dir, f"{event_type}_duration_{year}.tif")
        population_path = os.path.join(population_dir, f"population_{year}.tif")
        output_path = os.path.join(output_dir, f"exposure_{year}.tif")

        calculate_population_exposure(duration_path, population_path, output_path, output_nodata=np.nan)


def calculate_cmip6_exposure(base_dir, output_dir, default_population_path, population_map, models, scenarios, scenario_periods_dict, event_folders):
    """Calculate population exposure for CMIP6 warming periods."""

    for model in models:
        for scenario in scenarios:
            for period_name in scenario_periods_dict[scenario]:
                population_path = population_map.get((scenario, period_name), default_population_path)

                for event_type, folder_name in event_folders.items():
                    duration_path = os.path.join(base_dir, folder_name, model, scenario, period_name, f"{model}_{scenario}_{period_name}_{event_type}_yearmean.tif")

                    if not os.path.exists(duration_path):
                        print(f"Skip: {duration_path}")
                        continue

                    event_output_dir = os.path.join(output_dir, event_type, model, scenario, period_name)
                    output_path = os.path.join(event_output_dir, f"{model}_{scenario}_{period_name}_{event_type}_exposure.tif")

                    calculate_population_exposure(duration_path, population_path, output_path, invalid_below=-1, output_nodata=-9999.0)
