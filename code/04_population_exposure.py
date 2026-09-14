# -*- coding: utf-8 -*-
"""
Population exposure calculation.

Calculates total and age-specific population exposure to SEPH, SHEP, and CHEP
as event duration multiplied by population, with exposure expressed in
person-days.

Input:
    Event duration:
        GeoTIFF files containing annual or multi-year mean annual compound-event
        duration (days).

    Total population:
        GeoTIFF files containing population counts.

    Age-specific population:
        GeoTIFF files containing population counts for ages 0–14, 15–64,
        and >=65 years.

    Historical annual exposure uses annual event duration and population from
    the corresponding year. Future warming period exposure uses multi-year mean
    annual event duration and population from the representative year of the
    corresponding warming period.

Output:
    GeoTIFF files containing total or age-specific population exposure
    in person-days.
"""

import os
import rasterio
import numpy as np


def calculate_population_exposure(duration_tif, population_tif, output_tif, out_nodata=-9999.0):
    """Calculate population exposure from event duration and population."""

    with rasterio.open(duration_tif) as src_dur, rasterio.open(population_tif) as src_pop:
        if src_dur.width != src_pop.width or src_dur.height != src_pop.height or src_dur.transform != src_pop.transform or src_dur.crs != src_pop.crs:
            raise ValueError(f"Raster spatial information does not match:\nDuration: {duration_tif}\nPopulation: {population_tif}")

        duration = src_dur.read(1).astype("float32")
        population = src_pop.read(1).astype("float32")

        valid = np.ones(duration.shape, dtype=bool)

        if src_dur.nodata is not None:
            valid &= duration != src_dur.nodata
        if src_pop.nodata is not None:
            valid &= population != src_pop.nodata

        valid &= duration != -9999
        valid &= population != -9999
        valid &= np.isfinite(duration)
        valid &= np.isfinite(population)

        exposure = np.full(duration.shape, out_nodata, dtype="float32")
        exposure[valid] = duration[valid] * population[valid]

        profile = src_dur.profile.copy()
        profile.update(dtype="float32", nodata=out_nodata, compress="lzw", count=1)

        os.makedirs(os.path.dirname(output_tif), exist_ok=True)

        with rasterio.open(output_tif, "w", **profile) as dst:
            dst.write(exposure, 1)


def calculate_historical_exposure(duration_dir, population_dir, output_dir, event_type, years=range(1981, 2021)):
    """Calculate annual population exposure for the historical period."""

    os.makedirs(output_dir, exist_ok=True)

    for year in years:
        duration_tif = os.path.join(duration_dir, f"{event_type}_duration_{year}.tif")
        population_tif = os.path.join(population_dir, f"population_{year}.tif")
        output_tif = os.path.join(output_dir, f"{event_type}_exposure_{year}.tif")

        if not os.path.exists(duration_tif) or not os.path.exists(population_tif):
            print(f"Skip: {year}")
            continue

        calculate_population_exposure(duration_tif, population_tif, output_tif)


def calculate_future_exposure(duration_tif, population_tif, output_tif):
    """Calculate total population exposure for a future warming period."""

    calculate_population_exposure(duration_tif, population_tif, output_tif)


def calculate_age_specific_exposure(duration_tif, age_population_tif, output_tif):
    """Calculate age-specific population exposure."""

    calculate_population_exposure(duration_tif, age_population_tif, output_tif)
