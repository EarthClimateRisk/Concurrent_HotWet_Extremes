# Concurrent_HotWet_Extremes

## Overview

This repository contains Python scripts for identifying three types of compound hot–wet extremes and assessing differences in their future population exposure:

* **SEPH**: extreme precipitation followed by a heatwave
* **SHEP**: heatwave followed by extreme precipitation
* **CHEP**: concurrent heatwave and extreme precipitation

The workflow includes extreme-event threshold calculation, single and compound event identification, population exposure calculation, exposure attribution, and country-level integrated risk assessment.

## Requirements

Tested with Python 3.10.

Main Python packages:

```text
numpy
pandas
xarray
rasterio
geopandas
openpyxl
netCDF4
```

Install dependencies using:

```bash
pip install numpy pandas xarray rasterio geopandas openpyxl netCDF4
```

## Input Data

### 1. ERA5-Land climate data

The observational climate data cover 1981–2020 and are provided on a common 1° × 1° grid.

Inputs include:

* Daily maximum wet-bulb temperature
* Daily precipitation

The daily data are stored in monthly NetCDF files.

### 2. CMIP6/ISIMIP3b climate data

Future climate analyses use the multimodel ensemble mean derived from seven bias-adjusted CMIP6 models provided by ISIMIP3b:

* CanESM5
* EC-Earth3
* GFDL-ESM4
* IPSL-CM6A-LR
* MIROC6
* MPI-ESM1-2-HR
* MRI-ESM2-0

Analysis-ready daily maximum wet-bulb temperature and precipitation data are provided for:

* historical
* SSP2-4.5
* SSP5-8.5

All climate inputs are provided on a 1° × 1° grid.

The original bias-adjusted CMIP6 data are publicly available from the ISIMIP data repository:[ISIMIP data repository](https://data.isimip.org/)

### 3. Population data

Historical population data are annual gridded population counts for 1981–2020.

Future warming-level analyses use population data from representative years under SSP2 and SSP5.

### 4. Population age structure

Age-structure data include three age groups:

* 0–14 years
* 15–64 years
* ≥65 years

Age-specific population data are used to calculate age-specific population exposure.

### 5. Human Development Index

Country-level HDI data are used in the integrated risk analysis.

Social vulnerability is represented as:

```text
Social vulnerability = 1 - HDI
```

### 6. Spatial boundaries

The regional boundary dataset is used for regional aggregation.

### 7. Analysis mask

The mask excludes Greenland, Antarctica, and desert regions with annual precipitation below 100 mm.

## Script Description

### 1. `01_extreme_thresholds.py`

Calculates grid-cell-specific thresholds for heatwaves and extreme precipitation.

**Input**

* Monthly NetCDF files containing daily maximum wet-bulb temperature
* Monthly NetCDF files containing daily precipitation

**Output**

* Wet-bulb temperature threshold NetCDF files
* Precipitation threshold NetCDF files

**Key parameters**

```text
Heatwave threshold: P95 daily maximum wet-bulb temperature
Extreme precipitation threshold: P90 wet-day precipitation
Historical baseline: 1981–2010
```

### 2. `02_single_extreme_events.py`

Identifies individual heatwave and extreme-precipitation events.

**Input**

* Daily maximum wet-bulb temperature
* Daily precipitation
* Grid-cell-specific P95 and P90 thresholds
* 1° analysis mask

**Output**

* Daily event indicators
* Event-duration data, with duration recorded on the first day of each event

**Key parameter**

```text
Minimum heatwave duration: 3 consecutive days
```

### 3. `03_compound_hotwet_events.py`

Identifies SEPH, SHEP, and CHEP from the single-extreme event series.

**Input**

* Daily heatwave occurrence (0/1)
* Daily extreme-precipitation occurrence (0/1)
* Heatwave duration recorded on event start days
* Extreme-precipitation duration recorded on event start days

**Output**

* Compound-event frequency
* Compound-event duration
* Compound-event indicator arrays

**Key parameters**

```text
SEPH/SHEP maximum interval: 7 days
CHEP minimum overlap: 1 day
```

### 4. `04_population_exposure.py`

Calculates total and age-specific population exposure.

**Input**

* Compound-event duration GeoTIFF
* Total or age-specific population GeoTIFF

**Output**

* Population exposure GeoTIFF in person-days

Population exposure is calculated as:

```text
Population exposure = event duration × population
```

Historical annual exposure uses population from the corresponding year. Future warming-period exposure uses mean annual duration and population from the representative year.

### 5. `05_exposure_attribution.py`

Quantifies the drivers of changes in total and older-adult population exposure.

**Input**

For total population exposure:

* Global mean compound-event duration
* Global population

For older-adult exposure:

* Compound-event duration
* Total population
* Proportion of population aged ≥65 years

**Output**

Total population exposure contributions:

* Climate
* Population
* Climate × Population

Older-adult population exposure contributions:

* Climate
* Population
* Ageing
* Climate × Population
* Climate × Ageing
* Population × Ageing
* Climate × Population × Ageing

### 6. `06_integrated_risk.py`

Calculates country-level integrated risk for SEPH, SHEP, and CHEP.

**Input**

* Normalized country-level compound-event duration
* Normalized population
* HDI

**Output**

* Country-level integrated risk

Integrated risk is calculated as:

```text
Risk = duration × population × (1 - HDI)
```

## Workflow

```text
1. Extreme-event threshold calculation
        ↓
2. Heatwave and extreme-precipitation identification
        ↓
3. SEPH / SHEP / CHEP identification
        ↓
4. Population exposure calculation
        ↓
5. Exposure attribution
        ↓
6. Integrated risk calculation
```

## Output Data

### Thresholds

Grid-cell-specific P95 wet-bulb temperature and P90 wet-day precipitation thresholds.

### Compound hot–wet events

Annual GeoTIFF files containing:

* annual event frequency
* annual cumulative event duration (days)

for SEPH, SHEP, and CHEP from ERA5-Land and CMIP6/ISIMIP3b analyses.

### Population exposure

GeoTIFF files containing:

* total population exposure
* older-adult population exposure

Population exposure is expressed in person-days.

### Attribution

Tables containing percentage contributions of climate, population, ageing, and their interactions to changes in population exposure.

### Integrated risk

Country-level integrated risk results for the historical baseline and future warming periods.

## Data Availability

The analysis-ready climate input data are available on Figshare.

Climate data: https://doi.org/10.6084/m9.figshare.33734311
