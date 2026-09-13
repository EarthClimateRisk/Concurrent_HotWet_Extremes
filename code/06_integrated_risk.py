# -*- coding: utf-8 -*-
"""
Integrated risk calculation.

Calculates country-level integrated risk for SEPH, SHEP, and CHEP using
normalized compound-event duration, population, and HDI.

Input:
    Duration:
        Excel file containing normalized country-level compound-event duration
        for the historical baseline and future warming periods.

    Population:
        Excel file containing normalized country-level population for the
        historical baseline and representative future years.

    HDI:
        Excel file containing normalized country-level HDI for the historical
        baseline and future warming periods.

Output:
    Country-level integrated risk calculated as:

        risk = duration × population × (1 - HDI)

    Risk is calculated separately for SEPH, SHEP, and CHEP for the historical
    baseline and four future warming periods.
"""

import pandas as pd


def calculate_integrated_risk(duration_excel, population_excel, hdi_excel):
    """Calculate country-level integrated risk for SEPH, SHEP, and CHEP."""

    duration_df = pd.read_excel(duration_excel)
    pop_df = pd.read_excel(population_excel)
    hdi_df = pd.read_excel(hdi_excel)

    country = duration_df.iloc[:, 0]
    duration_data = duration_df.iloc[:, 1:]

    duration_groups = [
        duration_data.columns[0:5],
        duration_data.columns[5:10],
        duration_data.columns[10:15]
    ]

    pop_cols = ["pop2000", "pop2_2030", "pop2_2050", "pop5_2030", "pop5_2040"]
    hdi_cols = ["ssp245_hdi1", "ssp245_hdi2", "ssp245_hdi3", "ssp5_hdi2", "ssp5_hdi3"]

    results = {}

    for i, group in enumerate(duration_groups):
        risk_df = pd.DataFrame()
        risk_df["country"] = country

        for j in range(5):
            dur_col = group[j]
            pop_col = pop_cols[j]
            hdi_col = hdi_cols[j]
            risk_col = f"risk_{dur_col}"

            risk_df[risk_col] = duration_df[dur_col] * pop_df[pop_col] * (1 - hdi_df[hdi_col])

        results[f"group{i + 1}"] = risk_df

    return results
