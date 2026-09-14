# -*- coding: utf-8 -*-
"""
Population-exposure attribution.

Calculates the global attribution of changes in total population exposure and
older-adult population exposure to SEPH, SHEP, and CHEP.

Input:
    Total population exposure attribution:
        Excel file containing global mean compound-event duration for the
        historical baseline and future warming periods.
        Excel file containing global population for the historical baseline
        and representative future years.

    Older-adult exposure attribution:
        Shapefile containing grid-level compound-event duration, total
        population, and >=65 year population proportion for the historical
        and future periods.

Output:
    Total population exposure:
        Percentage contributions from climate change, population change,
        and their interaction.

    Older-adult population exposure:
        Percentage contributions from climate change, population change,
        ageing, climate-population interaction, climate-ageing interaction,
        population-ageing interaction, and the three-way interaction.
"""

import numpy as np
import pandas as pd
import geopandas as gpd


def calculate_total_exposure_attribution(duration_excel, population_excel):
    """Calculate global attribution of total population exposure."""

    duration_df = pd.read_excel(duration_excel)
    population_df = pd.read_excel(population_excel)

    duration_df.columns = duration_df.columns.str.replace("_mean", "", regex=False)

    duration_row = duration_df.iloc[0]
    population_row = population_df.iloc[0]

    group_cols = [duration_df.columns[i:i + 5] for i in range(1, duration_df.shape[1], 5)]
    population_columns = ["pop2_2030", "pop2_2050", "pop5_2030", "pop5_2040"]
    population_base = "pop2000"

    results = []

    for group_index, col_group in enumerate(group_cols):
        base_col = col_group[0]

        for i, future_col in enumerate(col_group[1:]):
            P0 = population_row[population_base]
            Pt = population_row[population_columns[i]]
            C0 = duration_row[base_col]
            Ct = duration_row[future_col]

            delta_P = Pt - P0
            delta_C = Ct - C0

            contribution_population = delta_P * C0
            contribution_climate = delta_C * P0
            contribution_interaction = delta_P * delta_C
            total = contribution_population + contribution_climate + contribution_interaction

            if total != 0:
                population_pct = contribution_population / total * 100
                climate_pct = contribution_climate / total * 100
                interaction_pct = contribution_interaction / total * 100
            else:
                population_pct = climate_pct = interaction_pct = np.nan

            results.append({
                "Group": f"Group{group_index + 1}",
                "Target": future_col,
                "Climate_Contribution": climate_pct,
                "Population_Contribution": population_pct,
                "Climate_Population_Interaction": interaction_pct
            })

    return pd.DataFrame(results)


def calculate_older_adult_simulations(df, event_code, du_base_id, du_target_id, pop_base, pop_target, age_base, age_target):
    """Calculate the eight simulation states used for older-adult attribution."""

    df = df.copy()

    du_base = f"DU{event_code}{du_base_id}"
    du_target = f"DU{event_code}{du_target_id}"

    df["sim_control"] = df[du_base] * df[pop_base] * df[age_base]
    df["sim_climate"] = df[du_target] * df[pop_base] * df[age_base]
    df["sim_population"] = df[du_base] * df[pop_target] * df[age_base]
    df["sim_ageing"] = df[du_base] * df[pop_base] * df[age_target]
    df["sim_population_ageing"] = df[du_base] * df[pop_target] * df[age_target]
    df["sim_climate_ageing"] = df[du_target] * df[pop_base] * df[age_target]
    df["sim_climate_population"] = df[du_target] * df[pop_target] * df[age_base]
    df["sim_all"] = df[du_target] * df[pop_target] * df[age_target]

    return df


def calculate_older_adult_contributions(sim):
    """Calculate the seven contributions from globally summed simulations."""

    control = sim["sim_control"]
    climate = sim["sim_climate"]
    population = sim["sim_population"]
    ageing = sim["sim_ageing"]
    population_ageing = sim["sim_population_ageing"]
    climate_ageing = sim["sim_climate_ageing"]
    climate_population = sim["sim_climate_population"]
    all_change = sim["sim_all"]

    C_climate = climate - control
    C_population = population - control
    C_ageing = ageing - control
    C_population_ageing = population_ageing - population - ageing + control
    C_climate_ageing = climate_ageing - climate - ageing + control
    C_climate_population = climate_population - climate - population + control
    C_climate_population_ageing = all_change - population_ageing - climate_ageing - climate_population + climate + population + ageing - control

    total = C_climate + C_population + C_ageing + C_population_ageing + C_climate_ageing + C_climate_population + C_climate_population_ageing

    if np.isclose(total, 0):
        return {
            "Climate_Contribution": np.nan,
            "Population_Contribution": np.nan,
            "Ageing_Contribution": np.nan,
            "Population_Ageing_Interaction": np.nan,
            "Climate_Ageing_Interaction": np.nan,
            "Climate_Population_Interaction": np.nan,
            "Climate_Population_Ageing_Interaction": np.nan
        }

    return {
        "Climate_Contribution": C_climate / total * 100,
        "Population_Contribution": C_population / total * 100,
        "Ageing_Contribution": C_ageing / total * 100,
        "Population_Ageing_Interaction": C_population_ageing / total * 100,
        "Climate_Ageing_Interaction": C_climate_ageing / total * 100,
        "Climate_Population_Interaction": C_climate_population / total * 100,
        "Climate_Population_Ageing_Interaction": C_climate_population_ageing / total * 100
    }


def calculate_older_adult_exposure_attribution(age65_shapefile):
    """Calculate global seven-factor attribution of >=65-year population exposure."""

    df = gpd.read_file(age65_shapefile)

    events = {
        "CHEP": "haf",
        "SHEP": "htf",
        "SEPH": "fth"
    }

    scenarios = {
        "SSP245_1.5C_vs_HIS": {
            "du_base_id": 1, "du_target_id": 2,
            "pop_base": "his20", "pop_target": "SSP230",
            "age_base": "HIS00", "age_target": "SSP230_1"
        },
        "SSP245_2.0C_vs_HIS": {
            "du_base_id": 1, "du_target_id": 3,
            "pop_base": "his20", "pop_target": "SSP250",
            "age_base": "HIS00", "age_target": "SSP250_1"
        },
        "SSP245_2.0C_vs_1.5C": {
            "du_base_id": 2, "du_target_id": 3,
            "pop_base": "SSP230", "pop_target": "SSP250",
            "age_base": "SSP230_1", "age_target": "SSP250_1"
        },
        "SSP585_1.5C_vs_HIS": {
            "du_base_id": 1, "du_target_id": 4,
            "pop_base": "his20", "pop_target": "SSP530",
            "age_base": "HIS00", "age_target": "SSP530_1"
        },
        "SSP585_2.0C_vs_HIS": {
            "du_base_id": 1, "du_target_id": 5,
            "pop_base": "his20", "pop_target": "SSP540",
            "age_base": "HIS00", "age_target": "SSP540_1"
        },
        "SSP585_2.0C_vs_1.5C": {
            "du_base_id": 4, "du_target_id": 5,
            "pop_base": "SSP530", "pop_target": "SSP540",
            "age_base": "SSP530_1", "age_target": "SSP540_1"
        }
    }

    sim_cols = [
        "sim_control",
        "sim_climate",
        "sim_population",
        "sim_ageing",
        "sim_population_ageing",
        "sim_climate_ageing",
        "sim_climate_population",
        "sim_all"
    ]

    results = []

    for event_name, event_code in events.items():
        for scenario_name, scenario in scenarios.items():
            df_sim = calculate_older_adult_simulations(
                df,
                event_code,
                scenario["du_base_id"],
                scenario["du_target_id"],
                scenario["pop_base"],
                scenario["pop_target"],
                scenario["age_base"],
                scenario["age_target"]
            )

            global_sim = df_sim[sim_cols].sum()
            contributions = calculate_older_adult_contributions(global_sim)

            result = {
                "Event": event_name,
                "Scenario": scenario_name
            }

            result.update(contributions)
            results.append(result)

    return pd.DataFrame(results)
