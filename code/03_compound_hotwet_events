# -*- coding: utf-8 -*-
"""
Compound hot-wet event identification.

Identifies three types of compound hot-wet events from daily extreme-event
and event-duration data stored in CSV format. The CSV files should be read
into arrays before calling the functions in this module.

Input:
    Heatwave event data: daily binary data, where 1 indicates a heatwave day
    and 0 indicates no heatwave.

    Extreme-precipitation event data: daily binary data, where 1 indicates an
    extreme-precipitation day and 0 indicates no extreme precipitation.

    Heatwave duration data: the duration of each heatwave is recorded only on
    the first day of the event, with all other days set to 0.

    Extreme-precipitation duration data: the duration of each precipitation
    event is recorded only on the first day of the event, with all other days
    set to 0.

Compound-event definitions:
    SEPH: extreme precipitation followed by a heatwave within 7 days after
    the precipitation event ends.

    SHEP: heatwave followed by extreme precipitation within 7 days after
    the heatwave ends.

    CHEP: heatwave and extreme precipitation overlap for at least 1 day.
"""

import numpy as np


def calculate_compound_events(heat_matrix, flood_matrix, heat_duration_mat, flood_duration_mat, event_type):
    """Calculate SEPH, SHEP, or CHEP events for all grid cells."""

    event_count = np.zeros(heat_matrix.shape[0], dtype=np.float32)
    event_duration = np.zeros(heat_matrix.shape[0], dtype=np.float32)
    compound_dura_event = np.zeros_like(heat_matrix, dtype=np.int32)

    if event_type in ["SEPH", "SHEP"]:
        compound_heat_event = np.zeros_like(heat_matrix, dtype=np.int32)
        compound_precipitation_event = np.zeros_like(heat_matrix, dtype=np.int32)

    for i in range(heat_matrix.shape[0]):
        heat = heat_matrix[i, :]
        flood = flood_matrix[i, :]
        heat_duration = heat_duration_mat[i, :]
        flood_duration = flood_duration_mat[i, :]

        if event_type == "SEPH":
            event_count[i], event_duration[i], compound_heat_event[i, :], compound_precipitation_event[i, :], compound_dura_event[i, :] = _calculate_FloodToHeat_for_point(heat, flood, heat_duration, flood_duration)

        elif event_type == "SHEP":
            event_count[i], event_duration[i], compound_heat_event[i, :], compound_precipitation_event[i, :], compound_dura_event[i, :] = _calculate_HeatToFlood_for_point(heat, flood, heat_duration, flood_duration)

        elif event_type == "CHEP":
            event_count[i], event_duration[i], compound_dura_event[i, :] = _calculate_HeatAndFlood_for_point(heat, flood, heat_duration, flood_duration)

        else:
            raise ValueError("event_type must be 'SEPH', 'SHEP', or 'CHEP'.")

    if event_type in ["SEPH", "SHEP"]:
        return event_count, event_duration, compound_heat_event, compound_precipitation_event, compound_dura_event

    return event_count, event_duration, compound_dura_event


# ============================================================
# SEPH: extreme precipitation followed by heatwave
# ============================================================

def _calculate_FloodToHeat_for_point(heatwave_matrix, flood_matrix, heatwave_duration_matrix, flood_duration_matrix):
    heat_duration = heatwave_duration_matrix.copy()
    flood_duration = flood_duration_matrix.copy()

    compound_event_tracker = np.zeros_like(heatwave_matrix, dtype=np.int32)
    compound_heat = np.zeros_like(heatwave_matrix, dtype=np.int32)
    compound_precipitation = np.zeros_like(heatwave_matrix, dtype=np.int32)

    A1_count = 0
    A2_count = 0
    internal_day = 0

    while internal_day < flood_matrix.shape[0]:
        if flood_duration[internal_day] >= 1:
            flood_duration_number = flood_duration[internal_day]
            flood_end_day = internal_day + flood_duration_number
            flood_end_day = min(flood_end_day, flood_matrix.shape[0] - 1)
            heat_duration_window = heat_duration[flood_end_day:flood_end_day + 7]

            if np.any(heat_duration_window > 1):
                A1_count = flood_duration[internal_day]
                compound_precipitation[internal_day:flood_end_day] = 1

                heat_indices = np.where(heat_duration_window > 1)[0]
                heat_start = heat_indices[0]
                heat_start_global = flood_end_day + heat_start
                heat_end_global = heat_start_global + heat_duration[heat_start_global]

                A2_count = heat_duration[heat_start_global]
                compound_heat[heat_start_global:heat_end_global] = 1

                compound_event_tracker[internal_day] += A1_count + A2_count
                internal_day = heat_end_global
            else:
                internal_day += 1
        else:
            internal_day += 1

    sum_across_time = np.sum(compound_event_tracker, axis=0)
    event_count = np.count_nonzero(compound_event_tracker > 0, axis=0)

    return event_count, sum_across_time, compound_heat, compound_precipitation, compound_event_tracker


# ============================================================
# SHEP: heatwave followed by extreme precipitation
# ============================================================

def _calculate_HeatToFlood_for_point(heatwave_matrix, flood_matrix, heatwave_duration_matrix, flood_duration_matrix):
    heat_duration = heatwave_duration_matrix.copy()
    flood_duration = flood_duration_matrix.copy()

    compound_event_tracker = np.zeros_like(heatwave_matrix, dtype=np.int32)
    compound_heat = np.zeros_like(heatwave_matrix, dtype=np.int32)
    compound_precipitation = np.zeros_like(heatwave_matrix, dtype=np.int32)

    A1_count = 0
    A2_count = 0
    internal_day = 0

    while internal_day < heatwave_matrix.shape[0]:
        if heat_duration[internal_day] > 1:
            heat_duration_number = heat_duration[internal_day]
            heatwave_end_day = internal_day + int(heat_duration_number)
            heatwave_end_day = min(heatwave_end_day, heatwave_matrix.shape[0] - 1)
            flood_duration_window = flood_duration[heatwave_end_day:heatwave_end_day + 7]

            if np.any(flood_duration_window >= 1):
                A1_count = heat_duration[internal_day]
                compound_heat[internal_day:heatwave_end_day] = 1

                flood_indices = np.where(flood_duration_window >= 1)[0]
                flood_start = flood_indices[0]
                flood_start_global = heatwave_end_day + flood_start
                flood_end_global = flood_start_global + flood_duration[flood_start_global]

                A2_count = flood_duration[flood_start_global]
                compound_precipitation[flood_start_global:flood_end_global] = 1

                compound_event_tracker[internal_day] += A1_count + A2_count
                internal_day = flood_end_global
            else:
                internal_day += 1
        else:
            internal_day += 1

    sum_across_time = np.sum(compound_event_tracker, axis=0)
    event_count = np.count_nonzero(compound_event_tracker > 0, axis=0)

    return event_count, sum_across_time, compound_heat, compound_precipitation, compound_event_tracker


# ============================================================
# CHEP: concurrent heatwave and extreme precipitation
# ============================================================

def _calculate_HeatAndFlood_for_point(heatwave_matrix, exflood_matrix, heatwave_duration_matrix, flood_duration_matrix):
    event_count = 0
    sum_across_time = 0

    exheat_duration = heatwave_duration_matrix.copy()
    exflood_duration = flood_duration_matrix.copy()
    event_matrix = exflood_matrix * 10 + heatwave_matrix

    compound_event_tracker = np.zeros_like(heatwave_matrix, dtype=np.int32)
    post_event_calculated_mask = np.zeros_like(heatwave_matrix)

    current_day = 0
    lookback_days = 10

    while current_day < heatwave_matrix.shape[0]:
        A1_count = A2_count = A3_count = 0
        last_flood_start = flood_start_global2 = flood_start_global3 = 0
        flood_end_global2 = flood_end_global3 = 0
        flood_start_global = flood_end_global = 0
        internal_day = current_day

        while internal_day < heatwave_matrix.shape[0]:
            if exheat_duration[internal_day] > 1:
                heat_duration_number = exheat_duration[internal_day]
                heatwave_end_day = min(internal_day + int(heat_duration_number) - 1, heatwave_matrix.shape[0] - 1)
                flood_duration_window = exflood_matrix[internal_day:heatwave_end_day + 1]

                if np.any(flood_duration_window >= 1):
                    A1_count = exheat_duration[internal_day]
                    post_event_calculated_mask[internal_day:heatwave_end_day + 1] = 1

                    flood_indices = np.where(flood_duration_window >= 1)[0]
                    flood_start = flood_indices[0]
                    flood_duration_end_index = flood_indices[-1]
                    flood_start_global = internal_day + flood_start
                    flood_end_global = internal_day + flood_duration_end_index

                    if np.any(exflood_matrix[flood_start_global:flood_end_global + 1] == 0):
                        for j in range(flood_end_global, max(flood_end_global - lookback_days, -1), -1):
                            if exflood_duration[j] > 0 and j + exflood_duration[j] - 1 >= flood_end_global:
                                last_flood_start = j
                                flood_end_global2 = j + exflood_duration[j] - 1
                                unprocessed = post_event_calculated_mask[last_flood_start:flood_end_global2 + 1] == 0
                                A2_count = np.sum((event_matrix[last_flood_start:flood_end_global2 + 1] > 0) & unprocessed)
                                post_event_calculated_mask[last_flood_start:flood_end_global2 + 1] = 1
                                break

                        if flood_start_global == internal_day:
                            search_start = max(0, flood_start_global - lookback_days)

                            for i in range(flood_start_global - 1, search_start - 1, -1):
                                if flood_duration_matrix[i] > 0 and i + flood_duration_matrix[i] - 1 >= flood_start_global:
                                    flood_start_global2 = i
                                    unprocessed = post_event_calculated_mask[flood_start_global2:flood_start_global + 1] == 0
                                    A3_count = np.sum((event_matrix[flood_start_global2:flood_start_global + 1] > 0) & unprocessed)
                                    post_event_calculated_mask[flood_start_global2:flood_start_global + 1] = 1
                                    break

                    else:
                        for j in range(flood_start_global, max(flood_start_global - lookback_days, -1), -1):
                            if exflood_duration[j] > 0 and j + exflood_duration[j] - 1 >= flood_end_global:
                                flood_start_global3 = j
                                flood_end_global3 = j + exflood_duration[j] - 1
                                unprocessed = post_event_calculated_mask[flood_start_global3:flood_end_global3 + 1] == 0
                                A2_count = np.sum((event_matrix[flood_start_global3:flood_end_global3 + 1] > 0) & unprocessed)
                                post_event_calculated_mask[flood_start_global3:flood_end_global3 + 1] = 1
                                break

                    compound_event_tracker[internal_day] += A1_count + A2_count + A3_count
                    latest_flood_end = max(flood_end_global2, flood_end_global3, flood_end_global)
                    current_day = latest_flood_end + 1
                    break
                else:
                    internal_day += 1
            else:
                internal_day += 1
        else:
            break

        sum_across_time = np.sum(compound_event_tracker, axis=0)
        event_count = np.count_nonzero(compound_event_tracker > 0, axis=0)

    return event_count, sum_across_time, compound_event_tracker.reshape(-1).astype(np.int32)
