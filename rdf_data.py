# read in data

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ScatteringFactors.ScatteringFactors import calculate_coherent_scattering_factor
from utility import normalize_elemental_abundances, clean_element_name


def read_rdf_data(filename):
    file = open(filename, "r")

    file_str = []
    for line in file:
        file_str.append(line)

    line_index = 2
    data = {}
    data['r'] = []
    pair_columns = []
    while (line_index < len(file_str)):
        line_str = file_str[line_index]
        if not line_str.startswith("  "):
            current_pair = "{} - {}".format(*line_str.split())
            pair_columns.append(current_pair)
            data[current_pair] = []
        else:
            line_str_split = line_str.split()
            if len(data['r']) == len(data[current_pair]):
                data['r'].append(float(line_str_split[0]))
            data[current_pair].append(float(line_str_split[1]))
        line_index += 1

    df = pd.DataFrame(data)
    df = df[["r"] + pair_columns]
    df = df.set_index("r")
    return df


def calculate_atomic_sq(df, atomic_density, max_q=20):
    q = np.arange(0.01, max_q, 0.01)
    r = df.index.values

    sq_df = pd.DataFrame()
    sq_df['q'] = q

    plt.figure()
    for col in df.columns:
        g_r = df[col]
        int_val = []
        delta_r = 0.5 * np.pi / np.max(r)
        modification = np.sin(delta_r) / delta_r
        for val in q:
            int_val.append(np.trapz(r * np.sin(val * r) * (g_r - 1) * modification, r))
        s_q = 1 + 4 * np.pi * atomic_density / q * np.array(int_val)
        sq_df[col] = s_q

    sq_df = sq_df.set_index('q')
    return sq_df


def calculate_combined_sq(s_q_df, concentrations):
    concentrations = normalize_elemental_abundances(concentrations)
    q = s_q_df.index.values

    # calculate denominator for weighting
    denom = np.zeros(q.shape)
    for element, c in concentrations.iteritems():
        f = calculate_coherent_scattering_factor(element, q)
        denom += f * c
    denom = denom ** 2

    sq_weighted_df = pd.DataFrame(index=q)
    sq_weighted_df.index.name = "q"
    s_q_sum = np.zeros(q.shape)

    for col in s_q_df.columns:
        # getting element names:
        elements_str = col.split(" - ")
        element1 = clean_element_name(elements_str[0])
        element2 = clean_element_name(elements_str[1])

        if element1 == element2:
            factor = 1.0
        else:
            factor = 2.0

        f_element1 = calculate_coherent_scattering_factor(element1, q)
        f_element2 = calculate_coherent_scattering_factor(element2, q)

        weight = factor * (concentrations[element1] * f_element1 * concentrations[element2] * f_element2) / denom
        s_q = s_q_df[col] * weight
        s_q_sum += s_q
        sq_weighted_df[col] = s_q

    sq_weighted_df["sum"] = s_q_sum
    return sq_weighted_df


def calculate_weighted_fr(sq_weighted_df, use_modification=False, q_max=None, r_limits=(0, 10)):
    r = np.arange(r_limits[0], r_limits[1], 0.01)

    q = sq_weighted_df.index.values

    if q_max is None:
        q_max = np.max(q)

    index = q <= q_max
    q = q[index]

    fr_weighted_df = pd.DataFrame(index=r)
    fr_weighted_df.index.name = "r"

    if use_modification:
        modification = np.sin(q * np.pi / np.max(q)) / (q * np.pi / np.max(q))
    else:
        modification = 1

    q = np.array(q)

    for col in sq_weighted_df.columns:
        sq = np.array(sq_weighted_df[col])[index]
        fr = 2.0 / np.pi * np.trapz(modification * q * (sq - 1) *
                                    np.array(np.sin(np.mat(q).T * np.mat(r))).T, q)
        fr_weighted_df[col] = fr

    return fr_weighted_df


def calculate_weighted_rdf(fr_weighted_df, atomic_density):
    r = fr_weighted_df.index.values
    rdf_weighted_df = pd.DataFrame(index=r)
    rdf_weighted_df.index.name = "r"
    for col in fr_weighted_df.columns:
        rdf_weighted_df[col] = 1 + fr_weighted_df[col] / (4.0 * np.pi * r * atomic_density)

    return rdf_weighted_df


if __name__ == "__main__":
    import os

    files = os.listdir('.')
    # for file in files:
    # if file.startswith("RDFDAT") and not file.endswith(".png"):
    rdf_df = read_rdf_data("RDFDAT-2")
    s_q_df = calculate_atomic_sq(rdf_df, 0.07)
    s_q_df.to_csv("RDFDAT-1-s_q.csv")
    s_q_df = pd.read_csv("RDFDAT-1-s_q.csv", index_col=0)
    concentrations = {
        "Mg": 2,
        "Si": 1,
        "O": 4
    }

    calculate_combined_sq(s_q_df, concentrations)

    plt.show()
