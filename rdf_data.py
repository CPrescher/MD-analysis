
#read in data
import numpy as np
import pandas as pd
from copy import copy
import matplotlib.pyplot as plt

from ScatteringFactors import calculate_coherent_scattering_factor

dpi = 100


def create_rdf_plot(filename):
    file = open(filename, "r")

    file_str = []
    for line in file:
        file_str.append(line)

    line_index = 2
    data = {}
    data['r'] = []
    pair_columns = []
    while(line_index < len(file_str)):
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

    plt.figure()
    # plot individual lines
    _, data_columns = df.shape
    for col in xrange(1, data_columns):
        plt.plot(df.iloc[:, 0], df.iloc[:, col], label=df.columns[col])

    plt.legend()
    plt.tight_layout()
    plt.savefig(filename + ".png", dpi=dpi)

    return df


def calculate_atomic_sq(df, number_density):

    q = np.linspace(0, 20, 2000)
    r = df["r"]

    s_q_df = pd.DataFrame()
    s_q_df['q'] = q

    _, data_columns = df.shape

    plt.figure()
    for col in xrange(1, data_columns):
        g_r = df.iloc[:, col]
        int_val = []
        delta_r = 0.5 * np.pi / np.max(r)
        modification = np.sin(delta_r) / delta_r
        for val in q:
            int_val.append(np.trapz(r * np.sin(val * r) * (g_r - 1) * modification, r))
        s_q = 1 + 4 * np.pi * number_density / q * np.array(int_val)
        s_q_df[df.columns[col]] = s_q

        plt.plot(q, s_q, label=df.columns[col])
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("s_q.png", dpi=dpi)
    return s_q_df


def calculate_combined_sq(s_q_df, concentrations):

    concentrations = normalize_elemental_abundances(concentrations)
    q = s_q_df["q"]

    # calculate denominator for weighting
    denom = np.zeros(q.shape)
    for element, c in concentrations.iteritems():
        f = calculate_coherent_scattering_factor(element, q)
        denom += f * c
    denom = denom ** 2

    s_q_weighted_df = pd.DataFrame()
    s_q_weighted_df["q"] = q
    s_q_sum = np.zeros(q.shape)

    plt.figure()

    for col in s_q_df.columns:
        if col == "q":
            continue

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
        s_q_weighted_df[col] = s_q

        plt.plot(q, s_q, label=col)

    s_q_weighted_df["sum"] = s_q_sum
    plt.plot(q, s_q_sum, label="sum")
    plt.legend(loc=4)
    plt.savefig("S_q_X.png", dpi=dpi)
    return s_q_weighted_df


################################
#### HELPER FUNCTIONS ##########
################################


def clean_element_name(element_str):
    element = element_str[:2]
    if element[1] == "_":
        element = element[0]
    else:
        element = element[0] + element[1].lower()
    return element


def normalize_elemental_abundances(elemental_abundances):
    """
    normalizes elemental abundances to 1
    :param elemental_abundances: dictionary with elements as key and abundances as relative numbers
    :return: normalized elemental abundances dictionary dictionary
    """
    sum = 0.0
    for key, val in elemental_abundances.iteritems():
        sum += val

    result = copy(elemental_abundances)

    for key in result:
        result[key] /= sum

    return result

if __name__ == "__main__":
    import os
    files = os.listdir('.')
    # for file in files:
        # if file.startswith("RDFDAT") and not file.endswith(".png"):
    rdf_df = create_rdf_plot("RDFDAT-2")
    s_q_df = calculate_atomic_sq(rdf_df,  0.07)
    s_q_df.to_csv("RDFDAT-1-s_q.csv")
    s_q_df = pd.read_csv("RDFDAT-1-s_q.csv", index_col=0)
    concentrations = {
        "Mg": 2,
        "Si": 1,
        "O": 4
    }

    calculate_combined_sq(s_q_df, concentrations)

    plt.show()
