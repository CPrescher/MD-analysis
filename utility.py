# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'

import os
from copy import copy

import matplotlib.pyplot as plt

from ScatteringFactors import ScatteringFactors


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


def convert_density_to_atoms_per_cubic_angstrom(elemental_abundances, density):
    """
    Converts densities in g/cm3 into atoms per A^3
    :param elemental_abundances: dictionary with elements as key and abundances as relative numbers
    :param density: density in g/cm^3
    :return: density in atoms/A^3
    """

    # get_smallest abundance
    norm_elemental_abundances = normalize_elemental_abundances(elemental_abundances)
    mean_z = 0.0
    for key, val in norm_elemental_abundances.iteritems():
        mean_z += val * ScatteringFactors.atomic_weights['AW'][key]
    return density / mean_z * .602214129


def convert_atoms_per_cubic_angstrom_to_density(elemental_abundances, atomic_density):
    """
    Converts densities in g/cm3 into atoms per A^3
    :param elemental_abundances: dictionary with elements as key and abundances as relative numbers
    :param density: density in g/cm^3
    :return: density in atoms/A^3
    """

    # get_smallest abundance
    norm_elemental_abundances = normalize_elemental_abundances(elemental_abundances)
    mean_z = 0.0
    for key, val in norm_elemental_abundances.iteritems():
        mean_z += val * ScatteringFactors.atomic_weights['AW'][key]
    return atomic_density * mean_z / .602214129


def create_folder(parent_folder, folder_name):
    folder_path = os.path.join(parent_folder, folder_name)
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    return folder_path


def create_df_plot(df, filename, dpi=100, columns=None):
    plt.figure()
    # plot individual lines
    _, data_columns = df.shape
    if columns is None:
        for col in df.columns:
            plt.plot(df.index.values, df[col], label=col)
    else:
        plt.plot(df.index.values, df[columns], label=columns)

    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(filename, dpi=dpi)
    plt.close()