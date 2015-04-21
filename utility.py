# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'

import os
from copy import copy

import matplotlib.pyplot as plt
import numpy as np

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

    plt.legend(loc=4)
    plt.tight_layout()
    plt.savefig(filename, dpi=dpi)
    plt.close()


def create_df_surface_plot(df, filename,
                           x_label= "", y_label="", z_label="",
                           level_limits=(0,2), level_bins=100,
                           start_index=0, skip=0,
                           x_limits=None,  dpi=100, fig_size=(14, 5.8)):
    X = []
    Y = []
    Z = []

    x = df.index.values
    if x_limits is not None:
        data_ind = np.logical_and(x>=x_limits[0], x<=x_limits[1])
    else:
        data_ind = np.ones(x.shape, dtype=bool)
    x = x[data_ind]

    for ind, col in enumerate(df.columns):
        if ind<start_index:
            continue
        if skip is not 0:
            if not (ind+start_index) % skip == 0:
                continue

        X.append(x)
        Y.append(np.ones(len(x)) * float(col))
        Z.append(np.array(df[col])[data_ind])

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)

    plt.figure(figsize=fig_size)

    levels = plt.MaxNLocator(nbins=100).tick_values(level_limits[0], level_limits[1])
    plt.contourf(X, Y, Z, levels=levels)
    plt.colorbar(label=z_label, ticks = np.linspace(level_limits[0], level_limits[1], 5))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.savefig(filename, dpi=dpi)
    plt.close()

def create_df_stack_plot(df, filename, sep=0.25, start_index=0, skip=0, x_limits=None,
                         x_label="", y_label="",
                         text_pos=(0.9, 0.8), dpi=100, fig_size=(8,10)):
    plt.figure(figsize=fig_size)
    # plot individual lines
    x = np.array(df.index.values)
    if x_limits is not None:
        data_ind = np.logical_and(x>=x_limits[0], x<=x_limits[1])
    else:
        data_ind = np.ones(x.shape, dtype=bool)

    x = x[data_ind]
    for ind, col in enumerate(df.columns):
        if ind<start_index:
            continue

        if skip is not 0:
            if not (ind+start_index) % skip == 0:
                continue
        plt.plot(x, np.array(df[col])[data_ind]+(ind-start_index) * sep, 'k-')
        plt.text((np.max(x)-np.min(x))*text_pos[0]+np.min(x), (ind-start_index)*sep+text_pos[1], "{:.2f}".format(float(col)))

    plt.ylim(-0.5, sep*(len(df.columns)-start_index)+2)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.savefig(filename, dpi=dpi)
    plt.close()
