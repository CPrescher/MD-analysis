# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'

from multiprocessing import Pool

import pandas as pd
import numpy as np


def read_history_data(filename, start_line=3):
    file = open(filename, "r")

    file_str = []
    for line in file:
        file_str.append(line)

    ind = start_line - 1

    positions = []
    number_of_atoms = []
    side_length = []

    while(ind<len(file_str)):
        if file_str[ind].split()[0] == "File":
            ind += 2
        time_step, atoms = read_time_and_atom_num(file_str[ind])
        number_of_atoms.append(atoms)
        side_length.append(read_side_length(file_str[ind + 1]))
        positions.append(read_atoms(file_str[ind+4:ind+4+2*atoms]))
        ind+=atoms*2+4

    return positions, side_length, number_of_atoms

def calculate_rdf_all(position_list, side_length, number_of_atoms, r_limits=(0.01, 10), r_step=0.01):
    # create a list of arguments for the Pool
    argument_list = []
    for ind, positions in enumerate(position_list):
        argument_list.append((positions, number_of_atoms[ind], side_length[ind], r_limits, r_step))

    pool = Pool(processes=16)
    rdf_dataframes = pool.starmap(calculate_rdf, argument_list)

    # create average for each step
    rdf_df = rdf_dataframes[0]

    #sum all data frames
    for i in range(1, len(rdf_dataframes)):
        for col in rdf_df.columns:
            rdf_df[col] += rdf_dataframes[i][col]

    #divide by the number of dataframes
    rdf_df = rdf_df / float(len(rdf_dataframes))

    return rdf_df


def calculate_rdf(positions, number_of_atoms, side_length, r_limits=(0.01, 10), r_step=0.01):
    real_r = np.arange(r_limits[0], r_limits[1], r_step)
    hist_r = np.arange(r_limits[0] - r_step / 2, r_limits[1] + r_step / 2, r_step)

    shell_volume = 4 * np.pi * hist_r ** 2 * r_step
    shell_volume = shell_volume[:-1]

    volume = side_length**3

    element_num = {}
    for key in positions.keys():
        element_num[key], _ = positions[key].shape

    elements = list(positions.keys())

    df = pd.DataFrame(index=real_r)
    df.index.name="r"

    for i in range(len(elements)):
        for j in range(i, len(elements)):
            distances = np.array([])
            for n in range(element_num[elements[i]]):
                distances = np.concatenate((distances, distance(positions[elements[i]][n, :], positions[elements[j]],
                                                                side_length)))
            hist, _ = np.histogram(distances, bins=hist_r)
            hist = np.array(hist)/float(element_num[elements[i]])
            partial_density = element_num[elements[j]]/volume
            df["{} - {}".format(elements[i], elements[j])] = hist / (shell_volume*partial_density)
    return df


def read_time_and_atom_num(line_str):
    line_str_splitted = line_str.split()
    time_step = int(line_str_splitted[1])
    number_of_atoms = int(line_str_splitted[2])
    return time_step, number_of_atoms


def read_side_length(line_str):
    line_str_splitted = line_str.split()
    side_length = float(line_str_splitted[0])
    return side_length


def read_atoms(file_str):
    positions = {}
    for line_index in range(0, len(file_str), 2):
        atom_type = file_str[line_index].split()[0]
        if atom_type not in positions:
            positions[atom_type] = []
        positions[atom_type].append(read_position(file_str[line_index + 1]))

    for key in positions.keys():
        positions[key] = np.array(positions[key])
    return positions


def distance(center_point, distant_points, dimensions):
    delta = np.abs(center_point - distant_points)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


def read_position(line_str):
    line_str_splitted = line_str.split()
    return [float(line_str_splitted[0]), float(line_str_splitted[1]), float(line_str_splitted[2])]