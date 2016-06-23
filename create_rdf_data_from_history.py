# -*- coding: utf8 -*-
import os
import datetime

__author__ = 'Clemens Prescher'

from lib.utility import create_folder, create_df_plot
from lib.history_data import read_history_data, calculate_rdf_all

#PARAMETER
folder_name = "../Simulation"
output_folder = "../RDF_history"

#os.mkdir(output_folder)


file_numbers = (1, 46)
r_limits = (0.1, 15)
r_step = 0.05

for id in range(file_numbers[0], file_numbers[1]+1):
    filename = "HISTORY-{}".format(id)
    print("{}: Analyzing {}".format(datetime.datetime.now(), filename))
    positions, side_length, number_of_atoms=read_history_data(os.path.join(folder_name, filename))
    rdf_df = calculate_rdf_all(positions, side_length, number_of_atoms, r_limits, r_step)
    rdf_df.to_csv(os.path.join(output_folder, "atomic_rdf-{}.csv".format(id)))
    create_df_plot(rdf_df, os.path.join(output_folder,"atomic_rdf-{}.eps".format(id)), legend_loc="best")


