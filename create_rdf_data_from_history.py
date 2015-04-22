# -*- coding: utf8 -*-
import os
import datetime

__author__ = 'Clemens Prescher'

from lib.utility import create_folder
from lib.history_data import read_history_data, calculate_rdf_all

#PARAMETER
folder_name = "../MD-4900at"
file_numbers = (1, 46)

output_folder = "../MD-4900at/results"

r_limits = (0.1, 15)
r_step = 0.05
##################


output_path = create_folder(output_folder, "atomic_rdf_own")

files = os.listdir(folder_name)

for id in range(file_numbers[0], file_numbers[1]+1):
    filename = "HISTORY-{}".format(id)
    print("{}: Analyzing {}".format(datetime.datetime.now(), filename))
    positions, side_length, number_of_atoms=read_history_data(os.path.join(folder_name, filename))
    rdf_df = calculate_rdf_all(positions, side_length, number_of_atoms, r_limits, r_step)
    rdf_df.to_csv(os.path.join(output_path, "atomic_rdf-{}".format(id)))


