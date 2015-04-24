# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'

import unittest
import os
import time

from lib.history_data import read_history_data
from lib.history_data import read_time_and_atom_num, read_side_length, read_position, read_atoms
from lib.history_data import calculate_rdf

test_directory = os.path.dirname(__file__)
data_directory = os.path.join(test_directory, 'data')
output_directory = os.path.join(test_directory, 'output')
hist_file_path = os.path.join(data_directory, "HISTORY")

class TestHistoryData(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_time_and_atom_num(self):
        example_str = "timestep      1000      4900         0         3    0.001000"
        time, atom_num = read_time_and_atom_num(example_str)
        self.assertEqual(time, 1000)
        self.assertEqual(atom_num, 4900)

    def test_read_side_length(self):
        example_str = "   35.82       0.000       0.000 "
        self.assertEqual(35.82, read_side_length(example_str))

    def test_read_position(self):
        example_str = " -1.3596E+01  8.4633E-01  5.0859E+00"
        self.assertEqual([-13.596, 0.84633, 5.0859], read_position(example_str))

    def test_read_atoms(self):
        example_str = """MG2+             1   24.305000    1.200000
 -1.3596E+01  8.4633E-01  5.0859E+00
MG2+             2   24.305000    1.200000
 -7.2111E+00 -1.6219E-02  1.6485E+00
SI4+             3   24.305000    1.200000
  1.6648E+01 -9.2389E+00  5.7718E+00
O_2+             4   24.305000    1.200000
  7.7851E+00  2.8863E+00  6.4479E+00"""
        example_str = example_str.split("\n")
        atoms = read_atoms(example_str)
        self.assertEqual(len(atoms), 3)
        self.assertEqual(len(atoms["MG2+"]), 2)
        self.assertEqual(len(atoms["SI4+"]), 1)

    def test_read_history_data(self):
        positions, side_length, number_of_atoms = read_history_data(hist_file_path)
        self.assertEqual(len(positions[0]), 3)
        self.assertEqual(number_of_atoms[0], 4900)
        self.assertEqual(side_length[0], 35.82)

    def test_calculate_rdf(self):
        positions, side_length, number_of_atoms = read_history_data(hist_file_path)
        t1 = time.time()
        df = calculate_rdf(positions[0], number_of_atoms[0], float(side_length[0]))
        df.to_csv(os.path.join(output_directory, "single_step.csv"))

        print("Calculation of one step took: {} s", time.time()-t1 )



