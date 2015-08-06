# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'
import unittest

from lib.output_data import convert_time_str

class HistoryDataTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_convert_time_str(self):
        self.assertEqual(convert_time_str("800.000f"), 0.8)
        self.assertEqual(convert_time_str("1.500p"), 1.5)

        self.assertEqual(convert_time_str("39.486m"), 2369.16)
        self.assertEqual(convert_time_str("19.200s"), 19.2)

