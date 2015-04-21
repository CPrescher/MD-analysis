# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'

import os
import pandas as pd

from utility import create_df_stack_plot, create_df_surface_plot


output_folder = "../MD-4900at/results"

sq_series_df = pd.read_csv(os.path.join(output_folder, "sq_series.csv"), index_col=0)
rdf_series_df = pd.read_csv(os.path.join(output_folder, "rdf_series.csv"), index_col=0)

create_df_stack_plot(sq_series_df, os.path.join(output_folder, "sq_series_stack.png"),
                     skip=2, x_limits=(0,12))
create_df_stack_plot(rdf_series_df, os.path.join(output_folder, "rdf_series_stack.png"), skip=2)

create_df_surface_plot(sq_series_df, os.path.join(output_folder, "sq_series_surface.png"), skip=2,
                       level_limits=(-0.2, 2))
create_df_surface_plot(rdf_series_df, os.path.join(output_folder, "rdf_series_surface.png"), skip=2,
                       level_limits=(-0.2, 2.3))