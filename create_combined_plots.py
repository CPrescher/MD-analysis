# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'

import os
import pandas as pd

from utility import create_df_stack_plot, create_df_surface_plot


output_folder = "../MD-4900at/results"

sq_series_df = pd.read_csv(os.path.join(output_folder, "sq_series.csv"), index_col=0)
rdf_series_df = pd.read_csv(os.path.join(output_folder, "rdf_series.csv"), index_col=0)

create_df_stack_plot(sq_series_df, os.path.join(output_folder, "sq_series_stack.png"),
                     start_index=3, skip=2, x_limits=(0,12), text_pos=(0.9, 1.08),
                     x_label="q $(\AA^{-1})$", y_label="S(q)")

create_df_stack_plot(rdf_series_df, os.path.join(output_folder, "rdf_series_stack.png"),
                     start_index=3, skip=2, x_label="r $(\AA)$", y_label="g(r)", text_pos=(0.9, 1.08))

create_df_surface_plot(sq_series_df, os.path.join(output_folder, "sq_series_surface.png"), start_index=2,
                       level_limits=(-0.2, 2))


create_df_surface_plot(rdf_series_df,
                       os.path.join(output_folder, "rdf_series_surface.png"),
                       start_index=3,
                       skip=2,
                       level_limits=(-0.2, 2.3),
                       x_label="r $(\AA})$",
                       y_label="Pressure (GPa)",
                       z_label="g(r)")
create_df_surface_plot(rdf_series_df,
                       os.path.join(output_folder, "rdf_series_surface_detail.png"),
                       start_index=3,
                       skip=2,
                       level_limits=(.4, 1.5),
                       x_limits=(2,6),
                       x_label="r $(\AA^{-1})$",
                       y_label="Pressure (GPa)",
                       z_label="g(r)")