# -*- coding: utf8 -*-
__author__ = 'Clemens Prescher'

import os
import pandas as pd

from lib.utility import create_df_stack_plot, create_df_surface_plot, create_folder
from config import output_folder, stack_index, plot_start_index, plot_skip, plot_gr_text_pos, plot_sq_text_pos


output_path = create_folder(output_folder, "combined_plots")


sq_series_df = pd.read_csv(os.path.join(output_folder, "sq_series.csv"), index_col=0)
rdf_series_df = pd.read_csv(os.path.join(output_folder, "rdf_series.csv"), index_col=0)

#############################################################################
### Stack Plots
#############################################################################

create_df_stack_plot(sq_series_df,
                     os.path.join(output_path, "sq_series_stack.png"),
                     start_index=plot_start_index,
                     skip=plot_skip,
                     x_limits=(0,12),
                     text_pos=plot_sq_text_pos,
                     x_label="q $(\AA^{-1})$",
                     y_label="S(q)")

create_df_stack_plot(rdf_series_df,
                     os.path.join(output_path, "rdf_series_stack.png"),
                     start_index=plot_start_index,
                     skip=plot_skip,
                     x_label="r $(\AA)$",
                     y_label="g(r)",
                     text_pos=plot_gr_text_pos)

#############################################################################
### Surface Plots
#############################################################################

create_df_surface_plot(sq_series_df,
                       os.path.join(output_path, "sq_series_surface.png"),
                       start_index=plot_start_index,
                       level_limits=(-0.2, 2),
                       x_label="r $(\AA^{-1})$",
                       y_label=stack_index,
                       skip=plot_skip
                       )

create_df_surface_plot(rdf_series_df,
                       os.path.join(output_path, "rdf_series_surface.png"),
                       start_index=plot_start_index,
                       skip=plot_skip,
                       level_limits=(-0.2, 2.3),
                       x_label="r $(\AA})$",
                       y_label=stack_index,
                       z_label="g(r)")

create_df_surface_plot(rdf_series_df,
                       os.path.join(output_path, "rdf_series_surface_detail.png"),
                       start_index=plot_start_index,
                       skip=plot_skip,
                       level_limits=(.4, 1.5),
                       x_limits=(2,6),
                       x_label="r $(\AA^{-1})$",
                       y_label=stack_index,
                       z_label="g(r)")