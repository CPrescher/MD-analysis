import os

import pandas as pd
import numpy as np

from lib.utility import create_folder, create_df_plot
from lib.rdf_data import read_rdf_data
from lib.rdf_data import calculate_atomic_sq, calculate_combined_sq
from lib.rdf_data import calculate_weighted_fr, calculate_weighted_rdf

from config import folder_name, output_folder, elemental_abundances, stack_index

# create output folders
atomic_rdf_path = create_folder(output_folder, "atomic_rdf")
weighted_fr_path = create_folder(output_folder, "weighted_fr")
weighted_rdf_path = create_folder(output_folder, "weighted_rdf")
atomic_sq_path = create_folder(output_folder, "atomic_sq")
weighted_sq_path = create_folder(output_folder, "combined_sq")

conditions_df = pd.read_csv(os.path.join(output_folder, "conditions.csv"),
                            index_col=0)

sq_weighted_series_df = pd.DataFrame()
rdf_weighted_series_df = pd.DataFrame()

files = os.listdir(folder_name)
output_path = create_folder(output_folder, "conditions")

for id, row in conditions_df.iterrows():
    filename = "RDFDAT-{}".format(id)
    print(filename)
    atomic_density = np.array(conditions_df["atomic_density"])[id-1]
    stack_value = np.array(conditions_df[stack_index])[id - 1]


    rdf_df = read_rdf_data(os.path.join(folder_name, filename))
    rdf_df.to_csv(os.path.join(atomic_rdf_path, filename + "_atomic_rdf.csv"))

    sq_df = calculate_atomic_sq(rdf_df, atomic_density=atomic_density)
    sq_df.to_csv(os.path.join(atomic_sq_path, filename + "_atomic_sq.csv"))

    sq_weighted_df = calculate_combined_sq(sq_df, elemental_abundances)
    sq_weighted_df.to_csv(os.path.join(weighted_sq_path, filename + "_combined_sq.csv"))

    fr_weighted_df = calculate_weighted_fr(sq_weighted_df, q_max=12)
    fr_weighted_df.to_csv(os.path.join(weighted_fr_path, filename + '_weighted_fr.csv'))

    rdf_weighted_df = calculate_weighted_rdf(fr_weighted_df, atomic_density=atomic_density)
    rdf_weighted_df.to_csv(os.path.join(weighted_rdf_path, filename + '_weighted_rdf.csv'))


    # creating the plots
    create_df_plot(rdf_df, os.path.join(atomic_rdf_path, filename + "_atomic_rdf.png"))
    create_df_plot(sq_df, os.path.join(atomic_sq_path, filename + "_atomic_sq.png"))
    create_df_plot(sq_weighted_df, os.path.join(weighted_sq_path, filename + "_weighted_sq.png"))
    create_df_plot(fr_weighted_df, os.path.join(weighted_fr_path, filename + "_weighted_fr.png"), columns="sum")
    create_df_plot(rdf_weighted_df, os.path.join(weighted_rdf_path, filename + "_weighted_rdf.png"), columns="sum")

    sq_weighted_series_df["q"] = np.array(sq_weighted_df.index.values)
    sq_weighted_series_df[str(stack_value)] = np.array(sq_weighted_df["sum"])

    rdf_weighted_series_df["r"] = rdf_weighted_df.index.values
    rdf_weighted_series_df[str(stack_value)] = np.array(rdf_weighted_df["sum"])


sq_weighted_series_df = sq_weighted_series_df.set_index("q")
rdf_weighted_series_df = rdf_weighted_series_df.set_index("r")

sq_weighted_series_df.to_csv(os.path.join(output_folder, "sq_series.csv"))
rdf_weighted_series_df.to_csv(os.path.join(output_folder, "rdf_series.csv"))






