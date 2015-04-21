import os

from utility import create_folder, create_df_plot
from rdf_data import read_rdf_data
from rdf_data import calculate_atomic_sq, calculate_combined_sq
from rdf_data import calculate_weighted_fr, calculate_weighted_rdf


# script parameter
folder_name = "../MD-4900at"
output_folder = "../MD-4900at/results"
minimum_time_for_average = 5
number_of_atoms = 4900
elemental_abundances = {"Mg": 2,
                        "Si": 1,
                        "O": 4}

files = os.listdir(folder_name)

# create output folders
atomic_rdf_path = create_folder(output_folder, "atomic_rdf")
weighted_fr_path = create_folder(output_folder, "weighted_fr")
weighted_rdf_path = create_folder(output_folder, "weighted_rdf")
atomic_sq_path = create_folder(output_folder, "atomic_sq")
weighted_sq_path = create_folder(output_folder, "combined_sq")


# for filename in files:
# if "RDFDAT-" in filename:
# print filename
# for file in files:
# if file.startswith("RDFDAT") and not file.endswith(".png"):


filename = "RDFDAT-1"
atomic_density = .05
# reading the data
rdf_df = read_rdf_data(os.path.join(folder_name, filename))
rdf_df.to_csv(os.path.join(atomic_rdf_path, filename + "_atomic_rdf.csv"))

sq_df = calculate_atomic_sq(rdf_df, atomic_density=atomic_density)
sq_df.to_csv(os.path.join(atomic_sq_path, filename + "_atomic_sq.csv"))

sq_weighted_df = calculate_combined_sq(sq_df, elemental_abundances)
sq_weighted_df.to_csv(os.path.join(weighted_sq_path, filename + "_combined_sq.csv"))

fr_weighted_df = calculate_weighted_fr(sq_weighted_df)
fr_weighted_df.to_csv(os.path.join(weighted_fr_path, filename + '_weighted_fr.csv'))

rdf_weighted_df = calculate_weighted_rdf(fr_weighted_df, atomic_density=atomic_density)
rdf_weighted_df.to_csv(os.path.join(weighted_rdf_path, filename + '_weighted_rdf.csv'))


# creating the plots
create_df_plot(rdf_df, os.path.join(atomic_rdf_path, filename + "_atomic_rdf.png"))
create_df_plot(sq_df, os.path.join(atomic_sq_path, filename + "_atomic_sq.png"))
create_df_plot(sq_weighted_df, os.path.join(weighted_sq_path, filename + "_weighted_sq.png"))
create_df_plot(fr_weighted_df, os.path.join(weighted_fr_path, filename + "_weighted_fr.png"), columns="sum")
create_df_plot(rdf_weighted_df, os.path.join(weighted_rdf_path, filename + "_weighted_rdf.png"), columns="sum")



