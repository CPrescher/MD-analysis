##
# C

import os

import pandas as pd
import numpy as np

from lib.utility import convert_atoms_per_cubic_angstrom_to_density, create_folder
from lib.output_data import analyze_output

from config import folder_name, output_folder, minimum_time_for_average, number_of_atoms, elemental_abundances

files = os.listdir(folder_name)
output_path = create_folder(output_folder, "conditions")


id = []
pressure = []
temperature = []
volume = []
atomic_density = []
real_density = []

for filename in files:
    if "OUTPUT-" in filename:
        print(filename)
        df = analyze_output(os.path.join(folder_name, filename),
                            show=False,
                            output_dir=output_path)

        df_subset = df[df["time (ps)"] > minimum_time_for_average]

        id.append(int(filename[7:]))
        pressure.append(np.mean(df_subset["pressure (kbar)"]))
        temperature.append(np.mean(df_subset["temperature (K)"]))
        volume.append(np.mean(df_subset["volume (A^3)"]))
        atomic_density.append(float(number_of_atoms)/np.mean(df_subset["volume (A^3)"]))

df = pd.DataFrame()
df["id"] = id
df["pressure (GPa)"] = np.array(pressure) / 10.0
df["temperature (K)"] = temperature
df["volume (A^3)"] = volume
df["atomic_density"] = atomic_density
df["density"] = convert_atoms_per_cubic_angstrom_to_density(elemental_abundances, np.array(atomic_density))
df = df.set_index(["id"])
df = df.sort_index()

df.to_csv(os.path.join(output_folder, 'conditions.csv'))
