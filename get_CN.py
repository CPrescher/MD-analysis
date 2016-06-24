import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import output_folder

rdf_folder = os.path.join(output_folder, "atomic_rdf")

conditions_df = pd.read_csv(os.path.join(output_folder, "conditions.csv"),
                            index_col=0)


pressure = []
CN = []
cutoff = 2.1
for id, row in conditions_df.iterrows():
    filename = "RDFDAT-{}_atomic_rdf.csv".format(id)
    file_path = os.path.join(rdf_folder, filename)
    print(file_path)
    rdf_df = pd.read_csv(file_path)

    plt.plot(rdf_df['r'], rdf_df['SI4+ - O_2-'])
    ind = np.where(rdf_df['r']<cutoff)[0]

    r_int = rdf_df['r'][ind]
    rdf_int = rdf_df['SI4+ - O_2-'][ind]

    CN.append(np.trapz(r_int**2*rdf_int, r_int)*4*np.pi*row['atomic_density']*2/3)
    pressure.append(row['pressure (GPa)'])


plt.figure()
plt.plot(pressure, CN, 'ko')
plt.show()
