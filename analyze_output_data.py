import os
from output_data import analyze_output

folder_name = "../MD-4900at"
output_folder= "../MD-4900at/results"

files = os.listdir(folder_name)

for filename in files:
    if "OUTPUT-" in filename:
        print filename
        analyze_output(os.path.join(folder_name, filename),
                       show=False,
                       output_dir=output_folder)
