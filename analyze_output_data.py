import os
from output_data import analyze_output

files = os.listdir(".")

for file in files:
    if "OUTPUT-" in file:
        num = int(file.split('-')[1])
        analyze_output(num, show=False)
