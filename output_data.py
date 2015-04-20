import matplotlib.pyplot as plt
import pandas as pd
import os


def is_separator(line_str):
    if " --------------------------" in line_str:
        return True
    return False


def is_empty(line_str):
    if line_str is "" or line_str is " ":
        return True
    return False


def analyze_output(filename, show=True, save_figure=True, save_csv=True, output_dir=None):

    if output_dir is None:
        output_dir = os.path.dirname(filename)

    file = open(filename, 'r')
    filename = os.path.basename(filename)
    # read in the file into an line array
    file_str = []
    for line in file:
        file_str.append(line)



    step = []
    time = []
    cpu_time = []
    energy = []
    volume = []
    temperature = []
    pressure = []

    line_index = 0

    while line_index < len(file_str):
        if is_separator(file_str[line_index]):
            line_index += 1
            if not is_empty(file_str[line_index]) and not is_separator(file_str[line_index]):
                split_line = file_str[line_index].split()
                if len(split_line):
                    step.append(int(split_line[0]))
                    energy.append(float(split_line[1]))
                    temperature.append(float(split_line[2]))
                else:
                    continue

                line_index += 1
                split_line = file_str[line_index].split()
                if len(split_line):
                    time.append(float(split_line[0]))

                line_index += 1
                split_line = file_str[line_index].rstrip().split()
                if len(split_line):
                    cpu_time.append(float(split_line[0]))
                    volume.append(float(split_line[1]))
                    pressure.append(float(split_line[9]))

        line_index += 1

    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(time, pressure)
    plt.xlabel("time (ps)")
    plt.ylabel("Pressure (kbar)")
    plt.subplot(2, 2, 2)
    plt.plot(time, temperature)
    plt.xlabel("time (ps)")
    plt.ylabel("Temperature (K)")
    plt.subplot(2,2,3)
    plt.plot(time, energy)
    plt.xlabel("time (ps)")
    plt.ylabel("Energy")
    plt.subplot(2,2,4)
    plt.plot(time, volume)
    plt.xlabel("time (ps)")
    plt.ylabel("Volume (A^3)")
    plt.tight_layout()

    if save_figure:
        plt.savefig(os.path.join(output_dir,filename+'.png'), dpi=300)

    if show:
        plt.show()

    plt.close()

    df = pd.DataFrame()
    df["step"] = step
    df["cpu_time (s)"] = cpu_time
    df["time (ps)"] = time
    df["temperature (K)"] = temperature
    df["pressure (kbar)"] = pressure
    df["volume (A^3)"] = volume
    df["energy"] = energy

    if save_csv:
        df.to_csv(os.path.join(output_dir,filename+'.csv'))

    return df


if __name__ == "__main__":
    print analyze_output("OUTPUT-2", show=False)