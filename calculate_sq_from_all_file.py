import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d

filename = "RDFDAT"

data = np.loadtxt(filename, skiprows=3)

number_density = 0.087

r = data[:, 0]
g_r = data[:, 1]

# plt.plot(r, g_r)
g_r = gaussian_filter1d(g_r, 3)
plt.plot(r, g_r)
plt.show()
q = np.linspace(0, 25, 1000)

int_val = []
for val in q:
    int_val.append(np.trapz(r * np.sin(val * r) * (g_r - 1), r))

s_q = 1 + 4 * np.pi * number_density / q * int_val

plt.figure()
plt.plot(q, s_q)
plt.tight_layout()
plt.show()
# plt.savefig("sq_output.png", dpi=300)
