import matplotlib.pyplot as plt
import numpy as np

from dft_functions import *


fpath = "../example/"
outpath = "../output/"

struc = ParseQEO()
struc.Process(fpath + "m1.scf.out")
Ef = struc.Fermi

bands = bands_from_gnu(fpath + "bands.m1.dat.gnu")

fig, ax = plt.subplots(figsize=(10, 8))
for i in bands:  # Here we plots the bands
    ax.plot(i[:, 0], i[:, 1] - Ef, color="black", linewidth=3)

vlines = np.linspace(0, 1, 12) * max(i[:, 0])

for vline in vlines:
    ax.axvline(x=vline, ymin=0, ymax=100, linewidth=1.2, color="black", alpha=0.5)

xlabeltext = [
    r"${\Gamma}$",
    "Y",
    "C",
    "Z",
    r"${\Gamma}$",
    "A",
    "E",
    "Z",
    r"${\Gamma}$",
    "B",
    "D",
    "Z",
]

ax.set_xticks(ticks=vlines, labels=xlabeltext)
ax.tick_params(axis="both", which="major", labelsize=18)
ax.axhline(0, linestyle="--", color="k")

ax.set_xlim(0, max(i[:, 0]))
ax.set_ylim(-1, 2.5)

ax.set_ylabel(r"E - E$_F$ (eV)", fontsize=18, labelpad=-10)

plt.savefig(outpath + "bands_m1.png", dpi=150)
plt.show()
