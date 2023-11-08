import matplotlib.pyplot as plt

from dft_functions import *


fpath = "../example/"
outpath = "../output/"

struc = ParseQEO()
struc.Process(fpath + "m1.scf.out")
Ef = struc.Fermi

pdosname = fpath + "pdos/m1"
outname = fpath + "m1.scf"
energy, pdos, pdos_tot = vo2_dict3(
    pdosname, outname, n_cells=1, E_fermi=Ef, spin=False, sd_reversal=True
)


fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10, 2 * 4))

lw = 1.2

axs[0].plot(
    energy,
    pdos["v_pdos_dx2y2"][0],
    linewidth=2,
    linestyle="-",
    color="k",
    label=r"V1 $d_{x^2 - y^2}$",
)
axs[0].plot(
    energy,
    pdos["v_pdos_dzx"][0],
    linewidth=2,
    linestyle="-",
    color="c",
    label=r"V1 $d_{xz}$",
)
axs[0].plot(
    energy,
    pdos["v_pdos_dzy"][0],
    linewidth=2,
    linestyle="-",
    color="m",
    label=r"V1 $d_{yz}$",
)
axs[0].plot(
    energy,
    pdos["v_pdos_dxy"][0],
    linewidth=2,
    linestyle="-",
    color="orange",
    label=r"V1 $d_{xy}$",
)
axs[0].plot(
    energy,
    pdos["v_pdos_dz2"][0],
    linewidth=2,
    linestyle="-",
    color="y",
    label=r"V1 $d_{z^2}$",
)

axs[1].plot(
    energy,
    pdos["v_pdos_dx2y2"][2],
    linewidth=2,
    linestyle="-",
    color="k",
    label=r"V2 $d_{x^2 - y^2}$",
)
axs[1].plot(
    energy,
    pdos["v_pdos_dzx"][2],
    linewidth=2,
    linestyle="-",
    color="c",
    label=r"V2 $d_{xz}$",
)
axs[1].plot(
    energy,
    pdos["v_pdos_dzy"][2],
    linewidth=2,
    linestyle="-",
    color="m",
    label=r"V2 $d_{yz}$",
)
axs[1].plot(
    energy,
    pdos["v_pdos_dxy"][2],
    linewidth=2,
    linestyle="-",
    color="orange",
    label=r"V2 $d_{xy}$",
)
axs[1].plot(
    energy,
    pdos["v_pdos_dz2"][2],
    linewidth=2,
    linestyle="-",
    color="y",
    label=r"V2 $d_{z^2}$",
)

for ax in axs:
    ax.set_yticks([])
    ax.set_ylabel("DOS", fontsize=18)
    ax.legend(loc=2, frameon=False, fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=16, width=2, length=6)
    ax.axvline(0, linestyle="--", color="k")


axs[-1].set_xlim(-8, 7)
axs[-1].set_ylim(ymax=2.5)
axs[-1].set_xlabel("Energy (eV)", fontsize=18)

plt.savefig(outpath + "pdos_m1.png", dpi=150)
plt.show()
