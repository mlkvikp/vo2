# %%
import numpy as np
import re
import json
import itertools
import shutil, os, glob
from subprocess import run
import numpy as np


# %%
class ParseQEO:
    def __init__(self):
        self.BOHRtoA = 0.529177249
        self.RYtoeV = 13.605698066
        self.program_version = ""
        self.volume = 0.0
        self.alat = 0.0
        self.natoms = 0
        self.nat = 0
        self.nelect = 0
        self.Ecut = 0.0
        self.RhoCut = 0.0
        self.Econv = 0.0
        self.beta = 0.0
        self.Exch = ""
        self.energy = 0.0
        self.natoms = 0
        self.bandgap = 0.0
        self.bands = 0
        self.lattice = {"a": np.zeros(3), "b": np.zeros(3), "c": np.zeros(3)}
        self.atoms = []
        self.norms = {"a": 0.0, "b": 0.0, "c": 0.0}
        self.angles = {"alpha": 0.0, "beta": 0.0, "gamma": 0.0}
        self.kpts = 0
        self.bnddiagram = np.zeros(0)
        self.symweights = np.zeros(0)
        self.FermiTest = False
        self.Fermi = 0.0
        self.U_dict = {}  # dictionary for Hubbard U calculations

    def Process(self, filestring, spin_polarize=False, Hubbard_U=False, Debug=False):
        f = open(filestring, "r")
        linenum = 0
        #'i' stands for the line in this script
        for i in f:
            # debug mode
            if Debug:
                print(i)

            if linenum < 1000:
                if "number of k points=" in i:
                    self.kpts = int(i.split()[4])
                    ##new code for DOS
                    self.symweights = np.zeros(self.kpts)
                    next(f)
                    for j in range(self.kpts):
                        line = next(f)
                        # print(line)
                        self.symweights[j] = float(line.split()[9])

                    next
                if "Program PWSCF" in i:
                    self.program_version = i.split()[2]
                    next
                if "lattice parameter (alat)" in i:
                    self.alat = float(i.split()[4]) * self.BOHRtoA
                    next
                if "number of Kohn-Sham states" in i:
                    self.bands = int(i.split()[4])
                if "unit-cell volume" in i and "new" not in i:
                    self.volume = float(i.split()[3]) * (self.BOHRtoA**3.0)
                    next
                if "number of atoms/cell" in i:
                    self.natoms = int(i.split()[4])
                    next
                if "number of atomic types" in i:
                    self.nat = int(i.split()[5])
                    next
                if "number of electrons" in i:
                    self.nelect = float(i.split()[4])
                    next
                if "kinetic-energy cutoff" in i:
                    self.Ecut = float(i.split()[3]) * self.RYtoeV
                    next
                if "charge density cutoff" in i:
                    self.RhoCut = float(i.split()[4]) * self.RYtoeV
                    next
                # if "convergence threshold" in i:
                #   self.Econv = float(i.split()[3])
                #   next
                if "mixing beta" in i:
                    self.beta = float(i.split()[3])
                    next
                if "Exchange-correlation" in i:
                    self.Exch = i.split()[0]  # should be [2]
                    next
                if "a(1) =" in i:
                    tmp = i.split()
                    for j in range(0, 3):
                        self.lattice["a"][j] = tmp[j + 3]
                    next
                if "a(2) =" in i:
                    tmp = i.split()
                    for j in range(0, 3):
                        self.lattice["b"][j] = tmp[j + 3]
                    next
                if "a(3) =" in i:
                    tmp = i.split()
                    for j in range(0, 3):
                        self.lattice["c"][j] = tmp[j + 3]
                    next
                if "site n.     atom                  positions (alat units)" in i:
                    for j in range(0, self.natoms):
                        line = next(f).split()
                        self.atoms.append(
                            [
                                line[1],
                                float(line[6]) * self.alat,
                                float(line[7]) * self.alat,
                                float(line[8]) * self.alat,
                            ]
                        )
                    next
            if "!" in i:
                self.energy = float(i.split()[4]) * self.RYtoeV
            if "new unit-cell volume" in i:
                self.volume = float(i.split()[4]) * (self.BOHRtoA**3)

            if "Begin final coordinates" in i:
                while "End final coordinates" not in line:
                    line = next(f)
                    if "CELL_PARAMETERS" in line:
                        for j in ["a", "b", "c"]:
                            line = next(f)
                            tmp = line.split()
                            for k in range(0, 3):
                                self.lattice[j][k] = float(tmp[k])
                    if "ATOMIC_POSITIONS" in line:
                        if "(crystal)" in line:
                            for j in range(0, self.natoms):
                                line = next(f).split()
                                self.atoms[j] = [
                                    line[0],
                                    float(line[1]),
                                    float(line[2]),
                                    float(line[3]),
                                ]
                        # if "angstrom" in line:
                        #   for j in range(0,self.natoms):
                        #     line = next(f).split()
                        #     self.atoms[j] = [line[0],float(line[1]),float(line[2]),float(line[3])]
            if Hubbard_U:
                if "Tr[ns(na)] (up, down, total)" in i:
                    this_line = i.split()
                    natom = int(this_line[1])
                    up_num = float(this_line[-3])
                    down_num = float(this_line[-2])
                    tot_num = float(this_line[-1])

                    self.U_dict[f"atom {natom}"] = {
                        "natom": natom,
                        "up_num": up_num,
                        "down_num": down_num,
                        "tot_num": tot_num,
                    }

                    for spin_idx in ["up", "down"]:
                        # print(f'Processing spin {spin_idx} for atom {natom}')
                        next(f), next(f)
                        line = next(f)
                        eigs = line.split()
                        self.U_dict[f"atom {natom}"][f"{spin_idx} occupations"] = eigs

                        for skipline in range(
                            len(eigs) * 1 + 2
                        ):  # temporary, skip all the eigenvectors and the occupation matrix plus the headers
                            next(f)

                        temp_mat = np.zeros((len(eigs), len(eigs)))
                        for mat_idx in range(len(eigs)):
                            line = next(f)
                            # print (line)
                            temp_mat[mat_idx] = line.split()
                        self.U_dict[f"atom {natom}"][
                            f"{spin_idx} occupation matrix"
                        ] = temp_mat

                    this_line = next(f).split()
                    magmom = float(this_line[-1])
                    self.U_dict[f"atom {natom}"]["local magnetic moment"] = magmom
                    self.U_dict[f"atom {natom}"]["d states"] = [
                        "r^2-3*z^2",
                        "xz",
                        "yz",
                        "xy",
                        "x^2-y^2",
                    ]

            """
      band_condition = (
        "End of self-consistent calculation" in i or
        "End of band structure calculation" in i 
      )
      

      if band_condition: #The band diagram is stored in lines of 8 entries
        if np.floor(self.bands/8.)*8. <= self.bands:
          numlines = int(np.floor(self.bands/8.) + 1)
          remainder = int(self.bands - np.floor(self.bands/8.)*8.)
        else: 
          numlines = int(np.floor(self.bands/8.))
          remainder = 0
        
        if spin_polarize == True:
          self.bnddiagram =  np.zeros((self.kpts,self.bands, 2))
          spin_iter = [0,1]
        else:
          spin_iter = [0]
          self.bnddiagram = np.zeros((self.kpts,self.bands,1))
        
        for polarization in spin_iter:
          counter = 0
          while counter < self.kpts:
            # if 'SPIN' in line:
              # print('Parsing: ' + line) 
            
            line = next(f)
            
            if "k =" in line:
              line = next(f)
              counter1 = 0
              for j in range(0,numlines):
                line = next(f)
                for k in range(0,len(line.split())):
                  self.bnddiagram[counter][counter1 + k][polarization] = float(line.split()[k])
                counter1 += 8
              counter += 1
             
          next
        next  
      """

            if "highest occupied, lowest unoccupied level (ev)" in i:
                self.bandgap = float(i.split()[7]) - float(i.split()[6])
                next
            if "the Fermi energy is" in i:
                self.Fermi = float(i.split()[4])
                self.FermiTest = True
                next
            linenum += 1
        f.close()
        for i in ["a", "b", "c"]:
            self.norms[i] = np.linalg.norm(self.lattice[i])
        self.angles["alpha"] = (
            np.arccos(
                np.dot(self.lattice["b"], self.lattice["c"])
                / (self.norms["c"] * self.norms["b"])
            )
            * 180.0
            / np.pi
        )
        self.angles["gamma"] = (
            np.arccos(
                np.dot(self.lattice["a"], self.lattice["b"])
                / (self.norms["a"] * self.norms["b"])
            )
            * 180.0
            / np.pi
        )
        self.angles["beta"] = (
            np.arccos(
                np.dot(self.lattice["a"], self.lattice["c"])
                / (self.norms["a"] * self.norms["c"])
            )
            * 180.0
            / np.pi
        )
        if self.FermiTest == True:  # The bandgap is now in the band diagram
            self.bnddiagram = np.subtract(self.bnddiagram, self.Fermi)
            emin = np.zeros(self.kpts)
            emax = np.zeros(self.kpts)
            counter = 0
            # NEED TO FIX, BANDGAP IS OBTAINABLE ONLY
            # FOR SINGLE SPIN POLARIZATION
            """
      for j in self.bnddiagram[:,:,0]:
        emin[counter] = j[np.searchsorted(j,  0.0,side='right')-1]
        emax[counter] = j[np.searchsorted(j,  0.0,side='right')]
        counter += 1
      self.bandgap = float(np.min(emax-emin))
      """

    def to_JSON(self, outfile):
        for i in self.lattice:
            self.lattice[i] = self.lattice[i].tolist()
        self.bnddiagram = self.bnddiagram.tolist()  # JSON doesnt like numpy arrays
        return json.dump(
            self, outfile, default=lambda o: o.__dict__, sort_keys=True, indent=4
        )


# %%


# %%
def data_loader(fname, spin=False, j=0):
    # need to implement spinup and spin down separation, rn just summing
    with open(fname) as fid:
        data = fid.readlines()

        energy = []
        pdos = []
        if spin:
            for row in range(len(data)):
                data_rows = data[row]
                if data_rows[0][0] != "#":
                    data_row = data_rows.split()
                    energy.append(float(data_row[0]))
                    pdos.append(float(data_row[3 + 2 * j]) + float(data_row[4 + 2 * j]))
        else:
            for row in range(len(data)):
                data_rows = data[row]
                if data_rows[0][0] != "#":
                    data_row = data_rows.split()
                    energy.append(float(data_row[0]))
                    pdos.append(float(data_row[2 + j]))

        energy = np.asarray(energy)
        pdos = np.asarray(pdos)

    return energy, pdos


# %%
def vo2_pdos(fname, n_cells, E_fermi, spin=False, sd_reversal=False):
    n_O = 4 * n_cells
    n_V = 2 * n_cells

    _, pdos_tot = data_loader(fname + ".pdos_tot", spin)
    l = len(pdos_tot)

    o_pdos_s = np.zeros(l)
    o_pdos_p = np.zeros(l)

    v_pdos_s = np.zeros(l)
    v_pdos_p = np.zeros(l)

    v_pdos_dz2 = np.zeros(l)
    v_pdos_dzx = np.zeros(l)
    v_pdos_dzy = np.zeros(l)
    v_pdos_dx2y2 = np.zeros(l)
    v_pdos_dxy = np.zeros(l)

    v_pdos_4s = np.zeros(l)

    for n in range(1, n_O + 1):
        o_energy, o_pdos_s_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(O)_wfc#1(s)", spin
        )
        _, o_pdos_p_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(O)_wfc#2(p)", spin
        )
        o_pdos_s += o_pdos_s_part
        o_pdos_p += o_pdos_p_part

    if sd_reversal:
        # for some reason different pp's put 4s below 3d in projwfc
        for n in range(n_O + 1, n_O + n_V + 1):
            v_energy, v_pdos_s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#1(s)", spin
            )
            _, v_pdos_p_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#2(p)", spin
            )

            _, v_pdos_dz2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin
            )
            _, v_pdos_dzx_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=1
            )
            _, v_pdos_dzy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=2
            )
            _, v_pdos_dx2y2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=3
            )
            _, v_pdos_dxy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=4
            )

            _, v_pdos_4s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(s)", spin
            )

            v_pdos_s += v_pdos_s_part
            v_pdos_p += v_pdos_p_part

            v_pdos_dz2 += v_pdos_dz2_part
            v_pdos_dzx += v_pdos_dzx_part
            v_pdos_dzy += v_pdos_dzy_part
            v_pdos_dx2y2 += v_pdos_dx2y2_part
            v_pdos_dxy += v_pdos_dxy_part

            v_pdos_4s += v_pdos_4s_part

    else:
        for n in range(n_O + 1, n_O + n_V + 1):
            v_energy, v_pdos_s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#1(s)", spin
            )
            _, v_pdos_p_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#2(p)", spin
            )

            _, v_pdos_dz2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin
            )
            _, v_pdos_dzx_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=1
            )
            _, v_pdos_dzy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=2
            )
            _, v_pdos_dx2y2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=3
            )
            _, v_pdos_dxy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=4
            )

            _, v_pdos_4s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(s)", spin
            )

            v_pdos_s += v_pdos_s_part
            v_pdos_p += v_pdos_p_part

            v_pdos_dz2 += v_pdos_dz2_part
            v_pdos_dzx += v_pdos_dzx_part
            v_pdos_dzy += v_pdos_dzy_part
            v_pdos_dx2y2 += v_pdos_dx2y2_part
            v_pdos_dxy += v_pdos_dxy_part

            v_pdos_4s += v_pdos_4s_part

    energy = v_energy - E_fermi
    pdos = [
        o_pdos_s,
        o_pdos_p,
        v_pdos_s,
        v_pdos_p,
        v_pdos_dz2,
        v_pdos_dzx,
        v_pdos_dzy,
        v_pdos_dx2y2,
        v_pdos_dxy,
        v_pdos_4s,
    ]
    return energy, pdos_tot, pdos


# %%
def vo2_pdos2(fname, n_cells, E_fermi, spin=False, sd_reversal=False):
    n_O = 4 * n_cells
    n_V = 2 * n_cells

    _, pdos_tot = data_loader(fname + ".pdos_tot", spin)
    l = len(pdos_tot)

    o_pdos_s = []
    o_pdos_p = []

    v_pdos_s = []
    v_pdos_p = []

    v_pdos_dz2 = []
    v_pdos_dzx = []
    v_pdos_dzy = []
    v_pdos_dx2y2 = []
    v_pdos_dxy = []

    v_pdos_4s = []

    for n in range(1, n_O + 1):
        o_energy, o_pdos_s_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(O)_wfc#1(s)", spin
        )
        _, o_pdos_p_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(O)_wfc#2(p)", spin
        )
        o_pdos_s.append(o_pdos_s_part)
        o_pdos_p.append(o_pdos_p_part)

    if sd_reversal:
        # for some reason different pp's put 4s below 3d in projwfc
        for n in range(n_O + 1, n_O + n_V + 1):
            v_energy, v_pdos_s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#1(s)", spin
            )
            _, v_pdos_p_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#2(p)", spin
            )

            _, v_pdos_dz2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin
            )
            _, v_pdos_dzx_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=1
            )
            _, v_pdos_dzy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=2
            )
            _, v_pdos_dx2y2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=3
            )
            _, v_pdos_dxy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=4
            )

            _, v_pdos_4s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(s)", spin
            )

            v_pdos_s.append(v_pdos_s_part)
            v_pdos_p.append(v_pdos_p_part)

            v_pdos_dz2.append(v_pdos_dz2_part)
            v_pdos_dzx.append(v_pdos_dzx_part)
            v_pdos_dzy.append(v_pdos_dzy_part)
            v_pdos_dx2y2.append(v_pdos_dx2y2_part)
            v_pdos_dxy.append(v_pdos_dxy_part)

            v_pdos_4s.append(v_pdos_4s_part)

    else:
        for n in range(n_O + 1, n_O + n_V + 1):
            v_energy, v_pdos_s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#1(s)", spin
            )
            _, v_pdos_p_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#2(p)", spin
            )

            _, v_pdos_dz2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin
            )
            _, v_pdos_dzx_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=1
            )
            _, v_pdos_dzy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=2
            )
            _, v_pdos_dx2y2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=3
            )
            _, v_pdos_dxy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=4
            )

            _, v_pdos_4s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(s)", spin
            )

            v_pdos_s.append(v_pdos_s_part)
            v_pdos_p.append(v_pdos_p_part)

            v_pdos_dz2.append(v_pdos_dz2_part)
            v_pdos_dzx.append(v_pdos_dzx_part)
            v_pdos_dzy.append(v_pdos_dzy_part)
            v_pdos_dx2y2.append(v_pdos_dx2y2_part)
            v_pdos_dxy.append(v_pdos_dxy_part)

            v_pdos_4s.append(v_pdos_4s_part)

    energy = v_energy - E_fermi
    pdos = [
        np.array(o_pdos_s),
        np.array(o_pdos_p),
        np.array(v_pdos_s),
        np.array(v_pdos_p),
        np.array(v_pdos_dz2),
        np.array(v_pdos_dzx),
        np.array(v_pdos_dzy),
        np.array(v_pdos_dx2y2),
        np.array(v_pdos_dxy),
        np.array(v_pdos_4s),
    ]
    return energy, pdos_tot, pdos


# %%
def vo2_pdos3(fname, outname, n_cells, E_fermi, spin=False, sd_reversal=False):
    outfname = outname + ".out"

    struc = ParseQEO()
    struc.Process(outfname)

    atoms_names = [j[0] for j in struc.atoms[:]]

    n_at = atoms_names

    _, pdos_tot = data_loader(fname + ".pdos_tot", spin)
    l = len(pdos_tot)

    o_pdos_s = []
    o_pdos_p = []

    v_pdos_s = []
    v_pdos_p = []

    v_pdos_dz2 = []
    v_pdos_dzx = []
    v_pdos_dzy = []
    v_pdos_dx2y2 = []
    v_pdos_dxy = []

    v_pdos_4s = []

    ge_pdos_d = []
    ge_pdos_s = []
    ge_pdos_p = []

    O_ind = [i for i, x in enumerate(atoms_names) if x == "O"]
    O_inds = [x + 1 for x in O_ind]
    V_ind = [i for i, x in enumerate(atoms_names) if x == "V"]
    V_inds = [x + 1 for x in V_ind]
    Ge_ind = [i for i, x in enumerate(atoms_names) if x == "Ge"]
    Ge_inds = [x + 1 for x in Ge_ind]

    for n in O_inds:
        o_energy, o_pdos_s_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(O)_wfc#1(s)", spin
        )
        _, o_pdos_p_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(O)_wfc#2(p)", spin
        )
        o_pdos_s.append(o_pdos_s_part)
        o_pdos_p.append(o_pdos_p_part)

    for n in Ge_inds:
        ge_energy, ge_pdos_d_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(Ge)_wfc#1(d)", spin
        )
        _, ge_pdos_s_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(Ge)_wfc#2(s)", spin
        )
        _, ge_pdos_p_part = data_loader(
            fname + ".pdos_atm#" + str(n) + "(Ge)_wfc#3(p)", spin
        )
        ge_pdos_d.append(ge_pdos_d_part)
        ge_pdos_s.append(ge_pdos_s_part)
        ge_pdos_p.append(ge_pdos_p_part)

    if sd_reversal:
        # for some reason different pp's put 4s below 3d in projwfc
        for n in V_inds:
            v_energy, v_pdos_s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#1(s)", spin
            )
            _, v_pdos_p_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#2(p)", spin
            )

            _, v_pdos_dz2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin
            )
            _, v_pdos_dzx_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=1
            )
            _, v_pdos_dzy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=2
            )
            _, v_pdos_dx2y2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=3
            )
            _, v_pdos_dxy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(d)", spin, j=4
            )

            _, v_pdos_4s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(s)", spin
            )

            v_pdos_s.append(v_pdos_s_part)
            v_pdos_p.append(v_pdos_p_part)

            v_pdos_dz2.append(v_pdos_dz2_part)
            v_pdos_dzx.append(v_pdos_dzx_part)
            v_pdos_dzy.append(v_pdos_dzy_part)
            v_pdos_dx2y2.append(v_pdos_dx2y2_part)
            v_pdos_dxy.append(v_pdos_dxy_part)

            v_pdos_4s.append(v_pdos_4s_part)

    else:
        for n in V_inds:
            v_energy, v_pdos_s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#1(s)", spin
            )
            _, v_pdos_p_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#2(p)", spin
            )

            _, v_pdos_dz2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin
            )
            _, v_pdos_dzx_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=1
            )
            _, v_pdos_dzy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=2
            )
            _, v_pdos_dx2y2_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=3
            )
            _, v_pdos_dxy_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#3(d)", spin, j=4
            )

            _, v_pdos_4s_part = data_loader(
                fname + ".pdos_atm#" + str(n) + "(V)_wfc#4(s)", spin
            )

            v_pdos_s.append(v_pdos_s_part)
            v_pdos_p.append(v_pdos_p_part)

            v_pdos_dz2.append(v_pdos_dz2_part)
            v_pdos_dzx.append(v_pdos_dzx_part)
            v_pdos_dzy.append(v_pdos_dzy_part)
            v_pdos_dx2y2.append(v_pdos_dx2y2_part)
            v_pdos_dxy.append(v_pdos_dxy_part)

            v_pdos_4s.append(v_pdos_4s_part)

    energy = v_energy - E_fermi
    pdos = [
        np.array(o_pdos_s),
        np.array(o_pdos_p),
        np.array(ge_pdos_d),
        np.array(ge_pdos_s),
        np.array(ge_pdos_p),
        np.array(v_pdos_s),
        np.array(v_pdos_p),
        np.array(v_pdos_dz2),
        np.array(v_pdos_dzx),
        np.array(v_pdos_dzy),
        np.array(v_pdos_dx2y2),
        np.array(v_pdos_dxy),
        np.array(v_pdos_4s),
    ]
    return energy, pdos_tot, pdos


# %%
def vo2_dict(fname, n_cells, E_fermi, spin=False, sd_reversal=False):
    pdos_names = [
        "o_pdos_s",
        "o_pdos_p",
        "v_pdos_s",
        "v_pdos_p",
        "v_pdos_dz2",
        "v_pdos_dzx",
        "v_pdos_dzy",
        "v_pdos_dx2y2",
        "v_pdos_dxy",
        "v_pdos_4s",
    ]
    extra_pdos_names = ["pdos_tot", "e_g", "t_2g", "v_pdos_d"]

    energy, pdos_tot, pdos_vals = vo2_pdos(fname, n_cells, E_fermi, spin, sd_reversal)

    pdos = dict(zip(pdos_names, pdos_vals))

    e_g = pdos["v_pdos_dz2"] + pdos["v_pdos_dxy"]
    t_2g = pdos["v_pdos_dx2y2"] + pdos["v_pdos_dzx"] + pdos["v_pdos_dzy"]
    v_pdos_d = (
        pdos["v_pdos_dz2"]
        + pdos["v_pdos_dxy"]
        + pdos["v_pdos_dx2y2"]
        + pdos["v_pdos_dzx"]
        + pdos["v_pdos_dzy"]
    )

    extra_pdos_vals = [pdos_tot, e_g, t_2g, v_pdos_d]

    extra_pdos = dict(zip(extra_pdos_names, extra_pdos_vals))

    all_pdos = dict(pdos, **extra_pdos)
    return energy, all_pdos


# %%
def vo2_dict2(fname, n_cells, E_fermi, spin=False, sd_reversal=False):
    pdos_names = [
        "o_pdos_s",
        "o_pdos_p",
        "v_pdos_s",
        "v_pdos_p",
        "v_pdos_dz2",
        "v_pdos_dzx",
        "v_pdos_dzy",
        "v_pdos_dx2y2",
        "v_pdos_dxy",
        "v_pdos_4s",
    ]
    extra_pdos_names = ["e_g", "t_2g", "v_pdos_d"]

    energy, pdos_tot, pdos_vals = vo2_pdos2(fname, n_cells, E_fermi, spin, sd_reversal)

    pdos = dict(zip(pdos_names, pdos_vals))

    e_g = pdos["v_pdos_dz2"] + pdos["v_pdos_dxy"]
    t_2g = pdos["v_pdos_dx2y2"] + pdos["v_pdos_dzx"] + pdos["v_pdos_dzy"]
    v_pdos_d = (
        pdos["v_pdos_dz2"]
        + pdos["v_pdos_dxy"]
        + pdos["v_pdos_dx2y2"]
        + pdos["v_pdos_dzx"]
        + pdos["v_pdos_dzy"]
    )

    extra_pdos_vals = [e_g, t_2g, v_pdos_d]

    extra_pdos = dict(zip(extra_pdos_names, extra_pdos_vals))

    all_pdos = dict(pdos, **extra_pdos)
    return energy, all_pdos, pdos_tot


# %%
def vo2_dict3(fname, outname, n_cells, E_fermi, spin=False, sd_reversal=False):
    pdos_names = [
        "o_pdos_s",
        "o_pdos_p",
        "ge_pdos_d",
        "ge_pdos_s",
        "ge_pdos_p",
        "v_pdos_s",
        "v_pdos_p",
        "v_pdos_dz2",
        "v_pdos_dzx",
        "v_pdos_dzy",
        "v_pdos_dx2y2",
        "v_pdos_dxy",
        "v_pdos_4s",
    ]
    extra_pdos_names = ["e_g", "t_2g", "v_pdos_d"]

    energy, pdos_tot, pdos_vals = vo2_pdos3(
        fname, outname, n_cells, E_fermi, spin, sd_reversal
    )

    pdos = dict(zip(pdos_names, pdos_vals))

    e_g = pdos["v_pdos_dz2"] + pdos["v_pdos_dxy"]
    t_2g = pdos["v_pdos_dx2y2"] + pdos["v_pdos_dzx"] + pdos["v_pdos_dzy"]
    v_pdos_d = (
        pdos["v_pdos_dz2"]
        + pdos["v_pdos_dxy"]
        + pdos["v_pdos_dx2y2"]
        + pdos["v_pdos_dzx"]
        + pdos["v_pdos_dzy"]
    )

    extra_pdos_vals = [e_g, t_2g, v_pdos_d]

    extra_pdos = dict(zip(extra_pdos_names, extra_pdos_vals))

    all_pdos = dict(pdos, **extra_pdos)
    return energy, all_pdos, pdos_tot


# %%
def parse_filband(feig, npl=10):
    # feig : filband in bands.x input file
    # npl : number per line, 10 for bands.x, 6 for phonon

    f = open(feig, "r")
    lines = f.readlines()

    header = lines[0].strip()
    line = header.strip("\n")
    shape = re.split("[,=/]", line)
    nbnd = int(shape[1])
    nks = int(shape[3])
    eig = np.zeros((nks, nbnd + 1), dtype=np.float32)

    dividend = nbnd
    divisor = npl
    div = nbnd // npl + 1 if nbnd % npl == 0 else nbnd // npl + 2
    kinfo = []
    for index, value in enumerate(lines[1:]):
        value = value.strip(" \n")
        quotient = index // div
        remainder = index % div

        if remainder == 0:
            kinfo.append(value)
        else:
            value = re.split("[ ]+", value)
            a = (remainder - 1) * npl
            b = a + len(value)
            eig[quotient][a:b] = value

    f.close()

    return eig, nbnd, nks, kinfo


# %%
def bands_from_gnu(datafile, scale=1):
    zb = np.loadtxt(datafile)
    xb = np.unique(zb[:, 0])
    bands = []
    bndl = len(zb[zb[:, 0] == xb[1]])

    for i in range(0, bndl):
        bands.append(np.zeros([len(xb), 2]))
    for i in range(0, len(xb)):
        sel = zb[zb[:, 0] == xb[i]]
        for j in range(0, bndl):
            bands[j][i][0] = xb[i]
            bands[j][i][1] = np.multiply(sel[j][1], scale)

    return bands


# %%


# %%
def read_cohp(fname, headlen, toplen=10, bondlen=False):
    cohp = np.loadtxt(fname, skiprows=headlen)
    cohp_abs = abs(cohp)
    cohp_sum = cohp_abs.sum(axis=0)

    cohp_sum = np.delete(cohp_sum, 0)
    cohp_sum = np.delete(cohp_sum, 0)
    cohp_sum = np.delete(cohp_sum, 0)
    cohp_sum = np.delete(cohp_sum, np.arange(1, cohp_sum.size, 2))

    topids = (-cohp_sum).argsort()[:toplen]
    topids = [i for i in topids]

    cohp_labels = np.genfromtxt(fname, dtype="str", skip_header=3, skip_footer=801)
    if not bondlen:
        cohp_labels = [
            label.replace("(" + label.split("(")[1], "") for label in cohp_labels
        ]
    cohp_labels = [
        label.replace(label.split(":")[0] + ":", "") for label in cohp_labels
    ]
    cohpT = [row for row in cohp.T]
    cohpT_en = cohpT[0]
    cohpT_ave = cohpT[1]
    cohpT = cohpT[3::2]
    return cohp_sum, topids, cohp_labels, cohpT_en, cohpT_ave, cohpT


# %%
def read_doe(fname, headlen=6):
    doe = np.loadtxt(fname, skiprows=headlen)
    doe_energy = doe[:, 0]
    doe_val = doe[:, 1]
    doe_intval = doe[:, 2]
    return doe_energy, doe_val, doe_intval


# %%


# %%
def bondlen(fname, vvlim=3.2, volim=2.5):
    scell = np.loadtxt(fname, unpack=True)
    scellT = [row[:4] for row in scell.T]
    vs_scell = [v for v in scellT if v[0] == 23.0]
    os_scell = [o for o in scellT if o[0] == 8.0]
    o_scell_atoms = np.array([j[1:] for j in os_scell[:]])
    v_scell_atoms = np.array([j[1:] for j in vs_scell[:]])

    scell_norms = []

    for pair in itertools.combinations(v_scell_atoms, 2):
        pos1 = pair[0]
        pos2 = pair[1]

        dist = np.linalg.norm(pos1 - pos2)
        scell_norms.append(dist)

    scell_norms.sort()
    scell_shortest = [x for x in scell_norms if x <= vvlim]

    scell_norms_vo = []

    for pair in itertools.product(v_scell_atoms, o_scell_atoms):
        pos1 = pair[0]
        pos2 = pair[1]

        dist = np.linalg.norm(pos1 - pos2)
        scell_norms_vo.append(dist)

    scell_norms_vo.sort()
    scell_vo_shortest = [x for x in scell_norms_vo if x <= volim]

    return scell_shortest, scell_vo_shortest


# %%
def bondlen_fromfile(
    fname, use_fname=True, vvlim=3.2, volim=2.5, vvhalf=2.78, BtoA_trans=True
):
    if BtoA_trans:
        BtoA = 0.529177
    else:
        BtoA = 1

    if use_fname:
        struc = ParseQEO()
        struc.Process(fname)
    else:
        pass

    vs = [v for v in struc.atoms if v[0] == "V"]
    v_atoms = np.array([j[1:] for j in vs[:]])

    os = [o for o in struc.atoms if o[0] == "O"]
    o_atoms = np.array([j[1:] for j in os[:]])

    lat = np.array(list(struc.lattice.values()))

    v_scell_atoms = np.array([np.dot(lat.T, atom) for atom in v_atoms])
    o_scell_atoms = np.array([np.dot(lat.T, atom) for atom in o_atoms])

    scell_norms = []

    for pair in itertools.combinations(v_scell_atoms, 2):
        pos1 = pair[0]
        pos2 = pair[1]

        dist = np.linalg.norm(pos1 - pos2)
        scell_norms.append(dist)

    scell_norms.sort()
    scell_shortest = [x * BtoA for x in scell_norms if x * BtoA <= vvlim]

    scell_norms_vo = []

    for pair in itertools.product(v_scell_atoms, o_scell_atoms):
        pos1 = pair[0]
        pos2 = pair[1]

        dist = np.linalg.norm(pos1 - pos2)
        scell_norms_vo.append(dist)

    scell_norms_vo.sort()
    scell_vo_shortest = [x * BtoA for x in scell_norms_vo if x * BtoA <= volim]

    return scell_shortest, scell_vo_shortest


# %%
def sscell(fname, only_vv=True):
    struc = ParseQEO()
    struc.Process(fname)

    vs = [v for v in struc.atoms if v[0] == "V"]
    v_atoms = np.array([j[1:] for j in vs[:]])

    lat = np.array(list(struc.lattice.values()))

    atoms100 = v_atoms + np.array([1, 0, 0])
    atoms100 = np.append(v_atoms, atoms100, axis=0)

    atoms010 = v_atoms + np.array([0, 1, 0])
    atoms010 = np.append(atoms100, atoms010, axis=0)

    atoms001 = v_atoms + np.array([0, 0, 1])
    atoms001 = np.append(atoms010, atoms001, axis=0)

    atoms110 = v_atoms + np.array([1, 1, 0])
    atoms110 = np.append(atoms001, atoms110, axis=0)

    atoms101 = v_atoms + np.array([1, 0, 1])
    atoms101 = np.append(atoms110, atoms101, axis=0)

    atoms011 = v_atoms + np.array([0, 1, 1])
    atoms011 = np.append(atoms101, atoms011, axis=0)

    atoms111 = v_atoms + np.array([1, 1, 1])
    atoms111 = np.append(atoms011, atoms111, axis=0)

    v_sscell_atoms = np.array([np.dot(lat.T, atom) for atom in atoms111])

    return v_sscell_atoms


# %%
def bondlen_from_sscell(sscell, vvlim=3.2, volim=2.5, vvhalf=2.78, BtoA_trans=True):
    if BtoA_trans:
        BtoA = 0.529177
    else:
        BtoA = 1

    scell_norms = []

    for pair in itertools.combinations(sscell, 2):
        pos1 = pair[0]
        pos2 = pair[1]

        dist = np.linalg.norm(pos1 - pos2)
        scell_norms.append(dist)

    scell_norms.sort()
    scell_shortest = [x * BtoA for x in scell_norms if x * BtoA <= vvlim]

    return scell_shortest


# %%
def angle_from_sscell(sscell, normal=np.array([1, 0, 0]), vvhalf=2.78, BtoA_trans=True):
    if BtoA_trans:
        BtoA = 0.529177
    else:
        BtoA = 1

    angles = []

    for pair in itertools.combinations(sscell, 2):
        pos1 = pair[0]
        pos2 = pair[1]

        vec = pos1 - pos2
        dist = np.linalg.norm(vec)

        if dist * BtoA <= vvhalf:
            c = (
                np.dot(vec, normal) / np.linalg.norm(vec) / np.linalg.norm(normal)
            )  # -> cosine of the angle
            angle = np.arccos(np.clip(c, -1, 1))

            if angle > np.pi / 2:
                angle = np.pi - angle

            angles.append(angle * 180 / np.pi)
    return angles


# %%


# %%
def percentile(a, b, percentile):
    return (1 - percentile) * a + percentile * b


# %%
def dz2_combo(dxy, dz2):
    return 0.5 * dxy + np.sqrt(3) * 0.5 * dz2


# %%
def get_indices(outfname):
    struc = ParseQEO()
    struc.Process(outfname)

    atoms_names = [j[0] for j in struc.atoms[:]]

    O_ind = [i for i, x in enumerate(atoms_names) if x == "O"]
    O_inds = [x + 1 for x in O_ind]
    V_ind = [i for i, x in enumerate(atoms_names) if x == "V"]
    V_inds = [x + 1 for x in V_ind]
    Ge_ind = [i for i, x in enumerate(atoms_names) if x == "Ge"]
    Ge_inds = [x + 1 for x in Ge_ind]

    return O_inds, V_inds, Ge_inds


# %%
import matplotlib.pyplot as plt

datafile1 = "../../dmft/data/m2/test3/v1/bands.m2.dat.gnu"

bands1 = bands_from_gnu(datafile1)

kscale = 1  # max(bands[0][:,0])/max(ea)

fig, ax = plt.subplots(figsize=(10, 8))

for i in bands1:  # Here we plots the bands
    ax.plot(i[:, 0], i[:, 1] - 10, color="black", linewidth=2)

# for i in bands2: #Here we plots the bands
#     k_rat=max(bands1[0][:,0])/max(bands2[0][:,0])
#     ax.plot(k_rat*i[:,0],i[:,1]-Ef2,color="red", linewidth=2)

# for i in bands3: #Here we plots the bands
#     k_rat=max(bands1[0][:,0])/max(bands3[0][:,0])
#     ax.plot(k_rat*i[:,0],i[:,1]-Ef3,color="blue", linewidth=2)


vlines = np.linspace(0, 1, 12) * max(bands1[0][:, 0])

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
ax.set_xticks(vlines, xlabeltext, fontsize=14)
ax.axhline(0, linestyle="--", color="k")

# ax.scatter(ea*kscale, bandsa, marker='.', c='r', alpha=.2)

ax.plot([], [], color="black", linewidth=2, label="exp M2")
# ax.plot([],[],color="red", linewidth=2, label='rlx M2')
# ax.plot([],[],color="blue", linewidth=2, label='vc-rlx M2')
ax.plot([], [], color="green", linewidth=2, label="symm M2")

ax.legend(fontsize=18, frameon=False, loc=4)
ax.tick_params(axis="both", which="major", labelsize=17)

ax.set_ylabel(r"E - E$_F$ (eV)", fontsize=18, labelpad=-10)

ax.set_xlim(0, max(bands1[0][:, 0]))
ax.set_ylim(-1, 2.5)

plt.plot()
# %%
