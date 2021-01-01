import numpy as np
from math import sqrt, pi, log, ceil

class Checkpt:
    def __init__(self):
        self.n_atoms = 0
        self.charge = 0.0
        self.n_elec = 0
        self.n_basis = 0
        self.coords = np.array([])
        self.nuclei = np.array([])
        self.shellTypes = np.array([])
        self.primitives = np.array([])
        self.atomMap = np.array([])
        self.exp = np.array([])
        self.coeff = np.array([])
        self.spCoeff = np.array([])
        self.shellCoords = np.array([])
        self.overlap = np.array([])
        self.MOcoeff = np.array([])
        self.densityMatrix = np.array([])
        self.coreHam = np.array([])
        pass

    def read_real(self, lines, i, perline=5):
        split = lines[i].split()
        n_elms = int(split[-1])
        n_lines = ceil(n_elms/perline)
        sData = lines[i+1:i+n_lines+1]
        data = np.zeros(n_elms)
        for i in range(n_lines):
            stuff = sData[i].split()
            data[i*perline:i*perline+len(stuff)] = np.array(stuff, dtype=float)
        return data

    def flat_to_matrix(self, data):
        dim = len(data)
        dim = int((sqrt(dim*8 + 1) - 1)/2)
        pMat = np.zeros((dim, dim))
        pMat[np.tril_indices(dim)] = data
        pMat += pMat.transpose()
        pMat[np.diag_indices(dim)] = pMat[np.diag_indices(dim)]*0.5
        return pMat

    def import_file(self, fileLoc):
        with open(fileLoc, 'r') as file:
            lines = file.readlines()
            for i in range(len(lines)):
                line = lines[i]
                split = line.split()
                if "Number of atoms" in line:
                    self.n_atoms = int(split[-1])
                elif "Charge" in line:
                    self.charge = float(split[-1])
                elif "Number of electrons" in line:
                    self.n_elec = int(split[-1])
                elif "Current cartesian coordinates" in line:
                    self.coords = self.read_real(lines, i)
                elif "Nuclear charges" in line:
                    self.nuclei = self.read_real(lines, i)
                elif "Number of basis functions" in line:
                    self.n_basis = int(split[-1])
                elif "Shell types" in line:
                    self.shellTypes = self.read_real(lines, i, perline=6).astype(int)
                elif "Number of primitives per shell" in line:
                    self.primitives = self.read_real(lines, i, perline=6).astype(int)
                elif "Shell to atom map" in line:
                    self.atomMap = self.read_real(lines, i, perline=6).astype(int)
                elif "Primitive exponents" in line:
                    self.exp = self.read_real(lines, i)
                elif "P(S=P) Contraction coefficients"in line:
                    self.spCoeff = self.read_real(lines, i)
                elif "Contraction coefficients" in line:
                    self.coeff = self.read_real(lines, i)
                elif "Coordinates of each shell" in line:
                    self.shellCoords = self.read_real(lines, i)
                elif "Overlap Matrix" in line:
                    self.overlap = self.read_real(lines, i)
                elif "Alpha MO coefficients" in line:
                    self.MOcoeff = self.read_real(lines, i)
                elif "Total SCF Density" in line:
                    self.densityMatrix = self.read_real(lines, i)
                elif "Core Hamiltonian Matrix" in line:
                    self.coreHam = self.read_real(lines, i)