import os
from cffi import FFI
import numpy as np
import scipy.linalg as sl
from ase.units import Ha, Bohr, Ry
from ase.io import read
from ase import Atoms
from TB2J.utils import kmesh_to_R, symbol_number
from TB2J.myTB import  AbstractTB
import matplotlib.pyplot as plt
from TB2J_OpenMX.cmod._scfout_parser import ffi, lib
import copy

## Create the dictionary mapping ctypes to np dtypes.
ctype2dtype = {'int': 'i4', 'double': 'f8'}

## Integer types
for prefix in ('int', 'uint'):
    for log_bytes in range(4):
        ctype = '%s%d_t' % (prefix, 8 * (2**log_bytes))
        dtype = '%s%d' % (prefix[0], 2**log_bytes)
        ctype2dtype[ctype] = np.dtype(dtype)

## Floating point types
ctype2dtype['float'] = np.dtype('f4')
ctype2dtype['double'] = np.dtype('f8')


def asarray(ffi, ptr, length):
    ## Get the canonical C type of the elements of ptr as a string.
    T = ffi.getctype(ffi.typeof(ptr).item)

    if T not in ctype2dtype:
        raise RuntimeError("Cannot create an array for element type: %s" % T)

    return np.frombuffer(ffi.buffer(ptr, length * ffi.sizeof(T)),
                         ctype2dtype[T])


class OpenmxWrapper(AbstractTB):
    def __init__(self, path, prefix='openmx'):
        self.is_siesta = False
        self.is_orthogonal = False
        xyz_fname = os.path.join(path, prefix + '.xyz')
        fname = os.path.join(path, prefix + '.scfout')
        self.fname = fname
        self.R2kfactor = 2.0j * np.pi
        self.parse_scfoutput()
        self.Rdict = dict()
        for i, R in enumerate(self.R):
            self.Rdict[tuple(R)] = i
        atoms = read(xyz_fname)
        self.atoms = Atoms(atoms.get_chemical_symbols(),
                           cell=self.cell,
                           positions=self.positions)
        self.norbs_to_basis(self.atoms, self.norbs)
        self.nspin = 2
        self.nbasis = self.nspin * self.norb
        self._name='OpenMX'

    def solve(self, k):
        phase = np.exp(self.R2kfactor * (self.R @ k))
        Hk = np.einsum('rij, r->ij', self.H, phase)
        Sk = np.einsum('rij, r->ij', self.S, phase)
        return sl.eigh(Hk, Sk)

    def solve_all(self, kpts):
        nk = len(kpts)
        evals = np.zeros((nk, self.norb))
        evecs = np.zeros((nk, self.norb, self.norb), dtype=complex)
        for ik, k in enumerate(kpts):
            evals[ik], evecs[ik] = self.solve(k)
        return evals, evecs

    def HSE_k(self, k, convention=2):
        phase = np.exp(self.R2kfactor * (self.R @ k))
        Hk = np.einsum('rij, r->ij', self.H, phase)
        Sk = np.einsum('rij, r->ij', self.S, phase)
        evalue, evec = sl.eigh(Hk, Sk)
        return Hk, Sk, evalue, evec


    def HS_and_eigen(self, kpts):
        nk = len(kpts)
        evals = np.zeros((nk, self.nbasis))
        evecs = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        Hk = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        Sk = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        for ik, k in enumerate(kpts):
            phase = np.exp(self.R2kfactor * (self.R @ k))
            Hk[ik] = np.einsum('rij, r->ij', self.H, phase)
            Sk[ik] = np.einsum('rij, r->ij', self.S, phase)
            evals[ik], evecs[ik] = sl.eigh(Hk[ik], Sk[ik])
        return Hk, Sk, evals, evecs

    def get_hamR(self, R):
        return self.H[self.Rdict[tuple(R)]]

    def norbs_to_basis(self, atoms, norbs):
        self.basis = []
        symbols = atoms.get_chemical_symbols()
        sn = list(symbol_number(symbols).keys())
        for i, n in enumerate(norbs):
            for x in range(n):
                self.basis .append((sn[i], f'orb{x+1}', 'up'))
                self.basis .append((sn[i], f'orb{x+1}', 'down'))
        return self.basis

    def parse_scfoutput(self):
        argv0 = ffi.new("char[]", b"")
        argv = ffi.new("char[]", bytes(self.fname, encoding='ascii'))
        lib.read_scfout([argv0, argv])
        lib.prepare_HSR()
        self.ncell = lib.TCpyCell + 1
        self.natom = lib.atomnum
        self.norbs = np.copy(
            asarray(ffi, lib.Total_NumOrbs, self.natom + 1)[1:])

        if lib.SpinP_switch == 3:
            self.non_collinear = True
        elif lib.SpinP_switch == 1:
            self.non_collinear = False
        else:
            raise ValueError(
                " The value of SpinP_switch is %s. Can only get J from collinear and non-collinear mode."%lib.SpinP_switch)

        fnan = asarray(ffi, lib.FNAN, self.natom + 1)

        natn = []
        for iatom in range(self.natom):
            natn.append(asarray(ffi, lib.natn[iatom + 1], fnan[iatom + 1] + 1))

        ncn = []
        for iatom in range(self.natom):
            ncn.append(asarray(ffi, lib.ncn[iatom + 1], fnan[iatom + 1] + 1))
        # atv
        #  x,y,and z-components of translation vector of
        #periodically copied cells
        #size: atv[TCpyCell+1][4];
        atv = []
        for icell in range(self.ncell):
            atv.append(asarray(ffi, lib.atv[icell], 4))
        atv = np.copy(np.array(atv))

        atv_ijk = []
        for icell in range(self.ncell):
            atv_ijk.append(asarray(ffi, lib.atv_ijk[icell], 4))
        atv_ijk = np.array(atv_ijk)
        self.R = atv_ijk[:, 1:]
        tv = []
        for i in range(4):
            tv.append(asarray(ffi, lib.tv[i], 4))
        tv = np.array(tv)
        self.cell = np.copy(tv[1:, 1:]) * Bohr

        rtv = []
        for i in range(4):
            rtv.append(asarray(ffi, lib.rtv[i], 4))
        rtv = np.array(rtv)
        self.rcell = np.copy(rtv[1:, 1:])

        Gxyz = []
        for iatom in range(self.natom):
            Gxyz.append(asarray(ffi, lib.Gxyz[iatom + 1], 60))
        self.positions = np.copy(np.array(Gxyz)[:, 1:4]) * Bohr

        self.MP = np.copy(asarray(ffi, lib.MP, self.natom + 1)[1:])

        self.norb = lib.T_NumOrbs
        norb = self.norb

        if self.non_collinear:
            HR = np.zeros([self.ncell, 4, lib.T_NumOrbs, lib.T_NumOrbs])
            for iR in range(0, self.ncell):
                for ispin in range(lib.SpinP_switch + 1):
                    for iorb in range(lib.T_NumOrbs):
                        HR[iR, ispin,
                           iorb, :] = asarray(ffi, lib.HR[iR][ispin][iorb],
                                              norb)

            HR_imag = np.zeros([self.ncell, 3, lib.T_NumOrbs, lib.T_NumOrbs])
            for iR in range(0, self.ncell):
                for ispin in range(3):
                    for iorb in range(lib.T_NumOrbs):
                        HR_imag[iR, ispin, iorb, :] = asarray(
                            ffi, lib.HR_imag[iR][ispin][iorb], norb)

            self.H = np.zeros(
                [self.ncell, lib.T_NumOrbs * 2, lib.T_NumOrbs * 2],
                dtype=complex)

            # up up
            for iR in range(self.ncell):
                self.H[iR, ::2, ::2] = HR[iR, 0, :, :] + 1j * HR_imag[iR, 0, :, :]
                # up down
                self.H[iR, ::2,
                   1::2] = HR[iR, 2, :, :] + 1j * (HR[iR, 3, :, :] +
                                                   HR_imag[iR, 2, :, :])
                # down up
                self.H[iR,
                   1::2, ::2] = HR[iR, 2, :, :] - 1j * (HR[iR, 3, :, :] +
                                                          HR_imag[iR, 2, :, :])
                # down down
                self.H[iR, 1::2, 1::2] = HR[iR, 1, :, :] + 1j * HR_imag[iR, 1, :, :]
        else:  # collinear

            HR = np.zeros([self.ncell, 4, lib.T_NumOrbs, lib.T_NumOrbs])
            for iR in range(0, self.ncell):
                for ispin in range(lib.SpinP_switch + 1):
                    for iorb in range(lib.T_NumOrbs):
                        HR[iR, ispin,
                           iorb, :] = asarray(ffi, lib.HR[iR][ispin][iorb],
                                              norb)

            self.H = np.zeros(
                [self.ncell, lib.T_NumOrbs * 2, lib.T_NumOrbs * 2],
                dtype=complex)

            # up up
            for iR in range(self.ncell):
                self.H[iR, ::2, ::2] = HR[iR, 0, :, :]
                self.H[iR, 1::2, 1::2] = HR[iR, 1, :, :]
        self.efermi = lib.ChemP * Ha
        self.H *= Ha

        SR = np.zeros([self.ncell, lib.T_NumOrbs, lib.T_NumOrbs])
        for iR in range(0, self.ncell):
            for iorb in range(lib.T_NumOrbs):
                SR[iR, iorb, :] = asarray(ffi, lib.SR[iR][iorb], norb)
        self.S = np.kron( SR, np.eye(2))
        lib.free_HSR()
        lib.free_scfout()
        print("Loading from scfout file OK!")

def test():
    openmx = OpenmxWrapper(
        path='/home/hexu/projects/TB2J_example/OPENMX/SrMnO3_FM_SOC/')
    #hsr = openmx.parse_scfoutput()
    #from banddownfolder.plot import plot_band
    #plot_band(hsr)
    #plt.savefig('band.pdf')
    #plt.show()


#test()
