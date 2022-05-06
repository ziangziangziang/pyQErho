# -*- coding: utf-8 -*-
from os import listdir, mkdir
from os.path import exists, join
import numpy as np
from utils import simpson_ndarray as simpson
from utils import realspace_grid_init, apostrophe_remover
from utils import constant_book, qe_cell, constant_book, qe_float
from utils import struc_fact, pwscf_parser
from upf_reader import UPF
from scipy.io import FortranFile
import time
from meshplot import mesh_plot

constants = constant_book()


def atomic_rho_g(ngm, ngl, gl, ntyp, upfs, cell, strf, nspina=1):
    '''
    adapted from PW/src/atomic_rho.f90
    ngm: int, number of g vectors
    ngl: int, number of g shells (length of g vectors)
    gl: array, length of g vectors
    upfs: dict, pseudopotentials
    ntyp: int, number of atom types
    strf: structure factor

    
    compute superposition of atomic charges in reciprocal space
    only support nspina = 1 in version 0.1
    '''
    eps8 = constants.eps8
    gl_uniq, gl_uniq_id = np.unique(gl, return_inverse=True)
    #print(gl_uniq_id.shape)
    #ind_sort = np.argsort(gl_uniq)
    #gl_uniq_sort = gl_uniq[ind_sort]
    #gl_uniq_id_sort = gl_uniq_id[ind_sort]
    

    rhocg = np.zeros(ngm, dtype=np.cdouble)
    tpiba = cell.tpiba
    msh = [i.msh for i in upfs]

    for nt in range(ntyp):
        gx = np.sqrt(gl_uniq) * tpiba
        factor = np.outer(gx, upfs[nt].pp_r[:msh[nt]])
        factor[0, :] += eps8
        factor[:, 0] += eps8
        ff = np.sin(factor)/(factor)
        aux = upfs[nt].rho_atom[:msh[nt]] * ff
        rhocgnt_uniq = simpson(msh[nt], aux, upfs[nt].pp_rab[:msh[nt]])*upfs[nt].rho_scale

        rhocgnt = rhocgnt_uniq[gl_uniq_id]

        rhocg += strf[:, nt] * rhocgnt/cell.omega
    
    return rhocg

   

class RhoInit():
    '''
    generate SAD charge density from pwscf.in file
    adapted from PW/src/run_pwscf.f90
    '''
    def __init__(self, pwscfin):
        self.pwscfin = pwscf_parser(pwscfin)
        self.nat = int(self.pwscfin['SYSTEM']['nat'])
        self.ntyp = int(self.pwscfin['SYSTEM']['ntyp'])
        self.atom_order = []
        self.tau = np.empty((self.nat, 3))
        self.type_atom_dict = {}        # {atom_name: [atom_order, ]}
        self.upf_files = {}              # {atom_name: upf}
        self.upfs = {}         # {atom_name: UPF class}
        self.upf_folder = apostrophe_remover(self.pwscfin['CONTROL']['pseudo_dir'])
        self.prefix = apostrophe_remover(self.pwscfin['CONTROL']['prefix'])
        self.save_folder = f'{self.prefix}.save'
        self.gamma_only = 0
        self.nspin = 1

        if not exists(self.save_folder):
            mkdir(self.save_folder)

        self.upf_order = {}
        self.upf_list = []
        for iele, ele in enumerate(self.pwscfin['ATOMIC_SPECIES']['text']):
            eleV = [i for i in ele.split(' ') if i]
            ele = eleV[0]
            V = eleV[2]
            self.upf_files[ele] = V
            self.upf_order[ele] = iele
            self.upf_list.append(ele)
            try:
                self.upfs[ele] = UPF(join(self.upf_folder, V))
            except FileNotFoundError:
                print('can not find upf file ' + join(self.upf_folder, V))
                quit()
            self.type_atom_dict[ele] = []
        self.ityp = [[] for i in self.type_atom_dict.keys()]
        cell_vector = np.array([[qe_float(j) for j in i.split(' ') if j] for i in self.pwscfin['CELL_PARAMETERS']['text']])
        self.cell = qe_cell(cell_vector)

        for iat, atom in enumerate(self.pwscfin['ATOMIC_POSITIONS']['text']):
            atpos = [i for i in atom.split(' ') if i]
            at = atpos[0]
            self.type_atom_dict[at].append(iat)         # record the order of the atom
            pos = np.array([qe_float(i) for i in atpos[1:4]])
            self.atom_order.append(at)
            if self.pwscfin['pos_angstrom']:
                self.tau[iat, :] = pos/self.cell.alat
            else:
                self.tau[iat, :] = np.matmul(pos, self.cell.at)
            self.ityp[self.upf_order[at]].append(iat)
            if at in self.type_atom_dict.keys():
                self.type_atom_dict[at].append(iat)
            else:
                self.type_atom_dict[at] = [iat]
        self.ecutwfc = qe_float(self.pwscfin['SYSTEM']['ecutwfc'])
        try:
            self.ecutrho = qe_float(self.pwscfin['SYSTEM']['ecutrho'])
        except KeyError:
            self.ecutrho = 4.0 * self.ecutwfc

        time_start = time.time()
        self.nr1, self.nr2, self.nr3, self.ngm, self.gvect, self.mill_g, self.ngl, self.gl = realspace_grid_init(self.cell.at, self.cell.bg, self.ecutrho/self.cell.tpiba2)

        time_start = time.time()
        self.strf = struc_fact(self.nat, self.tau, self.ntyp, self.ityp, self.nr1, self.nr2, self.nr3, self.cell.bg, self.ngm, self.gvect, self.mill_g)
        self.zp = np.array([self.upfs[self.upf_list[i]].zp*self.ityp[i].__len__() for i in range(self.ntyp)]).sum()

        time_start = time.time()
        self.rhog = atomic_rho_g(self.ngm, self.ngl, self.gl, self.ntyp, [self.upfs[i] for i in self.upf_list], self.cell, self.strf)

    def saverhog(self, chdens_path):
        rhog = np.empty((self.ngm,2), dtype=np.float)
        rhog[:, 0] = self.rhog.real
        rhog[:, 1] = self.rhog.imag
        rhog.reshape(2*self.ngm)

        f = FortranFile(chdens_path, 'w')
        f.write_record(np.array([self.gamma_only, self.ngm, self.nspin], dtype=np.int32))
        f.write_record(self.cell.bg*self.cell.tpiba)
        f.write_record((self.mill_g).reshape(self.ngm* 3))
        f.write_record(rhog)
        f.close()
    
    def rho_g2r(self):
        rhog = np.zeros((self.nr1, self.nr2, self.nr3), dtype=np.complex128)
        rhog[self.mill_g[:, 0], self.mill_g[:, 1], self.mill_g[:, 2]] = self.rhog
        rhor = np.fft.ifftn(rhog).real *self.cell.omega#/( self.nr1* self.nr2* self.nr3)
        # no need to devide the grid size
        return rhor

    def plot(self, output=None, filtering='none'):
        self.rhor = self.rho_g2r()
        mesh_plot(self.rhor, a=self.cell.at*self.cell.alat*constants.bohr_radius_angs, output_name=output, filtering=filtering)

    def read_charge(self, path):
        f = FortranFile(path, 'r')
        gamma_only, ngm_g, nspin = f.read_ints()

        bg = f.read_reals(np.float).reshape((3,3))
        miller_g = f.read_ints().reshape((ngm_g, 3))
        print(miller_g[:100])

        rhog_compact = f.read_reals(np.float).reshape([ngm_g, 2]).dot([1.0, 1.0j])
        f.close()
        rhog = np.zeros((self.nr1, self.nr2, self.nr3), dtype=np.complex128)
        rhog[miller_g[:, 0], miller_g[:, 1], miller_g[:, 2]] = rhog_compact
        rhor = np.fft.ifftn(rhog).real *self.cell.omega#/( self.nr1* self.nr2* self.nr3)
        return rhor


        


if __name__ == '__main__':
    start = time.time()
    pwscfin = RhoInit('pwscf.in')
    print(pwscfin.mill_g[:100])
    c = pwscfin.read_charge('chargedens-px.dat')
    pwscfin.saverhog(join(pwscfin.save_folder, 'charge-density.dat'))
 