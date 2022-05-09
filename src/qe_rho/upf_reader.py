from xml.etree import ElementTree as ET
from os.path import join
from os import listdir
import numpy as np
from .utils import simpson_ndarray as simpson


class UPF(dict):
    '''
    returns a dict of dict
    '''
    def __init__(self, upf_file_path:str):
        rcut = 10e0

        self.upf_file_path = upf_file_path
        try:
            tree = ET.parse(upf_file_path)
        except ET.ParseError: # upf file without the real root upf
            with open(upf_file_path, 'r') as upf_original:
                tree = ET.ElementTree(ET.fromstring('<UPF>\n' + upf_original.read() +'\n</UPF>'))

        root = tree.getroot()
        self.update(self.root2dict(root))
        self.rho_atom = self.str2array(self['UPF']['PP_RHOATOM']['body'])
        self.pp_r = self.str2array(self['UPF']['PP_MESH']['PP_R']['body'])
        self.pp_rab = self.str2array(self['UPF']['PP_MESH']['PP_RAB']['body'])
        #self.header = self.str2array(self['UPF']['PP_HEADER']['body'])
        
        self.mesh = len(self.rho_atom)
        for j in range(self.mesh):
            if self.pp_r[self.mesh-1-j] < rcut:
                break
        self.msh = (self.mesh - j)//2 * 2 - 1 # msh be odd
        #self.zp = simpson(self.mesh, self.rho_atom, self.pp_rab)
        if 'z_valence' in self['UPF']['PP_HEADER']['tag'].keys():
            self.zp = float(self['UPF']['PP_HEADER']['tag']['z_valence'])
        elif 'Z_valence' in self['UPF']['PP_HEADER']['tag'].keys():
            self.zp = float(self['UPF']['PP_HEADER']['tag']['Z_valence'])
        else:
            if 'valence' in self['UPF']['PP_HEADER']['body']:
                for i in self['UPF']['PP_HEADER']['body'].splitlines():
                    if 'Z valence' in i:
                        self.zp = float(i.replace('Z valence', '').strip())
                    elif 'z valence' in i:
                        self.zp = float(i.replace('z valence', '').strip())


        self.rho_scale = self.zp/(simpson(self.msh, self.rho_atom[:self.msh], self.pp_rab[:self.msh]))
        # the radial grid is defined up to r(mesh) but the auxiliary variable msh to limit the grid
        # is rcut = 10 a.u. to cutoff the numerical nose arising from the large-r tail
        # incases li

    def root2dict(self, root):
        temp_dict = {}
        temp_dict.update({root.tag: {'tag':root.attrib}})
        temp_dict[root.tag].update({'body':root.text})
        for item in root:
            #temp_dict.update({item.tag: self.root2dict(item)})
            temp_dict[root.tag].update(self.root2dict(item))
        return temp_dict

    def str2array(self, body:str):
        # first seperate the str to list by space, remove line break
        temp_list = [i.replace('\n', '') for i in body.split(' ')]
        # remove the empty items
        return np.array([float(i) for i in temp_list if i])
                    

class pseudo_upf():
    '''
    q-e/upflib/pseudo_types.f90
    '''
    def __init__(self):
        self.generated = '' # generator software
        self.author = '' # pseudopotential's author
        self.date = '' # generation date
        self.comment = '' # author's comment
        self.psd = '' # element label
        self.typ = '' # pseudo type: NC, SL, US, PAW, 1/r
    
class paw_in_upf():
    def __init__(self):
        self.ae_rho_atc = [] #  AE core charge (pseudo ccharge is already included in upf)
        self.pfunc = [] # Psi_i(r)*Psi_j(r)
        self.pfunc_rel = [] # Psi_i(r) * Psi_j(r) small component
        self.ftfunc = [] # as above, but for pseudo
        self.aewfc_rel = [] # as above, but for pseudo

        self.ae_vloc = [] # AE local potential
        self.oc = [] # starting occupation used to init

if __name__ == '__main__':
    UPF_FOLDER = 'UPFs/PSEUDOPOTENTIALS/'
    for i in listdir(UPF_FOLDER): 
        if i.startswith('Ti'):
            a = UPF(join(UPF_FOLDER, i))
            print(i, a.zp, simpson(len(a.pp_rab), a.rho_atom, a.pp_rab))
            print(a.rho_scale)
