import numpy as np

DEBUG = False

def struc_fact(nat, tau, ntyp, ityp, nr1, nr2, nr3, bg, ngm, g, mill):
    '''
    adapted from PW/src/struct_fact.f90
    encode the structure information
    variables:
        nat: number of atom in the unit cell
        tau(3, nat): positions of the atoms in the cell, crystal
        ntyp: type ofatoms
        ityp: for each atom gives the type (list)
        ngm: number of G vectors
        nr1, nr2, nr3: fft dimensions
        bg(3, 3): reciprocal crystal basis vectors
        g(3, ngm): coordiantes of the g vectors
        mill: miller index
    '''
    constants = constant_book()
    strf = np.empty((ngm, ntyp), dtype=np.cdouble)        # structure factor
    eigts1 = np.empty((2*nr1 + 1, nat), dtype=np.complex128) # phase in dimensiosn
    eigts2 = np.empty((2*nr2 + 1, nat), dtype=np.complex128)
    eigts3 = np.empty((2*nr3 + 1, nat), dtype=np.complex128)
    n1 = np.arange(2*nr1+1) - nr1
    n2 = np.arange(2*nr2+1) - nr2
    n3 = np.arange(2*nr3+1) - nr3

    bgtau = np.empty(3)                 # product of bg and tau

    for na in range(nat):
        bgtau = np.matmul(bg, tau[na])
        arg = constants.tpi * n1 * bgtau[0]
        eigts1[:, na] = np.exp(-1j * arg)
        arg = constants.tpi * n2 * bgtau[1]
        eigts2[:, na] = np.exp(-1j * arg)
        arg = constants.tpi * n3 * bgtau[2]
        eigts3[:, na] = np.exp(-1j * arg)

    for nt in range(ntyp):
        z = np.zeros(ngm, dtype=np.complex128)
        for na in range(nat):
            if na in ityp[nt]:
                z += eigts1[mill[:,0]+nr1, na] * eigts2[mill[:, 1]+nr2, na] * eigts3[mill[:, 2]+nr3, na]
        strf[:, nt] = z
    return strf
 
def pwscf_parser(file_path: str):
    '''
    read pwscf.in file in file_path
    '''
    pwscfin = {}
    with open(file_path, 'r') as f:
        pwscf_lines = f.readlines()

    for line in pwscf_lines:
        line = line.strip()
        if line.startswith('!') or line.startswith('/') or line=='':
            # skip comment line or namelist end mark
            pass
        else:
            # for namelists
            if line.startswith('&'):
                namelist = line[1:]
            elif line.startswith('ATOMIC_SPECIES'):
                namelist = 'ATOMIC_SPECIES'
            elif line.startswith('ATOMIC_POSITIONS'):
                namelist = 'ATOMIC_POSITIONS'
                if 'angstrom' in line:
                    pos_angstrom = True
                else:
                    pos_angstrom = False
            elif line.startswith('K_POINTS'):
                namelist = 'K_POINTS'
            elif line.startswith('CELL_PARAMETERS'):
                namelist = 'CELL_PARAMETERS'
            #for contents:
            if namelist in pwscfin.keys():
                pwscfin[namelist]['text'].append(line)
                if '=' in line:
                    key, value = line.split('=')
                    pwscfin[namelist][key.strip()]=value.strip('" ').strip("'")
            else:
                pwscfin[namelist] = {'text':[]}
    pwscfin['pos_angstrom'] = pos_angstrom

    return pwscfin


def simpson_ndarray(mesh:int, func:np.ndarray, rab:np.ndarray):
    '''
    adapted from q-e/upflib/simpsn.f90
    simpson rule integration. On input:
        mesh: the number of grid points ( should be odd_)
        func: function to be integrated
        rab:  dr
    output asum:
        \sum_i c_i f_i * rab_i = \int_0^\infty f(r) dr
        where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    '''
    funcr = func*rab
    first_term = funcr[..., :mesh-2:2] 
    second_term = 4*funcr[..., 1:mesh-1:2] 
    third_term = funcr[..., 2:mesh:2] 
    return (first_term+second_term+third_term).sum(axis=-1)/3


def simpson(mesh:int, func:np.ndarray, rab:np.ndarray):
    '''
    adapted from q-e/upflib/simpsn.f90
    simpson rule integration. On input:
        mesh: the number of grid points ( should be odd_)
        func: function to be integrated
        rab:  dr
    output asum:
        \sum_i c_i f_i * rab_i = \int_0^\infty f(r) dr
        where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    '''
    asum = 0
    r12 = 1/3
    f3 = func[0] * rab[0] * r12
    for i in range(1, mesh -2, 2):
        f1 = f3
        f2 = func[i] * rab[i] * r12
        f3 = func[i+1] * rab[i+1] * r12
        asum += f1 + 4*f2 + f3
    return asum

def hpsort_eps(ra,eps=1e-8):
    '''
    adapted from Modules/sort.f90
    to sort the miller indices by length of g-vector
    have to sort miller indices so QE can accept the charge density
    input: ra, 1-D array
    '''
    n = ra.shape[0]
    l = n//2 + 1
    ir = n
    ind = np.arange(n)
    while True:
        if l > 1:
            l = l-1
            rra = ra[l]
            iind = ind[l]

        else:
            # clear a s
            pass


class qe_cell():
    '''
    Modules/cell_base.f90
    '''
    def __init__(self, at):
        constants = constant_book()
        alat = np.sqrt(np.square(at[0]).sum()) 
        self.alat = alat / constants.bohr_radius_angs # alat in bohr
        self.at = at/alat
        self.tpiba = 2*constants.pi/self.alat
        self.tpiba2 = self.tpiba ** 2
        self.omega = cell_volume(self.at, self.alat)
        self.bg = cell_recips(self.at)

def cell_volume(a, alat=1):
    #omega = alat**3 * np.dot(a[0], np.cross(a[1], a[2]))
    omega = alat**3 * np.linalg.det(a)
    return abs(omega)

def cell_recips(at):
    '''
    Modules/recips.f90
    returns b in units of 2 pi / a
    '''
    den = np.linalg.det(at)
    b1 = np.cross(at[1], at[2])/den
    b2 = np.cross(at[2], at[0])/den
    b3 = np.cross(at[0], at[1])/den
    return np.array((b1, b2, b3))

class constant_book():
    '''
    Modules/constants.f90
    '''
    def __init__(self):
        # mathematical constants
        self.pi = np.float(3.14159265358979323846)
        self.tpi = self.pi * np.float(2)
        self.fpi = self.pi * np.float(4)
        self.sqrtpi = np.float(1.77245385090551602729)
        self.sqrtpm1 = np.float(1)/self.sqrtpi
        self.sqrt2 = np.float(1.41421356237309504880)

        # physical constants
        # physics.nist.gov/constants
        self.H_PLANCK_SI      = np.float(6.62606896e-34)   # J s
        self.K_BOLTZMANN_SI   = np.float(1.3806504e-23)    # J K^-1
        self.ELECTRON_SI      = np.float(1.602176487e-19)  # C
        self.ELECTRONVOLT_SI  = np.float(1.602176487e-19)  # J
        self.ELECTRONMASS_SI  = np.float(9.10938215e-31)   # Kg
        self.HARTREE_SI       = np.float(4.35974394e-18)   # J
        self.RYDBERG_SI       = self.HARTREE_SI/np.float(2.0)   # J
        self.BOHR_RADIUS_SI   = np.float(0.52917720859e-10)# m
        self.AMU_SI           = np.float(1.660538782e-27)  # Kg
        self.C_SI             = np.float(2.99792458e+8)    # m sec^-1
        self.MUNOUGHT_SI      = np.float(self.fpi*1.0e-7)       # N A^-2
        self.EPSNOUGHT_SI     = np.float(1)/(self.MUNOUGHT_SI * self.C_SI**2) # F m^-1
        self.angstrom2au      = np.float(1e-10)/self.BOHR_RADIUS_SI

        # Physical constants, atomic units:
        # AU for "Hartree" atomic units (e = m = hbar = 1)
        # RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
        self.K_BOLTZMANN_AU   = self.K_BOLTZMANN_SI / self.HARTREE_SI
        self.K_BOLTZMANN_RY   = self.K_BOLTZMANN_SI / self.RYDBERG_SI

        self.eps4  = np.float(1.0e-4)
        self.eps6  = np.float(1.0e-6)
        self.eps8  = np.float(1.0e-8)
        self.eps12 = np.float(1.0e-12)
        self.eps14 = np.float(1.0e-14)
        self.eps16 = np.float(1.0e-16)
        self.eps24 = np.float(1.0e-24)
        self.eps32 = np.float(1.0e-32)

        self.bohr_radius_cm = self.BOHR_RADIUS_SI * 100
        self.bohr_radius_angs = self.bohr_radius_cm * 1e8
        self.angstrom_au = 1/self.bohr_radius_angs


def good_fft_order(nr,primes=np.array([2, 3, 5, 7, 11])):
    '''
    FFTXlib/src/fft_support.f90
    an fft order is not good if not implemented (as on IBM with ESSL
    or emplemented butwith awful performances
    output:
        the closest higher number that is good
    '''
    nfftx = 2049 # max allowed fft dimension
    while nr < nfftx:
        mr, factors = primefactors(nr)
        if nr == mr * (primes ** factors).prod() and \
            factors[3]  == 0 and factors[4] == 0 and \
            mr == 1: # fftw, no factors 7 and 11
            break
        nr += 1
    assert nr < nfftx, NotImplementedError('fft box too large')
    return nr
    

def primefactors(nr, primes=np.array([2, 3, 5, 7, 11])):
    factors = np.zeros_like(primes)
    mr = nr
    for i, prime in enumerate(primes):
        maxpwr = int(np.log(mr)/np.log(prime) + .5) + 1
        for p in range(maxpwr):
            if mr == 1:
                return mr, factors
            elif mr % prime == 0:
                mr = mr / prime
                factors[i] += 1
            else:
                break
    return mr, factors

def realspace_grid_init(at, bg, gcutm):
    '''
    FFTXlib/fft_types.f90
         gcutm: ecutrho / (2*pi/a)^2, cut-off for G^2
    '''
    if DEBUG:
        print('realspace grid init input gcutm', gcutm)
    nr1 = int(np.sqrt(gcutm) * np.sqrt(np.square(at[0]).sum())) + 1
    nr2 = int(np.sqrt(gcutm) * np.sqrt(np.square(at[1]).sum())) + 1
    nr3 = int(np.sqrt(gcutm) * np.sqrt(np.square(at[2]).sum())) + 1
    if DEBUG:
        print('realspace_grid_init ', nr1, nr2, nr3)
    nr1, nr2, nr3, ngm, gvect, mill, ngl, gl = grid_set(nr1, nr2, nr3, bg, gcutm)
    #ngm, gvect, mill, ngl, gl= get_ngm(bg, nr1, nr2, nr3, gcutm)
    if DEBUG:
        print('before good fft order ', nr1, nr2, nr3)
    nr1 = good_fft_order(nr1)
    nr2 = good_fft_order(nr2)
    nr3 = good_fft_order(nr3)
    if DEBUG:
        print('after good fft order ', nr1, nr2, nr3)
    return nr1, nr2, nr3, ngm, gvect, mill, ngl, gl

def grid_set(nr1, nr2, nr3, bg, gcut, nproc=1, mype=0):
    '''
    FFTXlib/fft_types.f90
    returns in nr1, nr2, nr3 the minumal 3d real-space fft grid
    required to fit the G-vector spaere with G^2<=gcut
    '''
    g_indices = np.indices((2*nr1+1, 2*nr2+1, nr3+1)).transpose((1,2,3,0)) - np.array((nr1, nr2, 0))
    g_index = g_indices.transpose(3, 0, 1, 2)
    g = np.matmul(g_indices, bg)
    g_length = np.square(g).sum(-1)


    g_mask = np.where(g_length<=gcut, 1, 0)
    g_ix, g_iy, g_iz = np.abs(g_index * g_mask)

    ngm = np.count_nonzero(g_mask)+np.count_nonzero(g_iz)
    g_mask_nonzero = np.nonzero(g_mask)
    g_iz_nonzero = np.nonzero(g_iz)
    ngm = g_mask_nonzero[0].shape[0] + g_iz_nonzero[0].shape[0]

    gvect = np.empty((ngm, 3))
    mill = np.empty((ngm, 3), dtype=np.int32)
    gg = np.empty(ngm)
    gvect[:g_mask_nonzero[0].shape[0]] = g[g_mask_nonzero]
    gvect[g_mask_nonzero[0].shape[0]:] = g[g_iz_nonzero] * np.array((1, 1, -1))

    mill[:g_mask_nonzero[0].shape[0]] = g_indices[g_mask_nonzero]
    mill[g_mask_nonzero[0].shape[0]:] = g_indices[g_iz_nonzero]* np.array((1, 1, -1))
    gg[:g_mask_nonzero[0].shape[0]] = g_length[g_mask_nonzero]
    gg[g_mask_nonzero[0].shape[0]:] = g_length[g_iz_nonzero]

    # sort mill, gvect, gg
    #gg = np.flip(gg)
    #mill = np.flip(mill)
    #gvect = np.flip(gvect)
    ind_sort = np.argsort(gg, kind='heapsort')
    gg = gg[ind_sort]
    mill = mill[ind_sort]
    gvect = gvect[ind_sort]


    nb1 = g_ix.max()
    nb2 = g_iy.max()
    nb3 = g_iz.max()

    nr1 = 2*nb1 + 1
    nr2 = 2*nb2 + 1
    nr3 = 2*nb3 + 1
    return nr1, nr2, nr3, ngm, gvect, mill, ngm, gg
 
def apostrophe_remover(string):
    '''
    remove the "'" in the input string
    '''
    return string.replace("'", "").replace('"', '')

def qe_float(qe_number:str):
    '''
    float number was denoted in QE as like 1d+01, where in python is 1e+01
    '''
    try:
        num = float(qe_number)
    except ValueError:
        num = float(qe_number.replace('d', 'e'))
    return num
