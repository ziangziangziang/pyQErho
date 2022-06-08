# qeRho



The module is to generate a charge density file `charge-density.dat` for quantum espresso to read.

## to install

```
pip install qe_rho
```

## to use

```
# import superposition of atomic density
from qe_rho import SAD

# initialize SAD density
rho = SAD('tests/pwscf.in')

# save charge-density.dat to folder pwscf.save
rho.saverhog('pwscf.save')

# output real-space charge density in a 3D numpy array
rhor = rho.rho_g2r()
```

## TODO
- [ ] HDF5 support
- [ ] spin support



