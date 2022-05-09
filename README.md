# qeRho

The module is to generate the initial guess charge density for quantum-espresso _pwscf_ calculations. Currently the module supports `nspin=1`, binary output (no HDF5 supported). The output charge density should be sorted in next subversion.
## to install

```
pip install ./
```

## to use

```
# import superposition of atomic density
from qe_rho import SAD

# initialize SAD density
rho = SAD('tests/pwscf.in')

# output real-space charge density
rhor = rho.rho_g2r()

# save quantum espresso readable charge density
rho.saverhog('charge-density.dat')
```



