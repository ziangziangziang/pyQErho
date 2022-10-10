# qeRho

This Python module reads _QuantumEspresso_ input file and generate a electronic charge density by _super-positioning atom density_ method.

The generated charge density file `charge-density.dat` can be loaded by _QuantumEspresso_ directly.


## Get Starded

### How to install the module

```
pip install qe_rho
```

### How to generate charge density

```
from qe_rho import SAD # super-positioning of atom densities

# initialize SAD density
rho = SAD('tests/pwscf.in')

# save charge-density.dat to folder pwscf.save
rho.saverhog('pwscf.save')

# output real-space charge density in a 3D numpy array
rhor = rho.rho_g2r()
```

You can visualize the charge density with `tests/meshplot.py`:

![RESULT](/tests/my.png)



