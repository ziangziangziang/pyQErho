from qe_rho import SAD
import numpy as np
from meshplot import mesh_plot

a = SAD('pwscf.in')
#a.sort()
a.saverhog('pwscf.save')
rhor = a.rho_g2r()
#print('bg shape', a.cell.bg)
#print('bg shape', np.square(a.cell.bg).sum(-1))
mesh_plot(rhor, 'my.png', a=a.cell.at)
