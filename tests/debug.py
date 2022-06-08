from qe_rho import SAD
import numpy as np
from meshplot import mesh_plot

a = SAD('pwscf.in')
#a.sort()
a.saverhog('pwscf.save')
rhor = a.rho_g2r()
#print('bg shape', a.cell.bg)
#print('bg shape', np.square(a.cell.bg).sum(-1))
#mesh_plot(rhor, 'my.png', a=a.cell.at)
a.saverhog('pwscf.save')
print(rhor.sum())
#a.setrhor(rhor)
#print(a.rho_g2r().sum())

#mesh_plot(a.read_charge('/home/zhanz0j/Documents/pyQErho/chargedens-px.dat'), 'groundtrutu.png')
rhor = a.read_charge('pwscf.save/charge-density-pw.dat')
#mesh_plot(rhor, 'gt.png', a=a.cell.at)
print(rhor.sum())
