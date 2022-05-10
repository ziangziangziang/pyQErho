from qe_rho import SAD
from meshplot import mesh_plot

a = SAD('pwscf.in')
rhor = a.rho_g2r()
#mesh_plot(rhor, 'firstread.png')
print(rhor.sum())
a.setrhor(rhor)
print(a.rho_g2r().sum())
mesh_plot(a.rho_g2r(), 'secondread.png')

#mesh_plot(a.read_charge('/home/zhanz0j/Documents/pyQErho/chargedens-px.dat'), 'groundtrutu.png')
