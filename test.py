from actuator import *
import numpy as np
import math

# Choose an Air Gap
r_gap   = 0.001;            # [m]
# Choose the coil current
c_I     = 1;                # [Amp]
# Choose the wire properties
c_rho   = RHO_CU;           # [ohm m]
c_dw    = 0.0005         # [m]
# Choose the magnet Residual Flux Density
m_Br    = 1;                # [T] 
# Choose the volume of the magnet
m_Vm    = .009**2*np.pi*.01;  # [m^3]


act = Actuator(None, None, r_gap);

# choose the aspect ratio of the and coil and magnet
c_alpha = 2;            # 
m_beta  = 10/9.0;       #

# Calculate the magnet dimensions
m_rm    = (m_Vm/float(np.pi*m_beta))**(1/3.0);  # [m]
m_lm    = m_beta * m_rm;

# Calculate the coil dimensions
c_lw    = 2.513274
c_rc    = 0.01;     # [m]
c_lc    = 0.02;   #
c_Rc    = 0.0105; 
c_Nz    = 40;
c_Nr    = 1;


act.c = Coil(c_I, c_Rc, c_rc, c_lc, c_Nr, c_Nz, c_dw, c_lw, c_alpha);
act.m = Magnet(m_Br, m_lm, m_rm, m_Vm, m_beta);

print act.c
print act.m

for disp in np.linspace(0, 40, 81):
    print str(disp) +","+ str(act.filament_method(disp/1000.0, 10))+","+str(act.shell_method(disp/1000.0));

