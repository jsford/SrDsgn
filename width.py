from actuator import *
import numpy as np
import math

# Choose an Air Gap
r_gap   = 0.001;               # [m]
# Choose the coil current
c_I     = 1;                   # [Amp]
# Choose the wire properties
c_rho   = RHO_CU;              # [ohm m]
c_dw    = 0.0005               # [m]
# Choose the magnet Residual Flux Density
m_Br    = 1;                   # [T] 
# Choose the volume of the magnet
m_Vm    = .009**2*np.pi*.01;   # [m^3]


act = Actuator(None, None, r_gap);

# choose the aspect ratio of the and coil and magnet
c_alpha = 2;            # 
m_beta  = 10/9.0;       #

# Calculate the magnet dimensions
m_rm    = (m_Vm/float(np.pi*m_beta))**(1/3.0);  # [m]
m_lm    = m_beta * m_rm;                        # [m]

# Calculate the coil dimensions
c_lw    = 2.513274
c_rc    = 0.01;     # [m]
c_lc    = 0.02;     # [m]
c_Rc    = 0.0105;   # [m]
c_Nz    = 40;       #
c_Nr    = 1;        #


act.c = Coil(c_I, c_Rc, c_rc, c_lc, c_Nr, c_Nz, c_dw, c_lw, c_alpha);
act.m = Magnet(m_Br, m_lm, m_rm, m_Vm, m_beta);

pts_per_mm = 5;
diff = 20;
length_mm = 80;

length_pts = pts_per_mm*length_mm+1
displ = np.linspace(0, length_mm, length_pts);

f_s = np.array([0.0]*length_pts);

d = 0
for disp in displ:
    f_s[d] = act.shell_method(disp/1000.0);
    d+=1;
#    print str(disp) + "," + str(f_s[d]);
    
pad = np.array([0]*(pts_per_mm*np.abs(diff-length_mm)-1));

if(diff > length_mm and diff <= length_mm*2):
    f_s = np.append(f_s, pad);
    displ = np.linspace(0, diff, pts_per_mm*(diff));
elif(diff < length_mm):
    f_s = np.append(pad, f_s);
    displ = np.linspace(diff-length_mm, length_mm, (2*length_mm-diff)*pts_per_mm);
elif(diff == length_mm):
    displ = np.linspace(0, length_mm, pts_per_mm*length_mm+1);
else:
    print "ERROR: length_mm should be at least twice as long as diff!"

f_s_rev = f_s[::-1];
f_s_sum = f_s + f_s_rev;

import matplotlib.pyplot as plt

plt.xkcd();            # Uncomment for additional whimsy :D
plt.plot(displ, f_s_sum, 'b-', label='Shell Method');
plt.legend(loc='lower right', prop={'size':10});
plt.suptitle('Validating  the Filament and Shell Methods', fontsize=18);
plt.title('Jordan Ford  12/28/14', fontsize=12);
plt.xlabel('Displacement [mm]');
plt.ylabel('Force [N]');
plt.savefig('Width.pdf');
plt.show();
