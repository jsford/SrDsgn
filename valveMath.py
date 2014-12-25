# ===========================================================
# Jordan Ford
# 12/22/14
#
# Given lots of valve parameters, 
# find the ideal electromagnet shape and size. 
# ===========================================================


# ===========================================================
# Find the force required to actuate the valve.
# ===========================================================

# working area of valve [in^2]
# How much area disappears when the valve is open?
# Five slits. One inch tall. 1/11 inches wide.
valve_Area = 5*1/11.0*1.0; 

# working pressure [psi]
air_Pressure = 8;

# force exerted on the sliding grid [lbf]
plate_NormalForce = valve_Area * air_Pressure;

# static friction coefficient [unitless]
# teflon - teflon:  0.04
# teflon - steel:   0.05 - 0.20
valve_mu = 0.04;

# force required to actuate valve [lbf]
actuation_Force = valve_mu * plate_NormalForce;

print "Minimum force to actuate valve:";
print round(actuation_Force, 5), "[lbs]";
print round(actuation_Force*453.592, 5), "[grams]";

# ===========================================================
# Find a magnet configuration that can produce enough force
# to actuate the valve.
# ===========================================================

from actuator import *

# Choose a Coil Impedance
R = 4;  # [ohm]

c_I     = 1;                         # [Amp]
c_rc    = 0.01;                      # [m]
c_dw    = 0.00016;                   # [m]
c_lw    = R*np.pi*(0.5*c_dw)**2/float(RHO_CU);  # [m]
c_alpha = 2;                         #
m_Br    = 1;                         # [T]
m_Vm    = 3.218e-6                   # [m^3] .25 in. radius 1 in. height
m_beta  = 10/9.0;                    #
 
disp    = 0.001;

c = Coil(c_I, c_rc, c_dw, c_lw, c_alpha);
m = Magnet(m_Br, m_Vm, m_beta);
act = Actuator(m, c, disp);



print c.l_w
print c.N_r*c.N_z
















