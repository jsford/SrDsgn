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

#print "Minimum force to actuate valve:";
#print round(actuation_Force, 5), "[lbs]";
#print round(actuation_Force*453.592, 5), "[grams]";

# ===========================================================
# Find a magnet configuration that can produce enough force
# to actuate the valve.
# ===========================================================

from actuator import *
import numpy as np
import math

# Choose a Coil Impedance
R = 4;                      # [ohm]
# Choose an Air Gap
r_gap   = 0.0005;           # [m]
# Choose the coil current
c_I     = 1;                # [Amp]
# Choose the wire properties
c_rho   = RHO_CU;           # [ohm m]
c_dw    = 0.001;            # [m]
# Choose the magnet Residual Flux Density
m_Br    = 1;              # [T] 
# Choose the volume of the magnet
m_Vm    = .02**3;  #0.0254**3;        # [m^3]


act = Actuator(None, None, r_gap);

for alpha in np.linspace(1, 4, 13):
    for beta in np.linspace(0.5, 4, 15):
        # choose the aspect ratio of the and coil and magnet
        c_alpha = alpha;            # 
        m_beta  = beta;             #

        # Calculate the magnet dimensions
        m_rm    = (m_Vm/float(np.pi*m_beta))**(1/3.0);  # [m]
        m_lm    = m_beta * m_rm;

        # Calculate the coil dimensions
        c_lw    = R*np.pi*(0.5*c_dw)**2/float(c_rho);   # [m]
        c_rc    = m_rm + r_gap;     # [m]
        c_lc    = c_alpha * c_rc;   #
        c_Rc    = np.sqrt(c_lw*c_dw**2/float(np.pi*c_lc) + c_rc**2); 
        c_Nz    = int(math.ceil( c_lc / float(c_dw)));
        c_Nr    = int(math.ceil( (c_Rc - c_rc) / float(c_dw)));


        act.c = Coil(c_I, c_Rc, c_rc, c_lc, c_Nr, c_Nz, c_dw, c_lw, c_alpha);
        act.m = Magnet(m_Br, m_lm, m_rm, m_Vm, m_beta);

        import time
        t0 = time.time();
        f_s = act.shell_method(0.0);
        t1 = time.time();
        f_f = act.filament_method(0.0);
        t2 = time.time();
        print "SHELL TIME: " + str(t1 - t0);
        print "SHELL FORCE: " + str(f_s);
        print "FILAMENT TIME: " + str(t2 - t1);
        print "FILAMENT FORCE: " + str(f_f);

        max_force = 0;
        for disp in np.linspace(0, 40, 81):
            f = act.shell_method(disp/1000.0);
            if f > max_force:
                max_force = f;

        print ("COIL L/R: " + str(c_alpha) + "  \t" +
               "MAGNET L/R: " + str(m_beta) + "  \t" +
               "Max. Force: " + str(max_force));

