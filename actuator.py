# ===========================================================
# Jordan Ford
# 12/23/14
# Provide functions for calculating the axial force
# produced by a sleeve coil actuator for a given displacement.
#
# Equations used here come from the following paper:
# 
# Axial force between a thick coil and a cylindrical permanent magnet:
# optimising the geometry of an electromagnetic actuator.
# http://personal.mecheng.adelaide.edu.au/will.robertson/research/2012-magcoil.pdf
# ===========================================================

from scipy import special   # For elliptic integrals
import numpy as np

# Permeability of free space [m kg s^-2 A^-2]
MU_0 = 1.25663706 * 10**-6
# Resistivity of copper [ohm m]
RHO_CU = 1.7 * 10**-8

# Default number of windings used to convert 
# permanent magnets to equivalent coils.
N_M_DEFAULT = 10000;

class Magnet(object):
    def __init__(self, B_r, V_m, beta):
        self.B_r  = B_r;
        self.V_m  = V_m;
        self.beta = beta;
        self.R_m  = (V_m / (np.pi * self.beta) ) ** (1/3.0);
        self.l_m  = self.beta * self.R_m;

class Coil(object):
    def __init__(self, I, r_c, d_w, l_w, alpha, rho=RHO_CU):
        self.I     = I;     # Coil current
        self.d_w   = d_w;   # Wire diam.
        self.l_w   = l_w;   # Wire length
        self.rho   = rho;   # Wire resistivity
        self.alpha = alpha; # Coil aspect ratio
        self.r_c   = r_c;   # Coil inner radius
        self.l_c   = r_c * alpha;  # Coil length
        # Coil outer radius. Equation 21 in paper cited above.
        self.R_c   = np.sqrt(l_w*d_w**2 / (np.pi * self.l_c) + self.r_c**2);
        # Number of layers in the coil.
        self.N_r   = int(np.ceil((self.R_c - self.r_c) / self.d_w));
        # Number of windings per layer.
        self.N_z   = int(np.ceil(self.l_c / self.d_w));

class Actuator(object):
    def __init__(self, magnet, coil, r_gap):
        self.m     = magnet;
        self.c     = coil;
        self.r_gap = r_gap;

    def calc_axial_force(self, z, N_m=N_M_DEFAULT):
        force = 0;
        for n_m in range(1, N_m+1):
            for n_r in range(1, self.c.N_r+1):
                for n_z in range(1, self.c.N_z+1):
                    force += self.force_between_windings(
                             self.__radius_of_winding(n_r),
                             self.m.R_m,
                             z + self.__dist_between_coils(n_m, n_z, N_m),
                             N_m
                             );
        return force;
                    

    # Calculate the force between two coaxial windings 
    # of radius r1 and r2 separated by a distance z using
    # Equation 1 of [Robertson].
    # N_m is the number of windings used by the filament
    # method to discretize the permanent magnet.
    def force_between_windings(self, r1, r2, z, N_m):

        m = 4*r1*r2/float((r1+r2)**2+z**2);    # Argument of the Elliptic Integrals
        K = special.ellipk(m);            # Complete. First kind.
        E = special.ellipe(m);            # Complete. Second kind.

        # Current in the coil.
        I1 = self.c.I;

        # The magnet is treated as an equivalent coil.
        # I2 is the current in the equivalent coil.
        I2 = self.m.B_r * self.m.l_m / float( N_m * MU_0);   

        return ( MU_0 * I1 * I2 * z
                      * np.sqrt(m/float(4*r1*r2))
                      * (K - E*(m/2.0-1)/float(m-1))
               );

    # Equation 4 [Robertson]
    def __radius_of_winding(self, n_r):
        if self.c.N_r == 1:
            return self.c.r_c + self.c.d_w / 2.0;
        return self.c.r_c + (n_r-1)/float(self.c.N_r-1) * (self.c.R_c - self.c.r_c);

    # Equation 5 [Robertson]
    def __dist_between_coils(self, n_m, n_z, N_m):
        return (
                 -1/2.0 * (self.m.l_m + self.c.l_c)
                 + (n_z - 1)/float(self.c.N_z - 1) * self.c.l_c
                 + (n_m - 1)/float(N_m - 1) * self.m.l_m
               );
    
# Verify my implementation by reproducing the data in Fig. 4 of [Robertson]
def verify():
    c = Coil(1, .01, .0005, 2.576, 2);
    m = Magnet(1, np.pi*.009**2*.01, 10/9.0);
    actuateHer = Actuator(m, c, .001);

    print "Displacement [mm]  Force [N]"
    for d in np.linspace(0, 0.04 , 81, endpoint=True):
        print( str(d*1000) + 
               "\t\t   " +
               str(round(actuateHer.calc_axial_force(d, 1000), 6))
             );

