// unit: distance (Mpc/h), mass (Msun/h), time (second)

#ifndef AstroConstantWMAP5_H
#define AstroConstantWMAP5_H

#include <cmath>

namespace halo {
const double h    = 0.72; // times 100 km/s/Mpc
// cosmological model constant
const double omegaM=0.27;
const double omegaL=0.73;
const double omegaK=0.0;
const double omegaB=0.0462;
const double sig8=0.79; //0.81;
const double delc0=1.66;
const double f1 = pow(omegaM,5./9.);
const double f2 = 2.*pow(omegaM,6./11.);
const double ns = 0.96;

// constant number
const double Pi=M_PI;
const double h0=100.*h;
const double G= 6.67384*(1.989E-9)/3.086;
const double G_si=6.67384E-11; // m^3/(kg*s^2)
const double Gs=G_si*(1.9891E30)/(3.08567758E22)*1.E-6;  // = 4.25e-9;  // Mpc/Msol (km/s)^2
const double Gc2=1.4747551/(3.0856E19);
const double ckms = 3.e5; // c in units of 1 km/s
//const double rho_crit=2.775E11*h*h;
//const double rho0=omegaM*rho_crit;
const double AU=1.495E11;
const double pc=3.0856E16; // [m]
const double Mpc=pc*1.E6;
const double m_sun=-26.74;
const double Msol=1.9891e33; // solar mass [g]

// conversation factor
const double deg_rad=Pi/180.;
const double rad_deg=180/Pi;
}

#endif
