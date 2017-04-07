#ifndef COSMOGEOMETRY_H
#define COSMOGEOMETRY_H

#include "AstroConstant.h"
using namespace halo;

double angularDiameterDistance(double z1, double z2, int num=100)
{
   double sum = 0;
   double R2 = 1.0/(1.0+z1);
   double R1 = 1.0/(1.0+z2);
   double dR = (R2-R1)/((double)num);
   double R = R1;
   for (int i = 0 ; i < num ; i++, R += dR)
   {
      double term1 = omegaL*(R*R*R*R); // constant
      double term2 = omegaM*R;   // diluted matter

      double val1 = 1.0/sqrt(term1+term2);

      term1 = omegaL*(R+dR)*(R+dR)*(R+dR)*(R+dR);
      term2 = omegaM*(R+dR); 

      double val2 = 1.0/sqrt(term1+term2);

      sum += ((val1+val2)/2.0)*dR; // trapezium rule
      // note that da=a^2 dz; we multiply the a^2 into the square root terms
   }

   double result = sum*ckms/100./h/(1.0+z2); // 3000 MPc h^-1
   return result; // Mpc / rad  
}

double dVdtheta2dz(double z, int num=100)
// comoving volume [Mpc^3] per solid angle and redshift interval
{
   double da = angularDiameterDistance(0,z,num);

   return (ckms/100./h)*(1.+z)*(1.+z)*da*da/sqrt(omegaM*(1.+z)*(1.+z)*(1.+z)+omegaL);   

}


#endif
