//	BLHelicalUtils.hh
/*
This source file is part of G4beamline, http://g4beamline.muonsinc.com
Copyright (C) 2002-2013 by Tom Roberts, all rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
*/
//
//	A set of inline utility functions for helical dipoles

#ifndef BLHELICALUTILS_HH
#define BLHELICALUTILS_HH

#ifdef G4BL_GSL

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLElement.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLParam.hh"

#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/* how cylindrical coordinates stored */
#define cylPHI  0
#define cylRHO  1
#define cylZ    2

/* how rectangular coordinates stored */
#define rectX    0
#define rectY    1
#define rectZ    2

/* pitch convention */
#define  LEFT_HANDED_THREAD   (-1)
#define RIGHT_HANDED_THREAD   (+1)

#define MuonsInc_PITCH_CONVENTION    RIGHT_HANDED_THREAD

#define MuonsInc_Target_Radius_mm    159.15494

#define ESSENTIALLY_ZERO   1.E-33


G4ThreeVector CYLTOCARTESIAN(G4ThreeVector Bcyl,G4double phi)
{
	    const int PHI=cylPHI,RHO=cylRHO,Z=cylZ;                      /* notation */
            G4ThreeVector Bcart;

            Bcart[rectX] = (Bcyl[RHO] * cos(phi)- Bcyl[PHI]*sin(phi));
	    Bcart[rectY] = (Bcyl[RHO] * sin(phi)+ Bcyl[PHI]*cos(phi));
     	    Bcart[rectZ] = Bcyl[Z];
	    return Bcart;
}



G4ThreeVector SIMPLEFIELD(G4double b,G4double kz,G4double Bsolenoid,G4double model)
{
            const int X=rectX,Y=rectY,Z=rectZ;
            G4ThreeVector SFIELD;

	    SFIELD[X] = b * sin( MuonsInc_PITCH_CONVENTION * kz );
	    SFIELD[Y] = b * cos( MuonsInc_PITCH_CONVENTION * kz );
	    SFIELD[Z] = Bsolenoid + b; 

	    return SFIELD;
}

G4double ICOOLDiv(G4int n, G4double lambda)
{
  G4double k    =2*M_PI/lambda;
  G4double bip  =gsl_sf_bessel_In(n,n);
  G4double bipp =gsl_sf_bessel_In(n-1,n);
  return k*(n*bipp-(1+n)*bip);
}

G4double ICOOLFact(G4int n, G4double krad)
{
  return gsl_sf_gamma(n+1)*gsl_sf_pow_int(2.0/(n*krad),n);
}



G4ThreeVector ICOOLFIELD(G4int n, G4double rho, G4double psi, G4double k, G4double refrad)
{
  const int PHI=cylPHI, RHO=cylRHO, Z=cylZ;
  G4ThreeVector BICOOL;
  G4double bn,bnp;
  G4double kr,sz,cz,fac;
  G4double bref0=-1.0;

  kr  = k*rho;
  sz  = sin(n*psi);
  cz  = cos(n*psi);
  fac = ICOOLFact(n,k*refrad);

  bnp = 0.5*(gsl_sf_bessel_In(n+1,n*kr)+gsl_sf_bessel_In(n-1,n*kr));
  bn  = gsl_sf_bessel_In(n,n*kr);

  BICOOL[RHO]   =  bref0*refrad*fac*k*bnp*(-0.0*cz+1.0*sz);
  BICOOL[Z]     = -bref0*refrad*fac*k*bn*(0.0*sz+1.0*cz);

  if(rho > 0.0){
    BICOOL[PHI] = -BICOOL[Z]/kr;
  }else{
    BICOOL[PHI] = 0.0;
  }

  //  printf("n=%d Bfield Br=%e Bz=%e Bp=%e \n",
  //	 n,BICOOL[RHO],BICOOL[Z],BICOOL[PHI]);
  return BICOOL;
}

G4ThreeVector DIPOLF(G4double bd,G4double k,G4double rho,G4double psiangle)
{
	  const int PHI=cylPHI,RHO=cylRHO,Z=cylZ;        /* notation */
	  G4double t,I0,I1_t,psi;
          G4ThreeVector Bcyl;

	  t= k * rho;
	  I0   = gsl_sf_bessel_In(0,t);
	  I1_t = gsl_sf_bessel_In(1,t);

	  Bcyl[PHI]= 2 * bd * I1_t    * cos(psiangle);
	  Bcyl[RHO]= 2 * bd * (I0-I1_t) * sin(psiangle);

	  Bcyl[Z]=   - k * rho * Bcyl[PHI];
	  return Bcyl;
}


G4ThreeVector QUADF(G4double bprime,G4double k,G4double rho,G4double psiangle)
{
  const char *CoDE="QUADRUPOLEFIELD";
  const int PHI=cylPHI,RHO=cylRHO,Z=cylZ;                /* notation */
  static int messageCount=0;
  G4double t,I1,I2_t,Bprime;
  G4ThreeVector Bcyl;
  if( 0==messageCount++ )  printf("%s: version of 4/27/04\n", CoDE);
  Bprime= bprime/(tesla/meter) * (millimeter/meter) * (tesla/millimeter); 
  
  t= 2 * k * rho;
  I2_t = gsl_sf_bessel_In(2,t);
  I1   = gsl_sf_bessel_In(1,t);

  Bcyl[PHI] = 2/k * Bprime * I2_t * cos(2*psiangle); 
  Bcyl[RHO] = 1/k * Bprime * (I1 - 2*I2_t) * sin(2*psiangle);  // <== fixed! 
  Bcyl[Z] = - k * rho * Bcyl[PHI];

  return Bcyl;
}

#endif // G4BL_GSL
#endif // BLHELICALUTILS_HH
