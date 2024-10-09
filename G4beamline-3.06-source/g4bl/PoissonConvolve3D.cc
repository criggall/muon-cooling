//	PoissonConvolve3D.cpp - solve Poisson equation using MUDPACK, mud3sp.f
//
//	compile with "-DMAIN" to create a test program.

#ifdef G4BL_FFTW

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#include "PoissonConvolve3D.hh"		// includes fftw3.h

PoissonConvolve3D::PoissonConvolve3D(int _nx, int _ny, int _nz, float _xMin,
	     			float _xMax, float _yMin, float _yMax,
	     			float _zMin, float _zMax) : approx() {
	nx=_nx; ny=_ny; nz=_nz;
	nc = 4*nx*ny*(nz+1); // double size, except last dimension
	status = UNINITIALIZED;
	array = (float*)fftwf_malloc(sizeof(float)*8*nx*ny*nz);
	transf = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nc);
	greensFunction = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nc);
	sumX = 0.0; sumY = 0.0; sumZ = 0.0; totalCharge = 0.0;
	putMode = PRORATED;

	fftwf_set_timelimit(10.0); // otherwise can take up to 226 sec (!)

	gfPlan = fftwf_plan_dft_r2c_3d(2*nx,2*ny,2*nz,array,greensFunction,
							FFTW_ESTIMATE);
	plan1 = fftwf_plan_dft_r2c_3d(2*nx,2*ny,2*nz,array,transf,
							FFTW_ESTIMATE);
	plan2 = fftwf_plan_dft_c2r_3d(2*nx,2*ny,2*nz,transf,array,
							FFTW_ESTIMATE);

	updatePosition(_xMin,_xMax,_yMin,_yMax,_zMin,_zMax);
}

float PoissonConvolve3D::phi(float x, float y, float z)
{
	assert(status==PHI);
	int ix = (int)floor((x-xMin)/dx);
	int iy = (int)floor((y-yMin)/dy);
	int iz = (int)floor((z-zMin)/dz);
	if(ix < 0 || iy < 0 || iz < 0 || ix+1 >= nx || iy+1 >= ny || 
								iz+1 >= nz) {
		// outside grid, use the approximation
		float v=0.0;
		G4ThreeVector pos(x,y,z);
		for(unsigned i=0; i<approx.size(); ++i) {
			G4ThreeVector pt = approx[i].getPosition();
			pt -= pos;
			double r = pt.mag();
			v += approx[i].c / r / (4.0*M_PI);
		}
		return v;
	}

	// within grid, interpolate linearly among 8 surrounding grid points
	float fx = (x - (xMin+ix*dx))/dx;
	assert(-0.01 <= fx && fx <= 1.01);
	float fy = (y - (yMin+iy*dy))/dy;
	assert(-0.01 <= fy && fy <= 1.01);
	float fz = (z - (zMin+iz*dz))/dz;
	assert(-0.01 <= fz && fz <= 1.01);
	float gx=1.0-fx, gy=1.0-fy, gz=1.0-fz;

	float v =
		gx*gy*gz*array[index(ix,iy,iz)] +
		fx*gy*gz*array[index(ix+1,iy,iz)] +
		gx*fy*gz*array[index(ix,iy+1,iz)] +
		fx*fy*gz*array[index(ix+1,iy+1,iz)] +
		gx*gy*fz*array[index(ix,iy,iz+1)] +
		fx*gy*fz*array[index(ix+1,iy,iz+1)] +
		gx*fy*fz*array[index(ix,iy+1,iz+1)] +
		fx*fy*fz*array[index(ix+1,iy+1,iz+1)];

	return v;
}

bool PoissonConvolve3D::putCharge(float x, float y, float z, float c)
{
	assert(status==CHARGE);
	int ix = (int)floor((x-xMin)/dx);
	int iy = (int)floor((y-yMin)/dy);
	int iz = (int)floor((z-zMin)/dz);
	// (omit boundaries)
	if(ix<=0 || iy<=0 || iz<=0 || ix+1>=nx || iy+1>=ny || iz+1>=nz) {
		return false;
	}

	sumX += c*x;
	sumY += c*y;
	sumZ += c*z;
	totalCharge += c;

	c /= 4.0*M_PI;

	switch(putMode) {
	case LOWEST:
		array[index(ix,iy,iz)] += c;
		break;
	case NEAREST:
		if(x > xMin+ix*dx+dx/2.0) ++ix; // could be upper boundary
		if(y > yMin+iy*dy+dy/2.0) ++iy; // but is guaranteed to 
		if(z > zMin+iz*dz+dz/2.0) ++iz; // be in range.
		array[index(ix,iy,iz)] += c;
		break;
	case PRORATED:
		{ float fx = (x - (xMin+ix*dx))/dx;
		  assert(-0.01 <= fx && fx <= 1.01);
		  float fy = (y - (yMin+iy*dy))/dy;
		  assert(-0.01 <= fy && fy <= 1.01);
		  float fz = (z - (zMin+iz*dz))/dz;
		  assert(-0.01 <= fz && fz <= 1.01);
		  float gx=1.0-fx, gy=1.0-fy, gz=1.0-fz;
		  array[index(ix,iy,iz)] += gx*gy*gz*c;
		  array[index(ix+1,iy,iz)] += fx*gy*gz*c;
		  array[index(ix,iy+1,iz)] += gx*fy*gz*c;
		  array[index(ix+1,iy+1,iz)] += fx*fy*gz*c;
		  array[index(ix,iy,iz+1)] += gx*gy*fz*c;
		  array[index(ix+1,iy,iz+1)] += fx*gy*fz*c;
		  array[index(ix,iy+1,iz+1)] += gx*fy*fz*c;
		  array[index(ix+1,iy+1,iz+1)] += fx*fy*fz*c;
		}
		break;
	}

	return true;
}

void PoissonConvolve3D::zeroRhs()
{
	sumX = sumY = sumZ = totalCharge = 0;
	for(int i=0; i<2*nx; ++i) {
		for(int j=0; j<2*ny; ++j) {
			for(int k=0; k<2*nz; ++k) {
				array[index(i,j,k)] = 0.0;
			}
		}
	}

	// indicate CHARGE is now in array[]
	status = CHARGE;
}

G4ThreeVector PoissonConvolve3D::getE(float x, float y, float z)
{
	assert(status==PHI);
	int ix = (int)floor((x-xMin)/dx);
	int iy = (int)floor((y-yMin)/dy);
	int iz = (int)floor((z-zMin)/dz);
	if(ix < 0 || iy < 0 || iz < 0 || ix+1 >= nx || iy+1 >= ny || 
								iz+1 >= nz) {
		return approxE(x,y,z);
	}

#ifdef SIMPLE_LINEAR
	// Has the problem that the derivative is constant in each box.
	// variables used in mathematica notebook "interpolationOn3DGrid.nb"
	// -- it differentiates the 8-point interpolation formula.
	double x0=xMin+ix*dx, y0=yMin+iy*dy, z0=zMin+iz*dz;
	double x1=x0+dx, y1=y0+dy, z1=z0+dz;
	double v000=array[index(ix,iy,iz)];
	double v100=array[index(ix+1,iy,iz)];
	double v010=array[index(ix,iy+1,iz)];
	double v110=array[index(ix+1,iy+1,iz)];
	double v001=array[index(ix,iy,iz+1)];
	double v101=array[index(ix+1,iy,iz+1)];
	double v011=array[index(ix,iy+1,iz+1)];
	double v111=array[index(ix+1,iy+1,iz+1)];

	// E = - grad array[]
	double Ex=-(1.0/(dx*dy*dz))*(dy*(dz*(-v000+v100)+(v000-v001-v100+v101)
		*(z-z0))+(y-y0)*(dz*(v000-v010-v100+v110)-
		(v000-v001-v010+v011-v100+v101+v110-v111)*(z-z0)));
	double Ey=-(1.0/(dx*dy*dz))*(dx*(dz*(-v000+v010)+(v000-v001-v010+v011)
		*(z-z0))+(x-x0)*(dz*(v000-v010-v100+v110)-
		(v000-v001-v010+v011-v100+v101+v110-v111)*(z-z0)));
	double Ez=-(1.0/(dx*dy*dz))*(dx*(dy*(-v000+v001)+(v000-v001-v010+v011)
		*(y-y0))+(x-x0)*(dy*(v000-v001-v100+v101)-
		(v000-v001-v010+v011-v100+v101+v110-v111)*(y-y0)));
#else
	// linear interpolation between derivatives in adjacent cells.
	// Associate the derivative in each grid box with the center of the box;
	// interpolate the derivative between the current box and the closer
	// adjacent one, or the boundary if this box is an edge.
	float fx = (x - (xMin+ix*dx))/dx;
	assert(-0.01 <= fx && fx <= 1.01);
	float fy = (y - (yMin+iy*dy))/dy;
	assert(-0.01 <= fy && fy <= 1.01);
	float fz = (z - (zMin+iz*dz))/dz;
	assert(-0.01 <= fz && fz <= 1.01);
	float gx=1.0-fx, gy=1.0-fy, gz=1.0-fz;

	double v000=array[index(ix,iy,iz)];
	double v100=array[index(ix+1,iy,iz)];
	double v010=array[index(ix,iy+1,iz)];
	double v110=array[index(ix+1,iy+1,iz)];
	double v001=array[index(ix,iy,iz+1)];
	double v101=array[index(ix+1,iy,iz+1)];
	double v011=array[index(ix,iy+1,iz+1)];
	double v111=array[index(ix+1,iy+1,iz+1)];

	// partial d/dx
	double v0=gy*gz*v000 + fy*gz*v010 + gy*fz*v001 + fy*fz*v011;
	double v1=gy*gz*v100 + fy*gz*v110 + gy*fz*v101 + fy*fz*v111;
	double Ex=-(v1-v0)/dx; // for this grid box
	if(fx > 0.5) {
		if(ix+1 == nx-1) {
			// interpolate between center of box and xMax
			double Ex2=approxE(xMax,y,z).x();
			double f=2.0*fx-1.0;
			Ex = (1.0-f)*Ex + f*Ex2;
		} else {
			// interpolate between this box and next higher one
			double v200=array[index(ix+2,iy,iz)];
			double v210=array[index(ix+2,iy+1,iz)];
			double v201=array[index(ix+2,iy,iz+1)];
			double v211=array[index(ix+2,iy+1,iz+1)];
			double v2=gy*gz*v200 + fy*gz*v210 + gy*fz*v201 + 
								fy*fz*v211;
			double f=fx-0.5;
			Ex = (1.0-f)*Ex - f*(v2-v1)/dx;
		}
	} else {
		if(ix == 0) {
			// interpolate between center of box and xMin
			double Ex2=approxE(xMin,y,z).x();
			double f=fx*2.0;
			Ex = (1.0-f)*Ex2 + f*Ex;
		} else {
			// interpolate between this box and next lower one
			double v200=array[index(ix-1,iy,iz)];
			double v210=array[index(ix-1,iy+1,iz)];
			double v201=array[index(ix-1,iy,iz+1)];
			double v211=array[index(ix-1,iy+1,iz+1)];
			double v2=gy*gz*v200 + fy*gz*v210 + gy*fz*v201 + 
								fy*fz*v211;
			double f=fx+0.5;
			Ex = -(1.0-f)*(v0-v2)/dx + f*Ex;
		}
	}

	// partial d/dy
	v0=gx*gz*v000 + fx*gz*v100 + gx*fz*v001 + fx*fz*v101;
	v1=gx*gz*v010 + fx*gz*v110 + gx*fz*v011 + fx*fz*v111;
	double Ey=-(v1-v0)/dy; // for this grid box
	if(fy > 0.5) {
		if(iy+1 == ny-1) {
			// interpolate between center of box and yMax
			double Ey2=approxE(x,yMax,z).y();
			double f=2.0*fy-1.0;
			Ey = (1.0-f)*Ey + f*Ey2;
		} else {
			// interpolate between this box and next higher one
			double v020=array[index(ix,iy+2,iz)];
			double v120=array[index(ix+1,iy+2,iz)];
			double v021=array[index(ix,iy+2,iz+1)];
			double v121=array[index(ix+1,iy+2,iz+1)];
			double v2=gx*gz*v020 + fx*gz*v120 + gx*fz*v021 + 
								fx*fz*v121;
			double f=fy-0.5;
			Ey = (1.0-f)*Ey - f*(v2-v1)/dy;
		}
	} else {
		if(iy == 0) {
			// interpolate between center of box and yMin
			double Ey2=approxE(x,yMin,z).y();
			double f=fy*2.0;
			Ey = (1.0-f)*Ey2 + f*Ey;
		} else {
			// interpolate between this box and next lower one
			double v020=array[index(ix,iy-1,iz)];
			double v120=array[index(ix+1,iy-1,iz)];
			double v021=array[index(ix,iy-1,iz+1)];
			double v121=array[index(ix+1,iy-1,iz+1)];
			double v2=gx*gz*v020 + fx*gz*v120 + gx*fz*v021 + 
								fx*fz*v121;
			double f=fy+0.5;
			Ey = -(1.0-f)*(v0-v2)/dy + f*Ey;
		}
	}

	// partial d/dz
	v0=gx*gy*v000 + fx*gy*v100 + gx*fy*v010 + fx*fy*v110;
	v1=gx*gy*v001 + fx*gy*v101 + gx*fy*v011 + fx*fy*v111;
	double Ez=-(v1-v0)/dz; // for this grid box
	if(fz > 0.5) {
		if(iz+1 == nz-1) {
			// interpolate between center of box and zMax
			double Ez2=approxE(x,y,zMax).z();
			double f=2.0*fz-1.0;
			Ez = (1.0-f)*Ez + f*Ez2;
		} else {
			// interpolate between this box and next higher one
			double v002=array[index(ix,iy,iz+2)];
			double v102=array[index(ix+1,iy,iz+2)];
			double v012=array[index(ix,iy+1,iz+2)];
			double v112=array[index(ix+1,iy+1,iz+2)];
			double v2=gx*gy*v002 + fx*gy*v102 + gx*fy*v012 + 
								fx*fy*v112;
			double f=fz-0.5;
			Ez = (1.0-f)*Ez - f*(v2-v1)/dz;
		}
	} else {
		if(iz == 0) {
			// interpolate between center of box and zMin
			double Ez2=approxE(x,y,zMin).z();
			double f=fz*2.0;
			Ez = (1.0-f)*Ez2 + f*Ez;
		} else {
			// interpolate between this box and next lower one
			double v002=array[index(ix,iy,iz-1)];
			double v102=array[index(ix+1,iy,iz-1)];
			double v012=array[index(ix,iy+1,iz-1)];
			double v112=array[index(ix+1,iy+1,iz-1)];
			double v2=gx*gy*v002 + fx*gy*v102 + gx*fy*v012 + 
								fx*fy*v112;
			double f=fz+0.5;
			Ez = -(1.0-f)*(v0-v2)/dz + f*Ez;
		}
	}
#endif

	return G4ThreeVector(Ex,Ey,Ez);
}

G4ThreeVector PoissonConvolve3D::approxE(double x, double y, double z)
{
	double tol = ((dx<dy?dy:dx)<dz ? (dx<dy?dy:dx) : dz) * 0.05;
	G4ThreeVector v(0.0,0.0,0.0);
	G4ThreeVector pos(x,y,z);
	for(unsigned i=0; i<approx.size(); ++i) {
		G4ThreeVector p = approx[i].getPosition();
		p -= pos; // note - sign
		double r=p.mag();
		if(r < tol) continue; // should never happen
		double denom=4.0*M_PI*r*r*r;
		double c=approx[i].getCharge();
		G4ThreeVector E(c*p.x()/denom,c*p.y()/denom,c*p.z()/denom);
		v -= E; // note - sign above
	}
	return v;
}

float PoissonConvolve3D::approxPhi(double x, double y, double z)
{
	double tol = ((dx<dy?dy:dx)<dz ? (dx<dy?dy:dx) : dz) * 0.05;
	double phi=0.0;
	G4ThreeVector pos(x,y,z);
	for(unsigned i=0; i<approx.size(); ++i) {
		G4ThreeVector p = approx[i].getPosition();
		p -= pos; // note - sign
		double r=p.mag();
		if(r < tol) r = tol; // should never happen
		double c=approx[i].getCharge();
		phi += c/r;
	}
	return phi/(4.0*M_PI);
}

void PoissonConvolve3D::defaultApproximation()
{
	Charge c(sumX/totalCharge,sumY/totalCharge,sumZ/totalCharge,
								totalCharge);
	approx.clear();
	approx.push_back(c);
}

void PoissonConvolve3D::computeGreensFunction()
{
	// fill in the physical region
	// include the 1/N overall normalization
	// (Need one more sample than usual, for the symmetry.)
	float norm=0.125f/(float)(nx*ny*nz);
	for(int ix=0; ix<=nx; ++ix) {
		float x2=ix*dx*ix*dx;
		for(int iy=0; iy<=ny; ++iy) {
			float y2=iy*dy*iy*dy;
			for(int iz=0; iz<=nz; ++iz) {
				float z2=iz*dz*iz*dz;
				float r=sqrtf(x2+y2+z2);
				if(ix==0 && iy==0 && iz==0)
					r = (dx+dy+dz)/3.0;
				array[index(ix,iy,iz)] = norm/r;
			}
		}
	}

	// now fill in the un-physical region via symmetry
	for(int ix=0; ix<2*nx; ++ix) {
		int jx=ix;
		if(jx > nx) jx=2*nx-jx;
		assert(jx>=0 && jx<=nx);
		for(int iy=0; iy<2*ny; ++iy) {
			int jy=iy;
			if(jy > ny) jy=2*ny-jy;
			assert(jy>=0 && jy<=ny);
			for(int iz=0; iz<2*nz; ++iz) {
				if(ix<=nx && iy<=ny && iz<=nz) continue;
				int jz=iz;
				if(jz > nz) jz=2*nz-jz;
				assert(jz>=0 && jz<=nz);
				array[index(ix,iy,iz)] = array[index(jx,jy,jz)];
			}
		}
	}

	// compute the FFT of the Green's function
	fftwf_execute(gfPlan);

	// indicate we clobbered array[]
	status = UNINITIALIZED;
}

void PoissonConvolve3D::solve()
{
	assert(status==CHARGE);

	// compute the FFT of array[], which contains the charge,
	// into transf[].
	fftwf_execute(plan1);

	// multiply by the FFT of the Green's function
	for(int i=0; i<nc; ++i) {
		double tmp = transf[i][0]*greensFunction[i][0] -
					transf[i][1]*greensFunction[i][1];
		transf[i][1] = transf[i][1]*greensFunction[i][0] +
					transf[i][0]*greensFunction[i][1];
		transf[i][0] = tmp;
	}

	// inverse transform transf[] to get phi into array[]
	fftwf_execute(plan2); // transf[] is overwritten, but don't care

	// indicate a valid phi in array[]
	status = PHI;
}

void PoissonConvolve3D::printRhs()
{
	assert(status==CHARGE);
	for(int ix=0; ix<nx; ++ix) {
		for(int iy=0; iy<ny; ++iy) {
			for(int iz=0; iz<nz; ++iz) {
				int k=index(ix,iy,iz);
				if(array[k] == 0.0) continue;
				printf("rhs[%d] = %.3f  ix=%d iy=%d iz=%d\n",
					k,array[k],ix,iy,iz);
			}
		}
	}
}

void PoissonConvolve3D::printPhi()
{
	assert(status==PHI);
	for(int ix=0; ix<nx; ++ix) {
		for(int iy=0; iy<ny; ++iy) {
			for(int iz=0; iz<nz; ++iz) {
				int k=index(ix,iy,iz);
				printf("phi[%d] = %.3f  ix=%d iy=%d iz=%d\n",
					k,array[k],ix,iy,iz);
			}
		}
	}
}


#ifdef MAIN 
//	main program for testing

#include <sys/time.h>


/// gettimeus() returns the current time in microseconds.
static long long  gettimeus()
{
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_sec*1000000LL + tv.tv_usec;
}

/// gausran() is a simple Gaussian random number generator, limited to 6*sigma.
double gausran(double mean, double sigma)
{
	double v;
	v = drand48()+drand48()+drand48()+drand48()+drand48()+drand48()
	   +drand48()+drand48()+drand48()+drand48()+drand48()+drand48() - 6.0;
	return v*sigma+mean;
}

G4ThreeVector gaussianE(double x, double y, double z, double c, double sigma)
{
	double r=sqrt(x*x+y*y+z*z);
	c /= 4.0*M_PI;
	double Ex=-c*x*(exp(-(r*r/(2*sigma*sigma)))*sqrt(2/M_PI)/(r*r*sigma) -
			(erf(r/(sqrt(2.0)*sigma)))/(r*r*r));
	double Ey=-c*y*(exp(-(r*r/(2*sigma*sigma)))*sqrt(2/M_PI)/(r*r*sigma) -
			(erf(r/(sqrt(2.0)*sigma)))/(r*r*r));
	double Ez=-c*z*(exp(-(r*r/(2*sigma*sigma)))*sqrt(2/M_PI)/(r*r*sigma) -
			(erf(r/(sqrt(2.0)*sigma)))/(r*r*r));
	return G4ThreeVector(Ex,Ey,Ez);
}

int main(int argc, char *argv[])
{
	int N=129;
	double sigma=0.2;

	long long start;
	int ierror;

	printf("Test 0: single charge at the origin\n");
	for(int i=0; i<4; ++i) {
	//for(int i=0; i<5; ++i) {
	    switch(i) {
	    case 0:	N = 17;		break;
	    case 1:	N = 33;		break;
	    case 2:	N = 65;		break;
	    case 3:	N = 129;	break;
	    case 4:	N = 257;	break;
	    }
	    start=gettimeus();
	    PoissonConvolve3D poisson(N,N,N,-1.0,1.0,-1.0,1.0,-1.0,1.0);
	    printf("Constructor(N=%d) took %lld us\n",N,gettimeus()-start);
	    start=gettimeus();
	    poisson.zeroRhs();
	    double totalCharge=1000.0;
	    if(!poisson.putCharge(0.0,0.0,0.0,totalCharge)) {
	    	fprintf(stderr,"ERROR\n");
		exit(9);
	    }
	    poisson.defaultApproximation();
	    printf("Setup took %lld us\n",gettimeus()-start);
	    start = gettimeus();
	    poisson.solve();
	    printf("solve took %lld us;  N=%d\n", gettimeus()-start,N);
	    fflush(stdout);
	    if(true) {
		char fname[256];
		sprintf(fname,"Point,N=%d.txt",N);
		FILE *out = fopen(fname,"w");
		assert(out != 0);
		fprintf(out,"# %s\n",fname);
		fprintf(out,"#x y z r PHImodel PHIexact Emodel Eexact Eerror\n");
		double dx = 0.01;
		for(double r=dx; r<3.0; r+=dx) {
			double exact = poisson.approxPhi(r,0.0,0.0);
			for(int m=0; m<10; ++m) {
				double theta=M_PI*drand48();
				double phi=2.0*M_PI*drand48();
				double x=r*sin(theta)*cos(phi);
				double y=r*sin(theta)*sin(phi);
				double z=r*cos(theta);
assert(fabsf(sqrtf(x*x+y*y+z*z)-r) < 1e-6);
				double v = (double)poisson.phi(x,y,z);
				double denom = r*r*r*4.0*M_PI;
				G4ThreeVector Emodel=poisson.getE(x,y,z);
				G4ThreeVector Eexact=poisson.approxE(x,y,z);
//G4ThreeVector E2(totalCharge/4.0/M_PI*x/r/r/r,totalCharge/4.0/M_PI*y/r/r/r,totalCharge/4.0/M_PI*z/r/r/r);
//printf("r=%.2f Eexact=%.2f,%.2f,%.2f  E2=%.2f,%.2f,%.2f  Emodel=%.2f,%.2f,%.2f\n",r,Eexact.x(),Eexact.y(),Eexact.z(),E2.x(),E2.y(),E2.z(),Emodel.x(),Emodel.y(),Emodel.z());
				G4ThreeVector Eerror=Emodel;
				Eerror -= Eexact;
				fprintf(out,"%.4f %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f\n",
					x,y,z,r,v,exact,
					Emodel.mag(),Eexact.mag(),Eerror.mag());
			}
		}
		fclose(out);
		// scan along z, with constant x and constant y
		double x=0.9, y=0.0, dz=0.0017;
		sprintf(fname,"Point,N=%d,x=%.1f,y=%.1f.txt",N,x,y);
		out = fopen(fname,"w");
		assert(out != 0);
		fprintf(out,"# %s\n",fname);
		fprintf(out,"#x y z r PHImodel PHIexact Emodel Eexact Eerror Eexact_x Emodel_x\n");
		for(double z=-3.0; z<=3.0; z+=dz) {
			double exact = poisson.approxPhi(x,y,z);
			double r=sqrt(x*x+y*y+z*z);
			double v = (double)poisson.phi(x,y,z);
			G4ThreeVector Emodel=poisson.getE(x,y,z);
			G4ThreeVector Eexact=poisson.approxE(x,y,z);
			G4ThreeVector Eerror=Emodel;
			Eerror -= Eexact;
			fprintf(out,"%.4f %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
				x,y,z,r,v,exact,
				Emodel.mag(),Eexact.mag(),Eerror.mag(),
				Eexact.x(),Emodel.x());
		}
		fclose(out);
	    }
	}

	printf("Test 1: exact Gaussian charge distribution at the origin.\n");
	N = 65;
	PoissonConvolve3D poisson(N,N,N,-1.0,1.0,-1.0,1.0,-1.0,1.0);
	static double sig[8]={0.1,0.12,0.16,0.2,0.25,0.3,0.4,0.5};
	for(int isigma=0; isigma<8; ++isigma) {
	    sigma = sig[isigma];
	    start=gettimeus();
	    poisson.zeroRhs();
	    double totalCharge=0.0;
	    for(int ix=0; ix<N; ++ix) {
		for(int iy=0; iy<N; ++iy) {
			for(int iz=0; iz<N; ++iz) {
				G4ThreeVector p=poisson.getPosition(ix,iy,iz);
				double r=p.mag();
				double c=exp(-r*r/sigma/sigma/2.0);
				if(poisson.putCharge(p.x(),p.y(),p.z(),c))
					totalCharge += c;
			}
		}
	    }
	    poisson.defaultApproximation();
	    printf("Setup took %lld us\n",gettimeus()-start);
	    start = gettimeus();
	    poisson.solve();
	    printf("solve took %lld us;  N=%d\n", gettimeus()-start,N);
	    if(true) {
		char fname[64];
		sprintf(fname,"Gaussian,sigma=%.2f,N=%d.txt",sigma,N);
		FILE *out = fopen(fname,"w");
		assert(out != 0);
		fprintf(out,"#x y z r PHImodel PHIexact Emodel Eexact Eerror\n");
		double dx = 0.01;
		double sqrt2_sigma=sqrt(2.0)*sigma;
		for(double r=dx; r<3.0; r+=dx) {
			double exact = totalCharge/r * erf(r/sqrt2_sigma);
exact /= 4.0*M_PI;
			for(int m=0; m<5; ++m) {
				double theta=M_PI*drand48();
				double phi=2.0*M_PI*drand48();
				double x=r*sin(theta)*cos(phi);
				double y=r*sin(theta)*sin(phi);
				double z=r*cos(theta);
				double v = (double)poisson.phi(x,y,z);
				G4ThreeVector Emodel=poisson.getE(x,y,z);
				G4ThreeVector Eexact=gaussianE(x,y,z,
						totalCharge,sigma);
				G4ThreeVector Eerror=Emodel;
				Eerror -= Eexact;
				fprintf(out,"%.4f %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f\n",
					x,y,z,r,v,exact,
					Emodel.mag(),Eexact.mag(),Eerror.mag());
			}
		}
		fclose(out);
		if(sigma == 0.2) {
			static double points[16][3] = {
					{ 0.5, 0.5, 0.5},
					{ 0.5, 0.5,-0.5},
					{ 0.5,-0.5, 0.5},
					{ 0.5,-0.5,-0.5},
					{-0.5, 0.5, 0.5},
					{-0.5, 0.5,-0.5},
					{-0.5,-0.5, 0.5},
					{-0.5,-0.5,-0.5},
					{ 2.0, 2.0, 2.0},
					{ 2.0, 2.0,-2.0},
					{ 2.0,-2.0, 2.0},
					{ 2.0,-2.0,-2.0},
					{-2.0, 2.0, 2.0},
					{-2.0, 2.0,-2.0},
					{-2.0,-2.0, 2.0},
					{-2.0,-2.0,-2.0},
				};
			for(int i=0; i<16; ++i) {
				G4ThreeVector p(points[i][0],points[i][1],points[i][2]);
				G4ThreeVector E=poisson.getE(p.x(),p.y(),p.z());
				printf("Sign check: x=%4.1f y=%4.1f z=%4.1f: "
					" E: %7.0f %7.0f %7.0f\n",
					p.x(),p.y(),p.z(),E.x(),E.y(),E.z());
			}
		}
	    }
	}
exit(0);

	printf("Test 2: 1M charges, Gaussian distribution, PRORATED\n");
	sigma = 0.2;
	start=gettimeus();
	poisson.zeroRhs();
	poisson.setPutMode(PoissonConvolve3D::PRORATED);
	double totalCharge=0.0;
	for(int i=0; i<1000000; ++i) {
		double x = gausran(0.0,sigma);
		double y = gausran(0.0,sigma);
		double z = gausran(0.0,sigma);
		if(poisson.putCharge(x,y,z,0.00001)) totalCharge += 0.00001;
	}
	poisson.defaultApproximation();
	printf("Setup took %lld us\n",gettimeus()-start);
	start = gettimeus();
	poisson.solve();
	printf("solve took %lld us;  N=%d\n", gettimeus()-start,N);
	start = gettimeus();
	poisson.solve();
	printf("solve took %lld us\n",gettimeus()-start);
	if(true) {
		FILE *out = fopen("MillionCharges-PRORATED.txt","w");
		assert(out != 0);
		fprintf(out,"#x y z r PHImodel PHIexact Emodel Eexact Eerror\n");
		double dx = 1.0/(N-1)/3.0;
		double sqrt2_sigma=sqrt(2.0)*sigma;
		for(double r=dx; r<3.0; r+=dx) {
			double exact = totalCharge/r * erf(r/sqrt2_sigma);
			for(int m=0; m<5; ++m) {
				double theta=M_PI*drand48();
				double phi=2.0*M_PI*drand48();
				double x=r*sin(theta)*cos(phi);
				double y=r*sin(theta)*sin(phi);
				double z=r*cos(theta);
				double v = (double)poisson.phi(x,y,z);
				G4ThreeVector Emodel=poisson.getE(x,y,z);
				G4ThreeVector Eexact=gaussianE(x,y,z,
						totalCharge,sigma);
				G4ThreeVector Eerror=Emodel;
				Eerror -= Eexact;
				fprintf(out,"%.4f %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f\n",
					x,y,z,r,v,exact,
					Emodel.mag(),Eexact.mag(),Eerror.mag());
			}
		}
		fclose(out);
	}

	printf("Test 3: 1M charges, Gaussian distribution, NEAREST\n");
	sigma = 0.2;
	start=gettimeus();
	poisson.zeroRhs();
	poisson.setPutMode(PoissonConvolve3D::NEAREST);
	totalCharge=0.0;
	for(int i=0; i<1000000; ++i) {
		double x = gausran(0.0,sigma);
		double y = gausran(0.0,sigma);
		double z = gausran(0.0,sigma);
		if(poisson.putCharge(x,y,z,0.00001)) totalCharge += 0.00001;
	}
	poisson.defaultApproximation();
	printf("Setup took %lld us\n",gettimeus()-start);
	start = gettimeus();
	poisson.solve();
	printf("solve took %lld us;  N=%d\n", gettimeus()-start,N);
	start = gettimeus();
	poisson.solve();
	printf("solve took %lld us\n",gettimeus()-start);
	if(true) {
		FILE *out = fopen("MillionCharges-NEAREST.txt","w");
		assert(out != 0);
		fprintf(out,"#x y z r PHImodel PHIexact Emodel Eexact Eerror\n");
		double dx = 1.0/(N-1)/3.0;
		double sqrt2_sigma=sqrt(2.0)*sigma;
		for(double r=dx; r<3.0; r+=dx) {
			double exact = totalCharge/r * erf(r/sqrt2_sigma);
			for(int m=0; m<5; ++m) {
				double theta=M_PI*drand48();
				double phi=2.0*M_PI*drand48();
				double x=r*sin(theta)*cos(phi);
				double y=r*sin(theta)*sin(phi);
				double z=r*cos(theta);
				double v = (double)poisson.phi(x,y,z);
				G4ThreeVector Emodel=poisson.getE(x,y,z);
				G4ThreeVector Eexact=gaussianE(x,y,z,
						totalCharge,sigma);
				G4ThreeVector Eerror=Emodel;
				Eerror -= Eexact;
				fprintf(out,"%.4f %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f\n",
					x,y,z,r,v,exact,
					Emodel.mag(),Eexact.mag(),Eerror.mag());
			}
		}
		fclose(out);
	}

	printf("Test 4: timing for various grid sizes, 100k particles\n");
	sigma = 0.2;
	static int NN[4] = {33, 65, 129, 257};
	FILE *summary = fopen("summary.txt","w");
	assert(summary != 0);
	fprintf(summary,"#N tolerance nCycles time\n");
	for(int iNN=0; iNN<4; ++iNN) {
	    N = NN[iNN];
	    static float tol[4]={1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7};
	    for(int itol=0; itol<4; ++itol) {
		double tolerance = tol[itol];
		PoissonConvolve3D poisson(N,N,N,-1.0,1.0,-1.0,1.0,-1.0,1.0);
		printf("N=%d tolerance=%.1e\n",N,tolerance);
		start=gettimeus();
		poisson.zeroRhs();
		poisson.setPutMode(PoissonConvolve3D::PRORATED);
		totalCharge=0.0;
	        for(int ix=0; ix<N; ++ix) {
		    for(int iy=0; iy<N; ++iy) {
			for(int iz=0; iz<N; ++iz) {
				G4ThreeVector p=poisson.getPosition(ix,iy,iz);
				double r=p.mag();
				double c=exp(-r*r/sigma/sigma/2.0);
				if(poisson.putCharge(p.x(),p.y(),p.z(),c))
					totalCharge += c;
			}
		    }
	        }
		poisson.defaultApproximation();
		printf("Setup took %lld us\n",gettimeus()-start);
		start = gettimeus();
		poisson.solve();
	    	printf("solve took %lld us;  N=%d\n", gettimeus()-start,N);
		if(true) {
			fprintf(summary,"%d\t%.1e\t%d\t%.6f\n",N,tolerance,
				1,(double)start/1.0e6);
			char fname[64];
			sprintf(fname,"N=%d,tol=%.1e.txt",NN[iNN],tol[itol]);
			FILE *out = fopen(fname,"w");
			assert(out != 0);
			fprintf(out,"#x y z r PHImodel PHIexact Emodel Eexact Eerror\n");
			double dx = 1.0/(N-1)/3.0;
			double sqrt2_sigma=sqrt(2.0)*sigma;
			for(double r=dx; r<3.0; r+=dx) {
				double exact = totalCharge/r * erf(r/sqrt2_sigma);
				for(int m=0; m<5; ++m) {
					double theta=M_PI*drand48();
					double phi=2.0*M_PI*drand48();
					double x=r*sin(theta)*cos(phi);
					double y=r*sin(theta)*sin(phi);
					double z=r*cos(theta);
					double v = (double)poisson.phi(x,y,z);
					G4ThreeVector Emodel=poisson.getE(x,y,z);
					G4ThreeVector Eexact=gaussianE(x,y,z,
							totalCharge,sigma);
					G4ThreeVector Eerror=Emodel;
					Eerror -= Eexact;
					fprintf(out,"%.4f %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f\n",
						x,y,z,r,v,exact,
						Emodel.mag(),Eexact.mag(),Eerror.mag());
				}
			}
			fclose(out);
		}
	    }
	}
	fclose(summary);
	
	exit(0);
}

#endif // MAIN

#else // G4BL_FFTW
int dummy_poissonconvolve3d=0;
#endif // G4BL_FFTW
