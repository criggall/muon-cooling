//	PoissonConvolve3D.hh
/*
	NOTE: fftw must be installed with float libraries.
*/

#ifndef POISSONCONVOLVE3D_HH
#define POISSONCONVOLVE3D_HH

#include <assert.h>
#include <math.h>
#include <vector>

#include "fftw3.h"

#ifdef USE_DOUBLE
#define float double
#define sqrtf sqrt
#define fftwf_complex fftw_complex
#define fftwf_plan fftw_plan
#define fftwf_malloc fftw_malloc
#define fftwf_free fftw_free
#define fftwf_destroy_plan fftw_destroy_plan
#define fftwf_set_timelimit fftw_set_timelimit
#define fftwf_plan_dft_r2c_3d fftw_plan_dft_r2c_3d
#define fftwf_plan_dft_c2r_3d fftw_plan_dft_c2r_3d
#define fftwf_execute fftw_execute
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifdef MAIN
///@ testing - don't require Geant4 headers.
struct G4ThreeVector { 
	double xx, yy, zz; 
	G4ThreeVector(double x, double y, double z) { xx=x; yy=y; zz=z; }
	double x() { return xx; }
	double y() { return yy; }
	double z() { return zz; }
	double mag() { return sqrt(xx*xx+yy*yy+zz*zz); }
	G4ThreeVector& operator += (G4ThreeVector& r)
		{ xx += r.xx; yy += r.yy; zz += r.zz; return *this; }
	G4ThreeVector& operator -= (G4ThreeVector& r)
		{ xx -= r.xx; yy -= r.yy; zz -= r.zz; return *this; }
	G4ThreeVector& operator *= (double v)
		{ xx *= v; yy*= v; xx*= v; return *this; }
	G4ThreeVector& operator /= (double v)
		{ xx /= v; yy/= v; xx/= v; return *this; }
};
const double megavolt=1.0, meter=1.0;
#else
#include "G4ThreeVector.hh"
#endif //MAIN

/// Charge represents one point charge for the approximation outside the grid.
struct Charge {
	float x, y, z, c;
	Charge(float _x, float _y, float _z, float _c) 
		{ x=_x; y=_y; z=_z; c=_c; }
	G4ThreeVector getPosition() { return G4ThreeVector(x,y,z); }
	float getCharge() { return c; }
};

/**	class PoissonConvolve3D
 *
 *	A solver for Poisson's equation in 3-D Cartesian coordinates {x,y,z},
 *	using convolution of the Green's function (whch has infinite boundary 
 *	conditions). The convolution is sped up by using an FFT on a grid
 *	twice as large (in each dimension) as the solution (physical) grid.
 *
 *	This class uses the Poisson solver on a rectangular 3-d grid; inside
 *	the grid, it calculates E by applying the analytic derivative of the
 *	interpolation function to the values of phi at the surrounding 8
 *	grid points. Outside the grid, by default it uses the analytic value
 *	of E for a single point charge, with its charge equal to the total
 *	charge in the grid, and its position at the mean position of the
 *	charge in the grid. More complicated approximations, based on
 *	multiple point charges, can be used (not discussed here).
 *
 *	Units are determined by the code using this class. Normally the
 *	charge values will not be integers, but rather integers times the
 *	positron charge divided by epsilon_0 (the permittivity of free space),
 *	the grid sizes will have units of length; the result is that
 *	phi() and getE() will then have the appropriate units.
 *
 *	This class solves:
 *		del^2 phi = -4*pi*rho
 *	On a specified grid in {x,y,z}. rho(x,y,z) is the charge density.
 *	The -4*pi and the conversion from charge to charge density is handled
 *	internally in putCharge(); they are NOT handled in setRhs() (which
 *	is intended primarily for testing).
 *
 *	This class includes an approximation for use outside the grid. The
 *	default is simply to use the total charge put into the grid (via
 *	putCharge()), and the weighted mean position of that charge. Better
 *	approximations can be computed externally as a vector of point
 *	charges and used via setApproximation().
 *
 *	Note: array[] contains CHARGE between calls to zeroRhs() and solve();
 *	it contains PHI between calls to solve() and zeroRhs(); status
 *	reflects its contents, ans assert()-s are used to ensure accesses
 *	are valid..
 *
 *
 *	Accuracy for a 65x65x65 grid of unit size:
 *	... MORE TO COME ...
 *
 *
 *	The CPU time increases with increasing N, and can vary by as much as a
 *	factor of 3 for adjacent values. The following are good values:
 *	128, 120, 112, 100, 96, 90, 84, 80, 75, 70, 64, 60, 56, 48, 42, 36, 32.
 *	Below 28 the variations are small. The following are particularly SLOW
 *	values: 129, 127, 123, 107, 103, 101.
 **/
class PoissonConvolve3D {
public:
	/// enum PutMode determines how putCharge() places charge into
	/// the grid. PRORATED is the default.
	/// LOWEST places charge into the grid point lower than (x,y,z).
	/// NEAREST places charge into the nearest grid point to (x,y,z).
	/// PRORATED places charge into the 8 nearest grid points, pro-rated
	/// by proximity to (x,y,z).
	/// For 1M particles into a 129^3 grid, NEAREST and LOWER were less
	/// than 10% faster than PRORATED, giving slightly more ragged results.
	enum PutMode {LOWEST, NEAREST, PRORATED};
protected:
	enum Status {UNINITIALIZED, CHARGE, PHI};
	int nx, ny, nz;
	float xMin, xMax, yMin, yMax, zMin, zMax;
	Status status;
	float dx, dy, dz;	// grid spacings
	float *array;		// phi[index(ix,iy,iz)] or rho[index(ix,iy,iz)]
	fftwf_complex *transf;	// transform of array (working)
	fftwf_complex *greensFunction;
	std::vector<Charge> approx;
	float sumX, sumY, sumZ, totalCharge;	// for defaultApproximation()
	float vol;
	PutMode putMode;
	int nc;
	fftwf_plan gfPlan;
	fftwf_plan plan1;
	fftwf_plan plan2;
	int index(int ix, int iy, int iz) const {
		assert(ix>=0&&ix<2*nx && iy>=0&&iy<2*ny && iz>=0&&iz<2*nz);
		return iz+2*nz*(iy+2*ny*ix);  // C order
	}
public:
	/// Constructor.
	/// nx,ny,nz are the # grid points in the 3 dimensions -- NOTE:
	/// these are the size of the PHYSICAL grid, the actual grid is
	/// doubled in each dimension. They can be any integers >= 9, but
	/// powers of 2 have a performance advantage.
	/// Note that 512x512x512 allocates too much memory for a 32-bit
	/// program on an Intel processor.
	/// xMin,xMax,yMin,yMax,zMin,zMax are the grid limits in 3 dimensions.
	/// That implies array[index(0,iy,iz)] is the boundary at x=xMin; 
	/// array[index(nx-1,iy,iz)] is the boundary at x=xMax; etc.
	/// 1.0e-4 is appropriate.
	PoissonConvolve3D(int _nx, int _ny, int _nz, float _xMin, float _xMax,
			float _yMin, float _yMax, float _zMin, float _zMax);

	/// Destructor.
	virtual ~PoissonConvolve3D() {
		if(array) fftwf_free(array);
		if(transf) fftwf_free(transf);
		if(greensFunction) fftwf_free(greensFunction);
		fftwf_destroy_plan(gfPlan);
		fftwf_destroy_plan(plan1);
		fftwf_destroy_plan(plan2);
		approx.clear();
	}

	/// updatePosition() updates the position of the boundaries.
	/// Usually followed by a call to setApproximation(), with the
	/// approximation charge positions updated similarly, or by
	/// a series of calls to putCharge() and then defaultApproximation().
	void updatePosition(float _xMin, float _xMax,
			float _yMin, float _yMax, float _zMin, float _zMax) {
		xMin=_xMin; xMax=_xMax;
		yMin=_yMin; yMax=_yMax;
		zMin=_zMin; zMax=_zMax;
		dx = (xMax-xMin)/(nx-1);
		dy = (yMax-yMin)/(ny-1);
		dz = (zMax-zMin)/(nz-1);
		vol = dx*dy*dz;
		computeGreensFunction();
	}

	/// phi() works for all {x,y,z}, using the solution inside the grid
	/// and the approximation outside the grid.
	/// Note that phi() is not necessarily continuous or differentiable,
	/// but it should be approximately continuous.
	/// Do not differentiate this function, use getE() instead (it is
	/// MUCH more intelligent about the grid and its bounadries).
	/// Intended primarily for testing.
	float phi(float x, float y, float z);

	/// phiGrid() returns the value of the solution to Poisson's equation
	/// (i.e. the electrostatic potential) on grid points.
	/// ix,iy,iz MUST be in range: 0..nx-1, 0..ny-1, 0..nz-1
	/// Intended primarily for testing.
	float phiGrid(int ix, int iy, int iz) const
		{ assert(status==PHI);  return array[index(ix,iy,iz)]; }

	/// setPhi() will set a value in array[].
	/// Intended primarily for testing.
	void setPhi(int ix, int iy, int iz, float v)
		{ assert(status==PHI);  array[index(ix,iy,iz)] = v; }

	/// addPhi() will add v to a value in array[].
	/// Intended primarily for testing.
	void addPhi(int ix, int iy, int iz, float v)
		{ assert(status==PHI);  array[index(ix,iy,iz)] += v; }

	/// rhs() returns the charge at a grid point of the right-hand-side 
	/// of Poisson's equation.
	/// NOTE: ix,iy,iz MUST be in range: 0..nx-1, 0..ny-1, 0..nz-1
	float rhs(int ix, int iy, int iz) const
		{ assert(status==CHARGE);
		  return 4.0*M_PI*array[index(ix,iy,iz)]; }

	/// setRhs() will set an entry in rhs[]. Normally putCharge() is better.
	/// ix,iy,iz MUST be in range: 0..nx-1, 0..ny-1, 0..nz-1
	/// This directly sets the rhs array entry, so you probably need
	/// to multiply by -4*pi/(dx*dy*dz) to get the proper R.H.S. for
	/// Poisson's equation; putCharge() does that internally.
	/// Intended primarily for testing.
	void setRhs(int ix, int iy, int iz, float v)
		{ assert(status==CHARGE); array[index(ix,iy,iz)] = v; }

	/// putCharge() puts a charge into the rhs grid according to
	/// the current putMode. Does nothing if outside the grid.
	/// Multiplies by -1.0/(dx*dy*dz) to convert the charge to the R.H.S.
	/// of Poisson's equation, which is minus the charge density.
	/// NOTE: dividing by epsilon0 is up to the user's choice of units.
	/// Accumulates the total charge and mean position of charge actually
	/// placed into the grid (for use by defaultApproximation()).
	/// Returns true if the charge was placed within the grid, false if not.
	bool putCharge(float x, float y, float z, float c);
	bool putCharge(Charge &c) { return putCharge(c.x,c.y,c.z,c.c); }

	/// setPutMode() sets the mode of putCharge(). See PutMode above.
	void setPutMode(PutMode m) { putMode = m; }

	/// zeroRhs() will zero rhs.
	void zeroRhs();

	/// getE() returns the value of the E field at a given point.
	/// Uses the solution inside the grid and the approximation outside.
	/// It handles grid cells and boundaries intelligently.
	G4ThreeVector getE(float x, float y, float z);

	/// approxE() evaluates the approximation E field.
	/// {x,y,z} can be anywhere, but normally are outside the grid or on
	/// its boundary.
	G4ThreeVector approxE(double x, double y, double z);

	/// approxPhi() evaluates the approximation phi.
	/// {x,y,z} can be anywhere, but normally are outside the grid or on
	/// its boundary.
	float approxPhi(double x, double y, double z);

	/// getPosition() returns the position corresponding to grid point
	/// indexes. Intended primarily for testing.
	G4ThreeVector getPosition(int ix, int iy, int iz)
		{ return G4ThreeVector(xMin+ix*dx,yMin+iy*dy,zMin+iz*dz); }

	/// setApproximation() sets the approximation to use outside the
	/// grid, and on its boundary. v is a vector of point charges.
	/// Updates the boundary values of phi.
	void setApproximation(std::vector<Charge> v) 
		{ approx = v; }

	/// defaultApproximation() calls setApproximation() with the
	/// total charge and mean position in rhs (from putCharge()).
	/// Updates the boundary values of phi, and the interior if
	/// interior=true.
	void defaultApproximation();

	/// computeGreensFunction() does just that. Must be called whenever
	/// the physical size of the grid is changed.
	void computeGreensFunction();

	/// solve() will solve Poisson's equation 
	void solve();

	/// printRhs() prints the rhs.
	/// Intended primarily for testing.
	void printRhs();

	/// printPhi() prints phi.
	/// Intended primarily for testing.
	void printPhi();

	/// accessors.
	int getNx() const { return nx; }
	int getNy() const { return ny; }
	int getNz() const { return nz; }
};

#endif // POISSONCONVOLVE3D_HH
