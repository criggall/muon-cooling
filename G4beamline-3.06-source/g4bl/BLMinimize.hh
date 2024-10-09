//	BLMinimize.hh
#define USE_GSL_MULTIMIN
//#define USE_TMINUIT
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

#ifndef BLMINIMIZE_HH
#define BLMINIMIZE_HH

#ifdef G4BL_GSL

#include <stdio.h>
#include <vector>
#include <string>

#include "globals.hh"
#include "BLAssert.hh"

/**	class BLMinimizeFunction is the base class for a class defining a
 *	function to minimize.
 **/
class BLMinimizeFunction {
public:
	BLMinimizeFunction() { }
	virtual ~BLMinimizeFunction() { }

	/// operator() is the function to minimize. Note it is not const
	/// for flexibility in derived classes.
	/// Its argument is necessarily const.
	virtual double operator()(const std::vector<double> &x) = 0;
};

#ifdef USE_GSL_MULTIMIN
#ifndef G4BL_GSL
#error GSL not available
#endif

#include <gsl/gsl_multimin.h>

/**	class BLMinimize implements a general function minimizer using
 *	the GSL multimin functions.
 *
 *	Uses the the Simplex algorithm of Nelder and Mead.
 *
 *	Compilation for simple programs:
 *		g++ -I... test.cc -lgsl -lgslcblas -lm
 *
 *	NOTE: the parameters to the minimizer must all have scales that are
 *	comparable, or the minimizer will often converge to a wildly incorrect
 *	value. This applies to track fitting, where positions have scales of
 *	1 (millimeters), but angles have scales of 0.001 (milliradians).
 *	This class scales the parameters appropriately. it also permits the
 *	function to have arguments that are held constant during the 
 *	minimization by setting their scale to 0.0.
 *
 *	NOTE: only a single instance of BLMinimize can be active at a time
 *	(i.e. inside the minimize() function). Any number of instances can
 *	exist. This basically means your BLMinimizeFunction cannot itself
 *	use a minimizer internally.
 **/
class BLMinimize {
	std::vector<double> scale;
	std::vector<double> arg;
	std::vector<std::string> name;
	double val;
	double size;
	int iter;
	int print;
	double chisqMin;
	BLMinimizeFunction *func;
	gsl_multimin_function ff;
	gsl_multimin_fminimizer *m;
private:
	// f() is the function the gsl routines call (they are in C)
	static double f(const gsl_vector *x, void *param) {
		BLMinimize *p = (BLMinimize*)param;
		return p->call(x);
	}
	// call() will call the user's function, scaling the argument vector.
	double call(const gsl_vector *x) {
		if(print >= 2) printf("BLMinimize call args:");
		int j=0;
		for(unsigned i=0; i<scale.size(); ++i) {
			if(scale[i] != 0.0) // fixed args are already set
				arg[i] = gsl_vector_get(x,j++) * scale[i];
			if(print >= 2)
				printf(" %.4g",arg[i]);
		}
		if(print >= 2) printf("\n");
		double v = (*func)(arg);
		if(print >= 2) printf("BLMinimize call returns value=%.4g\n",v);
		return v;
	}
public:
	/// Constructor.
	BLMinimize(int _print=2) : scale(), arg(), name()
		{ val=GSL_NAN; size=GSL_NAN; iter=0; print=_print;
		  chisqMin=GSL_NAN; func=0; ff.f=0; ff.n=0; ff.params=0;  m=0;
		}

	virtual ~BLMinimize() { if(m) gsl_multimin_fminimizer_free(m); }

	/// setFunc() will set the function to minimize and the scale of its
	/// arguments. The number of non-zero entries in _scale[] determines
	/// the number of variables to minimize. names is a colon-separated
	/// list of parameter names.
	void setFunc(BLMinimizeFunction &f, const std::vector<double> &_scale,
						const std::string &names="") {
		func = &f;
		scale.clear();
		arg.clear();
		name.clear();
		std::vector<std::string> n = splitString(names);
		int nx = 0;
		for(unsigned i=0; i<_scale.size(); ++i) {
			double s=_scale[i];
			scale.push_back(s);
			arg.push_back(s);
			name.push_back(i<n.size() ? n[i] : "");
			if(s != 0.0) ++nx;
		}
		ff.f=BLMinimize::f;
		ff.n=nx;
		ff.params=this;
		if(m) gsl_multimin_fminimizer_free(m);
		m = 0;
		if(nx == 0) return;
		m = gsl_multimin_fminimizer_alloc(
		  	gsl_multimin_fminimizer_nmsimplex,nx);
	}

	/// setPrintDetail() sets the detail of printing per iteration:
	/// 0=no print, 1=print value and size only, 2=print value, size, and x.
	/// the default is 2.
	void setPrintDetail(int v=0) { print = v; }

	/// setChisqMin() sets an additional convergence criterium of chisq
	/// less than the given value. Can only be used when the function is
	/// positive definite, usually a chi-squared. Required for track-
	/// fitting of test tracks that have truly 0 chisq -- the normal 
	/// convergence often fails for such tracks.
	void setChisqMin(double v) { chisqMin = v; }

	/// minimize() finds the minimum, returning status (0=success).
	///   x		the initial and then the final argument values
	///   tolerance	the convergence limit on the simplex size
	///		NOTE: this is after scaling!
	///   maxIter	limits the number of iterations
	/// The initial step for x[i] is scale[i] (x[i] is scaled, so the
	/// step to the gsl routine is 1.0).
	int minimize(std::vector<double> &x, double tolerance=1e-3,
							int maxIter=100) {
		if(ff.n <= 0) return GSL_EINVAL;
		BLAssert(x.size()==scale.size());
		gsl_vector *init = gsl_vector_alloc(ff.n);
		gsl_vector *istep = gsl_vector_alloc(ff.n);
		unsigned j=0;
		for(unsigned i=0; i<scale.size(); ++i) {
			arg[i] = x[i]; // fixed args, variable ones overwritten
			if(scale[i] != 0.0) {
				gsl_vector_set(init,j,x[i]/scale[i]); 
				gsl_vector_set(istep,j,1.0);
				++j;
			}
		}
		BLAssert(j==ff.n); BLAssert(m!=0);
		gsl_multimin_fminimizer_set(m,&ff,init,istep);
		for(iter=1; iter<=maxIter; ++iter) {
			int status = gsl_multimin_fminimizer_iterate(m);
			if(status) {
				if(print >= 1)
					printf("BLMinimize FAILED status=%d\n",
									status);
				gsl_vector_free(init);
				gsl_vector_free(istep);
				return status;
			}
			val = m->fval;
			size = gsl_multimin_fminimizer_size(m);
			if(print >= 1)
			    printf("BLMinimize Iter %d: val=%.4g size=%.4g\n",
			    				iter,val,size);
			status = gsl_multimin_test_size(size,tolerance);
			if(status == GSL_SUCCESS || (!gsl_isnan(chisqMin) &&
			   val <= chisqMin)) {
				unsigned j=0;
				for(unsigned i=0; i<scale.size(); ++i) {
				    if(scale[i] != 0.0)
					x[i] = gsl_vector_get(m->x,j++) *
								scale[i];
				}
				gsl_vector_free(init);
				gsl_vector_free(istep);
				if(print >= 1) printf("BLMinimize SUCCESS\n");
				return 0;
			}
		}
		gsl_vector_free(init);
		gsl_vector_free(istep);
		if(print >= 1) printf("BLMinimize FAILED iter > %d\n",maxIter);
		return GSL_EMAXITER;
	}

	/// value() returns the function value from the previous minimize().
	double value() const { return val; }

	/// getSize() returns the Simplex size from the previous minimize().
	double getSize() { return size; }

	/// getIterations() returns the # iterations of the previous minimize().
	int getIterations() { return iter; }

private:
	std::vector<std::string> splitString(const std::string &s) {
		std::string tmp(s);
		std::vector<std::string> v;
		while(tmp.size() > 0) {
			unsigned i=tmp.find_first_of(":");
			if(i > tmp.size()) i = tmp.size();
			v.push_back(tmp.substr(0,i));
			tmp.erase(0,i+1);
		}
		return v;
	}
};

#endif // USE_GSL_MULTIMIN

#ifdef USE_TMINUIT

#include "TMinuit.h"

/**	class BLMinimize implements a general function minimizer using
 *	the TMinuit class from Root.
 *
 *	NOTE: only a single instance of BLMinimize can be active at a time
 *	(i.e. inside the minimize() function). Any number of instances can
 *	exist. This basically means your BLMinimizeFunction cannot itself
 *	use a minimizer internally.
 **/
class BLMinimize {
	std::vector<double> scale;
	std::vector<std::string> name;
	std::vector<double> arg;
	double val;
	int iter;
	int print;
	double chisqMin;
	BLMinimizeFunction *func;
	TMinuit *tm;
	static BLMinimize *current;
private:
	static void fcn(Int_t &npar, Double_t *grad, Double_t &val, Double_t *x,
								Int_t iflag) {
		BLAssert(current != 0);
		val = current->call(npar,x);
	}
	double call(int npar, double *x) {
		// npar seems to be 1 too small assert(npar == (int)arg.size());
		for(unsigned i=0; i<arg.size(); ++i)
			arg[i] = x[i];
		return (*func)(arg);
	}
public:
	/// Constructor.
	BLMinimize(int _print=2) : scale(), name(), arg()
		{ val=0.0; iter=0; print=_print; chisqMin=0.0; func=0; tm=0; }

	virtual ~BLMinimize() {  }

	/// setFunc() will set the function to minimize and the scale of its
	/// arguments. The number of non-zero entries in _scale[] determines
	/// the number of variables to minimize. names is a colon-separated
	/// list of parameter names.
	void setFunc(BLMinimizeFunction &f, const std::vector<double> &_scale,
						const std::string &names="") {
		func = &f;
		scale.clear();
		name.clear();
		arg.clear();
		std::vector<std::string> n = splitString(names);
		for(unsigned i=0; i<_scale.size(); ++i) {
			double s=_scale[i];
			scale.push_back(s);
			name.push_back(i<n.size() ? n[i] : "");
			arg.push_back(0.0);
		}
		
		if(tm) delete tm;
		tm = new TMinuit();
		tm->SetFCN(fcn);
		tm->SetPrintLevel(-1);
		tm->Command("SET PRINTOUT -1");
	}

	/// setPrintDetail() sets the detail of printing per iteration:
	/// 0=no print, 1=print value and size only, 2=print value, size, and x.
	/// the default is 2.
	void setPrintDetail(int v=0) { print = v; }

	/// setChisqMin() sets an additional convergence criterium of chisq
	/// less than the given value. Can only be used when the function is
	/// positive definite, usually a chi-squared. Required for track-
	/// fitting of test tracks that have truly 0 chisq -- the normal 
	/// convergence often fails for such tracks.
	void setChisqMin(double v) { chisqMin = v; }

	/// minimize() finds the minimum, returning status (0=success).
	///   x		the initial and then the final argument values
	///   tolerance	the convergence limit on the simplex size
	///		NOTE: this is after scaling!
	///   maxIter	limits the number of iterations
	/// The initial step for x[i] is scale[i].
	int minimize(std::vector<double> &x, double tolerance=1e-3,
							int maxIter=100) {
		BLAssert(tm != 0);
		BLAssert(x.size() == scale.size());
		for(unsigned i=0; i<scale.size(); ++i) {
			tm->DefineParameter(i,name[i].c_str(),x[i],scale[i],0.0,0.0);
		}
		//@ What do I do with tolerance????
		tm->SetMaxIterations(maxIter);
		current = this;
		int retval = tm->Migrad();
		current = 0;
		for(unsigned i=0; i<scale.size(); ++i) {
			double err;
			tm->GetParameter(i,x[i],err);
		}
		val = (*func)(x);
		return retval;
	}

	/// value() returns the function value from the previous minimize().
	double value() const { return val; }

	/// getSize() returns the Simplex size from the previous minimize().
	double getSize() { return 0.0; }

	/// getIterations() returns the # iterations of the previous minimize().
	int getIterations() { return iter; }

private:
	std::vector<std::string> splitString(const std::string &s) {
		std::string tmp(s);
		std::vector<std::string> v;
		while(tmp.size() > 0) {
			unsigned i=tmp.find_first_of(":");
			if(i > tmp.size()) i = tmp.size();
			v.push_back(tmp.substr(0,i));
			tmp.erase(0,i+1);
		}
		return v;
	}
};


#endif // USE_TMINUIT

#endif // G4BL_GSL
#endif // BLMINIMIZE_HH
