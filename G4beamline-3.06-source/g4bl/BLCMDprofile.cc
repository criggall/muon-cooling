//	BLCMDprofile.cc
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

#include <vector>
#include <stdio.h>

#include <G4ParticleDefinition.hh>

#include "BLManager.hh"
#include "BLParam.hh"
#include "BLEvaluator.hh"
#include "BLCoordinates.hh"
#include "BLWriteAsciiFile.hh"
#include "BLMPI.hh"

static double determinant(double **a,int n);	// at end of this file


/**	class BLCMDprofile will print profile information at a list of Z
 *	positions.
 **/
class BLCMDprofile : public BLCommand, public BLCallback {
	struct Entry : public BLManager::ZSteppingAction {
		G4double z;
		G4String require;
		// data for sumAllWorkers: 16+7+1=24 doubles
		G4double Mdata[16];
		G4double sumX, sumY, sumPx, sumPy, sumPz, sumPtot, sumW;
		G4double n;
		// end of data for sumAllWorkers
		BLEvaluator eval;
		G4double **M;
		G4ParticleDefinition *particle;
		BLCoordinateType coordinateType;
		Entry(G4double zpos, G4String& req, G4ParticleDefinition *pd,
				BLCoordinateType _coordinateType) :
							require(req), eval()
			{ z=zpos; sumX=sumY=sumPx=sumPy=sumPz=sumPtot=sumW=0.0;
			  n=0;
			  M=new double*[4]; 
			  M[0]=Mdata+0; M[1]=Mdata+4; M[2]=Mdata+8;
			  M[3]=Mdata+12;
			  for(int i=0; i<4; ++i) {
			  	for(int j=0; j<4; ++j)
			 		M[i][j]=0.0;
			  }
for(int k=0; k<12; ++k) BLAssert(Mdata[k] == 0.0);
			  particle = pd;
			  coordinateType = _coordinateType;
			}
		void UserZSteppingAction(const G4Track *track);
	};
	G4String z;
	G4String zloop;
	G4String require;
	G4String file;
	G4String particle;
	G4String coordinates;
	BLCoordinateType coordinateType;
	std::vector<Entry*> entries;
	unsigned place;
	G4ParticleDefinition *pd;
public:
	/// Constructor.
	BLCMDprofile();
	BLCMDprofile(BLCMDprofile& r);

	/// commandName() returns "profile".
	virtual G4String commandName() { return "profile"; }
	
	/// command() implements the profile command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// callback from BLCallback.
	virtual void callback(int type);
};

BLCMDprofile defaultProfile;

BLCMDprofile::BLCMDprofile() : BLCommand(), BLCallback(), entries()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("write beam profile information to a file");
	setDescription("This command accumulates the moments of the track\n"
		"distributions during the run, and at the end of\n"
		"run prints the mean, sigma, emittance (RMS), alpha, and\n"
		"beta (Twiss parameters) for the tracks. Each z\n"
		"position generates a line in the output file.\n\n"
		"Each value in z and zloop can be an expression using "
		"double constants and the usual C operators and "
		"functions.\n\n"
		"NOTE: This computation has two rather serious deficiencies:\n"
		"  1) it includes all tracks, with no sigma cut, so a single\n"
		"     outlier track can have a large effect on sigmas and\n"
		"     emittances.\n"
		"  2) its emittance does not include the vector potential.\n"
		"To avoid these deficiencies, it may be better to write "
		"NTuples and apply an external program such as ecalc9.\n\n"
		"This command is not placed into the geometry.");

	z = "";
	zloop = "";
	require = "";
	file = "-";
	particle = "mu+";
	coordinates = "Centerline";
	coordinateType = BLCOORD_CENTERLINE;
	place = 0;
	pd = 0;
}

BLCMDprofile::BLCMDprofile(BLCMDprofile& r) : BLCommand(r), 
		BLCallback(), entries()
{
	z = r.z;
	zloop = r.zloop;
	require = r.require;
	file = r.file;
	particle = r.particle;
	coordinates = r.coordinates;
	coordinateType = r.coordinateType;
	place = r.place;
	pd = 0;
}

int BLCMDprofile::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDprofile *p = new BLCMDprofile(defaultProfile);

	int retval = p->handleNamedArgs(namedArgs);

	p->coordinateType = BLCoordinates::getCoordinateType(p->coordinates);

	p->pd =  G4ParticleTable::GetParticleTable()->FindParticle(p->particle);
	if(!p->pd) printError("profile: unknown particle '%s'\n",p->particle.c_str());

	BLEvaluator eval;

	// register the z position(s) in the string z
	std::vector<G4double> vect = getList(p->z,",");
	if(vect.size() > 0) {
		for(unsigned i=0; i<vect.size(); ++i) {
			Entry *e = new Entry(vect[i],p->require,p->pd,
							p->coordinateType);
			BLManager::getObject()->registerZStep(vect[i],e,4);
			p->entries.push_back(e);
		}
	} else if(p->z.size() > 0) {
		printError("printf: Syntax error in z");
	}

	// register the z position(s) in the string zloop
	vect = getList(p->zloop,",:");
	if(vect.size() == 3) {
		if(vect[2] < 1.0*mm || vect[0] > vect[1]) {
			printError("printf: invalid zloop '%s'",
							p->zloop.c_str());
		} else {
			while(vect[0] <= vect[1]) {
			    Entry *e = new Entry(vect[0],p->require,p->pd,
							p->coordinateType);
			    BLManager::getObject()->registerZStep(vect[0],
								e,4);
			    p->entries.push_back(e);
			    vect[0] += vect[2];
			}
		}
	} else if(p->zloop.size() > 0) {
		printError("printf: invalid zloop '%s'",p->zloop.c_str());
	}

	BLManager::getObject()->registerCallback(p,2);

	p->print("");

	return retval;
}

void BLCMDprofile::defineNamedArgs()
{
	argString(z,"z","Comma-separated list of Z positions for profiling (mm)");
	argString(zloop,"zloop","Loop in z, first:last:incr (mm)");
	argString(require,"require","logical expression for cutting (default=true)");
	argString(particle,"particle","Name of particle to profile (default=mu+)");
	argString(file,"file","Output filename (default=stdout)");
	argString(file,"filename","Synonym for file");
	argString(coordinates,"coordinates","Coordinates: centerline or reference (default=c).");
}
void BLCMDprofile::Entry::UserZSteppingAction(const G4Track *track)
{
	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	if(track->GetDefinition() != particle) {
		return;
	}

	if(require.size() != 0) {
		eval.setTrackVariables(track,coordinateType); 
		double v = eval.evaluate(require);
		if(!eval.isOK()) {
		    BLCommand::printError("profile: invalid expression '%s'\n",
					require.c_str());
		    require = "0";
		    return;
		}
		if(v == 0.0)
			return;
	}

	G4double pos[4];
	G4ThreeVector momentum;
	G4double w=track->GetWeight();
	BLCoordinates *coord = (BLCoordinates*)track->GetUserInformation();
	if(!coord || !coord->isValid()) return;
	momentum = track->GetMomentum();
	coord->getCoords(coordinateType,pos);
	momentum = coord->getRotation() * momentum;
	sumW += w;
	sumX += w*pos[0];
	sumY += w*pos[1];
	sumPx += w*momentum[0];
	sumPy += w*momentum[1];
	sumPz += w*momentum[2];
	sumPtot += w*momentum.mag();
	M[0][0] += w*pos[0]*pos[0];
	M[0][1] += w*pos[0]*momentum[0];
	M[0][2] += w*pos[0]*pos[1];
	M[0][3] += w*pos[0]*momentum[1];
	M[1][0] += w*momentum[0]*pos[0];
	M[1][1] += w*momentum[0]*momentum[0];
	M[1][2] += w*momentum[0]*pos[1];
	M[1][3] += w*momentum[0]*momentum[1];
	M[2][0] += w*pos[1]*pos[0];
	M[2][1] += w*pos[1]*momentum[0];
	M[2][2] += w*pos[1]*pos[1];
	M[2][3] += w*pos[1]*momentum[1];
	M[3][0] += w*momentum[1]*pos[0];
	M[3][1] += w*momentum[1]*momentum[0];
	M[3][2] += w*momentum[1]*pos[1];
	M[3][3] += w*momentum[1]*momentum[1];
	++n;
}

void BLCMDprofile::callback(int type)
{
	// after beam tracking is complete, print the beam profile data
	if(type != 2) return;

	bool printing = (!BLMPI::isMPI() || BLMPI::isRank0());
	FILE *fd = 0;
	if(printing) {
		fd = stdout;
		if(file != "-") {
			fd = BLWriteAsciiFile::fopen(file);
			if(!fd) {
				G4Exception("profile",
					"Cannot write file - stdout used",
							JustWarning,file);
				fd = stdout;
			}
		}
		fprintf(fd,"# g4beamline profile\n");
		fprintf(fd,"#Z N meanX sigmaX meanY sigmaY emitX emitY emitTrans betaX betaY betaTrans alphaX alphaY alphaTrans meanP\n");
		fprintf(fd,"#mm - mm mm mm mm mm-rad mm-rad mm-rad mm mm mm - - - MeV/c\n");
	}
	
	// one message per entry is too many messages and overwhelms rank 0
	// so combine all entries into one call to sumAllWorkers() 
	if(BLMPI::isMPI()) {
		int ndata = entries.size()*24;
		double *data = new double[ndata];
		BLAssert(data != 0);
		double *d=data;
		for(unsigned i=0; i<entries.size(); ++i) {
			Entry *e = entries[i];
			for(int j=0; j<24; ++j)
				*d++ = e->Mdata[j];
		}
		BLMPI::sumAllWorkers(data,ndata);
		if(!printing) {
			return;
		}
		d=data;
		for(unsigned i=0; i<entries.size(); ++i) {
			Entry *e = entries[i];
			for(int j=0; j<24; ++j) {
				e->Mdata[j] = *d++;
			}
		}
	}

	for(unsigned i=0; i<entries.size(); ++i) {
		int n = entries[i]->n;
		G4double sumW = entries[i]->sumW;
		if(n <= 0 || fabs(sumW) < 1.0e-12) continue;
		G4double meanX = entries[i]->sumX/sumW;
		G4double meanY = entries[i]->sumY/sumW;
		G4double meanPx = entries[i]->sumPx/sumW;
		G4double meanPy = entries[i]->sumPy/sumW;
		G4double meanPz = entries[i]->sumPz/sumW;
		G4double meanPtot = entries[i]->sumPtot/sumW;
		G4double **M = entries[i]->M;
		M[0][0] = M[0][0]/sumW - meanX*meanX;
		M[0][1] = M[0][1]/sumW - meanX*meanPx;
		M[0][2] = M[0][2]/sumW - meanX*meanY;
		M[0][3] = M[0][3]/sumW - meanX*meanPy;
		M[1][0] = M[1][0]/sumW - meanPx*meanX;
		M[1][1] = M[1][1]/sumW - meanPx*meanPx;
		M[1][2] = M[1][2]/sumW - meanPx*meanY;
		M[1][3] = M[1][3]/sumW - meanPx*meanPy;
		M[2][0] = M[2][0]/sumW - meanY*meanX;
		M[2][1] = M[2][1]/sumW - meanY*meanPx;
		M[2][2] = M[2][2]/sumW - meanY*meanY;
		M[2][3] = M[2][3]/sumW - meanY*meanPy;
		M[3][0] = M[3][0]/sumW - meanPy*meanX;
		M[3][1] = M[3][1]/sumW - meanPy*meanPx;
		M[3][2] = M[3][2]/sumW - meanPy*meanY;
		M[3][3] = M[3][3]/sumW - meanPy*meanPy;
		G4double sigmaX = sqrt(M[0][0]);
		G4double sigmaY = sqrt(M[2][2]);
		G4double mass = pd->GetPDGMass();
		G4double emitX = sqrt(M[0][0]*M[1][1]-M[1][0]*M[0][1])/mass;
		G4double emitY = sqrt(M[2][2]*M[3][3]-M[3][2]*M[2][3])/mass;
		G4double emitTrans = sqrt(sqrt(determinant(M,4)))/mass;
		G4double alphaX = -M[0][1]/mass/emitTrans;
		G4double alphaY = -M[2][3]/mass/emitTrans;
		G4double alphaTrans = (alphaX+alphaY)/2.0;
		G4double betaX = M[0][0]*meanPz/mass/emitTrans;
		G4double betaY= M[2][2]*meanPz/mass/emitTrans;
		G4double betaTrans = (betaX+betaY)/2.0;
		fprintf(fd,"%.1f %d %.3f %.3f %.3f %.3f %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g\n",
			entries[i]->z,n,meanX,sigmaX,meanY,sigmaY,emitX,emitY,emitTrans,
			betaX,betaY,betaTrans,alphaX,alphaY,alphaTrans,meanPtot/MeV);
	}
	if(printing) {
		if(fd != stdout)
			BLWriteAsciiFile::fclose(fd);
		printf("profile: wrote '%s'\n",file.c_str());
	}
}

//==============================================================================
// Recursive definition of determinate using expansion by minors.
//
// Notes: 1) arguments:
//             a (double **) pointer to a pointer of an arbitrary square matrix
//             n (int) dimension of the square matrix
//
//        2) determinant is a recursive function, calling itself repeatedly
//           each time with a sub-matrix of the original till a terminal
//           2X2 matrix is achieved and a simple determinat can be computed.
//           As the recursion works backwards, cumulative determinants are
//           found till untimately, the final determinate is returned to the
//           initial function caller.
//
//        3) m is a matrix (4X4 in example)  and m13 is a minor of it.
//           A minor of m is a 3X3 in which a row and column of values
//           had been excluded.   Another minor of the submartix is also
//           possible etc.
//             m  a b c d   m13 . . . .
//                e f g h       e f . h     row 1 column 3 is elminated
//                i j k l       i j . l     creating a 3 X 3 sub martix
//                m n o p       m n . p
//
//        4) the following function finds the determinant of a matrix
//           by recursively minor-ing a row and column, each time reducing
//           the sub-matrix by one row/column.  When a 2X2 matrix is
//           obtained, the determinat is a simple calculation and the
//           process of unstacking previous recursive calls begins.
//
//                m n
//                o p  determinant = m*p - n*o
//
//        5) this function uses dynamic memory allocation on each call to
//           build a m X m matrix  this requires **  and * pointer variables
//           First memory allocation is ** and gets space for a list of other
//           pointers filled in by the second call to malloc.
//
//        6) C++ implements two dimensional arrays as an array of arrays
//           thus two dynamic malloc's are needed and have corresponsing
//           free() calles.
//
//        7) the final determinant value is the sum of sub determinants
//
//==============================================================================


static double determinant(double **a,int n)
{
    int i,j,j1,j2 ;                    // general loop and matrix subscripts
    double det = 0 ;                   // init determinant
    double **m = NULL ;                // pointer to pointers to implement 2d
                                       // square array

    if (n < 1)    {   }                // error condition, should never get here

    else if (n == 1) {                 // should not get here
        det = a[0][0] ;
        }

    else if (n == 2)  {                // basic 2X2 sub-matrix determinate
                                       // definition. When n==2, this ends the
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1] ;// the recursion series
        }


                                       // recursion continues, solve next sub-matrix
    else {                             // solve the next minor by building a
                                       // sub matrix
        det = 0 ;                      // initialize determinant of sub-matrix

                                           // for each column in sub-matrix
        for (j1 = 0 ; j1 < n ; j1++) {
                                           // get space for the pointer list
            m = (double **) malloc((n-1)* sizeof(double *)) ;

            for (i = 0 ; i < n-1 ; i++)
                m[i] = (double *) malloc((n-1)* sizeof(double)) ;
                       //     i[0][1][2][3]  first malloc
                       //  m -> +  +  +  +   space for 4 pointers
                       //       |  |  |  |          j  second malloc
                       //       |  |  |  +-> _ _ _ [0] pointers to
                       //       |  |  +----> _ _ _ [1] and memory for
                       //       |  +-------> _ a _ [2] 4 doubles
                       //       +----------> _ _ _ [3]
                       //
                       //                   a[1][2]
                      // build sub-matrix with minor elements excluded
            for (i = 1 ; i < n ; i++) {
                j2 = 0 ;               // start at first sum-matrix column position
                                       // loop to copy source matrix less one column
                for (j = 0 ; j < n ; j++) {
                    if (j == j1) continue ; // don't copy the minor column element

                    m[i-1][j2] = a[i][j] ;  // copy source element into new sub-matrix
                                            // i-1 because new sub-matrix is one row
                                            // (and column) smaller with excluded minors
                    j2++ ;                  // move to next sub-matrix column position
                    }
                }

            det += pow(-1.0,1.0 + j1 + 1.0) * a[0][j1] * determinant(m,n-1) ;
                                            // sum x raised to y power
                                            // recursively get determinant of next
                                            // sub-matrix which is now one
                                            // row & column smaller

            for (i = 0 ; i < n-1 ; i++) free(m[i]) ;// free the storage allocated to
                                            // to this minor's set of pointers
            free(m) ;                       // free the storage for the original
                                            // pointer to pointer
        }
    }
    return(det) ;
}

