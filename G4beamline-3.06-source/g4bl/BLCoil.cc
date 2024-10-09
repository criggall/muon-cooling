//	BLCoil.cc
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

See also the DISCLAIMER for the routine getSheetField() near the end of
this file.
*/

#ifdef G4BL_GSL

#include "Randomize.hh"

#include "BLAssert.hh"
#include "BLCoil.hh"
#include "BLParam.hh"
#include "BLTime.hh"

#include "gsl/gsl_sf_ellint.h"

std::map<G4String,BLCoil*> BLCoil::mapCoil;

BLCoil::BLCoil()
{
	// provide initial values for fields
	name = "";
	innerRadius = 0.0;
	outerRadius = 0.0;
	length = 0.0;
	nSheets = 0;
	material = "Cu";
	tolerance = 0.002;
	maxR = 0.0;
	minZ = 0.0;
	maxZ = 0.0;
	dR = 0.0;
	dZ = 0.0;
	nR = 0;
	nZ = 0;
	filename = "";
	exactComputation = 0;
	mapFile = "";
	mapBr = 0;
	mapBz = 0;
	sheetR0 = 0.0;
	sheetDR = 0.0;
	norm = 0.0;
	goodMap = false;
	reflectZ = true;
	currentOfMap = 1.0;
}

BLCoil::~BLCoil() 
{
	if(mapBr) delete[] mapBr;
	mapBr = 0;
	if(mapBz) delete[] mapBz;
	mapBz = 0;
}

BLCoil::BLCoil(const BLCoil& r)
{
	name = "";
	innerRadius = r.innerRadius;
	outerRadius = r.outerRadius;
	length = r.length;
	nSheets = r.nSheets;
	material = r.material;
	tolerance = r.tolerance;
	maxR = r.maxR;
	minZ = r.minZ;
	maxZ = r.maxZ;
	dR = r.dR;
	dZ = r.dZ;
	nR = r.nR;
	nZ = r.nZ;
	filename = r.filename;
	exactComputation = r.exactComputation;
	mapFile = r.mapFile;
	mapBr = 0;
	mapBz = 0;
	sheetR0 = r.sheetR0;
	sheetDR = r.sheetDR;
	norm = 0.0;
	goodMap = false;
	reflectZ = r.reflectZ;
	currentOfMap = r.currentOfMap;
}

void BLCoil::printCoil()
{
	printf("coil    %-7s ",name.c_str());
	printf("innerRadius=%.1f ",innerRadius);
	printf("outerRadius=%.1f ",outerRadius);
	printf("length=%.1f ",length);
	printf("material=%s ",material.c_str());
	printf("\n\t\t");
	if(mapFile == "") {
		printf("tolerance=%.3f ",tolerance);
		printf("nSheets=%d ",nSheets);
		printf("\n\t\t");
		printf("maxR=%.1f ",maxR);
		printf("maxZ=%.1f ",maxZ);
		printf("dR=%.1f ",dR);
		printf("dZ=%.1f ",dZ);
		printf("filename=%s ",filename.c_str());
		printf("\n");
	} else {
		printf("mapFile=%s\n",mapFile.c_str());
	}
}

BLCoil *BLCoil::find(G4String name)
{
	if(mapCoil.count(name) > 0)
		return mapCoil[name];
	return 0;
}

void BLCoil::addField(const G4double point[4], G4double field[6],
					G4double currentAmpPerMm2) const
{
	G4double z = point[2];
	G4double r = sqrt(point[0]*point[0] + point[1]*point[1]);

	if(goodMap && exactComputation == 0) {
		BLAssert(mapBr != 0 && mapBz != 0 && dR > 0.0 && dZ > 0.0);
		BLAssert(nR > 0 && nZ > 0);
		if(z < minZ+0.002*mm || z > maxZ-0.002*mm || r > maxR-0.002*mm)
			return;
		G4double sign = 1.0;
		if(reflectZ) {
			// map covers only z >= 0
			if(z < 0.0) {
				sign = -sign;
				z = -z;
			}
		} else {
			z -= minZ;
		}
		// 2-d linear average of the 4 surrounding points
		int i = (int)(r/dR);
		int j = (int)(z/dZ);
if(i >= 0 && i < nR-1 && j >= 0 && j < nZ-1) {
  ;
} else {
  printf("\ni=%d j=%d nR=%d nZ=%d r=%.3f z=%.3f dR=%.3f dZ=%.3f\n",i,j,nR,nZ,r,z,dR,dZ);
  printf("z=%.6f maxZ=%.6f r=%.6f maxR=%.5f\n",z,maxZ,r,maxR);
  fflush(stdout);
}
		BLAssert(i >= 0 && i < nR-1 && j >= 0 && j < nZ-1);
		float fr = 1.0 - (r - i*dR) / dR;
		BLAssert(fr >= 0.0 && fr <= 1.0);
		float fz = 1.0 - (z - j*dZ) / dZ;
		BLAssert(fz >= 0.0 && fz <= 1.0);
	
		G4double Bz = mapBz[j*nR+i]*fr*fz +
		   	mapBz[j*nR+i+1]*(1.0-fr)*fz +
		   	mapBz[(j+1)*nR+i]*fr*(1.0-fz) +
		   	mapBz[(j+1)*nR+i+1]*(1.0-fr)*(1.0-fz);
		G4double Br = mapBr[j*nR+i]*fr*fz +
		   	mapBr[j*nR+i+1]*(1.0-fr)*fz +
		   	mapBr[(j+1)*nR+i]*fr*(1.0-fz) +
		   	mapBr[(j+1)*nR+i+1]*(1.0-fr)*(1.0-fz);
		G4double phi = atan2(point[1], point[0]);
		currentAmpPerMm2 /= currentOfMap;
		field[0] += sign * currentAmpPerMm2 * Br * cos(phi);
		field[1] += sign * currentAmpPerMm2 * Br * sin(phi);
		field[2] += currentAmpPerMm2 * Bz;
		return;
	}

	// Computation for no map
if(nSheets > 0 && sheetDR > 0.0 && sheetR0 > 0.0); else printf("%s: nSheets=%d sheetR0=%.3f sheetDR=%.3f innerRadius=%.3f outerRadius=%.3f length=%.3f\n",name.c_str(),nSheets,sheetR0,sheetDR,innerRadius,outerRadius,length);
	BLAssert(nSheets > 0 && sheetDR > 0.0 && sheetR0 > 0.0);
	G4double Br=0.0, Bz=0.0;
	G4double rSheet = sheetR0;
	G4double currentAmpPerMm = currentAmpPerMm2 * sheetDR;
	for(int i=0; i<nSheets; ++i) {
		G4double br, bz;
		// avoid singularity when r==rSheet
		if(fabs(r-rSheet) < sheetDR/8.0) {
		    getSheetField(r+sheetDR/4.0,z,br,bz,rSheet,length,
		    					currentAmpPerMm);
		    G4double br1, bz1;
		    getSheetField(r-sheetDR/4.0,z,br1,bz1,rSheet,length,
		    					currentAmpPerMm);
		    br = (br+br1)/2.0;
		    bz = (bz+bz1)/2.0;
		} else {
		    getSheetField(r,z,br,bz,rSheet,length,currentAmpPerMm);
		}
		Br += br;
		Bz += bz;
		rSheet += sheetDR;
	}
	G4double phi = atan2(point[1], point[0]);
	field[0] += Br * cos(phi);
	field[1] += Br * sin(phi);
	field[2] += Bz;
}

void BLCoil::generateFieldMap()
{
	if(exactComputation != 0) {
		minZ = -DBL_MAX;
		maxZ = DBL_MAX;
		maxR = DBL_MAX;
		if(nSheets <= 0) nSheets = 1;
		sheetDR = (outerRadius-innerRadius) / nSheets;
		sheetR0 = innerRadius + sheetDR/2.0;
		return;
	}

	if(mapFile != "") {
		readMapFile();
		return;
	}

	struct FileHeader {
		long magic, nR, nZ, nSheets;
		double innerRadius, outerRadius, length, dR, dZ, error;
		// followed by nR*nZ  float-s of mapBr
		// followed by nR*nZ  float-s of mapBz
	};
	// Had serious problems with files written on a different OS that
	// passed the check but generated gross tracking errors and/or
	// exceptions in malloc/free. Probably different compilers padding
	// FileHeader differently.
#ifdef __linux
	const int MAGIC = 0x544A5232;	// 'TJR2'
#endif
#ifdef WIN32
	const int MAGIC = 0x544A5235;	// 'TJR5'
#endif
#ifdef __APPLE__
	const int MAGIC = 0x544A5234;	// 'TJR4'
#endif
#ifdef __CYGWIN__
	const int MAGIC = 0x544A5235;	// 'TJR5'
#endif

	printf("coilmap %-7s tolerance=%.5f [fraction of "
				"Bz(r=0,z=0)]\n",name.c_str(),tolerance);
	fflush(stdout);
	goodMap = false;

	// read FileHeader from filename (if it exists), and check for 
	// consistency with this coil
#ifdef WIN32
	FILE *in = fopen(filename.c_str(),"rb");
#else
	FILE *in = fopen(filename.c_str(),"r");
#endif
	FileHeader fh;
	if(in) {
		if(fread(&fh,sizeof(fh),1,in) != 1) {
			printf("BLCoil::generateFieldMap Cannot read header\n");
			if(in) fclose(in);
			in = 0;
		}
		if(fh.magic != MAGIC) {
			printf("BLCoil::generateFieldMap wrong magic\n");
			if(in) fclose(in);
			in = 0;
		}
		if(fh.error > tolerance) {
			printf("BLCoil::generateFieldMap error > tolerance\n");
			if(in) fclose(in);
			in = 0;
		}
		if(fabs(fh.innerRadius-innerRadius) > 0.001*mm) {
			printf("BLCoil::generateFieldMap wrong innerRadius\n");
			if(in) fclose(in);
			in = 0;
		}
		if(fabs(fh.outerRadius-outerRadius) > 0.001*mm) {
			printf("BLCoil::generateFieldMap wrong outerRadius\n");
			if(in) fclose(in);
			in = 0;
		}
		if(fabs(fh.length-length) > 0.001*mm) {
			printf("BLCoil::generateFieldMap wrong length\n");
			if(in) fclose(in);
			in = 0;
		}
		if((dR > 0.001*mm && fabs(fh.dR-dR) > 0.001*mm)) {
			printf("BLCoil::generateFieldMap wrong dR\n");
			if(in) fclose(in);
			in = 0;
		}
		if((dZ > 0.001*mm && fabs(fh.dZ-dZ) > 0.001*mm)) {
			printf("BLCoil::generateFieldMap wrong dZ\n");
			if(in) fclose(in);
			in = 0;
		}
		if((fh.nR-1)*fh.dR < maxR || (fh.nZ-1)*fh.dZ < maxZ) {
			printf("BLCoil::generateFieldMap wrong maxR or maxZ\n");
			if(in) fclose(in);
			in = 0;
		}
	}
	// setup data from file and read the maps from file
	if(in) {
		G4double orgMaxZ = maxZ;
		G4double orgMaxR = maxR;
		G4int orgNsheets = nSheets;
		nR = fh.nR;
		nZ = fh.nZ;
		dR = fh.dR;
		dZ = fh.dZ;
		nSheets = fh.nSheets;
		sheetDR=(outerRadius-innerRadius)/nSheets;
		sheetR0=innerRadius + sheetDR/2.0;
		norm = 0.0;
		maxR = (nR-1)*dR;
		maxZ = (nZ-1)*dZ;
		minZ = -maxZ;
		if(mapBr) delete[] mapBr;
		if(mapBz) delete[] mapBz;
		mapBr = new float[nZ*nR];
		mapBz = new float[nZ*nR];
		if(!mapBr || !mapBz) {
			G4Exception("coil","Out of memory",FatalException,"");
		}
		if(fread(mapBr,sizeof(float),nR*nZ,in) != (size_t)nR*nZ ||
		   fread(mapBz,sizeof(float),nR*nZ,in) != (size_t)nR*nZ) {
			fclose(in);
			in = 0;
			printf("BLCoil::generateFieldMap cannot read data\n");
			// restore invalid data that came from file
			maxZ = orgMaxZ;
			minZ = -maxZ;
			reflectZ = true;
			maxR = orgMaxR;
			nSheets = orgNsheets;
		}
	}
	// if file is correct, use it.
	if(in) {
		fclose(in);
		in = 0;
		goodMap = true;
		printf("coilmap %-7s read file '%s'  "
		    "dR=%.1f dZ=%.1f\n",name.c_str(),filename.c_str(),dR,dZ);
		minZ = -maxZ;
		reflectZ = true;
		currentOfMap = 1.0;
		return;
	}

	// file is not correct - generate map

	if(maxR <= 0.0)
		determineMaxR();
	if(maxZ <= 0.0)
		determineMaxZ();
	if(nSheets <= 0) {
		determineNsheets();
	} else {
		sheetDR = (outerRadius-innerRadius) / nSheets;
		sheetR0 = innerRadius + sheetDR/2.0;
	}

	// initial estimate of dR,dZ,nR,nZ
	if(dR <= 0.0)
		dR = maxR/10;
	nR = (int)floor(maxR/dR) + 2;
	//@@ for now, assume dZ = dR....
	if(dZ <= 0.0)
		dZ = dR;
	nZ = (int)floor(maxZ/dZ) + 2;

	// loop until we achieve the required tolerance. Limit to 8 doublings.
	const int MAX_ITER=8;
	int iter;
	G4double err = 1.0;
	for(iter=0; iter<MAX_ITER; ++iter) {
		// for simplicity, discard all previous computations
		if(mapBr) delete[] mapBr;
		if(mapBz) delete[] mapBz;
		mapBr = new float[nZ*nR];
		mapBz = new float[nZ*nR];
		if(!mapBr || !mapBz) {
			G4Exception("coil","Out of memory",FatalException,"");
		}
		err = estimateMapError();
		if(err < tolerance)
			break;
		dR /= 2.0;
		dZ /= 2.0;
		nR = (int)floor(maxR/dR) + 2;
		nZ = (int)floor(maxZ/dZ) + 2;
	}
	if(iter >= MAX_ITER) {
		G4Exception("coilmap","Failed to achieve required accuracy",
							FatalException, "");
	}
	printf("coilmap %-7s nR=%d nZ=%d dR=%.1f dZ=%.1f "
			"error=%.6f\n",name.c_str(),nR,nZ,dR,dZ,err);

	// fill the maps
	//@@ for simplicity, ignore all previous computations
	goodMap = false;
	time_t start = BLTime::time();
	for(int i=0; i<nR; ++i) {
		G4double r = i * dR;
		for(int j=0; j<nZ; ++j) {
			G4double z = j * dZ;
			int k = j * nR + i;
			BLAssert( k >= 0 && k < nZ*nR);
			G4double Br,Bz;
			getUnitField(r,z,Br,Bz);
			mapBr[k] = Br;
			mapBz[k] = Bz;
		}
//		printf("coilmap %-7s  %d%% %ld seconds\r",
//			name.c_str(),(100*(i+1))/nR,(long)(BLTime::time()-start));
//		fflush(stdout);
	}
	printf("coilmap %-7s  DONE %ld seconds\n",
		name.c_str(),(long)(BLTime::time()-start)); 
	fflush(stdout);

	// write the map
	if(filename != "") {
#ifdef WIN32
		FILE *out = fopen(filename.c_str(),"wb");
#else
		FILE *out = fopen(filename.c_str(),"w");
#endif
		if(out) {
			FileHeader fh;
			fh.magic = MAGIC;
			fh.nR = nR;
			fh.nZ = nZ;
			fh.nSheets = nSheets;
			fh.innerRadius = innerRadius;
			fh.outerRadius = outerRadius;
			fh.length = length;
			fh.dR = dR;
			fh.dZ = dZ;
			fh.error = err;
			if(fwrite(&fh,sizeof(fh),1,out) != 1 ||
			   fwrite(mapBr,sizeof(float),nR*nZ,out) != (size_t)nR*nZ ||
			   fwrite(mapBz,sizeof(float),nR*nZ,out) != (size_t)nR*nZ) {
				printf("coilmap %-7s ERROR "
					"writing file '%s'\n",
					name.c_str(),filename.c_str());
			}
			fclose(out);
			printf("coilmap %-7s wrote file '%s'   %lu bytes\n",
				name.c_str(),filename.c_str(),
				(long)(sizeof(fh)+sizeof(float)*2*nR*nZ));
		}
	}

	goodMap = true;
	minZ = -maxZ;
	reflectZ = true;
	currentOfMap = 1.0;
}

void BLCoil::readMapFile()
{
	FILE *in = fopen(mapFile.c_str(),"r");
	if(!in) {
		BLCommand::printError("Cannot open mapFile '%s'\n",mapFile.c_str());
		return;
error:
		BLCommand::printError("Syntax error in mapFile '%s'\n",mapFile.c_str());
		return;
	}

	char line[512];

	// read first line (parameters of the map)
	// units are millimeters, Tesla, and Amperes/mm^2.
	do {
		if(!fgets(line,sizeof(line),in)) goto error;
	} while(line[0] == '#');

	BLArgumentVector argv;
	BLArgumentMap argm;
	if(BLCommand::parseArgs(line,argv,argm)) goto error;

	// sample parameter line:
	// z0=0.0 dz=2.0 nz=1501 dr=5.0 nr=13 reflectZ=1 current=-101.1
	G4double z0 = atof(argm["z0"].c_str())*mm;
	dZ = atof(argm["dz"].c_str())*mm;
	nZ = atoi(argm["nz"].c_str());
	dR = atof(argm["dr"].c_str())*mm;
	nR = atoi(argm["nr"].c_str());
	reflectZ = (atoi(argm["reflectZ"].c_str()) != 0);
	currentOfMap = atof(argm["current"].c_str());
	maxR = (nR-1)*dR;
	maxZ = z0 + (nZ-1)*dZ;
	if(reflectZ) {
		minZ = -maxZ;
		if(z0 != 0.0) goto error;
	} else {
		minZ = z0;
	}
	if(argv.size() > 0 || nZ <= 0 || nR <= 0 || maxR <= 0.0 ||
					maxZ <= 0.0 || currentOfMap == 0.0)
		goto error;
	if(mapBr) delete[] mapBr;
	if(mapBz) delete[] mapBz;
	mapBr = new float[nZ*nR];
	mapBz = new float[nZ*nR];
	if(!mapBr || !mapBz) {
		G4Exception("coil","Out of memory",FatalException,"");
	}

	// read the map -- there are two maps, Bz and Br, each introduced
	// by a line with its name in column 1. Each map consists of nZ
	// lines of nR floats.
	int i=-1;
	float *map=0;
	while(fgets(line,sizeof(line),in)) {
		if(line[0] == '#') continue;
		if(strncmp(line,"Bz",2) == 0) {
			if(i != -1 && i != nZ) goto error;
			i = 0;
			map = mapBz;
			continue;
		} else if(strncmp(line,"Br",2) == 0) {
			if(i != -1 && i != nZ) goto error;
			i = 0;
			map = mapBr;
			continue;
		}
		if(i < 0 || i >= nZ || !map) goto error;
		char *p = line;
		for(int j=0; j<nR; ++j) {
			char *q=0;
			if(*p == '\0') goto error;
			*map++ = strtod(p,&q)*tesla;
			if(q == 0 || p == q) goto error;
			p = q;
		}
		++i;
	}
	if(i != nZ) goto error;
	fclose(in);
	goodMap = true; 
}

#ifdef STUB
void BLCoil::mapErrorsNTuple(G4String tuplename, G4String  category,
						G4int nevents)
{
	static const char *var[] = {"R","Z","Br","Bz","Er","Ez"};
	static int nvar = 6;
	float data[6];

	NTuple nt(tuplename,category,nvar,var);

	BLCoil best("temp",innerRadius,outerRadius,length,99,0.0,maxR,maxZ);
	if(norm == 0.0) {
		G4double Br,Bz;
		best.getUnitField(0.0,0.0,Br,Bz);
		norm = 1.0/Bz;
	}

	for(int i=0; i<nevents; ++i) {
		G4double r,z,Br0,Bz0,Br1,Bz1;
		r = G4UniformRand() * maxR;
		z = G4UniformRand() * 2.0 * maxZ - maxZ;
		best.getUnitField(r,z,Br0,Bz0);
		getUnitField(r,z,Br1,Bz1);
		data[0] = r;
		data[1] = z;
		data[2] = Br0;
		data[3] = Bz0;
		data[4] = fabs((Br0-Br1)*norm);
		data[5] = fabs((Bz0-Bz1)*norm);
		nt.append(data);
//		if(i%100 == 0) {
//			printf("BLCoil::mapErrorsNTuple '%s' %.0f%%\r",
//				tuplename.c_str(),(100.0*i)/nevents);
//			fflush(stdout);
//		}
		NTuple::update();
	}
	printf("BLCoil::mapErrorsNTuple '%s' DONE\n",tuplename.c_str());
	fflush(stdout);
}
#endif // STUB

void BLCoil::displayMapValues(G4double r, G4double z)
{
	if(!goodMap) {
		printf("BLCoil '%s' has no good map\n",name.c_str());
		return;
	}
	if(z < 0.0) {
		printf("z must be >= 0\n");
		return;
	}

	BLCoil best("temp",innerRadius,outerRadius,length,99,0.0,maxR,maxZ);
	G4double Br,Bz;
	getUnitField(r,z,Br,Bz);

	int ii = (int)(r/dR);
	int jj = (int)(z/dZ);
	if(ii+1 >= nR || jj+1 >= nZ) {
		printf("out of range\n");
		return;
	}
	printf("r=%.3f z=%.3f Bz=%.6f T  Br=%.6f T  Berr=%.6f T\n",
			r,z,Bz/tesla,Bz/tesla,
			computeError(r,z,best,false)/norm/tesla);
	G4double bzmin=Bz, bzmax=Bz, bzavg=0.0;
	G4double brmin=Br, brmax=Br, bravg=0.0;
	for(int i=ii; i<=ii+1; ++i) {
	    for(int j=jj; j<=jj+1; ++j) {
		if(bzmin > mapBz[j*nR+i]) bzmin = mapBz[j*nR+i];
		if(bzmax < mapBz[j*nR+i]) bzmax = mapBz[j*nR+i];
		if(brmin > mapBr[j*nR+i]) brmin = mapBr[j*nR+i];
		if(brmax < mapBr[j*nR+i]) brmax = mapBr[j*nR+i];
		bzavg += mapBz[j*nR+i]/4.0;
		bravg += mapBr[j*nR+i]/4.0;
		printf("Grid: r=%.3f z=%.3f Bz=%.6f T  Br=%.6f T\n",
			dR*i,dZ*j,mapBz[j*nR+i]/tesla,mapBr[j*nR+i]/tesla);
	    }
	}
	printf("BzVariation=%.4f BrVariation=%.4f (fractions of avg)\n",
		(bzmax-bzmin)/bzavg,(brmax-brmin)/bravg);
}

void BLCoil::getUnitField(G4double r, G4double z, G4double& Br, G4double& Bz)
{
	G4double point[4];
	point[0] = point[1] = point[2] = point[3] = 0.0;
	G4double f[6];
	f[0]=f[1]=f[2]=f[3]=f[4]=f[5]=0.0;
	point[0] = r;
	point[2] = z;
	addField(point,f,1.0);
	Br = f[0];
	Bz = f[2];
}

void BLCoil::determineMaxR()
{
	G4double Br,Bz;
	BLCoil best("temp",innerRadius,outerRadius,length,99,0.0,
					10.0*outerRadius,100.0*outerRadius);
	if(norm == 0.0) {
		best.getUnitField(0.0,0.0,Br,Bz);
		norm = 1.0/Bz;
	}

	// step out in maxR until 0.5*tolerance is achieved
	for(maxR=outerRadius+2.0*innerRadius;
	    maxR<=outerRadius+32.0*innerRadius;
	    maxR+=innerRadius) {
		best.getUnitField(maxR,0.0,Br,Bz);
		Bz = fabs(Bz);
		if(Bz*norm < 0.5*tolerance) break;
	}
	if(maxR <= outerRadius+32.0*innerRadius) {
		printf("BLCoil::determineMaxR '%s'  maxR=%.1f err=%.6f\n",
					name.c_str(),maxR,Bz*norm);
	} else {
		G4Exception("coilmap determineMaxR",
					"Failed to achieve required accuracy",
					FatalException, "");
	}
}

void BLCoil::determineMaxZ()
{
	G4double Br,Bz;
	BLCoil best("temp",innerRadius,outerRadius,length,99,0.0,maxR,
							100.0*outerRadius);
	if(norm == 0.0) {
		best.getUnitField(0.0,0.0,Br,Bz);
		norm = 1.0/Bz;
	}

	// step out in maxZ until 0.5*tolerance is achieved
	for(maxZ=length/2.0+4.0*innerRadius; maxZ<=length/2.0+32.0*innerRadius;
							maxZ+=innerRadius) {
		best.getUnitField(0.0,maxZ,Br,Bz);
		if(Bz*norm < 0.5*tolerance) break;
	}
	if(maxZ <= length/2.0+32.0*innerRadius) {
		printf("BLCoil::determineMaxZ '%s'  maxZ=%.1f err=%.6f\n",
					name.c_str(),maxZ,Bz*norm);
	} else {
		G4Exception("coilmap determineMaxZ",
					"Failed to achieve required accuracy",
					FatalException, "");
	}
}

void BLCoil::determineNsheets()
{
	BLCoil best("temp",innerRadius,outerRadius,length,99,0.0,maxR,maxZ);

	G4double exclude = 0.1 * innerRadius;
	G4double err = 0.0;
	const int MAX_SHEETS=40;
	for(nSheets=2; nSheets<=MAX_SHEETS; nSheets+=2) {
		sheetDR = (outerRadius-innerRadius) / nSheets;
		sheetR0 = innerRadius + sheetDR/2.0;
		err = 0.0;
		for(int i=0; i<50; ++i) {
			G4double r = innerRadius - exclude;
			G4double z = G4UniformRand() * length/2.0;
			if(r >= maxR) r = maxR;
			G4double e = computeError(r,z,best,false);
			if(e > err) err = e;
			r = innerRadius + (outerRadius-innerRadius)*G4UniformRand();
			if(r >= maxR) continue;
			z = length/2.0 + exclude;
			e = computeError(r,z,best,false);
			if(e > err) err = e;
			if(err > tolerance/4.0) break;
		}
		if(err <= tolerance/4.0) break;
	}

	if(nSheets <= MAX_SHEETS) {
		printf("BLCoil::determineNsheets '%s'  nSheets=%d err=%.6f\n",
					name.c_str(),nSheets,err);
	} else {
		G4Exception("coilmap determineNsheets",
					"Failed to achieve required accuracy",
					FatalException, "");
	}
}

G4double BLCoil::estimateMapError()
{
	BLCoil best("temp",innerRadius,outerRadius,length,99,0.0,maxR,maxZ);

	// generate 1000 random positions, "close"
	// to the coil in z, but exclude the coil (and nearby).
	G4double exclude = 0.1 * innerRadius;
	G4double err = 0.0;
	G4double genZ = length;
	if(length/2.0+exclude+10.0*mm > genZ) genZ = length/2.0+exclude+10.0*mm;
	for(int i=0; i<1000; ++i) {
		G4double r, z;
		do {
			r = G4UniformRand() * maxR;
			z = G4UniformRand() * genZ;
		} while(z < length/2+exclude &&
		        r > innerRadius-exclude &&
			r < outerRadius+exclude);
		G4double e = computeError(r,z,best,true);
		if(e > err) err = e;
		if(err > tolerance) {
			break;
		}
//		if(i%10 == 0)
//			printf("BLCoil::estimateMapError nR=%d nZ=%d %d%%\r",
//				nR,nZ,(100*i)/1000), fflush(stdout);
	}

	return err;
}

G4double BLCoil::computeError(G4double r, G4double z, BLCoil &best,
							G4bool makeMap)
{
	G4double Br,Bz;
	if(z < -maxZ || z > maxZ || r > maxR) return 0.0;
	G4bool originalGoodMap = goodMap;

	if(norm == 0.0) {
		best.getUnitField(0.0,0.0,Br,Bz);
		norm = 1.0/Bz;
	}

	if(makeMap) {
		// generate fieldMap for the 4 surrounding grid points
		BLAssert(mapBz != 0 && mapBr != 0 && nR > 0 && nZ > 0);
		BLAssert(z >= 0.0);
		goodMap = false;
		int i = (int)(r/dR);
		int j = (int)(z/dZ);
		BLAssert(i >= 0 && i < nR-1 && j >= 0 && j < nZ-1);
	
		getUnitField(i*dR,j*dZ,Br,Bz);
		mapBz[j*nR+i] = Bz;
		mapBr[j*nR+i] = Br;

		getUnitField((i+1)*dR,j*dZ,Br,Bz);
		mapBz[j*nR+(i+1)] = Bz;
		mapBr[j*nR+(i+1)] = Br;
	
		getUnitField(i*dR,(j+1)*dZ,Br,Bz);
		mapBz[(j+1)*nR+i] = Bz;
		mapBr[(j+1)*nR+i] = Br;

		getUnitField((i+1)*dR,(j+1)*dZ,Br,Bz);
		mapBz[(j+1)*nR+(i+1)] = Bz;
		mapBr[(j+1)*nR+(i+1)] = Br;

		goodMap = true;
	}

	G4double bestBr,bestBz;
	best.getUnitField(r,z,bestBr,bestBz);
	getUnitField(r,z,Br,Bz);
	G4double eBr = fabs(Br-bestBr)*norm;
	G4double eBz = fabs(Bz-bestBz)*norm;

	goodMap = originalGoodMap;
	return (eBr>eBz ? eBr : eBz);
}

// Originally this came from BTSheet.cc:
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of     *
// * FERMILAB.                                                        *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// BTSheet.cc
//
// Created: V.Daniel Elvira (4/19/00)
//                                      
//           
//   The BTSheet Class inherits from G4MagneticField. The class objects are
//   field maps generated by an infinitelly thin solenoidal current sheets.
//   The class data members are all the parameters necessary to generate
//   analytically a magnetic field in r,z space (there is phi symmetry).
//   No geometric volumes or materials are associated with BTSheets. 
//
// Modified: Tom Roberts (2003/02/28)
//
//   Changed from an object with obfuscated interface to a simple function.
//   Field value verified to be correct within 0.2% for a long solenoid and
//   for a pair of Helmholtz coils (at the geometric center of each).
//
// Modified: Tom Roberts 3/28/2003
//   SpecialFunctions/EllipticIntegrals.h => gsl/gsl_sf_ellint.h
//   This should be transparent, as SpecialFunctions/EllipticIntegrals
//   seems to be an interface to the gsl functions. I verified SINGLE accuracy
//   makes less than 0.0000001 difference in the values, except when ellpi 
//   gets large (>1000); in all cases 6 digit accuracy is maintained. But I
//   switched to DOUBLE accuracy because the time increase was only ~30% and
//   this is only a small part of the program total execution time.

void BLCoil::getSheetField(G4double r, G4double z, G4double& Br, G4double& Bz,
	G4double sheetR, G4double sheetLen, G4double currentAmpPerMm) const
{
  // TJR: within this routine the sheet extends from z=0 to z=sheetLen,
  //      externally it is centered at z=0.
  z += sheetLen/2.0;

//printf("BLCoil::getSheetField r=%.6g z=%.6g sheetR=%.6g sheetLen=%.6g\n",r,z,sheetR,sheetLen); fflush(stdout);
  G4double k; 

  G4double rposition = r/meter;
  G4double zposition = z/meter;

  G4double rad = sheetR/meter;
  G4double len = sheetLen/meter;
  G4double cur = currentAmpPerMm;

  //  Left-end of sheet

  G4double rho = rposition;
  z = -zposition;
  
  G4double r1 = sqrt( pow(rad + rho,2) + z*z );
  k = 2.0 / r1 * sqrt( rad * rho );
  G4double c2 = -4.0 * rad*rho / pow(rad+rho,2);

  G4double f1 = z * rad / ((rad + rho) * r1);
  G4double f2 = (rad - rho) / (2.0 * rad);

//printf("gsl_sf_ellint_E k=%.6g\n",k); fflush(stdout);
  G4double elle = gsl_sf_ellint_E(pi/2.0,k,GSL_PREC_DOUBLE);
//printf("gsl_sf_ellint_F k=%.6g\n",k); fflush(stdout);
  G4double ellf = gsl_sf_ellint_F(pi/2.0,k,GSL_PREC_DOUBLE);
//printf("gsl_sf_ellint_P k=%.6g c2=%.6g\n",k,c2); fflush(stdout);
  G4double ellpi = gsl_sf_ellint_P(pi/2.0,k,c2,GSL_PREC_DOUBLE);

  G4double bz1 = f1 * ( ellf + f2 * ( ellpi - ellf ) );

  G4double br1;

  if( rho != 0.0 ) {
    f1 = r1 / (4.*rho);
    br1 = f1 * ( 2. * ( ellf - elle ) - k * k * ellf );
  }
  else {
    br1 = 0.0;
  }

  G4double bz = -bz1;
  G4double br = -br1;      

  //  Right-end of sheet
      
  z = z + len;
  r1 = sqrt( pow(rad + rho,2) + z*z );
  k = 2.0 / r1 * sqrt( rad * rho );
  c2 = -4.0 * rad * rho / pow(rad + rho,2);

  f1 = z * rad / ( (rad + rho) * r1);
  f2 = (rad - rho) / (2.0 * rad);

//printf("gsl_sf_ellint_E k=%.6g\n",k); fflush(stdout);
  elle = gsl_sf_ellint_E(pi/2.0,k,GSL_PREC_DOUBLE);
//printf("gsl_sf_ellint_F k=%.6g\n",k); fflush(stdout);
  ellf = gsl_sf_ellint_F(pi/2.0,k,GSL_PREC_DOUBLE);
//printf("gsl_sf_ellint_P k=%.6g c2=%.6g\n",k,c2); fflush(stdout);
  ellpi = gsl_sf_ellint_P(pi/2.0,k,c2,GSL_PREC_DOUBLE);

  bz1 = f1 * ( ellf + f2 * ( ellpi - ellf ) );

  if ( rho != 0.0 ) {
    f1 = r1 / (4.0 * rho);
    br1 = f1 * ( 2.0 * ( ellf - elle ) - k * k * ellf );
  }      
  else {
    br1 = 0.0;
  }

  bz = bz + bz1;       // note change of sign
  br = br + br1;

  G4double permfactor = 4.0e-7;

  Br = permfactor * cur * br;
  Bz = permfactor * cur * bz;

}

#else // G4BL_GSL
int dummycoil=0;
#endif // G4BL_GSL
