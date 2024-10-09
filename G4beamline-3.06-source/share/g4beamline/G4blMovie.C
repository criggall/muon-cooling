//	G4blMovie.C - Root ACLiC macro to generate a movie
//
//	This macro must be compiled, as CINT cannot handle std::vector and
//	std::map correctly; the related root classes are too restrictive.

#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "TROOT.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TMarker.h"
#include "TColor.h"
#include "TText.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TSystem.h"

TPad *firstpad=0;	// convenience for command-line users

typedef std::map<TString,TString> Args;

class Setup *setup=0;
class Particle *particle=0;
class MoviePad;

/// TrajIdType is 64-bit trajectory identifier
typedef Long64_t TrajIdType;
TrajIdType getTrajId(int eventID, int trackID) {
	return (((TrajIdType)eventID)<<16) | (((TrajIdType)trackID)&0xFFFF);
}

static int nError=0;

/**	class Command is the base class for a command in the input file 
**/
class Command {
public:
	Command() { }
	virtual ~Command() { }
	bool handleArgs(TString cmd, Args args) {
		bool retval=true;
		for(Args::iterator i=args.begin(); i!=args.end(); ++i) {
			if(!defineArgs(i->first,i->second)) {
			    retval = false;
			    printf("Invalid argument to %s: '%s=%s'\n",
				cmd.Data(),i->first.Data(),i->second.Data());
			    ++nError;
			}
		}
		return retval;
	}
	virtual bool defineArgs(TString name, TString value) = 0;
	static void read(TString filename) {
		ifstream in(filename.Data());
		if(!in.good()) {
			printf("Cannot open '%s'\n",filename.Data());
			exit(1);
		}
		TString line;
		while(line.ReadLine(in).good()) {
			while(line.EndsWith("\\")) {
				line.Remove(TString::kTrailing,'\\');
				line.Append(" ");
				TString tmp;
				if(!tmp.ReadLine(in).good()) break;
				line.Append(tmp);
			}
			line.ReplaceAll("\t"," ");
			if(line.BeginsWith("#")) continue;
			printf("\n%s\n",line.Data());
			TObjArray *a=line.Tokenize(" ");
			a->SetOwner(1);
			int n=a->GetEntriesFast();
			Args args;
			TString cmd=((TObjString *)a->At(0))->String();
			for(int i=1; i<n; ++i) {
			    TString word=((TObjString *)a->At(i))->String();
			    int j=word.Index("=");
			    if(j > 0) {
			    	TString name = word(0,j);
				TString value = word(j+1,99999);
				args[name] = value;
			    } else {
			    	printf("Invalid argument '%s' to command %s\n",
					word.Data(),cmd.Data());
			        ++nError;
			    }
/***
			    TObjArray *b=word.Tokenize("=");
			    b->SetOwner(1);
			    if(b->GetEntriesFast() == 2) {
			    	TString name=((TObjString *)b->At(0))->String();
			    	TString val=((TObjString *)b->At(1))->String();
			    	args[name] = val;
			    } else {
			    	printf("Invalid argument '%s' to command %s\n",
					word.Data(),cmd.Data());
			        ++nError;
			    }
			    delete b;
***/
			}
			delete a;
			doCommand(cmd,args);
		}
	}
	static void doCommand(TString cmd, Args args);
};

/**	the particle command assigns colors to PDGid-s
**/
struct Particle : public Command {
	std::map<int,TString> colorMap;
public:
	Particle() : Command(), colorMap() { colorMap[0] = "#000000"; }
	bool defineArgs(TString name, TString value) {
		colorMap[name.Atoi()] = value;
		return true;
	}
	TString color(int pdgId) { 
		TString v = colorMap[pdgId];
		if(v == "") v = colorMap[0];
		return v;
	}
};

/**	the row command defines a new row in the movie.
**/
struct Row : public Command {
	float height;
	int firstPad;
	int lastPad;
	int pixelHeight;
public:
	Row(int n) : Command() {
		height = 1.0;
		firstPad =  n;
		lastPad = -99;
		pixelHeight = 0;
	}
	Row(Args args, int n) : Command() {
		height = 1.0;
		firstPad =  n;
		lastPad = -99;
		pixelHeight = 0;
		handleArgs("row",args);
	}
	bool defineArgs(TString name, TString value) {
		if(name == "height") { height = value.Atof(); return true; }
		return false;
	}
};

/**	the setup command specifies basic parameters for the movei
**/
struct Setup : public Command {
	TString borderColor;
	int borderSize;
	TString outputFile;
	TString pictureType;
	float tMin;
	float tMax;
	float duration;
	float frameRate;
	int marker;
	float markerSize;
	int windowWidth;
	int windowHeight;
	int textHeight;
	TString background;
	std::vector<Row*> rows;
	std::vector<MoviePad*> pads;
	TCanvas *canvas;
public:
	Setup() : Command(), rows(), pads() {
		borderColor="#C0C0C0";
		borderSize=5;
		outputFile = "m.mov";
		pictureType="jpg";
		tMin = 0.001;
		tMax = 20.0;
		duration = 20.0;
		frameRate = 24.0;
		marker = 21;
		markerSize = 1.0;
		windowWidth = 500;
		windowHeight = 500;
		textHeight = 15;
		background = "#F3F3F3";
		canvas = 0;

		rows.push_back(new Row(pads.size()));
	}
	bool defineArgs(TString name, TString value) {
	    if(name=="borderColor") { borderColor = value; return true; }
	    if(name=="borderSize") { borderSize = value.Atoi(); return true; }
	    if(name=="outputFile") { outputFile = value; return true; }
	    if(name=="pictureType") { pictureType = value; return true; }
	    if(name=="tMin") { tMin = value.Atof(); return true; }
	    if(name=="tMax") { tMax = value.Atof(); return true; }
	    if(name=="duration") { duration = value.Atof(); return true; }
	    if(name=="frameRate") { frameRate = value.Atof(); return true; }
	    if(name=="marker") { marker = value.Atoi(); return true; }
	    if(name=="markerSize") { markerSize = value.Atof(); return true; }
	    if(name=="windowWidth") { windowWidth = value.Atoi(); return true; }
	    if(name=="windowHeight") { windowHeight = value.Atoi(); return true; }
	    if(name=="textHeight") { textHeight = value.Atoi(); return true; }
	    if(name=="background") { background = value; return true; }
	    return false;
	}
	void layout();
	void drawOneFrame(double t);
	void drawMovie();
};

/**	class MoviePad is a base class for a pad that is displayed in the movie.
 *
 *	NOTE: the pad title must be set before MoviePad::prepareForDrawing()
 *	is called. This happens in MoviePad::defineArgs() for the title
 *	argument, so it only needs to be done if it is to be changed.
**/
class MoviePad : public Command, public TPad {
protected:
	TString cmdName;
	TString title;
	float xMin, xMax, yMin, yMax;	// coordinate limits
	TString background;
	int textHeight;
public:
	MoviePad(TString name) : Command(), TPad(), title("_") {
		cmdName = name;
		xMin = -100.0;
		xMax = 100.0;
		yMin = -100.0;
		yMax = 100.0;
		background = setup->background;
		textHeight = setup->textHeight;
		setup->pads.push_back(this);
	}
	virtual ~MoviePad() { }
	virtual bool defineArgs(TString name, TString value) {
if(name=="background") printf("MoviePad::defineArgs background=%s\n",value.Data());
	    if(name=="title") { title=value; SetTitle(title);  return true; }
	    if(name=="xMin") { xMin=value.Atof(); return true; }
	    if(name=="xMax") { xMax=value.Atof(); return true; }
	    if(name=="yMin") { yMin=value.Atof(); return true; }
	    if(name=="yMax") { yMax=value.Atof(); return true; }
	    if(name=="background") { background=value; return true; }
	    if(name=="textHeight") { textHeight = value.Atoi(); return true; }
	    return false;
	}
	virtual void init(Args args) {
		handleArgs(cmdName,args);
	}
	virtual void prepareForDrawing() {
		cd();
		Range(xMin,yMin,xMax,yMax);
		SetFillColor(TColor::GetColor(background));
		TString tmp=GetTitle();
		if(tmp != "") {
			TText *txt=new TText(0.5,0.97,tmp);
			txt->SetTextAlign(23); // h=centered, v=top
			txt->SetNDC();
			txt->SetTextSize(getTextSize());
			txt->Draw();
		}
		// draw border as a TPolyLine
		SetBorderMode(0);
		float cornerX[5],cornerY[5];
		cornerX[0] = cornerX[1] = cornerX[4] = 0.0001;
		cornerX[2] = cornerX[3] = 0.9999;
		cornerY[0] = cornerY[3] = cornerY[4] = 0.0001;
		cornerY[1] = cornerY[2] = 0.9999;
		TPolyLine *p = new TPolyLine(5,cornerX,cornerY);
		p->SetNDC();
		p->SetLineWidth(setup->borderSize);
		p->SetLineColor(TColor::GetColor(setup->borderColor));
		p->Draw();
	}
	virtual void drawPadContents(double t) { 
		(void)t;
	}
	double getTextSize() { return (double)textHeight /
					(double)(GetAbsHNDC()*GetWh()); }
	ClassDef(MoviePad,1);
};

/// Trajectory contains 1 trajectory, can interpolate between points
struct Trajectory {
	float tMin, tMax;	// values for existence of this trajectory
	float t, v1, v2;	// values of trajectory at most recent setT()
	struct TrajPoint { 
		float t;
		float v1, v2; 
		TrajPoint() { t=v1=v2=-7.0e30; }
		TrajPoint(float _t, float _v1, float _v2)
			{ t=_t; v1=_v1, v2=_v2; }
	};
	std::vector<TrajPoint> points;
	int index;
	bool bookends;	// bookends simplify the logic in setT (high usage)
public:
	Trajectory() : points() 
		{ v1=0; v2=0; tMin=1.0e30, tMax=-1.0e30; index=0; 
		  bookends=false;}
	void appendPoint(float _t, float _v1, float _v2)
		{ if(_t < tMin) tMin = _t;
		  if(_t > tMax) tMax = _t;
		  points.push_back(TrajPoint(_t,_v1,_v2)); }
	void setBookends() 
		{ if(bookends || points.size() < 2) return;
		  bookends = true;
		  tMin -= 1.0e-6;
		  points[0].t -= 2.0e-6;
		  points[points.size()-1].t += 1.0e-6; }
	// interpolate trajectory values; returns false if should not be drawn
	bool setT(float _t) { // interpolate trajecory to time _t
		setBookends();
		if(_t < tMin || _t > tMax || points.size() < 2) return false;
		if(_t < points[index].t) index = 0;
		while(index < (int)points.size()-2 && _t > points[index+1].t) 
			++index;
		assert(index>=0 && index+1<(int)points.size());
		TrajPoint *a=&points[index], *b=&points[index+1];
		assert(_t>=a->t && _t <= b->t);
		float f=(_t-a->t) / (b->t-a->t);
		t = _t;
		v1 = (1.0-f)*a->v1 + f*b->v1;
		v2 = (1.0-f)*a->v2 + f*b->v2;
		return true;
	}
};

///	ReferenceTrack contains the reference track; can interpolate 
///	between points.
struct ReferenceTrack {
	float t, Zref, Ptot;	// trajectory values at most recent setT()
	Trajectory traj;
public:
	ReferenceTrack(TNtuple *ntuple) : traj() {
		ntuple->SetBranchAddress("T",&t);
		ntuple->SetBranchAddress("Zref",&Zref);
		ntuple->SetBranchAddress("Ptot",&Ptot);
		long nEntries = ntuple->GetEntries();
		for(long i=0; i<nEntries; ++i) {
			ntuple->GetEntry(i);
			traj.appendPoint(t,Zref,Ptot);
		}
		printf("The Reference Track has %ld points, %.3f to %.3f ns\n",
			traj.points.size(),traj.tMin, traj.tMax);
	}
	bool setT(float _t) { // interpolate trajectory to time _t
		if(!traj.setT(_t)) return false;
		t = traj.t;
		Zref = traj.v1;
		Ptot = traj.v2;
		return true;
	}
};

///	TrackSet contains a set of Tracks; can interpolate between the
///	points of each track.
class TrackSet {
	std::map<TrajIdType,Trajectory> trajMap;
	std::map<TrajIdType,TMarker*> markerMap;
	std::map<TrajIdType,int> pdgIdMap;
public:
	TrackSet(TNtuple *nt, TString exprX, TString exprY)
				: trajMap(), markerMap(), pdgIdMap() {
		TTreeFormula *fx = new TTreeFormula("x",exprX,nt);
		if(fx->GetTree() == 0) {
			printf("Invalid formula for X '%s'\n",exprX.Data());
			exit(2);
		}
		TTreeFormula *fy = new TTreeFormula("y",exprY,nt);
		if(fy->GetTree() == 0) {
			printf("Invalid formula for Y '%s'\n",exprY.Data());
			exit(2);
		}
		Long64_t nEntries = nt->GetEntries();
		float t, EventID, TrackID, PDGid;
		nt->SetBranchAddress("t",&t);
		nt->SetBranchAddress("EventID",&EventID);
		nt->SetBranchAddress("TrackID",&TrackID);
		nt->SetBranchAddress("PDGid",&PDGid);
		float earliest=1.0e30, latest=-1.0e30;
		for(Long64_t i=0; i<nEntries; ++i) {
			if(i%1000 == 0)
				printf("Reading Movie/Tracks %lld/%lld\r",
					i,nEntries); fflush(stdout);
			if(nt->GetEntry(i) < 0) continue;
			float v1 = fx->EvalInstance();
			float v2 = fy->EvalInstance();
			TrajIdType id = getTrajId(EventID,TrackID);
			trajMap[id].appendPoint(t,v1,v2);
			pdgIdMap[id] = (int)PDGid;
			if(t < earliest) earliest = t;
			if(t > latest) latest = t;
		}
		printf("Found %ld tracks, from %.3f to %.3f ns\n",
			trajMap.size(),earliest,latest);
	}
	/// draws on the current TPad -- successive calls MUST use the same
	/// TPad.
	void draw(double _t) {
		std::map<TrajIdType,Trajectory>::iterator i;
		for(i=trajMap.begin(); i!=trajMap.end(); ++i) {
			TrajIdType id=i->first;
			Trajectory &traj=i->second;
			TMarker *m = markerMap[id];
			if(traj.setT(_t)) {
				if(!m) {
					markerMap[id] = m = 
						particleMarker(pdgIdMap[id]);
					m->Draw();
				}
				m->SetX(traj.v1);
				m->SetY(traj.v2);
			} else {
				if(m) gPad->RecursiveRemove(m);
				delete m;
				markerMap[id] = 0;
			}
		}
	}
	TMarker * particleMarker(int pdgId) {
		// mark=20 (solid circle) slows it down by an enormous factor
		// mark=21 is good
		int mark=setup->marker;

		// smaller marks make color hard to see, and increase artifacts
		float size=setup->markerSize;

		TMarker *m = new TMarker(0.0,0.0,mark);
		m->SetMarkerColor(TColor::GetColor(particle->color(pdgId)));
		m->SetMarkerSize(size);
		return m;
	}
};

/**	the plot command plots a scatter-plot of two expressions in a pad
	of the movie.
**/
class Plot : public MoviePad {
	TString filename;
	TString x;
	TString y;
	TFile *file;
	TrackSet *trackSet;
public:
	Plot() : MoviePad("plot") {
		filename = "";
		x = "x";
		y = "y";
		file = 0;
		trackSet = 0;
	}
	bool defineArgs(TString name, TString value) {
	    if(name=="filename") { filename = value; return true; }
	    if(name=="x") { x = value; return true; }
	    if(name=="y") { y = value; return true; }
	    return MoviePad::defineArgs(name,value);
	}
	void init(Args args) {
		handleArgs(cmdName,args);
		if(filename == "" && gFile != 0 && !gFile->IsZombie() &&
							gFile->IsOpen()) {
			file = gFile;
		} else {
			file = new TFile(filename);
			if(file->IsZombie()) {
				printf("Cannot open '%s'\n",file->GetPath());
				exit(3);
			}
		}
		file->cd();
		TNtuple *nt = (TNtuple *)file->Get("Movie/Tracks");
		if(!nt) {
			printf("cannot find NTuple "
				"Movie/Tracks in Root file\n");
			exit(4);
		}
		trackSet = new TrackSet(nt,x,y);
	}
	void prepareForDrawing() {
		if(title == "_") {
			TString tmp(y);
			tmp.Append(" vs. ").Append(x);
			SetTitle(tmp);
		}
		MoviePad::prepareForDrawing();
	}
	void drawPadContents(double t) {
		trackSet->draw(t);
	}
	ClassDef(Plot,1);
};

/**	the sideview command plots a side view of the apparatus with a moving
	marker indicating the current reference particle location.
**/
class SideView : public MoviePad {
	TString filename;
	TFile *file;
	std::vector<TPolyLine> list;
	TMarker *marker1, *marker2;
	ReferenceTrack *reference;
	float markerSize;
public:
	SideView() : MoviePad("sideview"), list() {
		filename = "";
		file = 0;
		marker1 = marker2 = 0;
		reference = 0;
		markerSize = 1.0;
	}
	bool defineArgs(TString name, TString value) {
	    if(name=="filename") { filename = value; return true; }
	    if(name=="markerSize") { markerSize = value.Atof(); return true; }
	    if(name=="zMin") { xMin=value.Atof(); return true; }
	    if(name=="zMax") { xMax=value.Atof(); return true; }
	    return MoviePad::defineArgs(name,value);
	}
	void init(Args args) {
		handleArgs(cmdName,args);
		if(filename == "" && gFile != 0 && !gFile->IsZombie() &&
							gFile->IsOpen()) {
			file = gFile;
		} else {
			file = new TFile(filename);
			if(file->IsZombie()) {
				printf("Cannot open '%s'\n",file->GetPath());
				exit(3);
			}
		}
		file->cd();
		TNtuple *nt = (TNtuple *)file->Get("Movie/Reference");
		if(!nt) {
			printf("cannot find NTuple "
				"Movie/Reference in Root file\n");
			exit(4);
		}
		reference = new ReferenceTrack(nt);
	}
	void prepareForDrawing() {
		MoviePad::prepareForDrawing();
		file->cd();
		TNtuple *nt = (TNtuple *)file->Get("Movie/Elements");
		if(!nt) {
			printf("cannot find NTuple "
				"Movie/Elements in Root file\n");
			exit(4);
		}
		float Zmin, Zmax, Ymin, Ymax, R, G, B;
		float extra = (xMax - xMin) / (GetWw()-1);
		nt->SetBranchAddress("Zmin",&Zmin);
		nt->SetBranchAddress("Zmax",&Zmax);
		nt->SetBranchAddress("Ymin",&Ymin);
		nt->SetBranchAddress("Ymax",&Ymax);
		nt->SetBranchAddress("R",&R);
		nt->SetBranchAddress("G",&G);
		nt->SetBranchAddress("B",&B);
		long nEntries = nt->GetEntries();
		for(long i=0; i<nEntries; ++i) {
			nt->GetEntry(i);
			float cornerX[5],cornerY[5];
			cornerX[0] = cornerX[1] = cornerX[4] = Zmin-extra;
			cornerX[2] = cornerX[3] = Zmax+extra;
			cornerY[0] = cornerY[3] = cornerY[4] = Ymin;
			cornerY[1] = cornerY[2] = Ymax;
			list.push_back(TPolyLine(5,cornerX,cornerY));
			list.back().SetFillColor(TColor::GetColor(R,G,B));
			list.back().SetFillStyle(4100);
			list.back().SetLineColor(TColor::GetColor(R,G,B));
			list.back().SetLineWidth(1);
		}
		for(unsigned i=0; i<list.size(); ++i) {
			list[i].Draw("f");
		}
		marker1 = new TMarker(0.0,0.0,20);
		marker1->SetMarkerColor(10); // white
		marker1->SetMarkerSize(markerSize);
		marker1->Draw();
		marker2 = new TMarker(0.0,0.0,28);
		marker2->SetMarkerColor(1); // black
		marker2->SetMarkerSize(markerSize);
		marker2->Draw();
	}
	void drawPadContents(double t) {
		reference->setT(t);
		marker1->SetX(reference->Zref);
		marker2->SetX(reference->Zref);
	}
	ClassDef(SideView,1);
};

/**	the space command displays a blank pad in the movie.
**/
class Space : public MoviePad {
public:
	Space() : MoviePad("space") { }
	bool defineArgs(TString name, TString value) 
		{ return MoviePad::defineArgs(name,value); }
	void init(Args args) { (void)args; }
	void prepareForDrawing() { 
		SetTitle("");
		MoviePad::prepareForDrawing();
		TText *text = new TText(0.5,0.5,title);
		text->SetNDC();
		text->SetTextAlign(22);	// centered H and V
		text->SetTextFont(102); // Courier New, high quality
		text->SetTextSize(getTextSize());
		text->Draw();
	}
	void drawPadContents(double t) {
		(void)t; 
	}
	ClassDef(Space,1);
};

/**	the time command displays the current time in a pad in the movie.
**/
class Time : public MoviePad {
	TText *text;
public:
	Time() : MoviePad("time") { }
	bool defineArgs(TString name, TString value) 
		{ return MoviePad::defineArgs(name,value); }
	void init(Args args) { (void)args; }
	void prepareForDrawing() { MoviePad::prepareForDrawing(); 
		text = new TText(0.5,0.5,"");
		text->SetNDC();
		text->SetTextAlign(22);	// centered H and V
		text->SetTextFont(102); // Courier New, high quality
		text->SetTextSize(getTextSize());
		text->Draw();
	}
	void drawPadContents(double t) { 
		text->SetText(0.5,0.5,TString::Format("t=%.3f ns",t));
	}
	ClassDef(Time,1);
};

/**	the position command displays the current Z position in a pad in
	the movie.
**/
class Position : public MoviePad {
	TString filename;
	TFile *file;
	TText *text;
	ReferenceTrack *reference;
public:
	Position() : MoviePad("position") {
		filename = "";
		file = 0;
		text = 0;
		reference = 0;
	}
	bool defineArgs(TString name, TString value) {
	    if(name=="filename") { filename = value; return true; }
	    return MoviePad::defineArgs(name,value);
	}
	void init(Args args) {
		handleArgs(cmdName,args);
		if(filename == "" && gFile != 0 && !gFile->IsZombie() &&
							gFile->IsOpen()) {
			file = gFile;
		} else {
			file = new TFile(filename);
			if(file->IsZombie()) {
				printf("Cannot open '%s'\n",file->GetPath());
				exit(3);
			}
		}
		file->cd();
		TNtuple *nt = (TNtuple *)file->Get("Movie/Reference");
		if(!nt) {
			printf("cannot find NTuple "
				"Movie/Reference in Root file\n");
			exit(4);
		}
		reference = new ReferenceTrack(nt);
	}
	void prepareForDrawing() { MoviePad::prepareForDrawing(); 
		text = new TText(0.5,0.5,"");
		text->SetNDC();
		text->SetTextAlign(22);	// centered H and V
		text->SetTextFont(102); // Courier New, high quality
		text->SetTextSize(getTextSize());
		text->Draw();
	}
	void drawPadContents(double t) { 
		reference->setT(t);
		text->SetText(0.5,0.5,TString::Format("z=%.1f mm",
							reference->Zref));
	}
	ClassDef(Position,1);
};

void Setup::layout() {
	if(pads.size() == 0) {
		printf("Nothing to layout\n");
		exit(1);
	}
	canvas = new TCanvas("Movie", "Movie", windowWidth, windowHeight);
	canvas->SetWindowSize(windowWidth + (windowWidth - canvas->GetWw()),
				windowHeight + (windowHeight - canvas->GetWh()));
	canvas->SetBorderMode(0);
	canvas->SetFillColor(TColor::GetColor(setup->borderColor));

	double totalHeight=0.0;
	for(unsigned i=0; i<rows.size(); ++i) {
		totalHeight += rows[i]->height;
		rows[i]->lastPad = 
			(i+1<rows.size() ? rows[i+1]->firstPad : pads.size())-1;
	}
	for(unsigned i=0; i<rows.size(); ++i) {
		rows[i]->pixelHeight = windowHeight*rows[i]->height/totalHeight;
		printf("Row %d: %d pads   %d pixels high\n",i+1,
		    rows[i]->lastPad-rows[i]->firstPad+1,rows[i]->pixelHeight);
	}
	unsigned row=0;
	double height=rows[row]->height/totalHeight;
	double bottom=1.0-height;

	double w=1.0/(double)(rows[row]->lastPad-rows[row]->firstPad+1);
	double x=0.0;
	for(unsigned i=0; i<pads.size(); ++i) {
		if(row+1 < rows.size() && (int)i > rows[row]->lastPad) {
		    ++row;
		    height=rows[row]->height/totalHeight;
		    bottom -= height;
		    w=1.0/(double)(rows[row]->lastPad-rows[row]->firstPad+1);
		    x = 0.0;
		}
		canvas->cd();
		pads[i]->Draw();
		pads[i]->SetPad(x,bottom,x+w,bottom+height);
		pads[i]->prepareForDrawing();
		x += w;
	}
}

void Setup::drawOneFrame(double t) {
	for(unsigned i=0; i<pads.size(); ++i) {
		pads[i]->cd();
		pads[i]->drawPadContents(t);
	}
}

void Setup::drawMovie()
{
	int nFrames = duration*frameRate;
	double dt=(tMax-tMin)/(double)nFrames;
	canvas->cd();
	for(int i=0; i<=nFrames; ++i) {
		double t = tMin + i*dt;
		printf("frame %d/%d  t=%.3f ns\r",i,nFrames,t); fflush(stdout);
		drawOneFrame(t);
		canvas->cd();
		char tmp[64];
		sprintf(tmp,"movie_%d.%s",i,pictureType.Data());
		canvas->Update();
		canvas->SaveAs(tmp);
	}
	printf("\n");
}

void Command::doCommand(TString cmd, Args args) {
	if(cmd == "setup") {
		setup->handleArgs("setup",args);
		return;
	}
	if(cmd == "particle") {
		particle->handleArgs("particle",args);
		return;
	}
	if(cmd == "pad") {
		MoviePad *p = new MoviePad("pad");
		p->init(args);
		return;
	}
	if(cmd == "plot") {
		Plot *p = new Plot();
		p->init(args);
		return;
	}
	if(cmd == "sideview") {
		SideView *p = new SideView();
		p->init(args);
		return;
	}
	if(cmd == "row") {
		setup->rows.push_back(new Row(args,setup->pads.size()));
		return;
	}
	if(cmd == "space") {
		(new Space())->handleArgs("space",args);
		return;
	}
	if(cmd == "time") {
		(new Time())->handleArgs("time",args);
		return;
	}
	if(cmd == "position") {
		(new Position())->init(args);
		return;
	}

	// unknown command - just print it and its args
	printf("Unknown cmd=%s\n",cmd.Data());
	for(Args::iterator i=args.begin(); i!=args.end(); ++i) {
		TString name=i->first;
		TString value=i->second;
		printf("    %s=%s\n",name.Data(),value.Data());
	}
	++nError;
}

///	G4blMovie() -- the macro function that gets called to create the movie.
void G4blMovie() {
	gROOT->GetStyle("Video");

	setup = new Setup();
	particle = new Particle();

	// the g4blmovie script puts its $1 into the environment MOVIE_IN
	TString movieIn=gSystem->Getenv("MOVIE_IN");
	printf("\nG4blMovie starting: Reading '%s'\n",movieIn.Data());
	Command::read(movieIn);

	if(nError > 0) {
		printf("%d errors found\n",nError);
		exit(1);
	}

	printf("\nRemoving old movie_*.%s and output file '%s'\n",
			setup->pictureType.Data(),setup->outputFile.Data());
	TString cmd="rm -f movie_*.";
	cmd.Append(setup->pictureType).Append(" ").Append(setup->outputFile);
	gSystem->Exec(cmd);

	printf("\nLaying out the pads into the movie\n");
	setup->layout();

	firstpad = setup->pads[0];

	printf("\nDrawing movie\n");
	setup->drawMovie();

	cmd = "ffmpeg -an -r ";
	cmd.Append(TString::Format("%.0f",setup->frameRate)).
			Append(" -i movie_%d.").Append(setup->pictureType).
			Append(" ").Append(setup->outputFile);
	printf("\nRunning %s\n",cmd.Data());
	gSystem->Exec(cmd);

	printf("\nRemoving old movie_*.%s\n",setup->pictureType.Data());
	cmd="rm -f movie_*.";
	cmd.Append(setup->pictureType);
	gSystem->Exec(cmd);
	printf("G4blMovie exiting.\n");
}
