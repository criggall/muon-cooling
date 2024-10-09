//	BLUserCode.cc
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

#include "G4ParticleTable.hh"

#include "BLManager.hh"
#include "BLUserCode.hh"

void BLUserCode::registerUserCode()
{
	BLManager::getObject()->registerUserCode(this);
}

int BLUserTrackFilter::getPDGid(G4Track *track)
{
	return track->GetDefinition()->GetPDGEncoding();
}

void BLUserTrackFilter::setMomentum(G4Track *track, G4ThreeVector momentum)
{
	G4double mass=track->GetDefinition()->GetPDGMass();
	G4double KE=sqrt(momentum.mag2()+mass*mass) - mass;
	track->SetMomentumDirection(momentum.unit());
	track->SetKineticEnergy(KE);
}


G4double BLUserTrackFilter::getMass(int PDGid)
{
	return G4ParticleTable::GetParticleTable()->FindParticle(PDGid)->
								GetPDGMass();
}

void BLUserTrackFilter::killTrack(G4Track *track)
{
	track->SetTrackStatus(fStopAndKill);
}

G4Track *BLUserTrackFilter::constructTrack(G4ThreeVector position,
	G4ThreeVector momentum, G4double time, int PDGid, G4double weight)
{
	G4ParticleDefinition *particle = 
		    G4ParticleTable::GetParticleTable()->FindParticle(PDGid);
	G4DynamicParticle *dyn = new G4DynamicParticle(particle,momentum);
	G4Track *track = new G4Track(dyn,time,position);
	track->SetWeight(weight);
	return track;
}

