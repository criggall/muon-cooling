//	BLMPI.hh -- MPI interface for parallelization
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

#ifndef BLMPI_HH
#define BLMPI_HH

#include <vector>
#include <globals.hh>

/**	class BLMPI implements MPI interface for parallelization.
 *
 *	BLMPI handles G4beamline-specific MPI stuff.
 *
 *	NOTE: No other class should call any MPI functions directly. This
 *	implies that the only #include <mpi.h> should be in BLMPI.cc.
 *
 *	NOTE: When running in MPI mode, threads are not permitted.
 *
 *	When running under MPI, all BLNTuples are physically written
 *	only by the rank0 instance; MPI is used by all other instances to
 *	send messages to rank0 to create BLNTuple-s and to append rows to
 *	BLNTuple-s.
 *
 *	While it sounds simple to have rank 0 parcel out events to the workers
 *	and collect NTuple data, in practice it is rather complex to get
 *	it to not waste time and closeup properly. A simple implementation
 *	can easily block in rank 0 sending to rank 1, while rank 1 is
 *	simulating events -- then rank 0 cannot service other ranks' requests.
 *	So there is a BeamTracksMsg for each worker, and the send to the
 *	worker is asynchronous.
 *
 *	This entire class is static.
 **/
class BLMPI {
	static int rank;
	static int nRanks;
	static int alarmTime;
	static bool processingFatalException;
	static void rank0AlarmHandler();
	static void mainRankZero();
	static void mainRankNonZero();
	static bool readyForBeam;
	static int nExiting;
	static bool inScan;
	static int inMPI;
	static double startupTime;
	static double startupWaitTime;
	static double computeTime;
	static double computeWaitTime;
public:
	/// init() will check whether or not we are running under MPI,
	/// and will initialize MPI if so. 
	/// MUST be called before any other BLMPI routine.
	/// NOTE: for all but rank0, sets stdout to /dev/null (stderr
	/// remains connected for all ranks.
	static void init(std::vector<G4String> &argV);

	/// isMPI() returns true if in MPI mode.
	static bool isMPI();

	/// isRank0() returns false in non-MPI mode, and returns true if
	/// this node is rank 0 in MPI.
	static bool isRank0() { return rank == 0; }

	/// isRank1() returns true in non-MPI mode, and returns true if
	/// this node is rank 1 in MPI. Used for things that should happen
	/// once per simulation, not for every worker.
	static bool isRank1() { return rank == 1; }

	/// getRank() returns the rank of this process.
	static int getRank() { return rank; }

	/// getNranks() returns the number of ranks. The number of workers
	/// is nRanks-1.
	static int getNranks() { return nRanks; }

	/// scan() will scan for and handle an incoming MPI message.
	/// If rank!=0 it returns immediately, unless an ABORT_MSG is
	/// received, in which case it calls closeupAndExit().
	/// In rank0:
	/// if wait=false, it will handle 0 or more messages and return.
	/// if wait=true, it handles all messages until all ExitingMsg
	/// messages have been received.
	static void scan(bool wait=false);

	/// sendExceptionMessage() will send an exception message to rank 0.
	/// Rank0 keeps track of all of them across all ranks.
	static void sendExceptionMessage(const char *origin, const char *code,
		const char *severity, const char *description,
		int eventID, int trackID, const char *particleName,
		G4double kineticEnergy, bool abortProgram);

	/// rank0BufferPrintf() will do a printf() in non-MPI mode or on rank
	/// 0, via a message on other ranks; the data are buffered until
	/// rank0BufferFlush() is called, or the buffer is full. The buffer
	/// is 32768 chars; each call is limited to 1024 chars.
	/// filename="-" prints to rank 0 stdout.
	static void rank0BufferPrintf(G4String filename, const char *format,
									...);

	/// rank0BufferFlush() will flush all buffers to rank 0.
	static void rank0BufferFlush();

	/// sumAllWorkers() will sum up values over all workers. Used in 
	/// after-beam callbacks. On rank 0, will set data[] to the sum from 
	/// ranks 1-N. On ranks 1-N, will send a SUM message to rank 0 with 
	/// the data. This is essentially a Reduce-sum operation, except that 
	/// no global synchronization is involved. This relies on the fact
	/// that calls to sumAllWorkers() will occur in the same order in
	/// all ranks, including 0. Note the 
	/// closeoutAndExit sequencing ensures that ranks 1-N do the callbacks 
	/// and exit, before rank 0 does the callbacks. Does nothing but return 
	/// true if not in MPI mode. Used in BLCMDprofile and BLCMDtotalenergy,
	/// and any other places that need a sum over all ranks.
	/// Returns true if data[] contains the sum (rank 0 or non-MPI mode).
	static bool sumAllWorkers(double data[], int nData);

	/// main will perform the computation in MPI mode.
	/// In non-MPI mode this is a no-op that returns immediately.
	/// In MPI mode this determines the collective/serial mode and then
	/// performs the computation, never returning.
	static void main();

	/// closeupAndExit() will do just that, for all ranks. This MUST be the
	/// final call, as it sequences the closeup, and then does
	/// MPI_Finalize() and _exit().
	static void closeupAndExit(int value);
};

#endif // BLMPI_HH
