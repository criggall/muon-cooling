#!/bin/bash
#	triplet.sh - demonstrate the tuning of a quad triplet
#

###
### check for required programs
###
EXIT=0
if ! type g4bl >/dev/null 2>&1
then
	echo The program g4beamline is required.
	echo It is available in Computer Programs at http://muonsinc.com .
	EXIT=1
fi
if ! type gminuit >/dev/null 2>&1
then
	echo The program gminuit is required.
	echo It is available in Computer Programs at http://muonsinc.com .
	EXIT=1
fi
if ! type wish >/dev/null 2>&1
then
	echo You need to have your system administrator install Tcl/Tk
	EXIT=1
fi
if ! type tclsh >/dev/null 2>&1
then
	echo You need to have your system administrator install Tcl/Tk
	EXIT=1
fi
if ! type gnuplot >/dev/null 2>&1
then
	echo You need to have your system administrator install gnuplot
	EXIT=1
fi
if test $EXIT != 0; then exit $EXIT; fi

###
### define temp files used
###
GNUPLOT=/tmp/gnuplot.$$
G4BL=/tmp/g4bl.$$
GMINUIT=/tmp/gminuit.$$
SCRIPT=/tmp/script.$$
PROFILE=/tmp/profile.$$
APERTURE=/tmp/aperture.$$
OUT=/tmp/out.$$
trap "kill %1 >/dev/null 2>&1; rm -f $GNUPLOT $G4BL $GMINUIT \
		$SCRIPT $PROFILE $APERTURE $OUT" 0 1 2 3

###
### start gnuplot from a pipe, so a single plot window is used
###
mkfifo $GNUPLOT
gnuplot <$GNUPLOT >/dev/null &
exec 3>$GNUPLOT   # keep the pipe open as long as we are alive

cat <<-! >$G4BL
# command-line parameters: P dP sigmaX sigmaY sigmaXp sigmaYp QF QD
physics QGSP_BIC
trackcuts keep=mu+
param worldMaterial=Vacuum
beam gaussian particle=mu+ nEvents=500 meanMomentum=\$P sigmaP=\$dP \
	sigmaX=\$sigmaX sigmaY=\$sigmaY sigmaXp=\$sigmaXp sigmaYp=\$sigmaYp
genericquad Quad fieldLength=100 ironLength=100 apertureRadius=100 \
	ironRadius=200 kill=1
box Box height=50 width=50 length=1
place Box z=0 rename=Beam color=0,1,0
place Quad z=2000 rename=Q1 gradient=\$QF
place Quad z=3000 rename=Q2 gradient=\$QD
place Quad z=4000 rename=Q3 gradient=\$QF
place Box z=6000 rename=End color=1,0,0
profile file=$PROFILE particle=mu+ z=0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,5900,6000,6100,6200
!

cat <<-! >$GMINUIT
param P 200.0 100.0 1000.0 0 fixed
param dP 0.0 0.0 100.0 0 fixed
param sigmaX 0.0 0.0 100 0 fixed
param sigmaY 0.0 0.0 100 0 fixed
param sigmaXp 0.0065 0.0 0.1 0 fixed
param sigmaYp 0.0065 0.0 0.1 0 fixed
# values for fringe=0: QF=5.912 QD=-7.933
# values for fringe enabled: QF=5.375 QD=-7.235
param QF 5.375 -15 15 0.5 limited
param QD -7.235 -15 15 0.5 limited
script $SCRIPT value 1
minuitcmd MIGRAD
!

cat <<-! >$APERTURE
1950 150
1950 100
2050 100
2050 150

2950 150
2950 100
3050 100
3050 150

3950 150
3950 100
4050 100
4050 150

1950 -150
1950 -100
2050 -100
2050 -150

2950 -150
2950 -100
3050 -100
3050 -150

3950 -150
3950 -100
4050 -100
4050 -150

0 0
6200 0
!

cat <<-! >$SCRIPT
echo >&2 Executing: $G4BL "\$@"
g4bl $G4BL "\$@" >$OUT 2>&1
echo >$GNUPLOT 'plot [0:6200][-120:120] "$PROFILE" using 1:(3*\$4) with lines title "3*sigmaX", "$APERTURE" with lines notitle, "$PROFILE" using 1:(-3*\$6) with lines title "-3*sigmaY"'
# chisq = sum of sigmaX^2 and sigmaY^2 at Z=6000
V=\`awk <$PROFILE '/^6000/ {print \$4*\$4+\$6*\$6}'\`
echo >&2 "Chisq=\$V   \$@"
echo "\$V"
!

chmod +x $SCRIPT
gminuit $GMINUIT

# omit the line indicating that gnuplot terminated
exec 1>/dev/null 2>&1
