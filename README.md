[Setup]

source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh

This repository aims to perform Tag and Probe studies starting from the code available here 
Bntuple: https://github.com/boundino/Bntuple.
You can find below the steps you have to follow to perform the full study.

1) Code/fitJPsi.C 

**Input**: This code gets as inputs ntuples in which all the variables related to the jspi candidates are 
           stored: e.g. Jpsi mass, pt, eta and all the variables related to the two muons used to build the 
           Jpsi candidates  (e.g. pt1, pt2). In particular all the variables used to define 
           tag, passing and passing probe are defined for the first and the second muon.

