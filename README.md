
This repository aims to perform Tag and Probe studies starting from the code available here 
Bntuple: https://github.com/boundino/Bntuple.
You can find below the steps you have to follow to perform the full study.

**Step 1) Code/fitJPsi.C** 

**Input**: This code gets as inputs ntuples in which all the variables related to the jspi candidates are 
           stored: e.g. Jpsi mass, pt, eta and all the variables related to the two muons used to build the 
           Jpsi candidates  (e.g. pt1, pt2). In addition all the variables used to define 
           tag, passing and passing probe are also stored in the input  for the first and the second muon.

**What it does**: it loops over the jspi candidates and it fill histograms of passing, failing and all probes 
                  in different pt and eta bins for the trigger, tracking, muonID studies.

**Output** : the output is a file called foutputData.root or foutputMC.root in which all the histos are stored as well 
             as TTree for trg, trk and MuonID separately that can be also used for fitting pourpose
             (the current code we are developing uses histograms and not TTree)

**How to run?** :

                cd Code
                root 
                .L fitJPsi.C+
                fitJPsi(true) for running on data
                fitJPsi(false) for running on MC

The input file name are UP TO DATE!
VERY IMPORTANT: YOU HAVE TO DEFINE WHERE TO PUT YOUR OUTPUT IN line 18 and 19.
              
**Step 2) Display/FitTnP_sample.C** 

**Input**: output of fitJPsi.C, please REMIND to update the input files name in the code. 
           If not, it will use as an input /afs/cern.ch/user/g/ginnocen/public/TnPResults/foutputData.root
           that is the current result for data. 
           (the current result for MC is /afs/cern.ch/user/g/ginnocen/public/TnPResults/foutputMC.root)

**What it does**: it does the simultaneous fit of passing, failing, all probes for Trg, Trk, MuonId 
                  studies vs pt and eta and it produces corresponding efficiencies.

**How to run?** MAKE sure of the input file names are ok. 
                then launch :
                
                cd Display
                source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh
                source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh

                root 
                .x FitTnP_sample.C
