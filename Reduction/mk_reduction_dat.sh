#!/bin/bash

. setup.sh
rm TechProp_2014/*.dat
./AnalysisScripts/mk_reduction_dat.py  /store/group/phys_higgs/cmshgg/processed/FixPFIsolationIssue/22AprFix/    ${storedir}/TechProp_2014   TechProp_2014.txt

wd=$PWD
cd AnalysisScripts
tar cf $wd/${version}.tar *.py $(find common reduction baseline massfac_mva_binned full_mva_binned jetanalysis photonjet -name \*.dat \
-or -name \*.py) aux common python
cd -

tar rf ${version}.tar JSON *.sh
gzip -f ${version}.tar


git tag -a ${version} -m "Tag used for reduction ${group}/${version}"
