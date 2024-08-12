#!/bin/bash

dirR2="MinBias_PbPb_5p36TeV_Hydjet_mini_R2/240312_211736"
dirR4="MinBias_PbPb_5p36TeV_Hydjet_mini_R4_Run2/240805_174123"
dirR6="MinBias_PbPb_5p36TeV_Hydjet_mini_R6/240312_181630"
dirR8="MinBias_PbPb_5p36TeV_Hydjet_mini_R8/240312_145733"

if [ ${1} -lt 10 ] 
then
	path="/eos/cms/store/group/phys_heavyions/clemahie/Forest_test/pythiaForest/${dirR4}/000${1}"
        echo "made it here"
else
	path="/eos/cms/store/group/phys_heavyions/clemahie/Forest_test/pythiaDigi/${dirR4}/00${1}"
fi


outPath="/eos/cms/store/group/phys_heavyions/clemahie/pythiaReco_miniAOD_${1}.root"

if test -d "${path}"; then
    echo "made it to running the script"
    root -b -l -q "/afs/cern.ch/user/c/clemahie/MC_MINIAOD_Forest/CMSSW_11_2_1_patch2/src/HeavyIonsAnalysis/Configuration/test/miniAOD_mcSkim.C(\"${path}/HiForestMiniAOD_*.root\", \"${outPath}\")" 
fi





