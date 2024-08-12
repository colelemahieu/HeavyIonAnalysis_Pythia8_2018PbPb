#!/bin/bash

dirR2="MinBias_PbPb_5p36TeV_Hydjet_mini_R2/240312_211736"
dirR4="MinBias_PbPb_5p36TeV_Hydjet_mini_R4_Run2/240805_174123"
#dirR4="MinBias_PbPb_5p36TeV_Hydjet_mini/231208_154825"
dirR6="MinBias_PbPb_5p36TeV_Hydjet_mini_R6/240312_181630"
dirR8="MinBias_PbPb_5p36TeV_Hydjet_mini_R8/240312_145733"

if [ ${1} -lt 10 ] 
then
	path="/eos/cms/store/group/phys_heavyions/clemahie/Forest_test/pythiaForest/${dirR4}/000${1}"
        #path="/eos/cms/store/group/phys_heavyions/clemahie/Forest/pythiaGenSim10000_2/${dirR4}/000${1}"
else
	path="/eos/cms/store/group/phys_heavyions/clemahie/Forest_test/pythiaForest/${dirR4}/00${1}"
        #path="/eos/cms/store/group/phys_heavyions/clemahie/Forest/pythiaGenSim10000_2/${dirR4}/000${1}"
fi


outPath="/eos/cms/store/group/phys_heavyions/clemahie/pythiaGen_jets_${1}.root"

if test -d "${path}"; then
    root -b -l -q "/afs/cern.ch/user/c/clemahie/MC_MINIAOD_Forest/CMSSW_11_2_1_patch2/src/HeavyIonsAnalysis/Configuration/test/miniAOD_genSkim.C(\"${path}/HiForestMiniAOD_*.root\", \"${outPath}\")" 
fi





