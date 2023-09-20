#ifndef WEABM_PARAMETERS_H
#define WEABM_PARAMETERS_H

using namespace std;

const int numTimeSteps=4000;

// const int xDim=251;
// const int yDim=3;
// const int zDim=100;

const int xDim=30;
const int yDim=30;
const int zDim=200;

const float fibroAgeThreshold=100;
const float fibroCol1Deposition=0.5;
const float fibroCol3Deposition=0.5;
const float myoFibroThreshold=150;
const float fibroActivationThreshold=0.1;
const float fibro_tnf_apoptosisFactor=0.1;
const float fibro_tgfb_apoptosisFactor=0.05;
const int fibrosPerStep=5;
const float fibro_tgfb_proliferationFactor=0.005;
const float macroNecroClearance=0.1;

const int pmnProliferationActivationThreshold=200;
const int macroProliferationActivationThreshold=400;

const float pmn_GCSF_proliferationFactor=0.05;
const int pmnsPerStep=5;

const int Th0sPerStep=5;

const int cellCapacity=5;

const int macrosPerStep=5;

const float k_myoblast_fusion_addition=0.5;

const float minimumCytokineValue=0.0001;

const float oxyHeal=0.1;

const float vascProlif=1.0;

const int outputWrite=1;
const int tissueWrite=1;

const string filename1="/home/chase/Desktop/Research/BETR/Matlab/Output1.csv";
const string filename2="/home/chase/Desktop/Research/BETR/Matlab/Output2.csv";

const float hgf_ssc_prolif_factor=1;
const float ssc_prolif_chance=0.2;
const float ssc_tgfb_age_factor=0.1;

const int numRules=53;
const int numRuleParams=31;

const float mrf_critValue=10;

const string ctInputFile="DogTest.csv";


#endif
