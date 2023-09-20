#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

using namespace std;

gridCell::gridCell(){}

gridCell::gridCell(int x, int y, int z, int index){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  id=index;
  innervation=1.0;
  collagen1=0.0;
  collagen3=0.0;
  viableMuscle=1.0;
  woundEdge=0;
  colEdge-0;
  satellitePool=false;
  life=1;
  ec_roll=0;
  ec_stick=0;
  necrosis=0;
  ec_activation=0;
  cellPopulation=0;
  vascularization=0;
  ec_migrate=0;
  contamination=0;
  tempSignal=0;
  counter=0;
}

endoCell::endoCell(){}

endoCell::endoCell(int x, int y, int z, int index){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  id=index;
}

epiCell::epiCell(){}

epiCell::epiCell(int x, int y, int z, int index){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  id=index;
}

satelliteSource::satelliteSource(){}

satelliteSource::satelliteSource(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  activation=0;
}


pmnSource::pmnSource(){}

pmnSource::pmnSource(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  activation=0;
}

macroSource::macroSource(){}

macroSource::macroSource(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  activation=0;
}

T_Source::T_Source(){}

T_Source::T_Source(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
}

fibroSource::fibroSource(){}

fibroSource::fibroSource(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  activation=0;
}


SimulationObject::SimulationObject(){
  xDimension=xDim;
  yDimension=yDim;
}

SimulationObject::SimulationObject(float internalParam[]){
  for (int i = 0; i < 1631; i++) {
    internalParameterization[i] = internalParam[i];
  }
  xDimension=xDim;
  yDimension=yDim;
};
