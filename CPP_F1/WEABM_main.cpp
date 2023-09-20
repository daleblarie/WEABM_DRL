#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

#include<chrono>

using namespace std;

// int xDim=1000; Fibers oriented along this axis
// int yDim=400;
// int zDim=300;



uniform_int_distribution<int>  distributionEdges(0,5);
uniform_int_distribution<int>  distributionX(0,xDim);
uniform_int_distribution<int>  distributionY(0,yDim);
uniform_int_distribution<int>  distributionZ(0,zDim);
uniform_int_distribution<int>  distribution3(0,2);
uniform_int_distribution<int>  distribution27(0,26);
uniform_real_distribution<float>  dist01(0.0,1.0);

// typedef std::chrono::high_resolution_clock myclock;
//   myclock::time_point beginning = myclock::now();


// SimulationObject::SimulationObject(int rseed, float* internalParameterization){
//   generator.seed(rseed);
//   outputSetup();
//   getRuleMatrix(internalParameterization,1260);
//   initialize();
//   cout<<"Initialization Complete "<<fibroCol3Deposition<<"\n";
//   injure_biopsy();
//   cout<<"Injury Applied\n";
//   setupBiopsySatellites();
// }

void SimulationObject::setupSimulation(int rseed){
  cout<<"RSEED="<<rseed<<"\n";
    generator.seed(rseed);
    outputSetup();
    getRuleMatrix(internalParameterization,1631);
    initialize();
    //initializeFromCT();
    cout<<"Initialization Complete "<<fibroCol3Deposition<<"\n";
    injure_biopsy();
    cout<<"Injury Applied\n";
    setupBiopsySatellites();
}

float* SimulationObject::sim_obj_main(){
  int step,i,j,k;
  for(step=0;step<numTimeSteps;step++){
    cout<<step<<" "<<bursts<<"\n";
    bursts=0;
    cellFunctions();
    evaporations();
    diffusions();
    gridOutput(step);
    aggregateData(step);
  }
  outputFinalize();
  k=0;
  for(j=0;j<43;j++){
    for(i=0;i<numTimeSteps;i++){
      allSignalsReturn[k]=allSignals[j][i];
      k++;
    }
  }
  cout<<"All Done\n";
  return allSignalsReturn;
}

void SimulationObject::singleStep(){
  cellFunctions();
  evaporations();
  diffusions();
  gridOutput(current_step);
  aggregateData(current_step);
}

float* SimulationObject::endSimulation(){
  outputFinalize();
  int i,j, k=0;
  for(j=0;j<43;j++){
    for(i=0;i<numTimeSteps;i++){
      allSignalsReturn[k]=allSignals[j][i];
      k++;
    }
  }
  cout<<"All Done\n";
  return allSignalsReturn;
}


void SimulationObject::cellFunctions(){
  int i,j,tx,ty,tz,tempid,length;
  current_step++;

  length=pmns.size();
  if(length>0){
    shuffle(pmns.begin(),pmns.end(),generator);
  }
  j=0;
  while(j<length){
    pmns[j].pmn_function(j, tissueGrid, cytoGrids, pmns, bursts, generator, RM);
    j++;
    length=pmns.size();
  }
  //  cout<<"Test 80\n";

  length=macrophages.size();
  if(length>0){
    shuffle(macrophages.begin(),macrophages.end(),generator);
  }
  j=0;
  while(j<length){
    macrophages[j].action(j, tissueGrid, cytoGrids, macrophages, generator, RM);
    j++;
    length=macrophages.size();
  }

  length=myoblasts.size();
  if(length>0){
    shuffle(myoblasts.begin(),myoblasts.end(),generator);
  }
  j=0;
  while(j<length){
    myoblasts[j].action(j, tissueGrid, cytoGrids, myoblasts,RM,generator);
    j++;
    length=myoblasts.size();
  }

  length=satellites.size();
  if(length>0){
    shuffle(satellites.begin(),satellites.end(),generator);
  }
  j=0;
  while(j<length){
    satellites[j].action(j, tissueGrid, cytoGrids, satellites, myoblasts);
    j++;
    length=satellites.size();
  }



    shuffle(gridIndexes.begin(),gridIndexes.end(),generator);
    for(i=0;i<gridIndexes.size();i++){
      tempid=gridIndexes[i];
      tz=floor(tempid/(xDim*yDim));
      ty=floor(tempid-tz*xDim*yDim)/xDim;
      tx=tempid-ty*xDim-tz*xDim*yDim;
      tissueGrid[tx][ty][tz].endothelialFunction(tissueGrid, cytoGrids, pmnSources, macroSources, T_Sources, fibroSources, satelliteSources,generator,RM);
    }

    length=fibroblasts.size();
    if(length>0){
      shuffle(fibroblasts.begin(),fibroblasts.end(),generator);}
      j=0;
      while(j<length){
        fibroblasts[j].action(j, tissueGrid, cytoGrids, fibroblasts, generator, current_step,RM);
        j++;
        length=fibroblasts.size();
      }






    //
    length=pmnSources.size();
    j=0;
    if(length>0){
      while(j<length){
        //    cout<<"Triggering PMN Source\n";
        pmnSources[j].proliferate(tissueGrid, cytoGrids, pmns, generator);
        j++;
        length=pmnSources.size();
      }
    }
    // cout<<"Test 220\n";
    //
    length=macroSources.size();
    j=0;
    if(length>0){
      while(j<length){
        macroSources[j].proliferate(tissueGrid, cytoGrids, macrophages, generator);
        j++;
        length=macroSources.size();
      }
    }

    length=fibroSources.size();
    j=0;
    if(length>0){
      while(j<length){
        fibroSources[j].proliferate(tissueGrid, cytoGrids, fibroblasts, generator);
        j++;
        length=fibroSources.size();
      }
    }

    length=satelliteSources.size();
    j=0;
    if(length>0){
      while(j<length){
        satelliteSources[j].proliferate(tissueGrid, satellites, myoblasts,generator, cytoGrids, TotalTNF);
        j++;
        length=satelliteSources.size();
      }
    }
    //
    length=T_Sources.size();
    for(j=0;j<length;j++){
      T_Sources[j].proliferate(tissueGrid, Th0s, generator, TotalTGFb);
    }
    //
    length=Th0s.size();
    if(length>0){
      shuffle(Th0s.begin(),Th0s.end(),generator);}
      j=0;
      while(j<length){
        Th0s[j].action(j, tissueGrid, cytoGrids, Th0s, Th1s, Th2s, Th17s, generator,RM);
        j++;
        length=Th0s.size();
      }

      // length=Th1s.size();
      // if(length>0){
      //   shuffle(Th1s.begin(),Th1s.end(),generator);}
      // j=0;
      // while(j<length){
      //   Th1s[j].action(j);
      //   j++;
      //   length=Th1s.size();
      // }
      // //
      // length=Th2s.size();
      // if(length>0){
      //   shuffle(Th2s.begin(),Th2s.end(),generator);}
      // j=0;
      // while(j<length){
      //   Th2s[j].action(j);
      //   j++;
      //   length=Th2s.size();
      // }
      //  cout<<"Test 160\n";
      //
      length=Th17s.size();
      if(length>0){
        shuffle(Th17s.begin(),Th17s.end(),generator);}
        j=0;
        while(j<length){
          Th17s[j].action(j, tissueGrid, cytoGrids, Th17s, generator,RM);
          j++;
          length=Th17s.size();
        }
      }


// Simulation Functions that are available outside of c++ via ctypes
extern "C" {
  void* simulationObjectAddress;
  SimulationObject* sim_obj;
  float* surfaceCytokineLevels;

  void* CreateInstance(float* internalParam){
    delete sim_obj;
    SimulationObject* sim_obj = new SimulationObject(internalParam);
    simulationObjectAddress = sim_obj;
    return simulationObjectAddress;
  }
  void InitializeSimulation(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->outputSetup();
    ref->getRuleMatrix(ref->internalParameterization,1631);
    ref->initialize();
//    ref->initializeFromCT();
    ref->updateOutputGrids();
    ref->getTopGrid();
  }
  void InitializeSimulationCT(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->outputSetup();
    ref->getRuleMatrix(ref->internalParameterization,1631);
    ref->initializeFromCT();
    ref->updateOutputGrids();
    ref->getTopGrid();
  }
  void DoInjury(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->injure_biopsy();
    ref->setupBiopsySatellites();
    ref->updateOutputGrids();
    ref->getTopGrid();
    ref->aggregateData(0);
  }
  void StepSimulation(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->cellFunctions();
    ref->evaporations();
    ref->diffusions();
    ref->gridOutput(ref->current_step);
    ref->aggregateData(ref->current_step);
    ref->getTopGrid();
    ref->updateOutputGrids();
    // cout<<ref->current_step<<endl;
  }

  float* getSurfaceCytokineLevels(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getSurfaceCytokines();
  }
  float* getTotalCytokineLevels(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getTotalCytokines();
  }
  float* getInnervation(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getInnervation();
  }
  float* getCollagen1(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getCollagen1();
  }
  float* getCollagen3(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getCollagen3();
  }
  float* getViableMuscle(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getViableMuscle();
  }
  float* getVascularization(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getVascularization();
  }
  float* getLife(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getLife();
  }
  float* getNecrosis(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getNecrosis();
  }
  float* getContamination(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->getContamination();
  }

  float* get_IL8(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL8();
  }
  float* get_MIP1b(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_MIP1b();
  }
  float* get_TNF(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_TNF();
  }
  float* get_HGF(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_HGF();
  }
  float* get_PAF(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_PAF();
  }
  float* get_IL1(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL1();
  }
  float* get_PDGF(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_PDGF();
  }
  float* get_IL13(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL13();
  }
  float* get_TGFb(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_TGFb();
  }
  float* get_VEGF(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_VEGF();
  }
  float* get_GCSF(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_GCSF();
  }
  float* get_IL10(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL10();
  }
  float* get_IL4(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL4();
  }
  float* get_IL17(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL17();
  }
  float* get_Myf5(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_Myf5();
  }
  float* get_IL6(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL6();
  }
  float* get_IL1ra(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL1ra();
  }
  float* get_sIL1r(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_sIL1r();
  }
  float* get_sTNFr(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_sTNFr();
  }
  float* get_endotoxin(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_endotoxin();
  }
  float* get_IFNg(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IFNg();
  }
  float* get_IL12(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IL12();
  }
  float* get_MRF4(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_MRF4();
  }
  float* get_cytotox(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_cytotox();
  }
  float* get_IGF(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_IGF();
  }
  float* get_DAMP(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_DAMP();
  }
  float* get_MCP1(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_MCP1();
  }
  float* get_FNE(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->get_FNE();
  }

  void apply_IL8(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL8(addedAmount);
  }
  void apply_MIP1b(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_MIP1b(addedAmount);
  }
  void apply_TNF(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_TNF(addedAmount);
  }
  void apply_HGF(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_HGF(addedAmount);
  }
  void apply_PAF(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_PAF(addedAmount);
  }
  void apply_IL1(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL1(addedAmount);
  }
  void apply_PDGF(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_PDGF(addedAmount);
  }
  void apply_IL13(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL13(addedAmount);
  }
  void apply_TGFb(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_TGFb(addedAmount);
  }
  void apply_VEGF(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_VEGF(addedAmount);
  }
  void apply_GCSF(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_GCSF(addedAmount);
  }
  void apply_IL10(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL10(addedAmount);
  }
  void apply_IL4(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL4(addedAmount);
  }
  void apply_IL17(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL17(addedAmount);
  }
  void apply_Myf5(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_Myf5(addedAmount);
  }
  void apply_IL6(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL6(addedAmount);
  }
  void apply_IL1ra(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL1ra(addedAmount);
  }
  void apply_sIL1r(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_sIL1r(addedAmount);
  }
  void apply_sTNFr(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_sTNFr(addedAmount);
  }
  void apply_endotoxin(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_endotoxin(addedAmount);
  }
  void apply_IFNg(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IFNg(addedAmount);
  }
  void apply_IL12(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IL12(addedAmount);
  }
  void apply_MRF4(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_MRF4(addedAmount);
  }
  void apply_cytotox(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_cytotox(addedAmount);
  }
  void apply_IGF(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_IGF(addedAmount);
  }
  void apply_DAMP(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_DAMP(addedAmount);
  }
  void apply_MCP1(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_MCP1(addedAmount);
  }
  void apply_FNE(void* ptr, float addedAmount){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    ref->apply_FNE(addedAmount);
  }
  int getNumPMN(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->pmns.size();
  }
  int getNumMacrophage(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->macrophages.size();
  }
  int getNumFibroblast(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->fibroblasts.size();
  }
  int getNumMyoblast(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->myoblasts.size();
  }
  int getNumSatellite(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->satellites.size();
  }
  int getNumTH0(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->Th0s.size();
  }
  int getNumTH1(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->Th1s.size();
  }
  int getNumTH2(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->Th2s.size();
  }
  int getNumTH17(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->Th17s.size();
  }
  int getNumPMNsource(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->pmnSources.size();
  }
  int getNumMacroSource(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->macroSources.size();
  }
  int getNumFibroSource(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->fibroSources.size();
  }

  int get_xDim(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->xDimension;
  }
  int get_yDim(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->yDimension;
  }
  int get_current_step(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->current_step;
  }

  // float* mainSimulation(int rseed, int numMatrixElements, float* internalParameterization){
  //   SimulationObject sim = SimulationObject(internalParameterization);
  //   cout<<"obj created"<<endl;
  //   sim.setupSimulation(rseed);
  //   float* result = sim.sim_obj_main();
  //   return result;
  // }
  float* endSimulation(void* ptr){
    SimulationObject* ref = reinterpret_cast<SimulationObject *>(ptr);
    return ref->endSimulation();
  }

  void setSeed(void* ptr, int newSeed){
      SimulationObject * ref = reinterpret_cast<SimulationObject *>(ptr);
      ref->setSeed(newSeed);
  }
}







// int main(int argc, char const *argv[]) {
//   /* code */
//   SimulationObject sim = SimulationObject();
//   float* result = sim.sim_obj_main();
//   return 0;
// }
