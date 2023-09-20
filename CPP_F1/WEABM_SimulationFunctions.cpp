#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

using namespace std;

float tempSignal[xDim][yDim][zDim];

extern uniform_real_distribution<float> dist01;

void wiggle(int* orientation, int* x, int* y);
void getAhead(int orient, int x, int y, int *xl, int *xm, int *xr, int *yl, int *ym, int *yr);
void adjustOrientation(int* orientation, int leftOrRight);
void move(int orient);
// void evaporate(float &cytokineGrid[xDim][yDim][zDim], float evaporationConstant);
// void gridOutput(int step);
float cytokineComboRule(int ruleRow, int cellIndex, float tnfr, float il1r);
float checkGridOccupancy(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]);
bool checkValidPoint(int x, int y, int z);
// void setupBiopsySatellites();

void SimulationObject::initialize(){
  int i,j,k,index=0;
  float temp;
  current_step=0;
//  endothelium.clear();
  gridIndexes.clear();
  fibroblasts.clear();
  pmnSources.clear();
  macroSources.clear();
  T_Sources.clear();
  fibroSources.clear();
  pmns.clear();
  macrophages.clear();

  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){
        temp=dist01(generator);
        tissueGrid[i][j][k].xLoc=i;
        tissueGrid[i][j][k].yLoc=j;
        tissueGrid[i][j][k].zLoc=k;
        tissueGrid[i][j][k].innervation=1.0;
        tissueGrid[i][j][k].collagen1=0.0;
        tissueGrid[i][j][k].collagen3=0.0;
        tissueGrid[i][j][k].viableMuscle=1.0;
        tissueGrid[i][j][k].vascularization=1.0;
        tissueGrid[i][j][k].life=100;
        tissueGrid[i][j][k].necrosis=0;
        tissueGrid[i][j][k].cellPopulation=0;
        tissueGrid[i][j][k].woundEdge=0;
        tissueGrid[i][j][k].ec_roll=0;
        tissueGrid[i][j][k].ec_stick=0;
        tissueGrid[i][j][k].ec_activation=0;
        tissueGrid[i][j][k].ec_migrate=0;
        tissueGrid[i][j][k].id=i+j*xDim+k*xDim*yDim;
//        cout<<"IDTEST="<<tissueGrid[i][j][k].id<<" "<<i<<" "<<j<<" "<<k<<"\n";
        //z=floor(id/(xDim*yDim));  y=floor((id-z*xDim*yDim)/xDim); x=id-y*xDim-z*xDim*yDim;
        gridIndexes.push_back(tissueGrid[i][j][k].id);
        if(temp<0.05){
          tissueGrid[i][j][k].satellitePool=true;
        } else {
          tissueGrid[i][j][k].satellitePool=false;
        }
        if(temp>0.05 && temp<0.07){
//          cout<<"Placing PMN Source\n";
          pmnSources.push_back(pmnSource(i,j,k));
        }
        if(temp>0.1 && temp<0.12){
          macroSources.push_back(macroSource(i,j,k));
        }
        if(temp>0.15 && temp<0.150001){
          T_Sources.push_back(T_Source(i,j,k));
        }
        if(temp>0.2 && temp<0.22){
          fibroSources.push_back(fibroSource(i,j,k));
        }
        index++;

        cytoGrids.IL8[i][j][k]=0.0;
        cytoGrids.MIP1b[i][j][k]=0.0;
        cytoGrids.TNF[i][j][k]=0.0;
        cytoGrids.PAF[i][j][k]=0.0;
        cytoGrids.IL1[i][j][k]=0.0;
        cytoGrids.PDGF[i][j][k]=0.0;
        cytoGrids.IL13[i][j][k]=0.0;
        cytoGrids.TGFb[i][j][k]=0.0;
        cytoGrids.VEGF[i][j][k]=0.0;
        cytoGrids.GCSF[i][j][k]=0.0;
        cytoGrids.IL10[i][j][k]=0.0;
        cytoGrids.IL4[i][j][k]=0.0;
        cytoGrids.IL17[i][j][k]=0.0;
        cytoGrids.Myf5[i][j][k]=0.0;
        cytoGrids.IL6[i][j][k]=0.0;
        cytoGrids.IL1ra[i][j][k]=0.0;
        cytoGrids.sIL1r[i][j][k]=0.0;
        cytoGrids.sTNFr[i][j][k]=0.0;
        cytoGrids.endotoxin[i][j][k]=0.0;
        cytoGrids.IFNg[i][j][k]=0.0;
        cytoGrids.DAMP[i][j][k]=0.0;
        cytoGrids.cytotox[i][j][k]=0.0;
        cytoGrids.FNE[i][j][k]=0.0;
      }
    }
  }
}

void SimulationObject::injure_biopsy(){
//half-ellipsoid injury: x2/a2 + y2/b2 + z2/c2 = 1; tissue us removed in the
//half-ellipsiod; there is a boundary around the volume of removed tissue
//which then has a partial endothelial injury, and kick-starts the healing process
//by releasing some IL8, PAF, and MIP1b
// cout<<"INJURING INSIDE C++"<<endl;
  int i,j,k,length;
//  cout<<"Test IB0 "<<tissueGrid[6][6][0].vascularization<<"\n";
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=zDim-1;k>10;k--){
//        cout<<"K="<<k<<"\n";
//        cout<<"Test IB0 "<<i<<","<<j<<","<<k<<","<<tissueGrid[6][6][0].vascularization<<"\n";
        tissueGrid[i][j][k].life=0;
        tissueGrid[i][j][k].innervation=0.0;
        tissueGrid[i][j][k].collagen1=0.0;
        tissueGrid[i][j][k].collagen3=0.0;
        tissueGrid[i][j][k].viableMuscle=0.0;
        tissueGrid[i][j][k].vascularization=0.0;
        tissueGrid[i][j][k].necrosis=0;
        tissueGrid[i][j][k].contamination=0;
        tissueGrid[i][j][k].satellitePool=false;
//        cout<<"Test IB0.5 "<<i<<","<<j<<","<<k<<","<<tissueGrid[6][6][0].vascularization<<"\n";
      }
    }
  }
//  cout<<"Test IB1 "<<tissueGrid[6][6][0].vascularization<<"\n";
  k=10;
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      tissueGrid[i][j][k].contamination=1;
    }
  }

  length=fibroSources.size();
  for(j=0;j<length;j++){
    if(fibroSources[j].zLoc>=10){
      fibroSources.erase(fibroSources.begin()+j);
      j--;
      length=fibroSources.size();
    }
  }
  length=T_Sources.size();
  for(j=0;j<length;j++){
    if(T_Sources[j].zLoc>=10){
      T_Sources.erase(T_Sources.begin()+j);
      j--;
      length=T_Sources.size();
    }
  }

  length=macroSources.size();
  for(j=0;j<length;j++){
    if(macroSources[j].zLoc>=10){
      macroSources.erase(macroSources.begin()+j);
      j--;
      length=macroSources.size();
    }
  }

  length=pmnSources.size();
  for(j=0;j<length;j++){
    if(pmnSources[j].zLoc>=10){
      pmnSources.erase(pmnSources.begin()+j);
      j--;
      length=pmnSources.size();
    }
  }

  // length=satelliteSources.size();
  // for(j=0;j<length;j++){
  //   if(satelliteSources[j].zLoc>=10){
  //     satelliteSources.erase(satelliteSources.begin()+j);
  //     j--;
  //     length=satelliteSources.size();
  //   }
  // }
//  cout<<"Test IB2 "<<tissueGrid[6][6][0].vascularization<<"\n";
  return;
}

void SimulationObject::setupBiopsySatellites(){ //setting up satellite sources along the edge of the simulation; satellite cells then move along the y-axis and remain in an x-column
  int i,j;
  for(i=0;i<xDim;i++){
//    for(j=1;j<zDim;j++){
  for(j=1;j<11;j++){
      satelliteSources.push_back(satelliteSource(i,0,j));
    }
  }
}

void SimulationObject::initializeFromCT(){
  ifstream CTfile;
  int i,j,k,index=0;
  vector<int> woundSurface,xs,ys,depths;
  string tempIN;
  int tempDepth;
  float temp,t2;
  int top,numElements,testCounter=0;
  cout<<"Starting CT Read\n";
  current_step=0;
//  endothelium.clear();
  gridIndexes.clear();
  fibroblasts.clear();
  pmnSources.clear();
  macroSources.clear();
  T_Sources.clear();
  fibroSources.clear();
  pmns.clear();
  macrophages.clear();
  cout<<"CT Test 10\n";

  CTfile.open(ctInputFile);
  cout<<"CT Test 11\n";
  while(!CTfile.eof()){
    CTfile>>tempIN;
    cout<<"tempin="<<tempIN<<"\n";
    t2=stof(tempIN);
    cout<<t2<<"\n";
    woundSurface.push_back(int(t2));
    testCounter++;
  }
  cout<<"CT Test 20\n";

  k=0;
  for(i=0; i<xDim;i++){
    for(j=0;j<yDim;j++){
      xs.push_back(i);
      ys.push_back(j);
      tempDepth=zDim-abs(woundSurface[k]);
      cout<<"TD="<<tempDepth<<"\n";
      depths.push_back(tempDepth);
      k=k+1;
    }
  }

  cout<<"CT Test 30\n";
  for(i=0;i<depths.size();i++){
    top=zDim-depths[i];
    for(j=0;j<top;j++){
      temp=dist01(generator);
      tissueGrid[xs[i]][ys[i]][j].xLoc=xs[i];
      tissueGrid[xs[i]][ys[i]][j].yLoc=ys[i];
      tissueGrid[xs[i]][ys[i]][j].zLoc=j;
      tissueGrid[xs[i]][ys[i]][j].innervation=1.0;
      tissueGrid[xs[i]][ys[i]][j].collagen1=0.0;
      tissueGrid[xs[i]][ys[i]][j].collagen3=0.0;
      tissueGrid[xs[i]][ys[i]][j].viableMuscle=1.0;
      tissueGrid[xs[i]][ys[i]][j].vascularization=1.0;
      tissueGrid[xs[i]][ys[i]][j].life=1.0;
      tissueGrid[xs[i]][ys[i]][j].necrosis=0;
      tissueGrid[xs[i]][ys[i]][j].cellPopulation=0;
      tissueGrid[xs[i]][ys[i]][j].woundEdge=0;
      tissueGrid[xs[i]][ys[i]][j].ec_roll=0;
      tissueGrid[xs[i]][ys[i]][j].ec_stick=0;
      tissueGrid[xs[i]][ys[i]][j].ec_activation=0;
      tissueGrid[xs[i]][ys[i]][j].ec_migrate=0;
      tissueGrid[xs[i]][ys[i]][j].id=xs[i]+ys[i]*xDim+j*xDim*yDim;
      gridIndexes.push_back(tissueGrid[xs[i]][ys[i]][j].id);
      if(temp<0.05){
        tissueGrid[xs[i]][ys[i]][j].satellitePool=true;
      } else {
        tissueGrid[xs[i]][ys[i]][j].satellitePool=false;
      }
      if(temp>0.05 && temp<0.07){
//          cout<<"Placing PMN Source\n";
        pmnSources.push_back(pmnSource(xs[i],ys[i],j));
      }
      if(temp>0.1 && temp<0.12){
        macroSources.push_back(macroSource(xs[i],ys[i],j));
      }
      if(temp>0.15 && temp<0.150001){
        T_Sources.push_back(T_Source(xs[i],ys[i],j));
      }
      if(temp>0.2 && temp<0.22){
        fibroSources.push_back(fibroSource(xs[i],ys[i],j));
      }
      index++;
      if(j==top-1){
        tissueGrid[xs[i]][ys[i]][j].contamination=1;
      }

      cytoGrids.IL8[xs[i]][ys[i]][j]=0.0;
      cytoGrids.MIP1b[xs[i]][ys[i]][j]=0.0;
      cytoGrids.TNF[xs[i]][ys[i]][j]=0.0;
      cytoGrids.PAF[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IL1[xs[i]][ys[i]][j]=0.0;
      cytoGrids.PDGF[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IL13[xs[i]][ys[i]][j]=0.0;
      cytoGrids.TGFb[xs[i]][ys[i]][j]=0.0;
      cytoGrids.VEGF[xs[i]][ys[i]][j]=0.0;
      cytoGrids.GCSF[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IL10[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IL4[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IL17[xs[i]][ys[i]][j]=0.0;
      cytoGrids.Myf5[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IL6[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IL1ra[xs[i]][ys[i]][j]=0.0;
      cytoGrids.sIL1r[xs[i]][ys[i]][j]=0.0;
      cytoGrids.sTNFr[xs[i]][ys[i]][j]=0.0;
      cytoGrids.endotoxin[xs[i]][ys[i]][j]=0.0;
      cytoGrids.IFNg[xs[i]][ys[i]][j]=0.0;
      cytoGrids.DAMP[xs[i]][ys[i]][j]=0.0;
      cytoGrids.cytotox[xs[i]][ys[i]][j]=0.0;
      cytoGrids.FNE[xs[i]][ys[i]][j]=0.0;
    }
  }
  CTfile.close();
  cout<<"Wound Successfulyl Read from CT\n";
}

int SimulationObject::checkForInitialNeighboringTissueVoid(int x, int y, int z){ //only to be used in initial surgery/injury
  int i,j,k,flag=0;

  for(i=x-1;i<=x+1;i++){
    for(j=y-1;j<=y+1;j++){
      for(k=z-1;k<=z+1;k++){
        if(tissueGrid[i][j][k].viableMuscle==0.0){
          flag=1;
          return flag;
        }
      }
    }
  }
  return flag;
}

void SimulationObject::diffuse(float (&cytokineGrid)[xDim][yDim][zDim], float diffusionConstant){
  int i,j,k,ii,jj,kk,iii,jjj,kkk;
  //float tempSignal[xDim][yDim][zDim];
  float nFactor,viability;
  nFactor=1.0/26.0*diffusionConstant;
//  cout<<"Test 510\n";
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){
        tempSignal[i][j][k]=0.0;
      }
    }
  }
//  cout<<"Test 520\n";
  //body diffusion: The temp signal grid cell that is not on the boundary of the simulation
  //gets a some cytokine from each of its neighbors
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=1;k<zDim;k++){

        // cout<<"Test 530 "<<i<<","<<j<<","<<k<<"\n";
        for(ii=i-1;ii<=i+1;ii++){
          for(jj=j-1;jj<=j+1;jj++){
            for(kk=k-1;kk<=k+1;kk++){
              // cout<<"Test 540 "<<ii<<","<<jj<<","<<kk<<"\n";
              // cout<<"Test 545 "<<i<<","<<j<<","<<k<<"\n";
              if(ii==i&&jj==j&&kk==k){
                continue;
              } else {
                iii=ii;
                jjj=jj;
                kkk=kk;
                if(ii<0){
                  iii=xDim-1;
                }
                if(ii>=xDim){
                  iii=0;
                }
                if(jj<0){
                  jjj=yDim-1;
                }
                if(jj>=yDim){
                  jjj=0;
                }
                if(kk<0){
                  continue;
                }
                if(kk>=zDim){
                  continue;
                }
                viability=tissueGrid[iii][jjj][kk].collagen1+tissueGrid[iii][jjj][kk].collagen3+tissueGrid[iii][jjj][kk].viableMuscle+tissueGrid[iii][jjj][kk].necrosis;
                if(viability<=0.0){
                  continue;
                }
                tempSignal[i][j][k]+=nFactor*cytokineGrid[iii][jjj][kkk];
    }}}}}}}
  //top/bottom diffusion: Cells on the edge will diffuse cytokines away out of the simulation,
  //but they will not receive any cytokine from outside the boundary of the sim.
 //main cytokine grid is updated with the temporary signal
// cout<<"Test 580\n";
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){
        cytokineGrid[i][j][k]-=cytokineGrid[i][j][k]*diffusionConstant;
        cytokineGrid[i][j][k]+=tempSignal[i][j][k];
      }
    }
  }
//  cout<<"Test 590\n";
}

void SimulationObject::evaporate(float (&cytokineGrid)[xDim][yDim][zDim], float evaporationConstant){
  int i,j,k;
  float viability;
//  cout<<"Test 610\n";
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){
        cytokineGrid[i][j][k]*=evaporationConstant;
        if(cytokineGrid[i][j][k]<minimumCytokineValue){
          cytokineGrid[i][j][k]=0.0;
        }
        viability=tissueGrid[i][j][k].collagen1+tissueGrid[i][j][k].collagen3+tissueGrid[i][j][k].viableMuscle+tissueGrid[i][j][k].necrosis;
        if(viability<=0){
          cytokineGrid[i][j][k]=0.0;
        }
      }
    }
  }
//  cout<<"Test 690\n";
}

void SimulationObject::outputSetup(){
  if(outputWrite==1){
	   outputFile.open(filename1);
  }
  if(tissueWrite==1){
    tissueFile.open(filename2);
  }
}

void SimulationObject::outputFinalize(){
  if(outputWrite==1){
	   outputFile.close();
  }
  if(tissueWrite==1){
    tissueFile.close();
  }
}

void SimulationObject::outputWriter(float cytokine[xDim][yDim][zDim], int step, int fixedDim){
    int i,k;
    for(i=0;i<xDim;i++){
      for(k=0;k<zDim;k++){
        outputFile<<step<<","<<i<<","<<fixedDim<<","<<k<<","
          <<tissueGrid[i][fixedDim][k].viableMuscle<<"\n";
    }
  }
}

void SimulationObject::binaryTissueOutput(){
//outputs total tissue volume in each voxel - to be used to generate images with 3dSlicer
  int i,j,k;
  float tissue;
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){
        tissue=tissueGrid[i][j][k].collagen1+tissueGrid[i][j][k].collagen3+tissueGrid[i][j][k].viableMuscle;
        tissueFile<<i<<","<<j<<","<<k<<","<<tissue<<"\n";
      }
    }
  }
}

void SimulationObject::aggregateOutput(){
  int i,j,k;
  float totalM,totalC1,totalC3;
  totalM=0;
  totalC1=0;
  totalC3=0;
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){

      }
    }
  }
  outputFile<<totalM<<","<<totalC1<<","<<totalC3<<","<<fibroblasts.size()<<"\n";
}

void SimulationObject::diffusions(){ //The constant is how much total gets transferred out of the voxel
  diffuse(cytoGrids.PAF,0.6);
  diffuse(cytoGrids.IL8,0.6);
  diffuse(cytoGrids.DAMP,0.45);
  diffuse(cytoGrids.cytotox,0.75);
  diffuse(cytoGrids.TNF,0.6);
  diffuse(cytoGrids.IL1,0.6);
  diffuse(cytoGrids.IL1ra,0.6);
  diffuse(cytoGrids.sIL1r,0.6);
  diffuse(cytoGrids.IL10,0.6);
  diffuse(cytoGrids.IFNg,0.6);
  diffuse(cytoGrids.TGFb,0.6);
  diffuse(cytoGrids.MCP1,0.6);
  diffuse(cytoGrids.IL4,0.4);
  diffuse(cytoGrids.MIP1b,0.6);
  diffuse(cytoGrids.HGF,0.6);
  // diffuse(cytoGrids.PDGF,0.6);
  // diffuse(cytoGrids.IL13,0.6);
  // diffuse(cytoGrids.VEGF,0.6);
  // diffuse(cytoGrids.GCSF,0.6);
  diffuse(cytoGrids.IL17,0.6);
  // diffuse(cytoGrids.Myf5,0.6);
  diffuse(cytoGrids.IL6,0.6);
  // diffuse(cytoGrids.sTNFr,0.6);
  // diffuse(cytoGrids.endotoxin,0.6);
  // diffuse(cytoGrids.IL12,0.6);
  diffuse(cytoGrids.MRF4,0.6);
  // diffuse(cytoGrids.IGF,0.6);
  diffuse(cytoGrids.FNE,0.5);


}

void SimulationObject::evaporations(){ //The constant is how much remains after evaporation
  evaporate(cytoGrids.PAF,0.7);
  //evaporate(cytoGrids.IL8,0.8);
  evaporate(cytoGrids.IL8,0.8);
  evaporate(cytoGrids.DAMP,0.6);
  evaporate(cytoGrids.cytotox,0.6);
  evaporate(cytoGrids.TNF,0.8);
  evaporate(cytoGrids.IL1,0.1);
  evaporate(cytoGrids.IL1ra,0.1);
  evaporate(cytoGrids.sIL1r,0.1);
  evaporate(cytoGrids.IL10,0.9);
  evaporate(cytoGrids.IFNg,0.8);
  evaporate(cytoGrids.TGFb,0.5);
  evaporate(cytoGrids.MCP1,0.5);
  evaporate(cytoGrids.IL4,0.9);
  evaporate(cytoGrids.MIP1b,0.1);
  evaporate(cytoGrids.HGF,0.1);
  evaporate(cytoGrids.PDGF,0.1);
  evaporate(cytoGrids.IL13,0.8);
  evaporate(cytoGrids.VEGF,0.1);
  evaporate(cytoGrids.GCSF,0.8);
  evaporate(cytoGrids.IL17,0.8);
  evaporate(cytoGrids.Myf5,0.1);
  evaporate(cytoGrids.IL6,0.9);
  evaporate(cytoGrids.sTNFr,0.1);
  evaporate(cytoGrids.endotoxin,0.1);
  evaporate(cytoGrids.IL12,0.8);
  evaporate(cytoGrids.MRF4,0.1);
  evaporate(cytoGrids.IGF,0.1);
  evaporate(cytoGrids.FNE,0.8);
}

void SimulationObject::aggregateData(int step){
  int i,j,k;

  TotalIL8=0.0;
  TotalMIP1b=0;
  TotalTNF=0.0;
  TotalPAF=0.0;
  TotalIL1=0.0;
  TotalPDGF=0;
  TotalIL13=0;
  TotalTGFb=0.0;
  TotalVEGF=0;
  TotalGCSF=0;
  TotalIL10=0;
  TotalIL4=0.0;
  TotalIL17=0.0;
  TotalMyf5=0.0;
  TotalIL6=0.0;
  TotalIL1ra=0.0;
  TotalsIL1r=0.0;
  TotalsTNFr=0.0;
  TotalIFNg=0.0;
  TotalMCP1=0.0;
  TotalHGF=0.0;
  TotalIL12=0.0;
  TotalMRF4=0.0;
  TotalIGF=0.0;
  TotalDAMP=0.0;

  TotalIL1r=0.0;
  TotalTNFr=0.0;

  TotalNecrosis=0.0;
  TotalC3=0.0;
  TotalC1=0;
  TotalViabMus=0.0;



  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){

        TotalIL8+=cytoGrids.IL8[i][j][k];
        TotalMIP1b+=cytoGrids.MIP1b[i][j][k];
        TotalTNF+=cytoGrids.TNF[i][j][k];
        TotalPAF+=cytoGrids.PAF[i][j][k];
        TotalIL1+=cytoGrids.IL1[i][j][k];
        TotalPDGF+=cytoGrids.PDGF[i][j][k];
        TotalIL13+=cytoGrids.IL13[i][j][k];
        TotalTGFb+=cytoGrids.TGFb[i][j][k];
        TotalVEGF+=cytoGrids.VEGF[i][j][k];
        TotalGCSF+=cytoGrids.GCSF[i][j][k];
        TotalIL10+=cytoGrids.IL10[i][j][k];
        TotalIL4+=cytoGrids.IL4[i][j][k];
        TotalMyf5+=cytoGrids.Myf5[i][j][k];
        TotalIL6+=cytoGrids.IL6[i][j][k];
        TotalIL1ra+=cytoGrids.IL1ra[i][j][k];
        TotalsIL1r+=cytoGrids.sIL1r[i][j][k];
        TotalsTNFr+=cytoGrids.sTNFr[i][j][k];
        TotalIFNg+=cytoGrids.IFNg[i][j][k];
        TotalMCP1+=cytoGrids.MCP1[i][j][k];
        TotalHGF+=cytoGrids.HGF[i][j][k];
        TotalIL12+=cytoGrids.IL12[i][j][k];
        TotalMRF4+=cytoGrids.MRF4[i][j][k];
        TotalIGF+=cytoGrids.IGF[i][j][k];
        TotalDAMP+=cytoGrids.DAMP[i][j][k];
        TotalNecrosis+=tissueGrid[i][j][k].necrosis;
        TotalC3+=tissueGrid[i][j][k].collagen3;
        // if(tissueGrid[i][j][k].collagen3>0){
        // cout<<i<<" "<<j<<" "<<k<<" "<<tissueGrid[i][j][k].collagen3<<"\n";}
        TotalC1+=tissueGrid[i][j][k].collagen1;
        TotalViabMus=TotalViabMus+tissueGrid[i][j][k].viableMuscle;
        TotalMCP1+=cytoGrids.MCP1[i][j][k];
      }
    }
  }
//  cout<<"Total DAMP="<<TotalDAMP<<"\n";
  for(i=0;i<macrophages.size();i++){
    TotalIL1r+=macrophages[i].IL1r;
    TotalTNFr+=macrophages[i].TNFr;
  }

  allSignals[0][step]=pmns.size();
  allSignals[1][step]=macrophages.size();
  allSignals[2][step]=fibroblasts.size();
  allSignals[3][step]=TotalIL8;
  allSignals[4][step]=TotalMIP1b;
  allSignals[5][step]=TotalTNF;
  allSignals[6][step]=TotalPAF;
  allSignals[7][step]=TotalIL1;
  allSignals[8][step]=TotalPDGF;
  allSignals[9][step]=TotalIL13;
  allSignals[10][step]=TotalTGFb;
  allSignals[11][step]=TotalVEGF;
  allSignals[12][step]=TotalGCSF;
  allSignals[13][step]=TotalIL10;
  allSignals[14][step]=TotalIL4;
  allSignals[15][step]=TotalIL17;
  allSignals[16][step]=TotalMyf5;
  allSignals[17][step]=TotalIL6;
  allSignals[18][step]=TotalIL1ra;
  allSignals[19][step]=TotalsIL1r;
  allSignals[20][step]=TotalsTNFr;
  allSignals[21][step]=TotalEndotoxin;
  allSignals[22][step]=TotalIFNg;
  allSignals[23][step]=TotalDAMP;
  allSignals[24][step]=TotalNecrosis;
  allSignals[25][step]=TotalC1;
  allSignals[26][step]=TotalC3;
  allSignals[27][step]=TotalViabMus;
  allSignals[28][step]=TotalMCP1;
  allSignals[29][step]=TotalHGF;
  allSignals[30][step]=TotalIL12;
  allSignals[31][step]=TotalMRF4;
  allSignals[32][step]=TotalIGF;
  allSignals[33][step]=Th0s.size();
  allSignals[34][step]=Th1s.size();
  allSignals[35][step]=Th2s.size();
  allSignals[36][step]=Th17s.size();
  allSignals[37][step]=TotalIL1r;
  allSignals[38][step]=TotalTNFr;
  allSignals[39][step]=myoblasts.size();
  allSignals[40][step]=satellites.size();
  allSignals[41][step]=macroSources.size();
  allSignals[42][step]=fibroSources.size();

//  cout<<"TEST "<<TotalIL8<<" "<<TotalIL8<<"\n";
  //  outputFile<<step<<","<<TotalDAMP<<","<<TotalNecrosis<<","<<pmns.size()<<","<<macrophages.size()<<","<<TotalTGFb<<","<<fibroblasts.size()<<","<<TotalViabMus<<","<<TotalC3<<"\n";

  // outputFile<<step<<","<<pmns.size()<<","<<macrophages.size()<<","<<fibroblasts.size()<<","
  // <<TotalIL8<<","<<TotalMIP1b<<","<<TotalTNF<<","<<TotalPAF<<","<<TotalIL1<<","
  // <<TotalPDGF<<","<<TotalIL13<<","<<TotalTGFb<<","<<TotalVEGF<<","<<TotalGCSF<<
  // ","<<TotalIL10<<","<<TotalIL4<<","<<TotalIL17<<","<<TotalMyf5<<","<<TotalIL6
  // <<","<<TotalIL1ra<<","<<TotalsIL1r<<","<<TotalsTNFr<<","<<TotalEndotoxin<<","
  // <<TotalIFNg<<","<<TotalDAMP<<","<<TotalNecrosis<<","<<TotalC1<<","<<TotalC3<<
  // ","<<TotalViabMus<<","<<TotalMCP1<<","<<TotalHGF<<","<<TotalIL12<<","<<TotalMRF4<<","<<TotalIGF<<
  // ","<<Th0s.size()<<","<<Th1s.size()<<","<<Th2s.size()<<","<<Th17s.size()<<","<<TotalIL1r<<","<<TotalTNFr<<","<<myoblasts.size()<<","<<satellites.size()<<"\n";

}

void SimulationObject::gridOutput(int step){
  int i,j,k;
  float temp1,temp2,temp3,temp4;

  for(i=0;i<xDim;i++){
//    for(j=0;j<yDim;j++){
    for(j=4;j<15;j+=10){
      for(k=0;k<zDim;k++){
        temp1=tissueGrid[i][j][k].viableMuscle;
        temp2=tissueGrid[i][j][k].collagen3;
        temp4=tissueGrid[i][j][k].collagen1;
        temp3=tissueGrid[i][j][k].life;
//        cout<<step<<","<<i<<","<<j<<","<<k<<","<<temp1<<","<<temp2<<"\n";
        tissueFile<<step<<","<<int(i)<<","<<int(j)<<","<<int(k)<<","<<float(temp1)<<","<<float(temp2)<<","<<float(temp3)<<","<<float(temp4)<<"\n";
      }
    }
  }
}

float checkGridOccupancy(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]){
  float answer,cont;
  answer=tissueGrid[x][y][z].viableMuscle+tissueGrid[x][y][z].collagen1+
    tissueGrid[x][y][z].collagen3+tissueGrid[x][y][z].necrosis;
  cont=tissueGrid[x][y][z].contamination;
//  if(cont>0){answer=100;}
  return answer;
}

float checkGridOccupancyC3(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]){
  float answer,cont;
  answer=tissueGrid[x][y][z].collagen3;
  return answer;
}

bool checkValidPoint(int x, int y, int z){
  bool answer=false;
  if(x>=0&&x<xDim){
    if(y>=0&&y<yDim){
      if(z>=0&&z<zDim){
        answer=true;
      }
    }
  }
  return answer;
}

void SimulationObject::getTopGrid(){
  int i,j,k;
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){
        if(checkGridOccupancy(i,j,k, tissueGrid)==0){
          topGrid[i][j]=k -1;
          break;
        }
      }
    }
  }
}

float cytokineProductionRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr){
  //Note that the result is contrained to be positive
  int i;
  float tempSum=0;
  float currentCyto[30];

  currentCyto[0]=cytoGrids.IL8[x][y][z];
  currentCyto[1]=cytoGrids.MIP1b[x][y][z];
  currentCyto[2]=cytoGrids.TNF[x][y][z];
  currentCyto[3]=cytoGrids.PAF[x][y][z];
  currentCyto[4]=cytoGrids.IL1[x][y][z];
  currentCyto[5]=cytoGrids.PDGF[x][y][z];
  currentCyto[6]=cytoGrids.IL13[x][y][z];
  currentCyto[7]=cytoGrids.TGFb[x][y][z];
  currentCyto[8]=cytoGrids.VEGF[x][y][z];
  currentCyto[9]=cytoGrids.GCSF[x][y][z];
  currentCyto[10]=cytoGrids.IL10[x][y][z];
  currentCyto[11]=cytoGrids.IL4[x][y][z];
  currentCyto[12]=cytoGrids.IL17[x][y][z];
  currentCyto[13]=cytoGrids.Myf5[x][y][z];
  currentCyto[14]=cytoGrids.IL6[x][y][z];
  currentCyto[15]=cytoGrids.IL1ra[x][y][z];
  currentCyto[16]=cytoGrids.sIL1r[x][y][z];
  currentCyto[17]=cytoGrids.sTNFr[x][y][z];
  currentCyto[18]=cytoGrids.endotoxin[x][y][z];
  currentCyto[19]=cytoGrids.IFNg[x][y][z];
  currentCyto[20]=cytoGrids.DAMP[x][y][z];
  currentCyto[21]=cytoGrids.cytotox[x][y][z];
  currentCyto[22]=cytoGrids.MCP1[x][y][z];
  currentCyto[23]=cytoGrids.HGF[x][y][z];
  currentCyto[24]=cytoGrids.IL12[x][y][z];
  currentCyto[25]=cytoGrids.MRF4[x][y][z];
  currentCyto[26]=cytoGrids.IGF[x][y][z];
  currentCyto[27]=il1r;
  currentCyto[28]=tnfr;
  currentCyto[29]=cytoGrids.FNE[x][y][z];

  // cout<<"############MATRIX DEBUG#########################################\n";
  // for(i=0;i<numRuleParams-1;i++){
  //   if(RM[ruleRow][i]!=0&&currentCyto[i]!=0){
  //     cout<<"Col="<<i<<" Coef="<<RM[ruleRow][i]<<" Cyto="<<currentCyto[i]<<"\n";
  //   }
  // }
  for(i=0;i<numRuleParams-1;i++){
    tempSum+=currentCyto[i]*RM[ruleRow][i];
  }
  tempSum+=float(RM[ruleRow][30]);
  if(tempSum<0){tempSum=0;}
  return tempSum;
}

float cytokineComboRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr){
  int i;
  float tempSum=0;
  float currentCyto[31];

  currentCyto[0]=cytoGrids.IL8[x][y][z];
  currentCyto[1]=cytoGrids.MIP1b[x][y][z];
  currentCyto[2]=cytoGrids.TNF[x][y][z];
  currentCyto[3]=cytoGrids.PAF[x][y][z];
  currentCyto[4]=cytoGrids.IL1[x][y][z];
  currentCyto[5]=cytoGrids.PDGF[x][y][z];
  currentCyto[6]=cytoGrids.IL13[x][y][z];
  currentCyto[7]=cytoGrids.TGFb[x][y][z];
  currentCyto[8]=cytoGrids.VEGF[x][y][z];
  currentCyto[9]=cytoGrids.GCSF[x][y][z];
  currentCyto[10]=cytoGrids.IL10[x][y][z];
  currentCyto[11]=cytoGrids.IL4[x][y][z];
  currentCyto[12]=cytoGrids.IL17[x][y][z];
  currentCyto[13]=cytoGrids.Myf5[x][y][z];
  currentCyto[14]=cytoGrids.IL6[x][y][z];
  currentCyto[15]=cytoGrids.IL1ra[x][y][z];
  currentCyto[16]=cytoGrids.sIL1r[x][y][z];
  currentCyto[17]=cytoGrids.sTNFr[x][y][z];
  currentCyto[18]=cytoGrids.endotoxin[x][y][z];
  currentCyto[19]=cytoGrids.IFNg[x][y][z];
  currentCyto[20]=cytoGrids.DAMP[x][y][z];
  currentCyto[21]=cytoGrids.cytotox[x][y][z];
  currentCyto[22]=cytoGrids.MCP1[x][y][z];
  currentCyto[23]=cytoGrids.HGF[x][y][z];
  currentCyto[24]=cytoGrids.IL12[x][y][z];
  currentCyto[25]=cytoGrids.MRF4[x][y][z];
  currentCyto[26]=cytoGrids.IGF[x][y][z];
  currentCyto[27]=il1r;
  currentCyto[28]=tnfr;
  currentCyto[29]=cytoGrids.FNE[x][y][z];

  // cout<<"Cytos RULE="<<MCP1[x][y][z]<<" "<<PAF[x][y][z]<<" "<<IFNg[x][y][z]<<" "<<IL10[x][y][z]<<"\n";
  //cout<<"############################################################################################\n";
  for(i=0;i<numRuleParams-1;i++){
    // if(RM[ruleRow][i]!=0){
    //   cout<<"Row="<<i<<"\n";
    // }
    tempSum+=currentCyto[i]*RM[ruleRow][i];
  }
//  cout<<"Const="<<RM[ruleRow][16]<<"\n";
  tempSum+=RM[ruleRow][30];
  return tempSum;
}

void SimulationObject::getRuleMatrix(float internalParam[], int numMatEls){
  int i,j,k=0;
  float input;
  for(i=0;i<numRules;i++){
    for(j=0;j<numRuleParams;j++){
      input=internalParam[k];
      if(abs(input)<0.0001){
        input=0.0;
      }
      RM[i][j]=input;
      k++;
    }
  }
  cout<<"NumRules="<<numRules<<" "<<"numCols="<<numRuleParams<<"\n";

  for(i=0;i<numRules;i++){
    for(j=0;j<numRuleParams;j++){
    cout<<RM[i][j]<<" ";
  }
  cout<<"\n";
  }
}

void matrixDebugTester(float test1, float test2, bool exact){
//Tester when implementing new rule rows to the rule matrix; not that the order
//in which things are added/subtracted matters in floating point arithmetic,
//which is the reason for the 'exact' flag.
  if(exact==true){
    if(test1!=test2){
      cout<<"Matrix Rule="<<test2<<"\n";
      cout<<"Coded Rule="<<test1<<"\n";
      exit(0);
    }
  } else {
    if(abs(test1-test2)>0.001){
      cout<<"Matrix Rule="<<test2<<"\n";
      cout<<"Coded Rule="<<test1<<"\n";
      exit(0);
    }
  }
}

float* SimulationObject::getSurfaceCytokines(){
  int i,j, height;

  for (i=0; i<28; i++){
    surfaceCytokines[i] = 0;
  }
  for (i=0; i<xDim;i++){
    for(j=0; j<yDim; j++){
      height = topGrid[i][j];
      surfaceCytokines[0] += cytoGrids.IL8[i][j][height];
      surfaceCytokines[1] += cytoGrids.MIP1b[i][j][height];
      surfaceCytokines[2] += cytoGrids.TNF[i][j][height];
      surfaceCytokines[3] += cytoGrids.HGF[i][j][height];
      surfaceCytokines[4] += cytoGrids.PAF[i][j][height];
      surfaceCytokines[5] += cytoGrids.IL1[i][j][height];
      surfaceCytokines[6] += cytoGrids.PDGF[i][j][height];
      surfaceCytokines[7] += cytoGrids.IL13[i][j][height];
      surfaceCytokines[8] += cytoGrids.TGFb[i][j][height];
      surfaceCytokines[9] += cytoGrids.VEGF[i][j][height];
      surfaceCytokines[10] += cytoGrids.GCSF[i][j][height];
      surfaceCytokines[11] += cytoGrids.IL10[i][j][height];
      surfaceCytokines[12] += cytoGrids.IL4[i][j][height];
      surfaceCytokines[13] += cytoGrids.IL17[i][j][height];
      surfaceCytokines[14] += cytoGrids.Myf5[i][j][height];
      surfaceCytokines[15] += cytoGrids.IL6[i][j][height];
      surfaceCytokines[16] += cytoGrids.IL1ra[i][j][height];
      surfaceCytokines[17] += cytoGrids.sIL1r[i][j][height];
      surfaceCytokines[18] += cytoGrids.sTNFr[i][j][height];
      surfaceCytokines[19] += cytoGrids.endotoxin[i][j][height];
      surfaceCytokines[20] += cytoGrids.IFNg[i][j][height];
      surfaceCytokines[21] += cytoGrids.IL12[i][j][height];
      surfaceCytokines[22] += cytoGrids.MRF4[i][j][height];
      surfaceCytokines[23] += cytoGrids.cytotox[i][j][height];
      surfaceCytokines[24] += cytoGrids.IGF[i][j][height];
      surfaceCytokines[25] += cytoGrids.DAMP[i][j][height];
      surfaceCytokines[26] += cytoGrids.MCP1[i][j][height];
      surfaceCytokines[27] += cytoGrids.FNE[i][j][height];
    }
  }
  return surfaceCytokines;
}

float* SimulationObject::getTotalCytokines(){
  int i,j,k;
  float sum=0.0;
  for (i=0; i<28; i++){
    totalCytokines[i] = 0;
  }

  for (i=0; i<xDim;i++){
    for(j=0; j<yDim; j++){
      for(k=0; k<zDim; k++){
        totalCytokines[0] += cytoGrids.IL8[i][j][k];
        // if(cytoGrids.IL8[i][j][k]!=0.0){
        //   cout<<"IJK="<<i<<" "<<j<<" "<<k<<" "<<cytoGrids.IL8[i][j][k]<<"\n";
        // }
//        sum+=cytoGrids.IL8[i][j][k];
        totalCytokines[1] += cytoGrids.MIP1b[i][j][k];
        totalCytokines[2] += cytoGrids.TNF[i][j][k];
        totalCytokines[3] += cytoGrids.HGF[i][j][k];
        totalCytokines[4] += cytoGrids.PAF[i][j][k];
        totalCytokines[5] += cytoGrids.IL1[i][j][k];
        totalCytokines[6] += cytoGrids.PDGF[i][j][k];
        totalCytokines[7] += cytoGrids.IL13[i][j][k];
        totalCytokines[8] += cytoGrids.TGFb[i][j][k];
        totalCytokines[9] += cytoGrids.VEGF[i][j][k];
        totalCytokines[10] += cytoGrids.GCSF[i][j][k];
        totalCytokines[11] += cytoGrids.IL10[i][j][k];
        totalCytokines[12] += cytoGrids.IL4[i][j][k];
        totalCytokines[13] += cytoGrids.IL17[i][j][k];
        totalCytokines[14] += cytoGrids.Myf5[i][j][k];
        totalCytokines[15] += cytoGrids.IL6[i][j][k];
        totalCytokines[16] += cytoGrids.IL1ra[i][j][k];
        totalCytokines[17] += cytoGrids.sIL1r[i][j][k];
        totalCytokines[18] += cytoGrids.sTNFr[i][j][k];
        totalCytokines[19] += cytoGrids.endotoxin[i][j][k];
        totalCytokines[20] += cytoGrids.IFNg[i][j][k];
        totalCytokines[21] += cytoGrids.IL12[i][j][k];
        totalCytokines[22] += cytoGrids.MRF4[i][j][k];
        totalCytokines[23] += cytoGrids.cytotox[i][j][k];
        totalCytokines[24] += cytoGrids.IGF[i][j][k];
        totalCytokines[25] += cytoGrids.DAMP[i][j][k];
        totalCytokines[26] += cytoGrids.MCP1[i][j][k];
        totalCytokines[27] += cytoGrids.FNE[i][j][k];
      }
    }
  }
//  cout<<"SUM="<<sum<<"\n";
  return totalCytokines;
}

void SimulationObject::updateOutputGrids(){
  int index;
  // cout<<"ENTERING UPDATE OUTPUT"<<endl;
  for (int i = 0; i < xDim; i++) {
    for (int j = 0; j < yDim; j++) {
      for (int k = 0; k < zDim; k++) {
        gridCell currentCell = tissueGrid[i][j][k];
        // if(tissueGrid[i][j][k].collagen3>0){
        // cout<<i<<" "<<j<<" "<<k<<" "<<tissueGrid[i][j][k].collagen3<<"\n";}
        index = i + xDim*j + xDim*yDim*k;

        // cout<<currentCell.collagen1<<", "
        // <<currentCell.collagen3<<", "
        // <<currentCell.viableMuscle<<", "
        // <<currentCell.vascularization<<", "
        // <<currentCell.life<<", "
        // <<currentCell.necrosis<<", "
        // <<currentCell.contamination<<endl;

        innervation_out[index] = currentCell.innervation;
        collagen1_out[index] = currentCell.collagen1;
        collagen3_out[index] = currentCell.collagen3;
        // if(collagen3_out[index]>0){cout<<"TEST="<<collagen3_out[index]<<"\n";}
        viableMuscle_out[index] = currentCell.viableMuscle;
        vascularization_out[index] = currentCell.vascularization;
        life_out[index] = currentCell.life;
        necrosis_out[index] = currentCell.necrosis;
        contamination_out[index] = currentCell.contamination;
      }
    }
  }
}

float* SimulationObject::getInnervation(){
  return innervation_out;
}
float* SimulationObject::getCollagen1(){
  return collagen1_out;
}
float* SimulationObject::getCollagen3(){
  return collagen3_out;
}
float* SimulationObject::getViableMuscle(){
  return viableMuscle_out;
}
float* SimulationObject::getVascularization(){
  return vascularization_out;
}
float* SimulationObject::getLife(){
  return life_out;
}
float* SimulationObject::getNecrosis(){
  return necrosis_out;
}
float* SimulationObject::getContamination(){
  return contamination_out;
}
float* SimulationObject::get_IL8(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL8[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_MIP1b(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.MIP1b[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_TNF(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.TNF[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_HGF(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.HGF[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_PAF(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.PAF[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL1(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL1[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_PDGF(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.PDGF[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL13(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL13[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_TGFb(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.TGFb[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_VEGF(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.VEGF[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_GCSF(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.GCSF[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL10(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL10[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL4(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL4[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL17(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL17[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_Myf5(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.Myf5[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL6(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL6[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL1ra(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL1ra[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_sIL1r(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.sIL1r[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_sTNFr(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.sTNFr[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_endotoxin(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.endotoxin[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IFNg(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IFNg[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IL12(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IL12[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_MRF4(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.MRF4[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_cytotox(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.cytotox[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_IGF(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.IGF[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_DAMP(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.DAMP[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_MCP1(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.MCP1[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}
float* SimulationObject::get_FNE(){
  for(int i=0;i<xDim;i++){
    for(int j=0;j<yDim;j++){
      for(int k=0;k<zDim;k++){
        float currentCyto = cytoGrids.FNE[i][j][k];
        int index = i + xDim*j + xDim*yDim*k;
        temp_out[index] = currentCyto;
      }
    }
  }
  return temp_out;
}

void SimulationObject::apply_IL8(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL8[i][j][height] += addedAmount;
      if (cytoGrids.IL8[i][j][height]<0){cytoGrids.IL8[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_MIP1b(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.MIP1b[i][j][height] += addedAmount;
      if (cytoGrids.MIP1b[i][j][height]<0){cytoGrids.MIP1b[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_TNF(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.TNF[i][j][height] += addedAmount;
      if (cytoGrids.TNF[i][j][height]<0){cytoGrids.TNF[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_HGF(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.HGF[i][j][height] += addedAmount;
      if (cytoGrids.HGF[i][j][height]<0){cytoGrids.HGF[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_PAF(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.PAF[i][j][height] += addedAmount;
      if (cytoGrids.PAF[i][j][height]<0){cytoGrids.PAF[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL1(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL1[i][j][height] += addedAmount;
      if (cytoGrids.IL1[i][j][height]<0){cytoGrids.IL1[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_PDGF(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.PDGF[i][j][height] += addedAmount;
      if (cytoGrids.PDGF[i][j][height]<0){cytoGrids.PDGF[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL13(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL13[i][j][height] += addedAmount;
      if (cytoGrids.IL13[i][j][height]<0){cytoGrids.IL13[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_TGFb(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.TGFb[i][j][height] += addedAmount;
      if (cytoGrids.TGFb[i][j][height]<0){cytoGrids.TGFb[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_VEGF(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.VEGF[i][j][height] += addedAmount;
      if (cytoGrids.VEGF[i][j][height] <0){cytoGrids.VEGF[i][j][height]=0;}
    }
  }
}
void SimulationObject::apply_GCSF(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.GCSF[i][j][height] += addedAmount;
      if (cytoGrids.GCSF[i][j][height]<0){cytoGrids.GCSF[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL10(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL10[i][j][height] += addedAmount;
      if (cytoGrids.IL10[i][j][height]<0){cytoGrids.IL10[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL4(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL4[i][j][height] += addedAmount;
      if (cytoGrids.IL4[i][j][height]<0){cytoGrids.IL4[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL17(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL17[i][j][height] += addedAmount;
      if (cytoGrids.IL17[i][j][height]<0){cytoGrids.IL17[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_Myf5(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.Myf5[i][j][height] += addedAmount;
      if (cytoGrids.Myf5[i][j][height]<0){cytoGrids.Myf5[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL6(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL6[i][j][height] += addedAmount;
      if (cytoGrids.IL6[i][j][height]<0){cytoGrids.IL6[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL1ra(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL1ra[i][j][height] += addedAmount;
      if (cytoGrids.IL1ra[i][j][height]<0){cytoGrids.IL1ra[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_sIL1r(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.sIL1r[i][j][height] += addedAmount;
      if (cytoGrids.sIL1r[i][j][height]<0){cytoGrids.sIL1r[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_sTNFr(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.sTNFr[i][j][height] += addedAmount;
      if (cytoGrids.sTNFr[i][j][height]<0){cytoGrids.sTNFr[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_endotoxin(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.endotoxin[i][j][height] += addedAmount;
      if (cytoGrids.endotoxin[i][j][height]<0){cytoGrids.endotoxin[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IFNg(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IFNg[i][j][height] += addedAmount;
      if (cytoGrids.IFNg[i][j][height]<0){cytoGrids.IFNg[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IL12(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IL12[i][j][height] += addedAmount;
      if (cytoGrids.IL12[i][j][height]<0){cytoGrids.IL12[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_MRF4(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.MRF4[i][j][height] += addedAmount;
      if (cytoGrids.MRF4[i][j][height]<0){cytoGrids.MRF4[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_cytotox(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.cytotox[i][j][height] += addedAmount;
      if (cytoGrids.cytotox[i][j][height]<0){cytoGrids.cytotox[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_IGF(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.IGF[i][j][height] += addedAmount;
      if (cytoGrids.IGF[i][j][height]<0){cytoGrids.IGF[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_DAMP(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.DAMP[i][j][height] += addedAmount;
      if (cytoGrids.DAMP[i][j][height]<0){cytoGrids.DAMP[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_MCP1(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.MCP1[i][j][height] += addedAmount;
      if (cytoGrids.MCP1[i][j][height]<0){cytoGrids.MCP1[i][j][height]=0;};
    }
  }
}
void SimulationObject::apply_FNE(float addedAmount){
  for (int i=0; i<xDim;i++){
    for(int j=0; j<yDim; j++){
      int height = topGrid[i][j];
      cytoGrids.FNE[i][j][height] += addedAmount;
      if (cytoGrids.FNE[i][j][height]<0){cytoGrids.FNE[i][j][height]=0;};
    }
  }
}
