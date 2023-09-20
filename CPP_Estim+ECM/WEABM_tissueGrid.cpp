#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

using namespace std;

void injuryPropagation(int x, int y, int z, float delta,
gridCell (&tissueGrid)[xDim][yDim][zDim]);
void tissueInjuryPropagation(int x, int y, int z, float delta);
void sourceGeneration(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim], vector<pmnSource> &pmnSources, vector<macroSource> &macroSources, vector<T_Source> &T_Sources, vector<fibroSource> &fibroSources,vector<satelliteSource> &satelliteSources, mt19937 &generator);

extern float checkGridOccupancy(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]);
extern bool checkValidPoint(int x, int y, int z);
extern void matrixDebugTester(float test1, float test2, bool exact);
extern uniform_real_distribution<float>  dist01;

extern float cytokineProductionRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr);
extern float cytokineComboRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr);


void gridCell::endothelium_activate(CytokineGrids &cytoGrids,
float (&RM)[numRules][numRuleParams]){
  int x,y,z;

  x=xLoc;
  y=yLoc;
  z=zLoc;
  ec_roll++;
  ec_stick++;
  cytoGrids.PAF[x][y][z]=cytokineProductionRule(RM, cytoGrids,24,x,y,z,0,0);
  if((collagen3+collagen1)<=viableMuscle){
    cytoGrids.IL8[x][y][z]=cytokineProductionRule(RM, cytoGrids,25,x,y,z,0,0);
  }
  if(woundEdge>0){
    cytoGrids.IFNg[x][y][z]=cytokineProductionRule(RM, cytoGrids,26,x,y,z,0,0);
  }
}

void gridCell::endothelialFunction(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids, vector<pmnSource> &pmnSources, vector<macroSource> &macroSources,
vector<T_Source> &T_Sources, vector<fibroSource> &fibroSources, vector<satelliteSource> &satelliteSources,mt19937 &generator,
float (&RM)[numRules][numRuleParams]){
  float totalTissue,test1,test2;
  int x,y,z,flag,i;
  x=xLoc;
  y=yLoc;
  z=zLoc;
//  cout<<"Test Endothelial Function 1\n";
  if(tissueGrid[x][y][z].contamination>0){
    tissueGrid[x][y][z].contamination=tissueGrid[x][y][z].contamination-cytoGrids.cytotox[x][y][z];
  }
  if(tissueGrid[x][y][z].contamination<0){
    tissueGrid[x][y][z].contamination=0;
  }
  if(tissueGrid[x][y][z].viableMuscle>0){
   tissueGrid[x][y][z].viableMuscle-=cytoGrids.cytotox[x][y][z];
//   cout<<"Inflammatory Muscle Damage\n";
   tissueGrid[x][y][z].necrosis+=cytoGrids.cytotox[x][y][z];
  }

  if(tissueGrid[x][y][z].viableMuscle<0){
    tissueGrid[x][y][z].viableMuscle=0;
  }

  if(tissueGrid[x][y][z].necrosis>1){
    tissueGrid[x][y][z].necrosis=1;
  }
  if(tissueGrid[x][y][z].necrosis<0){
    tissueGrid[x][y][z].necrosis=0;
  }
  if(tissueGrid[x][y][z].life>0){
    tissueGrid[x][y][z].life-=cytoGrids.cytotox[x][y][z]/5;
  }
  if(tissueGrid[x][y][z].collagen3>0){
    tissueGrid[x][y][z].collagen3-=cytoGrids.cytotox[x][y][z]/20;
  }

  if(tissueGrid[x][y][z].collagen3<0){
    tissueGrid[x][y][z].collagen3=0;
  }

  if(tissueGrid[x][y][z].collagen1>0){
    tissueGrid[x][y][z].collagen1-=cytoGrids.cytotox[x][y][z]/20;
  }

  if(tissueGrid[x][y][z].collagen1<0){
    tissueGrid[x][y][z].collagen1=0;
  }

  if(collagen1>0.5&&collagen3>0){
    collagen1+=0.05;
    collagen3-=0.05;
    if(collagen3<0){
      collagen3=0;
    }
  }

  totalTissue=collagen1+necrosis+viableMuscle+collagen3;
  if(totalTissue<=0){contamination=0;}
  if(totalTissue>=0.5&&collagen3>0){
    cytoGrids.FNE[xLoc][yLoc][zLoc]+=collagen3*10;
  }
  if(totalTissue>0){
    cytoGrids.IL6[x][y][z]=cytokineProductionRule(RM, cytoGrids,27,x,y,z,0,0);
    cytoGrids.TGFb[x][y][z]=cytokineProductionRule(RM, cytoGrids,28,x,y,z,0,0);
    if(contamination>0){
      cytoGrids.DAMP[xLoc][yLoc][zLoc]+=5;
    }
    colEdge=checkForMuscleCollagenEdge(tissueGrid);
    woundEdge=checkForWoundEdge(tissueGrid);
    vascularEdge=checkForVascularEdge(tissueGrid);
    muscleEdge=checkForMuscleEdge(tissueGrid);

    if(muscleEdge>0){
  //    cout<<"Muscle Edge at Height "<<z<<"\n";
      cytoGrids.IL4[x][y][z]=cytokineProductionRule(RM, cytoGrids,30,x,y,z,0,0);
    }
    if(woundEdge>0){
      ///ECM/////////
      cytoGrids.IL10[x][y][z]++;
      /////////////////////////
      cytoGrids.TGFb[x][y][z]=cytokineProductionRule(RM, cytoGrids,31,x,y,z,0,0);
      if(contamination>0){
        cytoGrids.IL8[x][y][z]=cytokineProductionRule(RM, cytoGrids,32,x,y,z,0,0);
      } else if (contamination<=0 && (collagen3+collagen1)<=viableMuscle){
        cytoGrids.IL8[x][y][z]=cytokineProductionRule(RM, cytoGrids,33,x,y,z,0,0);
      }
      else{
          if(viableMuscle<1){
                cytoGrids.MCP1[xLoc][yLoc][zLoc]+=0.5;
          }
      }
      if(viableMuscle>(collagen1+collagen3)){
        cytoGrids.TNF[x][y][z]=cytokineProductionRule(RM,cytoGrids,35,x,y,z,0,0);
//        cout<<"Muscle High\n";
      } else{
        cytoGrids.TNF[x][y][z]=cytokineProductionRule(RM,cytoGrids,36,x,y,z,0,0);
//        cout<<"Muscle Low\n";
      }
    }
    if(vascularEdge>0){
      proliferation(tissueGrid, cytoGrids, pmnSources, macroSources, T_Sources, fibroSources,satelliteSources, generator);
    }
    if(life<0.75 ||necrosis>0.25 || contamination>0){
      ec_activation=1;
    }
    if(ec_activation>=1){
      endothelium_activate(cytoGrids,RM);
    }
    injurySpread(oxyHeal, tissueGrid, cytoGrids);
  }
  flag=0;
  if(viableMuscle>=1){
    if(zLoc<zDim-5){
      for(i=1;i<=5;i++){
        if(tissueGrid[xLoc][yLoc][zLoc+1].viableMuscle<0.99){
          flag=1;
        }
      }
      if(flag==0){
        cytoGrids.IL10[xLoc][yLoc][zLoc]+=0.1;
      }

//       if(tissueGrid[xLoc][yLoc][zLoc+1].viableMuscle>=1){
// //        cout<<"IL10 SUPPLEMENTATION\n";
//         cytoGrids.IL10[xLoc][yLoc][zLoc]+=0.1;
//
//       }
    }
  }
}

void gridCell::injurySpread(float oxyHeal,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids){
  int x,y,z;
  x=xLoc;
  y=yLoc;
  z=zLoc;

  if(life<0){
    life=0;
  }
  //Three damage regimes: damaged but healing, ischemia, infarction
  if(life>=0.50){
    life=min(float(1),life+oxyHeal);
  }
  if((life<0.50)&&(life>0.30)){ //ischemia
    ec_roll++;
    life-=0.01;
    injuryPropagation(x,y,z,-0.005, tissueGrid);
  }
  if(life<=0.30){ //infarction
    ec_stick++;
    life-=.02;
    injuryPropagation(x,y,z,-0.01, tissueGrid);
  }
  if(life<=0){
    life=0;
  }

  if(life<necrosis){
    necrosis+=((necrosis-life)/2);
    viableMuscle-=((necrosis-life)/2);
  }
  if(necrosis>=0.60){
    necrosis+=0.02;
    viableMuscle-=0.02;
    cytoGrids.MCP1[x][y][z]+=4;
  }
  if((necrosis<0.60)&&(necrosis>0.25)){ //ischemia
    necrosis+=0.01;
    viableMuscle-=0.01;
    cytoGrids.MCP1[x][y][z]+=2;
  }
  if(necrosis>0 && necrosis<=0.25){ //infarction
    necrosis-=0.005;
    viableMuscle+=0.005;
    cytoGrids.MCP1[x][y][z]+=1;
  }
  if(necrosis>1){
    necrosis=1;
  }
  if(necrosis<0){
    necrosis=0;
  }
  if(viableMuscle<0){
    viableMuscle=0;
  }
  if(viableMuscle>1){
    viableMuscle=1;
  }
}

void injuryPropagation(int x, int y, int z, float delta,
  gridCell (&tissueGrid)[xDim][yDim][zDim]) {
  int i,j,k,ii,jj,kk;
  //  cout<<"Test 51-IP-1\n";
  for(i=x-1;i<=x+1;i++){
    for(j=y-1;j<=y+1;j++){
      for(k=z-1;k<=z+1;k++){
        if(i==x&&j==y&&k==z){
          continue;
        }
        else{
          ii=i;
          jj=j;
          kk=k;
          if(i<0){
            ii=xDim-1;
          }
          if(i>=xDim){
            ii=0;
          }
          if(j<0){
            jj=yDim-1;
          }
          if(j>=yDim){
            jj=0;
          }
          if(k<0){
            continue;
          }
          if(k>=zDim){
            continue;
          }
          tissueGrid[ii][jj][kk].life+=delta;
          if(tissueGrid[ii][jj][kk].life<0){
            tissueGrid[ii][jj][kk].life=0;
          }
        }
      }
    }
  }
}

int gridCell::checkForWoundEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]){
  int i,j,k,flag,ii,jj,kk;
  float cellContent,cellNec;
  flag=0;
  if(zLoc>=zDim-1){
//    cout<<"Re Epithelization Test\n";
    return flag;
  }
  for(i=xLoc-1;i<=xLoc+1;i++){
    for(j=yLoc-1;j<=yLoc+1;j++){
      for(k=zLoc-1;k<=zLoc+1;k++){
        if(i==xLoc&&j==yLoc&&k==zLoc){
          continue;
        }
        else{
          ii=i;
          jj=j;
          kk=k;
          if(i<0){
            ii=xDim-1;
          }
          if(i>=xDim){
            ii=0;
          }
          if(j<0){
            jj=yDim-1;
          }
          if(j>=yDim){
            jj=0;
          }
          if(k<0){
            continue;
          }
          if(k>=zDim){
            continue;
          }
          cellContent=tissueGrid[ii][jj][kk].viableMuscle+tissueGrid[ii][jj][kk].collagen1+tissueGrid[ii][jj][kk].collagen3+tissueGrid[ii][jj][kk].necrosis;
          cellNec=tissueGrid[ii][jj][kk].necrosis;
          if(cellContent<=0){
            flag=1;
            return flag;
          }
        }
      }
    }
  }
  return flag;
}

int gridCell::checkForMuscleCollagenEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]){  //returns 1 if you are next to a muscle voxel or collagen1 voxel
  int i,j,k,flag,ii,jj,kk;
  float cellContent,cellNec;
  flag=0;

  for(i=xLoc-1;i<=xLoc+1;i++){
    for(j=yLoc-1;j<=yLoc+1;j++){
      for(k=zLoc-1;k<=zLoc+1;k++){
        if(i==xLoc&&j==yLoc&&k==zLoc){
          continue;
        }
        else{
          ii=i;
          jj=j;
          kk=k;
          if(i<0){
            ii=xDim-1;
          }
          if(i>=xDim){
            ii=0;
          }
          if(j<0){
            jj=yDim-1;
          }
          if(j>=yDim){
            jj=0;
          }
          if(k<0){
            continue;
          }
          if(k>=zDim){
            continue;
          }
          cellContent=tissueGrid[ii][jj][kk].viableMuscle+tissueGrid[ii][jj][kk].collagen1;
          if(cellContent>=0){
            flag=1;
            return flag;
          }
        }
      }
    }
  }
  return flag;
}

int gridCell::checkForVascularEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]){
  int i,j,k,flag,ii,jj,kk;
  float cellContent,vasc;
  flag=0;

  for(i=xLoc-1;i<=xLoc+1;i++){
    for(j=yLoc-1;j<=yLoc+1;j++){
      for(k=zLoc-1;k<=zLoc+1;k++){
        if(i==xLoc&&j==yLoc&&k==zLoc){
          continue;
        }
        else{
          ii=i;
          jj=j;
          kk=k;
          if(i<0){
            ii=xDim-1;
          }
          if(i>=xDim){
            ii=0;
          }
          if(j<0){
            jj=yDim-1;
          }
          if(j>=yDim){
            jj=0;
          }
          if(k<0){
            continue;
          }
          if(k>=zDim){
            continue;
          }
          cellContent=tissueGrid[ii][jj][kk].viableMuscle+tissueGrid[ii][jj][kk].collagen1+tissueGrid[ii][jj][kk].collagen3+tissueGrid[ii][jj][kk].necrosis;
          vasc=tissueGrid[ii][jj][kk].life;

          if(cellContent>=0.9&&vasc<0.1){
            //            cout<<"Test44 "<<ii<<","<<jj<<","<<kk<<","<<cellContent<<"\n";
            flag=1;
            return flag;
          }
        }
      }
    }
  }
  //  cout<<"Test 46\n";
  return flag;
}

int gridCell::checkForMuscleEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]){
  int yp,flag=0;
  //
  //
  // //for BIOPSY ONLY:
  // if(yLoc==0&&checkGridOccupancy(xLoc,yLoc,zLoc,tissueGrid)<1){
  //   flag=1;
  //   return flag;
  // }
  //
  //
  // //FOR BIOPSY ONLY
  //
  // yp=yLoc+1;
  // if(yp>=yDim){
  //   yp=0;
  // }
  // if(tissueGrid[xLoc][yp][zLoc].viableMuscle>0&&
  //   tissueGrid[xLoc][yLoc][zLoc].viableMuscle<tissueGrid[xLoc][yp][zLoc].viableMuscle){
  //   flag=1;
  //   return flag;
  // }
  // yp=yLoc-1;
  // if(yp<0){
  //   yp=yDim-1;
  // }
  // if(tissueGrid[xLoc][yp][zLoc].viableMuscle>0&&
  //   tissueGrid[xLoc][yLoc][zLoc].viableMuscle<tissueGrid[xLoc][yp][zLoc].viableMuscle){
  //   flag=1;
  //   return flag;
  // }
  // return flag;
  if(tissueGrid[xLoc][yLoc][zLoc].viableMuscle>0&&checkForWoundEdge(tissueGrid)==1){
//    cout<<xLoc<<" "<<yLoc<<" "<<zLoc<<" "<<tissueGrid[xLoc][yLoc][zLoc].viableMuscle<<" "<<checkForWoundEdge(tissueGrid)<<"\n";
    flag=1;
    return flag;
  }
  return flag;
}

void gridCell::proliferation(gridCell (&tissueGrid)[xDim][yDim][zDim], CytokineGrids &cytoGrids, vector<pmnSource> &pmnSources, vector<macroSource> &macroSources, vector<T_Source> &T_Sources, vector<fibroSource> &fibroSources, vector<satelliteSource> &satelliteSources, mt19937 &generator){
  vector<int> indexes;
  int x,y,z,q,i,j,k,tx,ty,tz;
  int candidates[26][3];
  float gridOccupancy,prolif;

  x=xLoc;
  y=yLoc;
  z=zLoc;
  q=0;
  //
  //Finding the neighbors that need to be healed
  for(i=xLoc-1;i<=xLoc+1;i++){
    for(j=yLoc-1;j<=yLoc+1;j++){
      for(k=zLoc-1;k<=zLoc+1;k++){
        if(i==x&&j==y&&k==z){
          continue;
        } else {
          candidates[q][0]=i;
          candidates[q][1]=j;
          candidates[q][2]=k;
          if(checkValidPoint(i,j,k)){
            gridOccupancy=checkGridOccupancy(i,j,k, tissueGrid);
            if(tissueGrid[i][j][k].life<1&&gridOccupancy>0.7){
              indexes.push_back(q);
            }
          }
          q=q+1;
        }
      }
    }
  }
  //  cout<<"Index Size="<<indexes.size()<<"\n";
  if(indexes.size()>0){
    prolif=(cytoGrids.VEGF[x][y][z]+cytoGrids.IL17[x][y][z]-cytoGrids.IFNg[x][y][z]);
    //prolif=1;
    if(prolif>0){
      for(i=0;i<indexes.size();i++){
        tx=candidates[indexes[i]][0];
        ty=candidates[indexes[i]][1];
        tz=candidates[indexes[i]][2];
        if(tissueGrid[tx][ty][tz].life<1){
          tissueGrid[tx][ty][tz].life+=0.5;
          if(counter<1){
            sourceGeneration(tx,ty,tz, tissueGrid, pmnSources, macroSources, T_Sources, fibroSources, satelliteSources, generator);
            tissueGrid[tx][ty][tz].counter++;
          }
          if(tissueGrid[tx][ty][tz].life>1){
            tissueGrid[tx][ty][tz].life=1;
          }
        }
      }
    }
  }
}

void sourceGeneration(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim], vector<pmnSource> &pmnSources, vector<macroSource> &macroSources, vector<T_Source> &T_Sources, vector<fibroSource> &fibroSources, vector<satelliteSource> &satelliteSources, mt19937 &generator){
  float temp,vm;
  temp=dist01(generator);
//  vm=1.0;

  vm=tissueGrid[x][y][z].viableMuscle;
//  cout<<"SGTEST "<<x<<" "<<y<<" "<<z<<" "<<tissueGrid[x][y][z].viableMuscle<<" "<<tissueGrid[x][y][z].collagen3<<"\n";
  if(temp>0.05 && temp<0.0505){
    pmnSources.push_back(pmnSource(x,y,z));
  }
  if(temp>0.1 && temp<0.103){
    macroSources.push_back(macroSource(x,y,z));
  }
  if(temp>0.15 && temp<0.150001){
    T_Sources.push_back(T_Source(x,y,z));
  }
  if(temp>0.2 && temp<0.203){
    fibroSources.push_back(fibroSource(x,y,z));
  }
//  if(temp>0.203 && temp<0.206){
  if((temp>0.203 && temp<0.21)&&(vm>=0.99)){
    satelliteSources.push_back(satelliteSource(x,y,z));
  }
}
