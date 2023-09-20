#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

extern float cytokineProductionRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr);
extern float cytokineComboRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr);
extern uniform_int_distribution<int>  distributionEdges,distributionX,distributionY,distributionZ,distribution27;
extern uniform_real_distribution<float>  dist01;

extern void matrixDebugTester(float test1, float test2, bool exact);

float getTotalDAMP(gridCell(&tissueGrid)[xDim][yDim][zDim],CytokineGrids &cytoGrids);

macrophage::macrophage(){}

macrophage::macrophage(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  TNFr=0;
  IL1r=0;
  activation=0;
  age=20;
  wbc_stick=0;
  wbc_migrate=0;
  wbc_roll=0;
}

void macrophage::action(int index, gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
vector<macrophage>& macrophages,
mt19937 &generator,
float (&RM)[numRules][numRuleParams]
){
  int x,y,z;
  float test,col3factor,col1factor,test1,test2;

  x=xLoc;
  y=yLoc;
  z=zLoc;

  col3factor=1.0-tissueGrid[x][y][z].collagen3;
  col1factor=1.0-tissueGrid[x][y][z].collagen1;

  if(cytoGrids.sTNFr[x][y][z]<=100){
    TNFr=min(float(100),(cytoGrids.TNF[x][y][z]+cytoGrids.sTNFr[x][y][z]));
  }
  else{
    TNFr=min(float(100),max(float(0),cytoGrids.TNF[x][y][z]-cytoGrids.sTNFr[x][y][z]));
  }
  IL1r=min(float(1),max(float(0),cytokineProductionRule(RM,cytoGrids,0,x,y,z,IL1r,TNFr)));
  cytoGrids.IL1ra[x][y][z]=cytokineProductionRule(RM,cytoGrids,1,x,y,z,IL1r,TNFr);
  cytoGrids.sTNFr[x][y][z]=cytoGrids.sTNFr[x][y][z]+TNFr/2;
  cytoGrids.sIL1r[x][y][z]=cytokineProductionRule(RM,cytoGrids,2,x,y,z,IL1r,TNFr);
  //activation=cytoGrids.MCP1[x][y][z]+cytoGrids.PAF[x][y][z]+cytoGrids.IFNg[x][y][z]-cytoGrids.IL10[x][y][z];
  activation=cytokineComboRule(RM,cytoGrids,3,x,y,z,IL1r,TNFr);

  if(activation>0){
    if(tissueGrid[x][y][z].necrosis>0){
      tissueGrid[x][y][z].necrosis=tissueGrid[x][y][z].necrosis-macroNecroClearance;
    }
    if(tissueGrid[x][y][z].necrosis<=0){
      tissueGrid[x][y][z].necrosis=0;
    }
    cytoGrids.TGFb[x][y][z]=cytokineProductionRule(RM,cytoGrids,4,x,y,z,IL1r,TNFr);
    //cytoGrids.TGFb[x][y][z]=cytoGrids.TGFb[x][y][z]-cytoGrids.DAMP[x][y][z]+float(1)+cytoGrids.cytotox[x][y][z];
    if(cytoGrids.TGFb[x][y][z]<0){
      cytoGrids.TGFb[x][y][z]=0;
    }
    cytoGrids.GCSF[x][y][z]=cytokineProductionRule(RM, cytoGrids, 5,x,y,z,IL1r,TNFr);
    //cytoGrids.GCSF[x][y][z]=cytoGrids.GCSF[x][y][z]+cytoGrids.endotoxin[x][y][z]+cytoGrids.PAF[x][y][z]+cytoGrids.TNF[x][y][z]+cytoGrids.IFNg[x][y][z];
    cytoGrids.IL12[x][y][z]=cytokineProductionRule(RM,cytoGrids,6,x,y,z,IL1r,TNFr);
    //cytoGrids.IL12[x][y][z]=cytoGrids.IL12[x][y][z]+cytoGrids.TNF[x][y][z]+cytoGrids.IL1[x][y][z];
    cytoGrids.IL10[x][y][z]=cytokineProductionRule(RM,cytoGrids,7,x,y,z,IL1r,TNFr);
    //cytoGrids.IL10[x][y][z]=(cytoGrids.IL10[x][y][z]+cytoGrids.TNF[x][y][z]+cytoGrids.IL1[x][y][z]+cytoGrids.IL6[x][y][z]+10*cytoGrids.cytotox[x][y][z]);
    //cytoGrids.IL4[x][y][z]=cytokineProductionRule(RM,cytoGrids,8,x,y,z,IL1r,TNFr);
    //cytoGrids.IL4[x][y][z]=cytoGrids.IL4[x][y][z]+cytoGrids.MCP1[x][y][z];
    cytoGrids.IL1[x][y][z]=cytokineProductionRule(RM,cytoGrids,9,x,y,z,IL1r,TNFr);

    if(tissueGrid[x][y][z].collagen1+tissueGrid[x][y][z].collagen3<1){
      cytoGrids.TNF[x][y][z]=cytoGrids.TNF[x][y][z]+(cytokineProductionRule(RM, cytoGrids, 10,x,y,z,IL1r,TNFr))*0.1;
      cytoGrids.IL8[x][y][z]=cytoGrids.IL8[x][y][z]+(cytokineProductionRule(RM, cytoGrids, 11,x,y,z,IL1r,TNFr))*0.5;
      if(cytoGrids.TNF[x][y][z]<0){cytoGrids.TNF[x][y][z]=0;}
      if(cytoGrids.IL8[x][y][z]<0){cytoGrids.IL8[x][y][z]=0;}
    } else {
      cytoGrids.TNF[x][y][z]=cytoGrids.TNF[x][y][z]+(cytokineProductionRule(RM, cytoGrids, 10,x,y,z,IL1r,TNFr))*0.025;
      cytoGrids.IL8[x][y][z]=cytoGrids.IL8[x][y][z]+(cytokineProductionRule(RM, cytoGrids, 11,x,y,z,IL1r,TNFr))*0.05;
      if(cytoGrids.TNF[x][y][z]<0){cytoGrids.TNF[x][y][z]=0;}
      if(cytoGrids.IL8[x][y][z]<0){cytoGrids.IL8[x][y][z]=0;}
    }
//    if(tissueGrid[x][y][z].viableMuscle<0.99){
//      cout<<"MACRO MCP SECRETION "<<z<<" "<<tissueGrid[x][y][z].viableMuscle<<"\n";
      cytoGrids.MCP1[x][y][z]=cytokineProductionRule(RM, cytoGrids, 12,x,y,z,IL1r,TNFr);
//    }
    cytoGrids.PDGF[x][y][z]=cytokineProductionRule(RM, cytoGrids, 13,x,y,z,IL1r,TNFr);

    if((wbc_stick==1)&&(tissueGrid[x][y][z].ec_stick>=100)){
      wbc_migrate=1;
    }
    if(wbc_roll==1){
      wbc_stick=1;
    }
    wbc_roll=1;
  }
  if(activation<0){
    cytoGrids.IL10[x][y][z]=cytokineProductionRule(RM,cytoGrids,14,x,y,z,IL1r,TNFr);
    cytoGrids.TGFb[x][y][z]+=1;
  }
  if(wbc_roll==1){
    sniff(tissueGrid, cytoGrids, generator);
  }
  else{
    sniff(tissueGrid, cytoGrids, generator);
    sniff(tissueGrid, cytoGrids, generator);
  }
  age--;
  if(age<0){
    x=xLoc;
    y=yLoc;
    z=zLoc;
    macrophages.erase(macrophages.begin()+index);
    tissueGrid[x][y][z].cellPopulation--;
    return;
  }
  if(activation>20){
    activation=20;
  }
}

void macrophage::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
 CytokineGrids &cytoGrids,
 mt19937 &generator){
  int x,y,z,i,j,k,m,nextLocIndex,ii,jj,kk,temp,pop;
  float chemotaxis[27];
  float cytoMax=-1.0;
  int coords[27][3];
  x=xLoc;
  y=yLoc;
  z=zLoc;
  //sensingcytokine concentration in surrounding environment; motion outside the grid is prohibited by the '-1' value
  m=0;
  for(i=x-1;i<=x+1;i++){
    for(j=y-1;j<=y+1;j++){
      for(k=z-1;k<=z+1;k++){
        ii=i;
        jj=j;
        kk=k;
        if(ii<0){
          ii=xDim-1;
        }
        if(ii>=xDim){
          ii=0;
        }
        if(jj<0){
          jj=yDim-1;
        }
        if(jj>=yDim){
          jj=0;
        }
        if(kk<0){
          chemotaxis[m]=-1;
        } else if(kk>=zDim){
          chemotaxis[m]=-1;
        } else if(kk>=0&&kk<zDim){
          pop=tissueGrid[ii][jj][kk].cellPopulation;
          if(tissueGrid[ii][jj][kk].cellPopulation>=cellCapacity){
            chemotaxis[m]=-1;
          } else{
            chemotaxis[m]=(cytoGrids.PAF[ii][jj][kk]+cytoGrids.MCP1[ii][jj][kk])*(cellCapacity-pop);
          }
        }
        coords[m][0]=ii;
        coords[m][1]=jj;
        coords[m][2]=kk;
        m=m+1;
      }
    }
  }
  //finding maximum cytokine concentration for gradient following
  for(i=0;i<27;i++){
    if(chemotaxis[i]>cytoMax){
      cytoMax=chemotaxis[i];
      nextLocIndex=i;
    }
  }
  if(cytoMax>0){
    x=coords[nextLocIndex][0];
    y=coords[nextLocIndex][1];
    z=coords[nextLocIndex][2];
  } else {
    temp=distribution27(generator);
    x=coords[temp][0];
    y=coords[temp][1];
    z=coords[temp][2];
    if(z<0){z=0;}
    if(z>=zDim){z=zDim-1;}
  }
  if(tissueGrid[x][y][z].cellPopulation<cellCapacity){
    tissueGrid[x][y][z].cellPopulation++;
    tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
    xLoc=x;
    yLoc=y;
    zLoc=z;
  }
}

void macroSource::proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
vector<macrophage>& macrophages,
mt19937 &generator){
  int pop,i,n;
  float temp;
  float signal,baselimit,limit,ratio,iratio,totaldamp;

  baselimit=0.02;

  totaldamp=getTotalDAMP(tissueGrid,cytoGrids);
  if(totaldamp>1){
//    baselimit/=totaldamp*float(0.001);
  baselimit/=totaldamp;
  }
  //cout<<"baselimit="<<baselimit<<" "<<totaldamp<<"\n";
  if(tissueGrid[xLoc][yLoc][zLoc].collagen1>0.5){
    activation=0;
    return;
  }

  pop=tissueGrid[xLoc][yLoc][zLoc].cellPopulation;
  signal=cytoGrids.MCP1[xLoc][yLoc][zLoc];
  if(signal>0){
    activation++;
  } else {
    activation--;
  }
  if(activation<0){
    activation=0;
    return;
  }
  iratio=float(activation)/float(macroProliferationActivationThreshold);
  if(iratio>1){
    iratio=1;
  }
  ratio=4*(iratio)-2.0;
  if(pop<2){
    temp=dist01(generator);
    limit=baselimit*(0.5*(erf(ratio)+1)); //error function shifted by 2 gives a sigmoid-like curve from x=-2,2 with -1<y<1
    if(temp<limit){
      macrophages.push_back(macrophage(xLoc,yLoc,zLoc));
      tissueGrid[xLoc][yLoc][zLoc].cellPopulation++;
    }
  }
}

void macrophage::heal(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]){
  tissueGrid[x][y][z].life+=0.1;
  if(tissueGrid[x][y][z].life>1){
    tissueGrid[x][y][z].life=1;
  }
}

float getTotalDAMP(gridCell(&tissueGrid)[xDim][yDim][zDim],CytokineGrids &cytoGrids){
  int i,j,k;
  float sum;
  for(i=0;i<xDim;i++){
    for(j=0;j<yDim;j++){
      for(k=0;k<zDim;k++){
        sum+=cytoGrids.DAMP[i][j][k];
      }
    }
  }
  return sum;
}
