#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

extern uniform_int_distribution<int>  distributionEdges,distributionX,distributionY,distributionZ,distribution27;
extern uniform_real_distribution<float>  dist01;

extern void matrixDebugTester(float test1, float test2, bool exact);

extern float cytokineProductionRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr);

pmn::pmn(){}

pmn::pmn(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  wbc_migrate=0;
  wbc_stick=0;
  wbc_roll=0;
  age=20;
}

void pmn::pmn_function(int pmnID,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
  vector<pmn> &pmns,
  int &bursts,
  mt19937 &generator,
  float (&RM)[numRules][numRuleParams]
){
  int x,y,z,roll;
  float tnf,paf,il1,test;
  x=xLoc;
  y=yLoc;
  z=zLoc;
  sniff(tissueGrid, cytoGrids, generator);
  sniff(tissueGrid, cytoGrids, generator);
//  cout<<"MIGRATE="<<wbc_migrate<<"\n";
  //  if(wbc_migrate>0||tissueGrid[x][y][z].contamination>0.0
  if(wbc_migrate>0||tissueGrid[x][y][z].contamination>0.0){
  //if(wbc_migrate>0){
//    if(tissueGrid[x][y][z].viableMuscle<0.99){
//      cout<<tissueGrid[x][y][z].viableMuscle<<" "<<tissueGrid[x][y][z].necrosis<<"\n";
      // cout<<"PMN Coords="<<x<<" "<<y<<" "<<z<<" "<<wbc_migrate<<"\n";
      // cout<<"DAMP="<<cytoGrids.DAMP[x][y][z]<<"\n";
      // cout<<cytoGrids.TNF[x][y][z]<<" "<<cytoGrids.IL1[x][y][z]<<" "<<cytoGrids.IL10[x][y][z]<<" "<<tissueGrid[x][y][z].contamination<<"\n";
      pmn_burst(pmnID, tissueGrid, cytoGrids, pmns, bursts, RM);
//    }
  } else {
    roll=tissueGrid[x][y][z].ec_roll;
    if((roll>3)&&(wbc_roll==1)){
      sniff(tissueGrid, cytoGrids, generator);
    }
    else{
      sniff(tissueGrid, cytoGrids, generator);
      sniff(tissueGrid, cytoGrids, generator);
    }
    x=xLoc;
    y=yLoc;
    z=zLoc;
    tnf=cytoGrids.TNF[x][y][z];
    paf=cytoGrids.PAF[x][y][z];
    il1=cytoGrids.IL1[x][y][z];
    if(tnf+paf>1){
      //wbc_stick=tnf+paf;
      wbc_stick=cytokineProductionRule(RM,cytoGrids,37,x,y,z,0,0);
    }
    if(wbc_stick>=1&&tissueGrid[x][y][z].ec_stick>=1){
      //wbc_migrate=max(float(0),(cytoGrids.TNF[x][y][z]+cytoGrids.IL1[x][y][z]+cytoGrids.DAMP[x][y][z]-cytoGrids.IL10[x][y][z]));
      wbc_migrate=max(float(0),cytokineProductionRule(RM,cytoGrids,38,x,y,z,0,0));
    }
    age--;
    if(tnf+paf+il1==0){age--;}
    if(age<0){
      x=xLoc;
      y=yLoc;
      z=zLoc;
      pmns.erase(pmns.begin()+pmnID);
      tissueGrid[x][y][z].cellPopulation--;
      return;
    }
  }
}

void pmn::pmn_burst(int pmnID,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
  vector<pmn> &pmns,
  int &bursts,
  float (&RM)[numRules][numRuleParams]
){
  int x,y,z;
  float tnf,ifng,il10,gcsf,pmn_pcd,test1,test2;
  pmn_pcd=0;
  x=xLoc;
  y=yLoc;
  z=zLoc;
  tissueGrid[x][y][z].ec_stick=0;
  tissueGrid[x][y][z].ec_roll=0;
  tissueGrid[x][y][z].ec_migrate=0;

  bursts++;
  //cytoGrids.cytotox[x][y][z]=min(float(1),cytoGrids.TNF[x][y][z]+cytoGrids.DAMP[x][y][z]);  //cytotox can be treated as ROS
  cytoGrids.cytotox[x][y][z]=min(float(1),cytokineProductionRule(RM,cytoGrids,39,x,y,z,0,0));
  //cytoGrids.TGFb[x][y][z]=cytoGrids.TGFb[x][y][z]+float(1)+10*cytoGrids.cytotox[x][y][z]+cytoGrids.IL6[x][y][z];
  cytoGrids.TGFb[x][y][z]=cytokineProductionRule(RM,cytoGrids,40,x,y,z,0,0);
  //cytoGrids.IL10[x][y][z]=cytoGrids.IL10[x][y][z]+cytoGrids.TNF[x][y][z]+cytoGrids.IL1[x][y][z]+cytoGrids.IL6[x][y][z];
  cytoGrids.IL10[x][y][z]=cytokineProductionRule(RM,cytoGrids,41,x,y,z,0,0);
  //cytoGrids.TNF[x][y][z]+=4;
  cytoGrids.TNF[x][y][z]=cytokineProductionRule(RM,cytoGrids,42,x,y,z,0,0);
  //cytoGrids.IL1[x][y][z]=cytoGrids.IL1[x][y][z]+1;
  cytoGrids.IL1[x][y][z]=cytokineProductionRule(RM,cytoGrids,43,x,y,z,0,0);
  cytoGrids.IL6[x][y][z]+=0.5*tissueGrid[x][y][z].contamination;

  test1=cytokineProductionRule(RM, cytoGrids,44,x,y,z,0,0);
  test2=cytoGrids.MCP1[x][y][z]+2*cytoGrids.TNF[x][y][z]+cytoGrids.IL1[x][y][z];
    //cytoGrids.MCP1[x][y][z]+=2*cytoGrids.TNF[x][y][z]+cytoGrids.IL1[x][y][z];
  // if(test1!=test2){
  //   cout<<"Matrix Rule="<<test1<<"\n";
  //   cout<<"Coded Rule="<<test2<<"\n";
  //   cout<<"DEBUGS="<<cytoGrids.MCP1[xLoc][yLoc][zLoc]<<" "<<cytoGrids.TNF[xLoc][yLoc][zLoc]<<" "<<cytoGrids.IL1[xLoc][yLoc][zLoc]<<"\n";
  //   exit(0);
  // }
//  if(tissueGrid[x][y][z].viableMuscle<0.99){
//    cout<<"PMN MCP SECRETOION "<<z<<"\n";
    cytoGrids.MCP1[x][y][z]=cytokineProductionRule(RM,cytoGrids,44,x,y,z,0,0);
//  }
  tnf=cytoGrids.TNF[x][y][z];
  ifng=cytoGrids.IFNg[x][y][z];
  il10=cytoGrids.IL10[x][y][z];
  gcsf=cytoGrids.GCSF[x][y][z];
  //pmn_pcd+=max(float(0),(tnf+ifng+gcsf-il10)/100);
  pmn_pcd+=max(float(0),cytokineProductionRule(RM,cytoGrids,45,x,y,z,0,0)/100);
  //  cout<<"PCD="<<pmn_pcd<<","<<tnf<<","<<ifng<<","<<gcsf<<","<<il10<<","<<(tnf+ifng+gcsf-il10)/100<<","<<max(float(0),(tnf+ifng+gcsf-il10)/100)<<"\n";
  age+=pmn_pcd;
  age--;
  if(age<0){
    //    cout<<"Age Death\n";
    x=xLoc;
    y=yLoc;
    z=zLoc;
    pmns.erase(pmns.begin()+pmnID);
    tissueGrid[x][y][z].cellPopulation--;
  }
}

void pmn::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
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
          if(tissueGrid[ii][jj][kk].cellPopulation>=cellCapacity){
            chemotaxis[m]=-1;
          } else{
            pop=tissueGrid[ii][jj][kk].cellPopulation;
            chemotaxis[m]=(cytoGrids.DAMP[ii][jj][kk])*(cellCapacity-pop);
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
  //cout<<"Climb Test "<<zLoc<<","<<z<<"\n";
  if(z<zLoc && cytoMax>0){
    // cout<<"Climb Test"<<DAMP[xLoc][yLoc][zLoc]<<","<<DAMP[x][y][z]<<","<<tissueGrid[xLoc][yLoc][zLoc].necrosis<<","<<tissueGrid[x][y][z].necrosis<<"\n";
    // cout<<"CT2\n";
  }
  if(tissueGrid[x][y][z].cellPopulation<cellCapacity&&tissueGrid[x][y][z].life>=0){
    tissueGrid[x][y][z].cellPopulation++;
    tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
    xLoc=x;
    yLoc=y;
    zLoc=z;
  } else {
    //    cout<<"PMN Movement Error " <<tissueGrid[x][y][z].cellPopulation<<"    "<<tissueGrid[x][y][z].vascularization<<"\n";
  }
}

void pmnSource::proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
vector<pmn>& pmns,
mt19937 &generator){
  int pop,i,n;
  float temp,limit,signal,ratio,iratio;
  pop=tissueGrid[xLoc][yLoc][zLoc].cellPopulation;
  if(pop<cellCapacity){
    signal=cytoGrids.PAF[xLoc][yLoc][zLoc];
    if(signal>0){
      activation++;
    } else {
      activation--;
    }
    if(activation<0){
      activation=0;
    }
    iratio=float(activation)/float(pmnProliferationActivationThreshold);
    if(iratio>1){
      iratio=1;
    }
    ratio=4*(iratio)-2.0;
    if(cytoGrids.DAMP[xLoc][yLoc][zLoc]>0){
      //cout<<"Damp Proliferation\n";
      limit=0.1*(0.5*(erf(ratio)+1));
    } else  if (cytoGrids.TNF[xLoc][yLoc][zLoc]-2*cytoGrids.IL10[xLoc][yLoc][zLoc]>0){
//      cout<<"TNF Proliferation\n";
      limit=0.03*(0.5*(erf(ratio)+1));
    } else {
      limit=0.0001*(0.5*(erf(ratio)+1));
    }
    temp=dist01(generator);
    if(temp<limit){
      pmns.push_back(pmn(xLoc,yLoc,zLoc));
      tissueGrid[xLoc][yLoc][zLoc].cellPopulation++;
    }
  }
}
