#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

extern uniform_int_distribution<int>  distributionEdges,distributionX,distributionY,distributionZ,distribution27;
extern uniform_real_distribution<float>  dist01;

extern float checkGridOccupancy(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]);
extern bool checkValidPoint(int x, int y, int z);
extern float checkGridOccupancyC3(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]);

extern float cytokineProductionRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids,int ruleRow, int x, int y, int z, float il1r, float tnfr);
extern void matrixDebugTester(float test1, float test2, bool exact);

fibroblast::fibroblast(){}

fibroblast::fibroblast(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  isMyofibroblast=false;
  sma=0;
  age=50;
}

void fibroblast::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
mt19937 &generator
){
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
            if(isMyofibroblast==false){
              chemotaxis[m]=(cytoGrids.IL4[ii][jj][kk]+cytoGrids.IL1[ii][jj][kk]+cytoGrids.IFNg[ii][jj][kk])*(cellCapacity-pop)+cytoGrids.IL8[ii][jj][kk];
            } else {
              chemotaxis[m]=(cytoGrids.FNE[ii][jj][kk])*(cellCapacity-pop);
              //              cout<<"FNE="<<zLoc<<" "<<FNE[ii][jj][kk]<<"\n";
            }
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
  // updating the cell population grid with this pMN's new location
  if(tissueGrid[x][y][z].cellPopulation<cellCapacity){
    //	if(tissueGrid[x][y][z].cellPopulation<cellCapacity&&tissueGrid[x][y][z].vascularization>=0.3){
    tissueGrid[x][y][z].cellPopulation++;
    tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
    xLoc=x;
    yLoc=y;
    zLoc=z;
  }
}

//Fibroblasts undergo chemotaxis.  Next, they sense the state of their location,
//meaning determining how much muscle/collagen is present at the current location.
//If there is space, the fibroblast will deposit collagen.
//Differentiation state is controlled by the boolean, "isMyofibroblast,"
//and the fibroblast differentiates in response to TGFb.  Fibroblasts deposit
//collagen3, while myofibroblasts deposit collagen1 (scar)

//Possible alternative behavior: Deposit collagen  no matter what and have
//it diffuse

// void fibroblast::action(int index,
// gridCell (&tissueGrid)[xDim][yDim][zDim],
// CytokineGrids &cytoGrids,
// vector<fibroblast> &fibroblasts,
// mt19937 &generator,


void fibroblast::action(int index,
gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
vector<fibroblast> &fibroblasts,
mt19937 &generator,int step,
float (&RM)[numRules][numRuleParams]){
  float gridOccupancy,hgfadd,test,goc3;
  int x,y,z;
  int fa,fc,fc2;
  float test1, test2;

  fc=0;
  fc2=0;
  sniff(tissueGrid, cytoGrids, generator);
  x=xLoc;
  y=yLoc;
  z=zLoc;

  gridOccupancy=checkGridOccupancy(xLoc,yLoc,zLoc, tissueGrid);
  activation=cytokineProductionRule(RM,cytoGrids,48,x,y,z,0,0);
  //test1=2*cytoGrids.TGFb[xLoc][yLoc][zLoc]-cytoGrids.IL1[xLoc][yLoc][zLoc]-cytoGrids.IFNg[xLoc][yLoc][zLoc];
  test1=2*cytoGrids.TGFb[xLoc][yLoc][zLoc]-cytoGrids.IL1[xLoc][yLoc][zLoc]-1.5*cytoGrids.IFNg[xLoc][yLoc][zLoc];
  if(test1<0){
    test1=0;
  }
//  test2=cytokineProductionRule(RM,cytoGrids,49,x,y,z,0,0);
  //test2=cytokineProductionRule(RM,cytoGrids,49,x,y,z,0,0)-cytoGrids.IFNg[xLoc][yLoc][zLoc];

  // sma+=cytokineProductionRule(RM,cytoGrids,49,x,y,z,0,0);
    sma+=test1;
  if(tissueGrid[xLoc][yLoc][zLoc].collagen3>0.5){
    sma++;
  }
  if(sma<0){sma=0;}
  if(sma>myoFibroThreshold){
    isMyofibroblast=true;
  }
  if(activation>fibroActivationThreshold){
    cytoGrids.VEGF[xLoc][yLoc][zLoc]=cytokineProductionRule(RM,cytoGrids,51,x,y,z,0,0);
    //hgfadd=cytoGrids.IL1[xLoc][yLoc][zLoc]+cytoGrids.IL6[xLoc][yLoc][zLoc]+cytoGrids.TNF[xLoc][yLoc][zLoc]-cytoGrids.TGFb[xLoc][yLoc][zLoc];
    hgfadd=cytokineProductionRule(RM,cytoGrids,52,x,y,z,0,0);
    if(hgfadd<0){
      hgfadd=0;
    }
    cytoGrids.HGF[xLoc][yLoc][zLoc]+=hgfadd;
    if(isMyofibroblast==false){
      cytoGrids.TGFb[xLoc][yLoc][zLoc]=cytokineProductionRule(RM,cytoGrids,50,x,y,z,0,0);
      if(gridOccupancy<1){
        fillCollagen(gridOccupancy, tissueGrid, cytoGrids);
      } else {
        fa=fillAdjacentF(tissueGrid, cytoGrids);
        if(fa==0){
          while(fa==0){
            fc++;
            mod_sniff(tissueGrid, generator);
            fa=fillAdjacentF(tissueGrid, cytoGrids);
            if(fc==20){
              sma+=2;
              fa=420;
            }
          }
        }
      }
    } else { //it is is a myofibrolast
      cytoGrids.TGFb[xLoc][yLoc][zLoc]+=1;
      goc3=checkGridOccupancyC3(xLoc,yLoc,zLoc, tissueGrid);
      if(gridOccupancy>=1 && goc3>0){
        fillCollagen(gridOccupancy, tissueGrid, cytoGrids);
      } else {
       fa=fillAdjacentM(tissueGrid, cytoGrids);
        fc2++;
        if(fa==0){
          mod_sniff(tissueGrid, generator);
          fa=fillAdjacentM(tissueGrid, cytoGrids);
          if(fc2==10){
            age--;
            fa=420;
          }
        }
      }
    }
  } else {
//    cout<<"activation="<<activation<<" "<<xLoc<<" "<<yLoc<<" "<<zLoc<<" "<<checkGridOccupancy(xLoc,yLoc,zLoc,tissueGrid)<<"\n";
    sniff(tissueGrid, cytoGrids, generator);
  }
  death(index, tissueGrid, cytoGrids, fibroblasts);
}

int fibroblast::fillAdjacentF(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids){
  vector<int> indexes;
  int x,y,z,q,i,j,k,newX,newY,newZ;
  int candidates[26][3];
  float gridOccupancy;
  float tissueMin=1.0;
  float tempCol,maxCol,colDep;
  int indexMin=-1;
  int flag=0;

  x=xLoc;
  y=yLoc;
  z=zLoc;
  q=0;
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
            if(gridOccupancy<1){
              indexes.push_back(q);
              if(gridOccupancy<=tissueMin){
                tissueMin=gridOccupancy;
                indexMin=q;
              }
            }
          }
          q=q+1;
        }
      }
    }
  }
  if(indexes.size()>0){
    //////////////////////NEW METHOD - Collagen goes to least full grid cell
    newX=candidates[indexMin][0];
    newY=candidates[indexMin][1];
    newZ=candidates[indexMin][2];
    if(isMyofibroblast==false){
      tempCol=max(float(0),(5*cytoGrids.TGFb[xLoc][yLoc][zLoc]-cytoGrids.IL1[xLoc][yLoc][zLoc]))*fibroCol3Deposition;
      if(tempCol<1){
//      cout<<"tempCOlA="<<tempCol<<"\n";
      }
      tissueGrid[newX][newY][newZ].collagen3+=tempCol;
      gridOccupancy=checkGridOccupancy(newX,newY,newZ, tissueGrid);
      if(gridOccupancy>1){
        tissueGrid[newX][newY][newZ].collagen3-=(gridOccupancy-1);
      }
    } else { //it is a myifibroblase
      maxCol=tissueGrid[newX][newY][newZ].collagen3;
      if(maxCol>0){
        colDep=min(maxCol,fibroCol1Deposition);
        tissueGrid[newX][newY][newZ].collagen3-=colDep;
        if(tissueGrid[newX][newY][newZ].collagen3<0){
          tissueGrid[newX][newY][newZ].collagen3=0;
        }
        tissueGrid[newX][newY][newZ].collagen1+=colDep;
      }

      gridOccupancy=checkGridOccupancy(newX,newY,newZ, tissueGrid);
      if(gridOccupancy>1){
        tissueGrid[newX][newY][newZ].collagen1-=(gridOccupancy-1);
      }
    }
    if(tissueGrid[newX][newY][newZ].cellPopulation<cellCapacity){
      tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
      tissueGrid[newX][newY][newZ].cellPopulation++;
      xLoc=newX;
      yLoc=newY;
      zLoc=newZ;
    }
    flag=1;
  }
  indexes.clear();
  return flag;
}

int fibroblast::fillAdjacentM(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids){
  vector<int> indexes;
  int x,y,z,q,i,j,k,newX,newY,newZ;
  int candidates[26][3];
  float gridOccupancy,goc3;
  float tissueMin=1.0;
  float tempCol,maxCol,colDep;
  int indexMin=-1;
  int flag=0;

  x=xLoc;
  y=yLoc;
  z=zLoc;
  q=0;
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
            goc3=checkGridOccupancyC3(xLoc,yLoc,zLoc, tissueGrid);
            if(gridOccupancy>=1 || goc3>0){
              indexes.push_back(q);
              if(gridOccupancy<=tissueMin){
                tissueMin=gridOccupancy;
                indexMin=q;
              }
            }
          }
          q=q+1;
        }
      }
    }
  }
  if(indexes.size()>0 && indexMin>=0){
    //////////////////////NEW METHOD - Collagen goes to least full grid cell
    newX=candidates[indexMin][0];
    newY=candidates[indexMin][1];
    newZ=candidates[indexMin][2];
    if(isMyofibroblast==false){
      tempCol=max(float(0),(5*cytoGrids.TGFb[xLoc][yLoc][zLoc]-cytoGrids.IL1[xLoc][yLoc][zLoc]))*fibroCol3Deposition;
      if(tempCol<1){
//      cout<<"tempCOlA="<<tempCol<<"\n";
      }
      tissueGrid[newX][newY][newZ].collagen3+=tempCol;
      gridOccupancy=checkGridOccupancy(newX,newY,newZ, tissueGrid);
      if(gridOccupancy>1){
        tissueGrid[newX][newY][newZ].collagen3-=(gridOccupancy-1);
      }
    } else { //it is a myifibroblase
      maxCol=tissueGrid[newX][newY][newZ].collagen3;
      if(maxCol>0){
        colDep=min(maxCol,fibroCol1Deposition);
        tissueGrid[newX][newY][newZ].collagen3-=colDep;
        if(tissueGrid[newX][newY][newZ].collagen3<0){
          tissueGrid[newX][newY][newZ].collagen3=0;
        }
        tissueGrid[newX][newY][newZ].collagen1+=colDep;
//        cout<<"Adj Fill\n";
      }

      gridOccupancy=checkGridOccupancy(newX,newY,newZ, tissueGrid);
      if(gridOccupancy>1){
        tissueGrid[newX][newY][newZ].collagen1-=(gridOccupancy-1);
      }
    }
    if(tissueGrid[newX][newY][newZ].cellPopulation<cellCapacity){
      tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
      tissueGrid[newX][newY][newZ].cellPopulation++;
      xLoc=newX;
      yLoc=newY;
      zLoc=newZ;
    }
    flag=1;
  }
  indexes.clear();
  return flag;
}

//fibroblasts experience more rapid aging in the presence of TNF, but apoptosis
//is delayed by TGFb
void fibroblast::death(int index,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
vector<fibroblast> &fibroblasts){
  float factor;
  age--;
  //	age-=(TNF[xLoc][yLoc][zLoc]*fibro_tnf_apoptosisFactor-TGFb[xLoc][yLoc][zLoc]*fibro_tgfb_apoptosisFactor);
  factor=(cytoGrids.TNF[xLoc][yLoc][zLoc]*fibro_tnf_apoptosisFactor-cytoGrids.TGFb[xLoc][yLoc][zLoc]*fibro_tgfb_apoptosisFactor);
  if(index<100){
    //    cout<<"FACTOR="<<factor<<"\n";
  }
  if(age<=0){
    tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
    fibroblasts.erase(fibroblasts.begin()+index);
  }
}

void fibroblast::fillCollagen(float gridOccupancy,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids){
  float tempCol,maxCol,colDep;
  if(isMyofibroblast==false){
    tempCol=max(float(0),(5*cytoGrids.TGFb[xLoc][yLoc][zLoc]-cytoGrids.IL1[xLoc][yLoc][zLoc]))*fibroCol3Deposition;
    if(tempCol<1){
//    cout<<"tempCOl="<<tempCol<<"\n";
    }
    tissueGrid[xLoc][yLoc][zLoc].collagen3+=tempCol;
    gridOccupancy=checkGridOccupancy(xLoc,yLoc,zLoc, tissueGrid);
    if(gridOccupancy>1){
      tissueGrid[xLoc][yLoc][zLoc].collagen3-=(gridOccupancy-1);
    }
  } else {
    maxCol=tissueGrid[xLoc][yLoc][zLoc].collagen3;
    if(maxCol>0){
      colDep=min(maxCol,fibroCol1Deposition);
      tissueGrid[xLoc][yLoc][zLoc].collagen3-=colDep;
      if(tissueGrid[xLoc][yLoc][zLoc].collagen3<0){
        tissueGrid[xLoc][yLoc][zLoc].collagen3=0;
      }
      tissueGrid[xLoc][yLoc][zLoc].collagen1+=colDep;
//      cout<<"Regualr Fill\n";
    }

    gridOccupancy=checkGridOccupancy(xLoc,yLoc,zLoc, tissueGrid);
    if(gridOccupancy>1){
      tissueGrid[xLoc][yLoc][zLoc].collagen1-=(gridOccupancy-1);
    }
  }
}

void fibroSource::proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
vector<fibroblast> &fibroblasts,
mt19937 &generator
){
  int pop,i,n;
  float temp,ratio,iratio,signal,limit,baselimit,signal1,signal2;

  if(tissueGrid[xLoc][yLoc][zLoc].collagen1>0.5){
    activation=0;
//    cout<<"Activation Halted by Fibrosis\n";
    return;
  }

  baselimit=0.03;
  pop=tissueGrid[xLoc][yLoc][zLoc].cellPopulation;
  signal1=cytoGrids.TGFb[xLoc][yLoc][zLoc];
  signal2=cytoGrids.MCP1[xLoc][yLoc][zLoc];
  if(signal1>0 ||  signal2>0){
    activation++;
    // if(activation>50){
    //   activation=100;
    // }
  } else {
    //activation--;
    activation-=2;
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
    limit=baselimit*(0.5*(erf(ratio)+1));
    if(temp<limit){
      fibroblasts.push_back(fibroblast(xLoc,yLoc,zLoc));
      tissueGrid[xLoc][yLoc][zLoc].cellPopulation++;
    }
  }
}

void fibroblast::mod_sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
mt19937 &generator
){
  int i,j,ii,jj,x,y,z,flag,m;
  vector<int> indexes;
  bool val;
  int coords[9][2];

  flag=0;
  m=0;
  x=xLoc;
  y=yLoc;
  z=zLoc;

  for(i=x-1;i<=x+1;i++){
    for(j=y-1;j<=y+1;j++){
      ii=i;
      jj=j;
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
    }
  }
  if(checkValidPoint(ii,jj,z+1)){
    if(tissueGrid[ii][jj][z+1].cellPopulation<cellCapacity){
      coords[m][0]=ii;
      coords[m][1]=jj;
      indexes.push_back(m);
      m=m+1;
    }
  }
  if(m>0){
    shuffle(indexes.begin(),indexes.end(),generator);
    x=coords[indexes[0]][0];
    y=coords[indexes[0]][1];
    xLoc=x;
    yLoc=y;
    zLoc=z+1;
    tissueGrid[x][y][z+1].cellPopulation++;
    tissueGrid[x][y][z].cellPopulation--;
  } else {
    //    cout<<"M=0TEST\n";
  }
  indexes.clear();
}
