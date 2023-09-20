#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

extern uniform_int_distribution<int>  distributionEdges,distributionX,distributionY,distributionZ;
extern uniform_real_distribution<float>  dist01;

extern bool checkValidPoint(int x, int y, int z);
extern float checkGridOccupancy(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]);

extern float cytokineProductionRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids,int ruleRow, int x, int y, int z, float il1r, float tnfr);

extern void matrixDebugTester(float test1, float test2, bool exact);

satellite::satellite(){}

satellite::satellite(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  age=20;
}

myoblast::myoblast(){}

myoblast::myoblast(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
  readyToFuse=0;
  age=20;
}


void myoblast::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids){
  float cytoMax,temp;
  int i,newY;

  cytoMax=0;
  newY=yLoc;
  for(i=yLoc-1;i<=yLoc+1;i++){
    if(checkValidPoint(xLoc,i,zLoc)){
      temp=cytoGrids.IL4[xLoc][i][zLoc];
      if(temp>cytoMax){
        cytoMax=temp;
        newY=i;
      }
    }
  }
  if(newY!=yLoc){
    if(tissueGrid[xLoc][newY][zLoc].cellPopulation<cellCapacity){
      tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
      tissueGrid[xLoc][newY][zLoc].cellPopulation++;
      yLoc=newY;
    }
  }
}

void myoblast::sniff2(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids){
  int x,y,z,i,j,k,m,nextLocIndex,ii,jj,kk,temp,pop;
  float chemotaxis[27];
  float cytoMax=-1.0;
  int coords[27][3];
  x=xLoc;
  y=yLoc;
  z=zLoc;
  // for(i=0;i<27;i++){
  //   chemotaxis[i]=-1.0;
  // }
  nextLocIndex=-1;
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
            chemotaxis[m]=(cytoGrids.IL4[ii][jj][kk]);
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
//  cout<<"NLI="<<nextLocIndex<<"\n";
  if(cytoMax>0 && nextLocIndex>=0){
    x=coords[nextLocIndex][0];
    y=coords[nextLocIndex][1];
    z=coords[nextLocIndex][2];
    age++;
  }
  // updating the cell population grid with this pMN's new location
  if(tissueGrid[x][y][z].cellPopulation<cellCapacity&& nextLocIndex>=0){
    //	if(tissueGrid[x][y][z].cellPopulation<cellCapacity&&tissueGrid[x][y][z].vascularization>=0.3){
    tissueGrid[x][y][z].cellPopulation++;
    tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
    xLoc=x;
    yLoc=y;
    zLoc=z;
  } else {
      // cout<<"NLI="<<nextLocIndex<<" "<<xLoc<<" "<<yLoc<<" "<<zLoc<<"\n";
      // cout<<cytoGrids.IL4[xLoc][yLoc][zLoc]<<" "<<tissueGrid[xLoc][yLoc][zLoc].viableMuscle<<" "<<tissueGrid[xLoc][yLoc][zLoc].necrosis<<"\n";
      // cout<<tissueGrid[xLoc][yLoc][zLoc].muscleEdge<<" "<<tissueGrid[xLoc][yLoc][zLoc].woundEdge<<"\n";
//    sniff2(tissueGrid, cytoGrids);
  }
}

void myoblast::action(int index,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
  vector<myoblast> &myoblasts,
  float (&RM)[numRules][numRuleParams],mt19937 &generator
){
  int x,y,z,i,j,k,edgeCheck,ii,jj,kk;
  float gridOccupancy, myoblast_fusion_addition;
  x=xLoc;
  y=yLoc;
  z=zLoc;
  cytoGrids.MRF4[x][y][z]=cytokineProductionRule(RM,cytoGrids,47,x,y,z,0,0);
  //cytoGrids.MRF4[x][y][z]+=cytoGrids.Myf5[x][y][z];
//  cout<<"MRF4TEST="<<cytoGrids.MRF4[x][y][z]<<"\n";
  //cytoGrids.Myf5[x][y][z]+=cytoGrids.HGF[x][y][z];
  cytoGrids.Myf5[x][y][z]=cytokineProductionRule(RM,cytoGrids,46,x,y,z,0,0);


  age--;
  health-=cytoGrids.cytotox[x][y][z];
  if(age<0){
    tissueGrid[x][y][z].cellPopulation--;
//    cout<<"MYO DEATH\n";
    myoblasts.erase(myoblasts.begin()+index);
    return;
  }
  sniff2(tissueGrid, cytoGrids);
  myoblast_fusion_addition=k_myoblast_fusion_addition*(min(float(1.0),cytoGrids.HGF[x][y][z])-cytoGrids.cytotox[x][y][z]); //(-NFkB)?
 //myoblast_fusion_addition=k_myoblast_fusion_addition;
////////////////////////////TEST////////////////////////
//  myoblast_fusion_addition=0.5;
//////////////////TEST//////////////////////////////
  // cout<<"MFA TEST "<<k_myoblast_fusion_addition<<" "<<cytoGrids.HGF[x][y][z]<<" "<<(min(float(1.0),cytoGrids.HGF[x][y][z]))<<"\n";
  // cout<<"MFA BASE="<<myoblast_fusion_addition<<"\n";
  if(myoblast_fusion_addition<0){
    myoblast_fusion_addition=0;
  }
  edgeCheck=0;
  for(i=-1;i<=1;i+=2){
    for(j=-1;j<=1;j+=2){
      for(k=-1;k<=1;k+=2){
          ii=xLoc+i;
          jj=yLoc+j;
          kk=zLoc+k;
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
            kk=0;
          }
          if(kk>=zDim){
            kk=zDim-1;
          }
          if((tissueGrid[ii][jj][kk].viableMuscle>0) && (checkGridOccupancy(x,y,z, tissueGrid)<1)){
            edgeCheck=1;
          }
        }
    }
  }
//  if(tissueGrid[x][y][z].muscleEdge>0){
  if(tissueGrid[x][y][z].muscleEdge>0 || edgeCheck>0){
    readyToFuse=1;
//    cout<<" Ready To Fuse "<<x<<" "<<y<<" "<<z<<"\n";
  } else {
//    cout<<"EDGE TEST "<<tissueGrid[x][y][z].muscleEdge<<" "<<x<<" "<<y<<" "<<z<<" "<<checkGridOccupancy(x,y,z-1, tissueGrid)<<" "<<checkGridOccupancy(x,y,z, tissueGrid)<<"\n";
  }
  // if(edgeCheck==1){
  //   readyToFuse=1;
  // }
  x=xLoc;
  y=yLoc;
  z=zLoc;
  gridOccupancy=checkGridOccupancy(x,y,z, tissueGrid);
  //If there is room, myoblast fuses to end of muscle fiber
  if(readyToFuse==1){
    if(gridOccupancy<1){
//      cout<<"FUSING\n";
      tissueGrid[x][y][z].viableMuscle+=myoblast_fusion_addition;
      gridOccupancy=checkGridOccupancy(x,y,z, tissueGrid);
      if(gridOccupancy>1){
        tissueGrid[x][y][z].viableMuscle-=(gridOccupancy-1);
      }
  //    cout<<"Removing FUSED MYO "<<x<<" "<<y<<" "<<z<<"\n";
      myoblasts.erase(myoblasts.begin()+index);  //myoblast fuses into muscle fiber
      return;
    } else {
//      cout<<"Fusing ADJACENT\n";
      fuseAdjacent(tissueGrid,myoblast_fusion_addition,cytoGrids,generator);
      myoblasts.erase(myoblasts.begin()+index);
      return;
    }
  }
  return;
}

void myoblast::fuseAdjacent(gridCell (&tissueGrid)[xDim][yDim][zDim], float mfa, CytokineGrids &cytoGrids,mt19937 &generator){
  vector<int> indexes;
  int x,y,z,q,i,j,k,newX,newY,newZ;
  int candidates[2600][3];
  float gridOccupancy;
  float tissueMin=1.0;
  float tempCol,maxCol,colDep;
  int indexMin=100000;
  int flag=0;

  x=xLoc;
  y=yLoc;
  z=zLoc;
  q=0;
  //Finding the neighbors that need to be healed
  for(i=xLoc-1;i<=xLoc+1;i++){
//        for(j=yLoc-1;j<=yLoc+1;j++){
    for(j=yLoc-1;j<=yLoc+5;j++){
  //    cout<<"j="<<j<<"\n";
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
//              cout<<"GO "<<gridOccupancy<<" "<<i<<" "<<j<<" "<<k<<"\n";
              indexes.push_back(q);
              if(k<indexMin){
                indexMin=q;
              }
              q=q+1;
            }
          }
//          q=q+1;
        }
      }
    }
  }
  if(q==0){
//    cout<<"Q=0\n";
  }
  if(indexes.size()>0){
    newX=candidates[indexMin][0];
    newY=candidates[indexMin][1];
    newZ=candidates[indexMin][2];
//    cout<<newX<<" "<<newY<<" "<<newZ<<" "<<checkGridOccupancy(newX,newY,newZ, tissueGrid)<<"\n";
//    cout<<"MFA="<<mfa<<"\n";
    tissueGrid[newX][newY][newZ].viableMuscle+=mfa;
    gridOccupancy=checkGridOccupancy(newX,newY,newZ, tissueGrid);
    if(gridOccupancy>1){
      tissueGrid[newX][newY][newZ].viableMuscle-=(gridOccupancy-1);
    }
//    cout<<"Post Deposition"<<newX<<" "<<newY<<" "<<newZ<<" "<<checkGridOccupancy(newX,newY,newZ, tissueGrid)<<"\n";
  }
  //extending the IL4 gradient so more cells can fuse////////////////////////
  for(i=xLoc-1;i<=xLoc+1;i++){
    for(j=yLoc-1;j<=yLoc+1;j++){
      if(i==x&&j==y){
//        if(i==x&&j==y&&k==z){
          continue;
        } else {
          if(checkValidPoint(i,j,k)){
            gridOccupancy=checkGridOccupancy(i,j,k, tissueGrid);
            if(gridOccupancy<1){
              cytoGrids.IL4[i][j][k]=cytoGrids.IL4[xLoc][yLoc][zLoc]+dist01(generator);
              cytoGrids.IL4[xLoc][yLoc][zLoc]=0;
            }
          }
        }
      }
    }
  return;
}

void satellite::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids){
  float cytoMax,temp;
  int i,newY;

  cytoMax=0;
  newY=yLoc;
  for(i=yLoc-1;i<=yLoc+1;i++){
    if(checkValidPoint(xLoc,i,zLoc)){
      temp=cytoGrids.TNF[xLoc][i][zLoc]+cytoGrids.IL6[xLoc][i][zLoc];
      if(temp>cytoMax){
        cytoMax=temp;
        newY=i;
      }
    }
  }

  if(newY!=yLoc){
    if(tissueGrid[xLoc][newY][zLoc].cellPopulation<cellCapacity){
      tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
      tissueGrid[xLoc][newY][zLoc].cellPopulation++;
      yLoc=newY;
    }
  }
}

void satellite::action(int index,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
vector<satellite> &satellites,
vector<myoblast> &myoblasts){
  int x,y,z;
  x=xLoc;
  y=yLoc;
  z=zLoc;
  //cout<<"Satellite Action\n";
  sniff(tissueGrid, cytoGrids);
  //cout<<"Position="<<" "<<xLoc<<" "<<yLoc<<" "<<zLoc<<" "<<cytoGrids.TNF[x][y][z]<<" "<<cytoGrids.IL6[x][y][z]<<" "<<cytoGrids.IL10[x][y][z]<<" "<<cytoGrids.IL4[x][y][z]<<"\n";
  if(cytoGrids.TNF[x][y][z]+cytoGrids.IL6[x][y][z]<cytoGrids.IL10[x][y][z]+cytoGrids.IL4[x][y][z]){
    satellites.erase(satellites.begin()+index);
  //  cout<<"Differentiating into Myoblast\n";
    myoblasts.push_back(myoblast(x,y,z));
    return;
  }
  age--;
  if(age<=0){
    tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
    satellites.erase(satellites.begin()+index);
  }
  return;
}

void satelliteSource::proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
vector<satellite> &satellites,vector<myoblast> &myoblasts,
mt19937 &generator,CytokineGrids &cytoGrids,
float &TotalTNF){
  int pop,i,n;
  float temp,limit,active,gridOccupancy;
//  cout<<zLoc<<" "<<cytoGrids.IL4[xLoc][yLoc][zLoc]<<"\n";
  pop=tissueGrid[xLoc][yLoc][zLoc].cellPopulation;
  if(zLoc==0){
    activation=0;
  }
  if(cytoGrids.IL4[xLoc][yLoc][zLoc]>0) {
   // if(activation==0){cout<<"XXXXX ACTIVATING "<<zLoc<<" "<<cytoGrids.IL4[xLoc][yLoc][zLoc]<<" "<<cytoGrids.MCP1[xLoc][yLoc][zLoc]<<"\n";
   //   cout<<tissueGrid[xLoc][yLoc][zLoc].viableMuscle<<" "<<tissueGrid[xLoc][yLoc][zLoc].necrosis<<"\n";
   // }
    activation++;
  } else {
    if(activation>0){
//    cout<<"DEACTIVATING "<<zLoc<<"\n";
}
    activation-=2;
  }   // if(activation==0){cout<<"XXXXX ACTIVATING "<<zLoc<<" "<<cytoGrids.IL4[xLoc][yLoc][zLoc]<<" "<<cytoGrids.MCP1[xLoc][yLoc][zLoc]<<"\n";
   //   cout<<tissueGrid[xLoc][yLoc][zLoc].viableMuscle<<" "<<tissueGrid[xLoc][yLoc][zLoc].necrosis<<"\n";
   // }
  if(activation<0){
    activation=0;
  }
  if(activation>10){
    activation=10;
    tissueGrid[xLoc][yLoc][zLoc].viableMuscle+=0.5;
    gridOccupancy=checkGridOccupancy(xLoc,yLoc,zLoc,tissueGrid);
    if(gridOccupancy>1){
      if(tissueGrid[xLoc][yLoc][zLoc].collagen3>0){
        tissueGrid[xLoc][yLoc][zLoc].collagen3=0;
      } else {
        tissueGrid[xLoc][yLoc][zLoc].viableMuscle-=(gridOccupancy-1);
      }
    }
    if(pop<cellCapacity){
        limit=0.1;
        temp=dist01(generator);
        if(temp<limit){
          myoblasts.push_back(myoblast(xLoc,yLoc,zLoc));
          tissueGrid[xLoc][yLoc][zLoc].cellPopulation++;
        }
    }
  }
}
