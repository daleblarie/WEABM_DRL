#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "WEABM_agents.h"
#include "WEABM_Parameters.h"

extern uniform_real_distribution<float>  dist01;
extern uniform_int_distribution<int>  distribution3,distribution27;

extern uniform_int_distribution<int>  distributionEdges,distributionX,distributionY,distributionZ;

extern float cytokineProductionRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr);
extern float cytokineComboRule(float (&RM)[numRules][numRuleParams], CytokineGrids &cytoGrids, int ruleRow, int x, int y, int z, float il1r, float tnfr);
extern void matrixDebugTester(float test1, float test2, bool exact);

int getWiggleCoord(int x, int dimMax, mt19937 &generator);

Th0::Th0(){}

Th0::Th0(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
	age=50;
}

Th1::Th1(){}

Th1::Th1(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
	age=50;
}

Th2::Th2(){}

Th2::Th2(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
	age=50;
}

Th17::Th17(){}

Th17::Th17(int x, int y, int z){
  xLoc=x;
  yLoc=y;
  zLoc=z;
	age=50;
}

void Th0::action(int index,
gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
vector<Th0> &Th0s,
vector<Th1> &Th1s,
vector<Th2> &Th2s,
vector<Th17> &Th17s,
mt19937 &generator,
float (&RM)[numRules][numRuleParams]){

	int x,y,z;
	float normalizer,proTh1,proTh2,proTh17,temp;

  sniff(generator);
  x=xLoc;
  y=yLoc;
  z=zLoc;
  proTh1=cytokineProductionRule(RM,cytoGrids,15,x,y,z,0,0);
  proTh2=cytokineProductionRule(RM,cytoGrids,16,x,y,z,0,0);
  proTh17=cytokineProductionRule(RM,cytoGrids,17,x,y,z,0,0);
	normalizer=proTh1+proTh2+proTh17;
	proTh1/=normalizer;
	proTh2/=normalizer;
	proTh17/=normalizer;

	age--;

	if(age<0){
		tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
		Th0s.erase(Th0s.begin()+index);
		return;
  }

	temp=dist01(generator);

	if(temp<proTh1){
		Th1s.push_back(Th1(xLoc,yLoc,zLoc));
		Th0s.erase(Th0s.begin()+index);
    return;
	}
	if(temp>proTh1&&temp<proTh1+proTh2){
		Th2s.push_back(Th2(xLoc,yLoc,zLoc));
		Th0s.erase(Th0s.begin()+index);
    return;
	}
	if(temp>(proTh1+proTh2)){
		Th17s.push_back(Th17(xLoc,yLoc,zLoc));
		Th0s.erase(Th0s.begin()+index);
    return;
	}
}

void Th1::action(int index,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
  vector<Th1> &Th1s,
  mt19937 &generator,
  float (&RM)[numRules][numRuleParams]
){
	int id,x,y,z;

  sniff(tissueGrid, cytoGrids, generator);
  x=xLoc;
	y=yLoc;
	z=zLoc;
	id=yLoc*xDim+xLoc;
  //activation=cytokineProductionRule(RM,cytoGrids,18,x,y,z,0,0);
  //if(activation>0){
	if(cytoGrids.IL12[x][y][z]>0){
    cytoGrids.IFNg[x][y][z]=cytokineProductionRule(RM,cytoGrids,19,x,y,z,0,0);
	}
	age--;
	if(age<0){
		tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
		Th1s.erase(Th1s.begin()+index);
    return;
	}
}

void Th2::action(int index,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
  vector<Th2> &Th2s,
  mt19937 &generator,
  float (&RM)[numRules][numRuleParams]){
  int id,x,y,z;

  sniff(tissueGrid, cytoGrids, generator);
  x=xLoc;
	y=yLoc;
	z=zLoc;

  //activation=cytokineProductionRule(20,x,y,z,0,0);
  //if(activation>0){
	if(cytoGrids.IL10[x][y][z]>0){
    cytoGrids.IL4[x][y][z]=cytokineProductionRule(RM,cytoGrids,21,x,y,z,0,0);
    cytoGrids.IL10[x][y][z]=cytokineProductionRule(RM,cytoGrids,22,x,y,z,0,0);
	}
	age--;
	if(age<0){
		tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
		Th2s.erase(Th2s.begin()+index);
    return;
	}
}

void Th17::action(int index,
  gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids,
  vector<Th17> &Th17s,
  mt19937 &generator,
  float (&RM)[numRules][numRuleParams]){
  int id,x,y,z;

  sniff(tissueGrid, cytoGrids, generator);
  x=xLoc;
	y=yLoc;
	z=zLoc;

	id=yLoc*xDim+xLoc;

  cytoGrids.IL17[x][y][z]=cytokineProductionRule(RM,cytoGrids,23,x,y,z,0,0);
	age--;
	if(age<0){
		tissueGrid[xLoc][yLoc][zLoc].cellPopulation--;
		Th17s.erase(Th17s.begin()+index);
    return;
	}
}

void Th17::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
mt19937 &generator
){
	int x,y,z,i,j,k,m,nextLocIndex,ii,jj,kk,temp;
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
          ii=yDim-1;
        }
        if(jj>=yDim){
          jj=0;
        }
        if(kk<0){
          chemotaxis[m]=-1;
        } else if(kk>=zDim){
          chemotaxis[m]=-1;
        } else if(kk>=0&&kk<=zDim){
          if(tissueGrid[ii][jj][kk].cellPopulation>=cellCapacity){
            chemotaxis[m]=-1;
          }
        } else{
          chemotaxis[m]=cytoGrids.TGFb[ii][jj][kk]+cytoGrids.IL1[ii][jj][kk];
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

void Th1::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
mt19937 &generator){
	int x,y,z,i,j,k,m,nextLocIndex,ii,jj,kk,temp;
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
          ii=yDim-1;
        }
        if(jj>=yDim){
          jj=0;
        }
        if(kk<0){
          chemotaxis[m]=-1;
        } else if(kk>=zDim){
          chemotaxis[m]=-1;
        } else if(kk>=0&&kk<=zDim){
          if(tissueGrid[ii][jj][kk].cellPopulation>=cellCapacity){
            chemotaxis[m]=-1;
          }
        } else{
          chemotaxis[m]=cytoGrids.IL4[ii][jj][kk]+cytoGrids.IL10[ii][jj][kk];
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
//  cout<<"NLI="<<nextLocIndex<<"\n";
}

void Th2::sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
CytokineGrids &cytoGrids,
mt19937 &generator){
  int x,y,z,i,j,k,m,nextLocIndex,ii,jj,kk,temp;
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
          ii=yDim-1;
        }
        if(jj>=yDim){
          jj=0;
        }
        if(kk<0){
          chemotaxis[m]=-1;
        } else if(kk>=zDim){
          chemotaxis[m]=-1;
        } else if(kk>=0&&kk<=zDim){
          if(tissueGrid[ii][jj][kk].cellPopulation>=cellCapacity){
            chemotaxis[m]=-1;
          }
        } else{
          chemotaxis[m]=cytoGrids.TGFb[ii][jj][kk]+cytoGrids.IL1[ii][jj][kk];
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

void Th0::sniff(mt19937 &generator){
  int x,y,z;
  xLoc=getWiggleCoord(xLoc,xDim, generator);
  yLoc=getWiggleCoord(yLoc,yDim, generator);
  if(zLoc>0){
    zLoc=getWiggleCoord(zLoc,zDim, generator);
  } else {
    zLoc=0;
  }
}

int getWiggleCoord(int x, int dimMax, mt19937 &generator
){
  int temp;
  temp=distribution3(generator);
  if(x==dimMax){
    cout<<"DIM ERROR "<<x<<"\n";
  }
  if(temp==0){
    x=x-1;
    if(x<0){
      x=dimMax-1;
    }
  }
  if(temp==2){
    x=x+1;
    if(x>=dimMax){
      x=0;
    }
  }
  return x;
}

void T_Source::proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
vector<Th0> &Th0s,
mt19937 &generator,
float TotalTGFb
){
  int pop,i,n;
  float temp,limit;
  limit=0.01;
  pop=tissueGrid[xLoc][yLoc][zLoc].cellPopulation;
  if(pop<cellCapacity){
    n=int(1+TotalTGFb/100);
    for(i=0;i<n;i++){
      temp=dist01(generator);
      if(temp<limit){
        Th0s.push_back(Th0(xLoc,yLoc,zLoc));
        tissueGrid[xLoc][yLoc][zLoc].cellPopulation++;
      }
    }
  }
}
