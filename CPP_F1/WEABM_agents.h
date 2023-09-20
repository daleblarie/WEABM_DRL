#ifndef WEABM_AGENTS_H
#define WEABM_AGENTS_H
#include<vector>
#include <stdlib.h>
#include <iostream>
#include "WEABM_Parameters.h"
#include <fstream>

using namespace std;

// Forward delcaration of all the cell types so that it doesnt throw errors
// when other cell types are given as parameters to cell functions
struct gridCell;
struct pmn;
struct macrophage;
struct Th0;
struct Th1;
struct Th2;
struct Th17;
struct fibroblast;
struct satellite;
struct myoblast;
struct pmnSource;
struct macroSource;
struct T_Source;
struct fibroSource;
struct satelliteSource;




struct CytokineGrids{
public:
  float IL8[xDim][yDim][zDim];
  float MIP1b[xDim][yDim][zDim];
  float TNF[xDim][yDim][zDim];
  float HGF[xDim][yDim][zDim];
  float PAF[xDim][yDim][zDim];
  float IL1[xDim][yDim][zDim];
  float PDGF[xDim][yDim][zDim];
  float IL13[xDim][yDim][zDim];
  float TGFb[xDim][yDim][zDim];
  float VEGF[xDim][yDim][zDim];
  float GCSF[xDim][yDim][zDim];
  float IL10[xDim][yDim][zDim];
  float IL4[xDim][yDim][zDim];
  float IL17[xDim][yDim][zDim];
  float Myf5[xDim][yDim][zDim];
  float IL6[xDim][yDim][zDim];
  float IL1ra[xDim][yDim][zDim];
  float sIL1r[xDim][yDim][zDim];
  float sTNFr[xDim][yDim][zDim];
  float endotoxin[xDim][yDim][zDim];
  float IFNg[xDim][yDim][zDim];
  float IL12[xDim][yDim][zDim];
  float MRF4[xDim][yDim][zDim];
  float cytotox[xDim][yDim][zDim];
  float IGF[xDim][yDim][zDim];
  float DAMP[xDim][yDim][zDim];
  float MCP1[xDim][yDim][zDim];
  float FNE[xDim][yDim][zDim];
};
struct gridCell{
  int xLoc,yLoc,zLoc,id;   //ind(i,j,k) = k + j*nz + i*ny*nz
  int cellPopulation; //number of cells that exist at that location
  float innervation;  //functional innervation in tissue voxel
  float collagen1; //type 1 (permanent) collagen (scar)
  float collagen3; //wound initially fills with weaker type 3 collagen, comprising ECM/granulation tissue
  float viableMuscle; //viable muscle tissue
  float vascularization; //vascularization of tissue; 1 is max value representing fully vascularized
  float life;
  float necrosis;
  float contamination;

  int counter; //debugging variable

  //Endothelial Cytokines

  //Endothelial Cell: Neutrophil Controls
  int ec_activation,ec_roll,ec_stick,ec_migrate;
  float tempSignal;
  bool satellitePool;
  int woundEdge,colEdge,muscleEdge,vascularEdge;

  gridCell();
  gridCell(int x, int y, int z, int index);

  void endothelialFunction(gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids, vector<pmnSource> &pmnSources, vector<macroSource> &macroSources,
  vector<T_Source> &T_Sources, vector<fibroSource> &fibroSources, vector<satelliteSource> &satelliteSources, mt19937 &generator,
  float (&RM)[numRules][numRuleParams]);
  void endothelium_activate(CytokineGrids &cytoGrids,float (&RM)[numRules][numRuleParams]);
  void injurySpread(float oxyHeal,
    gridCell (&tissueGrid)[xDim][yDim][zDim],
  CytokineGrids &cytoGrids);
  int checkForVascularEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]);
  int checkForWoundEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]);
  int checkForMuscleEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]);
  int checkForMuscleCollagenEdge(gridCell (&tissueGrid)[xDim][yDim][zDim]);
//  void injuryPropagation(int x, int y, int z, float delta);
  void proliferation(gridCell (&tissueGrid)[xDim][yDim][zDim], CytokineGrids &cytoGrids, vector<pmnSource> &pmnSources, vector<macroSource> &macroSources, vector<T_Source> &T_Sources, vector<fibroSource> &fibroSources, vector<satelliteSource> &satelliteSources,mt19937 &generator);
};

struct endoCell{
   public:
	    int xLoc;
	    int yLoc;
      int zLoc;
	    int id;

  endoCell();
  endoCell(int x, int y, int z, int index);
};

struct epiCell{
  public:
    int xLoc;
    int yLoc;
    int zLoc;
    int id;

    epiCell();
    epiCell(int x, int y, int z, int index);
};

struct pmn{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int age;
     float wbc_roll; //selectins
     float wbc_stick; //integrins
     float wbc_migrate; //diapedesis

     pmn();
     pmn(int x, int y, int z);

     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
      CytokineGrids &cytoGrids,
      mt19937 &generator);
     void pmn_function(int pmnID,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
       vector<pmn> &pmns,
       int &bursts,
       mt19937 &generator,
       float (&RM)[numRules][numRuleParams]);
     void pmn_burst(int pmnID,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
       vector<pmn> &pmns,
       int &bursts,
       float (&RM)[numRules][numRuleParams]);
};

struct macrophage{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int age;
     float wbc_stick,wbc_migrate,wbc_roll;
     float TNFr,IL1r,activation;

     macrophage();
     macrophage(int x, int y, int z);

     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
      CytokineGrids &cytoGrids,
      mt19937 &generator);
     void action(int index, gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids,
     vector<macrophage>& macrophages,
     mt19937 &generator,
     float (&RM)[numRules][numRuleParams]
);
     void heal(int x, int y, int z, gridCell (&tissueGrid)[xDim][yDim][zDim]);
};

struct Th0{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int id;
     int age;

     Th0();
     Th0(int x, int y, int z);

     void sniff(mt19937 &generator);
     void action(int index,
     gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids,
     vector<Th0> &Th0s,
     vector<Th1> &Th1s,
     vector<Th2> &Th2s,
     vector<Th17> &Th17s,
     mt19937 &generator,
     float (&RM)[numRules][numRuleParams]);
};

struct Th1{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int id;
     int age;

     Th1();
     Th1(int x, int y, int z);

     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids,
     mt19937 &generator);
     void action(int index,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
       vector<Th1> &Th1s,
       mt19937 &generator,
       float (&RM)[numRules][numRuleParams]);
};

struct Th2{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int id;
     int age;

     Th2();
     Th2(int x, int y, int z);

     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids,
     mt19937 &generator);
     void action(int index,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
       vector<Th2> &Th2s,
       mt19937 &generator,
       float (&RM)[numRules][numRuleParams]);
};

struct Th17{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int id;
     int age;

     Th17();
     Th17(int x, int y, int z);

     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids,
     mt19937 &generator);
     void action(int index,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
       vector<Th17> &Th17s,
       mt19937 &generator,
       float (&RM)[numRules][numRuleParams]);
};

struct fibroblast{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int id;
     int age;
     float activation;
     float sma;
     bool isMyofibroblast;

     void death(int index,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
     vector<fibroblast> &fibroblasts);
     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids,
     mt19937 &generator);
     void mod_sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
     mt19937 &generator);

     fibroblast();
     fibroblast(int x, int y, int z);
     void action(int index,
     gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids,
     vector<fibroblast> &fibroblasts,
     mt19937 &generator,int step,
   float (&RM)[numRules][numRuleParams]);
     int fillAdjacentF(gridCell (&tissueGrid)[xDim][yDim][zDim],CytokineGrids &cytoGrids);
     int fillAdjacentM(gridCell (&tissueGrid)[xDim][yDim][zDim],CytokineGrids &cytoGrids);
     void fillCollagen(float gridOccupancy,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids);

};

struct satellite{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int age;

     bool isMobile;

     satellite();
     satellite(int x, int y, int z);

     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids);
     void action(int index,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
       vector<satellite> &satellites,
       vector<myoblast> &myoblasts);
     void proliferation();
};

struct satelliteSource{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int activation;

     satelliteSource();
     satelliteSource(int x, int y, int z);

     void proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
     vector<satellite> &satellites, vector<myoblast> &myoblasts,
     mt19937 &generator,CytokineGrids &cytoGrids,
    float &TotalTNF);
};

struct myoblast{
   public:
     int xLoc;
     int yLoc;
     int zLoc;
     int id;
     bool fusedBack,fusedFront;
     float mCadherin;

     int age;
     float health;
     float readyToFuse;

     myoblast();
     myoblast(int x, int y, int z);

     void sniff(gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids);
     void sniff2(gridCell (&tissueGrid)[xDim][yDim][zDim],
     CytokineGrids &cytoGrids);
     void action(int index,
       gridCell (&tissueGrid)[xDim][yDim][zDim],
       CytokineGrids &cytoGrids,
       vector<myoblast> &myoblasts,
       float (&RM)[numRules][numRuleParams],mt19937 &generator);
     void proliferation();
     void fuseAdjacent(gridCell (&tissueGrid)[xDim][yDim][zDim],float mfa,CytokineGrids &cytoGrids,mt19937 &generator);
};

struct pmnSource{
  public:
    int xLoc,yLoc,zLoc,activation;
    pmnSource();
    pmnSource(int x, int y, int z);
    void proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
    CytokineGrids &cytoGrids,
    vector<pmn>& pmns,
    mt19937 &generator);
};

struct macroSource{
  public:
    int xLoc,yLoc,zLoc,activation;
    macroSource();
    macroSource(int x, int y, int z);
    void proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
    CytokineGrids &cytoGrids,
    vector<macrophage>& macrophages,
    mt19937 &generator);
};

struct T_Source{
  public:
    int xLoc,yLoc,zLoc,activation;
    T_Source();
    T_Source(int x, int y, int z);
    void proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
    vector<Th0> &Th0s,
    mt19937 &generator,
    float TotalTGFb);
};

struct fibroSource{
  public:
    int xLoc,yLoc,zLoc,activation;
    fibroSource();
    fibroSource(int x, int y, int z);
    void proliferate(gridCell (&tissueGrid)[xDim][yDim][zDim],
    CytokineGrids &cytoGrids,
    vector<fibroblast> &fibroblasts,
    mt19937 &generator);
};

struct SimulationObject{
public:
  mt19937 generator;

  float RM[numRules][numRuleParams];
  float internalParameterization[1260];

  int xDimension, yDimension;

  CytokineGrids cytoGrids;

  vector<pmn> pmns;
  vector<macrophage> macrophages;
  vector<Th0> Th0s;
  vector<Th1> Th1s;
  vector<Th2> Th2s;
  vector<Th17> Th17s;
  vector<fibroblast> fibroblasts;
  vector<satellite> satellites;
  vector<myoblast> myoblasts;
  vector<int> gridIndexes;
  vector<pmnSource> pmnSources;
  vector<macroSource> macroSources;
  vector<T_Source> T_Sources;
  vector<fibroSource> fibroSources;
  vector<satelliteSource> satelliteSources;

  ofstream outputFile,tissueFile;

  gridCell tissueGrid[xDim][yDim][zDim];


  float TotalIL8,TotalMIP1b,TotalTNF,TotalPAF,TotalIL1,TotalPDGF,TotalIL13,TotalTGFb,
  TotalVEGF,TotalGCSF,TotalIL10,TotalIL4,TotalIL17,TotalMyf5,TotalIL6,TotalIL1ra,
  TotalsIL1r,TotalsTNFr,TotalEndotoxin,TotalIFNg,TotalDAMP,TotalNecrosis,TotalC1,
  TotalC3,TotalViabMus,TotalMCP1,TotalHGF,TotalIL12,TotalMRF4,TotalIGF,TotalIL1r,TotalTNFr;

  int bursts, current_step;

  float allSignals[43][numTimeSteps];
  float allSignalsReturn[43*numTimeSteps];

  void setupSimulation(int rseed);
  void singleStep();
  float* endSimulation();

  void getRuleMatrix(float internalParam[], int numMatEls);
  void initialize();
  void initializeFromCT();
  void injure_biopsy();
  void setupBiopsySatellites();

  void cellFunctions();
  void evaporations();
  void evaporate(float (&cytokineGrid)[xDim][yDim][zDim], float evaporationConstant);
  void diffusions();
  void diffuse(float (&cytokineGrid)[xDim][yDim][zDim], float diffusionConstant);
  void aggregateData(int step);


  int checkForInitialNeighboringTissueVoid(int x, int y, int z);

  void outputSetup();
  void outputFinalize();
  void outputWriter(float cytokine[xDim][yDim][zDim], int step, int fixedDim);
  void binaryTissueOutput();
  void aggregateOutput();
  void gridOutput(int step);

  // Cytokine Grids for ouput
  int topGrid[xDim][yDim];
  float surfaceCytokines[28];
  float totalCytokines[28];
  float temp_out[xDim*yDim*zDim];
  float innervation_out[xDim*yDim*zDim];
  float collagen1_out[xDim*yDim*zDim];
  float collagen3_out[xDim*yDim*zDim];
  float viableMuscle_out[xDim*yDim*zDim];
  float vascularization_out[xDim*yDim*zDim];

  float life_out[xDim*yDim*zDim];
  float necrosis_out[xDim*yDim*zDim];
  float contamination_out[xDim*yDim*zDim];


  // Functions used for ctypes
    void getTopGrid();

    void updateSurfaceCells();
    float* getSurfaceCytokines();
    float* getTotalCytokines();
    void updateOutputGrids();
    float* getInnervation();
    float* getCollagen1();
    float* getCollagen3();
    float* getViableMuscle();
    float* getVascularization();
    float* getLife();
    float* getNecrosis();
    float* getContamination();
    float* get_IL8();
    float* get_MIP1b();
    float* get_TNF();
    float* get_HGF();
    float* get_PAF();
    float* get_IL1();
    float* get_PDGF();
    float* get_IL13();
    float* get_TGFb();
    float* get_VEGF();
    float* get_GCSF();
    float* get_IL10();
    float* get_IL4();
    float* get_IL17();
    float* get_Myf5();
    float* get_IL6();
    float* get_IL1ra();
    float* get_sIL1r();
    float* get_sTNFr();
    float* get_endotoxin();
    float* get_IFNg();
    float* get_IL12();
    float* get_MRF4();
    float* get_cytotox();
    float* get_IGF();
    float* get_DAMP();
    float* get_MCP1();
    float* get_FNE();

    void apply_IL8(float addedAmount);
    void apply_MIP1b(float addedAmount);
    void apply_TNF(float addedAmount);
    void apply_HGF(float addedAmount);
    void apply_PAF(float addedAmount);
    void apply_IL1(float addedAmount);
    void apply_PDGF(float addedAmount);
    void apply_IL13(float addedAmount);
    void apply_TGFb(float addedAmount);
    void apply_VEGF(float addedAmount);
    void apply_GCSF(float addedAmount);
    void apply_IL10(float addedAmount);
    void apply_IL4(float addedAmount);
    void apply_IL17(float addedAmount);
    void apply_Myf5(float addedAmount);
    void apply_IL6(float addedAmount);
    void apply_IL1ra(float addedAmount);
    void apply_sIL1r(float addedAmount);
    void apply_sTNFr(float addedAmount);
    void apply_endotoxin(float addedAmount);
    void apply_IFNg(float addedAmount);
    void apply_IL12(float addedAmount);
    void apply_MRF4(float addedAmount);
    void apply_cytotox(float addedAmount);
    void apply_IGF(float addedAmount);
    void apply_DAMP(float addedAmount);
    void apply_MCP1(float addedAmount);
    void apply_FNE(float addedAmount);

    int seed;
    void setSeed(int newSeed){seed = newSeed; generator.seed(seed);}

  SimulationObject();
  SimulationObject(float internalParam[]);
  float* sim_obj_main();
};




#endif
