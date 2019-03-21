#ifndef HiggsToMuMuProcessor_h
#define HiggsToMuMuProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include <string>

using namespace lcio;
using namespace marlin;

class TFile;
class TTree;

const int NMAX_PFOS     = 500;
const int NMAX_TRACKS   = 500;
const int NMAX_CLUSTERS = 500;

namespace EVENT{
  class ReconstructedParticle;
  class Track;
  class Cluster;
}

class HiggsToMuMuProcessor : public Processor{
  
 public:
  
  virtual Processor* newProcessor() { return new HiggsToMuMuProcessor; }
  
  HiggsToMuMuProcessor();
  
  //Called at the begin of the job before anything is read.
  //Use to initialize the processor, e.g. book histograms.
  virtual void init();
  
  //Called for every run.
  virtual void processRunHeader( LCRunHeader* run );
  
  //Called for every event - the working horse.
  virtual void processEvent( LCEvent* evt ); 
  
  virtual void check( LCEvent* evt ); 

  //Called after data processing for clean up.
  virtual void end();
  
 protected:

  //Prepare NTuple to check the data
  void makeNTuple();
 
  //get energy deposit on calorimeter
  //[0]:Ecal energy, [1]:Hcal energy
  void getCalEnergy( ReconstructedParticle* pfo , float* cale );

  //Read process list and add to map
  void readProcessList(const char* fname);
 
  //get tau decay mode information
  int getTauDecayMode(const MCParticle *mcp);

  MCParticle* findMCParticleFromIsolep( ReconstructedParticle* rp, LCCollection* AllPFOs );

  //Process map
  typedef std::multimap < std::string,
                          std::pair < int, float > > TMapProcess;
  TMapProcess _mapProcess;
  bool _useProcessMap;

  //Input collection name.
  std::string _colAllPFOs;
  std::string _colIsoleps;
  std::string _colPFOsWithoutIsoleps;
  std::string _colISR;
  std::string _colMC;
  std::string _mcpfoRelation;
  LCCollection* _mcpCol;
  LCRelationNavigator* _navpfo;
  // output root file name
  std::string _rootfilename;
  // process list file
  std::string _processList;

  float _energetic;
  int _nRun;
  int _nEvt;

  //Stored data 
  struct RECO_DATA{
    //***********
    //event number
    //***********
    int evt;

    //************************
    //inserted by process file
    //************************
    int   processid;
    float xsec;

    //**************
    //number of PFOs
    //**************
    int npfos;

    //******************************************
    //read out the thrust and sphere information
    //******************************************
    float principalthrust;
    float majorthrust;
    float minorthrust;
    float costh_thrustaxis;
    //float thrustAxisX;
    //float thrustAxisY;
    //float thrustAxisZ;
    float sphericity;
    float aplanarity;

    //************
    //MC loop part
    //************
    float ISR1_E;
    float ISR2_E;
    float ISR1_costh;
    float ISR2_costh;

    int   type;
    int   tau_m;
    int   higgsdecay1;
    int   higgsdecay2;
    float MC_mumu_mass;
    float MC_mumu_costh;

    //counter for h->tautau and tau->mu
    int h_nt_m;

    //counter for h->WW*,ZZ* and their decay
    int h_nq;
    int h_ne;
    int h_nm;
    int h_nt;
    int h_nn;
    int h_nboson_m;

    //2f process determination
    int   forff1;
    int   forff2;
    float ff1_costh;
    float ff2_costh;
    float massff;
    int   ff_nt_m;

    //4f process determination
    int   forffff1;
    int   forffff2;
    int   forffff3;
    int   forffff4;
    int   nq;
    int   ne;
    int   nm;
    int   nt;
    int   nn;
    int   no;
    int   nt_m;
    float ffff1_costh;
    float ffff1_t_costh;
    float ffff2_costh;
    float ffff2_t_costh;
    float ffff3_costh;
    float ffff3_t_costh;
    float ffff4_costh;
    float ffff4_t_costh;
    float mass1;
    float mass2;

    int   MC_muminus_parent;
    int   MC_muplus_parent;
    float MC_mumupair_mass;

    float MC_mu_E;
    //float MC_Z_mumu_mass;

    int ntaumode;
    int taumode[NMAX_PFOS];

    //*********************************************
    //check the relation between MCParticle and PFO
    //*********************************************
    float esum_mc;
    float esum_charged_mc;
    float esum_neutral_mc;
    float overlay_charged_mc;
    float overlay_neutral_mc;
    float overlay_E_mc;
    float costh_missmom_mc;
    float Pt_all_mc;

    int   n_highPt_mc;
    float Pt_mc[NMAX_PFOS];
    float costh_mc[NMAX_PFOS];
    float pfoPt[NMAX_PFOS];
    float pfocosth[NMAX_PFOS];

    float mumu_mass_mc;
    float mumu_costh_mc;
    int   n_muminus_mc;
    int   n_muplus_mc;

    float mumu_mass_mc_smear_1_3;
    float mumu_mass_mc_smear_5_4;
    float mumu_mass_mc_smear_3_4;
    float mumu_mass_mc_smear_2_4;
    float mumu_mass_mc_smear_1_4;
    float mumu_mass_mc_smear_5_5;
    float mumu_mass_mc_smear_3_5;
    float mumu_mass_mc_smear_2_5;
    float mumu_mass_mc_smear_1_5;
    float mumu_mass_mc_smear_5_6;
    float mumu_mass_mc_smear_3_6;
    float mumu_mass_mc_smear_2_6;
    float mumu_mass_mc_smear_1_6;

    //*************
    //PFO loop part
    //*************
    int   ntrks;
    float esum;
    float psum_x;
    float psum_y;
    float psum_z;
    float maxPt;
    float Pt;
    float Pt_all;
    float mass_vis;
    float misse;
    float missmom_x;
    float missmom_y;
    float missmom_z;
    float missmom_mag;
    float costh_missmom;
    float mass_inv;

    float pfo_e[NMAX_PFOS];
    float pfo_px[NMAX_PFOS];
    float pfo_py[NMAX_PFOS];
    float pfo_pz[NMAX_PFOS];
    float pfo_mom_mag[NMAX_PFOS];
    float pfo_chrg[NMAX_PFOS];
    int   pfo_pdg[NMAX_PFOS];
    float pfo_costh[NMAX_PFOS];
    float pfo_Pt[NMAX_PFOS];
    int   pfo_ntrk[NMAX_PFOS];
    int   pfo_nclus[NMAX_PFOS];

    float pfo_d0[NMAX_PFOS];
    float pfo_z0[NMAX_PFOS];
    float pfo_phi[NMAX_PFOS];
    float pfo_omega[NMAX_PFOS];
    float pfo_tanlambda[NMAX_PFOS];
    float pfo_d0sig[NMAX_PFOS];
    float pfo_z0sig[NMAX_PFOS];
    float pfo_Chi2[NMAX_PFOS];
    int   pfo_Ndf[NMAX_PFOS];
    float pfo_Chi2Ndf[NMAX_PFOS];
    float pfo_dEdx[NMAX_PFOS];
    float pfo_dEdxError[NMAX_PFOS];
    float pfo_RadiusOfInnermostHit[NMAX_PFOS];
    int   pfo_VXD[NMAX_PFOS];
    int   pfo_FTD[NMAX_PFOS];
    int   pfo_SIT[NMAX_PFOS];
    int   pfo_TPC[NMAX_PFOS];
    int   pfo_SET[NMAX_PFOS];
    int   pfo_ETD[NMAX_PFOS];//always 0

    float pfo_ecal[NMAX_PFOS];
    float pfo_hcal[NMAX_PFOS];
    float pfo_yoke[NMAX_PFOS];
    float pfo_lcal[NMAX_PFOS];
    float pfo_lhcal[NMAX_PFOS];
    float pfo_bcal[NMAX_PFOS];
    float pfo_ecalfrac[NMAX_PFOS];
    float pfo_EbyP[NMAX_PFOS];

    int   pfo_pid_mc[NMAX_PFOS];
    float pfo_E_mc[NMAX_PFOS];
    float pfo_Px_mc[NMAX_PFOS];
    float pfo_Py_mc[NMAX_PFOS];
    float pfo_Pz_mc[NMAX_PFOS];
    //int   pfo_taudecaymode[NMAX_PFOS];
    //int   count_overlay;
    //int   charged_overlay;

    int   parent_pid_mc[NMAX_PFOS];
    float parent_E_mc[NMAX_PFOS];
    float parent_Px_mc[NMAX_PFOS];
    float parent_Py_mc[NMAX_PFOS];
    float parent_Pz_mc[NMAX_PFOS];
    int   pfo_mother_pid[NMAX_PFOS];

    //*************************
    //isolated lepton loop part
    //*************************
    int   n_isoleps;
    int   n_muplus;
    float muplus_E;
    float muplus_costh;
    float muplus_charge_costh;
    float muplus_azi;
    float muplus_Pt;
    float muplus_mom_mag;
    int   n_muminus;
    float muminus_E;
    float muminus_costh;
    float muminus_charge_costh;
    float muminus_azi;
    float muminus_Pt;
    float muminus_mom_mag;
    float sum_charge_costh;
    float leadingmu_E;
    float subleadingmu_E;
    float leadingmu_costh;
    float subleadingmu_costh;

    float muplus_d0;
    float muplus_z0;
    float muplus_phi;
    float muplus_omega;
    float muplus_tanlambda;
    float muplus_d0sig;
    float muplus_z0sig;
    float muplus_Chi2;
    int   muplus_Ndf;
    float muplus_Chi2Ndf;
    float muplus_dEdx;
    float muplus_dEdxError;
    float muplus_RadiusOfInnermostHit;
    float muminus_d0;
    float muminus_z0;
    float muminus_phi;
    float muminus_omega;
    float muminus_tanlambda;
    float muminus_d0sig;
    float muminus_z0sig;
    float muminus_Chi2;
    int   muminus_Ndf;
    float muminus_Chi2Ndf;
    float muminus_dEdx;
    float muminus_dEdxError;
    float muminus_RadiusOfInnermostHit;

    float muplus_ecal;
    float muplus_hcal;
    float muplus_yoke;
    float muplus_lcal;
    float muplus_lhcal;
    float muplus_bcal;
    float muplus_ecalfrac;
    float muplus_EbyP;
    int   muplus_VXD;
    int   muplus_FTD;
    int   muplus_SIT;
    int   muplus_TPC;
    int   muplus_SET;
    int   muplus_ETD;
    float muminus_ecal;
    float muminus_hcal;
    float muminus_yoke;
    float muminus_lcal;
    float muminus_lhcal;
    float muminus_bcal;
    float muminus_ecalfrac;
    float muminus_EbyP;
    int   muminus_VXD;
    int   muminus_FTD;
    int   muminus_SIT;
    int   muminus_TPC;
    int   muminus_SET;
    int   muminus_ETD;

    //MC matching of isolated lepton
    int mc_muplus_PDG;
    int parent_mc_muplus_PDG;
    int mc_muminus_PDG;
    int parent_mc_muminus_PDG;

    float mumu_E;
    float mumu_mass;
    float mumu_costh;
    float mumu_Pt;
    float mumu_acop;
    float sigma_mumu_mass;
    float recoilmass;
    float mumu_costh_tobeam;
    float mumu_mom_mag;

    float muminus_momres;
    float muplus_momres;
    float pseudomass;

    //float mumu_mass_mc;
    //float mumu_E_mc;
    //float mumu_costh_mc;

    //*************************************
    //loop of after isolated letpon tagging
    //*************************************
    int   n_pfoswo;
    float pfowo_Pt[NMAX_PFOS];
    float pfowo_costh[NMAX_PFOS];
    int   n_highPt;
    int   n_highPt_neutral;
    int   n_highPt_charged;

  };

  RECO_DATA _data;

  //ROOT
  TFile* _otfile;
  TTree* _dataTree;
};

#endif



