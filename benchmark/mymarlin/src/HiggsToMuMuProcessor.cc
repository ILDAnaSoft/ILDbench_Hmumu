#include "HiggsToMuMuProcessor.h"

//STL
#include <iostream>
#include <fstream>
#include <set>

#include <vector>
#include <stack>

//LCIO
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "UTIL/PIDHandler.h"
#include "LCIOSTLTypes.h"

//Vertex informaion
#include "VertexInfo.h"

//----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "TRandom.h"

using namespace lcio;
using namespace marlin;
using namespace std;

HiggsToMuMuProcessor aHiggsToMuMuProcessor;

HiggsToMuMuProcessor::HiggsToMuMuProcessor() : Processor( "HiggsToMuMuProcessor" ) {
  
  //modify processor description
  _description = "HiggsToMuMuProcessor does whatever it does ...";
  
  //register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputAllPFOsCollection", 
			   "Name of the PFOs collection",
			   _colAllPFOs,
			   std::string("PandoraPFOs") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputPrimaryVertexCollection", 
			   "Name of the primary vertex collection (from LCFIPlus)",
			   _colPVtx,
			   std::string("PrimaryVertex") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputIsolepsCollection",
			   "Name of the isolated lepton collection",
			   _colIsoleps,
			   std::string("Isoleps") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputPFOsWithoutIsolepsCollection", 
			   "Name of the PFOs after isolated lepton tagging collection",
			   _colPFOsWithoutIsoleps,
			   std::string("PFOsWithoutIsoleps") );

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection", 
			   "Name of the MC particle collection",
			   _colMC,
			   std::string("MCParticle") );

  registerInputCollection( LCIO::LCRELATION,
                           "MCPFORelation",
                           "Relation between MC and PFO particles",
                           _mcpfoRelation,
                           std::string("RecoMCTruthLink"));

  registerProcessorParameter( "EnergeticPt",
			      "For counting energetic charged particles in the collection after isolated lepton tagging",
			      _energetic,
			      float(5) );
  /*
  registerProcessorParameter( "ProcessListFile",
                              "Name of process list file",
                              _processList,
                              std::string("procid500.txt") ) ;
  */
  registerProcessorParameter( "RootFileName",
                              "Name of Root file (default: output.root)",
                              _rootfilename, 
                              std::string("output.root") );
}


void HiggsToMuMuProcessor::init() { 

  //streamlog_out(DEBUG) << "   init called  " << std::endl;
  std::cout << "init called" << std::endl;

  // usually a good idea to
  printParameters();

  // make Ntuple
  makeNTuple();

  // read process list
  //readProcessList( _processList.c_str() );

  _nRun = 0;
  _nEvt = 0;
}

void HiggsToMuMuProcessor::processRunHeader( LCRunHeader* run ) { 
  _nRun++;
} 

void HiggsToMuMuProcessor::processEvent( LCEvent * evt ) { 
  _nEvt++;
  if( _nEvt % 1000 == 0 ) std::cout << "processing event "<< _nEvt << std::endl;

  // this gets called for every event 
  // usually the working horse ...

  // Clear memory
  memset( &_data, 0, sizeof(_data) );
  
  // fill histogram from LCIO data :
  LCCollection* AllPFOs = evt->getCollection( _colAllPFOs ) ;
  LCCollection* PVtx = evt->getCollection( _colPVtx );
  LCCollection* Isoleps = 0;
  try{
    Isoleps = evt->getCollection( _colIsoleps );
  }
  catch(...){
    cout << "no isolated leptons in this event" << endl;
  }
  LCCollection* PFOsWithoutIsoleps = 0;
  try{
    PFOsWithoutIsoleps = evt->getCollection( _colPFOsWithoutIsoleps );
  }
  catch(...){
    cout << "no pfos after isolated lepton tagging in this event" << endl;
  }
  LCCollection* AllMC = evt->getCollection( _colMC ) ;
  _navpfo = new LCRelationNavigator( evt->getCollection( _mcpfoRelation ) );

  //******************************************
  //read out the thrust and sphere information
  //******************************************
  _data.principalthrust = AllPFOs->parameters().getFloatVal( "principleThrustValue" );
  _data.majorthrust     = AllPFOs->parameters().getFloatVal( "majorThrustValue" );
  _data.minorthrust     = AllPFOs->parameters().getFloatVal( "minorThrustValue" );
  FloatVec ThrustAxis;
  AllPFOs->parameters().getFloatVals( "principleThrustAxis" , ThrustAxis );
  //_data.thrustAxisX = ThrustAxis[0];
  //_data.thrustAxisY = ThrustAxis[1];
  //_data.thrustAxisZ = ThrustAxis[2];
  TVector3 thrustvec = TVector3( ThrustAxis[0], ThrustAxis[1], ThrustAxis[2] );
  _data.costh_thrustaxis = TVector3(0,0,1).Dot( thrustvec.Unit() );

  _data.sphericity = AllPFOs->parameters().getFloatVal( "sphericity" );
  _data.aplanarity = AllPFOs->parameters().getFloatVal( "aplanarity" );


  //************
  //MC loop part
  //************
  int nMC = AllMC->getNumberOfElements();
  int nq = 0, ne = 0, nm = 0, nt = 0, nn = 0, no = 0;
  TLorentzVector forff1_4mom(0,0,0,0), forff2_4mom(0,0,0,0);
  TLorentzVector for4f1_4mom(0,0,0,0), for4f2_4mom(0,0,0,0), for4f3_4mom(0,0,0,0), for4f4_4mom(0,0,0,0);
  int count_muminus = 0, count_muplus = 0;
  TLorentzVector MC_muminus_4mom(0,0,0,0), MC_muplus_4mom(0,0,0,0);

  for( int i = 0; i < nMC; i++ ){
    MCParticle* MC = dynamic_cast< MCParticle* >( AllMC->getElementAt(i) );

    //find ISR photons
    if( i == 0 && MC->getPDG() == 22 ){
      TVector3 mom = TVector3( MC->getMomentum() );
      _data.ISR1_E = MC->getEnergy();
      _data.ISR1_costh = TVector3(0,0,1).Dot( mom.Unit() );
    }
    if( i == 1 && MC->getPDG() == 22 ){
      TVector3 mom = TVector3( MC->getMomentum() );
      _data.ISR2_E = MC->getEnergy();
      _data.ISR2_costh = TVector3(0,0,1).Dot( mom.Unit() );
    }

    //find initial particles to distinguish qqh, nnh, etc.
    if( i == 2 ) _data.type = abs( MC->getPDG() );

    //count intial tau decay to muons
    if( _data.type == 15 && i == 2 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ) _data.tau_m++;
    }
    if( _data.type == 15 && i == 3 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ) _data.tau_m++;
    }

    //find Higgs decay
    if( abs( MC->getPDG() ) == 25 && MC->getDaughters().size() == 2 ){
      vector< MCParticle* > HiggsDaughters = MC->getDaughters();
      MCParticle* daughter1 = HiggsDaughters[0];
      MCParticle* daughter2 = HiggsDaughters[1];
      int daughter1PDG = daughter1->getPDG();
      int daughter2PDG = daughter2->getPDG();
      _data.higgsdecay1 = daughter1PDG;
      _data.higgsdecay2 = daughter2PDG;

      //compute generator-level muon pair mass
      if( abs( daughter1PDG ) == 13 && abs( daughter2PDG ) == 13 ){
	TLorentzVector daughter1_4mom = TLorentzVector( daughter1->getMomentum(), daughter1->getEnergy() );
	TLorentzVector daughter2_4mom = TLorentzVector( daughter2->getMomentum(), daughter2->getEnergy() );
	_data.MC_mumu_mass = ( daughter1_4mom + daughter2_4mom ).M();
	_data.MC_mumu_costh = ( daughter1_4mom.Vect().Unit() ).Dot( daughter2_4mom.Vect().Unit() );
      }

      //count number of tau decays to muon
      if( abs( daughter1PDG ) == 15 && abs( daughter2PDG ) == 15 ){
	int tau1decay = getTauDecayMode( daughter1 );
	int tau2decay = getTauDecayMode( daughter2 );
	if( abs( tau1decay ) == 2 ) _data.h_nt_m++;
	if( abs( tau2decay ) == 2 ) _data.h_nt_m++;
      }

      //check real final state when Higgs decays to W/Z
      if( abs( daughter1PDG ) == 23 || abs( daughter1PDG ) == 24 ){
	MCParticle* Bosondaughter1 = daughter1->getDaughters()[0];
	MCParticle* Bosondaughter2 = daughter1->getDaughters()[1];
	MCParticle* Bosondaughter3 = daughter2->getDaughters()[0];
	MCParticle* Bosondaughter4 = daughter2->getDaughters()[1];

	int PDG[4] = {0,0,0,0};
	PDG[0] = abs( Bosondaughter1->getPDG() );
	PDG[1] = abs( Bosondaughter2->getPDG() );
	PDG[2] = abs( Bosondaughter3->getPDG() );
	PDG[3] = abs( Bosondaughter4->getPDG() );

	for( int j = 0; j <= 3; j++ ){
	  if( PDG[j] >= 1 && PDG[j] <= 5 ) _data.h_nq++;
	  else if( PDG[j] == 11 ) _data.h_ne++;
	  else if( PDG[j] == 13 ) _data.h_nm++;
	  else if( PDG[j] == 15 ) _data.h_nt++;
	  else if( PDG[j] == 12 || PDG[j] == 14 || PDG[j] == 16 ) _data.h_nn++;
	}

	if( abs( PDG[0] ) == 15 ){
	  int decay = getTauDecayMode( Bosondaughter1 );
	  if( abs( decay ) == 2 ) _data.h_nboson_m++;
	}
	if( abs( PDG[1] ) == 15 ){
	  int decay = getTauDecayMode( Bosondaughter2 );
	  if( abs( decay ) == 2 ) _data.h_nboson_m++;
	}
	if( abs( PDG[2] ) == 15 ){
	  int decay = getTauDecayMode( Bosondaughter3 );
	  if( abs( decay ) == 2 ) _data.h_nboson_m++;
	}
	if( abs( PDG[3] ) == 15 ){
	  int decay = getTauDecayMode( Bosondaughter4 );
	  if( abs( decay ) == 2 ) _data.h_nboson_m++;
	}


      }

    }

    //find final state of 2f
    if( i == 6 ){
      _data.forff1 = MC->getPDG();
      forff1_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      _data.ff1_costh = TVector3(0,0,1).Dot( forff1_4mom.Vect().Unit() );
      if( abs( MC->getPDG() ) == 15 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ) _data.ff_nt_m++;
      }
    }
    if( i == 7 ){
      _data.forff2 = MC->getPDG();
      forff2_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      _data.ff2_costh = TVector3(0,0,1).Dot( forff2_4mom.Vect().Unit() );
      if( abs( MC->getPDG() ) == 15 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ) _data.ff_nt_m++;
      }
    }

    //find final state of 4f
    if( i == 6 ){
      _data.forffff1 = MC->getPDG();
      for4f1_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      _data.ffff1_costh = TVector3(0,0,1).Dot( for4f1_4mom.Vect().Unit() );
      if( abs( MC->getPDG() ) == 15 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ){
	  _data.nt_m++;
	  vector< MCParticle* >::const_iterator iter;
	  for( iter = MC->getDaughters().begin();
	       iter != MC->getDaughters().end();
	       ++iter ){
	    if( abs( (*iter)->getPDG() ) == 13 ){
	      TLorentzVector iter_4mom = TLorentzVector( (*iter)->getMomentum(), (*iter)->getEnergy() );
	      _data.ffff1_t_costh = TVector3(0,0,1).Dot( iter_4mom.Vect().Unit() );
	    }
	  }
	}
      }
    }
    if( i == 7 ){
      _data.forffff2 = MC->getPDG();
      for4f2_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      _data.ffff2_costh = TVector3(0,0,1).Dot( for4f2_4mom.Vect().Unit() );
      if( abs( MC->getPDG() ) == 15 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ){
	  _data.nt_m++;
	  vector< MCParticle* >::const_iterator iter;
	  for( iter = MC->getDaughters().begin();
	       iter != MC->getDaughters().end();
	       ++iter ){
	    if( abs( (*iter)->getPDG() ) == 13 ){
	      TLorentzVector iter_4mom = TLorentzVector( (*iter)->getMomentum(), (*iter)->getEnergy() );
	      _data.ffff2_t_costh = TVector3(0,0,1).Dot( iter_4mom.Vect().Unit() );
	    }
	  }
	}
      }
    }
    if( i == 8 ){
      _data.forffff3 = MC->getPDG();
      for4f3_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      _data.ffff3_costh = TVector3(0,0,1).Dot( for4f3_4mom.Vect().Unit() );
      if( abs( MC->getPDG() ) == 15 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ){
	  _data.nt_m++;
	  vector< MCParticle* >::const_iterator iter;
	  for( iter = MC->getDaughters().begin();
	       iter != MC->getDaughters().end();
	       ++iter ){
	    if( abs( (*iter)->getPDG() ) == 13 ){
	      TLorentzVector iter_4mom = TLorentzVector( (*iter)->getMomentum(), (*iter)->getEnergy() );
	      _data.ffff3_t_costh = TVector3(0,0,1).Dot( iter_4mom.Vect().Unit() );
	    }
	  }
	}
      }
    }
    if( i == 9 ){
      _data.forffff4 = MC->getPDG();
      for4f4_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      _data.ffff4_costh = TVector3(0,0,1).Dot( for4f4_4mom.Vect().Unit() );
      if( abs( MC->getPDG() ) == 15 ){
	int decay = getTauDecayMode( MC );
	if( abs( decay ) == 2 ){
	  _data.nt_m++;
	  vector< MCParticle* >::const_iterator iter;
	  for( iter = MC->getDaughters().begin();
	       iter != MC->getDaughters().end();
	       ++iter ){
	    if( abs( (*iter)->getPDG() ) == 13 ){
	      TLorentzVector iter_4mom = TLorentzVector( (*iter)->getMomentum(), (*iter)->getEnergy() );
	      _data.ffff4_t_costh = TVector3(0,0,1).Dot( iter_4mom.Vect().Unit() );
	    }
	  }
	}
      }
    }
    if( i >= 6 && i <= 9 ){
      int tempPDG = abs( MC->getPDG() ); 
      if( tempPDG >= 1 && tempPDG <= 5 ) nq++;
      else if( tempPDG == 11 ) ne++;
      else if( tempPDG == 13 ) nm++;
      else if( tempPDG == 15 ) nt++;
      else if( tempPDG == 12 || tempPDG == 14 || tempPDG == 16 ) nn++;
      else no++;
    }
    _data.nq = nq; _data.ne = ne; _data.nm = nm; _data.nt = nt; _data.nn = nn; _data.no = no;
    _data.massff = ( forff1_4mom + forff2_4mom ).M();
    _data.mass1 = ( for4f1_4mom + for4f2_4mom ).M();
    _data.mass2 = ( for4f3_4mom + for4f4_4mom ).M();

    //check MC muons and its pair mass
    if( MC->getPDG() == 13 && count_muminus == 0 ){
      MC_muminus_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      if( MC->getParents().size() > 0 ){
	MCParticle* parent = MC->getParents()[0];
	_data.MC_muminus_parent = parent->getPDG();
      }
      count_muminus++;
    }
    else if( MC->getPDG() == -13 && count_muplus == 0 ){
      MC_muplus_4mom = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
      if( MC->getParents().size() > 0 ){
	MCParticle* parent = MC->getParents()[0];
	_data.MC_muplus_parent = parent->getPDG();
      }
      count_muplus++;
    }

  }
  _data.MC_mumupair_mass = ( MC_muminus_4mom + MC_muplus_4mom ).M();

  //quick look
  for( int i = 0; i < nMC; i++ ){
    MCParticle *MC = dynamic_cast< MCParticle* >( AllMC->getElementAt(i) );
    if( abs(MC->getPDG()) == 13 ){
      _data.MC_mu_E = MC->getEnergy();
    }

  }

  /*
  //check initial Z boson decay
  //note: in DBD samples, there are no Z bosons in initial state, it means start from e+e- -> e+e-gamgam -> mumutautaugamgam something like this. Z boson is not appeared in MC information. I need to take care of it.
  TLorentzVector muminus_4mom_MC(0,0,0,0), muplus_4mom_MC(0,0,0,0);
  float MC_Z_mumu_mass = 0;
  int counter1 = 0, counter2 = 0;
  for( int i = 0; i < nMC; i++ ){
    MCParticle* MC = dynamic_cast< MCParticle* >( AllMC->getElementAt(i) );
    if( MC->getParents().size() > 0 && MC->getParents()[0]->getPDG() != 25 ){//select particles not from Higgs
      if( MC->getPDG() == 13 ){//mu-
	muminus_4mom_MC = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
	counter1++;
      }
      if( MC->getPDG() == -13 ){//mu+
	muplus_4mom_MC = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
	counter2++;
      }
    }
    if( counter1 == 1 && counter2 == 1 ) break;
  }
  MC_Z_mumu_mass = ( muminus_4mom_MC + muplus_4mom_MC ).M();
  _data.MC_Z_mumu_mass = MC_Z_mumu_mass;
  */
  /*
  // map to store tau decay modes
  set<int> mcTauDecayModes;
  for( int i = 0; i < nMC; i++ ){
    MCParticle* MC = dynamic_cast< MCParticle* >( AllMC->getElementAt(i) );
    int mode = getTauDecayMode( MC );
    if(mode != 0) mcTauDecayModes.insert(mode);
  }
  _data.ntaumode = 0;

  for( set<int>::iterator iter = mcTauDecayModes.begin();
       iter != mcTauDecayModes.end();
       ++iter ){
    _data.taumode[_data.ntaumode ] = *iter;
    ++_data.ntaumode;
  }
  */


  //*********************************************
  //Check the relation between MCParticle and PFO
  //*********************************************
  float esum_mc = 0, overlay_E_mc = 0;
  float esum_charged_mc = 0, esum_neutral_mc = 0;
  float overlay_charged_mc = 0, overlay_neutral_mc = 0;
  TVector3 psum_mc(0,0,0);
  int n_highPt_mc = 0;
  TLorentzVector muminus_mom4_mc(0,0,0,0), muplus_mom4_mc(0,0,0,0);

  for( int i = 0; i < nMC; i++ ){
    MCParticle *MC = dynamic_cast< MCParticle* >( AllMC->getElementAt(i) );
    double energyDiff = 9999;
    if( _navpfo->getRelatedFromObjects( MC ).size() > 0 && MC->getGeneratorStatus() == 1 ){//only select MCParticle with related PFO and stable MCParticle (require generator status == 1)
      ReconstructedParticle *pfo = 0;
      LCObjectVec vec = _navpfo->getRelatedFromObjects( MC );
      for( unsigned int a = 0; a < vec.size(); a++ ){
	ReconstructedParticle *recpfo = dynamic_cast< ReconstructedParticle* >( vec[a] );
	if( recpfo->getCharge() == MC->getCharge() && TMath::Abs( recpfo->getEnergy() - MC->getEnergy() ) < energyDiff ){
	  energyDiff = TMath::Abs( recpfo->getEnergy() - MC->getEnergy() );
	  pfo = recpfo;
	}
      }

      //check general variables
      esum_mc += MC->getEnergy();
      psum_mc += MC->getMomentum();
      if( MC->isOverlay() == true ) _data.overlay_mc[i] = 1;
      if( MC->isOverlay() == true ) overlay_E_mc += MC->getEnergy();
      if( MC->getCharge() == 0 ){//neutral
	esum_neutral_mc += MC->getEnergy();
	if( MC->isOverlay() == true ){//overlay
	  overlay_neutral_mc = MC->getEnergy();
	}
      }
      else if( MC->getCharge() != 0 ){//charged
	esum_charged_mc += MC->getEnergy();
	if( MC->isOverlay() == true ){//overlay
	  overlay_charged_mc = MC->getEnergy();
	}
      }

      //check without isolated muon collection
      if( PFOsWithoutIsoleps != 0 ){
	int npfowo = PFOsWithoutIsoleps->getNumberOfElements();
	for( int j = 0; j < npfowo; j++ ){
	  ReconstructedParticle *pfo_temp = dynamic_cast< ReconstructedParticle* >( PFOsWithoutIsoleps->getElementAt(j) );
	  if( pfo == pfo_temp ){
	    TVector3 mom = MC->getMomentum();
	    float Pt = TMath::Sqrt( mom[0]*mom[0] + mom[1]*mom[1] );
	    TVector3 pfomom = pfo_temp->getMomentum();
	    float pfoPt = TMath::Sqrt( pfomom[0]*pfomom[0] + pfomom[1]*pfomom[1] );
	    if( Pt > _energetic && MC->getCharge() != 0 && pfo_temp->getCharge() != 0 ){
	      n_highPt_mc++;
	      _data.Pt_mc[j] = Pt;
	      _data.costh_mc[j] = TVector3(0,0,1).Dot( mom.Unit() );
	      _data.pfoPt[j] = pfoPt;
	      _data.pfocosth[j] = TVector3(0,0,1).Dot( pfomom.Unit() );
	    }
	  }
	}
      }

      //check isolated muons
      if( Isoleps != 0 ){
	int nisolep = Isoleps->getNumberOfElements();
	for( int k = 0; k < nisolep; k++ ){
	  ReconstructedParticle *isolep_temp = dynamic_cast< ReconstructedParticle* >( Isoleps->getElementAt(k) );
	  if( pfo == isolep_temp ){
	    TLorentzVector mom4 = TLorentzVector( MC->getMomentum(), MC->getEnergy() );
	    if( MC->getCharge() > 0 ){
	      muplus_mom4_mc = mom4;
	      _data.n_muplus_mc++;
	    }
	    if( MC->getCharge() < 0 ){
	      muminus_mom4_mc = mom4;
	      _data.n_muminus_mc++;
	    }
	  }
	}
      }

    }
  }
  _data.esum_mc = esum_mc;
  _data.esum_charged_mc = esum_charged_mc;
  _data.esum_neutral_mc = esum_neutral_mc;
  _data.overlay_E_mc = overlay_E_mc;
  _data.overlay_charged_mc = overlay_charged_mc;
  _data.overlay_neutral_mc = overlay_neutral_mc;
  _data.costh_missmom_mc = TVector3(0,0,1).Dot( -psum_mc.Unit() );
  _data.Pt_all_mc = TMath::Sqrt( psum_mc[0]*psum_mc[0] + psum_mc[1]*psum_mc[1] );

  _data.n_highPt_mc = n_highPt_mc;

  _data.mumu_mass_mc = ( muminus_mom4_mc + muplus_mom4_mc ).M();
  _data.mumu_costh_mc = ( muminus_mom4_mc.Vect().Unit() ).Dot( muplus_mom4_mc.Vect().Unit() );

  //test for smearing
  float pT_muminus = TMath::Sqrt( muminus_mom4_mc[0]*muminus_mom4_mc[0]
				 +muminus_mom4_mc[1]*muminus_mom4_mc[1] );
  float pT_muplus  = TMath::Sqrt( muplus_mom4_mc[0]*muplus_mom4_mc[0]
				 +muplus_mom4_mc[1]*muplus_mom4_mc[1] );

  float resolution = 1E-3;
  float mom_res = TMath::Sqrt( resolution );
  TLorentzVector smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
						 gRandom->Gaus(0,mom_res*pT_muminus),
						 0,
						 0 );
  TLorentzVector smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
						gRandom->Gaus(0,mom_res*pT_muplus),
						0,
						0 );
  TLorentzVector muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  TLorentzVector muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_1_3 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 5E-4;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_5_4 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 3E-4;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_3_4 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 2E-4;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_2_4 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 1E-4;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_1_4 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 5E-5;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_5_5 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 3E-5;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_3_5 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 2E-5;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_2_5 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 1E-5;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_1_5 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 5E-6;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_5_6 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 3E-6;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_3_6 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 2E-6;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_2_6 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();

  resolution = 1E-6;
  mom_res = TMath::Sqrt( resolution );
  smear_muminus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muminus),
				  gRandom->Gaus(0,mom_res*pT_muminus),
				  0,
				  0 );
  smear_muplus = TLorentzVector( gRandom->Gaus(0,mom_res*pT_muplus),
				 gRandom->Gaus(0,mom_res*pT_muplus),
				 0,
				 0 );
  muminus_mom4_mc_smear = muminus_mom4_mc + smear_muminus;
  muplus_mom4_mc_smear  = muplus_mom4_mc + smear_muplus;
  _data.mumu_mass_mc_smear_1_6 = ( muminus_mom4_mc_smear + muplus_mom4_mc_smear ).M();


  //*************
  //PFO loop part
  //*************
  const float CM = 500.;
  const float initial_Px = CM * TMath::Sin(0.007);
  TLorentzVector CM_4mom(initial_Px,0,0,CM);
  int npfos = AllPFOs->getNumberOfElements();
  TVector3 calcPt(0,0,0);
  float maxPt = 0, Pt = 0;
  int ntrksum = 0;
  float esum = 0; TVector3 psum(0,0,0);
  for( int i = 0; i < npfos; i++ ){
    ReconstructedParticle* pfo = dynamic_cast< ReconstructedParticle* >( AllPFOs->getElementAt(i) );
    const EVENT::TrackVec & trkvec = pfo->getTracks();
    const EVENT::ClusterVec& clusvec = pfo->getClusters();

    //use for z0 correction due to z-position smearing
    Vertex *pvtx = dynamic_cast< Vertex* >( PVtx->getElementAt(0) );
    Double_t z_pvtx = pvtx->getPosition()[2];

    //number of tracks
    ntrksum += trkvec.size();

    //basic variables
    esum += pfo->getEnergy();
    psum += TVector3( pfo->getMomentum() );
    calcPt += TVector3( pfo->getMomentum() );

    //calculating Pt for charged particle, and find maximum Pt among them
    float Ptcal = 0;
    if( pfo->getCharge() != 0 ){
      float Px = pfo->getMomentum()[0];
      float Py = pfo->getMomentum()[1];
      Ptcal = TMath::Sqrt( Px*Px + Py*Py );
    }
    if( Ptcal > maxPt ) maxPt = Ptcal;

    //store information of PFO itself
    _data.pfo_e[i]    = pfo->getEnergy();
    _data.pfo_px[i]   = pfo->getMomentum()[0];
    _data.pfo_py[i]   = pfo->getMomentum()[1];
    _data.pfo_pz[i]   = pfo->getMomentum()[2];
    _data.pfo_mom_mag[i] = TVector3( pfo->getMomentum() ).Mag();
    _data.pfo_chrg[i] = pfo->getCharge();
    _data.pfo_pdg[i]  = pfo->getType();

    //calculating Pt for every particle and summing up
    float forPx = pfo->getMomentum()[0];
    float forPy = pfo->getMomentum()[1];
    _data.pfo_Pt[i] = TMath::Sqrt( forPx*forPx + forPy*forPy );
    Pt += TMath::Sqrt( forPx*forPx + forPy*forPy );
    
    TVector3 pfo_mom = TVector3( pfo->getMomentum() );
    _data.pfo_costh[i] = TVector3(0,0,1).Dot( pfo_mom.Unit() );
    _data.pfo_ntrk[i] = trkvec.size();
    _data.pfo_nclus[i] = clusvec.size();

    //fill first track info
    if( trkvec.size() > 0 ){
      const Track* trk = trkvec[0];
      _data.pfo_d0[i] = trk->getD0();
      //_data.pfo_z0[i] = trk->getZ0();
      _data.pfo_z0[i] = trk->getZ0() - z_pvtx;
      _data.pfo_phi[i] = trk->getPhi();
      _data.pfo_omega[i] = trk->getOmega();
      _data.pfo_tanlambda[i] = trk->getTanLambda();
      _data.pfo_d0sig[i] = trk->getD0()/sqrt(trk->getCovMatrix()[0]);
      //_data.pfo_z0sig[i] = trk->getZ0()/sqrt(trk->getCovMatrix()[9]);
      _data.pfo_z0sig[i] = ( trk->getZ0() - z_pvtx )/sqrt(trk->getCovMatrix()[9]);

      _data.pfo_Chi2[i] = trk->getChi2();
      _data.pfo_Ndf[i]  = trk->getNdf();
      if( trk->getNdf() == 0 ) _data.pfo_Chi2Ndf[i] = 0;
      else _data.pfo_Chi2Ndf[i] = trk->getChi2() / float( trk->getNdf() );
      _data.pfo_dEdx[i] = trk->getdEdx();
      _data.pfo_dEdxError[i] = trk->getdEdxError();
      _data.pfo_RadiusOfInnermostHit[i] = trk->getRadiusOfInnermostHit();

      _data.pfo_VXD[i] = trk->getSubdetectorHitNumbers()[0];
      _data.pfo_FTD[i] = trk->getSubdetectorHitNumbers()[2];
      _data.pfo_SIT[i] = trk->getSubdetectorHitNumbers()[4];
      _data.pfo_TPC[i] = trk->getSubdetectorHitNumbers()[6];
      _data.pfo_SET[i] = trk->getSubdetectorHitNumbers()[8];
      _data.pfo_ETD[i] = trk->getSubdetectorHitNumbers()[10];
    }

    //fill sum of cluster info
    if( clusvec.size() > 0 ){
      for( ClusterVec::const_iterator iCluster=clusvec.begin();
	   iCluster!=clusvec.end();
	   ++iCluster) {
	_data.pfo_ecal[i]  += (*iCluster)->getSubdetectorEnergies()[0];
	_data.pfo_hcal[i]  += (*iCluster)->getSubdetectorEnergies()[1];
	_data.pfo_yoke[i]  += (*iCluster)->getSubdetectorEnergies()[2];
	_data.pfo_lcal[i]  += (*iCluster)->getSubdetectorEnergies()[3];
	_data.pfo_lhcal[i] += (*iCluster)->getSubdetectorEnergies()[4];
	_data.pfo_bcal[i]  += (*iCluster)->getSubdetectorEnergies()[5];
	
	_data.pfo_ecalfrac[i] = _data.pfo_ecal[i] / ( _data.pfo_ecal[i] + _data.pfo_hcal[i] );
	_data.pfo_EbyP[i]     = ( _data.pfo_ecal[i] + _data.pfo_hcal[i] ) / ( pfo_mom.Mag() );
      }
    }
    
    //matching with MC
    MCParticle* pfo_mcp = 0; 
    if ( _navpfo->getRelatedToWeights( pfo ).size() > 0 ) {
      pfo_mcp = dynamic_cast< MCParticle* > ( _navpfo->getRelatedToObjects( pfo )[0] );

      _data.pfo_pid_mc[i] = pfo_mcp->getPDG();
      _data.pfo_E_mc[i]   = pfo_mcp->getEnergy();
      _data.pfo_Px_mc[i]  = pfo_mcp->getMomentum()[0];
      _data.pfo_Py_mc[i]  = pfo_mcp->getMomentum()[1];
      _data.pfo_Pz_mc[i]  = pfo_mcp->getMomentum()[2];

      //_data.pfo_taudecaymode[i] = getTauDecayMode( pfo_mcp );

      /*
      //check overlay
      if( pfo_mcp->isOverlay() ){
	_data.count_overlay++;
	if( pfo_mcp->getCharge() != 0 ){
	  _data.charged_overlay++;
	}
      }
      */

      //check its parent
      if( pfo_mcp->getParents().size() > 0 ){
	MCParticle *parent = pfo_mcp->getParents()[0];
	_data.parent_pid_mc[i] =  parent->getPDG();
	_data.parent_E_mc[i]   =  parent->getEnergy();
	_data.parent_Px_mc[i]  =  parent->getMomentum()[0];
	_data.parent_Py_mc[i]  =  parent->getMomentum()[1];
	_data.parent_Pz_mc[i]  =  parent->getMomentum()[2];

	//finding original mother
	while( parent->getParents().size() > 0 && parent->getPDG() != 25 ){
	  parent = parent->getParents()[0];
	}
	_data.pfo_mother_pid[i] = parent->getPDG();
      }
    }
    
  }
  _data.npfos = npfos;
  _data.ntrks = ntrksum;
  Vertex *vtx = dynamic_cast< Vertex* >( PVtx->getElementAt(0) );
  _data.vtx_x = vtx->getPosition()[0];
  _data.vtx_y = vtx->getPosition()[1];
  _data.vtx_z = vtx->getPosition()[2];

  _data.esum = esum;
  _data.psum_x = psum[0];
  _data.psum_y = psum[1];
  _data.psum_z = psum[2];
  _data.maxPt = maxPt;
  _data.Pt = Pt;
  _data.Pt_all = TMath::Sqrt( calcPt[0]*calcPt[0] + calcPt[1]*calcPt[1] );
  TLorentzVector visible_4mom = TLorentzVector( psum[0], psum[1], psum[2], esum );
  _data.mass_vis = ( visible_4mom ).M();
  _data.misse = CM - esum;
  _data.missmom_x = -psum[0];
  _data.missmom_y = -psum[1];
  _data.missmom_z = -psum[2];
  _data.missmom_mag = psum.Mag();
  _data.costh_missmom = TVector3(0,0,1).Dot( -psum.Unit() );
  TLorentzVector invisible_4mom = TLorentzVector( -psum[0], -psum[1], -psum[2], CM - esum );
  _data.mass_inv = ( invisible_4mom ).M();


  //*************************
  //Isolated lepton loop part
  //*************************
  int n_isoleps = 0, n_muminus = 0, n_muplus = 0;
  TLorentzVector muminus_4mom(0,0,0,0), muplus_4mom(0,0,0,0);
  TVector3 muminus_3mom(0,0,0), muplus_3mom(0,0,0);
  float muminus_E = 0, muplus_E = 0, muminus_Pt = 0, muplus_Pt = 0, muminus_costh = 0, muplus_costh = 0;
  float muminus_azi = 0, muplus_azi = 0, muminus_mom_mag = 0, muplus_mom_mag = 0;
  float muminus_charge = 0, muplus_charge = 0;
  //for covariance matrix in momenta space
  float muplus_covMat0 = 0, muplus_covMat1 = 0, muplus_covMat2 = 0, muplus_covMat3 = 0, muplus_covMat4 = 0;
  float muplus_covMat5 = 0, muplus_covMat6 = 0, muplus_covMat7 = 0, muplus_covMat8 = 0, muplus_covMat9 = 0;
  float muminus_covMat0 = 0, muminus_covMat1 = 0, muminus_covMat2 = 0, muminus_covMat3 = 0;
  float muminus_covMat4 = 0, muminus_covMat5 = 0, muminus_covMat6 = 0, muminus_covMat7 = 0;
  float muminus_covMat8 = 0, muminus_covMat9 = 0;
  TLorentzVector muminus_4mom_mc(0,0,0,0), muplus_4mom_mc(0,0,0,0);
  Track *muplus_track = 0, *muminus_track = 0;
  if( Isoleps != 0 ){
    IntVec types;
    IntVec isoleptypes = Isoleps->getParameters().getIntVals( "ISOLepType", types );
    n_isoleps = Isoleps->getNumberOfElements();

    for( int i = 0; i < n_isoleps; i++ ){
      ReconstructedParticle* isolep = dynamic_cast< ReconstructedParticle* >( Isoleps->getElementAt(i) );
      int leptype = isoleptypes[i];

      //for correction of z0 due to z-position smearing
      Vertex *pvtx = dynamic_cast< Vertex* >( PVtx->getElementAt(0) );
      Double_t z_pvtx = pvtx->getPosition()[2];
      
      //preparing for MC matching
      MCParticle* mc_isolep = 0;
      if( _navpfo->getRelatedToWeights( isolep ).size() > 0 )
        mc_isolep = dynamic_cast< MCParticle* >( _navpfo->getRelatedToObjects( isolep )[0] );
      
      if( leptype == 13 ){//only selecting isolated muons
        const EVENT::TrackVec & trkvec = isolep->getTracks();
        const EVENT::ClusterVec & clusvec = isolep->getClusters();
	
	if( isolep->getCharge() > 0 ){//mu+
          n_muplus++;
          muplus_4mom = TLorentzVector( isolep->getMomentum(), isolep->getEnergy() );
          muplus_3mom = TVector3( isolep->getMomentum() );
          muplus_charge = float( isolep->getCharge() );
          muplus_E = isolep->getEnergy();
          muplus_Pt = TMath::Sqrt( muplus_3mom[0]*muplus_3mom[0] + muplus_3mom[1]*muplus_3mom[1] );
          muplus_costh = TVector3(0,0,1).Dot( muplus_3mom.Unit() );
	  muplus_azi = muplus_3mom.Phi();
          muplus_mom_mag = muplus_3mom.Mag();

	  //covariance matrix in momenta space
          muplus_covMat0 = isolep->getCovMatrix()[0];
          muplus_covMat1 = isolep->getCovMatrix()[1];
          muplus_covMat2 = isolep->getCovMatrix()[2];
          muplus_covMat3 = isolep->getCovMatrix()[3];
          muplus_covMat4 = isolep->getCovMatrix()[4];
          muplus_covMat5 = isolep->getCovMatrix()[5];
          muplus_covMat6 = isolep->getCovMatrix()[6];
          muplus_covMat7 = isolep->getCovMatrix()[7];
          muplus_covMat8 = isolep->getCovMatrix()[8];
          muplus_covMat9 = isolep->getCovMatrix()[9];

	  if( trkvec.size() > 0 ){//track parameters of mu+
	    muplus_track = trkvec[0];
	    const Track* trk = trkvec[0];
            _data.muplus_d0        = trk->getD0();
            //_data.muplus_z0        = trk->getZ0();
            _data.muplus_z0        = trk->getZ0() - z_pvtx;
            _data.muplus_phi       = trk->getPhi();
            _data.muplus_omega     = trk->getOmega();
            _data.muplus_tanlambda = trk->getTanLambda();
            _data.muplus_d0sig     = trk->getD0() / ( TMath::Sqrt( trk->getCovMatrix()[0] ) );
            //_data.muplus_z0sig     = trk->getZ0() / ( TMath::Sqrt( trk->getCovMatrix()[9] ) );
            _data.muplus_z0sig     = ( trk->getZ0() - z_pvtx ) / ( TMath::Sqrt( trk->getCovMatrix()[9] ) );

            _data.muplus_Chi2      = trk->getChi2();
            _data.muplus_Ndf       = trk->getNdf();
            if( trk->getNdf() == 0 ) _data.muplus_Chi2Ndf = 0;
            else _data.muplus_Chi2Ndf = trk->getChi2() / float( trk->getNdf() );
            _data.muplus_dEdx      = trk->getdEdx();
            _data.muplus_dEdxError = trk->getdEdxError();
            _data.muplus_RadiusOfInnermostHit = trk->getRadiusOfInnermostHit();

            _data.muplus_VXD = trk->getSubdetectorHitNumbers()[0];
            _data.muplus_FTD = trk->getSubdetectorHitNumbers()[2];
            _data.muplus_SIT = trk->getSubdetectorHitNumbers()[4];
            _data.muplus_TPC = trk->getSubdetectorHitNumbers()[6];
            _data.muplus_SET = trk->getSubdetectorHitNumbers()[8];
            _data.muplus_ETD = trk->getSubdetectorHitNumbers()[10];
	  }

	  if( clusvec.size() > 0 ){//clsuter parameters of mu+
            for( ClusterVec::const_iterator iClus = clusvec.begin();
                 iClus != clusvec.end();
                 ++iClus ){
              _data.muplus_ecal  += ( *iClus )->getSubdetectorEnergies()[0];
              _data.muplus_hcal  += ( *iClus )->getSubdetectorEnergies()[1];
              _data.muplus_yoke  += ( *iClus )->getSubdetectorEnergies()[2];
              _data.muplus_lcal  += ( *iClus )->getSubdetectorEnergies()[3];
              _data.muplus_lhcal += ( *iClus )->getSubdetectorEnergies()[4];
              _data.muplus_bcal  += ( *iClus )->getSubdetectorEnergies()[5];
              _data.muplus_ecalfrac = _data.muplus_ecal / ( _data.muplus_ecal + _data.muplus_hcal );
              _data.muplus_EbyP     = ( _data.muplus_ecal + _data.muplus_hcal ) / ( muplus_4mom.Vect().Mag() );
	    }
	  }

	  //MC matching
          _data.mc_muplus_PDG = mc_isolep->getPDG();
	  muplus_4mom_mc = TLorentzVector( mc_isolep->getMomentum(), mc_isolep->getEnergy() );
          MCParticle* parenttest = 0;
          if( mc_isolep->getParents().size() > 0 ) parenttest = mc_isolep;
	  if( parenttest != 0 ){
	    while( parenttest->getParents().size() > 0 ){
	      parenttest = parenttest->getParents()[0];
	      int PDG = parenttest->getPDG();
	      _data.parent_mc_muplus_PDG = PDG;
	      if( ( abs( PDG ) == 15 ) || ( abs( PDG ) == 25 ) ) break;
	    }
	  }
	
	}
	else if( isolep->getCharge() < 0 ){//mu-
          n_muminus++;
          muminus_4mom = TLorentzVector( isolep->getMomentum(), isolep->getEnergy() );
          muminus_3mom = TVector3( isolep->getMomentum() );
          muminus_charge = float( isolep->getCharge() );
          muminus_E = isolep->getEnergy();
          muminus_Pt = TMath::Sqrt( muminus_3mom[0]*muminus_3mom[0] + muminus_3mom[1]*muminus_3mom[1] );
          muminus_costh = TVector3(0,0,1).Dot( muminus_3mom.Unit() );
	  muminus_azi = muminus_3mom.Phi();
          muminus_mom_mag = muminus_3mom.Mag();

	  //covariance matrix in momenta space
          muminus_covMat0 = isolep->getCovMatrix()[0];
          muminus_covMat1 = isolep->getCovMatrix()[1];
          muminus_covMat2 = isolep->getCovMatrix()[2];
          muminus_covMat3 = isolep->getCovMatrix()[3];
          muminus_covMat4 = isolep->getCovMatrix()[4];
          muminus_covMat5 = isolep->getCovMatrix()[5];
          muminus_covMat6 = isolep->getCovMatrix()[6];
          muminus_covMat7 = isolep->getCovMatrix()[7];
          muminus_covMat8 = isolep->getCovMatrix()[8];
          muminus_covMat9 = isolep->getCovMatrix()[9];

	  if( trkvec.size() > 0 ){//track parameters of mu-
	    muminus_track = trkvec[0];
            const Track* trk = trkvec[0];
            _data.muminus_d0        = trk->getD0();
            //_data.muminus_z0        = trk->getZ0();
            _data.muminus_z0        = trk->getZ0() - z_pvtx;
            _data.muminus_phi       = trk->getPhi();
            _data.muminus_omega     = trk->getOmega();
            _data.muminus_tanlambda = trk->getTanLambda();
            _data.muminus_d0sig     = trk->getD0() / ( TMath::Sqrt( trk->getCovMatrix()[0] ) );
            //_data.muminus_z0sig     = trk->getZ0() / ( TMath::Sqrt( trk->getCovMatrix()[9] ) );
            _data.muminus_z0sig     = ( trk->getZ0() - z_pvtx ) / ( TMath::Sqrt( trk->getCovMatrix()[9] ) );

            _data.muminus_Chi2      = trk->getChi2();
            _data.muminus_Ndf       = trk->getNdf();
            if( trk->getNdf() == 0 ) _data.muminus_Chi2Ndf = 0;
            else _data.muminus_Chi2Ndf   = trk->getChi2() / float( trk->getNdf() );
            _data.muminus_dEdx      = trk->getdEdx();
            _data.muminus_dEdxError = trk->getdEdxError();
            _data.muminus_RadiusOfInnermostHit = trk->getRadiusOfInnermostHit();

            _data.muminus_VXD = trk->getSubdetectorHitNumbers()[0];
            _data.muminus_FTD = trk->getSubdetectorHitNumbers()[2];
            _data.muminus_SIT = trk->getSubdetectorHitNumbers()[4];
            _data.muminus_TPC = trk->getSubdetectorHitNumbers()[6];
            _data.muminus_SET = trk->getSubdetectorHitNumbers()[8];
            _data.muminus_ETD = trk->getSubdetectorHitNumbers()[10];
	  }

	  if( clusvec.size() > 0 ){//cluster parameters of mu-
            for( ClusterVec::const_iterator iClus = clusvec.begin();
                 iClus != clusvec.end();
                 ++iClus ){
              _data.muminus_ecal  += ( *iClus )->getSubdetectorEnergies()[0];
              _data.muminus_hcal  += ( *iClus )->getSubdetectorEnergies()[1];
              _data.muminus_yoke  += ( *iClus )->getSubdetectorEnergies()[2];
              _data.muminus_lcal  += ( *iClus )->getSubdetectorEnergies()[3];
              _data.muminus_lhcal += ( *iClus )->getSubdetectorEnergies()[4];
              _data.muminus_bcal  += ( *iClus )->getSubdetectorEnergies()[5];
              _data.muminus_ecalfrac = _data.muminus_ecal / ( _data.muminus_ecal + _data.muminus_hcal );
              _data.muminus_EbyP     = ( _data.muminus_ecal + _data.muminus_hcal ) / ( muminus_4mom.Vect().Mag() );
            }
	  }

	  //MC matching
          _data.mc_muminus_PDG = mc_isolep->getPDG();
	  muminus_4mom_mc = TLorentzVector( mc_isolep->getMomentum(), mc_isolep->getEnergy() );
          MCParticle* parenttest = 0;
          if( mc_isolep->getParents().size() > 0 ) parenttest = mc_isolep;
	  if( parenttest != 0 ){
	    while( parenttest->getParents().size() > 0 ){
	      parenttest = parenttest->getParents()[0];
	      int PDG = parenttest->getPDG();
	      _data.parent_mc_muminus_PDG = PDG;
	      if( ( abs( PDG ) == 15 ) || ( abs( PDG ) == 25 ) ) break;
	    }
	  }

	}
      
      }
      
    }
  }
  _data.n_isoleps = n_isoleps; _data.n_muplus = n_muplus; _data.n_muminus = n_muminus;
  _data.muplus_E = muplus_E; _data.muminus_E = muminus_E;
  _data.muplus_Pt = muplus_Pt; _data.muminus_Pt = muminus_Pt;
  _data.muplus_costh = muplus_costh; _data.muminus_costh = muminus_costh;
  _data.muplus_azi = muplus_azi; _data.muminus_azi = muminus_azi;
  _data.muplus_mom_mag = muplus_mom_mag; _data.muminus_mom_mag = muminus_mom_mag;
  _data.muplus_charge_costh  = muplus_charge  * muplus_costh;
  _data.muminus_charge_costh = muminus_charge * muminus_costh;
  _data.sum_charge_costh = _data.muplus_charge_costh + _data.muminus_charge_costh;
  if( muplus_E > muminus_E ){
    _data.leadingmu_E = muplus_E; _data.subleadingmu_E = muminus_E;
    _data.leadingmu_costh = muplus_costh; _data.subleadingmu_costh = muminus_costh;
  }
  else{
    _data.leadingmu_E = muminus_E; _data.subleadingmu_E = muplus_E;
    _data.leadingmu_costh = muminus_costh; _data.subleadingmu_costh = muplus_costh;
  }

  TLorentzVector mupair_4mom = ( muplus_4mom + muminus_4mom );
  TVector3 mupair_3mom = ( muplus_3mom + muminus_3mom );
  _data.mumu_mass = mupair_4mom.M();
  _data.mumu_E = mupair_4mom.E();
  _data.mumu_Pt = TMath::Sqrt( mupair_3mom[0]*mupair_3mom[0] + mupair_3mom[1]*mupair_3mom[1] );
  _data.mumu_mom_mag = mupair_3mom.Mag();
  _data.mumu_costh = ( muplus_3mom.Unit() ).Dot( muminus_3mom.Unit() );
  _data.mumu_costh_tobeam = TVector3(0,0,1).Dot( mupair_3mom.Unit() );
  _data.mumu_acop = TMath::Cos( muminus_3mom.Phi() - muplus_3mom.Phi() - TMath::Pi() );
  _data.recoilmass = ( CM_4mom - mupair_4mom ).M();

  //_data.mumu_mass_mc = ( muplus_4mom_mc + muminus_4mom_mc ).M();
  //_data.mumu_E_mc = ( muplus_4mom_mc + muminus_4mom_mc ).E();
  //_data.mumu_costh_mc = ( muplus_4mom_mc.Vect().Unit() ).Dot( muminus_4mom_mc.Vect().Unit() );

  //weight of muon event by event computed from covariance matrix in momenta space
  //beware!!! [0]:Px, [1]:Py, [2]:Pz, [3]:E
  float term1st = 0, term2nd = 0;
  term1st = - muplus_4mom[0] * ( muplus_4mom[0]*muminus_covMat0 + muplus_4mom[1]*muminus_covMat1
			       + muplus_4mom[2]*muminus_covMat3 + muplus_4mom[3]*muminus_covMat6 )
            - muplus_4mom[1] * ( muplus_4mom[0]*muminus_covMat1 + muplus_4mom[1]*muminus_covMat2
			       + muplus_4mom[2]*muminus_covMat4 + muplus_4mom[3]*muminus_covMat7 )
            - muplus_4mom[2] * ( muplus_4mom[0]*muminus_covMat3 + muplus_4mom[1]*muminus_covMat4
			       + muplus_4mom[2]*muminus_covMat5 + muplus_4mom[3]*muminus_covMat8 )
            + muplus_4mom[3] * ( muplus_4mom[0]*muminus_covMat6 + muplus_4mom[1]*muminus_covMat7
			       + muplus_4mom[2]*muminus_covMat8 + muplus_4mom[3]*muminus_covMat9 );
  term2nd = - muminus_4mom[0] * ( muminus_4mom[0]*muplus_covMat0 + muminus_4mom[1]*muplus_covMat1
				+ muminus_4mom[2]*muplus_covMat3 + muminus_4mom[3]*muplus_covMat6 )
            - muminus_4mom[1] * ( muminus_4mom[0]*muplus_covMat1 + muminus_4mom[1]*muplus_covMat2
			        + muminus_4mom[2]*muplus_covMat4 + muminus_4mom[3]*muplus_covMat7 )
            - muminus_4mom[2] * ( muminus_4mom[0]*muplus_covMat3 + muminus_4mom[1]*muplus_covMat4
			        + muminus_4mom[2]*muplus_covMat5 + muminus_4mom[3]*muplus_covMat8 )
            + muminus_4mom[3] * ( muminus_4mom[0]*muplus_covMat6 + muminus_4mom[1]*muplus_covMat7
			        + muminus_4mom[2]*muplus_covMat8 + muminus_4mom[3]*muplus_covMat9 );
  if( _data.mumu_mass > 0 ){
    _data.sigma_mumu_mass = ( 1. / _data.mumu_mass ) * TMath::Sqrt( term1st + term2nd );
  }

  //vertex information
  if( muminus_track != 0 && muplus_track != 0 ){
    VertexInfo vv;
    //vv.setBField(3.5);
    vv.addTrack( muminus_track );
    vv.addTrack( muplus_track  );
    //beamspot constraint
    float xyz[3] = { 10.e-3, 10.e-6, 0.3 }; //mm
    vv.setBeamSpotSize( xyz );
    vv.useIPcon( true );

    TVector3 vtxPos = vv.getVertexPosition();
    _data.vtxPos_x = vtxPos[0];
    _data.vtxPos_y = vtxPos[1];
    _data.vtxPos_z = vtxPos[2];
    _data.vtxChisq = vv.getVertexChisq();
    //cout << "vtxPos_z = " << vtxPos[2] << endl;
    //cout << "chi2 = " << vv.getVertexChisq() << endl;
  }

  //momentum resolution sigma{p}
  float muminus_momres = 0, muplus_momres = 0;
  float muminus_bottom = muminus_4mom[0]*muminus_4mom[0] + muminus_4mom[1]*muminus_4mom[1] + muminus_4mom[2]*muminus_4mom[2];
  float muminus_top1st = muminus_4mom[0]*muminus_4mom[0]*muminus_covMat0
                       + muminus_4mom[1]*muminus_4mom[1]*muminus_covMat2
                       + muminus_4mom[2]*muminus_4mom[2]*muminus_covMat5;
  float muminus_top2nd = muminus_4mom[0]*muminus_4mom[1]*muminus_covMat1
                       + muminus_4mom[0]*muminus_4mom[2]*muminus_covMat3
                       + muminus_4mom[1]*muminus_4mom[2]*muminus_covMat4;
  muminus_momres = 0.5 * TMath::Sqrt( ( muminus_top1st + 2*muminus_top2nd ) / muminus_bottom );
  _data.muminus_momres = muminus_momres;
  float muplus_bottom = muplus_4mom[0]*muplus_4mom[0] + muplus_4mom[1]*muplus_4mom[1] + muplus_4mom[2]*muplus_4mom[2];
  float muplus_top1st = muplus_4mom[0]*muplus_4mom[0]*muplus_covMat0
                      + muplus_4mom[1]*muplus_4mom[1]*muplus_covMat2
                      + muplus_4mom[2]*muplus_4mom[2]*muplus_covMat5;
  float muplus_top2nd = muplus_4mom[0]*muplus_4mom[1]*muplus_covMat1
                      + muplus_4mom[0]*muplus_4mom[2]*muplus_covMat3
                      + muplus_4mom[1]*muplus_4mom[2]*muplus_covMat4;
  muplus_momres = 0.5 * TMath::Sqrt( ( muplus_top1st + 2*muplus_top2nd ) / muplus_bottom );
  _data.muplus_momres = muplus_momres;

  //pseudomass method used at OPAL
  TVector3 boostVec = -CM_4mom.BoostVector();
  TLorentzVector boost_muminus_4mom = muminus_4mom;
  TLorentzVector boost_muplus_4mom  = muplus_4mom;
  boost_muminus_4mom.Boost( boostVec );
  boost_muplus_4mom.Boost( boostVec );

  TVector3 l_3mom = boost_muminus_4mom.Vect();
  TVector3 L_3mom = boost_muplus_4mom.Vect();
  float    l_E    = boost_muminus_4mom.E();
  float    L_E    = boost_muplus_4mom.E();
  float    mlep   = boost_muminus_4mom.M();
  float    mLEP   = boost_muplus_4mom.M();
  float    beam   = CM / 2.;
  float    l2     = l_3mom.Mag2();
  float    L2     = L_3mom.Mag2();
  float    sum2   = (l_3mom + L_3mom).Mag2();
  float    lL     = l_3mom.Dot( L_3mom );
  TVector3 cross  = L_3mom.Cross( l_3mom );
  float    cross2 = cross.Mag2();

  float P = beam*L_E - L_E*L_E + 0.5*mLEP*mLEP;
  float Q = -beam*l_E - lL + 0.5*mlep*mlep;
  float in = sum2*(beam-L_E)*(beam-L_E) - (P+Q)*(P+Q);
  float M2 = 2 * ( P*l2 - Q*L2 + (P-Q)*lL + TMath::Sqrt(cross2*in) ) / ( sum2 );
  _data.pseudomass = TMath::Sqrt(M2);


  //*************************************
  //loop of after isolated lepton tagging
  //*************************************
  int n_pfoswo = 0, n_highPt = 0;
  int n_highPt_neutral = 0, n_highPt_charged = 0;
  if( PFOsWithoutIsoleps != 0 ){
    n_pfoswo = PFOsWithoutIsoleps->getNumberOfElements();
    for( int i = 0; i < n_pfoswo; i++ ){
      ReconstructedParticle* pfo_wo = dynamic_cast< ReconstructedParticle* >( PFOsWithoutIsoleps->getElementAt(i) );
      if( pfo_wo->getCharge() != 0 ){
	TVector3 mom = TVector3( pfo_wo->getMomentum() );
        float pfoPt = TMath::Sqrt(  mom[0]*mom[0] + mom[1]*mom[1] );
        if( pfoPt > _energetic ){
	  n_highPt++;
	  _data.pfowo_Pt[i] = pfoPt;
	  _data.pfowo_costh[i] = TVector3(0,0,1).Dot( mom.Unit() );
	  if( pfo_wo->getCharge() == 0 ) n_highPt_neutral++;
	  else n_highPt_charged++;
	}
      }

    }
  }
  _data.n_pfoswo = n_pfoswo;
  _data.n_highPt = n_highPt;
  _data.n_highPt_neutral = n_highPt_neutral;
  _data.n_highPt_charged = n_highPt_charged;


  //****************
  //read process map
  //****************
  std::string sprocess = evt->parameters().getStringVal( "Process" );
  float xsection = evt->parameters().getFloatVal( "crossSection" );
  _data.xsec = xsection;
  _data.processid = evt->parameters().getIntVal( "ProcessID" );
  _data.epol = evt->parameters().getFloatVal( "Pol0" );
  _data.ppol = evt->parameters().getFloatVal( "Pol1" );


    //cout << "before: " << sprocess << endl;


    //cout << "after: " << sprocess << endl;

    /*
    if ( _mapProcess.find( sprocess ) == _mapProcess.end() ) {
      streamlog_out(ERROR) << "Error: process not found in the process map!" << std::endl; 
    }
    else {
      TMapProcess::iterator it;
      for ( it  = _mapProcess.lower_bound( sprocess );
            it != _mapProcess.upper_bound( sprocess );
            it++ ) {
	
        float xsproc = it->second.second;
        if (fabs( ( xsproc - xsection ) / xsproc ) < 0.0001 ) {
           _data.processid = it->second.first; 
           break;
        }
      } 
    }
    */


  // Fill the event data 
  _data.evt = _nEvt;
  _dataTree->Fill();
  delete _navpfo;

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() 
		       << std::endl ;
  
}


void HiggsToMuMuProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void HiggsToMuMuProcessor::end(){ 
  
//   std::cout << "HiggsToMuMuProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

  _otfile->cd();
  _dataTree->Write();
  _otfile->Write();
  _otfile->Close();

  std::cout << "END OF PROCESSING" << std::endl;

}

void HiggsToMuMuProcessor::makeNTuple() {

  //Output root file
  _otfile = new TFile( _rootfilename.c_str() , "RECREATE" );
  RECO_DATA &d = _data;

  //Define root tree
  _dataTree  = new TTree( "dataTree" , "events" );
  //event number
  _dataTree->Branch( "evt", &d.evt, "evt/I" );

  //inserted by original slcio file
  _dataTree->Branch( "processid", &d.processid, "processid/I" );
  _dataTree->Branch( "xsec"     , &d.xsec     , "xsec"        );
  _dataTree->Branch( "epol"     , &d.epol     , "epol"        );
  _dataTree->Branch( "ppol"     , &d.ppol     , "ppol"        );

  //read out the thrust and sphere information
  _dataTree->Branch( "principalthrust" , &d.principalthrust , "principalthrust"  );
  _dataTree->Branch( "majorthrust"     , &d.majorthrust     , "majorthrust"      );
  _dataTree->Branch( "minorthrust"     , &d.minorthrust     , "minorthrust"      );
  _dataTree->Branch( "costh_thrustaxis", &d.costh_thrustaxis, "costh_thrustaxis" );
  //_dataTree->Branch( "thrustAxisX"     , &d.thrustAxisX     , "thrustAxisX"      );
  //_dataTree->Branch( "thrustAxisY"     , &d.thrustAxisY     , "thrustAxisY"      );
  //_dataTree->Branch( "thrustAxisZ"     , &d.thrustAxisZ     , "thrustAxisZ"      );
  _dataTree->Branch( "sphericity"      , &d.sphericity      , "sphericity"       );
  _dataTree->Branch( "aplanarity"      , &d.aplanarity      , "aplanarity"       );

  //number of PFOs
  _dataTree->Branch( "npfos", &d.npfos, "npfos/I" );

  //MC loop part
  _dataTree->Branch( "ISR1_E"    , &d.ISR1_E    , "ISR1_E"     );
  _dataTree->Branch( "ISR2_E"    , &d.ISR2_E    , "ISR2_E"     );
  _dataTree->Branch( "ISR1_costh", &d.ISR1_costh, "ISR1_costh" );
  _dataTree->Branch( "ISR2_costh", &d.ISR2_costh, "ISR2_costh" );

  _dataTree->Branch( "type"         , &d.type         , "type/I"        );
  _dataTree->Branch( "tau_m"        , &d.tau_m        , "tau_m/I"       );
  _dataTree->Branch( "higgsdecay1"  , &d.higgsdecay1  , "higgsdecay1/I" );
  _dataTree->Branch( "higgsdecay2"  , &d.higgsdecay2  , "higgsdecay2/I" );
  _dataTree->Branch( "MC_mumu_mass" , &d.MC_mumu_mass , "MC_mumu_mass"  );
  _dataTree->Branch( "MC_mumu_costh", &d.MC_mumu_costh, "MC_mumu_costh" );

  _dataTree->Branch( "h_nt_m", &d.h_nt_m, "h_nt_m/I" );

  _dataTree->Branch( "h_nq"      , &d.h_nq      , "h_nq/I"       );
  _dataTree->Branch( "h_ne"      , &d.h_ne      , "h_ne/I"       );
  _dataTree->Branch( "h_nm"      , &d.h_nm      , "h_nm/I"       );
  _dataTree->Branch( "h_nt"      , &d.h_nt      , "h_nt/I"       );
  _dataTree->Branch( "h_nn"      , &d.h_nn      , "h_nn/I"       );
  _dataTree->Branch( "h_nboson_m", &d.h_nboson_m, "h_nboson_m/I" );

  _dataTree->Branch( "forff1"   , &d.forff1   , "forff1/I"  );
  _dataTree->Branch( "forff2"   , &d.forff2   , "forff2/I"  );
  _dataTree->Branch( "ff1_costh", &d.ff1_costh, "ff1_costh" );
  _dataTree->Branch( "ff2_costh", &d.ff2_costh, "ff2_costh" );
  _dataTree->Branch( "massff"   , &d.massff   , "massff"    );
  _dataTree->Branch( "ff_nt_m"  , &d.ff_nt_m  , "ff_nt_m/I" );

  _dataTree->Branch( "forffff1"     , &d.forffff1     , "forffff1/I"    );
  _dataTree->Branch( "forffff2"     , &d.forffff2     , "forffff2/I"    );
  _dataTree->Branch( "forffff3"     , &d.forffff3     , "forffff3/I"    );
  _dataTree->Branch( "forffff4"     , &d.forffff4     , "forffff4/I"    );
  _dataTree->Branch( "mass1"        , &d.mass1        , "mass1"         );
  _dataTree->Branch( "mass2"        , &d.mass2        , "mass2"         );
  _dataTree->Branch( "nq"           , &d.nq           , "nq/I"          );
  _dataTree->Branch( "ne"           , &d.ne           , "ne/I"          );
  _dataTree->Branch( "nm"           , &d.nm           , "nm/I"          );
  _dataTree->Branch( "nt"           , &d.nt           , "nt/I"          );
  _dataTree->Branch( "nn"           , &d.nn           , "nn/I"          );
  _dataTree->Branch( "no"           , &d.no           , "no/I"          );
  _dataTree->Branch( "nt_m"         , &d.nt_m         , "nt_m/I"        );
  _dataTree->Branch( "ffff1_costh"  , &d.ffff1_costh  , "ffff1_costh"   );
  _dataTree->Branch( "ffff2_costh"  , &d.ffff2_costh  , "ffff2_costh"   );
  _dataTree->Branch( "ffff3_costh"  , &d.ffff3_costh  , "ffff3_costh"   );
  _dataTree->Branch( "ffff4_costh"  , &d.ffff4_costh  , "ffff4_costh"   );
  _dataTree->Branch( "ffff1_t_costh", &d.ffff1_t_costh, "ffff1_t_costh" );
  _dataTree->Branch( "ffff2_t_costh", &d.ffff2_t_costh, "ffff2_t_costh" );
  _dataTree->Branch( "ffff3_t_costh", &d.ffff3_t_costh, "ffff3_t_costh" );
  _dataTree->Branch( "ffff4_t_costh", &d.ffff4_t_costh, "ffff4_t_costh" );

  _dataTree->Branch( "MC_muminus_parent", &d.MC_muminus_parent, "MC_muminus_parent/I" );
  _dataTree->Branch( "MC_muplus_parent" , &d.MC_muplus_parent , "MC_muplus_parent/I"  );
  _dataTree->Branch( "MC_mumupair_mass" , &d.MC_mumupair_mass , "MC_mumupair_mass"    );

  //_dataTree->Branch( "MC_Z_mumu_mass", &d.MC_Z_mumu_mass, "MC_Z_mumu_mass" );

  _dataTree->Branch( "MC_mu_E", &d.MC_mu_E, "MC_mu_E" );

  _dataTree->Branch( "ntaumode", &d.ntaumode, "ntaumode/I"          );
  _dataTree->Branch( "taumode" , &d.taumode , "taumode[ntaumode]/I" );

  //check the relation between MCParticle and PFO
  _dataTree->Branch( "overlay_mc"        , &d.overlay_mc        , "overlay_mc/I"       );
  _dataTree->Branch( "esum_mc"           , &d.esum_mc           , "esum_mc"            );
  _dataTree->Branch( "esum_charged_mc"   , &d.esum_charged_mc   , "esum_charged_mc"    );
  _dataTree->Branch( "esum_neutral_mc"   , &d.esum_neutral_mc   , "esum_neutral_mc"    );
  _dataTree->Branch( "overlay_E_mc"      , &d.overlay_E_mc      , "overlay_E_mc"       );
  _dataTree->Branch( "overlay_charged_mc", &d.overlay_charged_mc, "overlay_charged_mc" );
  _dataTree->Branch( "overlay_neutral_mc", &d.overlay_neutral_mc, "overlay_neutral_mc" );
  _dataTree->Branch( "costh_missmom_mc"  , &d.costh_missmom_mc  , "costh_missmom_mc"   );
  _dataTree->Branch( "Pt_all_mc"         , &d.Pt_all_mc         , "Pt_all_mc"          );

  _dataTree->Branch( "n_highPt_mc", &d.n_highPt_mc, "n_highPt_mc/I"   );
  _dataTree->Branch( "Pt_mc"      , &d.Pt_mc      , "Pt_mc[npfos]"    );
  _dataTree->Branch( "costh_mc"   , &d.costh_mc   , "costh_mc[npfos]" );
  _dataTree->Branch( "pfoPt"      , &d.pfoPt      , "pfoPt[npfos]"    );
  _dataTree->Branch( "pfocosth"   , &d.pfocosth   , "pfocosth[npfos]" );

  _dataTree->Branch( "mumu_mass_mc" , &d.mumu_mass_mc , "mumu_mass_mc"   );
  _dataTree->Branch( "mumu_costh_mc", &d.mumu_costh_mc, "mumu_costh_mc"  );
  _dataTree->Branch( "n_muminus_mc" , &d.n_muminus_mc , "n_muminus_mc/I" );
  _dataTree->Branch( "n_muplus_mc"  , &d.n_muplus_mc  , "n_muplus_mc/I"  );

  _dataTree->Branch( "mumu_mass_mc_smear_1_3", &d.mumu_mass_mc_smear_1_3, "mumu_mass_mc_smear_1_3" );
  _dataTree->Branch( "mumu_mass_mc_smear_5_4", &d.mumu_mass_mc_smear_5_4, "mumu_mass_mc_smear_5_4" );
  _dataTree->Branch( "mumu_mass_mc_smear_3_4", &d.mumu_mass_mc_smear_3_4, "mumu_mass_mc_smear_3_4" );
  _dataTree->Branch( "mumu_mass_mc_smear_2_4", &d.mumu_mass_mc_smear_2_4, "mumu_mass_mc_smear_2_4" );
  _dataTree->Branch( "mumu_mass_mc_smear_1_4", &d.mumu_mass_mc_smear_1_4, "mumu_mass_mc_smear_1_4" );
  _dataTree->Branch( "mumu_mass_mc_smear_5_5", &d.mumu_mass_mc_smear_5_5, "mumu_mass_mc_smear_5_5" );
  _dataTree->Branch( "mumu_mass_mc_smear_3_5", &d.mumu_mass_mc_smear_3_5, "mumu_mass_mc_smear_3_5" );
  _dataTree->Branch( "mumu_mass_mc_smear_2_5", &d.mumu_mass_mc_smear_2_5, "mumu_mass_mc_smear_2_5" );
  _dataTree->Branch( "mumu_mass_mc_smear_1_5", &d.mumu_mass_mc_smear_1_5, "mumu_mass_mc_smear_1_5" );
  _dataTree->Branch( "mumu_mass_mc_smear_5_6", &d.mumu_mass_mc_smear_5_6, "mumu_mass_mc_smear_5_6" );
  _dataTree->Branch( "mumu_mass_mc_smear_3_6", &d.mumu_mass_mc_smear_3_6, "mumu_mass_mc_smear_3_6" );
  _dataTree->Branch( "mumu_mass_mc_smear_2_6", &d.mumu_mass_mc_smear_2_6, "mumu_mass_mc_smear_2_6" );
  _dataTree->Branch( "mumu_mass_mc_smear_1_6", &d.mumu_mass_mc_smear_1_6, "mumu_mass_mc_smear_1_6" );

  //PFO loop part
  _dataTree->Branch( "ntrks"        , &d.ntrks        , "ntrks/I"       );
  _dataTree->Branch( "vtx_x"        , &d.vtx_x        , "vtx_x"         );
  _dataTree->Branch( "vtx_y"        , &d.vtx_y        , "vtx_y"         );
  _dataTree->Branch( "vtx_z"        , &d.vtx_z        , "vtx_z"         );
  _dataTree->Branch( "esum"         , &d.esum         , "esum"          );
  _dataTree->Branch( "psum_x"       , &d.psum_x       , "psum_x"        );
  _dataTree->Branch( "psum_y"       , &d.psum_y       , "psum_y"        );
  _dataTree->Branch( "psum_z"       , &d.psum_z       , "psum_z"        );
  _dataTree->Branch( "maxPt"        , &d.maxPt        , "maxPt"         );
  _dataTree->Branch( "Pt"           , &d.Pt           , "Pt"            );
  _dataTree->Branch( "Pt_all"       , &d.Pt_all       , "Pt_all"        );
  _dataTree->Branch( "mass_vis"     , &d.mass_vis     , "mass_vis"      );
  _dataTree->Branch( "misse"        , &d.misse        , "misse"         );
  _dataTree->Branch( "missmom_x"    , &d.missmom_x    , "missmom_x"     );
  _dataTree->Branch( "missmom_y"    , &d.missmom_y    , "missmom_y"     );
  _dataTree->Branch( "missmom_z"    , &d.missmom_z    , "missmom_z"     );
  _dataTree->Branch( "missmom_mag"  , &d.missmom_mag  , "missmom_mag"   );
  _dataTree->Branch( "costh_missmom", &d.costh_missmom, "costh_missmom" );
  _dataTree->Branch( "mass_inv"     , &d.mass_inv     , "mass_inv"      );

  _dataTree->Branch( "pfo_e"      , &d.pfo_e      , "pfo_e[npfos]"       );
  _dataTree->Branch( "pfo_px"     , &d.pfo_px     , "pfo_px[npfos]"      );
  _dataTree->Branch( "pfo_py"     , &d.pfo_py     , "pfo_py[npfos]"      );
  _dataTree->Branch( "pfo_pz"     , &d.pfo_pz     , "pfo_pz[npfos]"      );
  _dataTree->Branch( "pfo_mom_mag", &d.pfo_mom_mag, "pfo_mom_mag[npfos]" );
  _dataTree->Branch( "pfo_chrg"   , &d.pfo_chrg   , "pfo_chrg[npfos]"    );
  _dataTree->Branch( "pfo_pdg"    , &d.pfo_pdg    , "pfo_pdg[npfos]/I"   );
  _dataTree->Branch( "pfo_costh"  , &d.pfo_costh  , "pfo_costh[npfos]"   );
  _dataTree->Branch( "pfo_Pt"     , &d.pfo_Pt     , "pfo_Pt[npfos]"      );
  _dataTree->Branch( "pfo_ntrk"   , &d.pfo_ntrk   , "pfo_ntrk[npfos]/I"  );
  _dataTree->Branch( "pfo_nclus"  , &d.pfo_nclus  , "pfo_nclus[npfos]/I" );

  _dataTree->Branch( "pfo_d0"       , &d.pfo_d0       , "pfo_d0[npfos]"        );
  _dataTree->Branch( "pfo_z0"       , &d.pfo_z0       , "pfo_z0[npfos]"        );
  _dataTree->Branch( "pfo_phi"      , &d.pfo_phi      , "pfo_phi[npfos]"       );
  _dataTree->Branch( "pfo_omega"    , &d.pfo_omega    , "pfo_omega[npfos]"     );
  _dataTree->Branch( "pfo_tanlambda", &d.pfo_tanlambda, "pfo_tanlambda[npfos]" );
  _dataTree->Branch( "pfo_d0sig"    , &d.pfo_d0sig    , "pfo_d0sig[npfos]"     );
  _dataTree->Branch( "pfo_z0sig"    , &d.pfo_z0sig    , "pfo_z0sig[npfos]"     );
  _dataTree->Branch( "pfo_Chi2"     , &d.pfo_Chi2     , "pfo_Chi2[npfos]"      );
  _dataTree->Branch( "pfo_Ndf"      , &d.pfo_Ndf      , "pfo_Ndf[npfos]/I"     );
  _dataTree->Branch( "pfo_Chi2Ndf"  , &d.pfo_Chi2Ndf  , "pfo_Chi2Ndf[npfos]"   );
  _dataTree->Branch( "pfo_dEdx"     , &d.pfo_dEdx     , "pfo_dEdx[npfos]"      );
  _dataTree->Branch( "pfo_dEdxError", &d.pfo_dEdxError, "pfo_dEdxError[npfos]" );
  _dataTree->Branch( "pfo_RadiusOfInnermostHit", &d.pfo_RadiusOfInnermostHit, "pfo_RadiusOfInnermostHit[npfos]" );

  _dataTree->Branch( "pfo_VXD", &d.pfo_VXD, "pfo_VXD[npfos]/I" );
  _dataTree->Branch( "pfo_FTD", &d.pfo_FTD, "pfo_FTD[npfos]/I" );
  _dataTree->Branch( "pfo_SIT", &d.pfo_SIT, "pfo_SIT[npfos]/I" );
  _dataTree->Branch( "pfo_TPC", &d.pfo_TPC, "pfo_TPC[npfos]/I" );
  _dataTree->Branch( "pfo_SET", &d.pfo_SET, "pfo_SET[npfos]/I" );
  _dataTree->Branch( "pfo_ETD", &d.pfo_ETD, "pfo_ETD[npfos]/I" );

  _dataTree->Branch( "pfo_ecal"    , &d.pfo_ecal    , "pfo_ecal[npfos]"     );
  _dataTree->Branch( "pfo_hcal"    , &d.pfo_hcal    , "pfo_hcal[npfos]"     );
  _dataTree->Branch( "pfo_yoke"    , &d.pfo_yoke    , "pfo_yoke[npfos]"     );
  _dataTree->Branch( "pfo_lcal"    , &d.pfo_lcal    , "pfo_lcal[npfos]"     );
  _dataTree->Branch( "pfo_lhcal"   , &d.pfo_lhcal   , "pfo_lhcal[npfos]"    );
  _dataTree->Branch( "pfo_bcal"    , &d.pfo_bcal    , "pfo_bcal[npfos]"     );
  _dataTree->Branch( "pfo_ecalfrac", &d.pfo_ecalfrac, "pfo_ecalfrac[npfos]" );
  _dataTree->Branch( "pfo_EbyP"    , &d.pfo_EbyP    , "pfo_EbyP[npfos]"     );

  _dataTree->Branch( "pfo_pid_mc"      , &d.pfo_pid_mc      , "pfo_pid_mc[npfos]/I"       );
  _dataTree->Branch( "pfo_E_mc"        , &d.pfo_E_mc        , "pfo_E_mc[npfos]"           );
  _dataTree->Branch( "pfo_Px_mc"       , &d.pfo_Px_mc       , "pfo_Px_mc[npfos]"          );
  _dataTree->Branch( "pfo_Py_mc"       , &d.pfo_Py_mc       , "pfo_Py_mc[npfos]"          );
  _dataTree->Branch( "pfo_Pz_mc"       , &d.pfo_Pz_mc       , "pfo_Pz_mc[npfos]"          );
  //_dataTree->Branch( "pfo_taudecaymode", &d.pfo_taudecaymode, "pfo_taudecaymode[npfos]/I" );
  //_dataTree->Branch( "count_overlay"   , &d.count_overlay   , "count_overlay/I"           );
  //_dataTree->Branch( "charged_overlay" , &d.charged_overlay , "charged_overlay/I"         );

  _dataTree->Branch( "parent_pid_mc" , &d.parent_pid_mc , "parent_pid_mc[npfos]/I"  );
  _dataTree->Branch( "parent_E_mc"   , &d.parent_E_mc   , "parent_E_mc[npfos]"      );
  _dataTree->Branch( "parent_Px_mc"  , &d.parent_Px_mc  , "parent_Px_mc[npfos]"     );
  _dataTree->Branch( "parent_Py_mc"  , &d.parent_Py_mc  , "parent_Py_mc[npfos]"     );
  _dataTree->Branch( "parent_Pz_mc"  , &d.parent_Pz_mc  , "parent_Pz_mc[npfos]"     );
  _dataTree->Branch( "pfo_mother_pid", &d.pfo_mother_pid, "pfo_mother_pid[npfos]/I" );

  //loop of isolated lepton finder
  _dataTree->Branch( "n_isoleps"           , &d.n_isoleps           , "n_isoleps/I"          );
  _dataTree->Branch( "n_muplus"            , &d.n_muplus            , "n_muplus/I"           );
  _dataTree->Branch( "muplus_E"            , &d.muplus_E            , "muplus_E"             );
  _dataTree->Branch( "muplus_costh"        , &d.muplus_costh        , "muplus_costh"         );
  _dataTree->Branch( "muplus_azi"          , &d.muplus_azi          , "muplus_azi"           );
  _dataTree->Branch( "muplus_charge_costh" , &d.muplus_charge_costh , "muplus_charge_costh"  );
  _dataTree->Branch( "muplus_Pt"           , &d.muplus_Pt           , "muplus_Pt"            );
  _dataTree->Branch( "muplus_mom_mag"      , &d.muplus_mom_mag      , "muplus_mom_mag"       );
  _dataTree->Branch( "n_muminus"           , &d.n_muminus           , "n_muminus/I"          );
  _dataTree->Branch( "muminus_E"           , &d.muminus_E           , "muminus_E"            );
  _dataTree->Branch( "muminus_costh"       , &d.muminus_costh       , "muminus_costh"        );
  _dataTree->Branch( "muminus_charge_costh", &d.muminus_charge_costh, "muminus_charge_costh" );
  _dataTree->Branch( "sum_charge_costh"    , &d.sum_charge_costh    , "sum_charge_costh"     );
  _dataTree->Branch( "muminus_azi"         , &d.muminus_azi         , "muminus_azi"          );
  _dataTree->Branch( "muminus_Pt"          , &d.muminus_Pt          , "muminus_Pt"           );
  _dataTree->Branch( "muminus_mom_mag"     , &d.muminus_mom_mag     , "muminus_mom_mag"      );
  _dataTree->Branch( "leadingmu_E"         , &d.leadingmu_E         , "leadingmu_E"          );
  _dataTree->Branch( "subleadingmu_E"      , &d.subleadingmu_E      , "subleadingmu_E"       );
  _dataTree->Branch( "leadingmu_costh"     , &d.leadingmu_costh     , "leadingmu_costh"      );
  _dataTree->Branch( "subleadingmu_costh"  , &d.subleadingmu_costh  , "subleadingmu_costh"   );

  _dataTree->Branch( "muplus_d0"       , &d.muplus_d0       , "muplus_d0"        );
  _dataTree->Branch( "muplus_z0"       , &d.muplus_z0       , "muplus_z0"        );
  _dataTree->Branch( "muplus_d0sig"    , &d.muplus_d0sig    , "muplus_d0sig"     );
  _dataTree->Branch( "muplus_z0sig"    , &d.muplus_z0sig    , "muplus_z0sig"     );
  _dataTree->Branch( "muplus_phi"      , &d.muplus_phi      , "muplus_phi"       );
  _dataTree->Branch( "muplus_omega"    , &d.muplus_omega    , "muplus_omega"     );
  _dataTree->Branch( "muplus_tanlambda", &d.muplus_tanlambda, "muplus_tanlambda" );
  _dataTree->Branch( "muplus_Chi2"     , &d.muplus_Chi2     , "muplus_Chi2"      );
  _dataTree->Branch( "muplus_Ndf"      , &d.muplus_Ndf      , "muplus_Ndf/I"     );
  _dataTree->Branch( "muplus_Chi2Ndf"  , &d.muplus_Chi2Ndf  , "muplus_Chi2Ndf"   );
  _dataTree->Branch( "muplus_dEdx"     , &d.muplus_dEdx     , "muplus_dEdx"      );
  _dataTree->Branch( "muplus_dEdxError", &d.muplus_dEdxError, "muplus_dEdxError" );
  _dataTree->Branch( "muplus_RadiusOfInnermostHit", &d.muplus_RadiusOfInnermostHit, "muplus_RadiusOfInnermostHit" );
  _dataTree->Branch( "muminus_d0"       , &d.muminus_d0       , "muminus_d0"        );
  _dataTree->Branch( "muminus_z0"       , &d.muminus_z0       , "muminus_z0"        );
  _dataTree->Branch( "muminus_d0sig"    , &d.muminus_d0sig    , "muminus_d0sig"     );
  _dataTree->Branch( "muminus_z0sig"    , &d.muminus_z0sig    , "muminus_z0sig"     );
  _dataTree->Branch( "muminus_phi"      , &d.muminus_phi      , "muminus_phi"       );
  _dataTree->Branch( "muminus_omega"    , &d.muminus_omega    , "muminus_omega"     );
  _dataTree->Branch( "muminus_tanlambda", &d.muminus_tanlambda, "muminus_tanlambda" );
  _dataTree->Branch( "muminus_Chi2"     , &d.muminus_Chi2     , "muminus_Chi2"      );
  _dataTree->Branch( "muminus_Ndf"      , &d.muminus_Ndf      , "muminus_Ndf/I"     );
  _dataTree->Branch( "muminus_Chi2Ndf"  , &d.muminus_Chi2Ndf  , "muminus_Chi2Ndf"   );
  _dataTree->Branch( "muminus_dEdx"     , &d.muminus_dEdx     , "muminus_dEdx"      );
  _dataTree->Branch( "muminus_dEdxError", &d.muminus_dEdxError, "muminus_dEdxError" );
  _dataTree->Branch( "muminus_RadiusOfInnermostHit", &d.muminus_RadiusOfInnermostHit, "muminus_RadiusOfInnermostHit" );

  _dataTree->Branch( "muplus_ecal"    , &d.muplus_ecal    , "muplus_ecal"     );
  _dataTree->Branch( "muplus_hcal"    , &d.muplus_hcal    , "muplus_hcal"     );
  _dataTree->Branch( "muplus_yoke"    , &d.muplus_yoke    , "muplus_yoke"     );
  _dataTree->Branch( "muplus_lcal"    , &d.muplus_lcal    , "muplus_lcal"     );
  _dataTree->Branch( "muplus_lhcal"   , &d.muplus_lhcal   , "muplus_lhcal"    );
  _dataTree->Branch( "muplus_bcal"    , &d.muplus_bcal    , "muplus_bcal"     );
  _dataTree->Branch( "muplus_ecalfrac", &d.muplus_ecalfrac, "muplus_ecalfrac" );
  _dataTree->Branch( "muplus_EbyP"    , &d.muplus_EbyP    , "muplus_EbyP"     );
  _dataTree->Branch( "muplus_VXD"     , &d.muplus_VXD     , "muplus_VXD/I"    );
  _dataTree->Branch( "muplus_FTD"     , &d.muplus_FTD     , "muplus_FTD/I"    );
  _dataTree->Branch( "muplus_SIT"     , &d.muplus_SIT     , "muplus_SIT/I"    );
  _dataTree->Branch( "muplus_TPC"     , &d.muplus_TPC     , "muplus_TPC/I"    );
  _dataTree->Branch( "muplus_SET"     , &d.muplus_SET     , "muplus_SET/I"    );
  _dataTree->Branch( "muplus_ETD"     , &d.muplus_ETD     , "muplus_ETD/I"    );

  _dataTree->Branch( "muminus_ecal"    , &d.muminus_ecal    , "muminus_ecal"     );
  _dataTree->Branch( "muminus_hcal"    , &d.muminus_hcal    , "muminus_hcal"     );
  _dataTree->Branch( "muminus_yoke"    , &d.muminus_yoke    , "muminus_yoke"     );
  _dataTree->Branch( "muminus_lcal"    , &d.muminus_lcal    , "muminus_lcal"     );
  _dataTree->Branch( "muminus_lhcal"   , &d.muminus_lhcal   , "muminus_lhcal"    );
  _dataTree->Branch( "muminus_bcal"    , &d.muminus_bcal    , "muminus_bcal"     );
  _dataTree->Branch( "muminus_ecalfrac", &d.muminus_ecalfrac, "muminus_ecalfrac" );
  _dataTree->Branch( "muminus_EbyP"    , &d.muminus_EbyP    , "muminus_EbyP"     );
  _dataTree->Branch( "muminus_VXD"     , &d.muminus_VXD     , "muminus_VXD/I"    );
  _dataTree->Branch( "muminus_FTD"     , &d.muminus_FTD     , "muminus_FTD/I"    );
  _dataTree->Branch( "muminus_SIT"     , &d.muminus_SIT     , "muminus_SIT/I"    );
  _dataTree->Branch( "muminus_TPC"     , &d.muminus_TPC     , "muminus_TPC/I"    );
  _dataTree->Branch( "muminus_SET"     , &d.muminus_SET     , "muminus_SET/I"    );
  _dataTree->Branch( "muminus_ETD"     , &d.muminus_ETD     , "muminus_ETD/I"    );

  //MC matching of isolated lepton
  _dataTree->Branch( "mc_muplus_PDG"        , &d.mc_muplus_PDG        , "mc_muplus_PDG/I"         );
  _dataTree->Branch( "parent_mc_muplus_PDG" , &d.parent_mc_muplus_PDG , "parent_mc_muplus_PDG/I"  );
  _dataTree->Branch( "mc_muminus_PDG"       , &d.mc_muminus_PDG       , "mc_muminus_PDG/I"        );
  _dataTree->Branch( "parent_mc_muminus_PDG", &d.parent_mc_muminus_PDG, "parent_mc_muminus_PDG/I" );

  _dataTree->Branch( "mumu_E"           , &d.mumu_E           , "mumu_E"            );
  _dataTree->Branch( "mumu_mass"        , &d.mumu_mass        , "mumu_mass"         );
  _dataTree->Branch( "mumu_costh"       , &d.mumu_costh       , "mumu_costh"        );
  _dataTree->Branch( "mumu_acop"        , &d.mumu_acop        , "mumu_acop"         );
  _dataTree->Branch( "mumu_Pt"          , &d.mumu_Pt          , "mumu_Pt"           );
  _dataTree->Branch( "sigma_mumu_mass"  , &d.sigma_mumu_mass  , "sigma_mumu_mass"   );
  _dataTree->Branch( "recoilmass"       , &d.recoilmass       , "recoilmass"        );
  _dataTree->Branch( "mumu_costh_tobeam", &d.mumu_costh_tobeam, "mumu_costh_tobeam" );
  _dataTree->Branch( "mumu_mom_mag"     , &d.mumu_mom_mag     , "mumu_mom_mag"      );

  _dataTree->Branch( "vtxPos_x", &d.vtxPos_x, "vtxPos_x" );
  _dataTree->Branch( "vtxPos_y", &d.vtxPos_y, "vtxPos_y" );
  _dataTree->Branch( "vtxPos_z", &d.vtxPos_z, "vtxPos_z" );
  _dataTree->Branch( "vtxChisq", &d.vtxChisq, "vtxChisq" );

  _dataTree->Branch( "muminus_momres", &d.muminus_momres, "muminus_momres" );
  _dataTree->Branch( "muplus_momres" , &d.muplus_momres , "muplus_momres"  );
  _dataTree->Branch( "pseudomass", &d.pseudomass, "pseudomass" );

  //_dataTree->Branch( "mumu_mass_mc" , &d.mumu_mass_mc , "mumu_mass_mc"  );
  //_dataTree->Branch( "mumu_E_mc"    , &d.mumu_E_mc    , "mumu_E_mc"     );
  //_dataTree->Branch( "mumu_costh_mc", &d.mumu_costh_mc, "mumu_costh_mc" );

  //loop pf after isolated lepton tagging
  _dataTree->Branch( "n_pfoswo"   , &d.n_pfoswo   , "n_pfoswo/I"         );
  _dataTree->Branch( "pfowo_Pt"   , &d.pfowo_Pt   , "pfowo_Pt[npfos]"    );
  _dataTree->Branch( "pfowo_costh", &d.pfowo_costh, "pfowo_costh[npfos]" );
  _dataTree->Branch( "n_highPt"   , &d.n_highPt   , "n_highPt/I"         );
  _dataTree->Branch( "n_highPt_neutral", &d.n_highPt_neutral, "n_highPt_neutral/I" );
  _dataTree->Branch( "n_highPt_charged", &d.n_highPt_charged, "n_highPt_charged/I" );

  return;
}

void HiggsToMuMuProcessor::getCalEnergy( ReconstructedParticle* pfo , float* cale ) {
  float ecal = 0;
  float hcal = 0;
  std::vector<lcio::Cluster*> clusters = pfo->getClusters();
  for ( std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();
        iCluster!=clusters.end();
        ++iCluster) {
	ecal += (*iCluster)->getSubdetectorEnergies()[0];
	hcal += (*iCluster)->getSubdetectorEnergies()[1];
  }
  cale[0] = ecal;
  cale[1] = hcal;
}

/*
void HiggsToMuMuProcessor::readProcessList( const char* fname ) {
  std::ifstream in( fname );
  if( in.fail() ){
    _useProcessMap = false;
    std::cout << "Process file not found: " << fname << std::endl;
  }
  else{
    _useProcessMap = true;
    while( !in.eof() ){
      std::string process;
      std::string pol;
      int   processid;
      float xsection;

      //in >> processid >> process >> xsection ;
      in >> processid >> process >> pol >> xsection ; //procid500.txt specific
      //in >> process >> pol >> processid >> xsection ; //procid250.txt specific

      if ( process == "" ) break;

      _mapProcess.insert( TMapProcess::value_type 
                           ( process, std::pair< int, float >( processid, xsection ) ) ) ;

    }
  }
}
*/

int HiggsToMuMuProcessor::getTauDecayMode( const MCParticle *mcp ){
  /*
    return value is signed: >0 for tau-, <0 for tau+
    1: 1-prong electron
    2: 1-prong muon
    3: 1-prong hadron
    4: 1-prong hadron + neutral(s)
    5: 3-prong
    6: 3-prong + neutral(s)
    7: 5 or more prong (+ neutrals)
    10: error, charge not conserved
   */
  //look for tau
  if(!mcp) return 0;
  while( mcp && abs( mcp->getPDG() ) != 15 && mcp->getParents().size() > 0 ) mcp = mcp->getParents()[0];
  if( abs( mcp->getPDG() ) != 15 ) return 0;

  for(;;){
    MCParticle* daughterTau(0);
    vector< MCParticle* >::const_iterator iter;
    for( iter = mcp->getDaughters().begin(); !daughterTau && iter != mcp->getDaughters().end(); ++iter ){
      if( abs( (*iter)->getPDG() ) == 15 ) daughterTau = *iter;
    }

    if( daughterTau ){
      if( daughterTau == mcp ) break; //safeguard
      mcp = daughterTau;
    }
    else break;
  }

  //tau found
  const MCParticle *tau = mcp;
  int sign = tau->getPDG() / abs( tau->getPDG() );

  //count number of charged / neutral stable daughters
  int ne(0);
  int nmu(0);
  int nphoton(0);
  int npi0(0);
  int npi(0);
  int nK0(0);
  int nK(0);

  //take care of overlapping parent/daughter relationships e.g. W decay
  vector< const MCParticle* > added;

  stack< const MCParticle* > q;
  q.push( tau );
  while( !q.empty() ){
    const MCParticle *p = q.top();
    q.pop();

    switch( abs( p->getPDG() ) ){
    case 11:
      ++ne; break;
    case 13:
      ++nmu; break;
    case 22:
      ++nphoton; break;
    case 111:
      ++npi0; break;
    case 211:
      ++npi; break;
    case 311:
    case 310:
    case 130:
      ++nK0; break;
    case 321:
      ++nK; break;
    default:
      vector< MCParticle* >::const_iterator iter;

      for( iter = p->getDaughters().begin(); iter != p->getDaughters().end(); ++iter ){
	if( (*iter)->isCreatedInSimulation() == false ){
	  if( find( added.begin(), added.end(), *iter ) == added.end() ){
	    q.push( *iter );
	    added.push_back( *iter );
	  }
	}
      }
      break;
    }
  }

  int nch = ne + nmu + npi + nK;
  //int nne = nphoton + npi0 + nK0;
  int nne = npi0 + nK0;

  if( nch == 1 && ne  == 1 ) return sign * 1;
  if( nch == 1 && nmu == 1 ) return sign * 2;
  if( nch == 1 && nne == 0 ) return sign * 3;
  if( nch == 1 && nne >  0 ) return sign * 4;
  if( nch == 3 && nne == 0 ) return sign * 5;
  if( nch == 3 && nne >  0 ) return sign * 6;
  if( ( nch % 2 ) == 0 )     return sign * 10;
  return sign * 7;
}


MCParticle* HiggsToMuMuProcessor::findMCParticleFromIsolep( ReconstructedParticle* rp, LCCollection* AllPFOs ){
  if( rp == 0 ) return 0;
  /*
  if( AllPFOs == 0 ) return 0;

  MCParticle* mcp(0);

  TVector3 vec( rp->getMomentum() );
  ReconstructedParticle* match(0);
  double delta(1e6);

  int npfo = AllPFOs->getNumberOfElements();
  for( int i = 0; i < npfo; i++ ){
    ReconstructedParticle* pfo = dynamic_cast< ReconstructedParticle* >( AllPFOs->getElementAt(i) );
    TVector3 vecpfo( rp->getMomentum() );
    double tmp = ( vec - vecpfo ).Mag();
    if( tmp < delta ){
      delta = tmp;
      match = pfo;
    }
  }

  if( match && _navpfo->getRelatedToWeights( match ).size() > 0 ){
    mcp = dynamic_cast< MCParticle* >( _navpfo->getRelatedToObjects( match )[0] );
  }
  return mcp;
  */

  MCParticle* mcp(0);
  const ReconstructedParticleVec& rps = rp->getParticles();
  if( rps.size() == 0 ) return 0;
  ReconstructedParticle* pfo = rps[0];
  if( _navpfo->getRelatedToWeights( pfo ).size() > 0 ){
    mcp = dynamic_cast< MCParticle* >( _navpfo->getRelatedToObjects( pfo )[0] );
  }
  return mcp;
}
