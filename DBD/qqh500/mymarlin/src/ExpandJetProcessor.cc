#include "ExpandJetProcessor.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include <marlin/VerbosityLevels.h>

using namespace lcio;
using namespace marlin;

ExpandJetProcessor aExpandJetProcessor;

ExpandJetProcessor::ExpandJetProcessor()
  :Processor("ExpandJetProcessor"){

  //processor description
  _description = "";

  //register steering parameters
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputCollection",
			   "Input collection of reconstructed jets after applying kT clustering",
			   _inputCollection,
			   std::string("JetsAfterGamGamRemoval") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputCollection",
			    "Output collection of pfos after applying kT clustering",
			    _outputCollection,
			    std::string("PFOsAfterGamGamRemoval") );
}

void ExpandJetProcessor::init(){
  printParameters();
}

void ExpandJetProcessor::processRunHeader( LCRunHeader* run ){
}

void ExpandJetProcessor::processEvent( LCEvent* evt ){
  LCCollection* jetCol = 0;
  try{
    jetCol = evt->getCollection( _inputCollection );
  }
  catch(...){
    std::cout << "no events in JetsAfterGamGamRemoval" << std::endl;
  }
  LCCollectionVec* outCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  outCol->setSubset( true );

  if( jetCol != 0){
    unsigned int njet = jetCol->getNumberOfElements();
    for( unsigned int i = 0; i < njet; i++ ){
      ReconstructedParticle* jet = dynamic_cast< ReconstructedParticle* >( jetCol->getElementAt(i) );
      const ReconstructedParticleVec& pfovec = jet->getParticles();
      for( unsigned int j = 0; j < pfovec.size(); j++ ){
	outCol->addElement( pfovec[j] );
      }
    }
  }

  evt->addCollection( outCol, _outputCollection.c_str() );
}

void ExpandJetProcessor::check( LCEvent* evt ){
}

void ExpandJetProcessor::end(){
}
