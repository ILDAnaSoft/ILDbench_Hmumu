#include "ISRFinder.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include <TVector3.h>

using namespace lcio;
using namespace marlin;

ISRFinder aISRFinder;

ISRFinder::ISRFinder()
  :Processor("ISRFinder"){
  //processor description
  _description = "";

  //register steering parameters
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputCollection",
			   "Input collection of PFOs w/o isolated leptons",
			   _inputCollection,
			   std::string("PFOsWithoutIsoleps") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputCollection",
			    "Output collection after removing ISR",
			    _outputCollection,
			    std::string("PFOsWithoutISR") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "ISRCollection",
			    "Output collection for ISR",
			    _isrCollection,
			    std::string("ISR") );

  registerProcessorParameter( "Energy",
			      "energy of PFO",
			      _energy,
			      float(10) );

  registerProcessorParameter( "ConeCosTheta",
			      "definition of cone for isolation",
			      _coneCosTheta,
			      float(0.95) );

  registerProcessorParameter( "Ratio",
			      "pfo_E divided by coneE for isolation",
			      _ratio,
			      float(0.05) );
}

void ISRFinder::init(){
}

void ISRFinder::processRunHeader( LCRunHeader *run){
}

void ISRFinder::processEvent( LCEvent *evt ){
  LCCollection *inputCol = 0;
  try{
    inputCol = evt->getCollection( _inputCollection );
  }
  catch(...){
    std::cout << "NO inputs" << std::endl;
  }
  LCCollectionVec *outputCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  outputCol->setSubset( true );
  LCCollectionVec *isrCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  isrCol->setSubset( true );

  if( inputCol != 0 ){
    int npfo = inputCol->getNumberOfElements();
    for( int i = 0; i < npfo; i++ ){
      ReconstructedParticle *pfo = dynamic_cast< ReconstructedParticle* >( inputCol->getElementAt(i) );
      if( pfo->getType() != 22 ) outputCol->addElement( pfo );

      if( pfo->getType() == 22 ){
	TVector3 pfo_3mom = TVector3( pfo->getMomentum() );
	float pfo_E = pfo->getEnergy();

	if( pfo_E > _energy ){
	  float coneE = 0;
	  for( int j = 0; j < npfo; j++ ){
	    //calcuate charged cone energy
	    ReconstructedParticle* pfo_j = dynamic_cast< ReconstructedParticle* >( inputCol->getElementAt(j) );
	    if( pfo == pfo_j ) continue;
	    if( pfo_j->getCharge() == 0 ) continue;
	    TVector3 pfo_j_3mom = TVector3( pfo_j->getMomentum() );
	    float conecosth = pfo_3mom.Unit().Dot( pfo_j_3mom.Unit() );
	    if( conecosth > _coneCosTheta ) coneE += pfo_j->getEnergy();
	  }
	  float ratioE = coneE / pfo_E;

	  if( ratioE < _ratio ){
	    //regarded as ISR
	    isrCol->addElement( pfo );
	  }
	  else{
	    //not ISR by ratio
	    outputCol->addElement( pfo );
	  }
	}
	else{
	  //not ISR by costh and energy
	  outputCol->addElement( pfo );
	}
      }

    }
  }

  evt->addCollection( outputCol, _outputCollection.c_str() );
  evt->addCollection( isrCol, _isrCollection.c_str() );
}

void ISRFinder::check( LCEvent *evt ){
}

void ISRFinder::end(){
}
