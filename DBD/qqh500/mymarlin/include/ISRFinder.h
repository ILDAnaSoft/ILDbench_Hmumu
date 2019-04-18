#ifndef ISRFinder_h
#define ISRFinder_h 1

#include <string>

#include <marlin/Processor.h>
#include <lcio.h>
#include <EVENT/ReconstructedParticle.h>

using namespace lcio;
using namespace marlin;

class ISRFinder : public Processor{
 public:
  virtual Processor* newProcessor() { return new ISRFinder; }
  ISRFinder();
  virtual void init();
  virtual void processRunHeader( LCRunHeader *run );
  virtual void processEvent( LCEvent* evt );
  virtual void check( LCEvent* evt );
  virtual void end();

 protected:
  std::string _inputCollection;
  std::string _outputCollection;
  std::string _isrCollection;

  float _energy;
  float _coneCosTheta;
  float _ratio;
};
#endif
