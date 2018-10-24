#ifndef _MYVTX_H_
#define _MYVTX_H_

#include "EVENT/Track.h"
#include "vertex_lcfi/inc/event.h"
#include "vertex_lcfi/inc/track.h"
#include "vertex_lcfi/util/inc/memorymanager.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "TVector3.h"
#include <marlinutil/GeometryUtil.h>

/*
some simple documentation by Daniel Jeans
1. make a VertexInfo object
VertexInfo vv();
3. add the LCIO::tracks you want to fit
vv.addTrack( trk1 );
vv.addTrack( trk2 );
4. if you want a beamspot constraint
float xyz[3]={10.e-3, 10.e-6, 0.3}; // in mm
vv.setBeamSpotSize(xyz);
vv.useIPcon(true);
5. do the fit, get the vtx position
TVector3 vtxPos = vv.getVertexPosition();
float vtxChisq = vv.getVertexChisq();
6. for more details, get the 3 eigenvectors and corresponding eigenvalues of the vertex ellipse
for (int i=0; i<3; i++) {
  TVector3 evec = vv.getEigenVector(int i);
  float eval = vv.getEigenValue(int i);
}
*/

class VertexInfo{
 public:
  VertexInfo() {
    _vtxValid=false; 
    _vtxPos.SetXYZ(999,999,999); 
    _chisq=9999;
    _useIPcons =false;
    _trimTracks=false;
    // bfield
    _bField=MarlinUtil::getBzAtOrigin();

    // IP region size in mm, used for IP constraint
    _ipSize[0] = 10.e-3;
    _ipSize[1] = 10.e-6;
    _ipSize[2] = 0.3;
    for (int i=0; i<3; i++) {
      _eigenValues[i]=-999;
      _eigenVectors[i].SetXYZ(999,999,999);
      _seed[i]=0;
    }
  }
  ~VertexInfo() {cleanup();}

  void setBeamSpotSize( float* xyz ) { 
    for ( int i=0; i<3; i++) _ipSize[i]=xyz[i];
  }

  float getBField() {return _bField;}

  void addTrack( EVENT::Track* trk );

  TVector3 getVertexPosition() {calculateVertexPosition(); return _vtxPos;}
  void calculateVertexPosition( );
  int getNtrack() { return lcfi_tracks.size() ; }
  float getVertexChisq() {calculateVertexPosition(); return _chisq;}

  TVector3 getEigenVector(int i) { assert (i>=0 && i<3); calculateVertexPosition(); return _eigenVectors[i];}
  float getEigenValue(int i) {assert (i>=0 && i<3); calculateVertexPosition(); return _eigenValues[i]; }

  float getVertexZ0(float x=0, float y=0);

  bool isValid() {return _vtxValid;}

  void trimTracks( bool b=true ) { _trimTracks=b; _vtxValid=false; }
  void useIPcon( bool b=true ) { _useIPcons=b; _vtxValid=false; }

  void setSeedPos( TVector3 seed ) {_seed = seed;}

 private:

  bool _useIPcons;
  bool _trimTracks;
  float _bField;
  TVector3 _seed;

  std::vector < vertex_lcfi::Track* > lcfi_tracks;
  std::vector < vertex_lcfi::TrackState* > lcfi_trackstates;

  TVector3 _vtxPos;
  bool _vtxValid;
  float _chisq;

  TVector3 _eigenVectors[3];
  float _eigenValues[3];

  float _ipSize[3];

  std::vector < IMPL::ReconstructedParticleImpl* > myrecoparts;

  void cleanup();

};
#endif
