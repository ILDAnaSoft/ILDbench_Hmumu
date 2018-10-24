#include "VertexInfo.h"
#include "HelixClass.h"

#include "vertex_lcfi/inc/lciointerface.h"
#include "vertex_lcfi/zvtop/include/VertexFitterKalman.h"
#include "vertex_lcfi/zvtop/include/interactionpoint.h"

#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <iostream>
#include <cassert>

using std::cout;
using std::endl;

void VertexInfo::addTrack( EVENT::Track* track ){
  // ensure that b field is nonzero
  assert ( fabs( _bField )>1e-5 );

  // convert to LCFIvertex Track object.
  // this need a recoParticle as input, with momentum, charge
  IMPL::ReconstructedParticleImpl* rp = new IMPL::ReconstructedParticleImpl();
  myrecoparts.push_back( rp );
  rp->addTrack( track );

  // calculate momentum of track at PCA
  HelixClass hel;
  hel.Initialize_Canonical( track->getPhi(),
                            track->getD0(),
                            track->getZ0(),
                            track->getOmega(),
                            track->getTanLambda(),
                            _bField );
  rp->setMomentum( hel.getMomentum() );
  rp->setCharge( hel.getCharge() );

  vertex_lcfi::Track* vtx_lcfi_trk = vertex_lcfi::trackFromLCIORP( NULL, rp );
  lcfi_tracks.push_back( vtx_lcfi_trk );

  _vtxValid=false;
  return;
}

void VertexInfo::calculateVertexPosition(){
  const float maxchisqcont=20;
  const float maxchisqdiff=10;

  if( _vtxValid ) return;
  std::vector<vertex_lcfi::TrackState *> trackStates;
  std::map<vertex_lcfi::TrackState *, const vertex_lcfi::Track *> trackMap;
  for( size_t i = 0; i < lcfi_tracks.size(); i++ ){
    vertex_lcfi::Track* ptr = lcfi_tracks[i];
    vertex_lcfi::TrackState *pstate = ptr->makeState();

    lcfi_trackstates.push_back( pstate );

    trackStates.push_back( pstate );
    trackMap[ pstate ] = lcfi_tracks[i];
  }
  vertex_lcfi::ZVTOP::VertexFitterKalman kalman;
  vertex_lcfi::util::Vector3 vtxpos;
  vertex_lcfi::util::Matrix3x3 vtxerr;
  double chi2fit;
  double chi2ip;
  std::map<vertex_lcfi::TrackState*, double> chi2map;
  vertex_lcfi::ZVTOP::InteractionPoint *ipConst = 0;

  if( _seed.Mag() > 0 ){
    vertex_lcfi::util::Vector3 vtxposSeed( _seed[0], _seed[1], _seed[2] );
    kalman.setSeed( vtxposSeed );
  }

  if( _useIPcons ){//use the IP constraint
    vertex_lcfi::util::Vector3 ippos(0,0,0);
    vertex_lcfi::util::SymMatrix3x3 iperr;
    iperr(0,0) = pow(_ipSize[0],2);
    iperr(1,1) = pow(_ipSize[1],2);
    iperr(2,2) = pow(_ipSize[2],2);
    iperr(1,0)=iperr(0,1)=0; //assume uncorrelated
    iperr(2,0)=iperr(0,2)=0;
    iperr(2,1)=iperr(1,2)=0;
    ipConst = new vertex_lcfi::ZVTOP::InteractionPoint( ippos, iperr );
  }

  if( !_trimTracks ){//just use all tracks
    kalman.fitVertex(trackStates,ipConst, vtxpos, vtxerr, chi2fit, chi2map, chi2ip);
  }
  else{//trim off bad tracks
    while( trackStates.size() > 1 ){
      kalman.fitVertex( trackStates,ipConst, vtxpos, vtxerr, chi2fit, chi2map, chi2ip );
      float largestchidqcont(0);
      float secondlargestchidqcont(0);
      vertex_lcfi::TrackState* worstTrack(0);
      for( std::map<vertex_lcfi::TrackState*, double>::iterator irr=chi2map.begin(); irr!=chi2map.end(); irr++ ){
	if( irr->second>largestchidqcont ){
	  secondlargestchidqcont=largestchidqcont;
	  largestchidqcont=irr->second;
	  worstTrack=irr->first;
	}
      }
      if( largestchidqcont>maxchisqcont || largestchidqcont-secondlargestchidqcont > maxchisqdiff ){
	if( find( trackStates.begin(), trackStates.end(), worstTrack )!=trackStates.end() ){
	  trackStates.erase( find( trackStates.begin(), trackStates.end(), worstTrack ) );//remove worst track
	}
	else{
	  cout << "WEIRD, shouldn;t happen!" << endl;
	}
      }
      else{// the track looks OK, so stop trimming
	break;
      }      
    }//while loop
  }

  double ee[9];
  int ii(0);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      ee[ii++]=vtxerr(i,j);
    }
  }

  TMatrixDSym m(3, ee);
  TMatrixDSymEigen me(m);

  TVectorD eigenval = me.GetEigenValues();
  TMatrixD eigenvec = me.GetEigenVectors();

  for(int i=0; i<3; i++){
    _eigenValues[i] = eigenval[i];
    _eigenVectors[i].SetXYZ(eigenvec(0,i), eigenvec(1,i), eigenvec(2,i) );
  }

  _vtxPos.SetXYZ( vtxpos.x(), vtxpos.y(), vtxpos.z() );
  _chisq = chi2fit;

  _vtxValid=true;

  if(ipConst) delete ipConst;
  return;
}

float VertexInfo::getVertexZ0(float x, float y){
  // calcaulte impact paramter of vertex wrt beamline @ x,y
  calculateVertexPosition();

  TVector2 pp = _eigenVectors[0].XYvector();
  TVector2 qq = TVector2(x,y) - _vtxPos.XYvector();
  float alpha = (pp.X()*qq.X() + pp.Y()*qq.Y())/pp.Mod2();

  float z0 = _vtxPos.Z() + alpha*_eigenVectors[0].Z();

  return z0;
}

void VertexInfo::cleanup(){
  for ( size_t i=0; i<myrecoparts.size(); i++) {
    delete myrecoparts[i];
  }
  myrecoparts.clear();

  for ( size_t i=0; i<lcfi_tracks.size(); i++) {
    delete lcfi_tracks[i];
  }
  lcfi_tracks.clear();


  for ( size_t i=0; i<lcfi_trackstates.size(); i++) {
    delete lcfi_trackstates[i];
  }
  lcfi_trackstates.clear();
}
