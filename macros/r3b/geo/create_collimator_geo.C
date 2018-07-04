#include <iomanip>
#include <iostream>
#include "TGeoManager.h"
#include "TMath.h"

// Create Matrix Unity
TGeoRotation *fGlobalRot = new TGeoRotation();

// Create a null translation
TGeoTranslation *fGlobalTrans = new TGeoTranslation();
TGeoRotation *fRefRot = NULL;

TGeoManager*  gGeoMan = NULL;

Double_t fThetaX = 0.;
Double_t fThetaY = 0.;
Double_t fThetaZ = 0.;
Double_t fPhi   = 0.;
Double_t fTheta = 0.;
Double_t fPsi   = 0.;
Double_t fX = 0.;
Double_t fY = 0.;
Double_t fZ = 0.;
Bool_t fLocalTrans = kFALSE;
Bool_t fLabTrans = kTRUE;


TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans *fRef);


void create_collimator_geo(const char* geoTag)
{

  fGlobalTrans->SetTranslation(0.0,0.0,0.0);

  // -------   Load media from media file   -----------------------------------
  FairGeoLoader*    geoLoad = new FairGeoLoader("TGeo","FairGeoLoader");
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  TString geoPath = gSystem->Getenv("VMCWORKDIR");
  TString medFile = geoPath + "/geometry/media_r3b.geo";
  geoFace->setMediaFile(medFile);
  geoFace->readMedia();
  gGeoMan = gGeoManager;
  // --------------------------------------------------------------------------



  // -------   Geometry file name (output)   ----------------------------------
  TString geoFileName = geoPath + "/geometry/collimator_";
  geoFileName = geoFileName + geoTag + ".geo.root";
  // --------------------------------------------------------------------------



  // -----------------   Get and create the required media    -----------------
  FairGeoMedia*   geoMedia = geoFace->getMedia();
  FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

  FairGeoMedium* mAir      = geoMedia->getMedium("Air");
  if ( ! mAir ) Fatal("Main", "FairMedium Air not found");
  geoBuild->createMedium(mAir);
  TGeoMedium* pMed2 = gGeoMan->GetMedium("Air");
  if ( ! pMed2 ) Fatal("Main", "Medium Air not found");

  FairGeoMedium* mVac      = geoMedia->getMedium("vacuum");
  if ( ! mVac ) Fatal("Main", "FairMedium vacuum not found");
  geoBuild->createMedium(mVac);
  TGeoMedium* pMed1 = gGeoMan->GetMedium("vacuum");
  if ( ! pMed1 ) Fatal("Main", "Medium vacuum not found");
  
  FairGeoMedium* mFe      = geoMedia->getMedium("iron");
  if ( ! mFe ) Fatal("Main", "FairMedium iron not found");
  geoBuild->createMedium(mFe);
  TGeoMedium* pMedFe = gGeoMan->GetMedium("iron");
  if ( ! pMedFe ) Fatal("Main", "Medium iron not found");
  // --------------------------------------------------------------------------



  // --------------   Create geometry and top volume  -------------------------
  gGeoMan = (TGeoManager*)gROOT->FindObject("FAIRGeom");
  gGeoMan->SetName("collimatorgeom");
  TGeoVolume* top = new TGeoVolumeAssembly("TOP");
  gGeoMan->SetTopVolume(top);
  // --------------------------------------------------------------------------



  // out-of-file geometry definition
  Double_t dx,dy,dz,rmin,rmax;
  Double_t a;
  Double_t thx, phx, thy, phy, thz, phz;
  Double_t z, density, radl, absl, w;
  Int_t nel, numed;
  
  
  //collimator 1
  dx = 0.0;  
  dy = 0.0;
  dz = -220.0;
  
  TGeoRotation *gRot1 = new TGeoRotation();
  gRot1->RotateX(0.);
  gRot1->RotateY(0.); 
  gRot1->RotateZ(0.);
  
  TGeoCombiTrans* pMatrix1 = new TGeoCombiTrans("", dx,dy,dz,gRot1);
  
  //---------------------------------------------------------------------------
  // Shape: CollimatorEnveloppe type: TGeoTube
  rmin = 0.;
  rmax = 5.;//5  
  dz   = 128.;//4
  TGeoShape *pCollimatorEnveloppe = new TGeoTube("CollimatorEnveloppe",rmin,rmax,dz);
  TGeoVolume *pCollimatorEnveloppe_log = new TGeoVolume("CollimatorEnveloppe",pCollimatorEnveloppe, pMed1);
  pCollimatorEnveloppe_log->SetVisLeaves(kTRUE);
  TGeoCombiTrans* pGlobal = GetGlobalPosition(pMatrix1);
  top->AddNode(pCollimatorEnveloppe_log, 0, pGlobal);
  
  //Collimator1
  // Combi transformation:
  dx = 0.000000;
  dy = 0.000000;
  dz = 0.000000;
  // Rotation:
  thx = 90.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *pMatrix7 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*
  pMatrix6 = new TGeoCombiTrans("", dx,dy,dz,pMatrix7);
  
  // Shape: Collimator1 type: TGeoTube
  rmin = 0.20000;
  rmax = 5.;
  dz   = 4.;
  TGeoShape *pCollimator1 = new TGeoTube("Collimator1",rmin,rmax,dz);
  TGeoVolume *pCollimator1_log = new TGeoVolume("Collimator1",pCollimator1, pMedFe);   
  pCollimator1_log->SetVisLeaves(kTRUE);
  pCollimatorEnveloppe_log->AddNode(pCollimator1_log, 0, pMatrix6);
  
  //Collimator2
  // Combi transformation:
  dx = 0.000000;
  dy = 0.000000;
  dz = 120.;
  // Rotation:
  thx = 90.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *pMatrix9 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*
  pMatrix8 = new TGeoCombiTrans("", dx,dy,dz,pMatrix9);
  
  // Shape: Collimator2 type: TGeoTube
  rmin = 1.;
  rmax = 5.;
  dz   = 4.;
  TGeoShape *pCollimator2 = new TGeoTube("Collimator2",rmin,rmax,dz);
  TGeoVolume *pCollimator2_log = new TGeoVolume("Collimator2",pCollimator2, pMedFe);   
  pCollimator2_log->SetVisLeaves(kTRUE);
  pCollimatorEnveloppe_log->AddNode(pCollimator2_log, 0, pMatrix8);
 
  
  // ---------------   Finish   -----------------------------------------------
  gGeoMan->CloseGeometry();
  gGeoMan->CheckOverlaps(0.001);
  gGeoMan->PrintOverlaps();
  gGeoMan->Test();

  TFile* geoFile = new TFile(geoFileName, "RECREATE");
  top->Write();
  geoFile->Close();
  // --------------------------------------------------------------------------
}



TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans *fRef)
{
  if (fLocalTrans == kTRUE ) {
    
    if ( ( fThetaX == 0 )  && ( fThetaY==0 )  && ( fThetaZ == 0 )
	 &&
	 ( fX == 0 ) && ( fY == 0 ) && ( fZ == 0 )
	 )  return fRef;
    
    
    // X axis
    Double_t xAxis[3] = { 1. , 0. , 0. };
    Double_t yAxis[3] = { 0. , 1. , 0. };
    Double_t zAxis[3] = { 0. , 0. , 1. };
    // Reference Rotation
    fRefRot = fRef->GetRotation();
    
    if (fRefRot) {
      Double_t mX[3] = {0.,0.,0.};
      Double_t mY[3] = {0.,0.,0.};
      Double_t mZ[3] = {0.,0.,0.};
      
      fRefRot->LocalToMasterVect(xAxis,mX);
      fRefRot->LocalToMasterVect(yAxis,mY);
      fRefRot->LocalToMasterVect(zAxis,mZ);
      
      Double_t a[4]={ mX[0],mX[1],mX[2], fThetaX };
      Double_t b[4]={ mY[0],mY[1],mY[2], fThetaY };
      Double_t c[4]={ mZ[0],mZ[1],mZ[2], fThetaZ };
      
      ROOT::Math::AxisAngle aX(a,a+4);
      ROOT::Math::AxisAngle aY(b,b+4);
      ROOT::Math::AxisAngle aZ(c,c+4);
      
      ROOT::Math::Rotation3D fMatX( aX );
      ROOT::Math::Rotation3D fMatY( aY );
      ROOT::Math::Rotation3D fMatZ( aZ );
      
      ROOT::Math::Rotation3D  fRotXYZ = (fMatZ * (fMatY * fMatX));
      
      //cout << fRotXYZ << endl;
      
      Double_t fRotable[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      fRotXYZ.GetComponents(
			    fRotable[0],fRotable[3],fRotable[6],
			    fRotable[1],fRotable[4],fRotable[7],
			    fRotable[2],fRotable[5],fRotable[8]
			    );
      TGeoRotation *pRot = new TGeoRotation();
      pRot->SetMatrix(fRotable);
      TGeoCombiTrans *pTmp = new TGeoCombiTrans(*fGlobalTrans,*pRot);
      
      // ne peut pas etre applique ici
      // il faut differencier trans et rot dans la multi.
      TGeoRotation rot_id;
      rot_id.SetAngles(0.0,0.0,0.0);
      
      TGeoCombiTrans c1;
      c1.SetRotation(rot_id);
      const Double_t *t = pTmp->GetTranslation();
      c1.SetTranslation(t[0],t[1],t[2]);
      
      TGeoCombiTrans c2;
      c2.SetRotation(rot_id);
      const Double_t *tt = fRefRot->GetTranslation();
      c2.SetTranslation(tt[0],tt[1],tt[2]);
      
      TGeoCombiTrans cc = c1 * c2 ;
      
      TGeoCombiTrans c3;
      c3.SetRotation(pTmp->GetRotation());
      TGeoCombiTrans c4;
      c4.SetRotation(fRefRot);
      
      TGeoCombiTrans ccc = c3 * c4;
      
      TGeoCombiTrans pGlobal;
      pGlobal.SetRotation(ccc.GetRotation());
      const Double_t *allt = cc.GetTranslation();
      pGlobal.SetTranslation(allt[0],allt[1],allt[2]);
      
      return  ( new TGeoCombiTrans( pGlobal ) );
      
    }else{
      
      cout << "-E- R3BDetector::GetGlobalPosition() \
	      No. Ref. Transformation defined ! " << endl;
      cout << "-E- R3BDetector::GetGlobalPosition() \
	      cannot create Local Transformation " << endl;
      return NULL;
    } //! fRefRot
    
  } else {
    // Lab Transf.
    if ( ( fPhi == 0 )  && ( fTheta==0 )  && ( fPsi == 0 )
	 &&
	 ( fX == 0 ) && ( fY == 0 ) && ( fZ == 0 )
	 )  return fRef;
    
    
    return ( new TGeoCombiTrans(*fGlobalTrans,*fGlobalRot) );
    
  }
}

