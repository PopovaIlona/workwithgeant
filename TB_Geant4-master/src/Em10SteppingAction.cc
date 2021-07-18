//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm10/src/Em10SteppingAction.cc
/// \brief Implementation of the Em10SteppingAction class
//
//
// $Id: Em10SteppingAction.cc 73033 2013-08-15 09:24:45Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4SystemOfUnits.hh"
#include  "Analysis.hh"
#include "Em10DetectorConstruction.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Em10SteppingAction.hh"
#include "Em10EventAction.hh"
#include "Em10RunAction.hh"
#include "G4Event.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"
#include <iomanip>
#include <G4Track.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10SteppingAction::Em10SteppingAction(Em10EventAction* EA, Em10RunAction* RA)
  :G4UserSteppingAction(),eventaction (EA),runaction (RA),
   IDold(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10SteppingAction::~Em10SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  G4double Theta,Thetaback,Ttrans,Tback,Tsec,Egamma,yend,zend,rend;

  G4int evno                   = eventaction->GetEventno(); //number of events
  const G4Track* track         = aStep->GetTrack();
  const G4StepPoint* prePoint  = aStep->GetPreStepPoint();  //prestep point
  const G4StepPoint* postPoint = aStep->GetPostStepPoint();	//poststep point
  G4int trackID                = track->GetTrackID();		//ID of the track
  G4int parentID               = track->GetParentID();
		
  const G4DynamicParticle* dynParticle = track->GetDynamicParticle();
  const G4ParticleDefinition* particle = dynParticle->GetDefinition();
  G4VPhysicalVolume* preVol            = prePoint->GetPhysicalVolume();
  G4VPhysicalVolume* postVol           = aStep->GetPostStepPoint()->GetPhysicalVolume();

  IDnow = evno+10000*trackID+100000000*parentID;               //WHAT IS IT??????????

 
  
		//energy deposition from step (not good!!!)
 if(preVol->GetName()=="Absorber"	&& trackID > 1 
									&& postPoint->GetStepStatus() == fGeomBoundary) 
	{ 
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  //G4double Egamma1;
	  //Egamma1 = aStep->GetTotalEnergyDeposit();
      //analysisManager->FillH1(5, Egamma1 / keV);
	}

  //transition X-ray spectrum incident && angular distribution (IN WINDOWR) 
  if(preVol->GetName()=="WindowRR" 		&& track->GetMomentumDirection().z()>0. 
										&& particle == G4Gamma::Gamma() 
										&& prePoint->GetStepStatus() == fGeomBoundary 
										&& preVol != postVol && trackID > 1) 
	{	
      G4double energ = aStep->GetTrack()->GetKineticEnergy();
      G4double tet = std::acos(aStep->GetTrack()->GetMomentumDirection().z());
   	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      analysisManager->FillH1(3, energ / keV);
	  analysisManager->FillH1(4, tet);
	  analysisManager->FillH2(2 , tet, energ / keV);
	}

  if(preVol->GetName()=="Box1" 		&& track->GetMomentumDirection().z()>0. 
										&& particle == G4Gamma::Gamma() 
										&& prePoint->GetStepStatus() == fGeomBoundary 
										&& preVol != postVol && trackID > 1) 
	{	
      G4double energ1 = aStep->GetTrack()->GetKineticEnergy();
      G4double tet1 = std::acos(aStep->GetTrack()->GetMomentumDirection().z());
   	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      analysisManager->FillH1(9, energ1 / keV);
	  analysisManager->FillH1(10, tet1);
	  //analysisManager->FillH2(2 , tet1, energ1 / keV);
	}
	
	G4Track* theTrack = aStep->GetTrack();	
	G4ThreeVector v3_unitX(0.0, 0.0, 1.0);
	G4double tetta = acos(((aStep->GetTrack()->GetMomentumDirection())*v3_unitX)/
								(sqrt((aStep->GetTrack()->GetMomentumDirection())*
			(aStep->GetTrack()->GetMomentumDirection()))*sqrt(v3_unitX*v3_unitX)));
 
	//just to kill particle

/* if(preVol->GetName() == "WindowR"  && particle != G4Gamma::Gamma())
 {
       theTrack -> SetTrackStatus(fStopAndKill);
 }
 */


  if(preVol->GetName()=="Absorber") {
    if(particle == G4Electron::Electron() ||
       particle == G4Positron::Positron()) {
      eventaction->CountStepsCharged();

    } else if(particle == G4Gamma::Gamma()) {
      eventaction->CountStepsNeutral();
		///////
    }

    if(prePoint->GetStepStatus() == fGeomBoundary &&
       preVol != postVol) {

      if(trackID == 1) {
        if(track->GetMomentumDirection().z()>0.) {

          eventaction->SetTr();
          Theta = std::acos(track->GetMomentumDirection().z());
          runaction->FillTh(Theta);

          Ttrans = track->GetKineticEnergy();
          runaction->FillTt(Ttrans);
		  
          yend= aStep->GetTrack()->GetPosition().y();
          zend= aStep->GetTrack()->GetPosition().x();
          rend = std::sqrt(yend*yend+zend*zend);
          runaction->FillR(rend);
		  
        } else {
          eventaction->SetRef();
          Thetaback = std::acos(aStep->GetTrack()->GetMomentumDirection().z());
          Thetaback -= 0.5*pi;
          runaction->FillThBack(Thetaback);
          Tback  = aStep->GetTrack()->GetKineticEnergy();
          runaction->FillTb(Tback);
        }
      }
      if(track->GetMomentumDirection().z()>0. &&
         particle == G4Gamma::Gamma()) 
	  {
			// G4double Egamma1;
		//Egamma1 = aStep->GetTotalEnergyDeposit();
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		//analysisManager->FillH1(3, Egamma1 / eV);
		  G4double tetta;
		  tetta = std::acos(aStep->GetTrack()->GetMomentumDirection().z());
        Egamma = aStep->GetTrack()->GetKineticEnergy();
        runaction->FillGammaSpectrum(Egamma);
		analysisManager->FillH1(1, Egamma / keV);
		analysisManager->FillH1(2, tetta);
   
	  
	 }
  }
}
  }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
