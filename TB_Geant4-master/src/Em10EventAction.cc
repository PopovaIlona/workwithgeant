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
/// \file electromagnetic/TestEm10/src/Em10EventAction.cc
/// \brief Implementation of the Em10EventAction class
//
//
// $Id: Em10EventAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em10EventAction.hh"

#include "Em10RunAction.hh"
#include "G4SystemOfUnits.hh"
#include "Em10CalorHit.hh"
#include "Em10EventActionMessenger.hh"
#include "Analysis.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <G4SDManager.hh>
#include <G4THitsMap.hh>
#include <G4SystemOfUnits.hh>
#include <G4Event.hh>
#include "EnergyTimeHit.hh"
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em10EventAction::Em10EventAction(Em10RunAction* Em10RA)
:G4UserEventAction(),calorimeterCollID(-1),eventMessenger(0),
 runaction(Em10RA),verboselevel(0),drawFlag("all"),printModulo(10000)
{
  eventMessenger = new Em10EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em10EventAction::~Em10EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
     
  if(verboselevel>1)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;
    
  if (calorimeterCollID==-1)
    {
     G4SDManager * SDman = G4SDManager::GetSDMpointer();
     calorimeterCollID = SDman->GetCollectionID("CalCollection");
    } 

  nstep = 0. ;
  nstepCharged = 0. ;
  nstepNeutral = 0. ;
  Nch = 0. ;
  Nne = 0. ;
  NE=0.;
  NP=0.;
  Transmitted=0.;
  Reflected  =0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::EndOfEventAction(const G4Event* evt)
{
	//Part for NEW SD
    G4SDManager* sdm = G4SDManager::GetSDMpointer();
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();
	G4HCofThisEvent* hcofevt = evt->GetHCofThisEvent();
 	G4int fabsorberETId = sdm->GetCollectionID("absorberET/energy_time");
	   
    // Hit collections IDs to be looped over
    vector<G4int> hitCollectionIds = {fabsorberETId};
    for (G4int collectionId : hitCollectionIds)
    {
	 if (collectionId == -1)
   		{
		 G4cout << "evtAction: scorer ID: " << collectionId << G4endl;
		 continue;
    	}
	EnergyTimeHitsCollection* hitCollection = 
				dynamic_cast<EnergyTimeHitsCollection*>(hcofevt->GetHC(collectionId));

      if (hitCollection)
        {
			G4double h_E_tmp=0;
			G4double ED = 0;
           for (auto hit: *hitCollection->GetVector())
            	{
				G4ThreeVector v3_unitX(0.0, 0.0, 1.0);
				G4double tetta = acos(((hit->GetMoment())*v3_unitX)/(sqrt((hit->GetMoment())*(hit->GetMoment()))*sqrt(v3_unitX*v3_unitX)));
              	
				analysis->FillNtupleDColumn(0, hit->GetTotalEnergyDeposit() / keV);
				analysis->FillNtupleDColumn(1, hit->GetTime() / (0.001* ns));
				analysis->FillNtupleDColumn(2, hit->GetPosition().x() / cm);
				analysis->FillNtupleDColumn(3, hit->GetPosition().y() / cm);
				analysis->FillNtupleDColumn(4, hit->GetPosition().z() / cm - 232.5);
				analysis->FillNtupleDColumn(5, hit->GetKinEnergy() / keV);
				analysis->FillNtupleDColumn(6, tetta/ rad);
				analysis->AddNtupleRow();
						
				analysis->FillH2(1 , hit->GetPosition().x() / cm, hit->GetPosition().y() / cm);
				//analysis->FillH2(3 , tetta/ rad, hit->GetTotalEnergyDeposit() / keV);
				//analysis->FillH2(2 , hit->GetPosition().x() / cm, hit->GetPosition().z() / cm);
				analysis->FillH3(1 , hit->GetPosition().x() / cm, hit->GetPosition().y() / cm, hit->GetPosition().z() / mm); //hit->GetKinEnergy() / keV);//hit->GetTime() / (ns)
				//analysis->FillH1(8,  hit->GetKinEnergy()/ eV);
		        analysis->FillH1(6, tetta/ rad);
				//analysis->FillH1(7, hit->GetMoment().z());
						 
				// if(hit->GetTotalEnergyDeposit()>0 * keV) h_E_tmp += hit->GetTotalEnergyDeposit() / keV;
				int n = rand()%10+1;
				// G4cout << "event" << evt->GetEventID() << " - " << n << G4endl;
				if((hit->GetTotalEnergyDeposit()>34.5 * keV && n > 9) || hit->GetTotalEnergyDeposit()< 34.5* keV){
				h_E_tmp += hit->GetTotalEnergyDeposit() / keV;
				ED = hit->GetTotalEnergyDeposit() / keV;}
				analysis->FillH1(7, ED);
				analysis->FillH2(3 , tetta/ rad, ED);

            
	  }
				analysis->FillH1(8, h_E_tmp);

        }
    }




  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  Em10CalorHitsCollection* CHC = 0;
  if (HCE)
      CHC = (Em10CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
   {
    int n_hit = CHC->entries();
   // if(verboselevel==2)
   // G4cout << "     " << n_hit
   //      << " hits are stored in Em10CalorHitsCollection." << G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4double totEAbs=0, totLAbs=0;
	G4ThreeVector mom;
    for (int i=0;i<n_hit;i++)
      { totEAbs += (*CHC)[i]->GetEdepAbs(); 
        totLAbs += (*CHC)[i]->GetTrakAbs();
      }
  if(verboselevel==2)
    G4cout
       << "   Absorber: total energy: " << std::setw(7) << 
                             G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << std::setw(7) <<
                             G4BestUnit(totLAbs,"Length")
       << G4endl;           

   // count event, add deposits to the sum ...
    runaction->CountEvent() ;
    runaction->AddTrackLength(totLAbs) ;
    runaction->AddnStepsCharged(nstepCharged) ;
    runaction->AddnStepsNeutral(nstepNeutral) ;
    if(verboselevel==2)
      G4cout << " Ncharged=" << Nch << "  ,   Nneutral=" << Nne << G4endl;
    runaction->CountParticles(Nch,Nne);
    runaction->AddEP(NE,NP);
    runaction->AddTrRef(Transmitted,Reflected) ;
    runaction->AddEdeps(totEAbs) ;
    runaction->FillEn(totEAbs) ;
	 
	  //analysisManager->FillNtupleDColumn(0, totEAbs / keV);
	  
	  //analysisManager->AddNtupleRow();
      analysisManager->FillH1(5, totEAbs / keV);
    nstep=nstepCharged+nstepNeutral ;
    runaction->FillNbOfSteps(nstep);
  }

  if(verboselevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
  
  
  //save rndm status
  if (runaction->GetRndmFreq() == 2)
    { 
     CLHEP::HepRandom::saveEngineStatus("endOfEvent.rndm");   
     G4int evtNb = evt->GetEventID();
     if (evtNb%printModulo == 0)
       { 
        G4cout << "\n---> End of Event: " << evtNb << G4endl;
        CLHEP::HepRandom::showEngineStatus();
       }
    }     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Em10EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::CountStepsCharged()
{
  nstepCharged += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::CountStepsNeutral()
{
  nstepNeutral += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::AddCharged() 
{
  Nch += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::AddNeutral() 
{
  Nne += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::AddE() 
{
  NE += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::AddP() 
{
  NP += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::SetTr()
{
  Transmitted = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10EventAction::SetRef()
{
  Reflected   = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  


