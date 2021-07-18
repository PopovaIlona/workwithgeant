#include "EnergyTimeSD.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4GammaConversion.hh"
#include "G4Electron.hh"
#include "G4OpticalPhoton.hh"
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
EnergyTimeSD::EnergyTimeSD(G4String name) :
  G4VSensitiveDetector(name)
{
    collectionName.insert("energy_time");
	
}


G4bool EnergyTimeSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
    //total energy deposit, global time and position from the step
	G4Track* track = aStep->GetTrack();
	G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();
	G4double edep = aStep->GetTotalEnergyDeposit();
	G4double deposit = aStep->GetTotalEnergyDeposit();
	
	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4ThreeVector posit = preStepPoint->GetPosition();
	const G4ParticleDefinition* particle = 
									track->GetDynamicParticle()->GetDefinition();
	G4double Time = preStepPoint->GetGlobalTime();
	G4int trackID                = track->GetTrackID();			 //ID of the track
	G4int parentID               = track->GetParentID();
	G4ThreeVector mom = track->GetMomentumDirection();
	G4double KinE = track->GetKineticEnergy();
		    		

   	EnergyTimeHit* hit = new EnergyTimeHit();
   	
	const G4VProcess* processDefinedTheStep
		= aStep->GetPostStepPoint()->GetProcessDefinedStep();
  
		
	if (particle == G4Gamma::Gamma()) //particle != G4Gamma::Gamma() && trackID == 1
	 {
	  hit->SetKinEnergy(KinE);
	  hit->SetMoment(mom);
	  hit->SetTime(Time);
      hit->SetPosition(posit);
	  hit->SetTotalEnergyDeposit(deposit);}
	  fHitsCollection ->insert(hit);
      return true;
  }


void EnergyTimeSD::Initialize(G4HCofThisEvent* hcof)
{
 fHitsCollection = new EnergyTimeHitsCollection(SensitiveDetectorName,
	                                                                collectionName[0]);
 if (fHitsCollectionId < 0)
  {
   fHitsCollectionId = G4SDManager
	           ::GetSDMpointer()->GetCollectionID(GetName() + "/" + collectionName[0]);
  }
   hcof->AddHitsCollection(fHitsCollectionId, fHitsCollection);
}
