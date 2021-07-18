#ifndef ENERGYTIMEHIT_HH
#define ENERGYTIMEHIT_HH

#include <G4VHit.hh>
#include <G4THitsMap.hh>
#include <G4ThreeVector.hh>

/**
  * Custom hit class used in task 4d.
  *
  * It holds infomation about energy deposits and position/time when
  * traversing a layer.
  */
class EnergyTimeHit : public G4VHit
{
public:
    // Memory allocation and de-allocation
    inline void* operator new(size_t);
    inline void  operator delete(void*);
    
	void SetKinEnergy(G4double KinE) { fKinEnergy = KinE; }
    void SetDeltaEnergy(G4double deltaE) { fDeltaEnergy = deltaE; }
	void SetTotalEnergyDeposit (G4double deposit) { fdeposit = deposit; }
    void SetTime(G4double time) { fTime = time; }
    void SetPosition(G4ThreeVector pos) { fPosition = pos; }
	void SetMoment(G4ThreeVector mom) { fMoment = mom; }

	G4double GetKinEnergy() const { return fKinEnergy; }
    G4double GetDeltaEnergy() const { return fDeltaEnergy; }
	G4double GetTotalEnergyDeposit() const { return fdeposit; }
    G4double GetTime() const { return fTime; }
    G4ThreeVector GetPosition() const { return fPosition; }
	G4ThreeVector GetMoment() const { return fMoment; }

private:
	G4double fKinEnergy;
    G4double fDeltaEnergy;
    G4double fTime;
	G4double fdeposit;
    G4ThreeVector fPosition;
	G4ThreeVector fMoment;
};

using EnergyTimeHitsCollection = G4THitsCollection<EnergyTimeHit>;

extern G4ThreadLocal G4Allocator<EnergyTimeHit> *hitAllocator;

inline void* EnergyTimeHit::operator new(size_t)
{
  if (!hitAllocator)
  {
      hitAllocator = new G4Allocator<EnergyTimeHit>;
  }
  return hitAllocator->MallocSingle();
}

inline void EnergyTimeHit::operator delete(void *aHit)
{
    if (!hitAllocator)
    {
        hitAllocator = new G4Allocator<EnergyTimeHit>;
    }
    hitAllocator->FreeSingle((EnergyTimeHit*) aHit);
}

#endif
