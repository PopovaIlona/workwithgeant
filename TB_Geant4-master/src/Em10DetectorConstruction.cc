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
/// \file electromagnetic/TestEm10/src/Em10DetectorConstruction.cc
/// \brief Implementation of the Em10DetectorConstruction class
//
//
// $Id: Em10DetectorConstruction.cc 73033 2013-08-15 09:24:45Z gcosmo $
//
//

#include "Em10DetectorConstruction.hh"
#include "Em10DetectorMessenger.hh"
#include "Em10CalorimeterSD.hh"
#include "Em10Materials.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"



#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4ProductionCuts.hh"

//#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <G4Orb.hh>
#include "G4UserLimits.hh"
#include <G4SDManager.hh>
#include <G4PSEnergyDeposit.hh>
#include "EnergyTimeSD.hh"
#include <G4GlobalMagFieldMessenger.hh>
#include <sstream>

внесение каких-то изменений для затеста GitHub

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10DetectorConstruction::Em10DetectorConstruction()
  :G4VUserDetectorConstruction(),
   fWorldChanged(false), fAbsorberMaterial(0), fGapMat(0), fSetUp("TRT1"),
   fWorldMaterial(0), fSolidWorld(0), fLogicWorld(0), fPhysicsWorld(0),
   fSolidRadiator(0),  fLogicRadiator(0),   fPhysicsRadiator(0),
   fRadiatorMat(0), fPipe(false), fPipeField(false),
   fSolidAbsorber(0),  fLogicAbsorber(0),   fPhysicsAbsorber(0),
   fMagField(0),       fCalorimeterSD(0),   fRegGasDet(0),
   fRadRegion(0), fMat(0), fsolidWinR(0), flogicWinR(0), fphysWinR(0),
   fsolidWin(0), flogicWin(0), fphysWin(0)
{
  fDetectorMessenger = new Em10DetectorMessenger(this);
  fMat               = new Em10Materials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10DetectorConstruction::~Em10DetectorConstruction()
{
  delete fDetectorMessenger;
  delete fMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Em10DetectorConstruction::Construct()
{
  return ConstructDetectorXTR();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Em10DetectorConstruction::ConstructDetectorXTR()
{
 // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  if( fSetUp == "TRT1" )//////ALICE --> TRT
  {
    return SimpleSetUpALICE();
  }
  else
  {
    G4cout <<
    "Experimental setup is unsupported. Check /XTRdetector/setup " <<G4endl;
    G4cout<<"Run default: TRT1 "<<G4endl;
    return SimpleSetUpALICE();

    //  return 0;
  }
}


		 	////////////////////	/////////////	//////////////////		 
					/////			/////	/////			/////						
					/////			/////////////			/////
					/////			////////				/////	   
  					/////			/////	////			/////		    
					/////			/////		////		///// 		


// Runs by : TestEm10 salice.mac

G4VPhysicalVolume* Em10DetectorConstruction::SimpleSetUpALICE()
{

	//Radiator, World, and detector parameters
  fWorldSizeZ        = 300.*cm;
  fWorldSizeR        = 20.*cm;  	
  fRadThickness      = 0.050*mm;//30.*um;//0.0215*mm;//	  //foil thickness
  fGasGap            = 3.00*mm;//2.*mm;//.5*mm;//		  //gap between foils
  fFoilNumber        = 30;//1;//						  //number of foils
  fAbsorberThickness = 30*mm;							  //detector thickness
  fAbsorberRadius    = 5.*cm;
  fWindowThick       = 35.*micrometer;					  //detector window thickness
  //fElectrodeThick    = 10.0*micrometer;
  //fGapThick          = 10.0*cm;  //fDetThickness      = 40.0*mm;
  //fDetLength         = 200.0*cm;
  //fDetGap            = 0.180*mm;
  fStartR            = 0.52*cm;
  fStartZ            = 200.0*mm -0.5*fWorldSizeZ;
  fModuleNumber      = 1;
  foilGasRatio       = fRadThickness/(fRadThickness+fGasGap); 
  G4int numberOfLayers = 1;
  G4double gapdetector = 210.*cm;
  G4double fGap = 2*mm;
  G4double numberOfPixelX = 10;
	G4double numberOfPixelY = 10;
	G4double pixelZ = 0.1*mm;
	G4double pixelY = 5*mm;
	G4double pixelX = 5*mm;


	 // Preparation of mixed radiator material
  G4double a, z, density, temperature, pressure;
  G4int nelements, natoms;
  G4NistManager* nist = G4NistManager::Instance();
 // G4Material* C2H4 = nist->FindOrBuildMaterial("G4_C2H4");
  G4Material*Galactic = new G4Material("Galactic", z = 1., a = 1.01*g / mole, 
  density = 1.e-25*g / cm3, kStateGas,temperature = 0.1*kelvin, 
                                                    pressure = 1.e-19*pascal);  
  G4Material* Xe    = new G4Material("Xenon",z=54., 
	                               a = 131.29*g/mole, density = 5.858*mg/cm3);
 // G4Material* CH2   = fMat->GetMaterial("CH2");
  G4Material* CO2   = fMat->GetMaterial("CO2");
  G4Material* PolyP = nist->FindOrBuildMaterial("G4_POLYPROPYLENE");

  G4Element* H  = new G4Element("Hydrogen","H" , z= 1., a =1.01*g/mole);  G4Element* C  = new G4Element("Carbon"  ,"C" , z= 6., a = 12.01*g/mole);
  G4Material* CH2low = new G4Material("CH2low"  , 0.46*g/cm3 ,2);  CH2low->AddElement(C, natoms=1);  CH2low->AddElement(H, natoms=2);
  G4Material* CH2 = new G4Material("CH2"  , 0.92*g/cm3 ,2);  CH2->AddElement(C, natoms=1);  CH2->AddElement(H, natoms=2);


  
  G4Material* Mylar = fMat->GetMaterial("Mylar");
  G4Material* Air   = fMat->GetMaterial("Air");
  //G4Material* Al    = fMat->GetMaterial("Al");
  G4Material* Al = new G4Material("Aluminum", z = 13.,  
	                               a = 26.98*g / mole, density = 2.7*g / cm3);
   G4Material* He   = fMat->GetMaterial("He");
  
  G4double foilDensity =   2.79*g / cm3; // Mylar //0.92*g/cm3; // 0.46*g/cm3;//C2H4 0.91*g/cm3;  // CH2 0.534*g/cm3; //Li
  G4double gasDensity  =  0.0012041*g/cm3; //0.00129*g/cm3;
	// Air // 1.977*mg/cm3; // CO2 0.178*mg/cm3; //He
  G4double totDensity   = foilDensity*foilGasRatio + gasDensity*(1.0-foilGasRatio);
  G4double fractionFoil = foilDensity*foilGasRatio/totDensity; 
  G4double fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;

  G4Material* radiatorMat = new G4Material("radiatorMat"  , totDensity,2);
  radiatorMat->AddMaterial( Mylar, fractionFoil);
  radiatorMat->AddMaterial( Air, fractionGas);
  G4Material* abs = new G4Material("abs", (0.2*0.0012041*g/cm3 + 0.8*5.858*mg/cm3) ,2);
  abs->AddMaterial( Xe, 0.8);
  abs->AddMaterial( CO2, 0.2);

  G4Material* pipe = new G4Material("pipe", (0.2*0.178*mg/cm3 + 0.8* 0.178*mg/cm3),2);
  pipe->AddMaterial( He, 0.8);
  pipe->AddMaterial( Air, 0.2);

	// default materials of the detector and TR radiator

  fRadiatorMat      = radiatorMat;//  Al; //
 // fRadiator1Mat     = radiator1Mat;
  fFoilMat          = Mylar;//CH2;//PolyP; //Al;// CH2;//Kapton;//Li;//CH2;//C2H4;
  fGasMat           = Air; //Galactic; //Air; // CO2; // He; //  
  fWindowMat        = Mylar ; //Galactic;//Air; // CH2;
  fAbsorberMaterial = abs; //Xe;//Galactic;//Air; //abs;//Galactic;// 
  fWorldMaterial    = Galactic;// Air; //He;// ;// Galactic;// CO2;
  
	//volumes definition
	//"World"

  G4Box* expHall_box1 = new G4Box("World", fWorldSizeR,fWorldSizeR,fWorldSizeZ/2.);

  G4LogicalVolume* fLogicWorld = new G4LogicalVolume(expHall_box1, 
												fWorldMaterial,  "World", 0, 0, 0);
  G4VPhysicalVolume* fPhysicsWorld
		= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0); 

 	// //color for vis
	// G4Colour color1(0, 0, 0,1);
	// G4VisAttributes* PlateVisAtt = new  G4VisAttributes(color1);
	// fLogicWorld->SetVisAttributes(PlateVisAtt);
	// //fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible);
	
	// TR radiator
  fRadThick = fFoilNumber*(fRadThickness + fGasGap) - fGasGap;// + fDetGap;
  fRadZ = fStartZ + 0.5*fRadThick;
  G4double thickness = fGap + fRadThick;

  fSolidRadiator   = new G4Box("Radiator", 
	                    1.1*fAbsorberRadius, 1.1*fAbsorberRadius,  0.5*fRadThick );
  fLogicRadiator   = new G4LogicalVolume(fSolidRadiator, fRadiatorMat, "Radiator"); 				
	
  vector<G4ThreeVector> radiatorPositions;
																	
 
    for (int i = 0; i < numberOfLayers; i++)
    {
		 radiatorPositions.push_back({0, 0, fRadZ + i*1 * thickness});
	}

    for (int i = 0; i < numberOfLayers; i++)
    {
        ostringstream aName; aName << "Radiator" << i;
																		
        new G4PVPlacement(nullptr, radiatorPositions[i], fLogicRadiator, 
			  aName.str(), fLogicWorld, 0, i);																									//bName.str(), fLogicWorld, 0, i);
        
    }
 
	//XTR region
  if( fRadRegion != 0 ) delete fRadRegion;
  if( fRadRegion == 0 ) fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);
  G4ThreeVector box1Position (0, 0, fRadZ + (numberOfLayers*0.5)* 
														thickness -fGap/2 + 0.0004*m);

  G4ThreeVector box4Position (0, 0, fRadZ + (numberOfLayers*0.5)* 
														  thickness -fGap/2 + 1.05*m);

  fsolidBox1 = new G4Box("Box1",fAbsorberRadius,fAbsorberRadius, 0.00005*m);
  flogicBox1 = new G4LogicalVolume(fsolidBox1,fWorldMaterial, "Box1", 0, 0, 0);
  fphysBox1  = new G4PVPlacement(0, box1Position,
										flogicBox1,"Box1",fLogicWorld,false,0);

 G4Box * fsolidBox4 = new G4Box("Box4",fAbsorberRadius,fAbsorberRadius, 1.0*m);
  G4LogicalVolume* flogicBox4 = new G4LogicalVolume(fsolidBox4,pipe, "Box4", 0, 0, 0);
  G4VPhysicalVolume * fphysBox4  = new G4PVPlacement(0, box4Position,			
											   flogicBox4,"Box4",fLogicWorld,false,0);

  G4Box *fsolidBox2 = new G4Box("Box2",fAbsorberRadius,fAbsorberRadius, 0.0075*mm);
  G4LogicalVolume* flogicBox2 = new G4LogicalVolume(fsolidBox2,PolyP, "Box2", 0, 0, 0);
  G4VPhysicalVolume* fphysBox2  = new G4PVPlacement(0, 
			G4ThreeVector(0,0, 1*m - 0.0075*mm), flogicBox2,"Box2",flogicBox4,false,0);

  G4Box *fsolidBox3 = new G4Box("Box3",fAbsorberRadius,fAbsorberRadius, 0.0075*mm);
  G4LogicalVolume* flogicBox3 = new G4LogicalVolume(fsolidBox3,PolyP, "Box3", 0, 0, 0);
  G4VPhysicalVolume* fphysBox3  = new G4PVPlacement(0,
		  G4ThreeVector(0,0, - 1*m + 0.0075*mm), flogicBox3,"Box3",flogicBox4,false,0);

  G4double box = fWindowThick/2 + fAbsorberThickness/2 + pixelZ/2;
  G4ThreeVector boxPosition (0, 0, fRadZ + (numberOfLayers*0.5)* 
	         thickness + gapdetector+fWindowThick + fAbsorberThickness/2.+pixelZ/2);
  G4ThreeVector WindowPosition (0, 0, -box + fWindowThick/2);
  G4ThreeVector AbsorberPosition (0, 0, -box + fWindowThick + fAbsorberThickness/2);
  
  G4RotationMatrix rotm  = G4RotationMatrix();
		rotm.rotateZ(0*deg); 	
		rotm.rotateX(0*deg); //15
		

  fsolidBox = new G4Box("Box",fAbsorberRadius,fAbsorberRadius, box);
  G4LogicalVolume* flogicBox = new G4LogicalVolume(fsolidBox,
													fWorldMaterial, "Box", 0, 0, 0);
  G4VPhysicalVolume* fphysBox  = new G4PVPlacement(G4Transform3D(rotm, boxPosition), 
											   flogicBox,"Box",fLogicWorld,false,0);
 
	//WindowR
  
  fsolidWinR = new G4Box("WindowR",fAbsorberRadius,fAbsorberRadius, fWindowThick);
  flogicWinR = new G4LogicalVolume(fsolidWinR,fWorldMaterial, "WindowR"); 
  fphysWinR  = new G4PVPlacement(0, WindowPosition, "WindowR",
													 flogicWinR,fphysBox,false,0);

  G4Box *fsolidWinRR = new G4Box("WindowRR",
								fAbsorberRadius,fAbsorberRadius, fWindowThick/4.);
  G4LogicalVolume* flogicWinRR = new G4LogicalVolume(fsolidWinRR,
											 fWorldMaterial, "WindowRR", 0, 0, 0);
  G4VPhysicalVolume* fphysWinRR  =  new G4PVPlacement(0, 
									 G4ThreeVector(0.,0.,-fWindowThick*0.5), 
									"WindowRR", flogicWinRR, fphysWinR, false, 0); 

	// window
  fsolidWin = new G4Box("Window",fAbsorberRadius, fAbsorberRadius, fWindowThick/2.); 
  flogicWin = new G4LogicalVolume(fsolidWin, fWindowMat, "Window");
  fphysWin  = new G4PVPlacement(0, G4ThreeVector(0.,0.,fWindowThick*0.5), "Window", 
													flogicWin, fphysWinR, false, 0); 

	// Absorber
  

  fSolidAbsorber   = new G4Box("Absorber", fAbsorberRadius, 
											fAbsorberRadius, fAbsorberThickness/2.);
  fLogicAbsorber   = new G4LogicalVolume(fSolidAbsorber, 
													 fAbsorberMaterial, "Absorber");
  fPhysicsAbsorber = new G4PVPlacement(0, AbsorberPosition,  "Absorber", 
											   fLogicAbsorber,fphysBox,  false,  0);

    vector<G4ThreeVector> PixelPositions;
	
	G4double startX = -(numberOfPixelX * pixelX)/2 + pixelX/2;
	G4double startY = -(numberOfPixelY * pixelY)/2 + pixelY/2;
	

	/*fsolidPixel = new G4Box("Pixel", pixelX/2, pixelY/2, pixelZ/2);

	flogicPixel = new G4LogicalVolume(fsolidPixel, fWindowMat, "Pixel");
	
	for (int p = 0; p < numberOfPixelX; p++){ 
		for (int u = 0; u < numberOfPixelY; u++)
		{ 
		//fphysPixel = new G4PVPlacement(0, 
		//G4ThreeVector(startX + p*pixelX, startY + u*pixelY, -box + fWindowThick + fAbsorberThickness + pixelZ/2),
		//"Pixel (i, u)", flogicPixel,fphysBox,  false,  0);
		}}*/
		
		
	//XTR detector region
  if( fRegGasDet != 0 ) delete fRegGasDet;
  if( fRegGasDet == 0 ) fRegGasDet = new G4Region("XTRdEdxDetector");
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);

	// Sensitive Detectors: Absorber 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );        //old sensitive detector
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  PrintGeometryParameters();

  //New SD
    EnergyTimeSD* absorberET = new EnergyTimeSD("absorberET");
	SDman->AddNewDetector(absorberET);
    fLogicAbsorber->SetSensitiveDetector(absorberET);


  return fPhysicsWorld;
}
void Em10DetectorConstruction::ConstructSDandField()
{
  
/*
  G4ThreeVector fieldValue = G4ThreeVector(0.0, 0.* tesla, 0.0);
  G4GlobalMagFieldMessenger* fMagFieldMessenger = 
					  new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);  // to set messenger verbosity
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::PrintGeometryParameters()
{
  G4cout << "\n The  WORLD   is made of "
       << fWorldSizeZ/mm << "mm of " << fWorldMaterial->GetName();
  G4cout << ", the transverse size (R) of the world is " << 
                                         fWorldSizeR/mm << " mm. " << G4endl;
  G4cout << " The ABSORBER is made of "
       << fAbsorberThickness/mm << "mm of " << fAbsorberMaterial->GetName();
  G4cout << ", the transverse size (R) is " << fAbsorberRadius/mm << 
            " mm. " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " 
         << fAbsorberZ/mm << "  mm." << G4endl;

  G4cout<<"fRadZ = "<<fRadZ/mm<<" mm"<<G4endl;
 
  G4cout<<"fStartZ = "<<fStartZ/mm<<" mm"<<G4endl;

  G4cout<<"fRadThick = "<<fRadThick/mm<<" mm"<<G4endl;
  G4cout<<"fFoilNumber = "<<fFoilNumber<<G4endl;
  G4cout<<"fRadiatorMat = "<<fRadiatorMat->GetName()<<G4endl;
  G4cout<<"WorldMaterial = "<<fWorldMaterial->GetName()<<G4endl;
  //  G4cout<<"fAbsorberZ = "<<fAbsorberZ/mm<<" mm"<<G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name
  G4Material* pttoMaterial;

  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
  {
    pttoMaterial = (*theMaterialTable)[J];
 
    if(pttoMaterial->GetName() == materialChoice)
    {
      fAbsorberMaterial = pttoMaterial;
      fLogicAbsorber->SetMaterial(pttoMaterial);
        // PrintCalorParameters();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetRadiatorMaterial(G4String materialChoice)
{
  // get the pointer to the material table

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name

  G4Material* pttoMaterial;
  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
  {
    pttoMaterial = (*theMaterialTable)[J];

    if(pttoMaterial->GetName() == materialChoice)
    {
      fRadiatorMat = pttoMaterial;
//      fLogicRadSlice->SetMaterial(pttoMaterial);
      // PrintCalorParameters();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name
  G4Material* pttoMaterial;

  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
  {
    pttoMaterial = (*theMaterialTable)[J];
 
    if(pttoMaterial->GetName() == materialChoice)
    {
      fWorldMaterial = pttoMaterial;
      fLogicWorld->SetMaterial(pttoMaterial);
       //  PrintCalorParameters();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  fAbsorberThickness = val;
  //  ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetRadiatorThickness(G4double val)
{
  // change XTR radiator thickness and recompute the calorimeter parameters
  fRadThickness = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetGasGapThickness(G4double val)
{
  // change XTR gas gap thickness and recompute the calorimeter parameters
  fGasGap = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberRadius(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  fAbsorberRadius = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldChanged=true;
  fWorldSizeZ = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetWorldSizeR(G4double val)
{
  fWorldChanged=true;
  fWorldSizeR = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberZpos(G4double val)
{
  fAbsorberZ  = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetMagField(G4double)
{
  //apply a global uniform magnetic field along X axis

  /* *********************************************************

  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;             //delete the existing magn field

  if(fieldValue!=0.)                        // create a new one if non null
  {
    magField = new G4UniformMagField(G4ThreeVector(fieldValue,0.,0.));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  }
  else
  {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }

  *************************************************************** */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetectorXTR());
}
