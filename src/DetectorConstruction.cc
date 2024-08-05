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
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"


#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

#include "Materials.hh"
//#include <Riostream.h>


DetectorConstruction::DetectorConstruction()
{
//  sizeXY = 8.7*cm;
  sizeXY = 550.0*cm;
  water_sizeZ = 30.*cm;
//  water_sizeXY = sizeXY;
//  water_sizeXY = 200*cm;
  water_sizeXY = 3.*cm;
  water_R = 15.*cm;
  sample_sizeXY = 2.0*cm;
//  sample_sizeXY = water_sizeXY;
//  sample_sizeXY = 6.*cm;
  sample_sizeZ = 5.*mm;
  alum_sizeZ = 0.4*cm;
//  alum_sizeZ = 1.835*cm; // for cylinder geometry
//  alum_sizeZ = 0.55*cm;  // from article
//!!!!!!!!!!!! change aluminium on glass illuminator
//  alum_sizeZ = 2.08*cm;
//  composite_sizeZ = 18.5*cm;
  composite_sizeZ = 2.0*cm;
  composite_thickness = 1.*cm;
  composite_Rmin = 7.7*cm;
  composite_Rmax = composite_Rmin+composite_thickness;
  cout<<"!!!!!!!!!!!!!!!!!rmin rmax "<<G4BestUnit(composite_Rmin,"Length")<<" "<<G4BestUnit(composite_Rmax,"Length")<<endl;
//  polyprop_sizeZ = 1.*um;
//  polyprop_sizeZ = 2.*cm;
//  polyprop_sizeZ = 0.5*cm; //new for 1.5g/cm2 0.5*cm
//  polyprop_sizeZ = 15.5*cm; //new for 15g/cm2 
//  polyprop_sizeZ = 0.5*cm; //new for 15g/cm2 
//  polyprop_sizeZ = 10.*cm; //new for 10g/cm2
  polyprop_sizeZ = 0.5*cm; // for test
//  polyprop_sizeZ1 = 10.*cm; //new for 10g/cm2 10.*cm
  polyprop_sizeZ1 = 10.*cm; //new for 1.5g/cm2 
//  polyprop_sizeZ1 = 15.5*cm; //new for 15g/cm2 
//  polyprop_sizeZ1 = 1.13*cm; //new for 10g/cm2 and Fe
//  polyprop_sizeZ = 11.7*cm;
//  lif_sizeX = 4.*cm;
//  lif_sizeY = 16.*cm;
//  lif_sizeZ = 0.5*cm;
  lif_sizeX = 5.*cm;
  lif_sizeY = 5.*cm;
  lif_sizeZ = 0.1*cm;
//  lif_R = 4.5*mm;
  lif_R = 0.8*cm;
  lif_H = 16.0*mm;
//  gap_sizeZ = 5.*mm;
  gap_sizeZ = 5.*mm;
//  gap_sizeZ = 15.*cm;
  gap_sizeZ1 = 3.*cm;
  gap_sizeZ2 = 2.*cm;
//  aluminium_R = 1.0*m;
//  aluminium_R = 1.45*m;  // for cylinder geometry
  aluminium_R = 1.5*m;  // for cylinder geometry
//  aluminium_lengthZ = 350.*cm;
  aluminium_lengthZ = 300.*cm;
  aluminium_positionZ = composite_sizeZ+gap_sizeZ+polyprop_sizeZ+0.5*alum_sizeZ;
  fCompDet = "yes";
  ComputeCalorParameters();
  DefineMaterials();
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fStepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
//  G4LogicalSkinSurface::CleanSurfaceTable();
//  G4LogicalBorderSurface::CleanSurfaceTable();
  ComputeCalorParameters();
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
//  env_mat = nist->FindOrBuildMaterial(env_mat_name);
  water_mat = nist->FindOrBuildMaterial("G4_WATER");
  polyethylene_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  polyprop_mat = nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
  cout<<"!!!!polyethylene name density "<<polyethylene_mat->GetName()<<" "<<polyethylene_mat->GetDensity()/(g/cm3)<<endl;
//  polyprop_mat = nist->FindOrBuildMaterial("G4_WATER");
//  polyprop_mat = nist->FindOrBuildMaterial("G4_Fe");
  cout<<"!!!!polypropilen name density "<<polyprop_mat->GetName()<<" "<<polyprop_mat->GetDensity()/(g/cm3)<<endl;
  alum_mat = nist->FindOrBuildMaterial("G4_Al");
//!!!!!!!!!!!!!!!!!  change aluminium on glass illuminator
//  alum_mat = nist->FindOrBuildMaterial("G4_GLASS_PLATE");
  cout<<"!!!!aluminium name density "<<alum_mat->GetName()<<" "<<alum_mat->GetDensity()/(g/cm3)<<endl;
  ferrum_mat = nist->FindOrBuildMaterial("G4_Fe");
  cout<<"!!!!ferrum name density "<<ferrum_mat->GetName()<<" "<<ferrum_mat->GetDensity()/(g/cm3)<<endl;
  world_mat = nist->FindOrBuildMaterial(world_mat_name);
  eye_mat = nist->FindOrBuildMaterial("G4_WATER");
  cout<<"!!!!eye name density "<<eye_mat->GetName()<<" "<<eye_mat->GetDensity()/(g/cm3)<<endl;
  cout<<"!!!!world name density "<<world_mat->GetName()<<" "<<world_mat->GetDensity()/(g/cm3)<<endl;
  testis_mat = nist->FindOrBuildMaterial("G4_WATER");
  spleen_mat = nist->FindOrBuildMaterial("G4_WATER");
  stomach_mat = nist->FindOrBuildMaterial("G4_WATER");
  brain_mat = nist->FindOrBuildMaterial("G4_WATER");

/* --------------- start sphere geometry ------------ */

  G4double world_Z = 1.03*aluminium_lengthZ;
  G4Tubs* solidWorld= new G4Tubs("World",            //its name
                       0, world_R, world_Z/2.0, 0, CLHEP::twopi);//size
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,        //its solid
                                       world_mat,   //its material
                                       "World");       //its name
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logicWorld,                            //its logical volume
                        "World",                                       //its name
                        0,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number

  G4double alum_Rmin = aluminium_R-alum_sizeZ;
  cout<<"!!!!!!! alum thickness "<<(aluminium_R-alum_Rmin)/cm<<endl;
  cout<<"!!!!!!! alum Rmin Rmax length "<<aluminium_R/cm<<" "<<alum_Rmin/cm<<" "<<aluminium_lengthZ/cm<<endl;
// aluminium wall
//  G4Sphere* solidAluminium =
//    new G4Sphere("Aluminium",                       //its name
//       alum_Rmin, aluminium_R, 0, 2.*pi, 0, pi);     //its size
//------------- start cylinder geometry
  G4Tubs* solidAluminium =
    new G4Tubs("Aluminium",                       //its name
       alum_Rmin, aluminium_R, aluminium_lengthZ/2.0, 0, CLHEP::twopi);     //its size
// subtraction illuminator
//  G4double ilum_R = 10.*cm;
//  G4double theta_illum = acos(ilum_R/aluminium_R);
//  G4Sphere* solidSubtr1 =
//    new G4Sphere("Subtr1",                       //its name
//       alum_Rmin, aluminium_R, 0, 2.*pi, pi-theta_illum, pi);     //its size
//  G4SubtractionSolid* solidAluminium1 = new G4SubtractionSolid("Aluminium1",
//                                        solidAluminium,
//                                        solidSubtr1);
// end subtraction
  logic_alum = new G4LogicalVolume(solidAluminium,        //its solid
                                       alum_mat,   //its material
                                       "Aluminium");       //its name
  physic_alum = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_alum,                            //its logical volume
                        "Aluminium",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_alum->CheckOverlaps();

  G4Tubs* solidAluminium1 =
    new G4Tubs("Aluminium1",                       //its name
       0, aluminium_R, alum_sizeZ/2.0, 0, CLHEP::twopi);     //its size
  G4LogicalVolume* logic_alum1 = new G4LogicalVolume(solidAluminium1,        //its solid
                                       alum_mat,   //its material
                                       "Aluminium1");       //its name
  G4VPhysicalVolume* physic_alum1 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0, 0, -(aluminium_lengthZ/2.0+alum_sizeZ/2.0) ),
                        logic_alum1,                            //its logical volume
                        "Aluminium1",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_alum1->CheckOverlaps();

  G4Tubs* solidAluminium2 =
    new G4Tubs("Aluminium2",                       //its name
       0, aluminium_R, alum_sizeZ/2.0, 0, CLHEP::twopi);     //its size
  G4LogicalVolume* logic_alum2 = new G4LogicalVolume(solidAluminium2,        //its solid
                                       alum_mat,   //its material
                                       "Aluminium2");       //its name
  G4VPhysicalVolume* physic_alum2 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0, 0, (aluminium_lengthZ/2.0+alum_sizeZ/2.0) ),
                        logic_alum2,                            //its logical volume
                        "Aluminium2",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_alum2->CheckOverlaps();

  G4double polyprop_Rmin = alum_Rmin-polyprop_sizeZ;
  G4double polyprop_Rmax = alum_Rmin;
  cout<<"!!!!!!! polyprop thickness "<<(polyprop_Rmax-polyprop_Rmin)/cm<<endl;
  G4Tubs* solidPolyprop =
    new G4Tubs("Polyprop",                       //its name
       polyprop_Rmin, polyprop_Rmax, (aluminium_lengthZ/4.0-polyprop_sizeZ), 0, CLHEP::twopi); //its size
  logic_polyprop = new G4LogicalVolume(solidPolyprop,        //its solid
                                       polyprop_mat,   //its material
                                       "Polyprop");       //its name
  physic_polyprop = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,-(aluminium_lengthZ/4.0-polyprop_sizeZ)),
                        logic_polyprop,                            //its logical volume
                        "Polyprop",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_polyprop->CheckOverlaps();

  G4Tubs* solidPolyprop01 =
    new G4Tubs("Polyprop01",                       //its name
       0, polyprop_Rmax, polyprop_sizeZ/2.0, 0, CLHEP::twopi);     //its size
  G4LogicalVolume* logic_polyprop01 = new G4LogicalVolume(solidPolyprop01,        //its solid
                                       polyprop_mat,   //its material
                                       "Polyprop01");       //its name
  G4VPhysicalVolume* physic_polyprop01 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0, 0, -(aluminium_lengthZ/2.0-polyprop_sizeZ/2.0) ),
                        logic_polyprop01,                            //its logical volume
                        "Polyprop01",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_polyprop01->CheckOverlaps();

  polyprop_Rmin = alum_Rmin-polyprop_sizeZ1;
  cout<<"!!!!!!! polyprop1 thickness "<<(polyprop_Rmax-polyprop_Rmin)/cm<<endl;
//  polyprop_Rmax = alum_Rmin;
  G4Tubs* solidPolyprop1 =
    new G4Tubs("Polyprop1",                       //its name
       polyprop_Rmin, polyprop_Rmax, (aluminium_lengthZ/4.0-polyprop_sizeZ1), 0, CLHEP::twopi); //its size
  logic_polyprop1 = new G4LogicalVolume(solidPolyprop1,        //its solid
                                       polyprop_mat,   //its material
                                       "Polyprop1");       //its name
  physic_polyprop1 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,(aluminium_lengthZ/4.0-polyprop_sizeZ1)),
                        logic_polyprop1,                            //its logical volume
                        "Polyprop1",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_polyprop1->CheckOverlaps();


  G4Tubs* solidPolyprop11 =
    new G4Tubs("Polyprop11",                       //its name
       0, polyprop_Rmax, polyprop_sizeZ1/2.0, 0, CLHEP::twopi);     //its size
  G4LogicalVolume* logic_polyprop11 = new G4LogicalVolume(solidPolyprop11,        //its solid
                                       polyprop_mat,   //its material
                                       "Polyprop11");       //its name
  G4VPhysicalVolume* physic_polyprop11 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0, 0, (aluminium_lengthZ/2.0-polyprop_sizeZ1/2.0) ),
                        logic_polyprop11,                            //its logical volume
                        "Polyprop11",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_polyprop11->CheckOverlaps();

// start man box
  G4double man_sizeX = 30.*cm;
  G4double man_sizeZ = 30.*cm;
  G4double man_sizeY = 150.*cm;
  cout<<"!!!!!!! man x y z "<<man_sizeX/cm<<" "<<man_sizeY/cm<<" "<<man_sizeZ/cm<<endl;
  G4Box* solidMan =
    new G4Box("Man",                       //its name
       0.5*man_sizeX, 0.5*man_sizeY, 0.5*man_sizeZ);     //its size
  G4LogicalVolume* logic_man = new G4LogicalVolume(solidMan,        //its solid
                                       water_mat,   //its material
                                       "Man");       //its name
//  G4VPhysicalVolume* physic_man = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,-(aluminium_lengthZ/2.0-polyprop_sizeZ-man_sizeZ) ),
//                        logic_man,                            //its logical volume
//                        "Man",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_man->CheckOverlaps();
// end man box


//------------- end cylinder geometry 

//    polypropylene wall
  
//  G4double polyprop_Rmin = alum_Rmin-polyprop_sizeZ;
//  G4double polyprop_Rmax = alum_Rmin;
//  cout<<"!!!!!!! polyprop thickness "<<(polyprop_Rmax-polyprop_Rmin)/cm<<endl;
//  G4Sphere* solidPolyprop =
//    new G4Sphere("Polyprop",                       //its name
//       polyprop_Rmin, polyprop_Rmax, 0, 2.*pi, pi/2., pi); //its size
// subtraction illuminator
//  G4Sphere* solidSubtr2 =
//    new G4Sphere("Subtr2",                       //its name
//       polyprop_Rmin, polyprop_Rmin, 0, 2.*pi, pi-theta_illum, pi);     //its size
//  G4SubtractionSolid* solidPolyprop1 = new G4SubtractionSolid("Polyprop1",
//                                        solidPolyprop,
//                                        solidSubtr2);
// end subtraction
//  logic_polyprop = new G4LogicalVolume(solidPolyprop,        //its solid
//                                       polyprop_mat,   //its material
//                                       "Polyprop");       //its name
//  physic_polyprop = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,0),
//                        logic_polyprop,                            //its logical volume
//                        "Polyprop",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_polyprop->CheckOverlaps();
//    polypropylene1 wall
//  G4double polyprop_Rmin1 = alum_Rmin-polyprop_sizeZ1;
//  G4double polyprop_Rmax1 = alum_Rmin;
//  G4Sphere* solidPolyprop1 =
//    new G4Sphere("Polyprop1",                       //its name
//       polyprop_Rmin1, polyprop_Rmax1, 0, 2.*pi, 0, pi/2.); //its size
//  logic_polyprop1 = new G4LogicalVolume(solidPolyprop1,        //its solid
//                                       polyprop_mat,   //its material
//                                       "Polyprop1");       //its name
//  physic_polyprop1 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,0),
//                        logic_polyprop1,                            //its logical volume
//                        "Polyprop1",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_polyprop1->CheckOverlaps();
  G4double Distance = 15.*cm;

//  composite tube
  G4double lifbox_sizeZ = lif_H + 2.0*mm;
  G4double lifbox_radius = lif_R + 1.0*mm;
  polyprop_Rmin = alum_Rmin-polyprop_sizeZ;
  if ( fCompDet == "yes" )
  {
// box for LiF sample
   polyprop_Rmin = alum_Rmin-polyprop_sizeZ;
   cout<<"!!!!!!!!!!!!!Lif box radius length "<<lifbox_radius/cm<<" "<<lifbox_sizeZ/cm<<endl;
   G4Tubs* solidBoxLiF =
    new G4Tubs("BLiF",                       //its name
       0., lifbox_radius, lifbox_sizeZ/2., 0., CLHEP::twopi);     //its size
   G4LogicalVolume* logic_blif = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF");       //its name
   G4VPhysicalVolume* physic_blif = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), -composite_sizeZ/2.),
                        G4ThreeVector(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), -aluminium_lengthZ/4.0),
                        logic_blif,                            //its logical volume
                        "BLiF",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif->CheckOverlaps();
// end box for LiF sample
// LiF sample
   G4Tubs* solidLiF =
    new G4Tubs("LiF",                       //its name
       0., lif_R, lif_H/2., 0., CLHEP::twopi);     //its size
   logic_lif = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF");       //its name
   physic_lif = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif,                            //its logical volume
                        "LiF",                                       //its name
                        logic_blif,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif->CheckOverlaps();

   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite1 =
    new G4Tubs("Composite1",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite1 = new G4LogicalVolume(solidComposite1,        //its solid
                                       composite_mat,   //its material
                                       "Composite1");       //its name
   G4ThreeVector rotAxis(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite1 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_composite1,                            //its logical volume
                        "Composite1",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite1->CheckOverlaps();
   G4LogicalVolume* logic_blif1 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF1");       //its name
   G4VPhysicalVolume* physic_blif1 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_blif1,                            //its logical volume
                        "BLiF1",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif1->CheckOverlaps();
   G4LogicalVolume* logic_lif1 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF1");       //its name
   G4VPhysicalVolume* physic_lif1 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif1,                            //its logical volume
                        "LiF1",                                       //its name
                        logic_blif1,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif1->CheckOverlaps();
   cout<<"lif1 name "<<physic_lif1->GetName()<<endl;
// disk1
   G4Tubs* solidDisk1 =
    new G4Tubs("Disk1",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk1 = new G4LogicalVolume(solidDisk1,        //its solid
                                       composite_mat,   //its material
                                       "Disk1");       //its name
   G4VPhysicalVolume* physic_disk1 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                     G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                     logic_disk1,                            //its logical volume
                     "Disk1",                                       //its name
                     logicWorld,                           //its mother  volume
                     false,                                  //no boolean operation
                     0);                                     //copy number
   physic_disk1->CheckOverlaps();
// disk2
   G4Tubs* solidDisk11 =
    new G4Tubs("Disk11",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk11 = new G4LogicalVolume(solidDisk11,        //its solid
                                       composite_mat,   //its material
                                       "Disk11");       //its name
   G4VPhysicalVolume* physic_disk11 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), (composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                      logic_disk11,                            //its logical volume
                      "Disk11",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk11->CheckOverlaps();
//--------------------------------------  multiple conteiners
   composite_thickness += 1.0*cm;
   composite_Rmax = composite_Rmin + composite_thickness;   
   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite2 =
    new G4Tubs("Composite2",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite2 = new G4LogicalVolume(solidComposite2,        //its solid
                                       composite_mat,   //its material
                                       "Composite2");       //its name
   rotAxis.set(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi*2.0/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite2 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_composite2,                            //its logical volume
                        "Composite2",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite2->CheckOverlaps();
   G4LogicalVolume* logic_blif2 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF2");       //its name
   G4VPhysicalVolume* physic_blif2 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_blif2,                            //its logical volume
                        "BLiF2",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif2->CheckOverlaps();
   G4LogicalVolume* logic_lif2 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF2");       //its name
   G4VPhysicalVolume* physic_lif2 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif2,                            //its logical volume
                        "LiF2",                                       //its name
                        logic_blif2,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif2->CheckOverlaps();
   cout<<"lif2 name "<<physic_lif2->GetName()<<endl;
// disk1
   G4Tubs* solidDisk2 =
    new G4Tubs("Disk2",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk2 = new G4LogicalVolume(solidDisk2,        //its solid
                                       composite_mat,   //its material
                                       "Disk2");       //its name
   G4VPhysicalVolume* physic_disk2 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                      logic_disk2,                            //its logical volume
                      "Disk2",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk2->CheckOverlaps();
// disk2
   G4Tubs* solidDisk22 =
    new G4Tubs("Disk22",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk22 = new G4LogicalVolume(solidDisk22,        //its solid
                                       composite_mat,   //its material
                                       "Disk22");       //its name
   G4VPhysicalVolume* physic_disk22 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), composite_thickness/2.0),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                      logic_disk22,                            //its logical volume
                      "Disk22",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk22->CheckOverlaps();

   composite_thickness += 1.0*cm;
   composite_Rmax = composite_Rmin + composite_thickness;   
   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite3 =
    new G4Tubs("Composite3",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite3 = new G4LogicalVolume(solidComposite3,        //its solid
                                       composite_mat,   //its material
                                       "Composite3");       //its name
   rotAxis.set(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi*3.0/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite3 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_composite3,                            //its logical volume
                        "Composite3",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite3->CheckOverlaps();
   G4LogicalVolume* logic_blif3 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF3");       //its name
   G4VPhysicalVolume* physic_blif3 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_blif3,                            //its logical volume
                        "BLiF3",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif3->CheckOverlaps();
   G4LogicalVolume* logic_lif3 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF3");       //its name
   G4VPhysicalVolume* physic_lif3 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif3,                            //its logical volume
                        "LiF3",                                       //its name
                        logic_blif3,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif3->CheckOverlaps();
   cout<<"lif3 name "<<physic_lif3->GetName()<<endl;
// disk1
   G4Tubs* solidDisk3 =
    new G4Tubs("Disk3",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk3 = new G4LogicalVolume(solidDisk3,        //its solid
                                       composite_mat,   //its material
                                       "Disk3");       //its name
   G4VPhysicalVolume* physic_disk3 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                      logic_disk3,                            //its logical volume
                      "Disk3",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk3->CheckOverlaps();
// disk2
   G4Tubs* solidDisk33 =
    new G4Tubs("Disk33",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk33 = new G4LogicalVolume(solidDisk33,        //its solid
                                       composite_mat,   //its material
                                       "Disk33");       //its name
   G4VPhysicalVolume* physic_disk33 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), composite_thickness/2.0),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                      logic_disk33,                            //its logical volume
                      "Disk33",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk33->CheckOverlaps();

   composite_thickness += 1.0*cm;
   composite_Rmax = composite_Rmin + composite_thickness;   
   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite4 =
    new G4Tubs("Composite4",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite4 = new G4LogicalVolume(solidComposite4,        //its solid
                                       composite_mat,   //its material
                                       "Composite4");       //its name
   rotAxis.set(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi*4.0/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite4 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_composite4,                            //its logical volume
                        "Composite4",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite4->CheckOverlaps();
   G4LogicalVolume* logic_blif4 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF4");       //its name
   G4VPhysicalVolume* physic_blif4 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -aluminium_lengthZ/4.),
                        logic_blif4,                            //its logical volume
                        "BLiF4",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif4->CheckOverlaps();
   G4LogicalVolume* logic_lif4 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF4");       //its name
   G4VPhysicalVolume* physic_lif4 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif4,                            //its logical volume
                        "LiF4",                                       //its name
                        logic_blif4,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif4->CheckOverlaps();
   cout<<"lif4 name "<<physic_lif4->GetName()<<endl;
// disk1
   G4Tubs* solidDisk4 =
    new G4Tubs("Disk4",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk4 = new G4LogicalVolume(solidDisk4,        //its solid
                                       composite_mat,   //its material
                                       "Disk4");       //its name
   G4VPhysicalVolume* physic_disk4 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                      logic_disk4,                            //its logical volume
                      "Disk4",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk4->CheckOverlaps();
//   cout<<"lif1 name "<<physic_lif1->GetName()<<endl;
// disk2
   G4Tubs* solidDisk44 =
    new G4Tubs("Disk44",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk44 = new G4LogicalVolume(solidDisk44,        //its solid
                                       composite_mat,   //its material
                                       "Disk44");       //its name
   G4VPhysicalVolume* physic_disk44 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), composite_thickness/2.0),
                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                        logic_disk44,                            //its logical volume
                        "Disk44",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_disk44->CheckOverlaps();

//--------------------------------- another part
   polyprop_Rmin = alum_Rmin-polyprop_sizeZ1;
   cout<<"!!!!!!!!!!!!!Lif box radius length "<<lifbox_radius/cm<<" "<<lifbox_sizeZ/cm<<endl;
   G4LogicalVolume* logic_blif5 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF5");       //its name
   G4VPhysicalVolume* physic_blif5 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), -composite_sizeZ/2.),
                        G4ThreeVector(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), aluminium_lengthZ/4.0),
                        logic_blif5,                            //its logical volume
                        "BLiF5",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif5->CheckOverlaps();
// end box for LiF sample
// LiF sample
   G4LogicalVolume* logic_lif5 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF5");       //its name
   G4VPhysicalVolume* physic_lif5 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif5,                            //its logical volume
                        "LiF5",                                       //its name
                        logic_blif5,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif5->CheckOverlaps();
   cout<<"lif5 name "<<physic_lif5->GetName()<<endl;

   composite_thickness = 1.0*cm;
   composite_Rmax = composite_Rmin + composite_thickness;   
   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite6 =
    new G4Tubs("Composite6",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite6 = new G4LogicalVolume(solidComposite6,        //its solid
                                       composite_mat,   //its material
                                       "Composite6");       //its name
   rotAxis.set(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi*1.0/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite6 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_composite6,                            //its logical volume
                        "Composite6",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite6->CheckOverlaps();
   G4LogicalVolume* logic_blif6 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF6");       //its name
   G4VPhysicalVolume* physic_blif6 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_blif6,                            //its logical volume
                        "BLiF6",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif6->CheckOverlaps();
   G4LogicalVolume* logic_lif6 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF6");       //its name
   G4VPhysicalVolume* physic_lif6 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif6,                            //its logical volume
                        "LiF6",                                       //its name
                        logic_blif6,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif6->CheckOverlaps();
   cout<<"lif6 name "<<physic_lif6->GetName()<<endl;
// disk1
   G4Tubs* solidDisk6 =
    new G4Tubs("Disk6",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk6 = new G4LogicalVolume(solidDisk6,        //its solid
                                       composite_mat,   //its material
                                       "Disk6");       //its name
   G4VPhysicalVolume* physic_disk6 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                      logic_disk6,                            //its logical volume
                      "Disk6",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk6->CheckOverlaps();
//   cout<<"lif1 name "<<physic_lif1->GetName()<<endl;
// disk2
   G4Tubs* solidDisk66 =
    new G4Tubs("Disk66",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk66 = new G4LogicalVolume(solidDisk66,        //its solid
                                       composite_mat,   //its material
                                       "Disk66");       //its name
   G4VPhysicalVolume* physic_disk66 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), composite_thickness/2.0),
                       G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                       logic_disk66,                            //its logical volume
                       "Disk66",                                       //its name
                       logicWorld,                           //its mother  volume
                       false,                                  //no boolean operation
                       0);                                     //copy number
   physic_disk66->CheckOverlaps();


   composite_thickness += 1.0*cm;
   composite_Rmax = composite_Rmin + composite_thickness;   
   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite7 =
    new G4Tubs("Composite7",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite7 = new G4LogicalVolume(solidComposite7,        //its solid
                                       composite_mat,   //its material
                                       "Composite7");       //its name
   rotAxis.set(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi*2.0/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite7 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_composite7,                            //its logical volume
                        "Composite7",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite7->CheckOverlaps();
   G4LogicalVolume* logic_blif7 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF7");       //its name
   G4VPhysicalVolume* physic_blif7 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_blif7,                            //its logical volume
                        "BLiF7",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif7->CheckOverlaps();
   G4LogicalVolume* logic_lif7 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF7");       //its name
   G4VPhysicalVolume* physic_lif7 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif7,                            //its logical volume
                        "LiF7",                                       //its name
                        logic_blif7,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif7->CheckOverlaps();
   cout<<"lif7 name "<<physic_lif7->GetName()<<endl;
// disk1
   G4Tubs* solidDisk7 =
    new G4Tubs("Disk7",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk7 = new G4LogicalVolume(solidDisk7,        //its solid
                                       composite_mat,   //its material
                                       "Disk7");       //its name
   G4VPhysicalVolume* physic_disk7 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                      logic_disk7,                            //its logical volume
                      "Disk7",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk7->CheckOverlaps();
//   cout<<"lif1 name "<<physic_lif1->GetName()<<endl;
// disk2
   G4Tubs* solidDisk77 =
    new G4Tubs("Disk77",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk77 = new G4LogicalVolume(solidDisk77,        //its solid
                                       composite_mat,   //its material
                                       "Disk77");       //its name
   G4VPhysicalVolume* physic_disk77 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), composite_thickness/2.0),
                       G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                       logic_disk77,                            //its logical volume
                       "Disk77",                                       //its name
                       logicWorld,                           //its mother  volume
                       false,                                  //no boolean operation
                       0);                                     //copy number
   physic_disk77->CheckOverlaps();


   composite_thickness += 1.0*cm;
   composite_Rmax = composite_Rmin + composite_thickness;   
   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite8 =
    new G4Tubs("Composite8",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite8 = new G4LogicalVolume(solidComposite8,        //its solid
                                       composite_mat,   //its material
                                       "Composite8");       //its name
   rotAxis.set(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi*3.0/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite8 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_composite8,                            //its logical volume
                        "Composite8",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite8->CheckOverlaps();
   G4LogicalVolume* logic_blif8 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF8");       //its name
   G4VPhysicalVolume* physic_blif8 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_blif8,                            //its logical volume
                        "BLiF8",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif8->CheckOverlaps();
   G4LogicalVolume* logic_lif8 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF8");       //its name
   G4VPhysicalVolume* physic_lif8 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif8,                            //its logical volume
                        "LiF8",                                       //its name
                        logic_blif8,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif8->CheckOverlaps();
   cout<<"lif8 name "<<physic_lif8->GetName()<<endl;
// disk1
   G4Tubs* solidDisk8 =
    new G4Tubs("Disk8",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk8 = new G4LogicalVolume(solidDisk8,        //its solid
                                       composite_mat,   //its material
                                       "Disk8");       //its name
   G4VPhysicalVolume* physic_disk8 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                      logic_disk8,                            //its logical volume
                      "Disk8",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk8->CheckOverlaps();
//   cout<<"lif1 name "<<physic_lif1->GetName()<<endl;
// disk2
   G4Tubs* solidDisk88 =
    new G4Tubs("Disk88",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk88 = new G4LogicalVolume(solidDisk88,        //its solid
                                       composite_mat,   //its material
                                       "Disk88");       //its name
   G4VPhysicalVolume* physic_disk88 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), composite_thickness/2.0),
                       G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                       logic_disk88,                            //its logical volume
                       "Disk88",                                       //its name
                       logicWorld,                           //its mother  volume
                       false,                                  //no boolean operation
                       0);                                     //copy number
   physic_disk88->CheckOverlaps();


   composite_thickness += 1.0*cm;
   composite_Rmax = composite_Rmin + composite_thickness;   
   cout<<"!!!!!!!!!!composite rmax rmin length thick "<<composite_Rmax/cm<<" "<<composite_Rmin/cm<<" "<<composite_sizeZ/cm
              <<" "<<composite_thickness/cm<<endl;
   G4Tubs* solidComposite9 =
    new G4Tubs("Composite9",                       //its name
       composite_Rmin, composite_Rmax, composite_sizeZ/2., 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_composite9 = new G4LogicalVolume(solidComposite9,        //its solid
                                       composite_mat,   //its material
                                       "Composite9");       //its name
   rotAxis.set(0.0, -(polyprop_Rmin-composite_Rmax-1.*mm), 0.0 );
   rotAxis.rotateZ(CLHEP::twopi*4.0/5.0);
   cout<<"!!!!!!!!rotAxis x y z "<<rotAxis.x()<<" "<<rotAxis.y()<<" "<<rotAxis.z()<<endl;;
   G4VPhysicalVolume* physic_composite9 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_composite9,                            //its logical volume
                        "Composite9",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite9->CheckOverlaps();
   G4LogicalVolume* logic_blif9 = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF9");       //its name
   G4VPhysicalVolume* physic_blif9 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(rotAxis.x(),rotAxis.y(), -composite_sizeZ/2.),
                        G4ThreeVector(rotAxis.x(),rotAxis.y(), aluminium_lengthZ/4.),
                        logic_blif9,                            //its logical volume
                        "BLiF9",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif9->CheckOverlaps();
   G4LogicalVolume* logic_lif9 = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF9");       //its name
   G4VPhysicalVolume* physic_lif9 = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif9,                            //its logical volume
                        "LiF9",                                       //its name
                        logic_blif9,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif9->CheckOverlaps();
   cout<<"lif9 name "<<physic_lif9->GetName()<<endl;
// disk1
   G4Tubs* solidDisk9 =
    new G4Tubs("Disk9",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk9 = new G4LogicalVolume(solidDisk9,        //its solid
                                       composite_mat,   //its material
                                       "Disk9");       //its name
   G4VPhysicalVolume* physic_disk9 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), -(composite_sizeZ+composite_thickness/2.0)),
                      G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.+composite_sizeZ/2.+composite_thickness/2.0)),
                      logic_disk9,                            //its logical volume
                      "Disk9",                                       //its name
                      logicWorld,                           //its mother  volume
                      false,                                  //no boolean operation
                      0);                                     //copy number
   physic_disk9->CheckOverlaps();
//   cout<<"lif1 name "<<physic_lif1->GetName()<<endl;
// disk2
   G4Tubs* solidDisk99 =
    new G4Tubs("Disk99",                       //its name
       0, composite_Rmax, composite_thickness/2.0, 0, CLHEP::twopi);     //its size
   G4LogicalVolume* logic_disk99 = new G4LogicalVolume(solidDisk99,        //its solid
                                       composite_mat,   //its material
                                       "Disk99");       //its name
   G4VPhysicalVolume* physic_disk99 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(rotAxis.x(), rotAxis.y(), composite_thickness/2.0),
                       G4ThreeVector(rotAxis.x(), rotAxis.y(), (aluminium_lengthZ/4.-composite_sizeZ/2.-composite_thickness/2.0)),
                       logic_disk99,                            //its logical volume
                       "Disk99",                                       //its name
                       logicWorld,                           //its mother  volume
                       false,                                  //no boolean operation
                       0);                                     //copy number
   physic_disk99->CheckOverlaps();


//----------------------------------  end of multiple containers

  }
  else 
  {

// box for LiF sample
   cout<<"Lif box radius length "<<lifbox_radius/cm<<" "<<lifbox_sizeZ/cm<<endl;
   G4Tubs* solidBoxLiF =
    new G4Tubs("BLiF",                       //its name
       0., lifbox_radius, lifbox_sizeZ/2., 0., CLHEP::twopi);     //its size
   G4LogicalVolume* logic_blif = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF");       //its name
   G4ThreeVector rotAxis(1.0, 0.0, 0.0 );
   G4RotationMatrix *rotMatrix1 = new G4RotationMatrix(rotAxis, -pi/2.);
//   G4VPhysicalVolume* physic_blif = new G4PVPlacement(rotMatrix1,  //rotation
//                        G4ThreeVector(0,0,-(polyprop_Rmin-Distance)),
   G4VPhysicalVolume* physic_blif = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(0, 0, 0),
                        G4ThreeVector(0, -(composite_Rmin-lifbox_radius-1.*mm), -lifbox_sizeZ/2.0),
                        logic_blif,                            //its logical volume
                        "BLiF",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_blif->CheckOverlaps();
// end box for LiF sample
// LiF sample
   G4Tubs* solidLiF =
    new G4Tubs("LiF",                       //its name
       0., lif_R, lif_H/2., 0., 2.0*pi);     //its size
   logic_lif = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF");       //its name
//  G4ThreeVector rotAxis(1.0, 0.0, 0.0 );
//  G4RotationMatrix *rotMatrix1 = new G4RotationMatrix(rotAxis, -pi/2.);
   physic_lif = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,-(world_R-25.*cm)),
                        G4ThreeVector(0,0,0),
                        logic_lif,                            //its logical volume
                        "LiF",                                       //its name
//                        logicWorld,                           //its mother  volume
                        logic_blif,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_lif->CheckOverlaps();
// end LiF
//   G4VPhysicalVolume* physic_blif2 = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(0, 0, 0),
//                        logic_blif,                            //its logical volume
//                        "BLiF2",                                       //its name
//                        solidComposite,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//   physic_blif2->CheckOverlaps();
//   G4VPhysicalVolume* physic_lif2 = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,0),
//                        logic_lif,                            //its logical volume
//                        "LiF2",                                       //its name
//                        logic_blif2,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//   physic_lif2->CheckOverlaps();


  }
// curtain with ro = 6.7g/cm2
  G4double curtain_sizeXY = 30.*cm;
  G4double curtain_sizeZ = 6.7*cm;
  G4Box* solidCurtain =
    new G4Box("Curtain",                       //its name
       0.5*curtain_sizeXY, 0.5*curtain_sizeXY, 0.5*curtain_sizeZ);     //its size
  G4LogicalVolume* logic_curtain = new G4LogicalVolume(solidCurtain,        //its solid
                                       water_mat,   //its material
                                       "Curtain");       //its name
//  G4VPhysicalVolume* physic_curtain = new G4PVPlacement(0,  //rotation
//                        G4ThreeVector(0,0,-(polyprop_Rmin-Distance+2.*cm+curtain_sizeZ)),
//                        logic_curtain,                            //its logical volume
//                        "Curtain",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//   physic_curtain->CheckOverlaps();
// end curtain
//  water box
  G4Box* solidWater =
    new G4Box("Water",                       //its name
       0.5*water_sizeXY, 0.5*water_sizeXY, 0.5*water_sizeZ);     //its size
  logic_water = new G4LogicalVolume(solidWater,        //its solid
                                       water_mat,   //its material
                                       "Water");       //its name
//  physic_water = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,0.5*water_sizeZ),
//                        logic_water,                            //its logical volume
//                        "Water",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_water->CheckOverlaps();
// end water box
// eye sample
  G4double eye_posZ = -water_sizeZ/2.+0.5*sample_sizeZ;
  G4Box* solidEye =
    new G4Box("Eye",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_eye = new G4LogicalVolume(solidEye,        //its solid
                                       eye_mat,   //its material
                                       "Eye");       //its name
//  physic_eye = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,eye_posZ),
//                        logic_eye,                            //its logical volume
//                        "Eye",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_eye->CheckOverlaps();
// end eye
// testis sample
  G4double testis_posZ = -water_sizeZ/2.+2.*cm+0.5*sample_sizeZ;
  G4Box* solidTestis =
    new G4Box("Testis",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_testis = new G4LogicalVolume(solidTestis,        //its solid
                                       testis_mat,   //its material
                                       "Testis");       //its name
//  physic_testis = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,testis_posZ),
//                        logic_testis,                            //its logical volume
//                        "Testis",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_testis->CheckOverlaps();
// end testis
// spleen sample
  G4double spleen_posZ = -water_sizeZ/2.+5.*cm+0.5*sample_sizeZ;
  G4Box* solidSpleen =
    new G4Box("Spleen",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_spleen = new G4LogicalVolume(solidSpleen,        //its solid
                                       spleen_mat,   //its material
                                       "Spleen");       //its name
//  physic_spleen = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,spleen_posZ),
//                        logic_spleen,                            //its logical volume
//                        "Spleen",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_spleen->CheckOverlaps();
// end spleen
// brain sample
  G4double brain_posZ = -water_sizeZ/2.+7.*cm+0.5*sample_sizeZ;
  G4Box* solidBrain =
    new G4Box("Brain",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_brain = new G4LogicalVolume(solidBrain,        //its solid
                                       brain_mat,   //its material
                                       "Brain");       //its name
//  physic_brain = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,brain_posZ),
//                        logic_brain,                            //its logical volume
//                        "Brain",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_brain->CheckOverlaps();
// end brain
// stomach sample
  G4double stomach_posZ = -water_sizeZ/2.+9.*cm+0.5*sample_sizeZ;
  G4Box* solidStomach =
    new G4Box("Stomach",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_stomach = new G4LogicalVolume(solidStomach,        //its solid
                                       stomach_mat,   //its material
                                       "Stomach");       //its name
//  physic_stomach = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,stomach_posZ),
//                        logic_stomach,                            //its logical volume
//                        "Stomach",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_stomach->CheckOverlaps();
// end stomach

/* --------------- end sphere geometry --------------------- */


  
/* ------------------ start -------------------------- *\
  G4Sphere* solidWorld= new G4Sphere("World",            //its name
                       0, world_R, 0, 2.*pi, 0, pi);//size
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,        //its solid
                                       world_mat,   //its material
                                       "World");       //its name
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logicWorld,                            //its logical volume
                        "World",                                       //its name
                        0,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
//  water sphere
  cout<<"water rmax "<<water_R/cm<<endl;
  G4Sphere* solidWater = new G4Sphere("Water",                       //its name
        0, water_R, 0, 2.*pi, 0, pi);     //its size
  logic_water = new G4LogicalVolume(solidWater,        //its solid
                                       water_mat,   //its material
                                       "Water");       //its name
  physic_water = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_water,                            //its logical volume
                        "Water",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_water->CheckOverlaps();

// eye sample
  G4double eye_Rmin = water_R-sample_sizeZ;
  cout<<"eye rmin rmax "<<eye_Rmin/cm<<" "<<water_R/cm<<endl;
  G4Sphere* solidEye =
    new G4Sphere("Eye",                       //its name
       eye_Rmin, water_R, 0, 2.*pi, 0, pi);     //its size
  logic_eye = new G4LogicalVolume(solidEye,        //its solid
                                       eye_mat,   //its material
                                       "Eye");       //its name
  physic_eye = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_eye,                            //its logical volume
                        "Eye",                                       //its name
                        logic_water,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_eye->CheckOverlaps();
// end eye sample

// testis sample
  G4double testis_Rmin = water_R-2.*cm-sample_sizeZ;
  G4double testis_Rmax = water_R-2.*cm;
  cout<<"testis rmin rmax "<<testis_Rmin/cm<<" "<<testis_Rmax/cm<<endl;
  G4Sphere* solidTestis =
    new G4Sphere("Testis",                       //its name
       testis_Rmin, testis_Rmax, 0, 2.*pi, 0, pi);     //its size
  logic_testis = new G4LogicalVolume(solidTestis,        //its solid
                                       testis_mat,   //its material
                                       "Testis");       //its name
  physic_testis = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_testis,                            //its logical volume
                        "Testis",                                       //its name
                        logic_water,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_testis->CheckOverlaps();
// end testis
// spleen sample
  G4double spleen_Rmin = water_R-5.*cm-sample_sizeZ;
  G4double spleen_Rmax = water_R-5.*cm;
  cout<<"spleen rmin rmax "<<spleen_Rmin/cm<<" "<<spleen_Rmax/cm<<endl;
  G4Sphere* solidSpleen =
    new G4Sphere("Spleen",                       //its name
       spleen_Rmin, spleen_Rmax, 0, 2.*pi, 0, pi);     //its size
  logic_spleen = new G4LogicalVolume(solidSpleen,        //its solid
                                       spleen_mat,   //its material
                                       "Spleen");       //its name
  physic_spleen = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_spleen,                            //its logical volume
                        "Spleen",                                       //its name
                        logic_water,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_spleen->CheckOverlaps();
// end spleen
// brain sample
  G4double brain_Rmin = water_R-7.*cm-sample_sizeZ;
  G4double brain_Rmax = water_R-7.*cm;
  cout<<"brain rmin rmax "<<brain_Rmin/cm<<" "<<brain_Rmax/cm<<endl;
  G4Sphere* solidBrain =
    new G4Sphere("Brain",                       //its name
       brain_Rmin, brain_Rmax, 0, 2.*pi, 0, pi);     //its size
  logic_brain = new G4LogicalVolume(solidBrain,        //its solid
                                       brain_mat,   //its material
                                       "Brain");       //its name
  physic_brain = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_brain,                            //its logical volume
                        "Brain",                                       //its name
                        logic_water,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_brain->CheckOverlaps();
// end brain
// stomach sample
  G4double stomach_Rmin = water_R-9.*cm-sample_sizeZ;
  G4double stomach_Rmax = water_R-9.*cm;
  cout<<"stomach rmin rmax "<<stomach_Rmin/cm<<" "<<stomach_Rmax/cm<<endl;
  G4Sphere* solidStomach =
    new G4Sphere("Stomach",                       //its name
       stomach_Rmin, stomach_Rmax, 0, 2.*pi, 0, pi);     //its size
  logic_stomach = new G4LogicalVolume(solidStomach,        //its solid
                                       stomach_mat,   //its material
                                       "Stomach");       //its name
  physic_stomach = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_stomach,                            //its logical volume
                        "Stomach",                                       //its name
                        logic_water,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_stomach->CheckOverlaps();
// end stomach

  G4double composite_Rmin = water_R+gap_sizeZ;
  G4double composite_Rmax = water_R+gap_sizeZ+composite_sizeZ;
  if ( fCompDet == "yes" )
  {
//  composite box
   G4Sphere* solidComposite =
    new G4Sphere("Composite",                       //its name
       composite_Rmin, composite_Rmax, 0, 2.*pi, 0, pi);     //its size
   logic_composite = new G4LogicalVolume(solidComposite,        //its solid
                                       composite_mat,   //its material
                                       "Composite");       //its name
   physic_composite = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_composite,                            //its logical volume
                        "Composite",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite->CheckOverlaps();
  }
//    polypropylene wall
  G4double polyprop_Rmin = composite_Rmax+gap_sizeZ;
  G4double polyprop_Rmax = composite_Rmax+gap_sizeZ+polyprop_sizeZ;
  G4Sphere* solidPolyprop =
    new G4Sphere("Polyprop",                       //its name
       polyprop_Rmin, polyprop_Rmax, 0, 2.*pi, 0, pi); //its size
  logic_polyprop = new G4LogicalVolume(solidPolyprop,        //its solid
                                       polyprop_mat,   //its material
                                       "Polyprop");       //its name
  physic_polyprop = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,-(composite_r_max+gap_sizeZ+0.5*polyprop_sizeZ)),
                        G4ThreeVector(0,0,0),
                        logic_polyprop,                            //its logical volume
                        "Polyprop",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_polyprop->CheckOverlaps();
//       aluminium wall
  G4double alum_Rmin = polyprop_Rmax;
  G4double alum_Rmax = polyprop_Rmax+alum_sizeZ;
  cout<<"!!!!!!!!!!!!!!!!!!!!aluminium R "<<alum_Rmax/cm<<endl;
  G4Sphere* solidAluminium =
    new G4Sphere("Aluminium",                       //its name
       alum_Rmin, alum_Rmax, 0, 2.*pi, 0, pi);     //its size
  logic_alum = new G4LogicalVolume(solidAluminium,        //its solid
                                       alum_mat,   //its material
                                       "Aluminium");       //its name
  physic_alum = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,-(composite_r_max+gap_sizeZ+polyprop_sizeZ+0.5*alum_sizeZ)),
//                        G4ThreeVector(0,0,-(composite_sizeZ+2*gap_sizeZ+polyprop_sizeZ+0.5*alum_sizeZ)),
                        G4ThreeVector(0,0,0),
                        logic_alum,                            //its logical volume
                        "Aluminium",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_alum->CheckOverlaps();
\* -------------------------- end ------------------- */

/* -------------------------- usual geometry start ------------ *\
  cout<<"world xy z "<<world_sizeXY<<" "<<world_sizeZ<<endl;
  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0);                     //copy number

//  G4double lifbox_sizeZ = lif_H + 2.0*mm;
//  G4double lifbox_radius = lif_R + 1.0*mm;
//  cout<<"Lif box radius length "<<lifbox_radius/cm<<" "<<lifbox_sizeZ/cm<<endl;
//  G4Tubs* solidBoxLiF =
//    new G4Tubs("BLiF",                       //its name
//       0., lifbox_radius, lifbox_sizeZ/2., 0., CLHEP::twopi);     //its size
//  G4LogicalVolume* logic_blif = new G4LogicalVolume(solidBoxLiF,        //its solid
//                                       polyprop_mat,   //its material
//                                       "BLiF");       //its name
//  G4ThreeVector rotAxis(1.0, 0.0, 0.0 );
//  G4RotationMatrix *rotMatrix1 = new G4RotationMatrix(rotAxis, -pi/2.);
//  G4VPhysicalVolume* physic_blif = new G4PVPlacement(rotMatrix1,  //rotation
//                        G4ThreeVector(0, 0, lifbox_radius),
//                        logic_blif,                            //its logical volume
//                        "BLiF",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//   physic_blif->CheckOverlaps();
// end box for LiF sample
// LiF sample
//   G4Tubs* solidLiF =
//    new G4Tubs("LiF",                       //its name
//       0., lif_R, lif_H/2., 0., 2.0*pi);     //its size
//   logic_lif = new G4LogicalVolume(solidLiF,        //its solid
//                                       lif_mat,   //its material
//                                       "LiF");       //its name
//   physic_lif = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,0),
//                        logic_lif,                            //its logical volume
//                        "LiF",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        logic_blif,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//   physic_lif->CheckOverlaps();
// end LiF

  G4double lifbox_sizeZ = sample_sizeZ + 2.0*mm;
  G4double lifbox_sizeXY = sample_sizeXY + 2.0*mm;
  cout<<"Lif box xy length "<<lifbox_sizeXY/cm<<" "<<lifbox_sizeZ/cm<<endl;
  G4Box* solidBoxLiF =
    new G4Box("BLiF",                       //its name
       0.5*lifbox_sizeXY, 0.5*lifbox_sizeXY, 0.5*lifbox_sizeZ);     //its size
  G4LogicalVolume* logic_blif = new G4LogicalVolume(solidBoxLiF,        //its solid
                                       polyprop_mat,   //its material
                                       "BLiF");       //its name
  G4VPhysicalVolume* physic_blif = new G4PVPlacement(0,  //rotation
                        G4ThreeVector(0, 0, lifbox_sizeZ/2.0),
                        logic_blif,                            //its logical volume
                        "BLiF",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number


  cout<<"!!!! LiF size x y z "<<sample_sizeXY<<" "<<sample_sizeXY<<" "<<sample_sizeZ<<endl;
  G4Box* solidLiF =
    new G4Box("LiF",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_lif = new G4LogicalVolume(solidLiF,        //its solid
                                       lif_mat,   //its material
                                       "LiF");       //its name
  physic_lif = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,0),
                        logic_lif,                            //its logical volume
                        "LiF",                                       //its name
                        logic_blif,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_lif->CheckOverlaps();
// end lif



//  water box
  cout<<"!!!! water size x y z "<<water_sizeXY<<" "<<water_sizeXY<<" "<<water_sizeZ<<endl;
  G4Box* solidWater =
    new G4Box("Water",                       //its name
       0.5*water_sizeXY, 0.5*water_sizeXY, 0.5*water_sizeZ);     //its size
  logic_water = new G4LogicalVolume(solidWater,        //its solid
                                       water_mat,   //its material
                                       "Water");       //its name
//  physic_water = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,0.5*water_sizeZ),
//                        logic_water,                            //its logical volume
//                        "Water",                                       //its name
//                        logicWorld,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_water->CheckOverlaps();

//  G4double sample_sizeZ = 0.5*cm;
// eye sample
//  cout<<"!!!! eye size x y z "<<sample_sizeXY<<" "<<sample_sizeXY<<" "<<sample_sizeZ<<endl;
  G4double eye_posZ = -water_sizeZ/2.+0.5*sample_sizeZ;
  G4Box* solidEye =
    new G4Box("Eye",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_eye = new G4LogicalVolume(solidEye,        //its solid
                                       eye_mat,   //its material
                                       "Eye");       //its name
//  physic_eye = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,eye_posZ),
//                        logic_eye,                            //its logical volume
//                        "Eye",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_eye->CheckOverlaps();
// end eye
// testis sample
  G4double testis_posZ = -water_sizeZ/2.+2.*cm+0.5*sample_sizeZ;
  G4Box* solidTestis =
    new G4Box("Testis",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_testis = new G4LogicalVolume(solidTestis,        //its solid
                                       testis_mat,   //its material
                                       "Testis");       //its name
//  physic_testis = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,testis_posZ),
//                        logic_testis,                            //its logical volume
//                        "Testis",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_testis->CheckOverlaps();
// end testis
// spleen sample
  G4double spleen_posZ = -water_sizeZ/2.+5.*cm+0.5*sample_sizeZ;
  G4Box* solidSpleen =
    new G4Box("Spleen",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_spleen = new G4LogicalVolume(solidSpleen,        //its solid
                                       spleen_mat,   //its material
                                       "Spleen");       //its name
//  physic_spleen = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,spleen_posZ),
//                        logic_spleen,                            //its logical volume
//                        "Spleen",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_spleen->CheckOverlaps();
// end spleen
// brain sample
  G4double brain_posZ = -water_sizeZ/2.+7.*cm+0.5*sample_sizeZ;
  G4Box* solidBrain =
    new G4Box("Brain",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_brain = new G4LogicalVolume(solidBrain,        //its solid
                                       brain_mat,   //its material
                                       "Brain");       //its name
//  physic_brain = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,brain_posZ),
//                        logic_brain,                            //its logical volume
//                        "Brain",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_brain->CheckOverlaps();
// end brain
// brain stomach
  G4double stomach_posZ = -water_sizeZ/2.+9.*cm+0.5*sample_sizeZ;
  G4Box* solidStomach =
    new G4Box("Stomach",                       //its name
       0.5*sample_sizeXY, 0.5*sample_sizeXY, 0.5*sample_sizeZ);     //its size
  logic_stomach = new G4LogicalVolume(solidStomach,        //its solid
                                       stomach_mat,   //its material
                                       "Stomach");       //its name
//  physic_stomach = new G4PVPlacement(0,                           //rotation
//                        G4ThreeVector(0,0,stomach_posZ),
//                        logic_stomach,                            //its logical volume
//                        "Stomach",                                       //its name
//                        logic_water,                           //its mother  volume
//                        false,                                  //no boolean operation
//                        0);                                     //copy number
//  physic_stomach->CheckOverlaps();
// end stomach


  if ( fCompDet == "yes" )
  {
//  composite box
   cout<<"composite z "<<-(gap_sizeZ)/cm<<endl;
   cout<<"composite xy thickness "<<sizeXY<<" "<<composite_thickness<<endl;
   G4Box* solidComposite =
    new G4Box("Composite",                       //its name
       0.5*sizeXY, 0.5*sizeXY, 0.5*composite_thickness);     //its size
   logic_composite = new G4LogicalVolume(solidComposite,        //its solid
                                       composite_mat,   //its material
                                       "Composite");       //its name
   physic_composite = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,-(0.5*composite_thickness+gap_sizeZ)),
                        logic_composite,                            //its logical volume
                        "Composite",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
   physic_composite->CheckOverlaps();
  }


//    polypropylene wall
  cout<<"!!!! polyprop size x y z "<<sizeXY<<" "<<sizeXY<<" "<<polyprop_sizeZ<<endl;
  cout<<"polypropilene z "<<-(composite_thickness+gap_sizeZ)/cm<<endl;
  G4Box* solidPolyprop =
    new G4Box("Polyprop",                       //its name
       0.5*sizeXY, 0.5*sizeXY, 0.5*polyprop_sizeZ);     //its size
  logic_polyprop = new G4LogicalVolume(solidPolyprop,        //its solid
                                       polyprop_mat,   //its material
                                       "Polyprop");       //its name
  physic_polyprop = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,-(composite_thickness+gap_sizeZ+0.5*polyprop_sizeZ)),
                        logic_polyprop,                            //its logical volume
                        "Polyprop",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_polyprop->CheckOverlaps();
//       aluminium wall
  aluminium_positionZ = -(composite_thickness+gap_sizeZ+polyprop_sizeZ+0.5*alum_sizeZ);
  cout<<"!!!! aluminium size x y z position "<<sizeXY<<" "<<sizeXY<<" "<<alum_sizeZ<<" "<<aluminium_positionZ<<endl;
  cout<<"aluminium z "<<-(composite_thickness+gap_sizeZ+polyprop_sizeZ)/cm<<endl;
  G4Box* solidAluminium =
    new G4Box("Aluminium",                       //its name
       0.5*sizeXY, 0.5*sizeXY, 0.5*alum_sizeZ);     //its size
  logic_alum = new G4LogicalVolume(solidAluminium,        //its solid
                                       alum_mat,   //its material
                                       "Aluminium");       //its name
  physic_alum = new G4PVPlacement(0,                           //rotation
                        G4ThreeVector(0,0,-(composite_thickness+gap_sizeZ+polyprop_sizeZ+0.5*alum_sizeZ)),
                        logic_alum,                            //its logical volume
                        "Aluminium",                                       //its name
                        logicWorld,                           //its mother  volume
                        false,                                  //no boolean operation
                        0);                                     //copy number
  physic_alum->CheckOverlaps();

\* --------------------- end usual geometry ---------------- */


  //
  // Envelope
  //
//  G4Box* solidEnv =
//    new G4Box("Envelope",                    //its name
//        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

//  logicEnv =
//    new G4LogicalVolume(solidEnv,            //its solid
//                        env_mat,             //its material
//                        "Envelope");         //its name

//  physEnv = new G4PVPlacement(0,                       //no rotation
//                    G4ThreeVector(),         //at (0,0,0)
//                    logicEnv,                //its logical volume
//                    "Envelope",              //its name
//                    logicWorld,              //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking






//  ConstructScoringVolumes();



// Set maximal step size
  //G4double maxStep = 0.1 * CLHEP::mm;
  //fStepLimit = new G4UserLimits(maxStep);
  //logicWorld->SetUserLimits(fStepLimit);
  //logicEnv->SetUserLimits(fStepLimit);
  //logictempl->SetUserLimits(fStepLimit);
  //logictemp_pl_plate->SetUserLimits(fStepLimit);
  //logicShape2->SetUserLimits(fStepLimit);
  //
  //always return the physical World
  //
  return physWorld;
}

void DetectorConstruction::SetWallThickness(G4double z)
{
 alum_sizeZ = z;
 aluminium_positionZ = (composite_thickness+gap_sizeZ+polyprop_sizeZ+0.5*alum_sizeZ);
// G4RunManager::GetRunManager()->GeometryHasBeenModified();
// G4ThreeVector xx = physic_alum->GetTranslation();
// G4double zz = dynamic_cast<G4Box*>(logic_alum->GetSolid())->GetZHalfLength()*2.0;
// cout<<"!!!!!!!!!!! alum positionz sizeZ "<<xx.z()<<" "<<zz<<endl;
}


void DetectorConstruction::DefineMaterials()
{
  materials = Materials::GetInstance();
  composite_mat = FindMaterial("Comp");
  lif_mat = FindMaterial("LiF");
//  eye_mat = FindMaterial("LiF");
//  world_mat = FindMaterial("WATER1");
//  water_mat = FindMaterial("G4_Water");
//  defaultMaterial  = FindMaterial("Galactic");
//  AlveolaMaterial = FindMaterial("AIR");
//  GrooveMaterial = FindMaterial("AIR");
//  FiberMaterial = FindMaterial("WLS");
//  Clad1Material = FindMaterial("CLAD1");
//  Clad2Material = FindMaterial("CLAD2");
//  PMTMaterial  = FindMaterial("GLASS");
//  CaloMaterial = FindMaterial("SCINT");
//  FeMaterial = FindMaterial("FE");
//  MirrorMaterial = FindMaterial("Mylar");
}

void DetectorConstruction::ConstructScoringVolumes(){


//  G4NistManager* nist = G4NistManager::Instance();

//========platic A-150-tissue  with incorporated Bi particles =======================
//  G4double density =  admix_c + (rho_tissue/rho_admix)*(rho_admix - admix_c);
//  G4double admix_fr = admix_c /density;// mass fractions
//  G4int ncomp = 2;

//  TissueWithAdmixture = new G4Material(TissueWithAdmixture_name, density, ncomp);
  //Tissue = nist->FindOrBuildMaterial(Tissue_name); //for standard NIST materials

//  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  ///material->DefineMaterials();
//  Tissue  = material->GetMaterial(Tissue_name);

//  Admixture = nist->FindOrBuildMaterial(Admixture_name);
//  TissueWithAdmixture->AddMaterial(Tissue, 1 - admix_fr);
//  TissueWithAdmixture->AddMaterial(Admixture, admix_fr);
  //============================================================================

  //


//  G4int  N = int(env_sizeZ/sc_vol_st);

//  G4Box* templatebox =  new G4Box("Layer", 0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*sc_vol_st);     //its size

//  logictempl = new G4LogicalVolume(templatebox,          //its solid
//                        nist->FindOrBuildMaterial(env_mat_name),           //its material
//                        "Layer");
//  G4String temp_str;
//  const G4String pref = "Layer";

//  G4VPhysicalVolume* phys_vol;

//  G4int Npos = int((pos2.z() + 0.5*env_sizeZ)/sc_vol_st);

// find position in local volume

//  G4double pos_log_z = pos2.z() + 0.5*env_sizeZ - Npos*sc_vol_st;

// construct template to pu Bi plate inside

//  logictemp_pl_plate = new G4LogicalVolume(templatebox,          //its solid
//                        nist->FindOrBuildMaterial(env_mat_name),           //its material
//                        "LayerPlBi");

/*
  for (G4int i = 0; i<N; i++){
      temp_str = pref + std::to_string(i);
      if (i!=Npos){
          phys_vol = new G4PVPlacement(0,                    //no rotation
                    G4ThreeVector(0, 0, -0.5*env_sizeZ + sc_vol_st/2 + sc_vol_st*i),
                    logictempl,                //its logical volume
                    temp_str,                     //its name
                    logicEnv,                     //its mother  volume
                    false,                        //no boolean operation
                    0,                            //copy number
                    checkOverlaps);                        //overlaps checking
          fScoringVolumes.push_back(phys_vol);
      }
      else {
          phys_vol = new G4PVPlacement(0,                    //no rotation
                  G4ThreeVector(0, 0, -0.5*env_sizeZ + sc_vol_st/2 + sc_vol_st*i),
                  logictemp_pl_plate,                //its logical volume
                  temp_str,                     //its name
                  logicEnv,                     //its mother  volume
                  false,                        //no boolean operation
                  0,                            //copy number
                  checkOverlaps);                        //overlaps checking
          fScoringVolumes.push_back(phys_vol);
      }
  };
*/

// place shape2 inside one of volumes

/*
  G4Box* solidShape2 =
    new G4Box("Shape2",                      //its name
              0.5*shape2_dxy,0.5*shape2_dxy, 0.5*shape2_dz
    ); //its size


//shape2_mat = nist->FindOrBuildMaterial(shape2_mat_name);
  logicShape2 =
    new G4LogicalVolume(solidShape2,         //its solid
                        TissueWithAdmixture,  //its material
                        "Shape2");           //its name 

// find mother volume

  G4VPhysicalVolume* bismutvol = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0, pos_log_z),                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logictemp_pl_plate,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  if (bismutvol->CheckOverlaps()) {
    std::cout << "overlaps seen, move Bismuth plate a little bit, abort run" << std::endl;
  }
*/
}

void DetectorConstruction::UpdateGeometry()
{
  cout<<"=================== Reinitialization ====================="<<endl;
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
//  G4LogicalSkinSurface::CleanSurfaceTable();
//  G4LogicalBorderSurface::CleanSurfaceTable();


//  ComputeCaloParameters();
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
//  G4RunManager::GetRunManager()->ReinitializeGeometry();
//  if ( physic_composite )
//   physic_composite->CheckOverlaps();
//  if ( physic_lif )
//   physic_lif->CheckOverlaps();
//  if ( physic_polyprop  )
//   physic_polyprop->CheckOverlaps();
//  if ( physic_alum )
//   physic_alum->CheckOverlaps();
//  if ( physic_water )
//   physic_water->CheckOverlaps();
/*
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
//  G4LogicalSkinSurface::CleanSurfaceTable();
//  G4LogicalBorderSurface::CleanSurfaceTable();
*/

  ComputeCalorParameters();
/*
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
*/
}

G4Material* DetectorConstruction::FindMaterial(G4String name)
{
    G4Material* material = G4Material::GetMaterial(name,true);
    if ( material )
     cout<<"material name density "<<material->GetName()<<" "<<material->GetDensity()/(g/cm3)<<endl;
    return material;
}

