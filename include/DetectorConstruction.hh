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


#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4UserLimits.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
//#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
//#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4PVPlacement.hh"
#include "G4HumanPhantomMaterial.hh"
//#include <Riostream.h>
using namespace std;

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Materials;
class DetectorMessenger;

/// Detector construction class to define materials and geometry.


class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void UpdateGeometry();
    void SetWallThickness(G4double);
//    const G4String GetLiFName() {return physic_lif->GetName();}; 
    G4VPhysicalVolume* GetLiF() {return physic_lif;}
    G4LogicalVolume* GetLiFLog() {return logic_lif;}
    G4VPhysicalVolume* GetWater() {return physic_water;}
    G4VPhysicalVolume* GetEye() const {return physic_eye;}
    G4VPhysicalVolume* GetTestis() const {return physic_testis;}
    G4VPhysicalVolume* GetSpleen() const {return physic_spleen;}
    G4VPhysicalVolume* GetBrain() const {return physic_brain;}
    G4VPhysicalVolume* GetStomach() const {return physic_stomach;}
    G4double GetSampleXY() const {return sample_sizeXY;}
    G4double GetSampleZ() const {return sample_sizeZ;}

//    std::vector<G4VPhysicalVolume*> GetScoringVolumes() const { return fScoringVolumes; }

// Detector parameters


//    const G4double shape2_dxy = 10*cm, shape2_dz = 1*cm;
//    const G4double env_sizeXY = 10*cm, env_sizeZ = 10*cm;
//    const G4double sc_vol_st = 1.0*cm;
//    const G4ThreeVector pos2 = G4ThreeVector(0, 0, 0.0*cm);
//    G4double world_sizeXY = 1.2*env_sizeXY;
//    G4double world_sizeZ  = 1.2*env_sizeZ;

//========platic A-150-tissue  with incorporated Bi particles =======================
    const G4String env_mat_name = "G4_WATER";
    const G4String world_mat_name = "G4_AIR";
    const G4String lif_mat_name = "LiF";
//    const G4String Admixture_name = "G4_Bi";
    //const G4String Tissue_name = "G4_A-150_TISSUE";
//    const G4String Tissue_name = "soft_tissue"; //soft tissues from G4HumanPhantomMaterial

//    const G4String TissueWithAdmixture_name = "TissueWithAdmixture";

    //G4double  rho_tissue= 1.00*g/cm3; //this is just water
//    G4double  rho_tissue= 0.9869*g/cm3; //soft tissue
//    G4double  rho_admix = 9.79*g/cm3;
//    G4double admix_c = 0*mg/L; //concetration in  mg per Litre ++ as per Kolobov++++
    //G4double admix_c = 3000000*mg/L; //half of the plate from Bi  ++++test++++
    //G4double admix_c = 9790000*mg/L; //100% of the plate from Bi  ++++test++++

//    G4Material* TissueWithAdmixture;
//    G4Material* Tissue;
//    G4Material* Admixture;
  //  G4Material* soft;
//==================================================================================


//

    //G4Material* shape2_mat;

//    G4Material* env_mat;
//    G4Material* world_mat;
//    G4LogicalVolume* logicEnv,*logicWorld,*logicShape2,*logictempl,*logictemp_pl_plate;
//    G4VPhysicalVolume* physEnv;
//    G4bool checkOverlaps = true;

    void DefineMaterials();
//    G4double GetZedge() {return (polyprop_sizeZ+composite_r_max+alum_sizeZ+1.*gap_sizeZ+0.5*cm);};
    G4double GetZedge() const {return (polyprop_sizeZ+composite_sizeZ+alum_sizeZ+1.*gap_sizeZ);};
//    G4double GetZedge() const {return (-alum_sizeZ/2.-13.0*cm);};
//    G4double GetZedge() const {return (composite_sizeZ+gap_sizeZ1+polyprop_sizeZ+alum_sizeZ);};
    G4double GetWorldR() const {return world_R;};
    G4double GetR() const {return ((water_R+gap_sizeZ+composite_sizeZ+gap_sizeZ+polyprop_sizeZ+alum_sizeZ));};
    G4double GetWaterR() const {return water_R;};
    void SetComposDet(G4String tt) { fCompDet = tt;}

  private:
    G4LogicalVolume* logic_alum, *logic_composite, *logic_water, *logic_polyprop, *logic_lif, *logic_polyprop1;
    G4LogicalVolume* logic_illuminator, *logic_eye, *logic_testis, *logic_spleen, *logic_stomach, *logic_brain;
    G4VPhysicalVolume* physic_alum, *physic_composite, *physic_water, *physic_polyprop, *physic_lif, *physic_polyprop1;
    G4VPhysicalVolume* physic_illuminator, *physic_eye, *physic_testis, *physic_spleen, *physic_stomach, *physic_brain;
    G4Material* composite_mat, *alum_mat, *water_mat, *polyprop_mat, *lif_mat, *world_mat, *polyethylene_mat, *ferrum_mat; 
    G4Material* illuminator_mat, *eye_mat, *testis_mat, *spleen_mat, *stomach_mat, *brain_mat;
    G4double world_sizeXY;
    G4double world_sizeZ;
    G4double world_R;
    G4double sizeXY;  
    G4double water_sizeZ;
    G4double water_R;
    G4double sample_sizeXY;
    G4double sample_sizeZ;
    G4double alum_sizeZ;
    G4double aluminium_lengthZ;
    G4double polyprop_sizeZ;
    G4double polyprop_sizeZ1;
    G4double composite_sizeZ;
    G4double composite_sizeXY;
    G4double composite_Rmax; 
    G4double composite_Rmin; 
    G4double composite_thickness;
    G4double gap_sizeZ;
    G4double gap_sizeZ1;
    G4double gap_sizeZ2;
    G4double lif_sizeZ;
    G4double lif_sizeX;
    G4double lif_sizeY;
    G4double lif_R;
    G4double lif_H;
    G4double water_sizeXY;
    G4double aluminium_R;
    G4UserLimits* fStepLimit = nullptr;
    Materials*         materials;
    void ConstructScoringVolumes();
    G4Material* FindMaterial(G4String);
    void ComputeCalorParameters();
    DetectorMessenger* detectorMessenger;  //pointer to the Messenger
    G4String fCompDet;
    G4double aluminium_positionZ;
  protected:
//    std::vector<G4VPhysicalVolume*> fScoringVolumes {};
//    G4int N_vol;
};

inline void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the world
     world_sizeXY = 2.1*sizeXY;
//     world_sizeZ = (polyprop_sizeZ+composite_r_max+alum_sizeZ+1.*gap_sizeZ)*2.5;
     world_sizeZ = (water_sizeZ)*2.5;
     world_R = (water_R+gap_sizeZ+composite_sizeZ+gap_sizeZ+polyprop_sizeZ+alum_sizeZ)*1.1;
     world_R = 1.02*aluminium_R;
     cout<<"world xy z R "<<world_sizeXY/cm<<" "<<world_sizeZ/cm<<" "<<world_R/cm<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
