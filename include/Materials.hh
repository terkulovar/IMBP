
#ifndef Materials_h
#define Materials_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class Materials
{
  public:

    ~Materials();
 
    static Materials* GetInstance();

    G4Material* GetMaterial(const G4String);
 
  private:
 
    Materials();

    void CreateMaterials();

  private:

    static Materials* instance;

    G4NistManager*     nistMan;

//    G4double           Clad1ReflIndex;
//    G4double           Clad2ReflIndex;
//    G4double           FiberReflIndex;
//    G4double           GlassReflIndex;


};

#endif /*Materials_h*/
