#ifndef MyRunAction_HH
#define MyRunAction_HH 1

#include "G4UserRunAction.hh"
//#include "G4Run.hh"
#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4Run;
class MyRunAction : public G4UserRunAction
{
public:
    MyRunAction();
    ~MyRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    void FillPerEvent(G4double, G4double, G4double);

    const G4double MinZ = 0.*CLHEP::cm;
    const G4double MaxZ = 100*CLHEP::cm;

    const G4double stepfordEdz = 1.* CLHEP::mm;
    const G4double stepforfluence = 1.* CLHEP::mm;

private:
    G4double Dose_MeV;
    G4double Dose_Gr;
    G4double Dose_Zi;

};

#endif
