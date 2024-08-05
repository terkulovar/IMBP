
#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(
                                           DetectorConstruction* Det)
:Detector(Det)
{ 
//  N03Dir = new G4UIdirectory("/N03/");
//  N03Dir->SetGuidance("UI commands of this example");
  
  detDir = new G4UIdirectory("/det/");
  detDir->SetGuidance("detector control");

  ComposDetCmd = new G4UIcmdWithAString("/det/setCompositeDet",this);
  ComposDetCmd->SetGuidance("put or not composite detector");
  ComposDetCmd->SetParameterName("choice",false);
  ComposDetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZCmd = new G4UIcmdWithADoubleAndUnit("/det/setWallZ",this);
  ZCmd->SetGuidance("Set thickness of wall");
  ZCmd->SetParameterName("thickness",false);
  ZCmd->SetUnitCategory("Length");
  ZCmd->SetRange("thickness>0");
  ZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete ComposDetCmd;   delete UpdateCmd; delete ZCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == ComposDetCmd )
   { Detector->SetComposDet(newValue);}

  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }

  if( command == ZCmd )
    Detector->SetWallThickness(ZCmd->GetNewDoubleValue(newValue));

//  if( command == MagFieldCmd )
//   { Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
