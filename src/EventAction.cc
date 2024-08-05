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
/// \file EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"
#include "MyRunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include <vector>

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(MyRunAction* runAction, PrimaryGeneratorAction* generatorAction)
: fRunAction(runAction),
  fGeneratorAction(generatorAction)
{
 Dose_MeV = 0;
 Dose_Gr = 0;
 Dose_Zi = 0;
 f1 = -1;
 f2 = -1;
 f3 = -1;
 f4 = -1;
 f5 = -1;
 iLiF = 0;
 dose.clear();
 edose.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdep(G4double energy, G4double energy1, G4String name)
{
 if ( f1 == 0 )
 {
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  analysisManager->FillH1(1,energy);
  f1 = 1;
  iLiF++;
 }
 myname = name;
 dose.push_back(energy);
 edose.push_back(energy1);
/*
 if ( name == "Eye" && f1 == 0 )
 {
//  analysisManager->FillNtupleDColumn(2,0,energy);
//  analysisManager->FillNtupleSColumn(2,1,name);
//  analysisManager->AddNtupleRow(2);
//  analysisManager->FillH1(7,energy);
//  cout<<"----------------------------- eye "<<energy<<endl;
  f1 = 1;
 }
 if ( name == "Testis" && f2 == 0 )
 {
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
//  analysisManager->FillH1(8,energy);
  f2 = 1;
 }
 if ( name == "Spleen" && f3 == 0 )
 {
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
//  analysisManager->FillH1(9,energy);
  f3 = 1;
 }
 if ( name == "Brain" && f4 == 0 )
 {
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
//  analysisManager->FillH1(10,energy);
  f4 = 1;
 }
 if ( name == "Stomach" && f5 == 0 )
 {
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
//  analysisManager->FillH1(16,energy);
  f5 = 1;
 }
*/
/*
 Dose_MeV += edep;
 G4double dd = ((edep/CLHEP::eV)*e_SI)/(mass/kg);
 Dose_Gr += dd;
// G4cout<<"Dose_MeV Dose_Gr "<<Dose_MeV<<" "<<Dose_Gr<<" "<<mass/kg<<" "<<((edep/CLHEP::eV)*e_SI)<<G4endl;
// G4cout<<"l "<<step/um<<G4endl;
 G4double L = (edep/CLHEP::keV)/(step/um);
 G4double a0 = 2.34e-4;
 G4double a1 = 2.14e-4;
 G4double a2 = 4.75e-5;
 G4double QF;
 if ( L > 0 )
 {
  QF = (1.-std::exp(-a1*L-a2*L*L))/(a0*L);
  Dose_Zi += QF*dd;
 }
// cout<<"!!!!!!!!!!!!!!! QF "<<QF<<" "<<Dose_Gr<<" "<<Dose_Zi<<" "<<dd<<endl;
// l=L(МэВ/см)/10 [кэв/мкм]
// QF = (1 - exp(-a01*l-a2*l^2))/(a0*l)  
// где a0 = 2.34×10-4 мкм/кэВ; a1 = 2.14×10-4 мкм/кэв; a2 = 4.75×10-5 (мкм/кэВ)^2
*/
}


void EventAction::BeginOfEventAction(const G4Event* evt)
{
  evtNb = evt->GetEventID();
  f1 = 0;
  f2 = 0;
  f3 = 0;
  f4 = 0;
  f5 = 0;
//  G4double My_phi = G4UniformRand()*CLHEP::twopi;
//  G4double My_theta = -1.0 + 2.0*G4UniformRand();
//  G4ThreeVector my_pos(cos(My_phi)*sin(My_theta),sin(My_phi)*sin(My_theta),cos(My_theta));
//  G4double R1 = Detector->GetWorldR();
//  cout<<"!!!!!!!!!! R  "<<R1/m<<endl;
//  my_pos.setMag(1.01*m);
  G4ThreeVector my_pos = fGeneratorAction->GetParticleGun()->GetParticlePosition();
  cout<<"my_pos x y z mag "<<my_pos.x()/cm<<" "<<my_pos.y()/cm<<" "<<my_pos.z()/cm<<" "<<my_pos.mag()/m<<endl;
  G4ThreeVector my_dir = fGeneratorAction->GetParticleGun()->GetParticleMomentumDirection();
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  analysisManager->FillH1(3,cos(my_pos.theta()));
  analysisManager->FillH1(4,cos(my_pos.angle(my_dir)));

//  distdEdx = 0.;

// get primary vertex
//  xprime = fGeneratorAction->GetParticleGun()->GetParticlePosition().x();
//  yprime = fGeneratorAction->GetParticleGun()->GetParticlePosition().y();
//  zprime = fGeneratorAction->GetParticleGun()->GetParticlePosition().z();

//  vdEdz = InitializeZVector(fRunAction->MinZ, fRunAction->MaxZ, fRunAction->stepfordEdz);
//  vEn = InitializeEnVector(fRunAction->MinZ, fRunAction->MaxZ, fRunAction->stepforfluence);
//  const DetectorConstruction* detConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//  fEdepV = Initialize_EinVol_Vector(detConstruction->GetScoringVolumes().size());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4double sum = 0;
  G4double sum1 = 0;
  const G4Run* aRun = G4RunManager::GetRunManager()->GetCurrentRun();
  G4int nevents = aRun->GetNumberOfEventToBeProcessed();
  if ( dose.size() > 0 )
  {
   for (int i = 0; i < dose.size(); i++)
   {
    G4double bb = dose.at(i);
    sum += bb;
    bb = edose.at(i);
    sum1 += bb;
   }
  }
  G4int evtNb = evt->GetEventID();
  G4double En = fGeneratorAction->GetParticleGun()->GetParticleEnergy();
  cout<<"!!!!!!!!!!!!!! event total "<<evtNb<<" "<<nevents<<endl;
  if ( evtNb%50000 == 0 || sum > 1.e-11 )
   cerr<<"events iLif dose edose en "<<evtNb<<" "<<iLiF<<" "<<sum<<" "<<sum1<<" "<<En<<" "<<myname<<endl;
  dose.clear();
  edose.clear();

//  const DetectorConstruction* detConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//  G4int evtNb = evt->GetEventID();
//  cout<<"!!!!!!!!!!!!!! event energy Dose_Gr Dose_Zi "<<evtNb<<" "<<En<<" "<<Dose_Gr<<" "<<Dose_Zi<<endl;
  G4ThreeVector kproton = fGeneratorAction->GetParticleGun()->GetParticleMomentumDirection();
  G4double X = fGeneratorAction->GetParticleGun()->GetParticlePosition().x()/cm;
  G4double Y = fGeneratorAction->GetParticleGun()->GetParticlePosition().y()/cm;
  G4double Z = fGeneratorAction->GetParticleGun()->GetParticlePosition().z()/cm;
  G4double phi = fGeneratorAction->GetParticleGun()->GetParticlePosition().phi();
  G4double theta = fGeneratorAction->GetParticleGun()->GetParticlePosition().theta();
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  analysisManager->FillNtupleDColumn(2,0,En);
  analysisManager->FillNtupleDColumn(2,1,kproton.x());
  analysisManager->FillNtupleDColumn(2,2,kproton.y());
  analysisManager->FillNtupleDColumn(2,3,kproton.z());
  analysisManager->FillNtupleDColumn(2,4,kproton.phi());
  analysisManager->FillNtupleDColumn(2,5,kproton.theta());
  analysisManager->FillNtupleDColumn(2,6,X);
  analysisManager->FillNtupleDColumn(2,7,Y);
  analysisManager->FillNtupleDColumn(2,8,Z);
  analysisManager->FillNtupleDColumn(2,9,phi);
  analysisManager->FillNtupleDColumn(2,10,theta);
//  analysisManager->AddNtupleRow(2);
  analysisManager->FillH1(0,En);

  if ( sum > 1.e-11 )
  {
   analysisManager->FillNtupleDColumn(4,0,sum/nevents);
   analysisManager->FillNtupleDColumn(4,1,sum1/nevents);
   analysisManager->AddNtupleRow(4);
  }

//  cout<<"phi theta "<<kproton.phi()<<" "<<kproton.theta()<<endl;
//  G4RunManager* runManager = G4RunManager::GetRunManager();
//  auto* pointerRun = (MyRunAction*)(runManager->GetUserRunAction());
//  pointerRun->FillPerEvent(Dose_MeV, Dose_Gr, Dose_Zi);
//  Dose_MeV = 0;
//  Dose_Gr = 0;
//  Dose_Zi = 0;
//  G4double dose;
//  G4double mass = 0.;

/*
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  for (uint i = 0; i<vdEdz.size();i++){
    if (vEn.at(i)>-0.1){
      man->FillNtupleDColumn(1,0,vdEdz.at(i));
      man->FillNtupleDColumn(1,1,fRunAction->stepfordEdz);
      man->FillNtupleDColumn(1,2,i*fRunAction->stepfordEdz);
      man->FillNtupleDColumn(1,3,vEn.at(i));
      man->FillNtupleIColumn(1,4,G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
      man->AddNtupleRow(1);
    }
  }
  // accumulate statistics in run action
  for (uint i=0;i<fEdepV.size();i++){
  // calculate mass to conver to grays
     if (fEdepV.at(i) > 0.) {
       if (detConstruction->GetScoringVolumes().at(i)->GetLogicalVolume()->GetName()== "Layer"){
         mass = detConstruction->GetScoringVolumes().at(i)->GetLogicalVolume()->GetMass();
       };
       if (detConstruction->GetScoringVolumes().at(i)->GetLogicalVolume()->GetName()== "LayerPlBi"){
         mass = detConstruction->GetScoringVolumes().at(i)->GetLogicalVolume()->GetMass() + detConstruction->logicShape2->GetMass();
       };
       //dose = (fEdepV.at(i)/CLHEP::eV)/(mass*e_SI);
        dose = ((fEdepV.at(i)/CLHEP::eV)*e_SI)/(mass/kg);
       man->FillNtupleDColumn(2,0,fEdepV.at(i));
       man->FillNtupleIColumn(2,1,i);
       man->FillNtupleIColumn(2,2,G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
       man->FillNtupleDColumn(2,3,dose);
//       man->FillNtupleDColumn(2,4,(detConstruction->sc_vol_st)*i + detConstruction->sc_vol_st/2);
       man->AddNtupleRow(2);
    }
  };
*/

}

/*
std::vector<G4double> EventAction::InitializeZVector(G4double Min_Z,G4double Max_Z, G4double step){
    G4double diff = Max_Z - Min_Z;
    std::vector<G4double> temp = {};
    for (int j = 0; j<int(diff/step); j++){
      temp.push_back(0);
    };
    return temp;
}

std::vector<G4double> EventAction::InitializeEnVector(G4double Min_Z,G4double Max_Z, G4double step){
    G4double diff = Max_Z - Min_Z;
    std::vector<G4double> temp = {};
    for (int j = 0; j<int(diff/step); j++){
      temp.push_back(-1.);
    };
    return temp;
}

std::vector<G4double> EventAction::Initialize_EinVol_Vector(G4int N){
  std::vector<G4double> temp = {};
  for (int j = 0; j<N; j++){
    temp.push_back(0.);
  };
  return temp;
}
*/
