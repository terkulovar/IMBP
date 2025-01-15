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
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "MyRunAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Run.hh"
#include "G4TouchableHistory.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MySteppingAction::MySteppingAction(DetectorConstruction* det, EventAction* eventAction)
:detector(det), fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MySteppingAction::~MySteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void MySteppingAction::UserSteppingAction(const G4Step* step)
{
  // collect energy deposited in this step

  G4VPhysicalVolume* volume11
  = step->GetPostStepPoint()->GetTouchable()->GetVolume();
  G4VPhysicalVolume* volume12
  = step->GetPreStepPoint()->GetTouchable()->GetVolume();
//  cout<<"history depth "<<step->GetPostStepPoint()->GetTouchable()->GetHistoryDepth()<<" "
//  <<step->GetPreStepPoint()->GetTouchable()->GetHistoryDepth()<<endl;

  G4double edepStep = step->GetTotalEnergyDeposit();
  G4double stepl = step->GetStepLength();
  G4Track* theTrack = step->GetTrack ();
  G4VPhysicalVolume* volume = theTrack->GetVolume();
//  G4VPhysicalVolume* volume1 = volume->GetLogicalVolume()->GetDaughter(1);
  G4ThreeVector vector = theTrack->GetPosition();
//  G4ThreeVector trans =
//    aStep->GetPreStepPoint()->GetTouchableHandle()->GetTranslation();
  G4ThreeVector moment = theTrack->GetMomentum();
  G4double kinetic_energy = theTrack->GetKineticEnergy();
  G4Material* material = volume->GetLogicalVolume()->GetMaterial();
  G4String name = theTrack->GetDefinition()->GetParticleName();
  const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4double Energy;
  G4ThreeVector my_dir, my_pos;
  if (generatorAction)
  {
   const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
   my_dir = generatorAction->GetParticleGun()->GetParticleMomentumDirection();
   my_pos = generatorAction->GetParticleGun()->GetParticlePosition();
   Energy = particleGun->GetParticleEnergy();
  }

  G4int evId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  G4double xi, yi, zi, xi_post, yi_post, zi_post;
  xi = step->GetPreStepPoint()->GetPosition().x();
  yi = step->GetPreStepPoint()->GetPosition().y();
  zi=  step->GetPreStepPoint()->GetPosition().z();
  xi_post = step->GetPostStepPoint()->GetPosition().x();
  yi_post = step->GetPostStepPoint()->GetPosition().y();
  zi_post = step->GetPostStepPoint()->GetPosition().z();
  G4double mystep = sqrt((xi_post-xi)*(xi_post-xi)+(yi_post-yi)*(yi_post-yi)+(zi_post-zi)*(zi_post-zi));
  const G4Run* aRun = G4RunManager::GetRunManager()->GetCurrentRun();
  G4int nevents = aRun->GetNumberOfEventToBeProcessed();
  G4double mass;
  G4double QF;
//  G4double L = (edepStep/CLHEP::MeV)/(stepl/cm);
  G4double L = (edepStep/CLHEP::MeV)/(mystep/cm);
  if ( L < 100 ) QF = 1;
  else if ( L > 100 && L <1000 ) QF = 0.32*L/10. - 2.2;
  else if (L > 1000 ) QF = 300*pow(L/10.,-0.5);
//  G4double zpos = (step->GetPreStepPoint()->GetPosition().z()+step->GetPostStepPoint()->GetPosition().z())/2.;
//   cout<<"event id energy L QF zpos "<<evId<<" "<<Energy<<" "<<L<<" "<<QF<<" "<<zpos<<endl;
//  G4String name1, name2;
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  if ( strcmp(volume->GetName(),"LiF")==0 || strcmp(volume->GetName(),"LiF1")==0 || strcmp(volume->GetName(),"LiF2")==0
     || strcmp(volume->GetName(),"LiF3")==0 || strcmp(volume->GetName(),"LiF4")==0 || strcmp(volume->GetName(),"LiF5")==0
     || strcmp(volume->GetName(),"LiF6")==0 || strcmp(volume->GetName(),"LiF7")==0 || strcmp(volume->GetName(),"LiF8")==0
     || strcmp(volume->GetName(),"LiF9")==0 )
  {
//   if ( name != "gamma" && name != "neutron" )
   {
//    cout<<"!!!!!!!!particle "<<name<<endl;
    const std::vector<G4Track*>* STracks = step->GetSecondary();
//    cout<<"!!!!!! sec size "<<STracks->size()<<endl;
    if ( STracks->size() > 0 )
    {
     for ( int i = 0; i < STracks->size(); i++ )
     {
      G4Track* t1 = STracks->at(i);
      G4String name2 = t1->GetDefinition()->GetParticleName();
      if ( name2 == "e-" )
      {
       if ( t1->GetTrackStatus() != fStopAndKill )
       {
        edepStep += t1->GetKineticEnergy();
        t1->SetTrackStatus(fStopAndKill);
//        cout<<"sec i energy name  "<<i<<" "<<t1->GetKineticEnergy()<<" "<<name2<<endl;
       }
      }
     }
    }
    L = (edepStep/CLHEP::MeV)/(stepl/cm);
    if ( L < 100 ) QF = 1;
    else if ( L > 100 && L <1000 ) QF = 0.32*L/10. - 2.2;
    else if (L > 1000 ) QF = 300*pow(L/10.,-0.5);
   } 
   mass = volume->GetLogicalVolume()->GetMass();
   G4double edep1 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
   man->FillNtupleIColumn(1,0,evId);
   man->FillNtupleDColumn(1,1,(Energy/CLHEP::MeV));
   man->FillNtupleSColumn(1,2,volume->GetName());
   man->FillNtupleSColumn(1,3,name);
   man->FillNtupleDColumn(1,4,vector.z()/cm);
   man->FillNtupleDColumn(1,5,(edepStep/CLHEP::MeV)/nevents);
   man->FillNtupleDColumn(1,6,edep1/nevents);
   man->FillNtupleDColumn(1,7,edep1*QF/nevents);
   man->FillNtupleDColumn(1,8,cos(my_pos.angle(my_dir)));
  //   man->FillNtupleDColumn(1,9,stepl/cm);
   man->FillNtupleDColumn(1,9,mystep/cm);
   man->FillNtupleDColumn(1,10,L);
   man->FillNtupleDColumn(1,11,QF);
   man->FillNtupleIColumn(1,12,particleGun->GetCurrentSourceIndex());
   man->AddNtupleRow(1);
   fEventAction->AddEdep( edep1, edep1*QF, volume->GetName() );
//   if ( edep1 > 1 )
   cout<<"z z edep edep name mass QF dose "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<" "<<edep1<<endl;
//   cerr<<"z z edep edep name mass QF dose "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<" "<<edep1<<endl;
  }



//   man->FillNtupleIColumn(1,0,evId);
//   man->FillNtupleDColumn(1,1,Energy);
//   man->FillNtupleDColumn(1,2,zpos);
//   man->FillNtupleDColumn(1,3,edep);
//   man->FillNtupleDColumn(1,4,edep1);
//   man->FillNtupleDColumn(1,5,stepl);
//   man->FillNtupleDColumn(1,6,QF);
//   man->AddNtupleRow(1);
//   G4double mass = detector->GetLiF()->GetLogicalVolume()->GetMass();
//   cout<<"mass l l "<<G4BestUnit(mass,"Mass")<<" "<<G4BestUnit(stepl,"Length")<<" "<<stepl<<endl;
//   fEventAction->AddEdep(edep, step, mass);   

}
