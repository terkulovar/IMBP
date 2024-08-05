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

//  if ( volume == detector->GetWater() )
//  {
//   my_total += edepStep/CLHEP::MeV;
//   G4double Energy;
//   const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
//    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
//   if (generatorAction)
//   {
//    const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
//    Energy = particleGun->GetParticleEnergy();
//   }
//   man->FillNtupleIColumn(1,0,evId);
//   if ( edepStep > 0 )
//   {
//    G4AnalysisManager *man = G4AnalysisManager::Instance();
//    man->FillNtupleDColumn(1,0,Energy/CLHEP::MeV);
//    man->FillNtupleDColumn(1,1,vector.z()/mm);
//    man->FillNtupleDColumn(1,2,edepStep/CLHEP::MeV);
//    man->FillNtupleDColumn(1,3,stepl/mm);
//    man->AddNtupleRow(1);
//   }
//   cout<<"z z edep edep e e edeptot "<<vector.z()<<" "<<G4BestUnit(vector.z(),"Length")<<" "<<vector.z()/cm<<" "<<volume->GetName()
//   <<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<endl;
//   G4double mass = volume->GetLogicalVolume()->GetDaughter(0)->GetLogicalVolume()->GetMass();
//   cout<<"daugter "<<volume->GetLogicalVolume()->GetNoDaughters()<<" "<<
//    volume->GetLogicalVolume()->GetDaughter(0)->GetName()<<" "<<
//    mass/kg<<endl;
//  }
//  if ( volume->GetName() == "BLiF" || volume->GetName() == "LiF"  )
//   cout<<"z z edep edep name post pre "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<
//   edepStep/CLHEP::MeV<<" "<<volume->GetName()<<" "<<volume11->GetLogicalVolume()->GetName()<<" "<<volume12->GetName()<<endl;
  G4int evId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
//  const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
//   G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
//  if (generatorAction)
//  {
//   const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
//   Energy = particleGun->GetParticleEnergy();
//  }
  const G4Run* aRun = G4RunManager::GetRunManager()->GetCurrentRun();
//   G4int nevents = aRun->GetNumberOfEvent();
  G4int nevents = aRun->GetNumberOfEventToBeProcessed();
//  cout<<"events "<<nevents<<endl;
  G4double mass;
  G4double QF;
  G4double L = (edepStep/CLHEP::MeV)/(stepl/cm);
  if ( L < 100 ) QF = 1;
  else if ( L > 100 && L <1000 ) QF = 0.32*L/10. - 2.2;
  else if (L > 1000 ) QF = 300*pow(L/10.,-0.5);
//  G4double zpos = (step->GetPreStepPoint()->GetPosition().z()+step->GetPostStepPoint()->GetPosition().z())/2.;
//   cout<<"event id energy L QF zpos "<<evId<<" "<<Energy<<" "<<L<<" "<<QF<<" "<<zpos<<endl;
//  G4String name1, name2;
  G4AnalysisManager *man = G4AnalysisManager::Instance();
//   cout<<"de z "<<edepStep/CLHEP::MeV<<" "<<vector.z()/cm<<" "<<volume->GetName()<<" "<<nevents<<endl;
//   man->FillH1(0,vector.z()/cm,(edepStep/CLHEP::MeV)/nevents);
//  const DetectorConstruction* detConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//  G4double XY = detConstruction->GetSampleXY();
//  G4double RW = detConstruction->GetWaterR();
//   cout<<"sample XY XY cm XY mm Z "<<XY<<" "<<XY/cm<<" "<<XY/mm<<" "<<detector->GetSampleZ()/cm<<endl;

  if ( volume == detector->GetEye() )
  {
   mass = volume->GetLogicalVolume()->GetMass();
   G4double edep1 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
//   man->FillH1(0,vector.z()/cm,edep1/nevents);
//   if ( stepl > 1.e-10 )
//    man->FillH1(1,vector.z()/cm,edep1*QF/nevents);
//   if ( name == "proton" )
//    man->FillH1(2,kinetic_energy);
//   if ( name == "neutron" )
//    man->FillH1(13,kinetic_energy);
   man->FillNtupleIColumn(1,0,evId);
   man->FillNtupleDColumn(1,1,(Energy/CLHEP::MeV)/nevents);
   man->FillNtupleSColumn(1,2,volume->GetName());
   man->FillNtupleSColumn(1,3,name);
   man->FillNtupleDColumn(1,4,vector.z()/cm);
   man->FillNtupleDColumn(1,5,(edepStep/CLHEP::MeV)/nevents);
   man->FillNtupleDColumn(1,6,edep1/nevents);
   man->FillNtupleDColumn(1,7,edep1*QF/nevents);
   man->AddNtupleRow(1);
//   fEventAction->AddEdep( edep1, edep1*QF, volume->GetName() );
//   cout<<"z z edep edep name mass QF name mass "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<" "<<nevents<<endl;
  }
  if ( volume == detector->GetTestis() )
  {
   mass = volume->GetLogicalVolume()->GetMass();
   G4double edep1 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
//   man->FillH1(0,vector.z()/cm,edep1/nevents);
//   if ( stepl > 1.e-10 )
//    man->FillH1(1,vector.z()/cm,edep1*QF/nevents);
//   if ( name == "proton" )
//    man->FillH1(3,kinetic_energy);
//   if ( name == "neutron" )
//    man->FillH1(14,kinetic_energy);
   man->FillNtupleIColumn(1,0,evId);
   man->FillNtupleDColumn(1,1,(Energy/CLHEP::MeV)/nevents);
   man->FillNtupleSColumn(1,2,volume->GetName());
   man->FillNtupleSColumn(1,3,name);
   man->FillNtupleDColumn(1,4,vector.z()/cm);
   man->FillNtupleDColumn(1,5,(edepStep/CLHEP::MeV)/nevents);
   man->FillNtupleDColumn(1,6,edep1/nevents);
   man->FillNtupleDColumn(1,7,edep1*QF/nevents);
   man->AddNtupleRow(1);
//   cout<<"z z edep edep name mass QF "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<endl;
  }
  if ( volume == detector->GetSpleen() )
  {
   mass = volume->GetLogicalVolume()->GetMass();
   G4double edep1 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
//   man->FillH1(0,vector.z()/cm,edep1/nevents);
//   if ( stepl > 1.e-10 )
//    man->FillH1(1,vector.z()/cm,edep1*QF/nevents);
//   if ( name == "proton" )
//    man->FillH1(4,kinetic_energy);
//   if ( name == "neutron" )
//    man->FillH1(15,kinetic_energy);
   man->FillNtupleIColumn(1,0,evId);
   man->FillNtupleDColumn(1,1,(Energy/CLHEP::MeV)/nevents);
   man->FillNtupleSColumn(1,2,volume->GetName());
   man->FillNtupleSColumn(1,3,name);
   man->FillNtupleDColumn(1,4,vector.z()/cm);
   man->FillNtupleDColumn(1,5,(edepStep/CLHEP::MeV)/nevents);
   man->FillNtupleDColumn(1,6,edep1/nevents);
   man->FillNtupleDColumn(1,7,edep1*QF/nevents);
   man->AddNtupleRow(1);
//   cout<<"z z edep edep name mass QF "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<endl;
  }
  if ( volume == detector->GetBrain() )
  {
   mass = volume->GetLogicalVolume()->GetMass();
   G4double edep1 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
//   man->FillH1(0,vector.z()/cm,edep1/nevents);
//   if ( stepl > 1.e-10 )
//    man->FillH1(1,vector.z()/cm,edep1*QF/nevents);
//   if ( name == "proton" )
//    man->FillH1(5,kinetic_energy);
//   if ( name == "neutron" )
//    man->FillH1(16,kinetic_energy);
   man->FillNtupleIColumn(1,0,evId);
   man->FillNtupleDColumn(1,1,(Energy/CLHEP::MeV)/nevents);
   man->FillNtupleSColumn(1,2,volume->GetName());
   man->FillNtupleSColumn(1,3,name);
   man->FillNtupleDColumn(1,4,vector.z()/cm);
   man->FillNtupleDColumn(1,5,(edepStep/CLHEP::MeV)/nevents);
   man->FillNtupleDColumn(1,6,edep1/nevents);
   man->FillNtupleDColumn(1,7,edep1*QF/nevents);
   man->AddNtupleRow(1);
//   cout<<"z z edep edep name mass QF "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<endl;
  }
  if ( volume == detector->GetStomach() )
  {
   mass = volume->GetLogicalVolume()->GetMass();
   G4double edep1 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
//   man->FillH1(0,vector.z()/cm,edep1/nevents);
//   if ( stepl > 1.e-10 )
//    man->FillH1(1,vector.z()/cm,edep1*QF/nevents);
//   if ( name == "proton" )
//    man->FillH1(6,kinetic_energy);
//   if ( name == "neutron" )
//    man->FillH1(17,kinetic_energy);
   man->FillNtupleIColumn(1,0,evId);
   man->FillNtupleDColumn(1,1,(Energy/CLHEP::MeV)/nevents);
   man->FillNtupleSColumn(1,2,volume->GetName());
   man->FillNtupleSColumn(1,3,name);
   man->FillNtupleDColumn(1,4,vector.z()/cm);
   man->FillNtupleDColumn(1,5,(edepStep/CLHEP::MeV)/nevents);
   man->FillNtupleDColumn(1,6,edep1/nevents);
   man->FillNtupleDColumn(1,7,edep1*QF/nevents);
   man->AddNtupleRow(1);
//   cout<<"z z edep edep name mass QF "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<endl;
  }
  if ( volume==detector->GetEye() || volume==detector->GetTestis() || volume==detector->GetSpleen() 
       || volume==detector->GetBrain() || volume==detector->GetStomach() || volume==detector->GetWater() )
  {
//   cout<<"z z edep edep name mass QF "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<endl;
//   fEventAction->AddEdep(Energy/CLHEP::MeV, volume->GetName() );
   man->FillNtupleDColumn(3,0,kinetic_energy/nevents);
   man->FillNtupleSColumn(3,1,name);
   man->FillNtupleDColumn(3,2,vector.z()/mm);
   man->FillNtupleDColumn(3,3,(edepStep/CLHEP::MeV)/nevents);
   man->FillNtupleDColumn(3,4,stepl/mm);
//   man->AddNtupleRow(3);
//   if ( name == "proton" )
   if ( stepl > 1.e-10 )
    man->FillH1(1,vector.z()/mm,(edepStep/CLHEP::MeV)/(stepl/mm)/nevents);
//    man->FillH1(1,vector.z()/mm,(edepStep/CLHEP::MeV)/nevents);
//   if ( name == "neutron" )
//   {
//    cout<<"z edep name de/x volume "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//    <<name<<" "<<(edepStep/CLHEP::MeV)/(stepl/mm)<<" "<<volume->GetName()<<endl;
//    man->FillH1(2,vector.z()/mm,(edepStep/CLHEP::MeV)/(stepl/mm));
//   }
  }  
//  if ( volume == detector->GetLiF() )
  if ( strcmp(volume->GetName(),"LiF")==0 || strcmp(volume->GetName(),"LiF1")==0 || strcmp(volume->GetName(),"LiF2")==0
     || strcmp(volume->GetName(),"LiF3")==0 || strcmp(volume->GetName(),"LiF4")==0 || strcmp(volume->GetName(),"LiF5")==0
     || strcmp(volume->GetName(),"LiF6")==0 || strcmp(volume->GetName(),"LiF7")==0 || strcmp(volume->GetName(),"LiF8")==0
     || strcmp(volume->GetName(),"LiF9")==0 )
  {
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
   man->AddNtupleRow(1);
   fEventAction->AddEdep( edep1, edep1*QF, volume->GetName() );
//   if ( edep1 > 1 )
   cout<<"z z edep edep name mass QF dose "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<" "<<edep1<<endl;
//   cerr<<"z z edep edep name mass QF dose "<<vector.z()<<" "<<vector.z()/cm<<" "<<edepStep<<" "<<edepStep/CLHEP::MeV<<" "
//   <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<" "<<edep1<<endl;
  }

/* ----------------- start ------------------- *\
  if ( volume == detector->GetLiF() ) //dosimeter
  {
   if ( kinetic_energy/CLHEP::MeV > 0.00 )
   { 
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    mass = volume->GetLogicalVolume()->GetMass();
    G4double edep11 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
//     mass = mass/kg;
    G4double edep1 = ((edepStep/CLHEP::eV)*e_SI)/(mass/kg);
    man->FillNtupleIColumn(1,0,evId);
    man->FillNtupleDColumn(1,1,Energy/CLHEP::MeV);
    man->FillNtupleSColumn(1,2,volume->GetName());
    man->FillNtupleDColumn(1,3,vector.z()/cm );
//     man->FillNtupleDColumn(1,4,edepStep/CLHEP::MeV);
    man->FillNtupleDColumn(1,4,edep1/nevents);
    man->FillNtupleDColumn(1,5,edep1*QF/nevents);
    man->AddNtupleRow(1);
//    cout<<"z z edep edep name mass QF name1 nev step en "<<vector.z()<<" "<<vector.z()/cm<<" "
//    <<edepStep<<" "<<edepStep/CLHEP::keV<<" "
//    <<name<<" "<<mass/kg<<" "<<QF<<" "<<volume->GetName()<<" "<<nevents<<" "<<stepl/cm<<" "<<kinetic_energy/CLHEP::keV<<endl;
   }
  }
\* -------------- end -------------------- */


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

  
/*
  G4double xi, yi, zi, xi_post, yi_post, zi_post, En, stepdist, xpr, ypr, zpr , zdist_fluence_pre, zdist_fluence_post;
  G4int i_z_dEdx,i_p_dEdx;
  G4int i_z_fluence,i_p_fluence;


  xi = step->GetPreStepPoint()->GetPosition().x();
  yi = step->GetPreStepPoint()->GetPosition().y();
  zi=  step->GetPreStepPoint()->GetPosition().z();

  xi_post = step->GetPostStepPoint()->GetPosition().x();
  yi_post = step->GetPostStepPoint()->GetPosition().y();
  zi_post = step->GetPostStepPoint()->GetPosition().z();

  xpr = fEventAction->xprime;
  ypr = fEventAction->yprime;
  zpr = fEventAction->zprime;


  zdist_fluence_pre = zi - zpr ;
  zdist_fluence_post = zi_post - zpr ;

  stepdist = sqrt( (xi-xi_post)*(xi-xi_post) + (yi-yi_post)*(yi-yi_post) + (zi-zi_post)*(zi-zi_post) );


  i_z_fluence = int(zdist_fluence_pre/fRunAction->stepforfluence);
  i_p_fluence = int(zdist_fluence_post/fRunAction->stepforfluence);
  En = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;

  if (mytrack->GetTrackID() == 1){
    if (fEventAction->distdEdx<fRunAction->MaxZ) {
        i_z_dEdx = int(fEventAction->distdEdx/fRunAction->stepfordEdz);
        fEventAction->vdEdz.at(i_z_dEdx) = fEventAction->vdEdz.at(i_z_dEdx) + edepStep;
        fEventAction->distdEdx+= stepdist;
        i_p_dEdx = int(fEventAction->distdEdx/fRunAction->stepfordEdz);
        if (i_p_dEdx != i_z_dEdx ) fEventAction->vEn.at(i_z_dEdx) = En;
    };
  };

  if (i_z_fluence != i_p_fluence){
    int ismaller,ibigger;
    if (G4RunManager::GetRunManager()->GetCurrentEvent()){
      if (i_z_fluence<i_p_fluence){
        ismaller = i_z_fluence;
        ibigger = i_p_fluence;
      } else
      {
        ismaller = i_p_fluence;
        ibigger = i_z_fluence;
      };
      for (int i=ismaller; i<ibigger;i++){
        if(En>0.00){ //keeping only fluences for particles with  Energy >10 keV
        man->FillNtupleIColumn(0,0,G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
        man->FillNtupleDColumn(0,1,xi);
        man->FillNtupleDColumn(0,2,yi);
        man->FillNtupleDColumn(0,3,i*fRunAction->stepforfluence);
        man->FillNtupleDColumn(0,4,En);
        man->FillNtupleIColumn(0,5,mytrack->GetDefinition()->GetPDGEncoding());
        man->FillNtupleSColumn(0,6,mytrack->GetDefinition()->GetParticleName());
        man->FillNtupleSColumn(0,7,mytrack->GetMaterial()->GetName());
        man->AddNtupleRow(0);
       }
      }
    }
  };
*/
/*
  G4VPhysicalVolume* volume  = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (fScoringVolumes.size() == 0) {
      const DetectorConstruction* detConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      fScoringVolumes = detConstruction->GetScoringVolumes();
  };

  for (uint i = 0; i<fScoringVolumes.size(); i++){
    if (volume == fScoringVolumes.at(i)){
        fEventAction->AddEdep(edepStep,i);
    };
    if (volume->GetMotherLogical()){
      if ((volume->GetMotherLogical()->GetName() != "World") && (volume->GetMotherLogical()->GetName() == fScoringVolumes.at(i)->GetLogicalVolume()->GetName())){
          fEventAction->AddEdep(edepStep,i);
      };
    };
  };
*/
}
