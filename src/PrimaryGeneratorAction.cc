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


#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"


using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
:Detector(DC)
{
//  G4int n_particle = 1;
//  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
//  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//  G4String particleName;
//  G4ParticleDefinition* particle
//    = particleTable->FindParticle(particleName="proton");
//  fParticleGun->SetParticleDefinition(particle);
//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  fParticleGun  = new G4GeneralParticleSource();
//  fParticleGun->SetParticleEnergy(100.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
/*
  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  G4double size = 0.3;
  G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double z0 = -0.5 * envSizeZ;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->GeneratePrimaryVertex(anEvent);
*/
//  fParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
//  G4ThreeVector xproton = fParticleGun->GetCurrentSource()->GetParticlePosition();
//  cout<<"before x y z "<<xproton.x()<<" "<<xproton.y()<<" "<<xproton.z()<<endl;;
  cout<<"----------------------------"<<endl;
//  G4double My_phi = G4UniformRand()*CLHEP::twopi;
//  G4double My_theta = -1.0 + 2.0*G4UniformRand();
//  G4ThreeVector my_pos(cos(My_phi)*sin(My_theta),sin(My_phi)*sin(My_theta),cos(My_theta));
//  G4double R1 = Detector->GetWorldR();
//  cout<<"!!!!!!!!!! R  "<<R1/m<<endl;
//  my_pos.setMag(R1);
//  fParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
//  cout<<"my_pos x y z mag "<<my_pos.x()/cm<<" "<<my_pos.y()/cm<<" "<<my_pos.z()/cm<<" "<<my_pos.mag()/m<<endl;
//  cout<<"pos get "<<fParticleGun->GetParticlePosition().x()<<" "<<fParticleGun->GetParticlePosition().y()
//  <<" "<<fParticleGun->GetParticlePosition().z()<<endl;
//  fParticleGun->GetCurrentSource()->GetPosDist()->SetPosDisType("Surface");
//  fParticleGun->GetCurrentSource()->GetPosDist()->SetPosDisShape("Sphere");
//  fParticleGun->GetCurrentSource()->GetPosDist()->SetRadius(Detector->GetWorldR());
//  cout<<"Radius "<<fParticleGun->GetCurrentSource()->GetPosDist()->GetRadius()/m<<endl;
//  fParticleGun->GetCurrentSource()->GetPosDist()->SetPosRot1(my_pos);
//  fParticleGun->GetCurrentSource()->GetPosDist()->SetPosRot2(my_pos);
//  fParticleGun->SetParticlePosition(my_pos);
//  G4double Alpha = fParticleGun->GetCurrentSource()->GetPosDist()->GetParAlpha();
//  G4double Theta = fParticleGun->GetCurrentSource()->GetPosDist()->GetParTheta();
//  G4double Phi = fParticleGun->GetCurrentSource()->GetPosDist()->GetParPhi();
//  cout<<"alpha theta phi "<<Alpha<<" "<<Theta<<" "<<Phi<<endl;
//  G4ThreeVector my_dir(0,0,-1);
//  my_dir.unit();
//  G4ThreeVector xproton = fParticleGun->GetParticlePosition();
//  cout<<"xproton mag x y z "<<xproton.mag()<<" "<<xproton.x()<<" "<<xproton.y()<<" "<<xproton.z()<<endl;
//  G4ThreeVector xproton1(xproton.x()/xproton.mag(),xproton.y()/xproton.mag(),xproton.z()/xproton.mag());
//  cout<<"xproton1 mag x y z "<<xproton1.mag()<<" "<<xproton1.x()<<" "<<xproton1.y()<<" "<<xproton1.z()<<endl;
//  xproton = fParticleGun->GetParticlePosition();
//  cout<<"xproton mag x y z "<<xproton.mag()<<" "<<xproton.x()<<" "<<xproton.y()<<" "<<xproton.z()<<endl;
//  G4ThreeVector my_dir1(0,0,1);
//  cout<<"angle 1 "<<my_dir.angle(my_dir1)<<endl;
//  G4ThreeVector cros = xproton1.cross(my_dir1);
//  cout<<"cros x y z "<<cros.x()<<" "<<cros.y()<<" "<<cros.z()<<endl;
//  double tt = xproton1.angle(my_dir1);
//  my_dir.rotate(-xproton1.theta(), cros);
//  my_dir.rotate(-tt, cros);
//  cout<<"angle 2 "<<my_dir.angle(xproton1)<<endl;
//  cout<<"xphot x y z my_dir x y z "<<xproton1.x()<<" "<<xproton1.y()<<" "<<xproton1.z()<<" "
//  <<my_dir.x()<<" "<<my_dir.y()<<" "<<my_dir.z()<<" "<<xproton1.mag()<<endl;
//  fParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(my_pos.unit());
  fParticleGun->GeneratePrimaryVertex(anEvent);
  G4ThreeVector kproton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector xproton = fParticleGun->GetCurrentSource()->GetPosDist()->GetParticlePos();
//  xproton = fParticleGun->GetParticlePosition();
  G4String name = fParticleGun->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
  cout<<"!!!!!!!!!!!!!!!!!!! Z edge Z "<<G4BestUnit(Detector->GetZedge(),"Length")<<" "<<Detector->GetZedge()<<endl;
  cout<<"energy "<<G4BestUnit(fParticleGun->GetParticleEnergy(),"Energy")<<" "<<name<<endl;
  cout<<"init particle dir Px Py Pz "<<name<<" "<<kproton.x()<<" "<<kproton.y()<<" "<<kproton.z()<<endl;
  cout<<"init particle theta phi mom "<<name<<" "<<kproton.theta()<<" "<<kproton.phi()<<endl;
  cout<<"init particle theta phi pos "<<name<<" "<<xproton.theta()<<" "<<xproton.phi()<<endl;
  cout<<"init particle x y z r perp "<<name<<" "<<G4BestUnit(xproton.x(),"Length")<<" "<<G4BestUnit(xproton.y(),"Length")
                            <<G4BestUnit(xproton.z(),"Length")<<" "<<G4BestUnit(xproton.mag(),"Length")<<" "
                            <<G4BestUnit(xproton.perp(),"Length")<<endl;
  cout<<"!!!!!!!!!! R World "<<G4BestUnit(Detector->GetWorldR(),"Length")<<endl;
  cout<<"!!!!!!!!!! R  "<<G4BestUnit(Detector->GetR(),"Length")<<endl;
}
