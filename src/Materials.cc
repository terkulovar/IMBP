
#include "Materials.hh"

Materials* Materials::instance = 0;

Materials::Materials()
{
  nistMan = G4NistManager::Instance();

  nistMan->SetVerbose(2);

  CreateMaterials();
}

Materials::~Materials()
{
}


Materials* Materials::GetInstance()
{
  if (instance == 0)
    {
      instance = new Materials();
    }
  return instance;
}

G4Material* Materials::GetMaterial(const G4String material)
{
  G4Material* mat =  nistMan->FindOrBuildMaterial(material);

  if (!mat) {mat = G4Material::GetMaterial(material);}
  if (!mat) {
     std::ostringstream o;
     o << "Material " << material << " not found!";
     G4Exception("Materials::GetMaterial","",
                 FatalException,o.str().c_str());
  }

  return mat;
}

void Materials::CreateMaterials()
{

 G4String symbol;             //a=mass of a mole;
 G4double a, z, density;      //z=mean number of protons;
                             // n=number of nucleons in an isotope;
  
 G4int ncomponents, natoms;
 G4double fractionmass;

//
// define Elements
//
  
 G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);   
 G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
 G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
 G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
 G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);
 G4Element* Lu = new G4Element("Luthethium", symbol="Lu", z= 71., a= 174.97*g/mole);
 G4Element* Yt = new G4Element("Yttrium", symbol="Yt", z= 39., a= 88.9*g/mole);
 G4Element* Pb = new G4Element("Lead", symbol="Pb", z= 82., a= 207.2*g/mole);
 G4Element* F = new G4Element("Fluorine", symbol="F", z= 9., a= 18.998*g/mole);
 G4Element* K = new G4Element("Potassium", symbol="K", z= 19., a= 39.098*g/mole);
 G4Element* Na = new G4Element("Sodium", symbol="Na", z= 11., a= 22.989*g/mole);
 G4Element* Ce = new G4Element("Cerium", symbol="Ce", z= 58., a= 140.115*g/mole);
 G4Element* Fe = new G4Element("Fe", symbol="Fe", z= 26., a= 55.85*g/mole);
 G4Element* Bi = new G4Element("Bismuth", symbol="Bi", z= 83., a= 208.98*g/mole);
 G4Element* W = new G4Element("Wolframium", symbol="W", z= 74., a= 183.84*g/mole);
 G4Element* Li = new G4Element("Lithium", symbol="Li", z= 3., a= 6.94*g/mole);
  



 G4Material* C2F4 =
  new G4Material("C2F4", density=2.2*g/cm3, ncomponents=2);
 C2F4->AddElement(C, natoms=2);
 C2F4->AddElement(F, natoms=4);

 G4Material* BI2O3 =
  new G4Material("BI@O#", density=8.9*g/cm3, ncomponents=2);
 BI2O3->AddElement(Bi, natoms=2);
 BI2O3->AddElement(O, natoms=3);

 G4Material* WC =
  new G4Material("WC", density=15.77*g/cm3, ncomponents=2);
 WC->AddElement(W, natoms=1);
 WC->AddElement(O, natoms=1);

// density = 0.385*2.2+0.585*8.9+0.03*15.77;
 G4Material* Comp =
  new G4Material("Comp", 4.05*g/cm3, ncomponents=3);
 Comp->AddMaterial(C2F4, 0.385);
 Comp->AddMaterial(BI2O3, 0.585);
 Comp->AddMaterial(WC, 0.03);

 G4Material* LiF =
  new G4Material("LiF", 2.635*g/cm3, ncomponents=2);
 LiF->AddElement(Li, natoms=1);
 LiF->AddElement(F, natoms=1);

 G4Material* WATER1 =
  new G4Material("WATER1", 0.06*g/cm3, ncomponents=2);
 WATER1->AddElement(H, natoms=2);
 WATER1->AddElement(O, natoms=1);


 G4cout << *(G4Material::GetMaterialTable()) << G4endl;



}
