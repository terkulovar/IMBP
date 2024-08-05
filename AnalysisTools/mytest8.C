#include <stdio.h>
#include <iostream>
//#include <string>
//#include <cstring>
void mytest8()
{
 Double_t pi = acos(0.0)*2.0;
 auto c1 = new TCanvas("c1","A Simple Graph Example",200,50,700, 600);
 c1->SetFillColor(0);
 c1->SetBorderMode(0);
// c1->SetGrid();
 gPad->SetFrameFillColor(0);
 gPad->SetFrameFillStyle(4000);
 gStyle->SetOptStat(0);
 auto hh = new TH1F("hh","",160,0,160);
 hh->SetMarkerStyle(20);
 hh->SetMarkerSize(0.8);
 hh->SetMarkerColor(2);
 hh->SetXTitle("depth in water [mm]");
 hh->SetYTitle("dE [MeV]");
 hh->SetTitle("water, proton 100MeV");
 cout<<"hh bins width center "<<hh->GetNbinsX()<<" "<<hh->GetBinWidth(1)<<" "<<hh->GetBinCenter(1)<<endl;
 auto f = new TFile("../build/output.root");
 auto ftree = (TTree*)f->Get("Dose");
 auto ftree1 = (TTree*)f->Get("Event");
 auto ftree2 = (TTree*)f->Get("Dose1");
 ftree->Print();
 Int_t evId;
 Double_t en;
 Double_t energy;
 Double_t edep, zpos, dose, ziv;
 Double_t dose1, ziv1;
 char name[20];
// char* name;
// std::string name;
 char name1[20];
 char name2[20];
// char* name1;
// char* name2;
 Double_t mass = 2.*2.*0.5/1000.;
 Double_t d0, d1, d2, d3, d4, d5, d6, d7, d8, d9;
 Double_t ed0, ed1, ed2, ed3, ed4, ed5, ed6, ed7, ed8, ed9;
 Double_t z0, z1, z2, z3, z4, z5, z6, z7, z8, z9;
 Double_t enav = 0;
 Double_t days = 1.;
 d0=0; d1=0; d2=0; d3=0; d4=0; d5=0; d6=0; d7=0; d8=0; d9=0;
 ed0=0; ed1=0; ed2=0; ed3=0; ed4=0; ed5=0; ed6=0; ed7=0; ed8=0; ed9=0;
 z0=0; z1=0; z2=0; z3=0; z4=0; z5=0; z6=0; z7=0; z8=0; z9=0;
 Double_t d11 = 0, z11 = 0;
// Double_t square = 9.;
// Double_t square = 4.0*pi*100.0*100.0;
// Double_t square = 4.0*pi*620.*620.;
 Double_t square = 2.0*pi*150.*(300.0+150.0)/2.0;
// Double_t square = 500.*500.*1.;
 Double_t spectrum_trapped_protons_max  = 2.27e+07;
 Double_t spectrum_trapped_protons_min  = 5.11e+07;
 Double_t spectrum_all_min  = 5.11e+07+4.89E+04+7.01E+03+1.98E+02+1.90E+02;
 Double_t spectrum_GCR_min  = 4.89E+04+7.01E+03+1.98E+02+1.90E+02;
 Double_t spenvis_trapped_protons_2010_min  = 789696;
 Double_t spenvis_trapped_protons_2010_368km_min  = 22982400;
 Double_t spenvis_trapped_protons_2010_368km_min_new  = 22302432;
 Double_t spenvis_trapped_protons_2010_368km_min_1Mev  = 4464288.;
 Double_t spenvis_trapped_protons_2010_368km_52degree_min_1Mev  = 4.03237e+06;
// Double_t spenvis_trapped_protons_2010_368km_52degree_min_1Mev  = 3.03237e+06;
 Double_t spenvis_trapped_protons_2010_min_1Mev  = 151099;
 Double_t spectrum_trapped_protons_min_1Mev  = 9.95E+06;
 Double_t spectrum_trapped_protons_max_1Mev  = 2.24E+06;
 Double_t spectrum_trapped_protons_min_10Mev  = 1.87E+06;
 Double_t spectrum_GCR_H_min  = 4.89E+04;
 Double_t spectrum_GCR_He_min  = 7.01E+03;
 Double_t spectrum_GCR_Li_min  = 3.03E+01;
 Double_t spenvis_trapped_protons_min  = 2.36313e+08;
 Double_t spenvis_trapped_protons_max  = 9.15764e+07;
 Double_t spenvis_GCR_min_368km_H = 47536.5;
 Double_t spenvis_GCR_min_368km_He = 6878.22;
 Double_t norm = spectrum_all_min;
// Double_t norm = spectrum_trapped_protons_min_10Mev;
 ftree1->SetBranchAddress("Energy",&energy);

 ftree2->SetBranchAddress("dose",&dose1);
 ftree2->SetBranchAddress("dose_zi",&ziv1);

 ftree->SetBranchAddress("Event",&evId);
 ftree->SetBranchAddress("En",&en);
 ftree->SetBranchAddress("name",&name);
 ftree->SetBranchAddress("particle",&name2);
 ftree->SetBranchAddress("zpos",&zpos);
 ftree->SetBranchAddress("Edep",&edep);
 ftree->SetBranchAddress("Dose_Gr",&dose);
 ftree->SetBranchAddress("Dose_Zi",&ziv);
 int nentries = ftree->GetEntriesFast();
 int nentries1 = ftree1->GetEntriesFast();
 int nentries2 = ftree2->GetEntriesFast();
 cout<<"entries2 "<<nentries2<<endl;
 Int_t iold = -1;
 Int_t iold1 = -1;
 Int_t iold2 = -1;
 Int_t iold3 = -1;
 Int_t iold4 = -1;
 Int_t iold5 = -1;
 Int_t iold6 = -1;
 Int_t iold7 = -1;
 Int_t iold8 = -1;
 Int_t iold9 = -1;
 Int_t itot = 0;
 Int_t i1 = 0;
 Int_t nevents = 0;
 Int_t nevents1 = 0;
 Int_t nevents2 = 0;
 Int_t nevents3 = 0;
 Int_t nevents4 = 0;
 Int_t nevents5 = 0;
 Int_t nevents6 = 0;
 Int_t nevents7 = 0;
 Int_t nevents8 = 0;
 Int_t nevents9 = 0;

 for ( Int_t i = 0; i < nentries1; i++ )
 {
  ftree1->GetEntry(i);
  enav += energy/nentries1;
//  cout<<"energy "<<energy<<" "<<name1<<endl;
 }
 for ( Int_t i = 0; i < nentries2; i++ )
 {
  ftree2->GetEntry(i);
  d11 += dose1;
  z11 += ziv1;
 }
 cout<<"d11 z11  "<<d11*norm*square*days<<" "<<z11*norm*square*days<<endl;



 for ( Int_t i = 0; i < nentries; i++ )
 {
  ftree->GetEntry(i);
  if ( i%10000 == 0 )
   cout<<"i dose ziv z "<<i<<" "<<dose<<" "<<edep<<" "<<ziv<<" "<<zpos<<" "<<evId<<" "<<name<<" "<<name2<<endl;
  if ( strcmp(name,"LiF" ) == 0 ) { d0 += dose; ed0 += edep; z0 += ziv; }
  if ( strcmp(name,"LiF1" ) == 0 ) { d1 += dose; ed1 += edep; z1 += ziv; }
  if ( strcmp(name,"LiF2" ) == 0 ) { d2 += dose; ed2 += edep; z2 += ziv; }
  if ( strcmp(name,"LiF3" ) == 0 ) { d3 += dose; ed3 += edep; z3 += ziv; }
  if ( strcmp(name,"LiF4" ) == 0 ) { d4 += dose; ed4 += edep; z4 += ziv; }
  if ( strcmp(name,"LiF5" ) == 0 ) { d5 += dose; ed5 += edep; z5 += ziv; }
  if ( strcmp(name,"LiF6" ) == 0 ) { d6 += dose; ed6 += edep; z6 += ziv; }
  if ( strcmp(name,"LiF7" ) == 0 ) { d7 += dose; ed7 += edep; z7 += ziv; }
  if ( strcmp(name,"LiF8" ) == 0 ) { d8 += dose; ed8 += edep; z8 += ziv; }
  if ( strcmp(name,"LiF9" ) == 0 ) { d9 += dose; ed9 += edep; z9 += ziv; }
  if ( strcmp(name,"LiF9" ) == 0 ) { cout<<"LiF9 "<<dose<<" "<<edep<<" "<<z0<<endl; }

  if ( evId != iold && strcmp(name,"LiF" ) == 0 )
  {
   cout<<"new event "<<evId<<endl;
   nevents++;
//   if ( i1 == 0 ) { hh->Fill(
   iold = evId;
  }
  if ( evId != iold1 && strcmp(name,"LiF1" ) == 0 )
  {
   nevents1++;
   iold1 = evId;
  }
  if ( evId != iold2 && strcmp(name,"LiF2" ) == 0 )
  {
   nevents2++;
   iold2 = evId;
  }
  if ( evId != iold3 && strcmp(name,"LiF3" ) == 0 )
  {
   nevents3++;
   iold3 = evId;
  }
  if ( evId != iold4 && strcmp(name,"LiF4" ) == 0 )
  {
   nevents4++;
   iold4 = evId;
  }
  if ( evId != iold5 && strcmp(name,"LiF5" ) == 0 )
  {
   nevents5++;
   iold5 = evId;
  }
  if ( evId != iold6 && strcmp(name,"LiF6" ) == 0 )
  {
   nevents6++;
   iold6 = evId;
  }
  if ( evId != iold7 && strcmp(name,"LiF7" ) == 0 )
  {
   nevents7++;
   iold7 = evId;
  }
  if ( evId != iold8 && strcmp(name,"LiF8" ) == 0 )
  {
   nevents8++;
   iold8 = evId;
  }
  if ( evId != iold9 && strcmp(name,"LiF9" ) == 0 )
  {
   nevents9++;
   iold9 = evId;
  }
//  enav += en;
//  hh->Fill(zpos,dose/1.);
//  hh1->Fill(zpos,ziv/1.);
 }
 cout<<"events "<<nevents<<endl;
 d0 = d0*norm*square*days;
 d1 = d1*norm*square*days;
 d2 = d2*norm*square*days;
 d3 = d3*norm*square*days;
 d4 = d4*norm*square*days;
 d5 = d5*norm*square*days;
 d6 = d6*norm*square*days;
 d7 = d7*norm*square*days;
 d8 = d8*norm*square*days;
 d9 = d9*norm*square*days;
 z0 = z0*norm*square*days;
 z1 = z1*norm*square*days;
 z2 = z2*norm*square*days;
 z3 = z3*norm*square*days;
 z4 = z4*norm*square*days;
 z5 = z5*norm*square*days;
 z6 = z6*norm*square*days;
 z7 = z7*norm*square*days;
 z8 = z8*norm*square*days;
 z9 = z9*norm*square*days;
// ed1 = ed1*norm*square*days;
// ed2 = ed2*norm*square*days;
// ed3 = ed3*norm*square*days;
// ed4 = ed4*norm*square*days;
// ed5 = ed5*norm*square*days;
 cout<<"average energy Mev "<<enav<<endl;
 cout<<"doses Gr "<<d0<<" "<<z0<<" "<<d1<<" "<<z1<<" "<<d2<<" "<<z2<<" "<<d3<<" "<<z3<<" "<<d4
                  <<" "<<z4<<" "<<d5<<" "<<z5<<" "<<d6<<" "<<z6<<" "<<d7<<" "<<z7<<endl;
 cout<<"edep Mev "<<ed0<<" "<<ed1<<" "<<ed2<<" "<<ed3<<" "<<ed4<<" "<<ed5<<" "<<ed6<<" "<<ed7<<endl;
 ofstream *out;
// out = new ofstream("spectrum_GCR_H_min_composite_water.txt", ios::app);
// out = new ofstream("spectrum_trapped_protons_min_composite_water_sphere_R_1m.txt", ios::app);
 out = new ofstream("spectrum_all_protons_He_C_O__min_no_yes_composite_iso_cylinder_L_300cm_0.5-10cm_polyprop.txt", ios::app);
 *out<<d0<<" "<<z0<<" "<<ed0<<" "<<d1<<" "<<z1<<" "<<ed1<<" "<<d2<<" "<<z2<<" "<<ed2<<" "<<d3<<" "<<z3<<" "<<ed3<<" "
 <<d4<<" "<<z4<<" "<<ed4<<" "<<d5<<" "<<z5<<" "<<ed5<<" "<<d6<<" "<<z6<<" "<<ed6<<" "<<d7<<" "<<z7<<" "<<ed7<<" "
 <<" "<<d8<<" "<<z8<<" "<<ed8<<" "" "<<d9<<" "<<z9<<" "<<ed9<<" "<<
  nevents<<" "<<nevents1<<" "<<nevents2<<" "<<nevents3<<" "<<nevents4<<" "<<nevents5<<" "<<nevents6<<" "<<nevents7
  <<" "<<nevents8<<" "<<nevents9<<endl;
 out->close();
 delete out;
}
