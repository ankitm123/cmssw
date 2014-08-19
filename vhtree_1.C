#define vhtree_cxx
#include "vhtree_1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <fstream>
void vhtree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L vhtree.C
//      Root > vhtree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   double delRMin=0.03;
   double delR=0.;
   double GenMuonPt=0.;                                                                                                            
   double delphi=0.;
   double deleta=0.;
   double RecoMuonPt=0.;
   double RecoMuonEta=0.;
   double inv_mass=0.;
   double ImpParam=0.;
   double MuonDZ=0.;
   double MuonChi2=0.;
   double SumMuonChi2=0.;
   double e1=0.;
   double e2=0.;
   double p1x=0.;
   double p2x=0.;
   double p1y=0.;
   double p2y=0.;
   double p1z=0.;
   double p2z=0.;
   double phi_1=0.;
   double phi_2=0;
   double theta_1=0;
   double theta_2=0;
   double metx=0.;
   double mety=0.;
   double p_1=0.;
   double p_2=0.;
   double x_1=0;
   double x_2=0;
   double pvis_1=0;
   double pvis_2=0;
   double inv_mass_vis=0;

   ofstream myfile("acceptance.file",ios::app);
   ofstream myfile2("acceptance2.file",ios::app);
   
   TH1F *h1PtReco = new TH1F("h1PtReco","Reco Muon Pt",100,0,300);
   TH1F *h1PtGen = new TH1F("h1PtGen","Gen Muon Pt",100,0,300);
   
   TH1F *h1PtLeading = new TH1F("h1PtLeading","Leading Muon Pt",100,0,100);
   TH1F *h1PtSubLeading = new TH1F("h1PtSubLeading","SubLeading Muon Pt",100,0,100);
   
   TH1F *h1PtLeadingi = new TH1F("h1PtLeadingi","Integrated Leading Muon Pt",100,0,100);
   TH1F *h1PtSubLeadingi = new TH1F("h1PtSubLeadingi","Integrated SubLeading Muon Pt",100,0,100);
   


   TH1F *h1EtaReco = new TH1F("h1EtaReco","Reco Muon Eta",100,-5,+5);
   //TH1F *h1EtaGen = new TH1F("h1EtaGen","Gen Muon Eta",100,-5,+5);
   
   TH1F *hPhiReco = new TH1F("hPhiReco","Difference in Phi between the two muons in multiples of #pi",100,-2,2);

   
   
   TH1F *hMCMatchEC = new TH1F("hMCMatchEC","Muon Pt resolution in End cap region",100,-3,3);
   TH1F *hMCMatchBa = new TH1F("hMCMatchBa","Muon Pt resolution in barrel region",100,-3,3);
   
   TH1F *hLeadMuon1 = new TH1F("hLeadMuon1","Leading Muon and Subleading Muon in EndCap",100,0,100);
   TH1F *hLeadMuon2 = new TH1F("hLeadMuon2","Leading Muon in Endcap, subleading muon in barrel (lead muon Pt)",100,0,100);
   TH1F *hLeadMuon3 = new TH1F("hLeadMuon3","Leading Muon in Barrel, subleading muon in endcap (lead muon Pt)",100,0,100);
   
   TH1F *hSubLeadMuon1 = new TH1F("hSubLeadMuon1","Leading Muon and Subleading Muon in EndCap",100,0,100);
   TH1F *hSubLeadMuon2 = new TH1F("hSubLeadMuon2","Leading Muon in Endcap, subleading muon in barrel (sublead muon Pt)",100,0,100);
   TH1F *hSubLeadMuon3 = new TH1F("hSubLeadMuon3","Leading Muon in Barrel, subleading muon in endcap (sublead muon Pt)",100,0,100);
   
   TH2F *h1 = new TH2F("h1","Correlation plot between Subleading and leading muon Pt",100,0,100,100,0,100);
   //TH3F *h2 = new TH3F("h2","acceptance plot",300,0,300,100,0,150,100,0,1);
   TH2F *h2 = new TH2F("h2","acceptance plot",100,0,100,100,0,100);
   //TProfile2D *h2 = new TProfile2D("h2","acceptance plot",100,0,300,100,0,300,0,1);  
   //TGraph2D *t2 = new TGraph2D();
   TH1F *hInvTau = new TH1F("hInvTau","Invariant Mass plot for Taus",100,0,300);
   TH1F *hInvMuon = new TH1F("hInvMuon","Invariant Mass plot for Muons(visible)",100,0,300);

   
   
   TH1F *hImPEC = new TH1F("hImPEC","Impact parameter plot End cap",100,-0.3,0.3);
   TH1F *hImPBa = new TH1F("hImPBa","Impact parameter plot barrel",100,-0.3,0.3);
   
   TH1F *hDZEC = new TH1F("hDZEC","Muon DZ plot End cap",100,-10,10);
   TH1F *hDZBa = new TH1F("hDZBa","Muon DZ plot barrel",100,-10,10);

   TH1F *hChiEC = new TH1F("hChiEC","Muon Chi plot End cap",100,0,16);
   TH1F *hChiBa = new TH1F("hChiBa","Muon Chi plot barrel",100,0,16);

   TH1F *hImPECEr = new TH1F("hImPECEr","Error in Impact parameter plot End cap",100,0,16);
   TH1F *hImPBaEr = new TH1F("hImPBaEr","Error in Impact parameter plot barrel",100,0,16);

   //2-D histograms

   TH2F *hEtaResEC = new TH2F("hEtaResEC","Eta vs MuonPt resolution in Endcap",100,-3,3,100,-5,5);
   //TH2F *hEtaResEC = new TH2F("hEtaResEC2","Eta vs MuonPt resolution in Endcap",100,-3,3,100,-5,5);
   TH2F *hEtaResBa = new TH2F("hEtaResBa","Eta vs MuonPt resolution in barrel",100,-3,3,100,-5,5);
   
   std::vector<double> leadmuonpt;
   std::vector<double> subleadmuonpt;
   
   


   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"Number of events "<<nentries<<endl;
   int dimuon=0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0;jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     if (nMuon!=2) continue;
     ++dimuon;
     //cout<<jentry<<" "<<nMET<<endl;  
     //cout<<nMuon<<" muons in event "<<jentry<<endl;
     //for (int i=0;i<nMuon;i++){
     phi_1=Muon_phi[0];
     phi_2=Muon_phi[1];
     theta_1=2*atan(exp(-1*Muon_eta[0]));
     theta_2=2*atan(exp(-1*Muon_eta[1]));
     metx=MET_metx[0];
     mety=MET_mety[0];
     //cout<<metx<<" "<<mety<<endl;
     p_1=(metx*sin(phi_2)-mety*cos(phi_2))/(sin(theta_1)*(cos(phi_1)-sin(phi_1)));
     p_2=(metx*sin(phi_1)-mety*cos(phi_1))/(sin(theta_2)*(cos(phi_2)-sin(phi_2)));
     pvis_1=sqrt(Muon_px[0]*Muon_px[0]+Muon_py[0]*Muon_py[0]+Muon_pz[0]*Muon_pz[0]);
     pvis_2=sqrt(Muon_px[1]*Muon_px[1]+Muon_py[1]*Muon_py[1]+Muon_pz[1]*Muon_pz[1]);
     //cout<<p_1<<endl;
     x_1=pvis_1/(pvis_1+p_1);
     x_2=pvis_2/(pvis_2+p_2);
     
     e1=Muon_energy[0];
     e2=Muon_energy[1];
     p1x=Muon_px[0];
     p2x=Muon_px[1];
     p1y=Muon_py[0];
     p2y=Muon_py[1];
     p1z=Muon_pz[0];
     p2z=Muon_pz[1];
     inv_mass_vis=sqrt((e1+e2)*(e1+e2)-((p1x+p2x)*(p1x+p2x)+(p1y+p2y)*(p1y+p2y)+(p1z+p2z)*(p1z+p2z)));
     inv_mass=inv_mass_vis/sqrt(fabs(x_1*x_2));
     hPhiReco->Fill(((Muon_phi[0]-Muon_phi[1])/(22.0/7.0)));
     //if ((fabs((Muon_phi[0]-Muon_phi[1])/(22.0/7.0))) >= 1.1 || (fabs((Muon_phi[0]-Muon_phi[1])/(22.0/7.0))) <= 0.9) cout<<(Muon_phi[0]-Muon_phi[1])/(22.0/7.0)<<" Pi"<<endl;
     if (inv_mass>0) {
       if ((fabs((Muon_phi[0]-Muon_phi[1])/(22.0/7.0))) >= 1.3 || (fabs((Muon_phi[0]-Muon_phi[1])/(22.0/7.0))) <= 0.7)hInvTau->Fill(inv_mass);
       hInvMuon->Fill(inv_mass_vis);
     }
       
     if (Muon_pt[0]>Muon_pt[1]){
       h1->Fill(Muon_pt[0],Muon_pt[1]);
       h1PtLeading->Fill(Muon_pt[0]);
       h1PtSubLeading->Fill(Muon_pt[1]);
       leadmuonpt.push_back(Muon_pt[0]);
       subleadmuonpt.push_back(Muon_pt[1]);
       
       
       //cout<<"first muon"<<endl;
       //cout<<Muon_pt[0]<<" "<<Muon_pt[1]<<endl;
     }
     
     else if (Muon_pt[0]<Muon_pt[1]){
       h1->Fill(Muon_pt[1],Muon_pt[0]);
       cout<<"second muon"<<endl;
       cout<<Muon_pt[1]<<" "<<Muon_pt[0]<<endl;
     }
     
     

     //Different cases
     if (Muon_pt[0]>Muon_pt[1]){
	 if (((fabs(Muon_eta[0]) > 1.6) && (fabs(Muon_eta[0]) < 2.2)) && ((fabs(Muon_eta[1]) > 1.6) && (fabs(Muon_eta[1]) < 2.2))){
	   hLeadMuon1->Fill(Muon_pt[0]);
	   hSubLeadMuon1->Fill(Muon_pt[1]);
	 }
	 if (((fabs(Muon_eta[0]) > 1.6) && (fabs(Muon_eta[0]) < 2.2)) && (fabs(Muon_eta[1]) < 1.6)){
	   hLeadMuon2->Fill(Muon_pt[0]);
	   hSubLeadMuon2->Fill(Muon_pt[1]);
	 }
	 if ((fabs(Muon_eta[0]) < 1.6) && ((fabs(Muon_eta[1]) > 1.6) && (fabs(Muon_eta[1]) < 2.2))){
	   hLeadMuon3->Fill(Muon_pt[0]);
           hSubLeadMuon3->Fill(Muon_pt[1]);
	 }
      
	   
	 //cout<<"Leading Muon Pt "<<Muon_pt[0]<<" "<<fabs(Muon_eta[0])<<endl;
	 //cout<<"Sub-Leading Muon Pt "<<Muon_pt[1]<<" "<<fabs(Muon_eta[1])<<endl;
	 
       } 
       else{
	 //cout<<"hi"<<endl;
	 //cout<<"Leading Muon Pt "<<Muon_pt[1]<<" "<<fabs(Muon_eta[1])<<endl;
	 //cout<<"Sub-Leading Muon Pt "<<Muon_pt[0]<<" "<<fabs(Muon_eta[0])<<endl;
	 if (((fabs(Muon_eta[1]) > 1.6) && (fabs(Muon_eta[1]) < 2.2)) && ((fabs(Muon_eta[0]) > 1.6) && (fabs(Muon_eta[0]) < 2.2))){
	   hLeadMuon1->Fill(Muon_pt[1]);
	   hSubLeadMuon1->Fill(Muon_pt[0]);
	 }
	 if (((fabs(Muon_eta[1]) > 1.6) && (fabs(Muon_eta[1]) < 2.2)) && (fabs(Muon_eta[0]) < 1.6)){
	   hLeadMuon2->Fill(Muon_pt[1]);
	   hSubLeadMuon2->Fill(Muon_pt[0]);
	 }
	 if ((fabs(Muon_eta[1]) < 1.6) && ((fabs(Muon_eta[0]) > 1.6) && (fabs(Muon_eta[0]) < 2.2))){
	   hLeadMuon3->Fill(Muon_pt[1]);
           hSubLeadMuon3->Fill(Muon_pt[0]);
	 }
      
       }
   }
   cout<<"Number of dimuon events "<<dimuon<<" BR of di tau to dimuon "<<((float)dimuon/(float)nentries)*100<<endl;
   

   Double_t weight=0.;
   h1PtLeading->ComputeIntegral();
   Double_t *integral = (h1PtLeading->GetIntegral());
   h1PtLeadingi->SetContent(integral);
   Int_t binsize=h1PtLeadingi->GetSize();
   for (int i=1;i<binsize-1;i++){
     Double_t entry=h1PtLeadingi->GetBinContent(i);
     h1PtLeadingi->SetBinContent(i,1-entry);
   }

   h1PtSubLeading->ComputeIntegral();
   Double_t *integral = (h1PtSubLeading->GetIntegral());
   h1PtSubLeadingi->SetContent(integral);
   binsize=h1PtSubLeadingi->GetSize();
  
   for (int i=1;i<binsize-1;i++){
     entry=h1PtSubLeadingi->GetBinContent(i);
     h1PtSubLeadingi->SetBinContent(i,1-entry);
   }
  //myfile<<"Leading muon pt"
   for (int i=1;i<leadmuonpt.size();i++){
     for (int j=1;j<subleadmuonpt.size();j++){
       weight=0;
       int k=h1PtLeadingi->GetXaxis()->FindBin(leadmuonpt[i]);
       int m=h1PtSubLeadingi->GetXaxis()->FindBin(subleadmuonpt[j]);
       if (k!=101 && m!=101){
	 if (h1PtLeadingi->GetBinContent(k) < h1PtSubLeadingi->GetBinContent(m))
	   weight=h1PtLeadingi->GetBinContent(k);
	 else if (h1PtLeadingi->GetBinContent(k) > h1PtSubLeadingi->GetBinContent(m))
	   weight=h1PtSubLeadingi->GetBinContent(m);
	 else
	   weight=h1PtSubLeadingi->GetBinContent(m);
	 //if (weight>0.8)myfile2<<leadmuonpt[i]<<" bin "<<k<<" "<<subleadmuonpt[j]<<" bin "<<m<<" "<<weight<<endl;
	 if (leadmuonpt[i]>subleadmuonpt[j]){
	   //myfile2<<leadmuonpt[i]<<" "<<subleadmuonpt[j]<<" "<<weight<<endl;
	   h2->Fill(leadmuonpt[i],subleadmuonpt[j],weight);
	 }
       }
     }
   }

   myfile2.close();


}



     /*
     for (int i=0;i<nMuon; i++){
       RecoMuonPt=Muon_pt[i];
       RecoMuonEta=Muon_eta[i];
       GenMuonPt=Muon_ptgen[i];
       ImpParam=Muon_trkD0[i];
       MuonDZ=Muon_trkDz[i];
       MuonChi2=Muon_globalChi2[i];
     
       if ((fabs(RecoMuonEta) > 1.6) && (fabs(RecoMuonEta) < 2.2)){
	 hMCMatchEC->Fill(RecoMuonPt-GenMuonPt);
	 hImPEC->Fill(ImpParam);
	 hDZEC->Fill(MuonDZ);
	 hChiEC->Fill(MuonChi2);
	 hEtaResEC->Fill(RecoMuonPt-GenMuonPt,RecoMuonEta);
       }
       if ((fabs(Muon_eta[i]) < 1.6)){
	 hMCMatchBa->Fill(RecoMuonPt-GenMuonPt);
	 hImPBa->Fill(ImpParam);
	 hDZBa->Fill(MuonDZ);
	 hChiBa->Fill(MuonChi2);
	 hEtaResBa->Fill(RecoMuonPt-GenMuonPt,RecoMuonEta);
       }
       h1PtGen->Fill(GenMuonPt);
       h1PtReco->Fill(RecoMuonPt);
       //h1EtaGen->Fill(GenParticle_eta[i]);
       h1EtaReco->Fill(RecoMuonEta);
     }

     for (int i = 0; i<nTau; i++){
       for (int j=0; j<nTau;j++){
	 if (Tau_charge[i]==1 && Tau_charge[j]==-1){
	   //cout<<Muon_charge[i]<<" "<<Muon_charge[j]<<endl;
	   metp1=0.;
	   metp2=0.;
	   for (int k=0; k<nMet; k++){
	     for (int l=0; l<nMet; l++){
	       if (k!=l){
		 metp1=Met_p[k];
		 metp2=Met_p[l];
	       }
	     }
	   } 
	   e1=Tau_energy[i];
	   e2=Tau_energy[j];
	   p1x=Tau_px[i];
	   p2x=Tau_px[j];
	   p1y=Tau_py[i];
	   p2y=Tau_py[j];
	   p1z=Tau_pz[i];
	   p2z=Tau_pz[j];
	   inv_mass=(e1+e2)*(e1+e2)-((p1x+p2x)*(p1x+p2x)+(p1y+p2y)*(p1y+p2y)+(p1z+p2z)*(p1z+p2z));
	   //cout<<e1<<" "<<e2<<" "<<" "<<p1x<<" "<<p2x<<" "<<inv_mass<<endl;
	   if (inv_mass>0) hInvMuon->Fill(sqrt(inv_mass));
	 }
       }
     }

   }
}
     */
