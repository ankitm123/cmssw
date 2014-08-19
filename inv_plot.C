void inv_plot(){
  hInvMuon->SetLineColor(1);
  hInvTau->SetLineColor(4);
  hInvMuon->SetXTitle(" Mass in GeV/c^{2}");
  hInvTau->SetXTitle("Mass in GeV/c^{2}");
  hInvMuon->SetYTitle("Number of events");
  hInvTau->SetYTitle("Number of events");
  
  hInvMuon->Draw();          //draw hist_2 first as it has a larger range
  hInvTau->Draw("same");
  
  leg_hist = new TLegend(0.5,0.6,0.79,0.79);
  leg_hist->SetHeader("Invariant Mass histograms");
  leg_hist->AddEntry(hInvMuon,"hInvMuon","l");
  leg_hist->AddEntry(hInvTau,"hInvTau","l");
  leg_hist->Draw();
}
