void plot_eta_pt(){
  hEtaPt10->SetLineColor(1);
  hEtaPt15->SetLineColor(2);
  hEtaPt20->SetLineColor(3);
  hEtaPt25->SetLineColor(4);
  hEtaPt30->SetLineColor(8);
  hEtaPt35->SetLineColor(6);

  
  hEtaPt10->SetStats(0);
  hEtaPt10->SetTitle("#eta as a function of p_{T} threshold");

  hEtaPt10->SetXTitle("#eta_{#mu}^{Reco}");
  hEtaPt10->SetYTitle("Events");
  
  hEtaPt10->Draw();
  hEtaPt15->Draw("SAME");
  hEtaPt20->Draw("SAME");
  hEtaPt25->Draw("SAME");
  hEtaPt30->Draw("SAME");
  hEtaPt35->Draw("SAME");
  
  
  leg_hist = new TLegend(0.6724138,0.684322,0.8520115,0.8728814);
  //leg_hist->SetHeader("Pt histograms");
  leg_hist->SetFillColor(0);
  leg_hist->AddEntry(hEtaPt10,"p_{T} > 10 GeV/c","l");
  leg_hist->AddEntry(hEtaPt15,"p_{T} > 15 GeV/c","l");
  leg_hist->AddEntry(hEtaPt20,"p_{T} > 20 GeV/c","l");
  leg_hist->AddEntry(hEtaPt25,"p_{T} > 25 GeV/c","l");
  leg_hist->AddEntry(hEtaPt30,"p_{T} > 30 GeV/c","l");
  leg_hist->AddEntry(hEtaPt35,"p_{T} > 35 GeV/c","l");
  leg_hist->Draw();
  c1->SaveAs("plot_eta_pt.png");
}
