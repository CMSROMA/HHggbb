void plotsPuId() {

  gStyle->SetOptStat(0);
  
  TFile myFile("../finalizedTrees_Radion/Radion_Radion_M-300_madgraph_default_CSV.root");
  
  TH1F *myEffEtaAss = (TH1F*)myFile.Get("myEffEtaAss");
  TH1F *myEffPtAss  = (TH1F*)myFile.Get("myEffPtAss");

  TH1F *myEffEtaBtagAss = (TH1F*)myFile.Get("myEffEtaBtagAss");
  TH1F *myEffPtBtagAss  = (TH1F*)myFile.Get("myEffPtBtagAss");

  TH1F *myNAssRMS   = (TH1F*)myFile.Get("myNAssRMS");
  TH1F *myAssRMS    = (TH1F*)myFile.Get("myAssRMS");
  TH1F *myNAssBetaS = (TH1F*)myFile.Get("myNAssBetaS");
  TH1F *myAssBetaS  = (TH1F*)myFile.Get("myAssBetaS");

  // cosmetics
  myNAssRMS   -> SetLineWidth(2);    myNAssRMS   -> SetLineColor(2);  
  myAssRMS    -> SetLineWidth(2);    myAssRMS    -> SetLineColor(4);  
  myNAssBetaS -> SetLineWidth(2);    myNAssBetaS -> SetLineColor(2);  
  myAssBetaS  -> SetLineWidth(2);    myAssBetaS  -> SetLineColor(4);  

  // titles
  myNAssRMS   -> SetTitle("");    myNAssRMS   -> GetXaxis()->SetTitle("RMS");
  myAssRMS    -> SetTitle("");    myAssRMS    -> GetXaxis()->SetTitle("RMS");
  myNAssBetaS -> SetTitle("");    myNAssBetaS -> GetXaxis()->SetTitle("#beta*");
  myAssBetaS  -> SetTitle("");    myAssBetaS  -> GetXaxis()->SetTitle("#beta*");

  // plots
  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0); 
  leg->SetBorderSize(0); 
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(myNAssRMS, "PU jet", "l");
  leg->AddEntry(myAssRMS,  "good jet", "l");

  TLegend *leg2;
  leg2 = new TLegend(0.6,0.6,0.85,0.85);
  leg2->SetFillStyle(0); 
  leg2->SetBorderSize(0); 
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(myEffEtaAss,     "jets", "l");
  leg2->AddEntry(myEffEtaBtagAss, "btagged jets", "l");

  TCanvas c2("c2","c2",1); 
  c2.cd();
  myAssRMS->DrawNormalized(); 
  myNAssRMS->DrawNormalized("same"); 
  leg->Draw();
  c2.SaveAs("myRMS.png");
  c2.SaveAs("myRMS.eps");

  TCanvas c3("c3","c3",1); 
  c3.cd();
  myAssBetaS->DrawNormalized(); 
  myNAssBetaS->DrawNormalized("same");
  leg->Draw();
  c3.SaveAs("myBetaS.png");
  c3.SaveAs("myBetaS.eps");

  TH2F *myHeffEta = new TH2F("myHeffEta","",100,-2.5,2.5,100,0.5,1.1);
  myHeffEta->GetXaxis()->SetTitle("#eta");
  myHeffEta->GetYaxis()->SetTitle("efficiency");

  TH2F *myHeffPt  = new TH2F("myHeffPt", "", 100,25.,100.,100,0.5,1.1);
  myHeffPt->GetXaxis()->SetTitle("p_{T} [GeV]");
  myHeffPt->GetYaxis()->SetTitle("efficiency");

  TCanvas c51("c51","c51",1);  
  myHeffEta   -> Draw();
  myEffEtaAss -> SetLineColor(4);
  myEffEtaAss -> SetLineWidth(2);
  myEffEtaAss -> Draw("same");
  myEffEtaBtagAss -> SetLineColor(2);
  myEffEtaBtagAss -> SetLineWidth(2);
  myEffEtaBtagAss -> Draw("same");
  leg2->Draw();
  c51.SaveAs("myEffEtaAss.png");
  c51.SaveAs("myEffEtaAss.root");

  TCanvas c61("c61","c61",1); 
  myHeffPt   -> Draw();
  myEffPtAss -> SetLineWidth(2);
  myEffPtAss -> SetLineColor(4);
  myEffPtAss -> Draw("same"); 
  myEffPtBtagAss -> SetLineWidth(2);
  myEffPtBtagAss -> SetLineColor(2);
  myEffPtBtagAss -> Draw("same");
  leg2->Draw();
  c61.SaveAs("myEffPtAss.png");
  c61.SaveAs("myEffPtAss.root");
}

