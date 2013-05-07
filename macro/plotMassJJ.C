#define plotMassJJ_cxx
#include "plotMassJJ.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

void plotMassJJ::Loop() {

  if (fChain == 0) return;

  // histos, 1 btag
  TH1F *H1_mjj_gen        = new TH1F("H1_mjj_gen",       "H1_mjj_gen",       25,0.,250.);
  TH1F *H1_mjj_pt         = new TH1F("H1_mjj_pt",        "H1_mjj_pt",        25,0.,250.);
  TH1F *H1_mjj_maxptjj    = new TH1F("H1_mjj_maxptjj",   "H1_mjj_maxptjj",   25,0.,250.);
  TH1F *H1_mjj_maxptmjj   = new TH1F("H1_mjj_maxptmjj",  "H1_mjj_maxptmjj",  25,0.,250.);
  TH1F *H1_mjj_mgg        = new TH1F("H1_mjj_mgg",       "H1_mjj_mgg",       25,0.,250.);
  // 
  TH1F *H1_mjj_btag_pt    = new TH1F("H1_mjj_btag_pt",   "H1_mjj_btag_pt",   25,0.,250.);
  TH1F *H1_mjj_btag_ptjj  = new TH1F("H1_mjj_btag_ptjj", "H1_mjj_btag_ptjj", 25,0.,250.);
  TH1F *H1_mjj_btag_ptmjj = new TH1F("H1_mjj_btag_ptmjj","H1_mjj_btag_ptmjj",25,0.,250.);
  TH1F *H1_mjj_btag_mgg   = new TH1F("H1_mjj_btag_mgg",  "H1_mjj_btag_mgg",  25,0.,250.);
  // 
  TH1F *H1_dMmin_btag     = new TH1F("H1_dMmin_btag",    "H1_dMmin_btag",    25,0.,300.);     
  TH1F *H1_dMmin          = new TH1F("H1_dMmin",         "H1_dMmin",         25,0.,300.);     

  // histos, 2 btag
  TH1F *H2_mjj_gen        = new TH1F("H2_mjj_gen",       "H2_mjj_gen",       25,0.,250.);
  TH1F *H2_mjj_pt         = new TH1F("H2_mjj_pt",        "H2_mjj_pt",        25,0.,250.);
  TH1F *H2_mjj_maxptjj    = new TH1F("H2_mjj_maxptjj",   "H2_mjj_maxptjj",   25,0.,250.);
  TH1F *H2_mjj_maxptmjj   = new TH1F("H2_mjj_maxptmjj",  "H2_mjj_maxptmjj",  25,0.,250.);
  TH1F *H2_mjj_mgg        = new TH1F("H2_mjj_mgg",       "H2_mjj_mgg",       25,0.,250.);
  //
  TH1F *H2_mjj_btag_pt    = new TH1F("H2_mjj_btag_pt",   "H2_mjj_btag_pt",   25,0.,250.);
  TH1F *H2_mjj_btag_ptjj  = new TH1F("H2_mjj_btag_ptjj", "H2_mjj_btag_ptjj", 25,0.,250.);
  TH1F *H2_mjj_btag_ptmjj = new TH1F("H2_mjj_btag_ptmjj","H2_mjj_btag_ptmjj",25,0.,250.);
  TH1F *H2_mjj_btag_mgg   = new TH1F("H2_mjj_btag_mgg",  "H2_mjj_btag_mgg",  25,0.,250.);
  // 
  TH1F *H2_dMmin_btag     = new TH1F("H2_dMmin_btag",    "H2_dMmin_btag",    25,0.,300.);     
  TH1F *H2_dMmin          = new TH1F("H2_dMmin",         "H2_dMmin",         25,0.,300.);     
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (btagCategory==1) {
      H1_mjj_gen        -> Fill(mjj_gen,        weight);
      H1_mjj_pt         -> Fill(mjj_pt,         weight);
      H1_mjj_maxptjj    -> Fill(mjj_maxptjj,    weight);
      H1_mjj_maxptmjj   -> Fill(mjj_maxptmjj,   weight);
      H1_mjj_mgg        -> Fill(mjj_mgg,        weight);
      // 
      H1_mjj_btag_pt    -> Fill(mjj_btag_pt,    weight);
      H1_mjj_btag_ptjj  -> Fill(mjj_btag_ptjj,  weight);
      H1_mjj_btag_ptmjj -> Fill(mjj_btag_ptmjj, weight);
      H1_mjj_btag_mgg   -> Fill(mjj_btag_mgg,   weight);
      //
      H1_dMmin          -> Fill(dMmin);
      H1_dMmin_btag     -> Fill(dMmin_btag);

    } else if (btagCategory==2) {
      H2_mjj_gen        -> Fill(mjj_gen,        weight);
      H2_mjj_pt         -> Fill(mjj_pt,         weight);
      H2_mjj_maxptjj    -> Fill(mjj_maxptjj,    weight);
      H2_mjj_maxptmjj   -> Fill(mjj_maxptmjj,   weight);
      H2_mjj_mgg        -> Fill(mjj_mgg,        weight);
      //
      H2_mjj_btag_pt    -> Fill(mjj_btag_pt,    weight);
      H2_mjj_btag_ptjj  -> Fill(mjj_btag_ptjj,  weight);
      H2_mjj_btag_ptmjj -> Fill(mjj_btag_ptmjj, weight);
      H2_mjj_btag_mgg   -> Fill(mjj_btag_mgg,   weight);
      // 
      H2_dMmin          -> Fill(dMmin);
      H2_dMmin_btag     -> Fill(dMmin_btag);
    }
  }

  // cosmetics
  H1_mjj_gen        -> SetLineColor(1);  
  H1_mjj_pt         -> SetLineColor(2);  
  H1_mjj_maxptjj    -> SetLineColor(3);
  H1_mjj_maxptmjj   -> SetLineColor(4); 
  H1_mjj_mgg        -> SetLineColor(6); 
  H1_mjj_btag_pt    -> SetLineColor(2);  H1_mjj_btag_pt    -> SetLineStyle(2); 
  H1_mjj_btag_ptjj  -> SetLineColor(3);  H1_mjj_btag_ptjj  -> SetLineStyle(2); 
  H1_mjj_btag_ptmjj -> SetLineColor(4);  H1_mjj_btag_ptmjj -> SetLineStyle(2); 
  H1_mjj_btag_mgg   -> SetLineColor(6);  H1_mjj_btag_mgg   -> SetLineStyle(2); 
  H1_dMmin          -> SetLineColor(1);    
  H1_dMmin_btag     -> SetLineColor(2);   
  //
  H1_mjj_gen        -> SetLineWidth(2);  
  H1_mjj_pt         -> SetLineWidth(2);  
  H1_mjj_maxptjj    -> SetLineWidth(2);
  H1_mjj_maxptmjj   -> SetLineWidth(2); 
  H1_mjj_mgg        -> SetLineWidth(2); 
  H1_mjj_btag_pt    -> SetLineWidth(2);  
  H1_mjj_btag_ptjj  -> SetLineWidth(2);  
  H1_mjj_btag_ptmjj -> SetLineWidth(2);  
  H1_mjj_btag_mgg   -> SetLineWidth(2);  
  H1_dMmin          -> SetLineWidth(2);    
  H1_dMmin_btag     -> SetLineWidth(2);   
  //
  H2_mjj_gen        -> SetLineColor(1);  
  H2_mjj_pt         -> SetLineColor(2);
  H2_mjj_maxptjj    -> SetLineColor(3);
  H2_mjj_maxptmjj   -> SetLineColor(4); 
  H2_mjj_mgg        -> SetLineColor(6); 
  H2_mjj_btag_pt    -> SetLineColor(2);  H2_mjj_btag_pt    -> SetLineStyle(2);
  H2_mjj_btag_ptjj  -> SetLineColor(3);  H2_mjj_btag_ptjj  -> SetLineStyle(2); 
  H2_mjj_btag_ptmjj -> SetLineColor(4);  H2_mjj_btag_ptmjj -> SetLineStyle(2);
  H2_mjj_btag_mgg   -> SetLineColor(6);  H2_mjj_btag_mgg   -> SetLineStyle(2);
  H2_dMmin          -> SetLineColor(1);   
  H2_dMmin_btag     -> SetLineColor(2);   
  //
  H2_mjj_gen        -> SetLineWidth(2);  
  H2_mjj_pt         -> SetLineWidth(2);  
  H2_mjj_maxptjj    -> SetLineWidth(2);
  H2_mjj_maxptmjj   -> SetLineWidth(2); 
  H2_mjj_mgg        -> SetLineWidth(2); 
  H2_mjj_btag_pt    -> SetLineWidth(2);  
  H2_mjj_btag_ptjj  -> SetLineWidth(2);  
  H2_mjj_btag_ptmjj -> SetLineWidth(2);  
  H2_mjj_btag_mgg   -> SetLineWidth(2);  
  H2_dMmin          -> SetLineWidth(2);    
  H2_dMmin_btag     -> SetLineWidth(2);   

  // titles
  H1_mjj_gen        -> SetTitle("");   H2_mjj_gen        -> SetTitle("");
  H1_mjj_pt         -> SetTitle("");   H2_mjj_pt         -> SetTitle("");
  H1_mjj_maxptjj    -> SetTitle("");   H2_mjj_maxptjj    -> SetTitle("");
  H1_mjj_maxptmjj   -> SetTitle("");   H2_mjj_maxptmjj   -> SetTitle(""); 
  H1_mjj_mgg        -> SetTitle("");   H2_mjj_mgg        -> SetTitle(""); 
  H1_mjj_btag_pt    -> SetTitle("");   H2_mjj_btag_pt    -> SetTitle("");
  H1_mjj_btag_ptjj  -> SetTitle("");   H2_mjj_btag_ptjj  -> SetTitle("");
  H1_mjj_btag_ptmjj -> SetTitle("");   H2_mjj_btag_ptmjj -> SetTitle("");
  H1_mjj_btag_mgg   -> SetTitle("");   H2_mjj_btag_mgg   -> SetTitle("");
  H1_dMmin          -> SetTitle("");   H2_dMmin          -> SetTitle("");   
  H1_dMmin_btag     -> SetTitle("");   H2_dMmin_btag     -> SetTitle("");   

  // histos
  TH1F *myHS = new TH1F("myHS","",100,0,250);
  myHS->SetMaximum(0.25);
  TH1F *myHB = new TH1F("myHB","",100,0,250);
  myHB->SetMaximum(0.15);    
  myHS->GetXaxis()->SetTitle("mjj [GeV]");
  myHB->GetXaxis()->SetTitle("mjj [GeV]");

  // legend
  TLegend *leg1;
  leg1 = new TLegend(0.6,0.6,0.85,0.85);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.05);
  leg1->SetFillColor(0);
  leg1->AddEntry(H2_mjj_gen,        "best MC truth match",   "l");  
  leg1->AddEntry(H2_mjj_pt,         "2 highest pt jets",   "l");  
  leg1->AddEntry(H2_mjj_maxptjj,    "highest pt(jj)",   "l");  
  leg1->AddEntry(H2_mjj_maxptmjj,   "highest pt(jj)/mjj",   "l");  
  leg1->AddEntry(H2_mjj_mgg,        "closest to mgg", "l");

  TLegend *leg2;
  leg2 = new TLegend(0.6,0.6,0.85,0.85);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(H2_mjj_gen,        "best MC truth match",   "l");  
  leg2->AddEntry(H2_mjj_btag_pt,    "btag, 2 highest pt jets",   "l");  
  leg2->AddEntry(H2_mjj_btag_ptjj,  "btag, highest pt(jj)",   "l");  
  leg2->AddEntry(H2_mjj_btag_ptmjj, "btag, highest pt(jj)/mjj",   "l");  
  leg2->AddEntry(H2_mjj_btag_mgg,   "btag, closest to mgg", "l");

  TLegend *leg3;
  leg3 = new TLegend(0.6,0.6,0.85,0.85);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(H1_dMmin,"min dM",   "l");  
  leg3->AddEntry(H1_dMmin_btag,"min dM, btag",   "l");  

  // plots
  TCanvas c1("c1","1 bjet, no btag",1);
  myHS->Draw();
  H1_mjj_gen        -> DrawNormalized("same");  
  H1_mjj_mgg        -> DrawNormalized("same");  
  H1_mjj_pt         -> DrawNormalized("same");
  H1_mjj_maxptjj    -> DrawNormalized("same"); 
  H1_mjj_maxptmjj   -> DrawNormalized("same"); 
  leg1->Draw();
  c1.SaveAs("H1_noBtag.root");

  TCanvas c11("c11","1 bjet, with btag",1);
  myHS->Draw();
  H1_mjj_gen        -> DrawNormalized("same");  
  H1_mjj_btag_mgg   -> DrawNormalized("same");  
  H1_mjj_btag_pt    -> DrawNormalized("same");  
  H1_mjj_btag_ptjj  -> DrawNormalized("same");  
  H1_mjj_btag_ptmjj -> DrawNormalized("same");  
  leg2->Draw();
  c11.SaveAs("H1_withBtag.root");

  TCanvas c2("c2","2 bjet, no btag",1);
  myHS->Draw();
  H2_mjj_gen        -> DrawNormalized("same");  
  H2_mjj_mgg        -> DrawNormalized("same");  
  H2_mjj_pt         -> DrawNormalized("same");
  H2_mjj_maxptjj    -> DrawNormalized("same"); 
  H2_mjj_maxptmjj   -> DrawNormalized("same"); 
  leg1->Draw();
  c2.SaveAs("H2_noBtag.root");

  TCanvas c22("c22","2 bjet, with btag",1);
  myHS->Draw();
  H2_mjj_gen        -> DrawNormalized("same");  
  H2_mjj_btag_mgg   -> DrawNormalized("same");  
  H2_mjj_btag_pt    -> DrawNormalized("same");  
  H2_mjj_btag_ptjj  -> DrawNormalized("same");  
  H2_mjj_btag_ptmjj -> DrawNormalized("same");  
  leg2->Draw();
  c22.SaveAs("H2_withBtag.root");

  /*
  TCanvas c11("c11","1 bjet",1);
  H1_dMmin      -> DrawNormalized();
  H1_dMmin_btag -> DrawNormalized("same");
  leg3->Draw();
  c11.SaveAs("dMmin1.root");

  TCanvas c21("c21","2 bjet",1);
  H2_dMmin      -> DrawNormalized();
  H2_dMmin_btag -> DrawNormalized("same");
  leg3->Draw();
  c21.SaveAs("dMmin2.root");
  */
}
