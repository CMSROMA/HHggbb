#include "bestJets_commonNtp.h"

using namespace std;

using namespace TMVA;

bestJets_commonNtp::bestJets_commonNtp( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer_commonNtp( "Radion", dataset ) {
  
  bTaggerType_ = bTaggerType;
  
  setSelectionType(selectionType);

  // jet regression                                                                                                   
  readerRegres = new Reader( "!Color:!Silent" );
  readerRegres->AddVariable( "hJet_pt",  &fRegr_pt);
  readerRegres->AddVariable( "hJet_eta", &fRegr_eta);
  readerRegres->AddVariable( "hJet_cef", &fRegr_cef);
  readerRegres->AddVariable( "hJet_nconstituents", &fRegr_nconst);
  readerRegres->AddVariable( "hJet_chf",    &fRegr_chf);
  readerRegres->AddVariable( "hJet_vtxPt",  &fRegr_vtxPt);
  readerRegres->AddVariable( "hJet_vtx3dL", &fRegr_vtx3dl);
  readerRegres->AddVariable( "MET", &fRegr_met);
  readerRegres->AddVariable( "hJet_dPhiMETJet", &fRegr_dPhiMet);
  readerRegres->BookMVA("BDTG method","data/regrWeights/TMVARegression_BDTG.weights.xml");
}

bestJets_commonNtp::~bestJets_commonNtp() {
  
  outFile_->Close();
  
  if (!tree_) return;
  delete tree_->GetCurrentFile();
}

double bestJets_commonNtp::delta_phi(double phi1, double phi2) {

  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}

void bestJets_commonNtp::finalize() {
  
  this->Init();
  
  std::string fullFlags = selectionType_ + "_";
  fullFlags+=bTaggerType_;
  this->set_flags(fullFlags); 
  this->createOutputFile();

  outFile_->cd();


  //-----------------------------------------------------------------------------
  // for the tree with selected events
  float invMassJJ_pt,      invMassJJ_maxptjj,   invMassJJ_maxptmjj,   invMassJJ_mgg; 
  float invMassJJ_btag_pt, invMassJJ_btag_ptjj, invMassJJ_btag_ptmjj, invMassJJ_btag_mgg, invMassJJ_btag_ptCSV;
  float invMassJJ_gen;
  int btagCategory_t;
  float dMmin_t, dMmin_btag_t;
  float minDrGenReco1_t, minDrGenReco2_t;
  int ngoodJets_t;
  float weight_t;

  TTree* myTrees = new TTree();
  myTrees->SetName("myTrees");
  myTrees->Branch( "mjj_pt",          &invMassJJ_pt,          "invMassJJ_pt/F" );
  myTrees->Branch( "mjj_maxptjj",     &invMassJJ_maxptjj,     "invMassJJ_maxptjj/F" );
  myTrees->Branch( "mjj_maxptmjj",    &invMassJJ_maxptmjj,    "invMassJJ_maxptmjj/F" );
  myTrees->Branch( "mjj_mgg",         &invMassJJ_mgg,         "invMassJJ_mgg/F" );  
  myTrees->Branch( "mjj_btag_pt",     &invMassJJ_btag_pt,     "invMassJJ_btag_pt/F" );  
  myTrees->Branch( "mjj_btag_ptCSV",  &invMassJJ_btag_ptCSV,  "invMassJJ_btag_ptCSV/F" );  
  myTrees->Branch( "mjj_btag_ptjj",   &invMassJJ_btag_ptjj,   "invMassJJ_btag_ptjj/F" );  
  myTrees->Branch( "mjj_btag_ptmjj",  &invMassJJ_btag_ptmjj,  "invMassJJ_btag_ptmjj/F" );  
  myTrees->Branch( "mjj_btag_mgg",    &invMassJJ_btag_mgg,    "invMassJJ_btag_mgg/F" );  
  myTrees->Branch( "mjj_gen",         &invMassJJ_gen,         "invMassJJ_gen/F" );
  myTrees->Branch( "dMmin",           &dMmin_t,               "dMmin_t/F" );    
  myTrees->Branch( "dMmin_btag",      &dMmin_btag_t,          "dMmin_btag_t/F" );    
  myTrees->Branch( "btagCategory",    &btagCategory_t,        "btagCategory_t/I" );
  myTrees->Branch( "ngoodJets",       &ngoodJets_t,           "ngoodJets_t/I" );
  myTrees->Branch( "ngoodJets",       &ngoodJets_t,           "ngoodJets_t/I" );
  myTrees->Branch( "minDrGenReco1",   &minDrGenReco1_t,       "minDrGenReco1_t/F" );
  myTrees->Branch( "minDrGenReco2",   &minDrGenReco2_t,       "minDrGenReco2_t/F" );


  // ------------------------------------------------------
  // analysis
  if (tree_ == 0) return;

  Long64_t nentries = tree_->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  // to check how often the choice is correct 
  float counter1 = 0.;
  float counter2 = 0.;
  //
  float counter_pt    = 0.;
  float counter_ptjj  = 0.;
  float counter_ptmjj = 0.;
  float counter_mgg   = 0.;
  // 
  float counter_btag_pt    = 0.;
  float counter_btag_ptjj  = 0.;
  float counter_btag_ptmjj = 0.;
  float counter_btag_mgg   = 0.;
  // 
  float counter_btagBasedOnPt  = 0.;
  float counter_btagBasedOnCSV = 0.;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = tree_->GetEntry(jentry);   nbytes += nb;
    
    // event weight for MC events: xsec, PU and smearings. It correspond to 19.62/fb , chiara
    weight_t = evweight;  


    // ---------------------------------------------------
    // gamma-gamma analysis

    // photons acceptance                                                                                                         
    if((TMath::Abs(ph1_SCEta)>1.4442&&TMath::Abs(ph1_SCEta)<1.566)||(TMath::Abs(ph2_SCEta)>1.4442&&TMath::Abs(ph2_SCEta)<1.566)
       || TMath::Abs(ph1_SCEta)>2.5 || TMath::Abs(ph2_SCEta)>2.5) continue;  // acceptance                                      

    // photon id                                                                                                                  
    bool idphot1(0), idphot2(0);
    idphot1 = (ph1_ciclevel >= cicselection);
    idphot2 = (ph2_ciclevel >= cicselection);
    if(cicselection>0) {
      if(!(idphot1)) continue;
      if(!(idphot2)) continue;
    }

    // extra cuts on photons: splitting events per photon class if needed                                                         
    int myEB=2;
    if( (fabs(ph1_SCEta)>1.4442 || fabs(ph2_SCEta)>1.4442) ) myEB=0;
    if( (fabs(ph1_SCEta)<1.4442 && fabs(ph2_SCEta)<1.4442) ) myEB=1;

    // the two selected photons for the analysis                                                                                  
    TLorentzVector t4phot1, t4phot2;
    t4phot1.SetPtEtaPhiM(ph1_pt,ph1_eta,ph1_phi,0.);
    t4phot2.SetPtEtaPhiM(ph2_pt,ph2_eta,ph2_phi,0.);
    TLorentzVector t4diPhot;
    t4diPhot.SetPtEtaPhiM( dipho_pt, dipho_eta, dipho_phi, PhotonsMass );
    TVector3 t3diPhot;
    t3diPhot.SetPtEtaPhi( dipho_pt, dipho_eta, dipho_phi );

    // photons pt cuts                                                                                                            
    if(ph1_pt<ptphot1cut * PhotonsMass/120.) continue;         // pt first photon                                                
    if(ph2_pt<ptphot2cut)                    continue;         // pt second photon     

    // invariant mass cut on photons                                                                                              
    if (PhotonsMass<100 || PhotonsMass>180) continue;


    // --------------------------------------------------------------------------------
    // preparing vectors with the infos used later on
    ecorrjet[0]     = j1_e;           ecorrjet[1]     = j2_e;           ecorrjet[2]     = j3_e;           ecorrjet[3]     = j4_e;
    ptcorrjet[0]    = j1_pt;          ptcorrjet[1]    = j2_pt;          ptcorrjet[2]    = j3_pt;          ptcorrjet[3]    = j4_pt;
    etajet[0]       = j1_eta;         etajet[1]       = j2_eta;         etajet[2]       = j3_eta;         etajet[3]       = j4_eta;
    phijet[0]       = j1_phi;         phijet[1]       = j2_phi;         phijet[2]       = j3_phi;         phijet[3]       = j4_phi;
    btagjprobjet[0] = j1_jetProbBtag; btagjprobjet[1] = j2_jetProbBtag; btagjprobjet[2] = j3_jetProbBtag; btagjprobjet[3] = j4_jetProbBtag;
    btagcsvjet[0]   = j1_csvBtag;     btagcsvjet[1]   = j2_csvBtag;     btagcsvjet[2]   = j3_csvBtag;     btagcsvjet[3]   = j4_csvBtag;
    cefjet[0]       = j1_emfrac;      cefjet[1]       = j2_emfrac;      cefjet[2]       = j3_emfrac;      cefjet[3]       = j4_emfrac;
    chfjet[0]       = j1_hadfrac;     chfjet[1]       = j2_hadfrac;     chfjet[2]       = j3_hadfrac;     chfjet[3]       = j4_hadfrac;
    vtxPtjet[0]     = j1_secVtxPt;    vtxPtjet[1]     = j2_secVtxPt;    vtxPtjet[2]     = j3_secVtxPt;    vtxPtjet[3]     = j4_secVtxPt;
    vtx3dljet[0]    = j1_secVtx3dL;   vtx3dljet[1]    = j2_secVtx3dL;   vtx3dljet[2]    = j3_secVtx3dL;   vtx3dljet[3]    = j4_secVtx3dL;
    nconstjet[0]    = (float)(j1_nNeutrals + j1_nCharged);   // + j1_ntk?? chiara!                                                                      
    nconstjet[1]    = (float)(j2_nNeutrals + j2_nCharged);
    nconstjet[2]    = (float)(j3_nNeutrals + j3_nCharged);
    nconstjet[3]    = (float)(j4_nNeutrals + j4_nCharged);

    // applying the jet regression
    TVector3 tempT3jet, t3met;
    // t3met.SetPtEtaPhi(met_corr_pfmet, 0, met_corr_phi_pfmet);   // chiara                                                                            
    t3met.SetPtEtaPhi(met_pfmet, 0, met_phi_pfmet);                // chiara  
    for (int ii=0; ii<4; ii++) {
      if ( ptcorrjet[ii]<-1) continue;

      fRegr_pt      = ptcorrjet[ii];
      fRegr_eta     = etajet[ii];
      fRegr_cef     = cefjet[ii];
      fRegr_nconst  = nconstjet[ii];
      fRegr_chf     = chfjet[ii];
      fRegr_vtxPt   = vtxPtjet[ii];
      fRegr_vtx3dl  = vtx3dljet[ii];
      // fRegr_met     = met_corr_pfmet;      // chiara                                                                                                 
      fRegr_met     = met_pfmet;              // chiara
      tempT3jet.SetPtEtaPhi(ptcorrjet[ii], etajet[ii], phijet[ii]);
      fRegr_dPhiMet = tempT3jet.DeltaPhi(t3met);
      float thePtCorr = readerRegres->EvaluateRegression("BDTG method")[0];

      // corrected pT and energy                                                                                                                        
      float correctionFactor = thePtCorr/ptcorrjet[ii];                                                                                              
      ptcorrjet[ii] = thePtCorr;                                                                                                                     
      ecorrjet[ii] *= correctionFactor;                                                                                                              
    }

    vector<int> v_puIdJets;
    for (int ij=0; ij<4; ij++) { 
      if ( ptcorrjet[ij]<-1)               continue;
      if ( btagcsvjet[ij]<=0)              continue;
      if ( ptcorrjet[ij] < ptjetacccut )   continue;
      if ( fabs(etajet[ij])>etajetacccut ) continue;
      if ( !passCutBasedJetId(ij) )        continue;
      v_puIdJets.push_back(ij);
    }

    // jets passing btagging ( + eta/pT/PUid cuts)                                                                                
    vector<int> v_looseCSV, v_mediumCSV, v_tightCSV;
    for (int ij=0; ij<int(v_puIdJets.size());ij++) {
      int index = v_puIdJets[ij];
      if (btagcsvjet[index]>0.244) v_looseCSV.push_back(index);
      if (btagcsvjet[index]>0.679) v_mediumCSV.push_back(index);
      if (btagcsvjet[index]>0.898) v_tightCSV.push_back(index);
    }

    // at least 2 preselected jets
    if (v_puIdJets.size()<2) continue;   
    std::pair<int,int> allJets_highPt = myTwoHighestPtJet(v_puIdJets);



    //  ----------------------------------------------------------------------
    // choice of analysis jets: no btag info used

    // 0) gen-level jets associated to the gen-level b
    TLorentzVector t4_gen1, t4_gen2;
    // t4_gen1.SetPtEtaPhiM(gr_j1_p4_pt, gr_j1_p4_eta, gr_j1_p4_phi, gr_j1_p4_mass);
    // t4_gen2.SetPtEtaPhiM(gr_j2_p4_pt, gr_j2_p4_eta, gr_j2_p4_phi, gr_j2_p4_mass);
    t4_gen1.SetPtEtaPhiM(gr_b1_p4_pt, gr_b1_p4_eta, gr_b1_p4_phi, gr_b1_p4_mass);
    t4_gen2.SetPtEtaPhiM(gr_b2_p4_pt, gr_b2_p4_eta, gr_b2_p4_phi, gr_b2_p4_mass);

    // reco jets closest to the above gen jets
    int jet1_genJ = -1;
    int jet2_genJ = -1;
    float minDrGenJ1 = 999.;
    float minDrGenJ2 = 999.;
    for (int jet=0; jet<(v_puIdJets.size()); jet++) {
      int index = v_puIdJets[jet];
      TLorentzVector t4jet;
      t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
      float theDr1 = t4_gen1.DeltaR(t4jet);
      float theDr2 = t4_gen2.DeltaR(t4jet);
      if ( theDr1<theDr2 && theDr1<minDrGenJ1 ) {
	minDrGenJ1 = theDr1;  
	jet1_genJ  = index;
      }
      if ( theDr2<theDr1 && theDr2<minDrGenJ2 ) {
	minDrGenJ2 = theDr2;  
	jet2_genJ  = index;
      }
    }      
    if (jet1_genJ==jet2_genJ) cout << "this can not happen" << endl;
    
    // if there is no good match with gen jet skip the event
    if (minDrGenJ1>0.07) continue;   
    if (minDrGenJ2>0.07) continue;   

    TLorentzVector t4jet1_genJ, t4jet2_genJ;
    t4jet1_genJ.SetPtEtaPhiE(ptcorrjet[jet1_genJ],etajet[jet1_genJ],phijet[jet1_genJ],ecorrjet[jet1_genJ]);
    t4jet2_genJ.SetPtEtaPhiE(ptcorrjet[jet2_genJ],etajet[jet2_genJ],phijet[jet2_genJ],ecorrjet[jet2_genJ]);
    TLorentzVector t4diJet_genJ = t4jet1_genJ + t4jet2_genJ;


    // -----------------------------------------
    // 1) highest pT jets
    int jet1_pt = allJets_highPt.first;
    int jet2_pt = allJets_highPt.second;    

    TLorentzVector t4jet1_pt, t4jet2_pt;
    t4jet1_pt.SetPtEtaPhiE(ptcorrjet[jet1_pt],etajet[jet1_pt],phijet[jet1_pt],ecorrjet[jet1_pt]);
    t4jet2_pt.SetPtEtaPhiE(ptcorrjet[jet2_pt],etajet[jet2_pt],phijet[jet2_pt],ecorrjet[jet2_pt]);
    TLorentzVector t4diJet_pt = t4jet1_pt + t4jet2_pt;


    // -----------------------------------------
    // 2) jets giving the highest pT(jj)
    int jet1_maxptjj = -1;
    int jet2_maxptjj = -1;
    float maxPt = -999.;
    for (int jetA=0; jetA<(v_puIdJets.size()-1); jetA++) {
      for (int jetB=jetA+1; jetB<v_puIdJets.size(); jetB++) {
	TLorentzVector t4jetA, t4jetB;
	int indexA = v_puIdJets[jetA];
	int indexB = v_puIdJets[jetB];
	t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	TLorentzVector t4AB = t4jetA + t4jetB;
	if ( t4AB.Pt() > maxPt) {
	  maxPt = t4AB.Pt();
	  jet1_maxptjj = indexA;
	  jet2_maxptjj = indexB;
	}
      }
    }

    TLorentzVector t4jet1_maxptjj, t4jet2_maxptjj;
    t4jet1_maxptjj.SetPtEtaPhiE(ptcorrjet[jet1_maxptjj],etajet[jet1_maxptjj],phijet[jet1_maxptjj],ecorrjet[jet1_maxptjj]);
    t4jet2_maxptjj.SetPtEtaPhiE(ptcorrjet[jet2_maxptjj],etajet[jet2_maxptjj],phijet[jet2_maxptjj],ecorrjet[jet2_maxptjj]);
    TLorentzVector t4diJet_maxptjj = t4jet1_maxptjj + t4jet2_maxptjj;


    // -----------------------------------------
    // 3) jets giving the highest pT(jj) / m(jj)
    int jet1_maxptmjj = -1;
    int jet2_maxptmjj = -1;
    float maxPtmjj = -999.;
    for (int jetA=0; jetA<(v_puIdJets.size()-1); jetA++) {
      for (int jetB=jetA+1; jetB<v_puIdJets.size(); jetB++) {
	TLorentzVector t4jetA, t4jetB;
	int indexA = v_puIdJets[jetA];
	int indexB = v_puIdJets[jetB];
	t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	TLorentzVector t4AB = t4jetA + t4jetB;
	if ( t4AB.Pt()/t4AB.M() > maxPtmjj) {
	  maxPtmjj = t4AB.Pt()/t4AB.M();
	  jet1_maxptmjj = indexA;
	  jet2_maxptmjj = indexB;
	}
      }
    }

    TLorentzVector t4jet1_maxptmjj, t4jet2_maxptmjj;
    t4jet1_maxptmjj.SetPtEtaPhiE(ptcorrjet[jet1_maxptmjj],etajet[jet1_maxptmjj],phijet[jet1_maxptmjj],ecorrjet[jet1_maxptmjj]);
    t4jet2_maxptmjj.SetPtEtaPhiE(ptcorrjet[jet2_maxptmjj],etajet[jet2_maxptmjj],phijet[jet2_maxptmjj],ecorrjet[jet2_maxptmjj]);
    TLorentzVector t4diJet_maxptmjj = t4jet1_maxptmjj + t4jet2_maxptmjj;


    // -----------------------------------------
    // 4) jets giving mjj closest to mgg
    int jet1_mgg  = -1;
    int jet2_mgg  = -1;
    float minDMgg = 999.;
    for (int jetA=0; jetA<(v_puIdJets.size()-1); jetA++) {
      for (int jetB=jetA+1; jetB<v_puIdJets.size(); jetB++) {
	TLorentzVector t4jetA, t4jetB;
	int indexA = v_puIdJets[jetA];
	int indexB = v_puIdJets[jetB];
	t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	TLorentzVector t4AB = t4jetA + t4jetB;
	float massJJ = t4AB.M();
	float massGG = t4diPhot.M();
	if ( fabs(massJJ-massGG) < minDMgg ) {
	  minDMgg  = fabs(massJJ-massGG); 
	  jet1_mgg = indexA;
	  jet2_mgg = indexB;
	}
      }
    }
    
    TLorentzVector t4jet1_mgg, t4jet2_mgg;
    t4jet1_mgg.SetPtEtaPhiE(ptcorrjet[jet1_mgg],etajet[jet1_mgg],phijet[jet1_mgg],ecorrjet[jet1_mgg]);
    t4jet2_mgg.SetPtEtaPhiE(ptcorrjet[jet2_mgg],etajet[jet2_mgg],phijet[jet2_mgg],ecorrjet[jet2_mgg]);
    TLorentzVector t4diJet_mgg = t4jet1_mgg + t4jet2_mgg;


    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    // choose jets looking to btag: highest pT ones
    int jet1btag_pt = -1;
    int jet2btag_pt = -1;
    if (v_mediumCSV.size()==1) { 
      jet1btag_pt = myHighestPtJet(v_mediumCSV);      
    } else if (v_mediumCSV.size()>1) {
      std::pair<int,int> mediumBtag_highPt = myTwoHighestPtJet(v_mediumCSV);
      jet1btag_pt = mediumBtag_highPt.first;
      jet2btag_pt = mediumBtag_highPt.second;
    }

    // choose jets looking to btag: highest CSV output
    int jet1btag_CSV = -1;
    int jet2btag_CSV = -1;
    if (v_mediumCSV.size()==1) { 
      jet1btag_CSV = myHighestBtagJet(v_mediumCSV);                          
    } else if (v_mediumCSV.size()>1) {
      std::pair<int,int> mediumBtag_csv = myTwoHighestBtagJet(v_mediumCSV);  
      jet1btag_CSV = mediumBtag_csv.first;
      jet2btag_CSV = mediumBtag_csv.second;
    }

    // -----------------------------------------

    // 1a) now apply the same criteria as above but giving priorities to the btagged jets
    //     bjets: pT ordered ; nobjets (btagCat=1): pT ordered
    int jet1_btag_pt = -1;
    int jet2_btag_pt = -1;
    if( v_mediumCSV.size()==1 ) {
      if (jet1_pt != jet1btag_pt && jet2_pt != jet1btag_pt) {
	jet2_btag_pt = jet1_pt;
	jet1_btag_pt = jet1btag_pt;
      } else if (jet1_pt == jet1btag_pt) {
	jet1_btag_pt = jet1btag_pt;
	jet2_btag_pt = jet2_pt;
      } else if (jet2_pt == jet1btag_pt) {
	jet1_btag_pt = jet1btag_pt;
	jet2_btag_pt = jet1_pt;
      }
    }
    if( v_mediumCSV.size()>1 ) {
      jet1_btag_pt = jet1btag_pt;
      jet2_btag_pt = jet2btag_pt;
    }

    TLorentzVector t4jet1_btag_pt, t4jet2_btag_pt;
    t4jet1_btag_pt.SetPtEtaPhiE(ptcorrjet[jet1_btag_pt],etajet[jet1_btag_pt],phijet[jet1_btag_pt],ecorrjet[jet1_btag_pt]);
    t4jet2_btag_pt.SetPtEtaPhiE(ptcorrjet[jet2_btag_pt],etajet[jet2_btag_pt],phijet[jet2_btag_pt],ecorrjet[jet2_btag_pt]);
    TLorentzVector t4diJet_btag_pt = t4jet1_btag_pt + t4jet2_btag_pt;

    // -----------------------------------------

    // 1b) now apply the same criteria as above but giving priorities to the btagged jets
    //     bjets: CSV ordered ; nobjets (btagCat=1): pT ordered
    int jet1_btag_ptCSV = -1;
    int jet2_btag_ptCSV = -1;
    if( v_mediumCSV.size()==1 ) {
      if (jet1_pt != jet1btag_CSV && jet2_pt != jet1btag_CSV) {
	jet2_btag_ptCSV = jet1_pt;
	jet1_btag_ptCSV = jet1btag_CSV;
      } else if (jet1_pt == jet1btag_CSV) {
	jet1_btag_ptCSV = jet1btag_CSV;
	jet2_btag_ptCSV = jet2_pt;
      } else if (jet2_pt == jet1btag_CSV) {
	jet1_btag_ptCSV = jet1btag_CSV;
	jet2_btag_ptCSV = jet1_pt;
      }
    }
    if( v_mediumCSV.size()>1 ) {
      jet1_btag_ptCSV = jet1btag_CSV;
      jet2_btag_ptCSV = jet2btag_CSV;
    }

    TLorentzVector t4jet1_btag_ptCSV, t4jet2_btag_ptCSV;
    t4jet1_btag_ptCSV.SetPtEtaPhiE(ptcorrjet[jet1_btag_ptCSV],etajet[jet1_btag_ptCSV],phijet[jet1_btag_ptCSV],ecorrjet[jet1_btag_ptCSV]);
    t4jet2_btag_ptCSV.SetPtEtaPhiE(ptcorrjet[jet2_btag_ptCSV],etajet[jet2_btag_ptCSV],phijet[jet2_btag_ptCSV],ecorrjet[jet2_btag_ptCSV]);
    TLorentzVector t4diJet_btag_ptCSV = t4jet1_btag_ptCSV + t4jet2_btag_ptCSV;


    // -----------------------------------------
    // 2a) now apply the same criteria as above but giving priorities to the btagged jets: 
    //     pair giving the highest pT(jj)/mjj, of which at least 1 is btag
    int jet1_btag_ptjj = -1;
    int jet2_btag_ptjj = -1;
    float maxPt2 = -999.;
    if( v_mediumCSV.size()==1 ) {
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
	int index = v_puIdJets[jet];
	if (index==jet1btag_pt) continue;
	TLorentzVector t4jet, t4bjet;
	t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
	t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag_pt],etajet[jet1btag_pt],phijet[jet1btag_pt],ecorrjet[jet1btag_pt]);      
	TLorentzVector t4_bnotb = t4jet + t4bjet; 
	if ( t4_bnotb.Pt() > maxPt2) {
	  maxPt2 = t4_bnotb.Pt();
	  jet1_btag_ptjj = jet1btag_pt;
	  jet2_btag_ptjj = index;
	}
      }
    }
    else if( v_mediumCSV.size()>1 ) {
      for (int jetA=0; jetA<(v_mediumCSV.size()-1); jetA++) { 
	for (int jetB=jetA+1; jetB<v_mediumCSV.size(); jetB++) { 
	  TLorentzVector t4jetA, t4jetB;
	  int indexA = v_mediumCSV[jetA];
	  int indexB = v_mediumCSV[jetB];
	  t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	  t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	  TLorentzVector t4AB = t4jetA + t4jetB;
	  if ( t4AB.Pt() > maxPt2) {
	    maxPt2 = t4AB.Pt();
	    jet1_btag_ptjj = indexA;
	    jet2_btag_ptjj = indexB;
	  }
	}
      }
    }

    TLorentzVector t4jet1_btag_ptjj, t4jet2_btag_ptjj;
    t4jet1_btag_ptjj.SetPtEtaPhiE(ptcorrjet[jet1_btag_ptjj],etajet[jet1_btag_ptjj],phijet[jet1_btag_ptjj],ecorrjet[jet1_btag_ptjj]);
    t4jet2_btag_ptjj.SetPtEtaPhiE(ptcorrjet[jet2_btag_ptjj],etajet[jet2_btag_ptjj],phijet[jet2_btag_ptjj],ecorrjet[jet2_btag_ptjj]);
    TLorentzVector t4diJet_btag_ptjj = t4jet1_btag_ptjj + t4jet2_btag_ptjj;

    // -----------------------------------------

    // 3a) now apply the same criteria as above but giving priorities to the btagged jets: 
    //     pair giving the highest pT(jj)/m(jj), of which at least 1 is btag
    int jet1_btag_ptmjj = -1;
    int jet2_btag_ptmjj = -1;
    float maxPtMjj2 = -999.;
    if( v_mediumCSV.size()==1 ) {
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
	int index = v_puIdJets[jet];
	if (index==jet1btag_pt) continue;
	TLorentzVector t4jet, t4bjet;
	t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
	t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag_pt],etajet[jet1btag_pt],phijet[jet1btag_pt],ecorrjet[jet1btag_pt]);      
	TLorentzVector t4_bnotb = t4jet + t4bjet; 
	if ( (t4_bnotb.Pt()/t4_bnotb.M()) > maxPtMjj2) {
	  maxPtMjj2 = t4_bnotb.Pt()/t4_bnotb.M();
	  jet1_btag_ptmjj = jet1btag_pt;
	  jet2_btag_ptmjj = index;
	}
      }
    }
    else if( v_mediumCSV.size()>1 ) {
      for (int jetA=0; jetA<(v_mediumCSV.size()-1); jetA++) { 
	for (int jetB=jetA+1; jetB<v_mediumCSV.size(); jetB++) { 
	  TLorentzVector t4jetA, t4jetB;
	  int indexA = v_mediumCSV[jetA];
	  int indexB = v_mediumCSV[jetB];
	  t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	  t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	  TLorentzVector t4AB = t4jetA + t4jetB;
	  if ( (t4AB.Pt()/t4AB.M()) > maxPtMjj2) {
	    maxPtMjj2 = t4AB.Pt()/t4AB.M();
	    jet1_btag_ptmjj = indexA;
	    jet2_btag_ptmjj = indexB;
	  }
	}
      }
    }

    TLorentzVector t4jet1_btag_ptmjj, t4jet2_btag_ptmjj;
    t4jet1_btag_ptmjj.SetPtEtaPhiE(ptcorrjet[jet1_btag_ptmjj],etajet[jet1_btag_ptmjj],phijet[jet1_btag_ptmjj],ecorrjet[jet1_btag_ptmjj]);
    t4jet2_btag_ptmjj.SetPtEtaPhiE(ptcorrjet[jet2_btag_ptmjj],etajet[jet2_btag_ptmjj],phijet[jet2_btag_ptmjj],ecorrjet[jet2_btag_ptmjj]);
    TLorentzVector t4diJet_btag_ptmjj = t4jet1_btag_ptmjj + t4jet2_btag_ptmjj;

    // -----------------------------------------

    // 4a) now apply the same criteria as above but giving priorities to the btagged jets: 
    //     pair giving the closest mgg - m(jj), of which at least 1 is btag
    int jet1_btag_mgg = -1;
    int jet2_btag_mgg = -1;
    float minDMgg2    = 999.;
    if( v_mediumCSV.size()==1 ) {
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
	int index = v_puIdJets[jet];
	if (index==jet1btag_pt) continue;
	TLorentzVector t4jet, t4bjet;
	t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
	t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag_pt],etajet[jet1btag_pt],phijet[jet1btag_pt],ecorrjet[jet1btag_pt]);      
	TLorentzVector t4_bnotb = t4jet + t4bjet; 
	float massJJ = t4_bnotb.M();
	float massGG = t4diPhot.M();
	if ( fabs(massJJ-massGG) < minDMgg2 ) {
	  minDMgg2 = fabs(massJJ-massGG);
	  jet1_btag_mgg = jet1btag_pt;
	  jet2_btag_mgg = index;
	}
      }
    } else if (v_mediumCSV.size()>1 ) {
      for (int jetA=0; jetA<(v_mediumCSV.size()-1); jetA++) { 
	for (int jetB=jetA+1; jetB<v_mediumCSV.size(); jetB++) { 
	  TLorentzVector t4jetA, t4jetB;
	  int indexA = v_mediumCSV[jetA];
	  int indexB = v_mediumCSV[jetB];
	  t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	  t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	  TLorentzVector t4AB = t4jetA + t4jetB;
	  float massJJ = t4AB.M();
	  float massGG = t4diPhot.M();
	  if ( fabs(massJJ-massGG) < minDMgg2 ) {
	    minDMgg2 = fabs(massJJ-massGG);
	    jet1_btag_mgg = indexA;
	    jet2_btag_mgg = indexB;
	  }
	}
      }
    }
    
    TLorentzVector t4jet1_btag_mgg, t4jet2_btag_mgg;
    t4jet1_btag_mgg.SetPtEtaPhiE(ptcorrjet[jet1_btag_mgg],etajet[jet1_btag_mgg],phijet[jet1_btag_mgg],ecorrjet[jet1_btag_mgg]);
    t4jet2_btag_mgg.SetPtEtaPhiE(ptcorrjet[jet2_btag_mgg],etajet[jet2_btag_mgg],phijet[jet2_btag_mgg],ecorrjet[jet2_btag_mgg]);
    TLorentzVector t4diJet_btag_mgg = t4jet1_btag_mgg + t4jet2_btag_mgg;

    
    // several invariant masses combinations
    invMassJJ_pt         = t4diJet_pt.M();
    invMassJJ_maxptjj    = t4diJet_maxptjj.M();
    invMassJJ_maxptmjj   = t4diJet_maxptmjj.M();
    if ( minDMgg<999 ) 
      invMassJJ_mgg      = t4diJet_mgg.M();
    else 
      invMassJJ_mgg      = -1.; 
    invMassJJ_btag_pt    = t4diJet_btag_pt.M();
    invMassJJ_btag_ptCSV = t4diJet_btag_ptCSV.M();
    invMassJJ_btag_ptjj  = t4diJet_btag_ptjj.M();
    invMassJJ_btag_ptmjj = t4diJet_btag_ptmjj.M();
    if ( minDMgg2<999 )
      invMassJJ_btag_mgg   = t4diJet_btag_mgg.M();
    else
      invMassJJ_btag_mgg = -1.;
    invMassJJ_gen        = t4diJet_genJ.M();


    // counters: best criterium with/wo btag 
    if ( jet1_genJ>=0 && jet2_genJ>=0 &&
	 jet1_pt>=0 && jet2_pt>=0 && 
	 jet1_maxptjj>=0 && jet2_maxptjj>=0 && 
	 jet1_maxptmjj>=0 && jet2_maxptmjj>=0 && 
	 jet1_mgg>=0 && jet2_mgg>=0 &&
	 jet1_btag_pt>=0 && jet2_btag_pt>=0 && 
	 jet1_btag_ptjj>=0 && jet2_btag_ptjj>=0 && 
	 jet1_btag_ptmjj>=0 && jet2_btag_ptmjj>=0 && 
	 jet1_btag_mgg>=0 && jet2_btag_mgg>=0) {

      counter1++; 
      if ( (jet1_pt==jet1_genJ       && jet2_pt==jet2_genJ)       || (jet1_pt==jet2_genJ       && jet2_pt==jet1_genJ) )       counter_pt++;
      if ( (jet1_maxptjj==jet1_genJ  && jet2_maxptjj==jet2_genJ)  || (jet1_maxptjj==jet2_genJ  && jet2_maxptjj==jet1_genJ) )  counter_ptjj++;
      if ( (jet1_maxptmjj==jet1_genJ && jet2_maxptmjj==jet2_genJ) || (jet1_maxptmjj==jet2_genJ && jet2_maxptmjj==jet1_genJ) ) counter_ptmjj++;
      if ( (jet1_mgg==jet1_genJ      && jet2_mgg==jet2_genJ)      || (jet1_mgg==jet2_genJ      && jet2_mgg==jet1_genJ) )      counter_mgg++;
      if ( (jet1_btag_pt==jet1_genJ    && jet2_btag_pt==jet2_genJ)    || (jet1_btag_pt==jet2_genJ    && jet2_btag_pt==jet1_genJ) )    counter_btag_pt++;
      if ( (jet1_btag_ptjj==jet1_genJ  && jet2_btag_ptjj==jet2_genJ)  || (jet1_btag_ptjj==jet2_genJ  && jet2_btag_ptjj==jet1_genJ) )  counter_btag_ptjj++;
      if ( (jet1_btag_ptmjj==jet1_genJ && jet2_btag_ptmjj==jet2_genJ) || (jet1_btag_ptmjj==jet2_genJ && jet2_btag_ptmjj==jet1_genJ) ) counter_btag_ptmjj++;
      if ( (jet1_btag_mgg==jet1_genJ   && jet2_btag_mgg==jet2_genJ)   || (jet1_btag_mgg==jet2_genJ   && jet2_btag_mgg==jet1_genJ) )   counter_btag_mgg++;
    }

    // best criterium for btagging: pt vs CSV output
    if (jet1_genJ>0 && jet2_genJ>0 && jet1btag_CSV>0 && jet2btag_CSV>0 && jet1btag_pt>0 && jet2btag_pt>0) {
      counter2++;
      if ( (jet1btag_CSV==jet1_genJ && jet2btag_CSV==jet2_genJ) || (jet1btag_CSV==jet2_genJ && jet2btag_CSV==jet1_genJ) ) counter_btagBasedOnCSV++;
      if ( (jet1btag_pt==jet1_genJ  && jet2btag_pt==jet2_genJ)  || (jet1btag_pt==jet2_genJ  && jet2btag_pt==jet1_genJ) )  counter_btagBasedOnPt++;
    }


    // ---------------------------------------------------------------
    // counting the number of bjets to categorize the events
    int btagCategory = -1;
    if (bTaggerType_=="CSV") { btagCategory = (v_mediumCSV.size()<=2) ? v_mediumCSV.size() : 2; }
    else cout << "this btag algo does not exist" << endl;


    // --------------------------------------------------------------
    // here we apply further cuts according to the btag category
    if (btagCategory<1 || btagCategory>2) continue;


    // filling the tree for selected events 
    btagCategory_t  = btagCategory;
    ngoodJets_t     = v_puIdJets.size(); 
    dMmin_t         = minDMgg;
    dMmin_btag_t    = minDMgg2;
    minDrGenReco1_t = minDrGenJ1;
    minDrGenReco2_t = minDrGenJ2;

    myTrees->Fill();

  } // loop over entries 


  // summary:
  cout << "max pt: "    << counter_pt/counter1    << endl; 
  cout << "max ptjj: "  << counter_ptjj/counter1  << endl; 
  cout << "max ptmjj: " << counter_ptmjj/counter1 << endl; 
  cout << "min dm gg: " << counter_mgg/counter1   << endl; 
  cout << "btag, max pt: "    << counter_btag_pt/counter1    << endl; 
  cout << "btag, max ptjj: "  << counter_btag_ptjj/counter1  << endl; 
  cout << "btag, max ptmjj: " << counter_btag_ptmjj/counter1 << endl; 
  cout << "btag, min dm gg: " << counter_btag_mgg/counter1   << endl; 
  cout << endl;
  cout << "best btag choice:" << endl;
  cout << "CSV: " << counter_btagBasedOnCSV/counter2 << endl;
  cout << "Pt: "  << counter_btagBasedOnPt/counter2  << endl;

  outFile_->cd();
  myTrees->Write();

} // finalize

Int_t bestJets_commonNtp::GetEntry(Long64_t entry) {

  if (!tree_) return 0;
  return tree_->GetEntry(entry);
}

Long64_t bestJets_commonNtp::LoadTree(Long64_t entry) {

  if (!tree_) return -5;
  Long64_t centry = tree_->LoadTree(entry);
  if (centry < 0) return centry;
  if (tree_->GetTreeNumber() != fCurrent) {
    fCurrent = tree_->GetTreeNumber();
  }
  return centry;
}

void bestJets_commonNtp::Init() {

  // The Init() function is called when the selector needs to initialize                                                                      
  // a new tree or chain. Typically here the branch addresses and branch                                                                      
  // pointers of the tree will be set.                                                                                                        
  // It is normally not necessary to make changes to the generated                                                                            
  // code, but the routine can be extended by the user if needed.                                                                             
  // Init() will be called many times when running on PROOF                                                                                   
  // (once per file to be processed).                                                                                                         

  // Set branch addresses and branch pointers                                                                          
  fCurrent = -1;
  tree_->SetMakeClass(1);

  tree_->SetBranchAddress("itype", &itype, &b_itype);
  tree_->SetBranchAddress("run", &run, &b_run);
  tree_->SetBranchAddress("lumis", &lumis, &b_lumis);
  tree_->SetBranchAddress("event", &event, &b_event);
  tree_->SetBranchAddress("weight", &weight, &b_weight);
  tree_->SetBranchAddress("evweight", &evweight, &b_evweight);
  tree_->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  tree_->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
  tree_->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree_->SetBranchAddress("rho", &rho, &b_rho);
  tree_->SetBranchAddress("met_pfmet", &met_pfmet, &b_met_pfmet);
  tree_->SetBranchAddress("met_phi_pfmet", &met_phi_pfmet, &b_met_phi_pfmet);
  tree_->SetBranchAddress("met_corr_pfmet", &met_corr_pfmet, &b_met_corr_pfmet);
  tree_->SetBranchAddress("met_corr_phi_pfmet", &met_corr_phi_pfmet, &b_met_corr_phi_pfmet);
  tree_->SetBranchAddress("met_corr_eta_pfmet", &met_corr_eta_pfmet, &b_met_corr_eta_pfmet);
  tree_->SetBranchAddress("met_corr_e_pfmet", &met_corr_e_pfmet, &b_met_corr_e_pfmet);
  tree_->SetBranchAddress("category", &category, &b_category);
  tree_->SetBranchAddress("ph1_e", &ph1_e, &b_ph1_e);
  tree_->SetBranchAddress("ph2_e", &ph2_e, &b_ph2_e);
  tree_->SetBranchAddress("ph1_pt", &ph1_pt, &b_ph1_pt);
  tree_->SetBranchAddress("ph2_pt", &ph2_pt, &b_ph2_pt);
  tree_->SetBranchAddress("ph1_phi", &ph1_phi, &b_ph1_phi);
  tree_->SetBranchAddress("ph2_phi", &ph2_phi, &b_ph2_phi);
  tree_->SetBranchAddress("ph1_eta", &ph1_eta, &b_ph1_eta);
  tree_->SetBranchAddress("ph2_eta", &ph2_eta, &b_ph2_eta);
  tree_->SetBranchAddress("ph1_r9", &ph1_r9, &b_ph1_r9);
  tree_->SetBranchAddress("ph2_r9", &ph2_r9, &b_ph2_r9);
  tree_->SetBranchAddress("ph1_isPrompt", &ph1_isPrompt, &b_ph1_isPrompt);
  tree_->SetBranchAddress("ph2_isPrompt", &ph2_isPrompt, &b_ph2_isPrompt);
  tree_->SetBranchAddress("ph1_SCEta", &ph1_SCEta, &b_ph1_SCEta);
  tree_->SetBranchAddress("ph2_SCEta", &ph2_SCEta, &b_ph2_SCEta);
  tree_->SetBranchAddress("ph1_SCPhi", &ph1_SCPhi, &b_ph1_SCPhi);
  tree_->SetBranchAddress("ph2_SCPhi", &ph2_SCPhi, &b_ph2_SCPhi);
  tree_->SetBranchAddress("ph1_hoe", &ph1_hoe, &b_ph1_hoe);
  tree_->SetBranchAddress("ph2_hoe", &ph2_hoe, &b_ph2_hoe);
  tree_->SetBranchAddress("ph1_sieie", &ph1_sieie, &b_ph1_sieie);
  tree_->SetBranchAddress("ph2_sieie", &ph2_sieie, &b_ph2_sieie);
  tree_->SetBranchAddress("ph1_pfchargedisogood03", &ph1_pfchargedisogood03, &b_ph1_pfchargedisogood03);
  tree_->SetBranchAddress("ph2_pfchargedisogood03", &ph2_pfchargedisogood03, &b_ph2_pfchargedisogood03);
  tree_->SetBranchAddress("ph1_pfchargedisobad04", &ph1_pfchargedisobad04, &b_ph1_pfchargedisobad04);
  tree_->SetBranchAddress("ph2_pfchargedisobad04", &ph2_pfchargedisobad04, &b_ph2_pfchargedisobad04);
  tree_->SetBranchAddress("ph1_etawidth", &ph1_etawidth, &b_ph1_etawidth);
  tree_->SetBranchAddress("ph2_etawidth", &ph2_etawidth, &b_ph2_etawidth);
  tree_->SetBranchAddress("ph1_phiwidth", &ph1_phiwidth, &b_ph1_phiwidth);
  tree_->SetBranchAddress("ph2_phiwidth", &ph2_phiwidth, &b_ph2_phiwidth);
  tree_->SetBranchAddress("ph1_eseffssqrt", &ph1_eseffssqrt, &b_ph1_eseffssqrt);
  tree_->SetBranchAddress("ph2_eseffssqrt", &ph2_eseffssqrt, &b_ph2_eseffssqrt);
  tree_->SetBranchAddress("ph1_pfchargedisobad03", &ph1_pfchargedisobad03, &b_ph1_pfchargedisobad03);
  tree_->SetBranchAddress("ph2_pfchargedisobad03", &ph2_pfchargedisobad03, &b_ph2_pfchargedisobad03);
  tree_->SetBranchAddress("ph1_sieip", &ph1_sieip, &b_ph1_sieip);
  tree_->SetBranchAddress("ph2_sieip", &ph2_sieip, &b_ph2_sieip);
  tree_->SetBranchAddress("ph1_sipip", &ph1_sipip, &b_ph1_sipip);
  tree_->SetBranchAddress("ph2_sipip", &ph2_sipip, &b_ph2_sipip);
  tree_->SetBranchAddress("ph1_ecaliso", &ph1_ecaliso, &b_ph1_ecaliso);
  tree_->SetBranchAddress("ph2_ecaliso", &ph2_ecaliso, &b_ph2_ecaliso);
  tree_->SetBranchAddress("ph1_ecalisobad", &ph1_ecalisobad, &b_ph1_ecalisobad);
  tree_->SetBranchAddress("ph2_ecalisobad", &ph2_ecalisobad, &b_ph2_ecalisobad);
  tree_->SetBranchAddress("ph1_badvtx_Et", &ph1_badvtx_Et, &b_ph1_badvtx_Et);
  tree_->SetBranchAddress("ph2_badvtx_Et", &ph2_badvtx_Et, &b_ph2_badvtx_Et);
  tree_->SetBranchAddress("ph1_isconv", &ph1_isconv, &b_ph1_isconv);
  tree_->SetBranchAddress("ph2_isconv", &ph2_isconv, &b_ph2_isconv);
  tree_->SetBranchAddress("ph1_ciclevel", &ph1_ciclevel, &b_ph1_ciclevel);
  tree_->SetBranchAddress("ph2_ciclevel", &ph2_ciclevel, &b_ph2_ciclevel);
  tree_->SetBranchAddress("ph1_sigmaEoE", &ph1_sigmaEoE, &b_ph1_sigmaEoE);
  tree_->SetBranchAddress("ph2_sigmaEoE", &ph2_sigmaEoE, &b_ph2_sigmaEoE);
  tree_->SetBranchAddress("ph1_ptoM", &ph1_ptoM, &b_ph1_ptoM);
  tree_->SetBranchAddress("ph2_ptoM", &ph2_ptoM, &b_ph2_ptoM);
  tree_->SetBranchAddress("ph1_isEB", &ph1_isEB, &b_ph1_isEB);
  tree_->SetBranchAddress("ph2_isEB", &ph2_isEB, &b_ph2_isEB);
  tree_->SetBranchAddress("ph1_s4ratio", &ph1_s4ratio, &b_ph1_s4ratio);
  tree_->SetBranchAddress("ph2_s4ratio", &ph2_s4ratio, &b_ph2_s4ratio);
  tree_->SetBranchAddress("ph1_e3x3", &ph1_e3x3, &b_ph1_e3x3);
  tree_->SetBranchAddress("ph2_e3x3", &ph2_e3x3, &b_ph2_e3x3);
  tree_->SetBranchAddress("ph1_e5x5", &ph1_e5x5, &b_ph1_e5x5);
  tree_->SetBranchAddress("ph2_e5x5", &ph2_e5x5, &b_ph2_e5x5);
  tree_->SetBranchAddress("PhotonsMass", &PhotonsMass, &b_PhotonsMass);
  tree_->SetBranchAddress("dipho_E", &dipho_E, &b_dipho_E);
  tree_->SetBranchAddress("dipho_pt", &dipho_pt, &b_dipho_pt);
  tree_->SetBranchAddress("dipho_eta", &dipho_eta, &b_dipho_eta);
  tree_->SetBranchAddress("dipho_phi", &dipho_phi, &b_dipho_phi);
  tree_->SetBranchAddress("dipho_cosThetaStar_CS", &dipho_cosThetaStar_CS, &b_dipho_cosThetaStar_CS);
  tree_->SetBranchAddress("dipho_tanhYStar", &dipho_tanhYStar, &b_dipho_tanhYStar);
  tree_->SetBranchAddress("dipho_Y", &dipho_Y, &b_dipho_Y);
  tree_->SetBranchAddress("vtx_ind", &vtx_ind, &b_vtx_ind);
  tree_->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
  tree_->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
  tree_->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
  tree_->SetBranchAddress("vtx_mva", &vtx_mva, &b_vtx_mva);
  tree_->SetBranchAddress("vtx_mva_2", &vtx_mva_2, &b_vtx_mva_2);
  tree_->SetBranchAddress("vtx_mva_3", &vtx_mva_3, &b_vtx_mva_3);
  tree_->SetBranchAddress("vtx_ptbal", &vtx_ptbal, &b_vtx_ptbal);
  tree_->SetBranchAddress("vtx_ptasym", &vtx_ptasym, &b_vtx_ptasym);
  tree_->SetBranchAddress("vtx_logsumpt2", &vtx_logsumpt2, &b_vtx_logsumpt2);
  tree_->SetBranchAddress("vtx_pulltoconv", &vtx_pulltoconv, &b_vtx_pulltoconv);
  tree_->SetBranchAddress("vtx_prob", &vtx_prob, &b_vtx_prob);
  tree_->SetBranchAddress("njets_passing_kLooseID", &njets_passing_kLooseID, &b_njets_passing_kLooseID);
  tree_->SetBranchAddress("j1_e", &j1_e, &b_j1_e);
  tree_->SetBranchAddress("j1_pt", &j1_pt, &b_j1_pt);
  tree_->SetBranchAddress("j1_phi", &j1_phi, &b_j1_phi);
  tree_->SetBranchAddress("j1_eta", &j1_eta, &b_j1_eta);
  tree_->SetBranchAddress("j1_beta", &j1_beta, &b_j1_beta);
  tree_->SetBranchAddress("j1_betaStar", &j1_betaStar, &b_j1_betaStar);
  tree_->SetBranchAddress("j1_betaStarClassic", &j1_betaStarClassic, &b_j1_betaStarClassic);
  tree_->SetBranchAddress("j1_dR2Mean", &j1_dR2Mean, &b_j1_dR2Mean);
  tree_->SetBranchAddress("j1_csvBtag", &j1_csvBtag, &b_j1_csvBtag);
  tree_->SetBranchAddress("j1_csvMvaBtag", &j1_csvMvaBtag, &b_j1_csvMvaBtag);
  tree_->SetBranchAddress("j1_jetProbBtag", &j1_jetProbBtag, &b_j1_jetProbBtag);
  tree_->SetBranchAddress("j1_tcheBtag", &j1_tcheBtag, &b_j1_tcheBtag);
  tree_->SetBranchAddress("j1_radionMatched", &j1_radionMatched, &b_j1_radionMatched);
  tree_->SetBranchAddress("j1_ptD", &j1_ptD, &b_j1_ptD);
  tree_->SetBranchAddress("j1_nSecondaryVertices", &j1_nSecondaryVertices, &b_j1_nSecondaryVertices);
  tree_->SetBranchAddress("j1_secVtxPt", &j1_secVtxPt, &b_j1_secVtxPt);
  tree_->SetBranchAddress("j1_secVtx3dL", &j1_secVtx3dL, &b_j1_secVtx3dL);
  tree_->SetBranchAddress("j1_secVtx3deL", &j1_secVtx3deL, &b_j1_secVtx3deL);
  tree_->SetBranchAddress("j1_emfrac", &j1_emfrac, &b_j1_emfrac);
  tree_->SetBranchAddress("j1_hadfrac", &j1_hadfrac, &b_j1_hadfrac);
  tree_->SetBranchAddress("j1_ntk", &j1_ntk, &b_j1_ntk);
  tree_->SetBranchAddress("j1_nNeutrals", &j1_nNeutrals, &b_j1_nNeutrals);
  tree_->SetBranchAddress("j1_nCharged", &j1_nCharged, &b_j1_nCharged);
  tree_->SetBranchAddress("j1_axis1", &j1_axis1, &b_j1_axis1);
  tree_->SetBranchAddress("j1_axis2", &j1_axis2, &b_j1_axis2);
  tree_->SetBranchAddress("j1_pull", &j1_pull, &b_j1_pull);
  tree_->SetBranchAddress("j2_e", &j2_e, &b_j2_e);
  tree_->SetBranchAddress("j2_pt", &j2_pt, &b_j2_pt);
  tree_->SetBranchAddress("j2_phi", &j2_phi, &b_j2_phi);
  tree_->SetBranchAddress("j2_eta", &j2_eta, &b_j2_eta);
  tree_->SetBranchAddress("j2_beta", &j2_beta, &b_j2_beta);
  tree_->SetBranchAddress("j2_betaStar", &j2_betaStar, &b_j2_betaStar);
  tree_->SetBranchAddress("j2_betaStarClassic", &j2_betaStarClassic, &b_j2_betaStarClassic);
  tree_->SetBranchAddress("j2_dR2Mean", &j2_dR2Mean, &b_j2_dR2Mean);
  tree_->SetBranchAddress("j2_csvBtag", &j2_csvBtag, &b_j2_csvBtag);
  tree_->SetBranchAddress("j2_csvMvaBtag", &j2_csvMvaBtag, &b_j2_csvMvaBtag);
  tree_->SetBranchAddress("j2_jetProbBtag", &j2_jetProbBtag, &b_j2_jetProbBtag);
  tree_->SetBranchAddress("j2_tcheBtag", &j2_tcheBtag, &b_j2_tcheBtag);
  tree_->SetBranchAddress("j2_radionMatched", &j2_radionMatched, &b_j2_radionMatched);
  tree_->SetBranchAddress("j2_ptD", &j2_ptD, &b_j2_ptD);
  tree_->SetBranchAddress("j2_nSecondaryVertices", &j2_nSecondaryVertices, &b_j2_nSecondaryVertices);
  tree_->SetBranchAddress("j2_secVtxPt", &j2_secVtxPt, &b_j2_secVtxPt);
  tree_->SetBranchAddress("j2_secVtx3dL", &j2_secVtx3dL, &b_j2_secVtx3dL);
  tree_->SetBranchAddress("j2_secVtx3deL", &j2_secVtx3deL, &b_j2_secVtx3deL);
  tree_->SetBranchAddress("j2_emfrac", &j2_emfrac, &b_j2_emfrac);
  tree_->SetBranchAddress("j2_hadfrac", &j2_hadfrac, &b_j2_hadfrac);
  tree_->SetBranchAddress("j2_ntk", &j2_ntk, &b_j2_ntk);
  tree_->SetBranchAddress("j2_nNeutrals", &j2_nNeutrals, &b_j2_nNeutrals);
  tree_->SetBranchAddress("j2_nCharged", &j2_nCharged, &b_j2_nCharged);
  tree_->SetBranchAddress("j2_axis1", &j2_axis1, &b_j2_axis1);
  tree_->SetBranchAddress("j2_axis2", &j2_axis2, &b_j2_axis2);
  tree_->SetBranchAddress("j2_pull", &j2_pull, &b_j2_pull);
  tree_->SetBranchAddress("j3_e", &j3_e, &b_j3_e);
  tree_->SetBranchAddress("j3_pt", &j3_pt, &b_j3_pt);
  tree_->SetBranchAddress("j3_phi", &j3_phi, &b_j3_phi);
  tree_->SetBranchAddress("j3_eta", &j3_eta, &b_j3_eta);
  tree_->SetBranchAddress("j3_beta", &j3_beta, &b_j3_beta);
  tree_->SetBranchAddress("j3_betaStar", &j3_betaStar, &b_j3_betaStar);
  tree_->SetBranchAddress("j3_betaStarClassic", &j3_betaStarClassic, &b_j3_betaStarClassic);
  tree_->SetBranchAddress("j3_dR2Mean", &j3_dR2Mean, &b_j3_dR2Mean);
  tree_->SetBranchAddress("j3_csvBtag", &j3_csvBtag, &b_j3_csvBtag);
  tree_->SetBranchAddress("j3_csvMvaBtag", &j3_csvMvaBtag, &b_j3_csvMvaBtag);
  tree_->SetBranchAddress("j3_jetProbBtag", &j3_jetProbBtag, &b_j3_jetProbBtag);
  tree_->SetBranchAddress("j3_tcheBtag", &j3_tcheBtag, &b_j3_tcheBtag);
  tree_->SetBranchAddress("j3_radionMatched", &j3_radionMatched, &b_j3_radionMatched);
  tree_->SetBranchAddress("j3_ptD", &j3_ptD, &b_j3_ptD);
  tree_->SetBranchAddress("j3_nSecondaryVertices", &j3_nSecondaryVertices, &b_j3_nSecondaryVertices);
  tree_->SetBranchAddress("j3_secVtxPt", &j3_secVtxPt, &b_j3_secVtxPt);
  tree_->SetBranchAddress("j3_secVtx3dL", &j3_secVtx3dL, &b_j3_secVtx3dL);
  tree_->SetBranchAddress("j3_secVtx3deL", &j3_secVtx3deL, &b_j3_secVtx3deL);
  tree_->SetBranchAddress("j3_emfrac", &j3_emfrac, &b_j3_emfrac);
  tree_->SetBranchAddress("j3_hadfrac", &j3_hadfrac, &b_j3_hadfrac);
  tree_->SetBranchAddress("j3_ntk", &j3_ntk, &b_j3_ntk);
  tree_->SetBranchAddress("j3_nNeutrals", &j3_nNeutrals, &b_j3_nNeutrals);
  tree_->SetBranchAddress("j3_nCharged", &j3_nCharged, &b_j3_nCharged);
  tree_->SetBranchAddress("j3_axis1", &j3_axis1, &b_j3_axis1);
  tree_->SetBranchAddress("j3_axis2", &j3_axis2, &b_j3_axis2);
  tree_->SetBranchAddress("j3_pull", &j3_pull, &b_j3_pull);
  tree_->SetBranchAddress("j4_e", &j4_e, &b_j4_e);
  tree_->SetBranchAddress("j4_pt", &j4_pt, &b_j4_pt);
  tree_->SetBranchAddress("j4_phi", &j4_phi, &b_j4_phi);
  tree_->SetBranchAddress("j4_eta", &j4_eta, &b_j4_eta);
  tree_->SetBranchAddress("j4_beta", &j4_beta, &b_j4_beta);
  tree_->SetBranchAddress("j4_betaStar", &j4_betaStar, &b_j4_betaStar);
  tree_->SetBranchAddress("j4_betaStarClassic", &j4_betaStarClassic, &b_j4_betaStarClassic);
  tree_->SetBranchAddress("j4_dR2Mean", &j4_dR2Mean, &b_j4_dR2Mean);
  tree_->SetBranchAddress("j4_csvBtag", &j4_csvBtag, &b_j4_csvBtag);
  tree_->SetBranchAddress("j4_csvMvaBtag", &j4_csvMvaBtag, &b_j4_csvMvaBtag);
  tree_->SetBranchAddress("j4_jetProbBtag", &j4_jetProbBtag, &b_j4_jetProbBtag);
  tree_->SetBranchAddress("j4_tcheBtag", &j4_tcheBtag, &b_j4_tcheBtag);
  tree_->SetBranchAddress("j4_radionMatched", &j4_radionMatched, &b_j4_radionMatched);
  tree_->SetBranchAddress("j4_ptD", &j4_ptD, &b_j4_ptD);
  tree_->SetBranchAddress("j4_nSecondaryVertices", &j4_nSecondaryVertices, &b_j4_nSecondaryVertices);
  tree_->SetBranchAddress("j4_secVtxPt", &j4_secVtxPt, &b_j4_secVtxPt);
  tree_->SetBranchAddress("j4_secVtx3dL", &j4_secVtx3dL, &b_j4_secVtx3dL);
  tree_->SetBranchAddress("j4_secVtx3deL", &j4_secVtx3deL, &b_j4_secVtx3deL);
  tree_->SetBranchAddress("j4_emfrac", &j4_emfrac, &b_j4_emfrac);
  tree_->SetBranchAddress("j4_hadfrac", &j4_hadfrac, &b_j4_hadfrac);
  tree_->SetBranchAddress("j4_ntk", &j4_ntk, &b_j4_ntk);
  tree_->SetBranchAddress("j4_nNeutrals", &j4_nNeutrals, &b_j4_nNeutrals);
  tree_->SetBranchAddress("j4_nCharged", &j4_nCharged, &b_j4_nCharged);
  tree_->SetBranchAddress("j4_axis1", &j4_axis1, &b_j4_axis1);
  tree_->SetBranchAddress("j4_axis2", &j4_axis2, &b_j4_axis2);
  tree_->SetBranchAddress("j4_pull", &j4_pull, &b_j4_pull);
  tree_->SetBranchAddress("JetsMass", &JetsMass, &b_JetsMass);
  tree_->SetBranchAddress("dijet_E", &dijet_E, &b_dijet_E);
  tree_->SetBranchAddress("dijet_Pt", &dijet_Pt, &b_dijet_Pt);
  tree_->SetBranchAddress("dijet_Eta", &dijet_Eta, &b_dijet_Eta);
  tree_->SetBranchAddress("dijet_Phi", &dijet_Phi, &b_dijet_Phi);
  tree_->SetBranchAddress("RadMass", &RadMass, &b_RadMass);
  tree_->SetBranchAddress("radion_E", &radion_E, &b_radion_E);
  tree_->SetBranchAddress("radion_Pt", &radion_Pt, &b_radion_Pt);
  tree_->SetBranchAddress("radion_Eta", &radion_Eta, &b_radion_Eta);
  tree_->SetBranchAddress("radion_Phi", &radion_Phi, &b_radion_Phi);
  tree_->SetBranchAddress("gr_radion_p4_pt", &gr_radion_p4_pt, &b_gr_radion_p4_pt);
  tree_->SetBranchAddress("gr_radion_p4_eta", &gr_radion_p4_eta, &b_gr_radion_p4_eta);
  tree_->SetBranchAddress("gr_radion_p4_phi", &gr_radion_p4_phi, &b_gr_radion_p4_phi);
  tree_->SetBranchAddress("gr_radion_p4_mass", &gr_radion_p4_mass, &b_gr_radion_p4_mass);
  tree_->SetBranchAddress("gr_hgg_p4_pt", &gr_hgg_p4_pt, &b_gr_hgg_p4_pt);
  tree_->SetBranchAddress("gr_hgg_p4_eta", &gr_hgg_p4_eta, &b_gr_hgg_p4_eta);
  tree_->SetBranchAddress("gr_hgg_p4_phi", &gr_hgg_p4_phi, &b_gr_hgg_p4_phi);
  tree_->SetBranchAddress("gr_hgg_p4_mass", &gr_hgg_p4_mass, &b_gr_hgg_p4_mass);
  tree_->SetBranchAddress("gr_hbb_p4_pt", &gr_hbb_p4_pt, &b_gr_hbb_p4_pt);
  tree_->SetBranchAddress("gr_hbb_p4_eta", &gr_hbb_p4_eta, &b_gr_hbb_p4_eta);
  tree_->SetBranchAddress("gr_hbb_p4_phi", &gr_hbb_p4_phi, &b_gr_hbb_p4_phi);
  tree_->SetBranchAddress("gr_hbb_p4_mass", &gr_hbb_p4_mass, &b_gr_hbb_p4_mass);
  tree_->SetBranchAddress("gr_g1_p4_pt", &gr_g1_p4_pt, &b_gr_g1_p4_pt);
  tree_->SetBranchAddress("gr_g1_p4_eta", &gr_g1_p4_eta, &b_gr_g1_p4_eta);
  tree_->SetBranchAddress("gr_g1_p4_phi", &gr_g1_p4_phi, &b_gr_g1_p4_phi);
  tree_->SetBranchAddress("gr_g1_p4_mass", &gr_g1_p4_mass, &b_gr_g1_p4_mass);
  tree_->SetBranchAddress("gr_g2_p4_pt", &gr_g2_p4_pt, &b_gr_g2_p4_pt);
  tree_->SetBranchAddress("gr_g2_p4_eta", &gr_g2_p4_eta, &b_gr_g2_p4_eta);
  tree_->SetBranchAddress("gr_g2_p4_phi", &gr_g2_p4_phi, &b_gr_g2_p4_phi);
  tree_->SetBranchAddress("gr_g2_p4_mass", &gr_g2_p4_mass, &b_gr_g2_p4_mass);
  tree_->SetBranchAddress("gr_b1_p4_pt", &gr_b1_p4_pt, &b_gr_b1_p4_pt);
  tree_->SetBranchAddress("gr_b1_p4_eta", &gr_b1_p4_eta, &b_gr_b1_p4_eta);
  tree_->SetBranchAddress("gr_b1_p4_phi", &gr_b1_p4_phi, &b_gr_b1_p4_phi);
  tree_->SetBranchAddress("gr_b1_p4_mass", &gr_b1_p4_mass, &b_gr_b1_p4_mass);
  tree_->SetBranchAddress("gr_b2_p4_pt", &gr_b2_p4_pt, &b_gr_b2_p4_pt);
  tree_->SetBranchAddress("gr_b2_p4_eta", &gr_b2_p4_eta, &b_gr_b2_p4_eta);
  tree_->SetBranchAddress("gr_b2_p4_phi", &gr_b2_p4_phi, &b_gr_b2_p4_phi);
  tree_->SetBranchAddress("gr_b2_p4_mass", &gr_b2_p4_mass, &b_gr_b2_p4_mass);
  tree_->SetBranchAddress("gr_j1_p4_pt", &gr_j1_p4_pt, &b_gr_j1_p4_pt);
  tree_->SetBranchAddress("gr_j1_p4_eta", &gr_j1_p4_eta, &b_gr_j1_p4_eta);
  tree_->SetBranchAddress("gr_j1_p4_phi", &gr_j1_p4_phi, &b_gr_j1_p4_phi);
  tree_->SetBranchAddress("gr_j1_p4_mass", &gr_j1_p4_mass, &b_gr_j1_p4_mass);
  tree_->SetBranchAddress("gr_j2_p4_pt", &gr_j2_p4_pt, &b_gr_j2_p4_pt);
  tree_->SetBranchAddress("gr_j2_p4_eta", &gr_j2_p4_eta, &b_gr_j2_p4_eta);
  tree_->SetBranchAddress("gr_j2_p4_phi", &gr_j2_p4_phi, &b_gr_j2_p4_phi);
  tree_->SetBranchAddress("gr_j2_p4_mass", &gr_j2_p4_mass, &b_gr_j2_p4_mass);
}

void bestJets_commonNtp::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;
  
  // default values                                                                                                     
  dopureeventWeight_ = true;
  cicselection = 4;

  ptphot1cut = 40.;
  ptphot2cut = 25.;

  ptjetacccut  = 25.;
  etajetacccut = 2.5;  
}

std::pair<int,int> bestJets_commonNtp::myTwoHighestPtJet(std::vector<int> jets) {

  float secondJpt = -999.;
  float firstJpt  = -998.;
  int secondJ     = -999;
  int firstJ      = -999;
  
  for (int ij=0; ij<int(jets.size());ij++) {
    int jetIndex = jets[ij]; 

    if (ptcorrjet[jetIndex]>firstJpt) {
      secondJpt = firstJpt;
      firstJpt  = ptcorrjet[jetIndex];
      secondJ   = firstJ;
      firstJ    = jetIndex;
    } else if (ptcorrjet[ij]>secondJpt) {
      secondJpt = ptcorrjet[jetIndex];
      secondJ   = jetIndex;
    }
  }

  return make_pair(firstJ,secondJ);
}

std::pair<int,int> bestJets_commonNtp::myTwoHighestBtagJet(std::vector<int> jets) {

  float secondJbt = -999.;
  float firstJbt  = -998.;
  int secondJ     = -999;
  int firstJ      = -999;
  
  for (int ij=0; ij<int(jets.size());ij++) {
    int jetIndex = jets[ij]; 
    if (btagcsvjet[jetIndex]>firstJbt) {
      secondJbt = firstJbt;
      firstJbt  = btagcsvjet[jetIndex];
      secondJ   = firstJ;
      firstJ    = jetIndex;
    } else if (btagcsvjet[ij]>secondJbt) {
      secondJbt = btagcsvjet[jetIndex];
      secondJ   = jetIndex;
    }
  }

  return make_pair(firstJ,secondJ);
}

int bestJets_commonNtp::myHighestPtJet(std::vector<int> jets) {

  float firstJpt = -998.;
  int firstJ     = -999;
  
  for (int ij=0; ij<int(jets.size());ij++) {
    int jetIndex = jets[ij]; 
    if (ptcorrjet[jetIndex]>firstJpt) {
      firstJpt  = ptcorrjet[jetIndex];
      firstJ    = jetIndex;
    } 
  }

  return firstJ;
}

int bestJets_commonNtp::myHighestBtagJet(std::vector<int> jets) {

  float firstJbt = -998.;
  int firstJ     = -999;
  
  for (int ij=0; ij<int(jets.size());ij++) {
    int jetIndex = jets[ij]; 
    if (btagcsvjet[jetIndex]>firstJbt) {
      firstJbt  = btagcsvjet[jetIndex];
      firstJ    = jetIndex;
    } 
  }

  return firstJ;
}

bool bestJets_commonNtp::passCutBasedJetId(int jet) {

  bool isGood = true;

  float thebetastarjet[4], thermsjet[4];  
  thebetastarjet[0] = j1_betaStarClassic;
  thebetastarjet[1] = j2_betaStarClassic;
  thebetastarjet[2] = j3_betaStarClassic;
  thebetastarjet[3] = j4_betaStarClassic;
  thermsjet[0]      = j1_dR2Mean;
  thermsjet[1]      = j2_dR2Mean;
  thermsjet[2]      = j3_dR2Mean;
  thermsjet[3]      = j4_dR2Mean;

  if ( ptcorrjet[jet]<-900) isGood = false;

  if ( fabs(etajet[jet]) < 2.5 ) {
    if ( thebetastarjet[jet] > 0.2 * log( nvtx - 0.64) )  isGood = false;
    if (thermsjet[jet] > 0.06)                            isGood = false;
  } else if (fabs(etajet[jet]) < 3.){
    if ( thermsjet[jet] > 0.05)  isGood =false;
  } else {
    if ( thermsjet[jet] > 0.055) isGood =false;
  }

  return isGood;
}
