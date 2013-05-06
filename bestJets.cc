#include "bestJets.h"

using namespace std;

bestJets::bestJets( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer( "Radion", dataset ) {
  
  bTaggerType_ = bTaggerType;
  
  setSelectionType(selectionType);
}

bestJets::~bestJets() {
  
  outFile_->Close();
  
  if (!tree_) return;
  delete tree_->GetCurrentFile();
}

double bestJets::delta_phi(double phi1, double phi2) {

  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}

void bestJets::finalize() {
  
  this->Init();
  
  std::string fullFlags = selectionType_ + "_";
  fullFlags+=bTaggerType_;
  this->set_flags(fullFlags); 
  this->createOutputFile();

  outFile_->cd();


  //-----------------------------------------------------------------------------
  // for the tree with selected events
  float invMassJJ_pt, invMassJJ_maxptjj, invMassJJ_maxptmjj;  
  float invMassJJ_mgg, invMassJJ_btag_pt, invMassJJ_btag_ptCSV;
  float invMassJJ_btag_ptjj, invMassJJ_btag_ptmjj, invMassJJ_btag_mgg;
  int btagCategory_t;
  float weight;
  float dMmin_t, dMmin_btag_t;
  int ngoodJets_t;

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
  myTrees->Branch( "dMmin",           &dMmin_t,               "dMmin_t/F" );    
  myTrees->Branch( "dMmin_btag",      &dMmin_btag_t,          "dMmin_btag_t/F" );    
  myTrees->Branch( "btagCategory",    &btagCategory_t,        "btagCategory_t/I" );
  myTrees->Branch( "ngoodJets",       &ngoodJets_t,           "ngoodJets_t/I" );
  myTrees->Branch( "weight", &weight, "weight/F" );

  // ------------------------------------------------------
  // analysis
  if (tree_ == 0) return;

  Long64_t nentries = tree_->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  if( xSection_!=1. && xSection_!=0. ) {
    std::cout << std::endl;
    std::cout << "-> Cross-Section: "      << xSection_             << " pb" << std::endl;
    std::cout << "-> # Generated Events: " << nGenEvents_           << std::endl;
    std::cout << "-> Event Weight: "       << xSection_/nGenEvents_ << std::endl << std::endl;
  } 

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = tree_->GetEntry(jentry);   nbytes += nb;
    
    // event weight for MC events: xsec and PU
    bool isMC = ( run<5 );
    if( isMC ) {
      if( nGenEvents_<=0 ) {
	std::cout << std::endl << "-> There must be a problem: nGenEvents is not positive!! Exiting" << std::endl;
	exit(99);
      }
      weight = xSection_ / nGenEvents_ ;

      // pu reWeighting
      if(dopureeventWeight_) weight *= pu_weight;
    }


    // ---------------------------------------------------
    // gamma-gamma analysis

    // photons acceptance
    if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
       || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance
    
    // photon id
    bool idphot1(0), idphot2(0), pxlphot1(1), pxlphot2(1);
    idphot1 = (idcicpfphot1 >= cicselection);
    idphot2 = (idcicpfphot2 >= cicselection);
    if(!(idphot1)) continue;
    if(!(idphot2)) continue;
    
    // extra cuts on photons: splitting events per photon class if needed
    int myEB=2;
    if( (fabs(etascphot1)>1.4442 || fabs(etascphot2)>1.4442) ) myEB=0;
    if( (fabs(etascphot1)<1.4442 && fabs(etascphot2)<1.4442) ) myEB=1;

    // the two selected photons for the analysis
    TLorentzVector t4phot1, t4phot2;  
    t4phot1.SetPtEtaPhiM(ptphot1,etaphot1,phiphot1,0.);
    t4phot2.SetPtEtaPhiM(ptphot2,etaphot2,phiphot2,0.);
    TLorentzVector t4diPhot;
    t4diPhot.SetPtEtaPhiM( ptggnewvtx, etagg, phigg, massggnewvtx );
    TVector3 t3diPhot;
    t3diPhot.SetPtEtaPhi( ptggnewvtx, etagg, phigg );

    // further cuts on photons -----------------------

    // invariant mass cut on photons
    if (t4diPhot.M()<100 || t4diPhot.M()>180) continue;
    
    // photons pt cuts 
    if(ptphot1<ptphot1cut* massggnewvtx/120.) continue;         // pt first photon
    if(ptphot2<ptphot2cut)                    continue;         // pt second photon


    // --------------------------------------------------------------------------------
    // jets, no btagging, passing cut based jetID, restricting to |eta|<2.5
    vector<int> v_puIdJets;
    for (int ij=0; ij<njets; ij++) { 
      if ( ptcorrjet[ij] < ptjetacccut )   continue;   
      if ( fabs(etajet[ij])>etajetacccut ) continue;   
      if ( !passCutBasedJetId(ij) )        continue;
      v_puIdJets.push_back(ij);
    }

    // jets passing btagging ( + eta/pT/PUid cuts)
    vector<int> v_looseJP,  v_mediumJP,  v_tightJP; 
    vector<int> v_looseCSV, v_mediumCSV, v_tightCSV;
    for (int ij=0; ij<int(v_puIdJets.size());ij++) {   
      int index = v_puIdJets[ij];
      if (btagjprobjet[index]>0.275) v_looseJP.push_back(index);
      if (btagjprobjet[index]>0.545) v_mediumJP.push_back(index);
      if (btagjprobjet[index]>0.790) v_tightJP.push_back(index);
      if (btagcsvjet[index]>0.244)   v_looseCSV.push_back(index);
      if (btagcsvjet[index]>0.679)   v_mediumCSV.push_back(index);
      if (btagcsvjet[index]>0.898)   v_tightCSV.push_back(index);
    }

    // at least 2 preselected jets
    if (v_puIdJets.size()<2) continue;   
    std::pair<int,int> allJets_highPt = myTwoHighestPtJet(v_puIdJets);




    // choice of analysis jets ---------------------

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

    // -----------------------------------------

    // choose jets looking to btag: pT
    int jet1btag_pt = -1;
    int jet2btag_pt = -1;
    if (v_looseCSV.size()==1) { 
      jet1btag_pt = myHighestPtJet(v_looseCSV);      
    } else if (v_looseCSV.size()>1) {
      std::pair<int,int> looseBtag_highPt = myTwoHighestPtJet(v_looseCSV);
      jet1btag_pt = looseBtag_highPt.first;
      jet2btag_pt = looseBtag_highPt.second;
    }

    // choose jets looking to btag: CSV
    int jet1btag_CSV = -1;
    int jet2btag_CSV = -1;
    if (v_looseCSV.size()==1) { 
      jet1btag_CSV = myHighestBtagJet(v_looseCSV);                          
    } else if (v_looseCSV.size()>1) {
      std::pair<int,int> looseBtag_csv = myTwoHighestBtagJet(v_looseCSV);  
      jet1btag_CSV = looseBtag_csv.first;
      jet2btag_CSV = looseBtag_csv.second;
    }

    // -----------------------------------------

    // 1a) now apply the same criteria as above but giving priorities to the btagged jets
    //     bjets: pT ordered ; nobjets (btagCat=1): pT ordered
    int jet1_btag_pt = -1;
    int jet2_btag_pt = -1;
    if( v_looseCSV.size()==1 ) {
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
    if( v_looseCSV.size()>1 ) {
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
    if( v_looseCSV.size()==1 ) {
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
    if( v_looseCSV.size()>1 ) {
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
    if( v_looseCSV.size()==1 ) {
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
	int index = v_puIdJets[jet];
	if (index==jet1btag_CSV) continue;
	TLorentzVector t4jet, t4bjet;
	t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
	t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag_CSV],etajet[jet1btag_CSV],phijet[jet1btag_CSV],ecorrjet[jet1btag_CSV]);      
	TLorentzVector t4_bnotb = t4jet + t4bjet; 
	if ( t4_bnotb.Pt() > maxPt2) {
	  maxPt2 = t4_bnotb.Pt();
	  jet1_btag_ptjj = jet1btag_CSV;
	  jet2_btag_ptjj = index;
	}
      }
    }
    else if( v_looseCSV.size()>1 ) {
      for (int jetA=0; jetA<(v_looseCSV.size()-1); jetA++) { 
	for (int jetB=jetA+1; jetB<v_looseCSV.size(); jetB++) { 
	  TLorentzVector t4jetA, t4jetB;
	  int indexA = v_looseCSV[jetA];
	  int indexB = v_looseCSV[jetB];
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
    if( v_looseCSV.size()==1 ) {
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
	int index = v_puIdJets[jet];
	if (index==jet1btag_CSV) continue;
	TLorentzVector t4jet, t4bjet;
	t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
	t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag_CSV],etajet[jet1btag_CSV],phijet[jet1btag_CSV],ecorrjet[jet1btag_CSV]);      
	TLorentzVector t4_bnotb = t4jet + t4bjet; 
	if ( (t4_bnotb.Pt()/t4_bnotb.M()) > maxPtMjj2) {
	  maxPtMjj2 = t4_bnotb.Pt()/t4_bnotb.M();
	  jet1_btag_ptmjj = jet1btag_CSV;
	  jet2_btag_ptmjj = index;
	}
      }
    }
    else if( v_looseCSV.size()>1 ) {
      for (int jetA=0; jetA<(v_looseCSV.size()-1); jetA++) { 
	for (int jetB=jetA+1; jetB<v_looseCSV.size(); jetB++) { 
	  TLorentzVector t4jetA, t4jetB;
	  int indexA = v_looseCSV[jetA];
	  int indexB = v_looseCSV[jetB];
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
    if( v_looseCSV.size()==1 ) {
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
	int index = v_puIdJets[jet];
	if (index==jet1btag_CSV) continue;
	TLorentzVector t4jet, t4bjet;
	t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
	t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag_CSV],etajet[jet1btag_CSV],phijet[jet1btag_CSV],ecorrjet[jet1btag_CSV]);      
	TLorentzVector t4_bnotb = t4jet + t4bjet; 
	float massJJ = t4_bnotb.M();
	float massGG = t4diPhot.M();
	if ( fabs(massJJ-massGG) < minDMgg2 ) {
	  minDMgg2 = fabs(massJJ-massGG);
	  jet1_btag_mgg = jet1btag_CSV;
	  jet2_btag_mgg = index;
	}
      }
    } else if (v_looseCSV.size()>1 ) {
      for (int jetA=0; jetA<(v_looseCSV.size()-1); jetA++) { 
	for (int jetB=jetA+1; jetB<v_looseCSV.size(); jetB++) { 
	  TLorentzVector t4jetA, t4jetB;
	  int indexA = v_looseCSV[jetA];
	  int indexB = v_looseCSV[jetB];
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

    
    // several invariant masses cut on jets 
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

    // ---------------------------------------------------------------
    // counting the number of bjets to categorize the events
    int btagCategory = -1;
    if (bTaggerType_=="JP")       { btagCategory = (v_looseJP.size()<=2) ? v_looseJP.size() : 2; }
    else if (bTaggerType_=="CSV") { btagCategory = (v_looseCSV.size()<=2) ? v_looseCSV.size() : 2; }
    else cout << "this btag algo does not exist" << endl;


    // --------------------------------------------------------------
    // here we apply further cuts according to the btag category
    if (btagCategory<1 || btagCategory>2) continue;


    // filling the tree for selected events 
    btagCategory_t = btagCategory;
    ngoodJets_t    = v_puIdJets.size(); 
    dMmin_t        = minDMgg;
    dMmin_btag_t   = minDMgg2;

    myTrees->Fill();

  } // loop over entries 


  outFile_->cd();
  myTrees->Write();

} // finalize

Int_t bestJets::GetEntry(Long64_t entry) {

  if (!tree_) return 0;
  return tree_->GetEntry(entry);
}

Long64_t bestJets::LoadTree(Long64_t entry) {

  if (!tree_) return -5;
  Long64_t centry = tree_->LoadTree(entry);
  if (centry < 0) return centry;
  if (tree_->GetTreeNumber() != fCurrent) {
    fCurrent = tree_->GetTreeNumber();
  }
  return centry;
}

void bestJets::Init() {

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

  tree_->SetBranchAddress("run", &run, &b_run);
  tree_->SetBranchAddress("event", &event, &b_event);
  tree_->SetBranchAddress("lumi", &lumi, &b_lumi);
  tree_->SetBranchAddress("H_event", &H_event, &b_H_event);
  tree_->SetBranchAddress("V_event", &V_event, &b_V_event);
  tree_->SetBranchAddress("WH_event", &WH_event, &b_WH_event);
  tree_->SetBranchAddress("ZH_event", &ZH_event, &b_ZH_event);
  tree_->SetBranchAddress("Zbb_event", &Zbb_event, &b_Zbb_event);
  tree_->SetBranchAddress("Vqq_event", &Vqq_event, &b_Vqq_event);
  tree_->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
  tree_->SetBranchAddress("massgg", &massgg, &b_massgg);
  tree_->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
  tree_->SetBranchAddress("ptggnewvtx", &ptggnewvtx, &b_ptggnewvtx);
  tree_->SetBranchAddress("phigg", &phigg, &b_phigg);
  tree_->SetBranchAddress("etagg", &etagg, &b_etagg);
  tree_->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx);
  tree_->SetBranchAddress("ptphot1", &ptphot1, &b_ptphot1);
  tree_->SetBranchAddress("ptphot2", &ptphot2, &b_ptphot2);
  tree_->SetBranchAddress("deltaRToTrackphot1", &deltaRToTrackphot1, &b_deltaRToTrackphot1);
  tree_->SetBranchAddress("deltaRToTrackphot2", &deltaRToTrackphot2, &b_deltaRToTrackphot2);
  tree_->SetBranchAddress("timephot1", &timephot1, &b_timephot1);
  tree_->SetBranchAddress("timephot2", &timephot2, &b_timephot2);
  tree_->SetBranchAddress("etaphot1", &etaphot1, &b_etaphot1);
  tree_->SetBranchAddress("etaphot2", &etaphot2, &b_etaphot2);
  tree_->SetBranchAddress("phiphot1", &phiphot1, &b_phiphot1);
  tree_->SetBranchAddress("phiphot2", &phiphot2, &b_phiphot2);
  tree_->SetBranchAddress("etascphot1", &etascphot1, &b_etascphot1);
  tree_->SetBranchAddress("etascphot2", &etascphot2, &b_etascphot2);
  tree_->SetBranchAddress("phiscphot1", &phiscphot1, &b_phiscphot1);
  tree_->SetBranchAddress("phiscphot2", &phiscphot2, &b_phiscphot2);
  tree_->SetBranchAddress("E1phot1", &E1phot1, &b_E1phot1);
  tree_->SetBranchAddress("E1phot2", &E1phot2, &b_E1phot2);
  tree_->SetBranchAddress("E9phot1", &E9phot1, &b_E9phot1);
  tree_->SetBranchAddress("E9phot2", &E9phot2, &b_E9phot2);
  tree_->SetBranchAddress("energyErrphot1", &energyErrphot1, &b_energyErrphot1);
  tree_->SetBranchAddress("energyErrphot2", &energyErrphot2, &b_energyErrphot2);
  tree_->SetBranchAddress("energySmearingphot1", &energySmearingphot1, &b_energySmearingphot1);
  tree_->SetBranchAddress("energySmearingphot2", &energySmearingphot2, &b_energySmearingphot2);
  tree_->SetBranchAddress("r9phot1", &r9phot1, &b_r9phot1);
  tree_->SetBranchAddress("r9phot2", &r9phot2, &b_r9phot2);
  tree_->SetBranchAddress("isemEGphot1", &isemEGphot1, &b_isemEGphot1);
  tree_->SetBranchAddress("isemEGphot2", &isemEGphot2, &b_isemEGphot2);
  tree_->SetBranchAddress("promptGamma", &promptGamma, &b_promptGamma);
  tree_->SetBranchAddress("LOGamma", &LOGamma, &b_LOGamma);
  tree_->SetBranchAddress("ISRGamma", &ISRGamma, &b_ISRGamma);
  tree_->SetBranchAddress("FSRGamma", &FSRGamma, &b_FSRGamma);
  tree_->SetBranchAddress("idmvaphot1", &idmvaphot1, &b_idmvaphot1);
  tree_->SetBranchAddress("idmvaphot2", &idmvaphot2, &b_idmvaphot2);
  tree_->SetBranchAddress("idcicphot1", &idcicphot1, &b_idcicphot1);
  tree_->SetBranchAddress("idcicphot2", &idcicphot2, &b_idcicphot2);
  tree_->SetBranchAddress("idcicnoelvetophot1", &idcicnoelvetophot1, &b_idcicnoelvetophot1);
  tree_->SetBranchAddress("idcicnoelvetophot2", &idcicnoelvetophot2, &b_idcicnoelvetophot2);
  tree_->SetBranchAddress("idcicpfphot1", &idcicpfphot1, &b_idcicpfphot1);
  tree_->SetBranchAddress("idcicpfphot2", &idcicpfphot2, &b_idcicpfphot2);
  tree_->SetBranchAddress("idcicpfnoelvetophot1", &idcicpfnoelvetophot1, &b_idcicpfnoelvetophot1);
  tree_->SetBranchAddress("idcicpfnoelvetophot2", &idcicpfnoelvetophot2, &b_idcicpfnoelvetophot2);
  tree_->SetBranchAddress("idelephot1", &idelephot1, &b_idelephot1);
  tree_->SetBranchAddress("idelephot2", &idelephot2, &b_idelephot2);
  tree_->SetBranchAddress("pid_isEMphot1", &pid_isEMphot1, &b_pid_isEMphot1);
  tree_->SetBranchAddress("pid_isEMphot2", &pid_isEMphot2, &b_pid_isEMphot2);
  tree_->SetBranchAddress("pid_haspixelseedphot1", &pid_haspixelseedphot1, &b_pid_haspixelseedphot1);
  tree_->SetBranchAddress("pid_haspixelseedphot2", &pid_haspixelseedphot2, &b_pid_haspixelseedphot2);
  tree_->SetBranchAddress("pid_jurECALphot1", &pid_jurECALphot1, &b_pid_jurECALphot1);
  tree_->SetBranchAddress("pid_jurECALphot2", &pid_jurECALphot2, &b_pid_jurECALphot2);
  tree_->SetBranchAddress("pid_twrHCALphot1", &pid_twrHCALphot1, &b_pid_twrHCALphot1);
  tree_->SetBranchAddress("pid_twrHCALphot2", &pid_twrHCALphot2, &b_pid_twrHCALphot2);
  tree_->SetBranchAddress("pid_HoverEphot1", &pid_HoverEphot1, &b_pid_HoverEphot1);
  tree_->SetBranchAddress("pid_HoverEphot2", &pid_HoverEphot2, &b_pid_HoverEphot2);
  tree_->SetBranchAddress("pid_hlwTrackphot1", &pid_hlwTrackphot1, &b_pid_hlwTrackphot1);
  tree_->SetBranchAddress("pid_hlwTrackphot2", &pid_hlwTrackphot2, &b_pid_hlwTrackphot2);
  tree_->SetBranchAddress("pid_etawidphot1", &pid_etawidphot1, &b_pid_etawidphot1);
  tree_->SetBranchAddress("pid_etawidphot2", &pid_etawidphot2, &b_pid_etawidphot2);
  tree_->SetBranchAddress("pid_hlwTrackNoDzphot1", &pid_hlwTrackNoDzphot1, &b_pid_hlwTrackNoDzphot1);
  tree_->SetBranchAddress("pid_hlwTrackNoDzphot2", &pid_hlwTrackNoDzphot2, &b_pid_hlwTrackNoDzphot2);
  tree_->SetBranchAddress("pid_hasMatchedConvphot1", &pid_hasMatchedConvphot1, &b_pid_hasMatchedConvphot1);
  tree_->SetBranchAddress("pid_hasMatchedConvphot2", &pid_hasMatchedConvphot2, &b_pid_hasMatchedConvphot2);
  tree_->SetBranchAddress("pid_hasMatchedPromptElephot1", &pid_hasMatchedPromptElephot1, &b_pid_hasMatchedPromptElephot1);
  tree_->SetBranchAddress("pid_hasMatchedPromptElephot2", &pid_hasMatchedPromptElephot2, &b_pid_hasMatchedPromptElephot2);
  tree_->SetBranchAddress("pid_sminphot1", &pid_sminphot1, &b_pid_sminphot1);
  tree_->SetBranchAddress("pid_sminphot2", &pid_sminphot2, &b_pid_sminphot2);
  tree_->SetBranchAddress("pid_smajphot1", &pid_smajphot1, &b_pid_smajphot1);
  tree_->SetBranchAddress("pid_smajphot2", &pid_smajphot2, &b_pid_smajphot2);
  tree_->SetBranchAddress("pid_ntrkphot1", &pid_ntrkphot1, &b_pid_ntrkphot1);
  tree_->SetBranchAddress("pid_ntrkphot2", &pid_ntrkphot2, &b_pid_ntrkphot2);
  tree_->SetBranchAddress("pid_ptisophot1", &pid_ptisophot1, &b_pid_ptisophot1);
  tree_->SetBranchAddress("pid_ptisophot2", &pid_ptisophot2, &b_pid_ptisophot2);
  tree_->SetBranchAddress("pid_ntrkcsphot1", &pid_ntrkcsphot1, &b_pid_ntrkcsphot1);
  tree_->SetBranchAddress("pid_ntrkcsphot2", &pid_ntrkcsphot2, &b_pid_ntrkcsphot2);
  tree_->SetBranchAddress("pid_ptisocsphot1", &pid_ptisocsphot1, &b_pid_ptisocsphot1);
  tree_->SetBranchAddress("pid_ptisocsphot2", &pid_ptisocsphot2, &b_pid_ptisocsphot2);
  tree_->SetBranchAddress("pid_ecalisophot1", &pid_ecalisophot1, &b_pid_ecalisophot1);
  tree_->SetBranchAddress("pid_ecalisophot2", &pid_ecalisophot2, &b_pid_ecalisophot2);
  tree_->SetBranchAddress("pid_hcalisophot1", &pid_hcalisophot1, &b_pid_hcalisophot1);
  tree_->SetBranchAddress("pid_hcalisophot2", &pid_hcalisophot2, &b_pid_hcalisophot2);
  tree_->SetBranchAddress("njets", &njets, &b_njets);
  tree_->SetBranchAddress("ecorrjet", ecorrjet, &b_ecorrjet);
  tree_->SetBranchAddress("ptjet", ptjet, &b_ptjet);
  tree_->SetBranchAddress("ptcorrjet", ptcorrjet, &b_ptcorrjet);
  tree_->SetBranchAddress("etajet", etajet, &b_etajet);
  tree_->SetBranchAddress("phijet", phijet, &b_phijet);
  tree_->SetBranchAddress("betajet", betajet, &b_betajet);
  tree_->SetBranchAddress("betastarjet", betastarjet, &b_betastarjet);
  tree_->SetBranchAddress("btagvtxjet", btagvtxjet, &b_btagvtxjet);
  tree_->SetBranchAddress("btagcsvjet", btagcsvjet, &b_btagcsvjet);
  tree_->SetBranchAddress("btagjprobjet", btagjprobjet, &b_btagjprobjet);
  tree_->SetBranchAddress("btagcsvmvajet", btagcsvmvajet, &b_btagcsvmvajet);
  tree_->SetBranchAddress("btagbjprobjet", btagbjprobjet, &b_btagbjprobjet);
  tree_->SetBranchAddress("btagvtxhighpurjet", btagvtxhighpurjet, &b_btagvtxhighpurjet);
  tree_->SetBranchAddress("btagtchejet", btagtchejet, &b_btagtchejet);
  tree_->SetBranchAddress("btagtchpjet", btagtchpjet, &b_btagtchpjet);
  tree_->SetBranchAddress("ptDjet", ptDjet, &b_ptDjet);
  tree_->SetBranchAddress("ptD_QCjet", ptD_QCjet, &b_ptD_QCjet);
  tree_->SetBranchAddress("axis2_QCjet", axis2_QCjet, &b_axis2_QCjet);
  tree_->SetBranchAddress("rmsjet", rmsjet, &b_rmsjet);
  tree_->SetBranchAddress("ntrkjet", ntrkjet, &b_ntrkjet);
  tree_->SetBranchAddress("nneutjet", nneutjet, &b_nneutjet);
  tree_->SetBranchAddress("nChg_QCjet", nChg_QCjet, &b_nChg_QCjet);
  tree_->SetBranchAddress("nNeutral_ptCutjet", nNeutral_ptCutjet, &b_nNeutral_ptCutjet);
  tree_->SetBranchAddress("jetIdSimple_mvajet", jetIdSimple_mvajet, &b_jetIdSimple_mvajet);
  tree_->SetBranchAddress("jetIdFull_mvajet", jetIdFull_mvajet, &b_jetIdFull_mvajet);
  tree_->SetBranchAddress("jetId_dR2Meanjet", jetId_dR2Meanjet, &b_jetId_dR2Meanjet);
  tree_->SetBranchAddress("jetId_betaStarClassicjet", jetId_betaStarClassicjet, &b_jetId_betaStarClassicjet);
  tree_->SetBranchAddress("jetId_frac01jet", jetId_frac01jet, &b_jetId_frac01jet);
  tree_->SetBranchAddress("jetId_frac02jet", jetId_frac02jet, &b_jetId_frac02jet);
  tree_->SetBranchAddress("jetId_frac03jet", jetId_frac03jet, &b_jetId_frac03jet);
  tree_->SetBranchAddress("jetId_frac04jet", jetId_frac04jet, &b_jetId_frac04jet);
  tree_->SetBranchAddress("jetId_frac05jet", jetId_frac05jet, &b_jetId_frac05jet);
  tree_->SetBranchAddress("jetId_betajet", jetId_betajet, &b_jetId_betajet);
  tree_->SetBranchAddress("jetId_betaStarjet", jetId_betaStarjet, &b_jetId_betaStarjet);
  tree_->SetBranchAddress("jetIdCutBased_wpjet", jetIdCutBased_wpjet, &b_jetIdCutBased_wpjet);
  tree_->SetBranchAddress("jetIdSimple_wpjet", jetIdSimple_wpjet, &b_jetIdSimple_wpjet);
  tree_->SetBranchAddress("jetIdFull_wpjet", jetIdFull_wpjet, &b_jetIdFull_wpjet);
  tree_->SetBranchAddress("assjet", assjet, &b_assjet);
  tree_->SetBranchAddress("partPdgIDjet", partPdgIDjet, &b_partPdgIDjet);
  tree_->SetBranchAddress("partMomPdgIDjet", partMomPdgIDjet, &b_partMomPdgIDjet);
  tree_->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree_->SetBranchAddress("vtxId", &vtxId, &b_vtxId);
  tree_->SetBranchAddress("vtxPos_x", &vtxPos_x, &b_vtxPos_x);
  tree_->SetBranchAddress("vtxPos_y", &vtxPos_y, &b_vtxPos_y);
  tree_->SetBranchAddress("vtxPos_z", &vtxPos_z, &b_vtxPos_z);
  tree_->SetBranchAddress("vtxIdMVA", &vtxIdMVA, &b_vtxIdMVA);
  tree_->SetBranchAddress("vtxIdEvtProb", &vtxIdEvtProb, &b_vtxIdEvtProb);
  tree_->SetBranchAddress("diPhotMVA", &diPhotMVA, &b_diPhotMVA);
  tree_->SetBranchAddress("diPhotMVA_vtx0", &diPhotMVA_vtx0, &b_diPhotMVA_vtx0);
  tree_->SetBranchAddress("diPhotMVA_vtxPair", &diPhotMVA_vtxPair, &b_diPhotMVA_vtxPair);
  tree_->SetBranchAddress("preselPairId", &preselPairId, &b_preselPairId);
  tree_->SetBranchAddress("tmva_dipho_MIT_dmom", &tmva_dipho_MIT_dmom, &b_tmva_dipho_MIT_dmom);
  tree_->SetBranchAddress("tmva_dipho_MIT_dmom_wrong_vtx", &tmva_dipho_MIT_dmom_wrong_vtx, &b_tmva_dipho_MIT_dmom_wrong_vtx);
  tree_->SetBranchAddress("tmva_dipho_MIT_vtxprob", &tmva_dipho_MIT_vtxprob, &b_tmva_dipho_MIT_vtxprob);
  tree_->SetBranchAddress("tmva_dipho_MIT_ptom1", &tmva_dipho_MIT_ptom1, &b_tmva_dipho_MIT_ptom1);
  tree_->SetBranchAddress("tmva_dipho_MIT_ptom2", &tmva_dipho_MIT_ptom2, &b_tmva_dipho_MIT_ptom2);
  tree_->SetBranchAddress("tmva_dipho_MIT_eta1", &tmva_dipho_MIT_eta1, &b_tmva_dipho_MIT_eta1);
  tree_->SetBranchAddress("tmva_dipho_MIT_eta2", &tmva_dipho_MIT_eta2, &b_tmva_dipho_MIT_eta2);
  tree_->SetBranchAddress("tmva_dipho_MIT_dphi", &tmva_dipho_MIT_dphi, &b_tmva_dipho_MIT_dphi);
  tree_->SetBranchAddress("tmva_dipho_MIT_ph1mva", &tmva_dipho_MIT_ph1mva, &b_tmva_dipho_MIT_ph1mva);
  tree_->SetBranchAddress("tmva_dipho_MIT_ph2mva", &tmva_dipho_MIT_ph2mva, &b_tmva_dipho_MIT_ph2mva);
  tree_->SetBranchAddress("ePUMet", &ePUMet, &b_ePUMet);
  tree_->SetBranchAddress("ePUMet2", &ePUMet2, &b_ePUMet2);
  tree_->SetBranchAddress("ePUMet3", &ePUMet3, &b_ePUMet3);
  tree_->SetBranchAddress("ePUMet4", &ePUMet4, &b_ePUMet4);
  tree_->SetBranchAddress("ePUMet5", &ePUMet5, &b_ePUMet5);
  tree_->SetBranchAddress("ecorrPUMet5", &ecorrPUMet5, &b_ecorrPUMet5);
  tree_->SetBranchAddress("phiPUMet", &phiPUMet, &b_phiPUMet);
  tree_->SetBranchAddress("phiPUMet2", &phiPUMet2, &b_phiPUMet2);
  tree_->SetBranchAddress("phiPUMet3", &phiPUMet3, &b_phiPUMet3);
  tree_->SetBranchAddress("phiPUMet4", &phiPUMet4, &b_phiPUMet4);
  tree_->SetBranchAddress("phiPUMet5", &phiPUMet5, &b_phiPUMet5);
  tree_->SetBranchAddress("phiCorrPUMet5", &phiCorrPUMet5, &b_phiCorrPUMet5);
  tree_->SetBranchAddress("phot1Metx", &phot1Metx, &b_phot1Metx);
  tree_->SetBranchAddress("phot2Metx", &phot2Metx, &b_phot2Metx);
  tree_->SetBranchAddress("leptonsMetx", &leptonsMetx, &b_leptonsMetx);
  tree_->SetBranchAddress("part_in_jetsMetx", &part_in_jetsMetx, &b_part_in_jetsMetx);
  tree_->SetBranchAddress("chg_vtx_unclMetx", &chg_vtx_unclMetx, &b_chg_vtx_unclMetx);
  tree_->SetBranchAddress("chg_novtx_unclMetx", &chg_novtx_unclMetx, &b_chg_novtx_unclMetx);
  tree_->SetBranchAddress("neutrals_unclMetx", &neutrals_unclMetx, &b_neutrals_unclMetx);
  tree_->SetBranchAddress("part_fail_puidMetx", &part_fail_puidMetx, &b_part_fail_puidMetx);
  tree_->SetBranchAddress("phot1Mety", &phot1Mety, &b_phot1Mety);
  tree_->SetBranchAddress("phot2Mety", &phot2Mety, &b_phot2Mety);
  tree_->SetBranchAddress("leptonsMety", &leptonsMety, &b_leptonsMety);
  tree_->SetBranchAddress("part_in_jetsMety", &part_in_jetsMety, &b_part_in_jetsMety);
  tree_->SetBranchAddress("chg_vtx_unclMety", &chg_vtx_unclMety, &b_chg_vtx_unclMety);
  tree_->SetBranchAddress("chg_novtx_unclMety", &chg_novtx_unclMety, &b_chg_novtx_unclMety);
  tree_->SetBranchAddress("neutrals_unclMety", &neutrals_unclMety, &b_neutrals_unclMety);
  tree_->SetBranchAddress("part_fail_puidMety", &part_fail_puidMety, &b_part_fail_puidMety);
  tree_->SetBranchAddress("scaling", &scaling, &b_scaling);
  tree_->SetBranchAddress("sMet", &sMet, &b_sMet);
  tree_->SetBranchAddress("eMet", &eMet, &b_eMet);
  tree_->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
  tree_->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
  tree_->SetBranchAddress("eSmearedMet", &eSmearedMet, &b_eSmearedMet);
  tree_->SetBranchAddress("phiSmearedMet", &phiSmearedMet, &b_phiSmearedMet);
  tree_->SetBranchAddress("eShiftedMet", &eShiftedMet, &b_eShiftedMet);
  tree_->SetBranchAddress("phiShiftedMet", &phiShiftedMet, &b_phiShiftedMet);
  tree_->SetBranchAddress("eShiftedScaledMet", &eShiftedScaledMet, &b_eShiftedScaledMet);
  tree_->SetBranchAddress("phiShiftedScaledMet", &phiShiftedScaledMet, &b_phiShiftedScaledMet);
  tree_->SetBranchAddress("eSmearedShiftedMet", &eSmearedShiftedMet, &b_eSmearedShiftedMet);
  tree_->SetBranchAddress("phiSmearedShiftedMet", &phiSmearedShiftedMet, &b_phiSmearedShiftedMet);
  tree_->SetBranchAddress("eShiftedScaledMetPUcorr", &eShiftedScaledMetPUcorr, &b_eShiftedScaledMetPUcorr);
  tree_->SetBranchAddress("phiShiftedScaledMetPUcorr", &phiShiftedScaledMetPUcorr, &b_phiShiftedScaledMetPUcorr);
  tree_->SetBranchAddress("eSmearedShiftedMePUcorrt", &eSmearedShiftedMePUcorrt, &b_eSmearedShiftedMetPUcorr);
  tree_->SetBranchAddress("phiSmearedShiftedMetPUcorr", &phiSmearedShiftedMetPUcorr, &b_phiSmearedShiftedMetPUcorr);
  tree_->SetBranchAddress("sCorrMet", &sCorrMet, &b_sCorrMet);
  tree_->SetBranchAddress("eCorrMet", &eCorrMet, &b_eCorrMet);
  tree_->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
  tree_->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
  tree_->SetBranchAddress("smuCorrMet", &smuCorrMet, &b_smuCorrMet);
  tree_->SetBranchAddress("emuCorrMet", &emuCorrMet, &b_emuCorrMet);
  tree_->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
  tree_->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
  tree_->SetBranchAddress("sNoHFMet", &sNoHFMet, &b_sNoHFMet);
  tree_->SetBranchAddress("eNoHFMet", &eNoHFMet, &b_eNoHFMet);
  tree_->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
  tree_->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
  tree_->SetBranchAddress("stcMet", &stcMet, &b_stcMet);
  tree_->SetBranchAddress("etcMet", &etcMet, &b_etcMet);
  tree_->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
  tree_->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
  tree_->SetBranchAddress("sglobalPfMet", &sglobalPfMet, &b_sglobalPfMet);
  tree_->SetBranchAddress("eglobalPfMet", &eglobalPfMet, &b_eglobalPfMet);
  tree_->SetBranchAddress("phiglobalPfMet", &phiglobalPfMet, &b_phiglobalPfMet);
  tree_->SetBranchAddress("signifglobalPfMet", &signifglobalPfMet, &b_signifglobalPfMet);
  tree_->SetBranchAddress("scentralPfMet", &scentralPfMet, &b_scentralPfMet);
  tree_->SetBranchAddress("ecentralPfMet", &ecentralPfMet, &b_ecentralPfMet);
  tree_->SetBranchAddress("phicentralPfMet", &phicentralPfMet, &b_phicentralPfMet);
  tree_->SetBranchAddress("signifcentralPfMet", &signifcentralPfMet, &b_signifcentralPfMet);
  tree_->SetBranchAddress("eassocPfMet", &eassocPfMet, &b_eassocPfMet);
  tree_->SetBranchAddress("phiassocPfMet", &phiassocPfMet, &b_phiassocPfMet);
  tree_->SetBranchAddress("signifassocPfMet", &signifassocPfMet, &b_signifassocPfMet);
  tree_->SetBranchAddress("eassocOtherVtxPfMet", &eassocOtherVtxPfMet, &b_eassocOtherVtxPfMet);
  tree_->SetBranchAddress("phiassocOtherVtxPfMet", &phiassocOtherVtxPfMet, &b_phiassocOtherVtxPfMet);
  tree_->SetBranchAddress("signifassocOtherVtxPfMet", &signifassocOtherVtxPfMet, &b_signifassocOtherVtxPfMet);
  tree_->SetBranchAddress("etrkPfMet", &etrkPfMet, &b_etrkPfMet);
  tree_->SetBranchAddress("phitrkPfMet", &phitrkPfMet, &b_phitrkPfMet);
  tree_->SetBranchAddress("signiftrkPfMet", &signiftrkPfMet, &b_signiftrkPfMet);
  tree_->SetBranchAddress("ecleanPfMet", &ecleanPfMet, &b_ecleanPfMet);
  tree_->SetBranchAddress("phicleanPfMet", &phicleanPfMet, &b_phicleanPfMet);
  tree_->SetBranchAddress("signifcleanPfMet", &signifcleanPfMet, &b_signifcleanPfMet);
  tree_->SetBranchAddress("ecleanedSaclayPfMet", &ecleanedSaclayPfMet, &b_ecleanedSaclayPfMet);
  tree_->SetBranchAddress("phicleanedSaclayPfMet", &phicleanedSaclayPfMet, &b_phicleanedSaclayPfMet);
  tree_->SetBranchAddress("signifcleanedSaclayPfMet", &signifcleanedSaclayPfMet, &b_signifcleanedSaclayPfMet);
  tree_->SetBranchAddress("eminTypeICleanSaclayPfMet", &eminTypeICleanSaclayPfMet, &b_eminTypeICleanSaclayPfMet);
  tree_->SetBranchAddress("phiminTypeICleanSaclayPfMet", &phiminTypeICleanSaclayPfMet, &b_phiminTypeICleanSaclayPfMet);
  tree_->SetBranchAddress("signifminTypeICleanSaclayPfMet", &signifminTypeICleanSaclayPfMet, &b_signifminTypeICleanSaclayPfMet);
  tree_->SetBranchAddress("globalPfSums", &globalPfSums, &b_globalPfSums);
  tree_->SetBranchAddress("spfMet", &spfMet, &b_spfMet);
  tree_->SetBranchAddress("epfMet", &epfMet, &b_epfMet);
  tree_->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
  tree_->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
  tree_->SetBranchAddress("spfMetType1", &spfMetType1, &b_spfMetType1);
  tree_->SetBranchAddress("epfMetType1", &epfMetType1, &b_epfMetType1);
  tree_->SetBranchAddress("phipfMetType1", &phipfMetType1, &b_phipfMetType1);
  tree_->SetBranchAddress("signifpfMetType1", &signifpfMetType1, &b_signifpfMetType1);
  tree_->SetBranchAddress("sMetGen", &sMetGen, &b_sMetGen);
  tree_->SetBranchAddress("eMetGen", &eMetGen, &b_eMetGen);
  tree_->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
  tree_->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
  tree_->SetBranchAddress("sMetGen2", &sMetGen2, &b_sMetGen2);
  tree_->SetBranchAddress("eMetGen2", &eMetGen2, &b_eMetGen2);
  tree_->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
  tree_->SetBranchAddress("npu", &npu, &b_npu);
  tree_->SetBranchAddress("NtotEvents", &NtotEvents, &b_NtotEvents);
  tree_->SetBranchAddress("xsection", &xsection, &b_xsection);
  tree_->SetBranchAddress("EquivLumi", &EquivLumi, &b_EquivLumi);
  tree_->SetBranchAddress("SampleID", &SampleID, &b_SampleID);
  tree_->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  tree_->SetBranchAddress("pt_weight", &pt_weight, &b_pt_weight);
  tree_->SetBranchAddress("gen_custom_processId", &gen_custom_processId, &b_gen_custom_processId);
  tree_->SetBranchAddress("gen_pt_gamma1", &gen_pt_gamma1, &b_gen_pt_gamma1);
  tree_->SetBranchAddress("gen_pt_gamma2", &gen_pt_gamma2, &b_gen_pt_gamma2);
  tree_->SetBranchAddress("gen_eta_gamma1", &gen_eta_gamma1, &b_gen_eta_gamma1);
  tree_->SetBranchAddress("gen_eta_gamma2", &gen_eta_gamma2, &b_gen_eta_gamma2);
  tree_->SetBranchAddress("gen_phi_gamma1", &gen_phi_gamma1, &b_gen_phi_gamma1);
  tree_->SetBranchAddress("gen_phi_gamma2", &gen_phi_gamma2, &b_gen_phi_gamma2);
  tree_->SetBranchAddress("gen_mass_diphoton", &gen_mass_diphoton, &b_gen_mass_diphoton);
  tree_->SetBranchAddress("gen_pt_diphoton", &gen_pt_diphoton, &b_gen_pt_diphoton);
  tree_->SetBranchAddress("gen_eta_diphoton", &gen_eta_diphoton, &b_gen_eta_diphoton);
  tree_->SetBranchAddress("gen_phi_diphoton", &gen_phi_diphoton, &b_gen_phi_diphoton);
  tree_->SetBranchAddress("gen_mass_dijet", &gen_mass_dijet, &b_gen_mass_dijet);
  tree_->SetBranchAddress("gen_pt_dijet", &gen_pt_dijet, &b_gen_pt_dijet);
  tree_->SetBranchAddress("gen_eta_dijet", &gen_eta_dijet, &b_gen_eta_dijet);
  tree_->SetBranchAddress("gen_phi_dijet", &gen_phi_dijet, &b_gen_phi_dijet);
  tree_->SetBranchAddress("gen_zeppenfeld", &gen_zeppenfeld, &b_gen_zeppenfeld);
  tree_->SetBranchAddress("gen_pt_lep1", &gen_pt_lep1, &b_gen_pt_lep1);
  tree_->SetBranchAddress("gen_pt_lep2", &gen_pt_lep2, &b_gen_pt_lep2);
  tree_->SetBranchAddress("gen_eta_lep1", &gen_eta_lep1, &b_gen_eta_lep1);
  tree_->SetBranchAddress("gen_eta_lep2", &gen_eta_lep2, &b_gen_eta_lep2);
  tree_->SetBranchAddress("gen_phi_lep1", &gen_phi_lep1, &b_gen_phi_lep1);
  tree_->SetBranchAddress("gen_phi_lep2", &gen_phi_lep2, &b_gen_phi_lep2);
  tree_->SetBranchAddress("gen_pid_lep1", &gen_pid_lep1, &b_gen_pid_lep1);
  tree_->SetBranchAddress("gen_pid_lep2", &gen_pid_lep2, &b_gen_pid_lep2);
  tree_->SetBranchAddress("ptele1", &ptele1, &b_ptele1);
  tree_->SetBranchAddress("ptele2", &ptele2, &b_ptele2);
  tree_->SetBranchAddress("etaele1", &etaele1, &b_etaele1);
  tree_->SetBranchAddress("etaele2", &etaele2, &b_etaele2);
  tree_->SetBranchAddress("phiele1", &phiele1, &b_phiele1);
  tree_->SetBranchAddress("phiele2", &phiele2, &b_phiele2);
  tree_->SetBranchAddress("eneele1", &eneele1, &b_eneele1);
  tree_->SetBranchAddress("eneele2", &eneele2, &b_eneele2);
  tree_->SetBranchAddress("sIeIeele1", &sIeIeele1, &b_sIeIeele1);
  tree_->SetBranchAddress("sIeIeele2", &sIeIeele2, &b_sIeIeele2);
  tree_->SetBranchAddress("dphiele1", &dphiele1, &b_dphiele1);
  tree_->SetBranchAddress("dphiele2", &dphiele2, &b_dphiele2);
  tree_->SetBranchAddress("detaele1", &detaele1, &b_detaele1);
  tree_->SetBranchAddress("detaele2", &detaele2, &b_detaele2);
  tree_->SetBranchAddress("hoeele1", &hoeele1, &b_hoeele1);
  tree_->SetBranchAddress("hoeele2", &hoeele2, &b_hoeele2);
  tree_->SetBranchAddress("mhitsele1", &mhitsele1, &b_mhitsele1);
  tree_->SetBranchAddress("mhitsele2", &mhitsele2, &b_mhitsele2);
  tree_->SetBranchAddress("d0ele1", &d0ele1, &b_d0ele1);
  tree_->SetBranchAddress("d0ele2", &d0ele2, &b_d0ele2);
  tree_->SetBranchAddress("dzele1", &dzele1, &b_dzele1);
  tree_->SetBranchAddress("dzele2", &dzele2, &b_dzele2);
  tree_->SetBranchAddress("invMassele1g1", &invMassele1g1, &b_invMassele1g1);
  tree_->SetBranchAddress("invMassele1g2", &invMassele1g2, &b_invMassele1g2);
  tree_->SetBranchAddress("invMassele2g1", &invMassele2g1, &b_invMassele2g1);
  tree_->SetBranchAddress("invMassele2g2", &invMassele2g2, &b_invMassele2g2);
  tree_->SetBranchAddress("oEmoPele1", &oEmoPele1, &b_oEmoPele1);
  tree_->SetBranchAddress("oEmoPele2", &oEmoPele2, &b_oEmoPele2);
  tree_->SetBranchAddress("mvanotrigele1", &mvanotrigele1, &b_mvanotrigele1);
  tree_->SetBranchAddress("mvanotrigele2", &mvanotrigele2, &b_mvanotrigele2);
  tree_->SetBranchAddress("mvatrigele1", &mvatrigele1, &b_mvatrigele1);
  tree_->SetBranchAddress("mvatrigele2", &mvatrigele2, &b_mvatrigele2);
  tree_->SetBranchAddress("matchconvele1", &matchconvele1, &b_matchconvele1);
  tree_->SetBranchAddress("matchconvele2", &matchconvele2, &b_matchconvele2);
  tree_->SetBranchAddress("chHadIso03ele1", &chHadIso03ele1, &b_chHadIso03ele1);
  tree_->SetBranchAddress("chHadIso03ele2", &chHadIso03ele2, &b_chHadIso03ele2);
  tree_->SetBranchAddress("nHadIso03ele1", &nHadIso03ele1, &b_nHadIso03ele1);
  tree_->SetBranchAddress("nHadIso03ele2", &nHadIso03ele2, &b_nHadIso03ele2);
  tree_->SetBranchAddress("photIso03ele1", &photIso03ele1, &b_photIso03ele1);
  tree_->SetBranchAddress("photIso03ele2", &photIso03ele2, &b_photIso03ele2);
  tree_->SetBranchAddress("pteleloose1", &pteleloose1, &b_pteleloose1);
  tree_->SetBranchAddress("pteleloose2", &pteleloose2, &b_pteleloose2);
  tree_->SetBranchAddress("etaeleloose1", &etaeleloose1, &b_etaeleloose1);
  tree_->SetBranchAddress("etaeleloose2", &etaeleloose2, &b_etaeleloose2);
  tree_->SetBranchAddress("phieleloose1", &phieleloose1, &b_phieleloose1);
  tree_->SetBranchAddress("phieleloose2", &phieleloose2, &b_phieleloose2);
  tree_->SetBranchAddress("eneeleloose1", &eneeleloose1, &b_eneeleloose1);
  tree_->SetBranchAddress("eneeleloose2", &eneeleloose2, &b_eneeleloose2);
  tree_->SetBranchAddress("sIeIeeleloose1", &sIeIeeleloose1, &b_sIeIeeleloose1);
  tree_->SetBranchAddress("sIeIeeleloose2", &sIeIeeleloose2, &b_sIeIeeleloose2);
  tree_->SetBranchAddress("dphieleloose1", &dphieleloose1, &b_dphieleloose1);
  tree_->SetBranchAddress("dphieleloose2", &dphieleloose2, &b_dphieleloose2);
  tree_->SetBranchAddress("detaeleloose1", &detaeleloose1, &b_detaeleloose1);
  tree_->SetBranchAddress("detaeleloose2", &detaeleloose2, &b_detaeleloose2);
  tree_->SetBranchAddress("hoeeleloose1", &hoeeleloose1, &b_hoeeleloose1);
  tree_->SetBranchAddress("hoeeleloose2", &hoeeleloose2, &b_hoeeleloose2);
  tree_->SetBranchAddress("mhitseleloose1", &mhitseleloose1, &b_mhitseleloose1);
  tree_->SetBranchAddress("mhitseleloose2", &mhitseleloose2, &b_mhitseleloose2);
  tree_->SetBranchAddress("d0eleloose1", &d0eleloose1, &b_d0eleloose1);
  tree_->SetBranchAddress("d0eleloose2", &d0eleloose2, &b_d0eleloose2);
  tree_->SetBranchAddress("dzeleloose1", &dzeleloose1, &b_dzeleloose1);
  tree_->SetBranchAddress("dzeleloose2", &dzeleloose2, &b_dzeleloose2);
  tree_->SetBranchAddress("invMasseleloose1g1", &invMasseleloose1g1, &b_invMasseleloose1g1);
  tree_->SetBranchAddress("invMasseleloose1g2", &invMasseleloose1g2, &b_invMasseleloose1g2);
  tree_->SetBranchAddress("invMasseleloose2g1", &invMasseleloose2g1, &b_invMasseleloose2g1);
  tree_->SetBranchAddress("invMasseleloose2g2", &invMasseleloose2g2, &b_invMasseleloose2g2);
  tree_->SetBranchAddress("oEmoPeleloose1", &oEmoPeleloose1, &b_oEmoPeleloose1);
  tree_->SetBranchAddress("oEmoPeleloose2", &oEmoPeleloose2, &b_oEmoPeleloose2);
  tree_->SetBranchAddress("mvanotrigeleloose1", &mvanotrigeleloose1, &b_mvanotrigeleloose1);
  tree_->SetBranchAddress("mvanotrigeleloose2", &mvanotrigeleloose2, &b_mvanotrigeleloose2);
  tree_->SetBranchAddress("mvatrigeleloose1", &mvatrigeleloose1, &b_mvatrigeleloose1);
  tree_->SetBranchAddress("mvatrigeleloose2", &mvatrigeleloose2, &b_mvatrigeleloose2);
  tree_->SetBranchAddress("matchconveleloose1", &matchconveleloose1, &b_matchconveleloose1);
  tree_->SetBranchAddress("matchconveleloose2", &matchconveleloose2, &b_matchconveleloose2);
  tree_->SetBranchAddress("chHadIso03eleloose1", &chHadIso03eleloose1, &b_chHadIso03eleloose1);
  tree_->SetBranchAddress("chHadIso03eleloose2", &chHadIso03eleloose2, &b_chHadIso03eleloose2);
  tree_->SetBranchAddress("nHadIso03eleloose1", &nHadIso03eleloose1, &b_nHadIso03eleloose1);
  tree_->SetBranchAddress("nHadIso03eleloose2", &nHadIso03eleloose2, &b_nHadIso03eleloose2);
  tree_->SetBranchAddress("photIso03eleloose1", &photIso03eleloose1, &b_photIso03eleloose1);
  tree_->SetBranchAddress("photIso03eleloose2", &photIso03eleloose2, &b_photIso03eleloose2);
  tree_->SetBranchAddress("ptelenontr801", &ptelenontr801, &b_ptelenontr801);
  tree_->SetBranchAddress("ptelenontr802", &ptelenontr802, &b_ptelenontr802);
  tree_->SetBranchAddress("etaelenontr801", &etaelenontr801, &b_etaelenontr801);
  tree_->SetBranchAddress("etaelenontr802", &etaelenontr802, &b_etaelenontr802);
  tree_->SetBranchAddress("phielenontr801", &phielenontr801, &b_phielenontr801);
  tree_->SetBranchAddress("phielenontr802", &phielenontr802, &b_phielenontr802);
  tree_->SetBranchAddress("eneelenontr801", &eneelenontr801, &b_eneelenontr801);
  tree_->SetBranchAddress("eneelenontr802", &eneelenontr802, &b_eneelenontr802);
  tree_->SetBranchAddress("ptelenontr901", &ptelenontr901, &b_ptelenontr901);
  tree_->SetBranchAddress("ptelenontr902", &ptelenontr902, &b_ptelenontr902);
  tree_->SetBranchAddress("etaelenontr901", &etaelenontr901, &b_etaelenontr901);
  tree_->SetBranchAddress("etaelenontr902", &etaelenontr902, &b_etaelenontr902);
  tree_->SetBranchAddress("phielenontr901", &phielenontr901, &b_phielenontr901);
  tree_->SetBranchAddress("phielenontr902", &phielenontr902, &b_phielenontr902);
  tree_->SetBranchAddress("eneelenontr901", &eneelenontr901, &b_eneelenontr901);
  tree_->SetBranchAddress("eneelenontr902", &eneelenontr902, &b_eneelenontr902);
  tree_->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
  tree_->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
  tree_->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
  tree_->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
  tree_->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
  tree_->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
  tree_->SetBranchAddress("enemu1", &enemu1, &b_enemu1);
  tree_->SetBranchAddress("enemu2", &enemu2, &b_enemu2);
  tree_->SetBranchAddress("pixhitsmu1", &pixhitsmu1, &b_pixhitsmu1);
  tree_->SetBranchAddress("pixhitsmu2", &pixhitsmu2, &b_pixhitsmu2);
  tree_->SetBranchAddress("trkhitsmu1", &trkhitsmu1, &b_trkhitsmu1);
  tree_->SetBranchAddress("trkhitsmu2", &trkhitsmu2, &b_trkhitsmu2);
  tree_->SetBranchAddress("hitsmu1", &hitsmu1, &b_hitsmu1);
  tree_->SetBranchAddress("hitsmu2", &hitsmu2, &b_hitsmu2);
  tree_->SetBranchAddress("chi2mu1", &chi2mu1, &b_chi2mu1);
  tree_->SetBranchAddress("chi2mu2", &chi2mu2, &b_chi2mu2);
  tree_->SetBranchAddress("matchmu1", &matchmu1, &b_matchmu1);
  tree_->SetBranchAddress("matchmu2", &matchmu2, &b_matchmu2);
  tree_->SetBranchAddress("d0mu1", &d0mu1, &b_d0mu1);
  tree_->SetBranchAddress("d0mu2", &d0mu2, &b_d0mu2);
  tree_->SetBranchAddress("dzmu1", &dzmu1, &b_dzmu1);
  tree_->SetBranchAddress("dzmu2", &dzmu2, &b_dzmu2);
  tree_->SetBranchAddress("chHadmu1", &chHadmu1, &b_chHadmu1);
  tree_->SetBranchAddress("chHadmu2", &chHadmu2, &b_chHadmu2);
  tree_->SetBranchAddress("nHadmu1", &nHadmu1, &b_nHadmu1);
  tree_->SetBranchAddress("nHadmu2", &nHadmu2, &b_nHadmu2);
  tree_->SetBranchAddress("photmu1", &photmu1, &b_photmu1);
  tree_->SetBranchAddress("photmu2", &photmu2, &b_photmu2);
  tree_->SetBranchAddress("puptmu1", &puptmu1, &b_puptmu1);
  tree_->SetBranchAddress("puptmu2", &puptmu2, &b_puptmu2);
  tree_->SetBranchAddress("ptmuloose1", &ptmuloose1, &b_ptmuloose1);
  tree_->SetBranchAddress("ptmuloose2", &ptmuloose2, &b_ptmuloose2);
  tree_->SetBranchAddress("etamuloose1", &etamuloose1, &b_etamuloose1);
  tree_->SetBranchAddress("etamuloose2", &etamuloose2, &b_etamuloose2);
  tree_->SetBranchAddress("phimuloose1", &phimuloose1, &b_phimuloose1);
  tree_->SetBranchAddress("phimuloose2", &phimuloose2, &b_phimuloose2);
  tree_->SetBranchAddress("enemuloose1", &enemuloose1, &b_enemuloose1);
  tree_->SetBranchAddress("enemuloose2", &enemuloose2, &b_enemuloose2);
  tree_->SetBranchAddress("pixhitsmuloose1", &pixhitsmuloose1, &b_pixhitsmuloose1);
  tree_->SetBranchAddress("pixhitsmuloose2", &pixhitsmuloose2, &b_pixhitsmuloose2);
  tree_->SetBranchAddress("trkhitsmuloose1", &trkhitsmuloose1, &b_trkhitsmuloose1);
  tree_->SetBranchAddress("trkhitsmuloose2", &trkhitsmuloose2, &b_trkhitsmuloose2);
  tree_->SetBranchAddress("hitsmuloose1", &hitsmuloose1, &b_hitsmuloose1);
  tree_->SetBranchAddress("hitsmuloose2", &hitsmuloose2, &b_hitsmuloose2);
  tree_->SetBranchAddress("chi2muloose1", &chi2muloose1, &b_chi2muloose1);
  tree_->SetBranchAddress("chi2muloose2", &chi2muloose2, &b_chi2muloose2);
  tree_->SetBranchAddress("matchmuloose1", &matchmuloose1, &b_matchmuloose1);
  tree_->SetBranchAddress("matchmuloose2", &matchmuloose2, &b_matchmuloose2);
  tree_->SetBranchAddress("d0muloose1", &d0muloose1, &b_d0muloose1);
  tree_->SetBranchAddress("d0muloose2", &d0muloose2, &b_d0muloose2);
  tree_->SetBranchAddress("dzmuloose1", &dzmuloose1, &b_dzmuloose1);
  tree_->SetBranchAddress("dzmuloose2", &dzmuloose2, &b_dzmuloose2);
  tree_->SetBranchAddress("ptmuvloose1", &ptmuvloose1, &b_ptmuvloose1);
  tree_->SetBranchAddress("ptmuvloose2", &ptmuvloose2, &b_ptmuvloose2);
  tree_->SetBranchAddress("etamuvloose1", &etamuvloose1, &b_etamuvloose1);
  tree_->SetBranchAddress("etamuvloose2", &etamuvloose2, &b_etamuvloose2);
  tree_->SetBranchAddress("phimuvloose1", &phimuvloose1, &b_phimuvloose1);
  tree_->SetBranchAddress("phimuvloose2", &phimuvloose2, &b_phimuvloose2);
  tree_->SetBranchAddress("enemuvloose1", &enemuvloose1, &b_enemuvloose1);
  tree_->SetBranchAddress("enemuvloose2", &enemuvloose2, &b_enemuvloose2);
  tree_->SetBranchAddress("pixhitsmuvloose1", &pixhitsmuvloose1, &b_pixhitsmuvloose1);
  tree_->SetBranchAddress("pixhitsmuvloose2", &pixhitsmuvloose2, &b_pixhitsmuvloose2);
  tree_->SetBranchAddress("trkhitsmuvloose1", &trkhitsmuvloose1, &b_trkhitsmuvloose1);
  tree_->SetBranchAddress("trkhitsmuvloose2", &trkhitsmuvloose2, &b_trkhitsmuvloose2);
  tree_->SetBranchAddress("hitsmuvloose1", &hitsmuvloose1, &b_hitsmuvloose1);
  tree_->SetBranchAddress("hitsmuvloose2", &hitsmuvloose2, &b_hitsmuvloose2);
  tree_->SetBranchAddress("chi2muvloose1", &chi2muvloose1, &b_chi2muvloose1);
  tree_->SetBranchAddress("chi2muvloose2", &chi2muvloose2, &b_chi2muvloose2);
  tree_->SetBranchAddress("matchmuvloose1", &matchmuvloose1, &b_matchmuvloose1);
  tree_->SetBranchAddress("matchmuvloose2", &matchmuvloose2, &b_matchmuvloose2);
  tree_->SetBranchAddress("d0muvloose1", &d0muvloose1, &b_d0muvloose1);
  tree_->SetBranchAddress("d0muvloose2", &d0muvloose2, &b_d0muvloose2);
  tree_->SetBranchAddress("dzmuvloose1", &dzmuvloose1, &b_dzmuvloose1);
  tree_->SetBranchAddress("dzmuvloose2", &dzmuvloose2, &b_dzmuvloose2);
  tree_->SetBranchAddress("hasPassedSinglePhot", &hasPassedSinglePhot, &b_hasPassedSinglePhot);
  tree_->SetBranchAddress("hasPassedDoublePhot", &hasPassedDoublePhot, &b_hasPassedDoublePhot);
  tree_->SetBranchAddress("chHadmuloose1", &chHadmuloose1, &b_chHadmuloose1);
  tree_->SetBranchAddress("chHadmuloose2", &chHadmuloose2, &b_chHadmuloose2);
  tree_->SetBranchAddress("nHadmuloose1", &nHadmuloose1, &b_nHadmuloose1);
  tree_->SetBranchAddress("nHadmuloose2", &nHadmuloose2, &b_nHadmuloose2);
  tree_->SetBranchAddress("photmuloose1", &photmuloose1, &b_photmuloose1);
  tree_->SetBranchAddress("photmuloose2", &photmuloose2, &b_photmuloose2);
  tree_->SetBranchAddress("puptmuloose1", &puptmuloose1, &b_puptmuloose1);
  tree_->SetBranchAddress("puptmuloose2", &puptmuloose2, &b_puptmuloose2);
  tree_->SetBranchAddress("chHadmuvloose1", &chHadmuvloose1, &b_chHadmuvloose1);
  tree_->SetBranchAddress("chHadmuvloose2", &chHadmuvloose2, &b_chHadmuvloose2);
  tree_->SetBranchAddress("nHadmuvloose1", &nHadmuvloose1, &b_nHadmuvloose1);
  tree_->SetBranchAddress("nHadmuvloose2", &nHadmuvloose2, &b_nHadmuvloose2);
  tree_->SetBranchAddress("photmuvloose1", &photmuvloose1, &b_photmuvloose1);
  tree_->SetBranchAddress("photmuvloose2", &photmuvloose2, &b_photmuvloose2);
  tree_->SetBranchAddress("puptmuvloose1", &puptmuvloose1, &b_puptmuvloose1);
  tree_->SetBranchAddress("puptmuvloose2", &puptmuvloose2, &b_puptmuvloose2);
}

void bestJets::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;
  
  // default values                                                                                                     
  dopureeventWeight_ = true;
  cicselection = 4;

  ptphot1cut = 40.;
  ptphot2cut = 25.;

  ptjetacccut  = 25.;
  etajetacccut = 2.5;  
}

std::pair<int,int> bestJets::myTwoHighestPtJet(std::vector<int> jets) {

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

std::pair<int,int> bestJets::myTwoHighestBtagJet(std::vector<int> jets) {

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

int bestJets::myHighestPtJet(std::vector<int> jets) {

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

int bestJets::myHighestBtagJet(std::vector<int> jets) {

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

bool bestJets::passCutBasedJetId(int jet) {

  bool isGood = true;

  if ( fabs(etajet[jet]) < 2.5 ) {
    if ( betastarjet[jet] > 0.2 * log( nvtx - 0.64 ) ) isGood = false;
    if ( rmsjet[jet] > 0.06)                           isGood = false;
  } else if (fabs(etajet[jet]) < 3.){
    if ( rmsjet[jet] > 0.05)  isGood =false;
  } else {
    if ( rmsjet[jet] > 0.055) isGood =false;
  }
  
  return isGood;
}
