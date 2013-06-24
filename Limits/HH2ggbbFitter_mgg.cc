/** \macro H2GGFitter.cc
 *
 * The analysis root trees produced in a simple format 
 *
 *     TFile file(filename,"RECREATE", "X->jj input tree for unbinned maximum-likelihood fit");
 *     TTree* outTree  = new TTree("XTojj","X->jj input tree for unbinned maximum-likelihood fit");
 *     Float_t mass;
 *     Int_t CAT3;
 *     Float_t weight;
 *
 *     outTree->Branch("mass",&mass,"mass/F");
 *     outTree->Branch("weight",&weight,"weight/F");
 *     outTree->Branch("CAT4",&CAT4,"CAT4/I");
 *     {
 *       .............
 *       outTree->Fill();
 *     }
 *
 *     file.Write();
 *     file.Close();
 *     delete outTree;
 *
 * are used as input files. They have to be produced for 
 * data and Monte Carlo signal and background data sets 
 * after all analysis selections to be applied.  *
 */
// Loading:  .L HH2ggbbFitter.cc
// Running:  runfits("hgg120-shapes-combined-Unbinned.root")  
//                

using namespace RooFit;
using namespace RooStats ;

static const Int_t NCAT = 2;  // chiara
Double_t MMIN = 100;        
Double_t MMAX = 180;       
std::string filePOSTfix="";
//double signalScaler=0.005;  // chiara
double signalScaler=1.00;

void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*);
void SigModelFit(RooWorkspace*, Float_t);
RooFitResult*  BkgModelFitBernstein(RooWorkspace*, Bool_t);
RooFitResult*  BkgModelFitExpo(RooWorkspace*, Bool_t);
void MakePlots(RooWorkspace*, Float_t, RooFitResult*);
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeBkgWS(RooWorkspace* w, const char* filename);
void SetConstantParams(const RooArgSet* params);

RooArgSet* defineVariables() {

  // define variables of the input ntuple
  RooRealVar* massggnewvtx  = new RooRealVar("massggnewvtx", "M(gg)",MMIN,MMAX,"GeV");
  RooRealVar* mjj = new RooRealVar("mjj", "M(jj)",10,300,"GeV");
  RooRealVar* mggjj = new RooRealVar("mggjj", "M(ggjj)",10,1500,"GeV");
  RooRealVar* runptPhot1 = new RooRealVar("runptPhot1","runptPhot1",0,1000,"GeV");
  RooRealVar* runptCorrJet1 = new RooRealVar("runptCorrJet1","runptCorrJet1",0,1000,"GeV");
  RooRealVar* absCosThetaStar = new RooRealVar("absCosThetaStar","absCosThetaStar",0,1,"");
  // RooRealVar* HT_jet = new RooRealVar("HT_jet","HT_jet",0,1000,"GeV");
  // RooRealVar* njets = new RooRealVar("njets", "# jets",-10,10,"");
  RooRealVar* weight = new RooRealVar("weight","Reweightings",0,10,"");
  RooRealVar* btagCategory = new RooRealVar("btagCategory","event category",0.9,2.1,"") ;
  
  RooArgSet* ntplVars = new RooArgSet(*massggnewvtx, *mjj, *mggjj, *runptPhot1, *runptCorrJet1, *absCosThetaStar, *btagCategory, *weight);
  
  return ntplVars;
}

void runfits(const Float_t mass=2000, Bool_t dobands = false) {

  //******************************************************************//
  //  Running mode  corresponds to the following cases
  //         - full run set:
  //         - create signal and background data sets 
  //         - make and fit signal and background  models 
  //         - write signal and background workspaces in root files
  //         - write data card
  //*******************************************************************//

  TString fileBaseName(TString::Format("Xhh.mX%.1f", mass));    
  TString fileBkgName(TString::Format("Xhh.inputbkg_8TeV", mass));
  
  TString card_name("Xhhgg_mgg_models_Bkg_8TeV_test.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  RooFitResult* fitresults;

  w->var("massggnewvtx")->setMin(MMIN);
  w->var("massggnewvtx")->setMax(MMAX);

  // Add data to the workspace
  cout << endl; cout << "Now AddSigData" << endl;
  AddSigData(w, mass);   

  cout << endl; cout << "Now AddBkgData" << endl;
  AddBkgData(w);         
  
  // Add the signal and background models to the workspace.
  // Inside this function you will find a discription our model.
  cout << endl; cout << "Now SigModelFit" << endl;
  SigModelFit(w, mass);      

  cout << endl; cout << "Now BkgModelFit" << endl;    
  // fitresults = BkgModelFitBernstein(w, dobands);   
  fitresults = BkgModelFitExpo(w, dobands);   
   
  // Make statistical treatment. Setup the limit on X production
  cout << endl; cout << "Now make signal workspace" << endl;
  MakeSigWS(w, fileBaseName+"_8TeV");     
  cout << endl; cout << "Now make background workspace" << endl;
  MakeBkgWS(w, fileBkgName);              

  cout << endl; cout << "Now prepare datacards" << endl;
  int ncat = NCAT;
  for (int c=0; c<ncat; c++) MakeDataCard_1Channel(w, fileBaseName, fileBkgName, c);   

  // Make plots for data and fit results
  cout << endl; cout << "Preparing final plots" << endl;
  MakePlots(w, mass, fitresults);   
   
  return;
}

// Signal Data Set
void AddSigData(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  TString inDir = "../finalizedTrees_Radion_presel/";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  int iMass = abs(mass);       
  TFile sigFile1(inDir+TString(Form("Radion_Radion_M-%d_madgraph_default_CSV.root", iMass)));  
  TTree* sigTree1 = (TTree*) sigFile1.Get("myTrees");

  // common preselection cut
  // TString mainCut("massggnewvtx>=100 && massggnewvtx<=180 && mjj>98.5 && mjj<151.5 && mggjj>270 && mggjj<330");   // chiara
  // TString mainCut("massggnewvtx>=100 && massggnewvtx<=180 && ( (btagCategory==1 && mjj<150 && mjj>90 && mggjj>260 && mggjj<335) || (btagCategory==2 && mjj<140 && mjj>95 && mggjj<320 && mggjj>255) )");   // chiara
  TString mainCut("massggnewvtx>=100 && massggnewvtx<=180 && ( (btagCategory==1 && mjj<150 && mjj>90 && mggjj>260 && mggjj<335 && runptCorrJet1>50 && runptPhot1>55 && absCosThetaStar<0.8) || (btagCategory==2 && mjj<140 && mjj>95 && mggjj<320 && mggjj>255 && runptPhot1>50) )");   // chiara

  // Create signal dataset composed with different productions, the weight is already applied in our ntuples
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"weight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
  RooDataSet* signal[NCAT];
  for (int c=0; c<ncat; ++c) {

    // 0) chiara: 1cat only
    // signal[c] =  (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut);   //chiara, for 1 cat only

    // a) chiara: EBEB and !EBEB
    // if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( abs(etaPhot1)<1.5 && abs(etaPhot2)<1.5 )"));
    // if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( abs(etaPhot1)>1.5 || abs(etaPhot2)>1.5 )"));

    // b) chiara: altoR9 and !altoR9
    // if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1>0.94 && r9Phot2>0.94 )"));
    // if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1<0.94 || r9Phot2<0.94 )"));

    // c) chiara: altoR9, EBEB and NOT
    // if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1>0.94 && r9Phot2>0.94 && abs(etaPhot1)<1.5 && abs(etaPhot2)<1.5)"));
    // if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1<0.94 || r9Phot2<0.94 || abs(etaPhot1)>1.5 || abs(etaPhot2)>1.5)"));
    
    // d) chiara: 1-btag, 2-btag 
    if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==1"));
    if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==2"));

    // e) chiara: 2-btag, 1-btag high R9, 1-btag low R9
    // if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==1 && r9Phot1>0.94 && r9Phot2>0.94"));
    // if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==1 && (r9Phot1<0.94 || r9Phot2<0.94)"));
    // if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==2"));
 

    w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
    
    cout << "cat " << c << ", signal[c]: " << endl;
    signal[c]->Print("v");
    cout << "---- for category " << c << ", nX for signal[c]:  " << signal[c]->sumEntries() << endl; 
    cout << endl;
  }

  // Create full weighted signal data set without categorization
  RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut);
  w->import(*signalAll, Rename("SigWeight"));
  cout << "now signalAll" << endl;
  signalAll->Print("v");
  cout << "---- nX for signalAll:  " << signalAll->sumEntries() << endl; 
  cout << endl;
}

// Data dataset
void AddBkgData(RooWorkspace* w) {

  // initializations
  Int_t ncat = NCAT;
  Float_t minMassFit(MMIN),maxMassFit(MMAX); 

  // retrieve the data tree; no common preselection cut applied yet; 
  TString inDir = "../finalizedTrees_Radion_presel/";
  TFile dataFile(inDir+"Radion_Data2012_default_CSV.root");
  TTree* dataTree = (TTree*) dataFile.Get("myTrees");

  // Variables
  RooArgSet* ntplVars = defineVariables();

  // common preselection cut
  // TString mainCut("massggnewvtx>=100 && massggnewvtx<=180 && mjj>98.5 && mjj<151.5 && mggjj>270 && mggjj<330");  // chiara
  // TString mainCut("massggnewvtx>=100 && massggnewvtx<=180 && ( (btagCategory==1 && mjj<150 && mjj>90 && mggjj>260 && mggjj<335) || (btagCategory==2 && mjj<140 && mjj>95 && mggjj<320 && mggjj>255) )");   // chiara
  TString mainCut("massggnewvtx>=100 && massggnewvtx<=180 && ( (btagCategory==1 && mjj<150 && mjj>90 && mggjj>260 && mggjj<335 && runptCorrJet1>50 && runptPhot1>55 && absCosThetaStar<0.8) || (btagCategory==2 && mjj<140 && mjj>95 && mggjj<320 && mggjj>255 && runptPhot1>50) )");   // chiara

  // Create dataset
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,mainCut,"weight");
  cout << endl;
  cout << "Data, everything: " << endl;
  Data.Print("v");
  cout << "---- nX:  " << Data.sumEntries() << endl;
  cout << endl;

  // split into NCAT categories
  RooDataSet* dataToFit[NCAT];  
  for (int c=0; c<ncat; ++c) {
    int theCat = c+1;

    // 0) chiara: 1 cat only
    // dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut);  

    // a) chiara: EBEB and !EBEB
    // if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( abs(etaPhot1)<1.5 && abs(etaPhot2)<1.5 )"));
    // if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( abs(etaPhot1)>1.5 || abs(etaPhot2)>1.5 )"));

    // b) chiara: altoR9 and !altoR9
    // if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1>0.94 && r9Phot2>0.94 )"));
    // if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1<0.94 || r9Phot2<0.94 )"));

    // c) chiara: EBEB and highR9 or NOT
    // if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1>0.94 && r9Phot2>0.94 && abs(etaPhot1)<1.5 && abs(etaPhot2)<1.5)"));
    // if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && ( r9Phot1<0.94 || r9Phot2<0.94 || abs(etaPhot1)>1.5 || abs(etaPhot2)>1.5)"));

    // d) chiara: 1-2 btag
    if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==1"));
    if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==2"));

    // e) chiara: 2-btag, 1-btag high R9, 1-btag low R9
    // if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==1 && r9Phot1>0.94 && r9Phot2>0.94"));
    // if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==1 && (r9Phot1<0.94 || r9Phot2<0.94)"));
    // if (c==2) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut+TString::Format(" && btagCategory==2"));

    cout << endl; cout << "for category = " << c << endl;
    dataToFit[c]->Print("v");
    cout << "---- nX:  " << dataToFit[c]->sumEntries() << endl;

    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }

  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("massggnewvtx"),mainCut);
  w->import(*data, Rename("Data"));
  cout << endl;
  cout << "data, no split" << endl;
  data->Print("v");
  cout << "---- nX:  " << data->sumEntries() << endl; 
}

// Fit signal with model pdfs
void SigModelFit(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  RooAbsPdf* MassggnewvtxSig[NCAT];

  Float_t MASS(mass);  
  Float_t minMassFit(MMIN),maxMassFit(MMAX); 

  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    MassggnewvtxSig[c] = (RooAbsPdf*)  w->pdf("MassggnewvtxSig"+TString::Format("_cat%d",c));
    ((RooRealVar*) w->var("massggnewvtx"+TString::Format("_sig_m0_cat%d",c)))->setVal(MASS);

    cout << "---------------- Peak Val = " 
    	 << w->var("massggnewvtx"+TString::Format("_sig_m0_cat%d",c))->getVal() 
      	 << ", Mass = " << MASS << endl;
      
    MassggnewvtxSig[c] ->fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE));

    // Plot to verify everything is ok
    RooPlot* plotMassggnewvtxAll = w->var("massggnewvtx")->frame(Range(100,180),Bins(19));
    sigToFit[c]->plotOn(plotMassggnewvtxAll);
    MassggnewvtxSig[c]->plotOn(plotMassggnewvtxAll);
    if (c==0) {
      MassggnewvtxSig[c]->plotOn(plotMassggnewvtxAll,Components("MassggnewvtxGaussSig_cat0"),LineStyle(kDashed),LineColor(kGreen));
      MassggnewvtxSig[c]->plotOn(plotMassggnewvtxAll,Components("MassggnewvtxCBSig_cat0"),LineStyle(kDashed),LineColor(kRed));
    } else if (c==1) {
      MassggnewvtxSig[c]->plotOn(plotMassggnewvtxAll,Components("MassggnewvtxGaussSig_cat1"),LineStyle(kDashed),LineColor(kGreen));
      MassggnewvtxSig[c]->plotOn(plotMassggnewvtxAll,Components("MassggnewvtxCBSig_cat1"),LineStyle(kDashed),LineColor(kRed));
    }
    MassggnewvtxSig[c]->paramOn(plotMassggnewvtxAll, ShowConstants(true), Layout(0.15,0.55,0.9), Format("NEU",AutoPrecision(2)));
    TCanvas* c1 = new TCanvas("c1","Massggnewvtx",0,0,500,500);
    c1->cd(1);
    plotMassggnewvtxAll->Draw();  
    c1->SaveAs("prelimSignal"+TString::Format("_cat%d.png",c));
    c1->SaveAs("prelimSignal"+TString::Format("_cat%d.root",c));

    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("SigPdfParam_cat%d",c), RooArgSet(*w->var("massggnewvtx"+TString::Format("_sig_m0_cat%d",c)),
								   *w->var("massggnewvtx"+TString::Format("_sig_sigma0_cat%d",c)),
								   *w->var("massggnewvtx"+TString::Format("_sig_sigma1_cat%d",c)),
								   *w->var("massggnewvtx"+TString::Format("_sig_alpha_cat%d",c)),
								   *w->var("massggnewvtx"+TString::Format("_sig_n_cat%d",c)), 
								   *w->var("massggnewvtx"+TString::Format("_sig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));
  }
}

RooFitResult* BkgModelFitExpo(RooWorkspace* w, Bool_t dobands) {

  Int_t ncat = NCAT;
  RooDataSet* data[NCAT];
  RooAbsPdf* MassggnewvtxBkg[NCAT];

  RooFitResult* fitresult[NCAT];;
  RooPlot* plotMassggnewvtxBkg[NCAT];

  Float_t minMassFit(MMIN),maxMassFit(MMAX); 
  
  // Fit data with background pdf for data limit
  RooRealVar* massggnewvtx = w->var("massggnewvtx");  
  massggnewvtx->setUnit("GeV");
  
  for (int c=0; c<ncat; ++c) {
    cout << "---------- category = " << c << endl;
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    MassggnewvtxBkg[c] = (RooAbsPdf*)  w->pdf("MassggnewvtxBkg"+TString::Format("_cat%d",c));
    fitresult[c] = MassggnewvtxBkg[c]->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));

    // Plot to verify everything is ok 
    plotMassggnewvtxBkg[c] = w->var("massggnewvtx")->frame(Range(100,180),Bins(19));
    data[c]->plotOn(plotMassggnewvtxBkg[c]);
    MassggnewvtxBkg[c]->plotOn(plotMassggnewvtxBkg[c]);
    MassggnewvtxBkg[c]->paramOn(plotMassggnewvtxBkg[c], ShowConstants(true), Layout(0.15,0.55,0.9), Format("NEU",AutoPrecision(2)));
    TCanvas* c1 = new TCanvas("c1","Massggnewvtx",0,0,500,500);
    c1->cd(1);
    plotMassggnewvtxBkg[c]->Draw();
    c1->SaveAs("prelimBkg"+TString::Format("_cat%d.png",c));
    c1->SaveAs("prelimBkg"+TString::Format("_cat%d.root",c));

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = MassggnewvtxBkg[c];
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotMassggnewvtxBkg[c]->getObject(1));
      
      for (int i=1; i<(plotMassggnewvtxBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotMassggnewvtxBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotMassggnewvtxBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotMassggnewvtxBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	massggnewvtx->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      massggnewvtx->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotMassggnewvtxBkg[c]->Draw("SAME"); 
    }
  }
  return fitresult;
}

RooFitResult* BkgModelFitBernstein(RooWorkspace* w, Bool_t dobands) {

  Int_t ncat = NCAT;
  
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
  RooBernstein* MassggnewvtxBkg[NCAT];
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotMassggnewvtxBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];
  RooAbsPdf* MassggnewvtxSig[NCAT];
  
  Float_t minMassFit(MMIN),maxMassFit(MMAX); 
  
  // Fit data with background pdf for data limit
  RooRealVar* massggnewvtx = w->var("massggnewvtx");  
  massggnewvtx->setUnit("GeV");
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
    // fit a la dijets
    ((RooRealVar*) w->var(TString::Format("massggnewvtx_bkg_8TeV_slope3_cat%d",c)))->setConstant(true);
    cout << "---------------- Parameter 3 set to const" << endl;
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("p1mod_cat%d",c),"","@0",*w->var(TString::Format("massggnewvtx_bkg_8TeV_slope1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("p2mod_cat%d",c),"","@0",*w->var(TString::Format("massggnewvtx_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("p3mod_cat%d",c),"","@0",*w->var(TString::Format("massggnewvtx_bkg_8TeV_slope3_cat%d",c)));
    RooFormulaVar *sqrtS = new RooFormulaVar(TString::Format("sqrtS_cat%d",c),"","@0",*w->var("sqrtS"));
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("x_cat%d",c),"","@0/@1",RooArgList(*massggnewvtx, *sqrtS));

    RooAbsPdf* MassggnewvtxBkgTmp0 = new RooGenericPdf(TString::Format("DijetBackground_%d",c), "pow(1-@0, @1)/pow(@0, @2+@3*log(@0))", RooArgList(*x, *p1mod, *p2mod, *p3mod));
    w->factory(TString::Format("massggnewvtx_bkg_8TeV_norm_cat%d[4000.0,0.0,10000000]",c));
    
    RooExtendPdf MassggnewvtxBkgTmp(TString::Format("MassggnewvtxBkg_cat%d",c),"",*MassggnewvtxBkgTmp0,*w->var(TString::Format("massggnewvtx_bkg_8TeV_norm_cat%d",c)));
    
    fitresult[c] = MassggnewvtxBkgTmp.fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));
    w->import(MassggnewvtxBkgTmp);

    //************************************************
    // Plot Massggnewvtx background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","Massggnewvtx Background Categories",0,0,500,500);
    Int_t nBinsMass(19);
    plotMassggnewvtxBkg[c] = massggnewvtx->frame(nBinsMass);
    data[c]->plotOn(plotMassggnewvtxBkg[c],LineColor(kWhite),MarkerColor(kWhite));    
    
    MassggnewvtxBkgTmp.plotOn(plotMassggnewvtxBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
    data[c]->plotOn(plotMassggnewvtxBkg[c]);    
    plotMassggnewvtxBkg[c]->Draw();  


    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = MassggnewvtxBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotMassggnewvtxBkg[c]->getObject(1));
      
      for (int i=1; i<(plotMassggnewvtxBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotMassggnewvtxBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotMassggnewvtxBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotMassggnewvtxBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	massggnewvtx->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      massggnewvtx->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotMassggnewvtxBkg[c]->Draw("SAME"); 
    }
  }
  return fitresult;
}

void SetConstantParams(const RooArgSet* params) {

  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}

void MakePlots(RooWorkspace* w, Float_t mass, RooFitResult* fitresults) {

  Int_t ncat = NCAT;
  
  cout << endl; cout << "Retreive everything:" << endl; 
  w->Print();

  // retrieve data sets from the workspace
  RooDataSet* dataAll   = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll = (RooDataSet*) w->data("SigWeight");
  
  // maximum 9 cat...
  RooDataSet* data[9];  
  RooDataSet* signal[9];
  RooAbsPdf*  MassggnewvtxGaussSig[9];
  RooAbsPdf*  MassggnewvtxCBSig[9];
  RooAbsPdf*  MassggnewvtxSig[9];
  RooAbsPdf*  MassggnewvtxBkg[9];  
  
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    MassggnewvtxGaussSig[c]  = (RooAbsPdf*)  w->pdf(TString::Format("MassggnewvtxGaussSig_cat%d",c));
    MassggnewvtxCBSig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("MassggnewvtxCBSig_cat%d",c));
    MassggnewvtxSig[c]       = (RooAbsPdf*)  w->pdf("MassggnewvtxSig"+TString::Format("_cat%d",c));
    MassggnewvtxBkg[c]       = (RooAbsPdf*)  w->pdf(TString::Format("MassggnewvtxBkg_cat%d",c));
  }
  
  // retrieve mass observable from the workspace
  RooRealVar* massggnewvtx = w->var("massggnewvtx");  
  massggnewvtx->setUnit("GeV");
  
  // retrieve pdfs after the fits
  RooAbsPdf* MassggnewvtxGaussSigAll = w->pdf("MassggnewvtxGaussSig");
  RooAbsPdf* MassggnewvtxCBSigAll    = w->pdf("MassggnewvtxCBSig");
  RooAbsPdf* MassggnewvtxSigAll      = w->pdf("MassggnewvtxSig");
  RooAbsPdf* MassggnewvtxBkgAll      = w->pdf("MassggnewvtxBkgAll");

  Float_t minMassFit(MMIN),maxMassFit(MMAX); 
  Float_t MASS(mass);
  int iMass = abs(mass);
  
  Int_t nBinsMass(19);  // chiara


  // ****************************
  cout << endl; cout << "Progress plotting: signal" << endl;     

  RooPlot* plotMassggnewvtxAll = massggnewvtx->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
  signalAll->plotOn(plotMassggnewvtxAll);
  
  gStyle->SetOptTitle(0);
  MassggnewvtxSigAll->plotOn(plotMassggnewvtxAll);
  MassggnewvtxSigAll->plotOn(plotMassggnewvtxAll,Components("MassggnewvtxGaussSig"),LineStyle(kDashed),LineColor(kGreen));
  MassggnewvtxSigAll->plotOn(plotMassggnewvtxAll,Components("MassggnewvtxCBSig"),LineStyle(kDashed),LineColor(kRed));
  MassggnewvtxSigAll->paramOn(plotMassggnewvtxAll, ShowConstants(true), Layout(0.15,0.55,0.9), Format("NEU",AutoPrecision(2)));
  plotMassggnewvtxAll->getAttText()->SetTextSize(0.03);
  
  TCanvas* c1 = new TCanvas("c1","Massggnewvtx",0,0,500,500);
  c1->cd(1);
  plotMassggnewvtxAll->Draw();  
  c1->SaveAs("plots/sigmodel_"+TString::Format("%d.png", iMass));
  c1->SaveAs("plots/sigmodel_"+TString::Format("%d.root", iMass));


  // ****************************
  cout << endl; cout << "Progress plotting: signal per categories" << endl;     
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  RooPlot* plotMassggnewvtx[9];
  for (int c=0; c<ncat; ++c) {
    plotMassggnewvtx[c] = massggnewvtx->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
    signal[c]->plotOn(plotMassggnewvtx[c],LineColor(kWhite),MarkerColor(kWhite));    
    MassggnewvtxSig[c]  ->plotOn(plotMassggnewvtx[c]);
    MassggnewvtxSig[c]  ->plotOn(plotMassggnewvtx[c],Components("MassggnewvtxGaussSig"+TString::Format("_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
    MassggnewvtxSig[c]  ->plotOn(plotMassggnewvtx[c],Components("MassggnewvtxCBSig"+TString::Format("_cat%d",c)),LineStyle(kDashed),LineColor(kRed));
    MassggnewvtxSig[c]  ->paramOn(plotMassggnewvtx[c], ShowConstants(true), Layout(0.15,0.55,0.9), Format("NEU",AutoPrecision(2)));
    signal[c]  ->plotOn(plotMassggnewvtx[c]);

    TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F("hist", "hist", 400, minMassFit, maxMassFit);
 
    plotMassggnewvtx[c]->SetTitle("");      
    plotMassggnewvtx[c]->SetMinimum(0.0);
    plotMassggnewvtx[c]->SetMaximum(1.40*plotMassggnewvtx[c]->GetMaximum());
    plotMassggnewvtx[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");

    TCanvas* ctmp = new TCanvas("ctmp","Massggnewvtx Background Categories",0,0,500,500);
    plotMassggnewvtx[c]->Draw();  
    
    plotMassggnewvtx[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
    legmc->AddEntry(plotMassggnewvtx[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotMassggnewvtx[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotMassggnewvtx[c]->getObject(3),"Crystal Ball component","L");
    legmc->AddEntry(plotMassggnewvtx[c]->getObject(2),"Gaussian Outliers","L");
    
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    
    TLatex *lat  = new TLatex(minMassFit+1.5,0.85*plotMassggnewvtx[c]->GetMaximum(),"#scale[1.0]{CMS Preliminary}");
    lat->Draw();
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d.png", iMass, c));
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d.root", iMass, c));
  }


  // ****************************
  cout << endl; cout << "Progress plotting: background" << endl;     
  TCanvas* c4 = new TCanvas("c4","Massggnewvtx Background Categories",0,0,400,400);
  RooPlot* plotMassggnewvtxBkg[9];
  for (int c=0; c<ncat; ++c) {
    plotMassggnewvtxBkg[c] = massggnewvtx->frame(nBinsMass);
    data[c]->plotOn(plotMassggnewvtxBkg[c],LineColor(kWhite),MarkerColor(kWhite));    
    MassggnewvtxBkg[c]->plotOn(plotMassggnewvtxBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
    data[c]->plotOn(plotMassggnewvtxBkg[c]);    
    MassggnewvtxBkg[c]->paramOn(plotMassggnewvtxBkg[c], ShowConstants(true), Layout(0.65,0.9,0.9), Format("NEU",AutoPrecision(4)));
    plotMassggnewvtxBkg[c]->getAttText()->SetTextSize(0.03);
    plotMassggnewvtxBkg[c]->Draw();  
    gPad->SetLogy(1);
    plotMassggnewvtxBkg[c]->SetAxisRange(0.1,plotMassggnewvtxBkg[c]->GetMaximum()*1.5,"Y");
    c4->SaveAs("plots/backgrounds_log_"+TString::Format("cat%d.png", c));
    c4->SaveAs("plots/backgrounds_log_"+TString::Format("cat%d.root", c));
  }
  
  TCanvas* c5 = new TCanvas("c5","Massggnewvtx Background Categories",0,0,400,400);
  RooPlot* plotMassggnewvtxBkg[9];
  for (int c=0; c<ncat; ++c) {
    plotMassggnewvtxBkg[c] = massggnewvtx->frame(nBinsMass);
    data[c]->plotOn(plotMassggnewvtxBkg[c],LineColor(kWhite),MarkerColor(kWhite));    
    MassggnewvtxBkg[c]->plotOn(plotMassggnewvtxBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
    data[c]->plotOn(plotMassggnewvtxBkg[c]);    
    MassggnewvtxBkg[c]->paramOn(plotMassggnewvtxBkg[c], ShowConstants(true), Layout(0.65,0.9,0.9), Format("NEU",AutoPrecision(4)));
    plotMassggnewvtxBkg[c]->getAttText()->SetTextSize(0.03);
    plotMassggnewvtxBkg[c]->Draw();  
    c5->SaveAs("plots/backgrounds_"+TString::Format("cat%d.png", c));
    c5->SaveAs("plots/backgrounds_"+TString::Format("cat%d.root", c));
  }
}

// Write signal pdfs and datasets into the workspace 
void MakeSigWS(RooWorkspace* w, const char* fileBaseName){
  
  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");  

  //********************************//
  // Retrieve P.D.F.s
  RooAbsPdf* MassggnewvtxSigPdf[NCAT];
  for (int c=0; c<ncat; ++c) {
    MassggnewvtxSigPdf[c] = (RooAbsPdf*)  w->pdf("MassggnewvtxSig"+TString::Format("_cat%d",c));
    wAll->import(*w->pdf("MassggnewvtxSig"+TString::Format("_cat%d",c)));
  }
  std::cout << "done with importing signal pdfs" << std::endl;

  // (2) Systematics on energy scale and resolution // chiara: per ora tutte le sistematiche non hanno senso
  // wAll->factory("CMS_hgg_sig_m0_absShift[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat0(massggnewvtx_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat1(massggnewvtx_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");

  // (3) Systematics on resolution: create new sigmas
  // wAll->factory("CMS_hgg_sig_sigmaScale[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat0(massggnewvtx_sig_sigma0_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(massggnewvtx_sig_sigma1_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat1(massggnewvtx_sig_sigma0_cat1, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(massggnewvtx_sig_sigma1_cat1, CMS_hgg_sig_sigmaScale)")

  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  
  return;
}

// Write background pdfs and datasets into the workspace 
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName) {

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;  

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  
  //********************************//
  // Retrieve the datasets and PDFs
  RooDataSet* data[NCAT];
  RooAbsPdf* MassggnewvtxBkgPdf[NCAT];
  for (int c=0; c<ncat; ++c) {

    /*
    cout << "For category " << c << endl;
    data[c]      = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    ((RooRealVar*) data[c]->get()->find("massggnewvtx"))->setBins(MMAX-MMIN) ;
    RooDataHist* dataBinned = data[c]->binnedClone();
    MassggnewvtxBkgPdf[c] = (RooExtendPdf*)  w->pdf(TString::Format("MassggnewvtxBkg_cat%d",c));
    //   wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));
    wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c)));
    wAll->import(*w->pdf(TString::Format("MassggnewvtxBkg_cat%d",c)));

    double mean = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_norm_cat%d",c))->getVal();
    double min = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_norm_cat%d",c))->getMin();
    double max = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_norm_cat%d",c))->getMax();
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_cat%d_norm[%g,%g,%g]", c, mean, min, max));
    
    double mean = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope1_cat%d",c))->getVal();
    double min = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope1_cat%d",c))->getMin();
    double max = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope1_cat%d",c))->getMax();
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_slope1_cat%d[%g,%g,%g]", c, mean, min, max));
    
    double mean = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope2_cat%d",c))->getVal();
    double min = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope2_cat%d",c))->getMin();
    double max = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope2_cat%d",c))->getMax();
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_slope2_cat%d[%g,%g,%g]", c, mean, min, max));
    
    double mean = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope3_cat%d",c))->getVal();
    
    double min = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope3_cat%d",c))->getMin();
    double max = wAll->var(TString::Format("massggnewvtx_bkg_8TeV_slope3_cat%d",c))->getMax();
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_slope3_cat%d[%g,%g,%g]", c, mean, mean, mean));
    */

    MassggnewvtxBkgPdf[c] = (RooAbsPdf*) w->pdf("MassggnewvtxBkg"+TString::Format("_cat%d",c));

    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    ((RooRealVar*) data[c]->get()->find("massggnewvtx"))->setBins(MMAX-MMIN) ;
    RooDataHist* dataBinned = data[c]->binnedClone();
    
    wAll->import(*w->pdf("MassggnewvtxBkg"+TString::Format("_cat%d",c)));
    wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c)));
  }
  std::cout << "done with importing background pdfs" << std::endl;
  
  // (2) do reparametrization of background
  // for (int c=0; c<ncat; ++c) {
  //wAll->factory(                       
    //	  TString::Format("EDIT::CMS_hgg_bkg_8TeV_cat%d(MassggnewvtxBkg_cat%d,",c,c) +
    //	  TString::Format(" massggnewvtx_bkg_exp_cat%d=CMS_hgg_bkg_8TeV_cat%d_norm,", c,c)
    //	  );
  //} 

  TString filename(wsDir+TString(fileBaseName)+".root");
  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;

  std::cout << std::endl; 
  std::cout << "observation:" << std::endl;
  for (int c=0; c<ncat; ++c) {
    std::cout << "cat " << c << ", " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries() << endl;
    wAll->data(TString::Format("data_obs_cat%d",c))->Print();
  }
  std::cout << std::endl;
  
  return;
}

Double_t effSigma(TH1 *hist) {
  
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

// preparing datacards
void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, int iChan) {

  TString cardDir = "datacards/"+filePOSTfix;
  Int_t ncat = NCAT;
  TString wsDir   = "../workspaces/"+filePOSTfix;

  // **********************
  // Retrieve the datasets
  cout << "Start retrieving dataset" << endl;
  
  RooDataSet* data[9];
  RooDataSet* signal[9];
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
  }

  RooRealVar*  lumi = w->var("lumi");

  // *****************************
  // Print Expected event yields
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data: " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("SigWeight")->sumEntries()  << endl;
  Float_t siglikeErr[6];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << signal[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*signal[c]->sumEntries();
  }
  cout << "====================================================" << endl;  


  // *****************************
  // Print Data Card int file
  TString filename(cardDir+TString(fileBaseName)+"_"+"_8TeV"+Form("_channel%d.txt",iChan));
  ofstream outFile(filename);

  outFile << "#CMS-HGGBB DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d datacardName.txt -U -m 125 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax *" << endl;
  outFile << "jmax *" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  outFile << "shapes data_obs * " << wsDir+TString(fileBkgName)+".root" << " w_all:data_obs_$CHANNEL" << endl;
  outFile << "shapes sig * "      << wsDir+TString(fileBaseName)+"_8TeV"+".inputsig.root" << " w_all:MassggnewvtxSig_$CHANNEL" << endl;
  outFile << "shapes bkg * "      << wsDir+TString(fileBkgName)+".root" << " w_all:MassggnewvtxBkg_$CHANNEL" << endl;

  outFile << "---------------" << endl;
  outFile << Form("bin          cat%d", iChan) << endl;
  outFile <<  "observation   "  <<  Form("%.10lg",data[iChan]->sumEntries()) << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                 " << Form("cat%d      cat%d", iChan, iChan) << endl;
  outFile << "process                 sig      bkg" << endl;
  outFile << "process                   0        1" << endl;
  // if(signalScaler==1.)
  // signalScaler=1./signal[2]->sumEntries()*20;
  outFile << "rate                   " 
	  << signal[iChan]->sumEntries()*signalScaler << " " << data[iChan]->sumEntries() << endl;
  outFile << "--------------------------------" << endl;
  outFile << "# signal scaled by " << signalScaler << endl;

  outFile << "lumi_8TeV       lnN  0.950/1.050  - " << endl;
  // outFile << "CMS_VV_eff_g         lnN  0.8/1.20      - # Signal Efficiency" << endl;
  // outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  // outFile << Form("CMS_hgg_sig_m0_absShift    param   1   0.0125   # displacement of the mean w.r.t. nominal in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("CMS_hgg_sig_sigmaScale     param   1   0.1   # multiplicative correction to sigmas in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_cat%d_norm           flatParam  # Normalization uncertainty on background slope",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope2_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope3_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // if (iChan != 2 )  outFile << Form("CMS_hgg_bkg_8TeV_slope1_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
}

Double_t effSigma(TH1 *hist) {

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

void R2JJFitter(double mass, std::string postfix="")
{
    filePOSTfix=postfix;
    if(postfix!="")
    {
      // for optimization studies
      MMIN=1000;
      if(mass==1000)
         signalScaler=0.034246;
      if((mass>1000)&&(mass<2000))
         signalScaler=0.02469;
      if(mass==2000)
         signalScaler=2.0;
    };
    runfits(mass, 1);
    if(postfix!="")
    {
      // for optimization studies
      MMIN=1000;
      if(mass==1000)
         signalScaler=0.033500;
      if((mass>1000)&&(mass<2000))
         signalScaler=0.02016;
      if(mass==2000)
         signalScaler=2.22222;
    };
    runfits(mass, 0);
    runfits(mass, 2);
}
