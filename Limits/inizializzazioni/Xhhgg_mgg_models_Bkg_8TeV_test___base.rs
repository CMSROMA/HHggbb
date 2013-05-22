massggnewvtx[0,180];

massggnewvtx_sig_m0[125, 120.0, 130.0];
massggnewvtx_sig_sigma0[2.9457, 0.0, 5.0];
massggnewvtx_sig_sigma1[1.2181, 0.0, 5.0];
massggnewvtx_sig_alpha[3.6665, 0., 5.]; 
massggnewvtx_sig_n[0.0000099936, 0., 0.5]; 
massggnewvtx_sig_frac[0.44303, 0, 1.0];

MassggnewvtxGaussSig = Gaussian(massggnewvtx, massggnewvtx_sig_m0, massggnewvtx_sig_sigma0);
MassggnewvtxCBSig    = CBShape(massggnewvtx, massggnewvtx_sig_m0, massggnewvtx_sig_sigma1, massggnewvtx_sig_alpha, massggnewvtx_sig_n);
MassggnewvtxSig      = AddPdf(MassggnewvtxGaussSig, MassggnewvtxCBSig, massggnewvtx_sig_frac);

massggnewvtx_sig_m0_cat0[124.965, 120.0, 130.0];
massggnewvtx_sig_sigma0_cat0[3.5431, 0.0, 5.0];
massggnewvtx_sig_sigma1_cat0[1.5412, 0.0, 5.0];
massggnewvtx_sig_alpha_cat0[4.98665, 0., 5.]; 
massggnewvtx_sig_n_cat0[0.000006687, 0., 0.5]; 
massggnewvtx_sig_frac_cat0[0.244542, 0, 1.0];

MassggnewvtxGaussSig_cat0 = Gaussian(massggnewvtx, massggnewvtx_sig_m0_cat0, massggnewvtx_sig_sigma0_cat0);
MassggnewvtxCBSig_cat0    = CBShape(massggnewvtx, massggnewvtx_sig_m0_cat0, massggnewvtx_sig_sigma1_cat0, massggnewvtx_sig_alpha_cat0, massggnewvtx_sig_n_cat0);
MassggnewvtxSig_cat0      = AddPdf(MassggnewvtxGaussSig_cat0, MassggnewvtxCBSig_cat0, massggnewvtx_sig_frac_cat0);

massggnewvtx_sig_m0_cat1[300.0, 100.0, 1000.0];
massggnewvtx_sig_sigma_cat1[35., 5.0, 200.0];
massggnewvtx_sig_alpha_cat1[ 0.8, 0.0, 3.0]; 
massggnewvtx_sig_n_cat1[130, 0.00001, 1000.0]; 
massggnewvtx_sig_gsigma_cat1[100, 50.0, 200.0];
massggnewvtx_sig_frac_cat1[0.5, 0, 1.0];

MassggnewvtxGaussSig_cat1 = Gaussian(massggnewvtx, massggnewvtx_sig_m0_cat1, massggnewvtx_sig_gsigma_cat1);
MassggnewvtxCBSig_cat1    = CBShape(massggnewvtx, massggnewvtx_sig_m0_cat1, massggnewvtx_sig_sigma_cat1, massggnewvtx_sig_alpha_cat1, massggnewvtx_sig_n_cat1);
MassggnewvtxSig_cat1      = AddPdf(MassggnewvtxGaussSig_cat1, MassggnewvtxCBSig_cat1, massggnewvtx_sig_frac_cat1);

massggnewvtx_bkg_8TeV_slope[1.0,0, 1];
massggnewvtx_bkg_8TeV_slope1[7,0.0, 100.0];
massggnewvtx_bkg_8TeV_slope2[5,0.0, 100.0];
massggnewvtx_bkg_8TeV_slope3[0.,0.0, 0.0];

massggnewvtx_bkg_8TeV_slope_cat0[1000.0,0, 10000000];
massggnewvtx_bkg_8TeV_slope1_cat0[10., -100.0, 500.0];
massggnewvtx_bkg_8TeV_slope2_cat0[5.,0.0, 100.0];
massggnewvtx_bkg_8TeV_slope3_cat0[0.,-10.0, 10.0];

massggnewvtx_bkg_8TeV_slope_cat1[1000.0,0, 10000000];
massggnewvtx_bkg_8TeV_slope1_cat1[10., -100.0, 500.0];
massggnewvtx_bkg_8TeV_slope2_cat1[5.,0.0, 100.0];
massggnewvtx_bkg_8TeV_slope3_cat1[0.,-10.0, 10.0];

massggnewvtx_bkg_exp_cat0[-0.022,-0.05, 0.];
MassggnewvtxBkg_cat0 = Exponential(massggnewvtx, massggnewvtx_bkg_exp_cat0);

wei[1,0,10];

sqrtS[8000., 8000., 8000.]
