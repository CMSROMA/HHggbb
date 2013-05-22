massggnewvtx[0,180];

massggnewvtx_sig_m0[124.965, 120.0, 130.0];
massggnewvtx_sig_sigma0[3.5431, 0.0, 5.0];
massggnewvtx_sig_sigma1[1.5412, 0.0, 5.0];
massggnewvtx_sig_alpha[4.98665, 0., 5.]; 
massggnewvtx_sig_n[0.000006687, 0., 0.5]; 
massggnewvtx_sig_frac[0.244542, 0, 1.0];

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

massggnewvtx_bkg_exp[-0.022,-0.05, 0.];
MassggnewvtxBkg = Exponential(massggnewvtx, massggnewvtx_bkg_exp);

massggnewvtx_bkg_exp_cat0[-0.022,-0.05, 0.];
MassggnewvtxBkg_cat0 = Exponential(massggnewvtx, massggnewvtx_bkg_exp_cat0);

wei[1,0,10];

sqrtS[8000., 8000., 8000.]
