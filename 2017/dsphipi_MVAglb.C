#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace RooFit;

//gErrorIgnoreLevel = kFatal;

void dsphipi_MVAglb() 
{
    bool doPUrew = true;
    
    gErrorIgnoreLevel = kInfo;

    //open root files where to store all plots
    TFile *fout = new TFile("plots_MVAcut/dsphipi_MVAglb_2017_output.root", "RECREATE");
    fout->cd();
    
    int which_etabin = 2;//1: barrel, 2:endcap
    
    double var1[] = {0.1,0.4825939929473169, 0.5414023593110158, 0.5700989959305268, 0.589402134098577, 0.6045252383372394, 0.6178331517033439, 0.6310203868218345, 0.6466055250200856, 0.8};
    double var2[] = {0.1,0.43649083638325026, 0.5165698177711587, 0.557519695439444, 0.5844797312847071, 0.605730110146498, 0.6247331238530781, 0.6435834258706714, 0.665673969188595, 0.8};
    double var3[] = {0.0,1.2, 2.4};
    
    //TFile *fout0 = new TFile("dsphipi_MVAglb2_2017_output.root", "RECREATE");
    //TH1D *SF_mva1 = new TH1D("SF_mva1", "SF vs MVA", 9, var1 ) ;
    //TH1D *SF_mva2 = new TH1D("SF_mva2", "SF vs MVA", 9, var2 ) ;
    TFile *fout0 = new TFile("dsphipi_MVAglb2_2017_output.root", "UPDATE");
    TH1D *SF_mva1 = ((TH1D*) fout0->Get("SF_mva1"));
    TH1D *SF_mva2 = ((TH1D*) fout0->Get("SF_mva2"));
    
    TString bdt_cutnos1[] = {"0.1","0.4825939929473169", "0.5414023593110158", "0.5700989959305268", "0.589402134098577", "0.6045252383372394", "0.6178331517033439", "0.6310203868218345", "0.6466055250200856", "0.8"};
    TString bdt_cutnos2[] = {"0.1","0.43649083638325026", "0.5165698177711587", "0.557519695439444", "0.5844797312847071", "0.605730110146498", "0.6247331238530781", "0.6435834258706714", "0.665673969188595", "0.8"};
    
    int no_bin = (which_etabin==1)?11:5;
    
    TString bdt_cutnos[10];
    TString eta_cut;
    
    if(which_etabin==1){
      copy(begin(bdt_cutnos1), end(bdt_cutnos1), begin(bdt_cutnos));
      eta_cut = " && abs(mu2_eta)<1.2";
    }
    else if(which_etabin==2){
      copy(begin(bdt_cutnos2), end(bdt_cutnos2), begin(bdt_cutnos));
      eta_cut = " && abs(mu2_eta)>=1.2";
    }
    
    //set here list of bdt cuts
    
    
    TString bdt_cutlist[] = {
                              "MVA_glbmu2< "+bdt_cutnos[0],
                              "MVA_glbmu2>="+bdt_cutnos[0]+" && MVA_glbmu2<"+bdt_cutnos[1],
                              "MVA_glbmu2>="+bdt_cutnos[1]+" && MVA_glbmu2<"+bdt_cutnos[2],
                              "MVA_glbmu2>="+bdt_cutnos[2]+" && MVA_glbmu2<"+bdt_cutnos[3],
                              "MVA_glbmu2>="+bdt_cutnos[3]+" && MVA_glbmu2<"+bdt_cutnos[4],
                              "MVA_glbmu2>="+bdt_cutnos[4]+" && MVA_glbmu2<"+bdt_cutnos[5],
                              "MVA_glbmu2>="+bdt_cutnos[5]+" && MVA_glbmu2<"+bdt_cutnos[6],
                              "MVA_glbmu2>="+bdt_cutnos[6]+" && MVA_glbmu2<"+bdt_cutnos[7],
                              "MVA_glbmu2>="+bdt_cutnos[7]+" && MVA_glbmu2<"+bdt_cutnos[8],
                              "MVA_glbmu2>="+bdt_cutnos[8]+" && MVA_glbmu2<"+bdt_cutnos[9],
                              "MVA_glbmu2>="+bdt_cutnos[9],
                             };
                             
    //TString bdt_cutlist00[] = bdt_cutlist;
    
    /*
    TString bdt_cutlist2[] = {
                              "MVA_glbmu2< "+bdt_cutnos[0],
                              "MVA_glbmu2>="+bdt_cutnos[0]+" && MVA_glbmu2<"+bdt_cutnos[3],
                              "MVA_glbmu2>="+bdt_cutnos[3]+" && MVA_glbmu2<"+bdt_cutnos[6],
                              "MVA_glbmu2>="+bdt_cutnos[6]+" && MVA_glbmu2<"+bdt_cutnos[9],
                              "MVA_glbmu2>="+bdt_cutnos[9],
                             };
                             
    TString bdt_cutlist[5];
    if(which_etabin==1){
      copy(begin(bdt_cutlist1), end(bdt_cutlist1), begin(bdt_cutlist));
    }
    else if(which_etabin==2){
      copy(begin(bdt_cutlist2), end(bdt_cutlist2), begin(bdt_cutlist));
    }
    */
    
    int ncut = sizeof(bdt_cutlist)/sizeof(bdt_cutlist[0]);

    cout<<"Size of ncut: "<<ncut<<endl;
    
    TString filename_text = "dsphipi_yield_MVAglb_2017.txt";
    //open text file where to write Ds yields
    ofstream fout_yield(filename_text);
    
    Double_t Dsyield_data[ncut];
    Double_t Dsyield_data_err[ncut];
    Double_t Dsyield_MC[ncut];
    Double_t Dsyield_MC_err[ncut];
    Double_t Bkgyield_data[ncut];
    Double_t Bkgyield_data_err[ncut];
    
    TFile *fin = new TFile("t3mminitree_xgb_2017_13may22_ptmu3_control.root", "READ");
    cout<<"opened input file t3mminitree_xgb_2017_13may22_ptmu3_control.root"<<endl;
    TTree *tin = (TTree*)fin->Get("OutputTree");
    
    TH1F *h_tripletmass[ncut];//all data range with mva cut
    TH1F *h_tripletmass_mc[ncut];//peak MC range with mva cut

    TH1F *h_tripletmass_full;
    TH1F *h_tripletmass_mc_full;
	
    Double_t events_SB = 0;
    TString binning_mass = "(62, 1.72, 2.01)";
    Int_t n_bins = 62;
	
    Double_t lumi_full = 37.9; //fb
    TString run_lable = "2017";
    Double_t xsection_mc = 1.06e10; //Ds Production Cross section
    int N_MC = 1479233;  //Total number of events in MC sample 2964234 1479233
    Double_t BR = 1.29e-5;  //Branching ratio Ds to Phi Pi

    TString common_cut = " && bs_sv_d2Dsig>2.0 && mu3_pt > 1.2 && " 
                         "((mu1_pt>3.5 && mu1_eta<1.2) || (mu1_pt>2.0 && mu1_eta>=1.2 && mu1_eta<=2.4)) && "
                         "((mu2_pt>3.5 && mu2_eta<1.2) || (mu2_pt>2.0 && mu2_eta>=1.2 && mu2_eta<=2.4)) && "
                         "abs(phiMass-1.02)<0.045 && "
                         "!(l1double_DoubleMu4_fired && !l1double_DoubleMu0_fired)";

    TString invmass_all_data = "";
    TString invmass_peak_MC = "";
    
	
    if(doPUrew){
        invmass_all_data  = "puFactor*(tripletMass<2.02 && tripletMass>1.62"+common_cut+eta_cut+" && isMC==0";
        invmass_peak_MC = "puFactor*(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==1";
    } else {
        invmass_all_data  = "(tripletMass<2.02 && tripletMass>1.62"+common_cut+eta_cut+" && isMC==0";
        invmass_peak_MC = "(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==1";
    }

    //Fill histograms for each cut
    for(int i=0; i<ncut; i++){
        TString s = std::to_string(i-1);
        
        TString bdt_cut = bdt_cutlist[i];
        TString bdt_cut_label = bdt_cutlist[i];
        bdt_cut_label = bdt_cut_label.ReplaceAll(".", "p");
        bdt_cut_label = bdt_cut_label.ReplaceAll(">", "_");
        bdt_cut_label = bdt_cut_label.ReplaceAll("<", "_");
    		
        tin->Draw("tripletMass>>h_tripletmass["+s+"]"+binning_mass, invmass_all_data+"&&"+bdt_cut+")");
        h_tripletmass[i]     = (TH1F *)gDirectory->Get("h_tripletmass["+s+"]");
        
        tin->Draw("tripletMass>>h_tripletmass_mc["+s+"]"+binning_mass, invmass_peak_MC+"&&"+bdt_cut+")");
        h_tripletmass_mc[i]     = (TH1F *)gDirectory->Get("h_tripletmass_mc["+s+"]");
        
        cout<<"Events passing selections "<<bdt_cutlist[i]<<" = "<<h_tripletmass[i]->GetEntries()<<endl;
    }
    //Fill histograms for full statistics
    tin->Draw("tripletMass>>h_tripletmass_full"+binning_mass, invmass_all_data+")");
    h_tripletmass_full     = (TH1F *)gDirectory->Get("h_tripletmass_full");
    
    tin->Draw("tripletMass>>h_tripletmass_mc_full"+binning_mass, invmass_peak_MC+")");
    h_tripletmass_mc_full   = (TH1F *)gDirectory->Get("h_tripletmass_mc_full");

    //Define roofit categories for simultaneous fit	
    std::map<std::string, TH1 *> hmap;
    RooCategory c("c", "c");
    for(int i = 0; i<ncut; i++){
        std::stringstream category;
        TString s = std::to_string(i);
        category << "cat" <<i;
        //define category
        c.defineType(category.str().c_str());
        //map category - TH1
        hmap[category.str()] = h_tripletmass[i];
    }
    RooRealVar x("x","3glb inv. mass (GeV)",1.72,2.01);
    RooDataHist dh("dh", "dh", x, Index(c), Import(hmap));

    // Construct a simultaneous pdf using category "c" as index
    RooSimultaneous simPdf("simPdf", "simultaneous pdf", c);
    
    
    //set ranges
    /*
    x.setRange("R1",1.825,1.89); //first peak D+(1.87GeV)
    x.setRange("R2",1.93,2.01); //second peak Ds(1.97)
    x.setRange("R3",1.72,1.815); //background    
    x.setRange("R4",1.90,1.925); //background    
    x.setRange("R5",1.995,2.01); //background    
    x.setRange("R6",1.72,2.01); //full range   
    

    RooRealVar mass2("mass2","Central value of Gaussian",1.965,1.94,2.0);
    RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.06);
    RooRealVar mass1("mass1","Central value of Gaussian",1.87,1.85,1.89);
    RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0003,0.3);
    */
    
    
    
    
    
    /*

    // Fit a Crystal ball p.d.f to the data
    RooRealVar alpha2("alpha2","alpha value CB",1.5,-20,20);
    RooRealVar n2("n2", "n2", 2, 0, 20);
    RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
    
    // Fit a Crystal ball p.d.f to the data
    RooRealVar alpha1("alpha1","alpha value CB",1.5,-20,20);
    RooRealVar n1("n1", "n1", 2, 0, 10);
    RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
       
    // Fit an exponential to the background
    RooRealVar a("a", "a", -5, -20, 0.0);
    RooExponential bg_exp("bg_exp", "bg_exp", x, a);
    bg_exp.fitTo(dh, Range("R3,R4"));
    */
    
    //x.setRange("R3",1.72,1.79); //background  
    //x.setRange("R4",1.906,1.919); //background 
    //RooRealVar mass1("mass1","Central value of Gaussian",1.86774,1.80,1.905);
    //RooRealVar mass2("mass2","Central value of Gaussian",1.96724,1.920,2.01);
    
    x.setRange("R3",1.72,1.82); //background    
    x.setRange("R4",1.90,1.925); //background    
    x.setRange("R5",2.0,2.01); //background 
    
    // Fit a Crystal ball p.d.f to the data
    RooRealVar mass1("mass1","Central value of Gaussian",1.86774,1.75,2.01);
    //RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0001,0.3);
    RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0001,0.3);
    
    RooRealVar alpha1("alpha1","alpha value CB",1.9,0.1,20);
    RooRealVar n1("n1", "n1", 2, 0, 10);
    RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
    
    // Fit a Crystal ball p.d.f to the data
    RooRealVar mass2("mass2","Central value of Gaussian",1.96724,1.85,2.51);
    //RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.3);
    RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.3);
    
    RooRealVar alpha2("alpha2","alpha value CB",1.9,0.1,20);
    RooRealVar n2("n2", "n2", 2, 0, 20);
    RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
       
    // Fit an exponential to the background
    RooRealVar a("a", "a", -5, -30, 5.0);
    RooExponential bg_exp("bg_exp", "bg_exp", x, a);
    //bg_exp.fitTo(dh, Range("R3,R4,R5"));
    bg_exp.fitTo(dh, Range("R3"));

    // Combine the models
    // relative contributions of background and signal will depend on category
    RooRealVar *nsig2[ncut];
    RooRealVar *nsig1[ncut];
    RooRealVar *nbkg[ncut];
    RooAddPdf *model[ncut];

    // Associate model with the categories. In our case, model is the same for all categories
    for(int i = 0; i<ncut; i++){
        TString category = "cat"+std::to_string(i);
        Int_t entries = h_tripletmass[i]->GetEntries(); //reference for normalisation variables
        nsig2[i] = new RooRealVar("nsig2_"+category,"#signal2 events",int(entries/4),5,int(entries/2)); 
        nsig1[i] = new RooRealVar("nsig1_"+category,"#signal1 events",int(entries/4),5,int(entries/2)); 
        nbkg[i] = new RooRealVar("nbkg_"+category,"#background events",int(entries/2),5,entries);
        model[i] = new RooAddPdf("model_"+category,"g+a",RooArgList(signal_CB1,signal_CB2,bg_exp),RooArgList(*nsig1[i],*nsig2[i],*nbkg[i]));
        simPdf.addPdf(*model[i], category);
    }
    // P e r f o r m   a   s i m u l t a n e o u s   f i t
    // ---------------------------------------------------
    // Perform simultaneous fit of model to data and model_ctl to data_ctl
    RooFitResult * r = simPdf.fitTo(dh, Save(true));
    cout<<"Fit results"<<endl;
    r->floatParsFinal().Print("s");

    TCanvas *c5 = new TCanvas("c5","c5",150,10,800,800);
    c5->SetLeftMargin(0.15);
    c5->Update();

    /*
    x.setRange("signal",1.93,2.01);
    x.setRange("sideband",1.72,1.8);
    */
    x.setRange("signal",1.93,2.01);
    x.setRange("sideband",1.72,1.8);

    //for each category: do plot and store integrals:
    for(int i = 0; i<ncut; i++){
        cout<<"Starting step "<< i <<endl;
        
        TString category = "cat"+std::to_string(i);
        // Make plot of binned dataset showing Poisson error bars (RooFit default)
        RooPlot* frame = x.frame(Title(" "));
        dh.plotOn(frame, Cut("c==c::"+category), DataError(RooAbsData::SumW2), Name("data_"+category));
        simPdf.plotOn(frame, Slice(c, category), ProjWData(c, dh), Name("model_"+category), Range("chi2"));
        simPdf.plotOn(frame, Slice(c, category), Components(bg_exp), LineColor(kGreen), LineStyle(kDashed), ProjWData(c, dh));
        //simPdf.plotOn(frame, Slice(c, category), Components(signal_CB2, signal_CB1), LineColor(kRed), LineStyle(kDashed), ProjWData(c, dh));
        //simPdf.plotOn(frame, Slice(c, category), Components(RooArgSet(signal_CB2, signal_CB1)), LineColor(kRed), LineStyle(kDashed), ProjWData(c, dh));
        simPdf.plotOn(frame, Slice(c, category), Components(RooArgSet(signal_CB1)), LineColor(kRed), LineStyle(kDashed), ProjWData(c, dh));
        simPdf.plotOn(frame, Slice(c, category), Components(RooArgSet(signal_CB2)), LineColor(kPink), LineStyle(kDashed), ProjWData(c, dh));
        //simPdf.paramOn(frame,Layout(0.12, 0.4, 0.9));
        frame->Draw();
        
        //add Lumi and Chi2 to plot
        std::stringstream stream_lumi;
        stream_lumi << std::fixed << std::setprecision(1) << lumi_full;
        TString strLumi = stream_lumi.str();
        TLatex* text_lumi = new TLatex(0.10,0.91, "\n\\text{data }"+run_lable+"\n\\text{    }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
        text_lumi->SetTextSize(0.05);
        text_lumi->SetNDC(kTRUE);
        text_lumi->Draw("same");
        
        cout<<"frame->chiSquare() "<<frame->chiSquare("model_"+category, "data_"+category, r->floatParsFinal().getSize())<<endl;
        
        Double_t Chi2 = frame->chiSquare("model_"+category, "data_"+category, r->floatParsFinal().getSize());
        std::stringstream stream_chi2;
        stream_chi2 << std::fixed << std::setprecision(2) << Chi2;
        std::string strChi2 = stream_chi2.str();
        TString chi2tstring = "\\chi^{2}\\text{/NDOF} = "+strChi2;
        TLatex* text_chi2 = new TLatex(0.20,0.74, chi2tstring);
        text_chi2->SetTextSize(0.04);
        text_chi2->SetNDC(kTRUE);
        text_chi2->Draw("same");
        

        TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");

        fout->WriteObject(c5,category);
        c5->SaveAs("dsphipi_fit_perMVA_new"+category+".png");

        //fraction of total events in 1.93,2.01 (n_signal_region_events/n_total_events)
        RooAbsReal* fsigregion_model = model[i]->createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fs = fsigregion_model->getVal();
        Double_t fs_err = fsigregion_model->getPropagatedError(*r);

        //fraction of background events in 1.93,2.01
        RooAbsReal* fsigregion_bkg = ((RooAbsPdf*)model[i]->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fb = fsigregion_bkg->getVal();
        Double_t fb_err = fsigregion_bkg->getPropagatedError(*r);

        Dsyield_data[i] = fs * (nsig2[i]->getVal()+nsig1[i]->getVal()+nbkg[i]->getVal()) - fb*nbkg[i]->getVal();;
        Dsyield_data_err[i] = pow( pow(fs_err,2) * pow(nsig2[i]->getVal()+nsig1[i]->getVal()+nbkg[i]->getVal(),2)  + ( pow(nsig2[i]->getPropagatedError(*r),2)+pow(nsig1[i]->getPropagatedError(*r),2)+pow(nbkg[i]->getPropagatedError(*r),2)) * pow(fs,2) + pow(fb_err,2) * pow(nbkg[i]->getVal(),2) + pow(nbkg[i]->getPropagatedError(*r),2)*pow(fb,2) , 0.5);;

        Bkgyield_data[i] = fb*nbkg[i]->getVal();
        Bkgyield_data_err[i] = sqrt( pow(fb_err, 2.0) * pow(nbkg[i]->getVal(), 2.0) + pow(fb, 2.0) * pow(nbkg[i]->getPropagatedError(*r), 2.0));
    }

    //drawing triplet mass in MC for region mass selection and integral
    TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
    h_tripletmass_mc_full->Draw();
    auto f1  = new TF1("f1","gaus",1.93,2.01);
    h_tripletmass_mc_full->Fit("f1", "R");
    f1->Draw("same");
    c3->Update();
    fout->WriteObject(c3,"3glb_invmass_mc");
    //Integrals MC
    Double_t n_mc_peak = f1->Integral(1.93, 2.01) / h_tripletmass_mc_full->Integral(h_tripletmass_mc_full->FindFixBin(1.93),h_tripletmass_mc_full->FindFixBin(2.01),"width") * h_tripletmass_mc_full->Integral(h_tripletmass_mc_full->FindFixBin(1.93),h_tripletmass_mc_full->FindFixBin(2.01));

    //After fitting, set parameters to constant
    mass1.setConstant(kTRUE);
    mass2.setConstant(kTRUE);
    sigma1.setConstant(kTRUE);
    sigma2.setConstant(kTRUE);
    n1.setConstant(kTRUE);
    n2.setConstant(kTRUE);
    alpha1.setConstant(kTRUE);
    alpha2.setConstant(kTRUE);
    a.setConstant(kTRUE);
    //create model for full dataset
    Int_t entries = h_tripletmass_full->GetEntries(); //reference for normalisation variables
    RooRealVar* nsig2_full = new RooRealVar("nsig2_full","#signal2 events",int(entries/4),5,int(entries/2));
    RooRealVar* nsig1_full = new RooRealVar("nsig1_full","#signal1 events",int(entries/4),5,int(entries/2));
    RooRealVar* nbkg_full = new RooRealVar("nbkg_full","#background events",int(entries/2),5,entries);
    RooAddPdf* model_full = new RooAddPdf("model_full","model_full",RooArgList(signal_CB1,signal_CB2,bg_exp),RooArgList(*nsig1_full,*nsig2_full,*nbkg_full));
    //fit full dataset
    RooDataHist dh_full("dh_full", "dh_full", x, Import(*h_tripletmass_full));
    RooFitResult * r_full = model_full->fitTo(dh_full, Save(true));
    
    // Make plot of binned dataset showing Poisson error bars (RooFit default)
    c5->cd();
    RooPlot* frame_full = x.frame(Title(" "));
    dh_full.plotOn(frame_full, Name("dh_full"));
    model_full->plotOn(frame_full, Name("model_full"));
    model_full->plotOn(frame_full, Components(bg_exp), LineColor(kGreen), LineStyle(kDashed));
    //model_full->plotOn(frame_full, Components(RooArgSet(signal_CB2, signal_CB1)), LineColor(kRed), LineStyle(kDashed) );
    model_full->plotOn(frame_full, Components(RooArgSet(signal_CB1)), LineColor(kRed), LineStyle(kDashed) );
    model_full->plotOn(frame_full, Components(RooArgSet(signal_CB2)), LineColor(kPink), LineStyle(kDashed) );
    frame_full->Draw();
    //Chi2
    cout<<"frame_full->chiSquare() "<<frame_full->chiSquare("model_full", "dh_full", r_full->floatParsFinal().getSize())<<endl;
    RooChi2Var chi2("chi2","chi2",*model_full,dh_full);
    int NDOF = n_bins-r_full->floatParsFinal().getSize();
    cout << "NDOF = " << n_bins << "-" <<r_full->floatParsFinal().getSize()<<" = "<<NDOF<<endl;
    cout << "chi2.getVal()/NDOF = " << chi2.getVal()/NDOF << endl ;
        
    //add Lumi and Chi2 to plot
    std::stringstream stream_lumi;
    stream_lumi << std::fixed << std::setprecision(1) << lumi_full;
    TString strLumi = stream_lumi.str();
    TLatex* text_lumi = new TLatex(0.10,0.91, "\n\\text{data }"+run_lable+"\n\\text{    }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
    text_lumi->SetTextSize(0.05);
    text_lumi->SetNDC(kTRUE);
    text_lumi->Draw("same");

    Double_t Chi2 = chi2.getVal()/NDOF;
    std::stringstream stream_chi2;
    stream_chi2 << std::fixed << std::setprecision(2) << Chi2;
    std::string strChi2 = stream_chi2.str();
    TString chi2tstring = "\\chi^{2}\\text{/NDOF} = "+strChi2;
    TLatex* text_chi2 = new TLatex(0.20,0.74, chi2tstring);

    TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
    cmslabel->SetNDC(kTRUE);
    cmslabel->Draw("same");

    fout->WriteObject(c5,"full");
    c5->SaveAs("dsphipi_fit_perMVA_new_full.png");

    //Integrals data
    //fraction of total events in 1.93,2.01 (n_signal_region_events/n_total_events)
    RooAbsReal* fsigregion_model = model_full->createIntegral(x,NormSet(x),Range("signal")); 
    Double_t fs = fsigregion_model->getVal();
    Double_t fs_err = fsigregion_model->getPropagatedError(*r_full);

    //fraction of background events in 1.93,2.01
    RooAbsReal* fsigregion_bkg = ((RooAbsPdf*)model_full->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("signal")); 
    Double_t fb = fsigregion_bkg->getVal();
    Double_t fb_err = fsigregion_bkg->getPropagatedError(*r);

    Double_t Dsyield_data_full = fs * (nsig2_full->getVal()+nsig1_full->getVal()+nbkg_full->getVal()) - fb*nbkg_full->getVal();;
    Double_t Dsyield_data_err_full = pow( pow(fs_err,2) * pow(nsig2_full->getVal()+nsig1_full->getVal()+nbkg_full->getVal(),2)  + ( pow(nsig2_full->getPropagatedError(*r),2)+pow(nsig1_full->getPropagatedError(*r),2)+pow(nbkg_full->getPropagatedError(*r),2)) * pow(fs,2) + pow(fb_err,2) * pow(nbkg_full->getVal(),2) + pow(nbkg_full->getPropagatedError(*r),2)*pow(fb,2) , 0.5);;

    Double_t Bkgyield_data_full = fb*nbkg_full->getVal();
    Double_t Bkgyield_data_err_full = sqrt( pow(fb_err, 2.0) * pow(nbkg_full->getVal(), 2.0) + pow(fb, 2.0) * pow(nbkg_full->getPropagatedError(*r), 2.0));
    //cout<<"n_mc_peak "<<n_mc_peak<<endl;
    Double_t Dsyield_MC_full = n_mc_peak*lumi_full*xsection_mc*BR/N_MC;
    Double_t Dsyield_MC_err_full = sqrt(n_mc_peak)*(lumi_full*xsection_mc*BR/N_MC); 

    cout<<"\n"<<run_lable<<" lumi="<<lumi_full<<endl;
    cout<<"=================\noverall data/MC scale factor:"<<endl;
    cout<<"data: "<<Dsyield_data_full<<" +- "<<Dsyield_data_err_full<<"\n";
    cout<<"MC: "<<Dsyield_MC_full<<" +- "<<Dsyield_MC_err_full<<endl;
    cout<<"scale factor: "<<Dsyield_data_full/Dsyield_MC_full<<" +- "<<sqrt(pow((Dsyield_data_err_full/Dsyield_MC_full),2.0) + pow((Dsyield_data_full/(Dsyield_MC_full*Dsyield_MC_full))*Dsyield_MC_err_full,2.0))<<endl;
   
    //overall SF to be applied
    Double_t SF_total = Dsyield_data_full/Dsyield_MC_full; 

    for(int i = 0; i<ncut; i++){
        Dsyield_MC[i] = SF_total*h_tripletmass_mc[i]->GetEntries()*lumi_full*xsection_mc*BR/N_MC;
        Dsyield_MC_err[i] = sqrt(SF_total*h_tripletmass_mc[i]->GetEntries())*(lumi_full*xsection_mc*BR/N_MC); 

        fout_yield<<Dsyield_data[i]<<"\t"<<Dsyield_data_err[i]<<"\n";

        cout<<"\n"<<bdt_cutlist[i]<<" lumi="<<lumi_full<<endl;
        cout<<"=================\noverall data/MC scale factor:"<<endl;
        cout<<"data: "<<Dsyield_data[i]<<" +- "<<Dsyield_data_err[i]<<"\n";
        cout<<"MC: "<<Dsyield_MC[i]<<" +- "<<Dsyield_MC_err[i]<<endl;
        double SF = Dsyield_data[i]/Dsyield_MC[i];
        double SF_err = sqrt(pow((Dsyield_data_err[i]/Dsyield_MC[i]),2.0) + pow((Dsyield_data[i]/(Dsyield_MC[i]*Dsyield_MC[i]))*Dsyield_MC_err[i],2.0));
        cout<<"scale factor: "<<SF<<" +- "<<SF_err<<endl;
        if(i>0&&i<10&&which_etabin==1){
          SF_mva1->SetBinContent(i,SF);
          SF_mva1->SetBinError(i,SF_err);
        }
        if(i>0&&i<10&&which_etabin==2){
          SF_mva2->SetBinContent(i,SF);
          SF_mva2->SetBinError(i,SF_err);
        }
    }
    fout->Close();
    
    fout0->cd();
    SF_mva1->Write("",TObject::kOverwrite);
    SF_mva2->Write("",TObject::kOverwrite);
    //SF_eta_mva->Write();
    fout0->Close();
    return;
 
}
