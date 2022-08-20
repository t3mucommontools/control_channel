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

struct sfMuon{
   float value=1.0;
   float error=0.0;
};

sfMuon GetMuonSF(const TH2F* _h, const double x, const double eta ){
    sfMuon val;
    if( x < _h->GetXaxis()->GetXmax() && x > _h->GetXaxis()->GetXmin() &&
        eta < _h->GetYaxis()->GetXmax() && eta > _h->GetYaxis()->GetXmin()){
        int ix = _h->GetXaxis()->FindBin(x);
        int ieta = _h->GetYaxis()->FindBin(std::abs(eta));
        val.value = _h->GetBinContent(ix,ieta);
        val.error = _h->GetBinError(ix,ieta);
        //cout<<"x="<<x<<" y="<<std::abs(eta)<<" ix="<<ix<<" iy="<<ieta<<endl;
    }else{
        val.value = 1;
        val.error = 0;
    }
    return val;
}

void dsphipi_MVAglb() 
{
    gErrorIgnoreLevel = kInfo;

    //open input rootfile
    TFile *fin = new TFile("/eos/user/f/fsimone/input_files/t3mminitree_xgb_2017_13may22_ptmu3_control.root", "READ");
    cout<<"opened input file t3mminitree_xgb_2017_13may22_ptmu3_control.root"<<endl;
    TTree *tin = (TTree*)fin->Get("OutputTree");
    TString var = "MVA_glbmu2";

    //IMPORTANT: event selection 
    TString common_cut = "bs_sv_d2Dsig>3.75 && mu3_pt > 1.2 && " 
                         "((mu1_pt>3.5 && mu1_eta<1.2) || (mu1_pt>2.0 && mu1_eta>=1.2 && mu1_eta<=2.4)) && "
                         "((mu2_pt>3.5 && mu2_eta<1.2) || (mu2_pt>2.0 && mu2_eta>=1.2 && mu2_eta<=2.4)) && "
                         "abs(phiMass-1.02)<0.045 && "
                         "l1double_DoubleMu0_fired";

    //first step: set optimal binning for MVA score to cut on
    TH1F *h_MVAmu2_prelim;
    TString binning_MVA = "(80, 0.1, 0.8)";
    tin->Draw("MVA_glbmu2>>h_MVAmu2_prelim"+binning_MVA, "("+common_cut+")");
    h_MVAmu2_prelim = (TH1F *)gDirectory->Get("h_MVAmu2_prelim");
    Int_t nbins = 15;
    Double_t edges[nbins+1];
    TString edges_s[nbins+1];
    Double_t proba[nbins+1];
    //compute x values corresponding to evenly distributed probabilities
    for (int i=0; i<nbins+1; i++) proba[i] = float(i)/float(nbins);
    h_MVAmu2_prelim->GetQuantiles(nbins+1, edges, proba);

    //edges to string with desired precision
    for (int i=0; i<nbins+1; i++){
        stringstream ss;
        ss << setprecision(3)<< edges[i];
        edges_s[i] = ss.str();
    }
    //set here list of MVA cuts
    TString bdt_cutlist[nbins];
    TString bdt_cutlabels[nbins];
    for (int i=0; i<nbins; i++){
        bdt_cutlist[i] = var+"<="+edges[i+1]+"&&"+var+">"+edges[i];
        cout<<bdt_cutlist[i]<<endl;
        bdt_cutlabels[i] = "MVAglb_mu2 ["+edges_s[i+1]+"; "+edges_s[i]+"]";
        cout<<bdt_cutlabels[i]<<endl;
    }
    int ncut = sizeof(bdt_cutlist)/sizeof(bdt_cutlist[0]);
    cout<<"Size of ncut: "<<ncut<<endl;
 
    //set here list of eta cuts
    TString eta_cutlist[2]; //0: barrel, 1:endcap
    eta_cutlist[0] = " && abs(mu2_eta)<1.2";
    eta_cutlist[1] = " && abs(mu2_eta)>=1.2";
    TString eta_regions[] = {"barrel", "endcap"};

    //create output rootfile
    TFile *fout = new TFile("plots_MVAcut/dsphipi_MVAglb_2017_output.root", "RECREATE");
    fout->cd();
    //to do: add boundaries to MVAmu2 plot
    h_MVAmu2_prelim->Write();
    //output histograms
    TH1F *SF_mva1 = new TH1F("SF_mva1", "SF vs MVA", ncut, edges) ;
    TH1F *SF_mva2 = new TH1F("SF_mva2", "SF vs MVA", ncut, edges) ;
    TH2F *SF_mva  = new TH2F("SF_mva",  "SF vs MVA-eta", ncut, edges, 2, 0, 2.4) ;

    //yields     
    Double_t Dsyield_data[ncut];
    Double_t Dsyield_data_err[ncut];
    Double_t Dsyield_MC[ncut];
    Double_t Dsyield_MC_err[ncut];
    Double_t Bkgyield_data[ncut];
    Double_t Bkgyield_data_err[ncut];
    
    //histograms     
    TH1F *h_tripletmass[ncut];//all data range with mva cut
    TH1F *h_tripletmass_mc[ncut];//peak MC range with mva cut

    TH1F *h_tripletmass_full;
    TH1F *h_tripletmass_mc_full;

    //binning and normalisation quantities
    Double_t events_SB = 0;
    TString binning_mass = "(58, 1.72, 2.01)";
    Int_t n_bins = 58;
	
    Double_t lumi_full = 37.9; //fb
    TString run_lable = "2017";
    Double_t xsection_mc = 1.06e10; //Ds Production Cross section
    int N_MC = 2440762;  //Total number of events in MC sample
    Double_t BR = 1.29e-5;  //Branching ratio Ds to Phi Pi

    //loop on eta regions
    for(int j=0; j<2; j++){
        TString eta_cut = eta_cutlist[j];
        TString eta_region = eta_regions[j];

        //event selection, mass cuts and eta region
        TString invmass_all_data  = "puFactor*(tripletMass<2.01 && tripletMass>1.72 && "+common_cut+eta_cut+" && isMC==0";
        TString invmass_SB_data  = "puFactor*(tripletMass<1.82 && tripletMass>1.72 && "+common_cut+eta_cut+" && isMC==0";
        TString invmass_peak_data  = "puFactor*(tripletMass<2.01 && tripletMass>1.93 && "+common_cut+eta_cut+" && isMC==0";
        TString invmass_peak_MC   = "puFactor*(tripletMass<2.01 && tripletMass>1.93 && "+common_cut+eta_cut+" && isMC==1";

        //Fill histograms for each MVA cut
        for(int i=0; i<ncut; i++){
            TString s = std::to_string(i);
            
            TString bdt_cut = bdt_cutlist[i];
        		
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
       
        x.setRange("R1",1.83,1.89); //first peak D+(1.87GeV)
        x.setRange("R2",1.93,2.01); //second peak Ds(1.97) 
        x.setRange("R3",1.72,1.82); //background    
        x.setRange("R4",1.90,1.92); //background    
        x.setRange("R5",2.0,2.01); //background 
        
        // Fit a Crystal ball p.d.f to the data
        RooRealVar mass1("mass1","Central value of Gaussian",1.87,1.86,1.88);
        RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0001,0.6);
        
        RooRealVar alpha1("alpha1","alpha value CB",1.0,0.01,10);
        RooRealVar n1("n1", "n1", 2.0, 0.5, 10);
        RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
        RooGaussian signal_G1("g1", "The signal distribution", x, mass1, sigma1); 
        
        // Fit a Crystal ball p.d.f to the data
        RooRealVar mass2("mass2","Central value of Gaussian",1.965,1.955,1.975);
        RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.6);
        
        RooRealVar alpha2("alpha2","alpha value CB",1.0,0.1,10);
        RooRealVar n2("n2", "n2", 2.0, 0.5, 10);
        RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
        RooGaussian signal_G2("g2", "The signal distribution", x, mass2, sigma2); 
           
        // Fit an exponential to the background
        RooRealVar a("a", "a", -2, -10, 0.0);
        RooExponential bg_exp("bg_exp", "bg_exp", x, a);
        //bg_exp.fitTo(dh, Range("R3,R4"));
        bg_exp.fitTo(dh, Range("R3,R4,R5"));
        //bg_exp.fitTo(dh, Range("R3"));

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
            nsig2[i] = new RooRealVar("nsig2_"+category,"#signal2 events",int(entries/6),5,int(entries/2)); 
            nsig1[i] = new RooRealVar("nsig1_"+category,"#signal1 events",int(entries/6),5,int(entries/2)); 
            nbkg[i] = new RooRealVar("nbkg_"+category,"#background events",int(entries/2),5,entries);
            //model[i] = new RooAddPdf("model_"+category,"g+a",RooArgList(signal_CB1,signal_CB2,bg_exp),RooArgList(*nsig1[i],*nsig2[i],*nbkg[i]));
            model[i] = new RooAddPdf("model_"+category,"g+a",RooArgList(signal_G1,signal_G2,bg_exp),RooArgList(*nsig1[i],*nsig2[i],*nbkg[i]));
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

        x.setRange("signal",1.93,2.01);
        x.setRange("sideband",1.72,1.8);

        //for each bin: do plot and store integrals:
        for(int i = 0; i<ncut; i++){
            TString category = "cat"+std::to_string(i);
            // Make plot of binned dataset showing Poisson error bars (RooFit default)
            RooPlot* frame = x.frame(Title(" "));
            dh.plotOn(frame, Cut("c==c::"+category), DataError(RooAbsData::SumW2), Name("data_"+category));
            simPdf.plotOn(frame, Slice(c, category), ProjWData(c, dh), Name("model_"+category), Range("chi2"));
            simPdf.plotOn(frame, Slice(c, category), Components(bg_exp), LineColor(kGreen), LineStyle(kDashed), ProjWData(c, dh));
            simPdf.plotOn(frame, Slice(c, category), Components(RooArgSet(signal_G1)), LineColor(kRed), LineStyle(kDashed), ProjWData(c, dh));
            simPdf.plotOn(frame, Slice(c, category), Components(RooArgSet(signal_G2)), LineColor(kRed), LineStyle(kDashed), ProjWData(c, dh));
            //simPdf.paramOn(frame,Layout(0.12, 0.4, 0.9));
            frame->Draw();
            
            //add Lumi to plot
            std::stringstream stream_lumi;
            stream_lumi << std::fixed << std::setprecision(1) << lumi_full;
            TString strLumi = stream_lumi.str();
            TLatex* text_lumi = new TLatex(0.10,0.91, "\n\\text{data }"+run_lable+"\n\\text{    }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
            text_lumi->SetTextSize(0.05);
            text_lumi->SetNDC(kTRUE);
            text_lumi->Draw("same");
            
            //add Chi2 to plot
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

            //add bin to plot
            TLatex* bin_label = new TLatex(0.20,0.38, "#bf{"+bdt_cutlabels[i]+"}");
            TLatex* eta_label = new TLatex(0.20,0.34, "#bf{"+eta_region+"}");
            bin_label->SetNDC(kTRUE);
            bin_label->SetTextSize(0.04);
            bin_label->Draw("same");
            eta_label->SetNDC(kTRUE);
            eta_label->SetTextSize(0.04);
            eta_label->Draw("same");
            
            //add CMS label to plot
            TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
            cmslabel->SetNDC(kTRUE);
            cmslabel->Draw("same");

            //save plot
            fout->WriteObject(c5,category+"_"+eta_region);
            c5->SaveAs("plots_MVAcut/dsphipi_fit_perMVA_"+category+"_"+eta_region+".png");

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
        fout->WriteObject(c3,"mc_"+eta_region);
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
            
        //add Lumi to plot
        std::stringstream stream_lumi;
        stream_lumi << std::fixed << std::setprecision(1) << lumi_full;
        TString strLumi = stream_lumi.str();
        TLatex* text_lumi = new TLatex(0.10,0.91, "\n\\text{data }"+run_lable+"\n\\text{    }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
        text_lumi->SetTextSize(0.05);
        text_lumi->SetNDC(kTRUE);
        text_lumi->Draw("same");

        //add Chi2 to plot
        Double_t Chi2 = chi2.getVal()/NDOF;
        std::stringstream stream_chi2;
        stream_chi2 << std::fixed << std::setprecision(2) << Chi2;
        std::string strChi2 = stream_chi2.str();
        TString chi2tstring = "\\chi^{2}\\text{/NDOF} = "+strChi2;
        TLatex* text_chi2 = new TLatex(0.20,0.74, chi2tstring);

        //add bin to plot
        TLatex* bin_label = new TLatex(0.20,0.64, eta_region);
        bin_label->SetTextSize(0.05);
        bin_label->SetNDC(kTRUE);
        bin_label->Draw("same");

        //add CMS label to plot
        TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");

        fout->WriteObject(c5,"full_"+eta_region);
        c5->SaveAs("plots_MVAcut/dsphipi_fit_perMVA_full_"+eta_region+".png");

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

            cout<<"\n"<<bdt_cutlabels[i]<<" lumi="<<lumi_full<<endl;
            cout<<"=================\noverall data/MC scale factor:"<<endl;
            cout<<"data: "<<Dsyield_data[i]<<" +- "<<Dsyield_data_err[i]<<"\n";
            cout<<"MC: "<<Dsyield_MC[i]<<" +- "<<Dsyield_MC_err[i]<<endl;
            double SF = Dsyield_data[i]/Dsyield_MC[i];
            double SF_err = sqrt(pow((Dsyield_data_err[i]/Dsyield_MC[i]),2.0) + pow((Dsyield_data[i]/(Dsyield_MC[i]*Dsyield_MC[i]))*Dsyield_MC_err[i],2.0));
            cout<<"scale factor: "<<SF<<" +- "<<SF_err<<endl;
            SF_mva->SetBinContent(i+1,j+1,SF);
            SF_mva->SetBinError(i+1,j+1,SF_err);
            if(j==0){ //barrel
              SF_mva1->SetBinContent(i+1,SF);
              SF_mva1->SetBinError(i+1,SF_err);
            }
            else if(j==1){ //endcap
              SF_mva2->SetBinContent(i+1,SF);
              SF_mva2->SetBinError(i+1,SF_err);
            }
        }
   
        //Control plot
        TH1F *hdata_bkg      =   new TH1F("hdata_bkg", "hdata_bkg", 80, 0.1, 0.8); //ncut, edges for variable binning
        TH1F *hdata_bkg_plus =   new TH1F("hdata_bkg_plus", "hdata_bkg_plus", 80, 0.1, 0.8);
        TH1F *hdata_bkg_minus =  new TH1F("hdata_bkg_minus", "hdata_bkg_minus", 80, 0.1, 0.8);
        TH1F *hdata_sgn =        new TH1F("hdata_sgn", "hdata_sgn", 80, 0.1, 0.8);
        TH1F *hdata_sgn_plus =   new TH1F("hdata_sgn_plus", "hdata_sgn_plus", 80, 0.1, 0.8);
        TH1F *hdata_sgn_minus =  new TH1F("hdata_sgn_minus", "hdata_sgn_minus", 80, 0.1, 0.8);
        TH1F *hmc_sgn =          new TH1F("hmc_sgn", "hmc_sgn", 80, 0.1, 0.8);
        TH1F *hmc_sgn_weighted = new TH1F("hmc_weighted", "hmc_sgn_weighted", 80, 0.1, 0.8);

        tin->Draw(var+">>hdata_bkg", invmass_SB_data+")");
        tin->Draw(var+">>hdata_bkg_plus", invmass_SB_data+")");
        tin->Draw(var+">>hdata_bkg_minus", invmass_SB_data+")");
        tin->Draw(var+">>hdata_sgn", invmass_peak_data+")");
        tin->Draw(var+">>hdata_sgn_plus", invmass_peak_data+")");
        tin->Draw(var+">>hdata_sgn_minus", invmass_peak_data+")");

        tin->Draw(var+">>hmc_sgn", invmass_peak_MC+")");

        hdata_bkg = (TH1F *)gDirectory->Get("hdata_bkg");
        hdata_bkg_plus = (TH1F *)gDirectory->Get("hdata_bkg_plus");
        hdata_bkg_minus = (TH1F *)gDirectory->Get("hdata_bkg_minus");
        hdata_sgn = (TH1F *)gDirectory->Get("hdata_sgn");
        hdata_sgn_plus = (TH1F *)gDirectory->Get("hdata_sgn_plus");
        hdata_sgn_minus = (TH1F *)gDirectory->Get("hdata_sgn_minus");
        hmc_sgn = (TH1F *)gDirectory->Get("hmc_sgn");

        TCanvas *c2 = new TCanvas("c2","c2",150,10,800,800);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);

        Double_t normMC = hmc_sgn->GetEntries();
        //Normalizing Monte Carlo 
        Double_t wNorm = lumi_full*xsection_mc*BR/N_MC;
        //Double_t wNorm = lumi_full*xsection_mc*BR/N_MC  *  n_mc_peak/hmc_sgn->GetEntries();
        cout<<"wNorm = lumi_full*xsection_mc*BR/N_MC = "<<wNorm<<endl;
        hmc_sgn->Scale(wNorm);

        //scaling the SB distribution to number of background events in 1.93,2.01
        Double_t normSB = hdata_bkg->GetEntries();
        hdata_bkg->Scale(Bkgyield_data_full/normSB);
        hdata_bkg_plus->Scale( (Bkgyield_data_full/normSB) * 1.10);
        hdata_bkg_minus->Scale( (Bkgyield_data_full/normSB) * 0.90);

        cout<<"Entries in  hdata_sgn before SB subtraction "<<hdata_sgn->GetEntries()<<endl;
        hdata_sgn->Add(hdata_bkg,-1); //subtract h2 from h1 : h1->Add(h2,-1)
        hdata_sgn_plus->Add(hdata_bkg_plus,-1); //subtract h2 from h1 : h1->Add(h2,-1)
        hdata_sgn_minus->Add(hdata_bkg_minus,-1); //subtract h2 from h1 : h1->Add(h2,-1)

        //Rescaling to same integral
        hmc_sgn->Scale( 1.0 / hmc_sgn->Integral());
        hdata_sgn->Scale( hmc_sgn->Integral()/hdata_sgn->Integral() );
        hdata_sgn_plus->Scale( hmc_sgn->Integral()/hdata_sgn_plus->Integral() );
        hdata_sgn_minus->Scale( hmc_sgn->Integral()/hdata_sgn_minus->Integral() );

        //plot makeup
        double Y_max = std::max(hmc_sgn->GetMaximum(), hdata_sgn->GetMaximum());
        Y_max = Y_max*1.4;
        hmc_sgn->GetYaxis()->SetRangeUser(0, Y_max);
    
        hmc_sgn->GetYaxis()->SetTitle("a.u.");
        hmc_sgn->GetYaxis()->SetTitleSize(22);
        hmc_sgn->GetYaxis()->SetTitleFont(43);
        hmc_sgn->GetYaxis()->SetTitleOffset(1.7);
    
        hmc_sgn->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        hmc_sgn->GetXaxis()->SetTitleOffset(1.5);
        hmc_sgn->GetXaxis()->SetTitle(var);
    
        hmc_sgn->SetLineColor(kBlue);
        hmc_sgn->SetLineWidth(3);
        hmc_sgn->SetFillStyle(3004);
        hmc_sgn->SetFillColor(kBlue);
        hdata_sgn->SetLineColor(kRed);
        hdata_sgn->SetLineWidth(3);
        hdata_sgn->SetFillStyle(3005);
        hdata_sgn->SetFillColor(kRed);
    
        hdata_sgn_plus->SetLineColor(kBlack);
        hdata_sgn_plus->SetFillStyle(3001);
        hdata_sgn_plus->SetLineWidth(2);
        hdata_sgn_plus->SetLineStyle(2);
    
        hdata_sgn_minus->SetLineColor(kBlack);
        hdata_sgn_minus->SetFillStyle(3001);
        hdata_sgn_minus->SetLineWidth(2);
        hdata_sgn_minus->SetLineStyle(3);
 
        hmc_sgn->SetTitle(""); // Remove the title
        hmc_sgn->SetStats(0);

        hmc_sgn->Draw("hist");
        hdata_sgn->Draw("hist same");
        hdata_sgn_plus->Draw("hist same");
        hdata_sgn_minus->Draw("hist same");

        //Double_t x_low = 0.45, x_high = 0.89, y_low = 0.50, y_high = 0.89; //top right 
        Double_t x_low = 0.12, x_high = 0.45, y_low = 0.45, y_high = 0.89; //top left 
        TLegend*leg = new TLegend(x_low, y_low, x_high, y_high);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.040);
        //leg->SetHeader(varlable, "L");
        leg->AddEntry(hmc_sgn,"D_{s}#rightarrow#phi(#mu#mu)#pi MC","f");
        leg->AddEntry(hdata_sgn,"data "+run_lable+" (SB subtracted)","f");
        leg->AddEntry(hdata_sgn_plus,"data "+run_lable+" (SB +10\% subt.)","f");
        leg->AddEntry(hdata_sgn_minus,"data "+run_lable+" (SB -10\% subt.)","f");
        leg->Draw();
        bin_label->Draw();
        c2->Update();
    
        fout->cd();
        fout->WriteObject(c2,var+"_controlplot_before_"+eta_region);

        //weighted MC
        hmc_sgn_weighted = (TH1F*)hmc_sgn->Clone("hmc_sgn_weighted");
        int nbinsx = hmc_sgn->GetXaxis()->GetNbins();
        for(int b=0; b<nbinsx; b++){
            Double_t mva_h = hmc_sgn->GetBinContent(b+1);
            Double_t mva_x = ((TAxis*)hmc_sgn->GetXaxis())->GetBinCenter(b+1);
            sfMuon SF;
            if(j==0) SF = GetMuonSF(SF_mva, mva_x, 0.6); //barrel
            else if(j==1) SF = GetMuonSF(SF_mva, mva_x, 1.8); //endcap
            hmc_sgn_weighted->SetBinContent(b+1, mva_h*SF.value);
            //cout<<"bin "<<b+1<<" mva_x "<<mva_x<<" SF "<<SF.value<<"+-"<<SF.error<<endl;
        }
        //if(j==0) hmc_sgn_weighted->Multiply(SF_mva1);            
        //else if(j==1) hmc_sgn_weighted->Multiply(SF_mva2);            

        //Rescaling to same integral
        hmc_sgn_weighted->Scale( 1.0 / hmc_sgn_weighted->Integral());
        hdata_sgn->Scale( hmc_sgn_weighted->Integral()/hdata_sgn->Integral() );
        hdata_sgn_plus->Scale( hmc_sgn_weighted->Integral()/hdata_sgn_plus->Integral() );
        hdata_sgn_minus->Scale( hmc_sgn_weighted->Integral()/hdata_sgn_minus->Integral() );
        hmc_sgn_weighted->GetYaxis()->SetRangeUser(0, Y_max);

        TCanvas *c6 = new TCanvas("c6","c6",150,10,800,800);

        hmc_sgn_weighted->Draw("hist");
        hdata_sgn->Draw("hist same");
        hdata_sgn_plus->Draw("hist same");
        hdata_sgn_minus->Draw("hist same");
        c6->Update();
        leg->Draw();
        bin_label->Draw();
        c6->Update();
    
        hmc_sgn_weighted->SetStats(0);
        fout->cd();
        fout->WriteObject(c6,var+"_controlplot_after_"+eta_region);

    }
    fout->cd();
    fout->WriteObject(SF_mva1,"SF_mva1");
    fout->WriteObject(SF_mva2,"SF_mva2");
    fout->WriteObject(SF_mva,"SF_mva");

    fout->Close();
    return;
}
