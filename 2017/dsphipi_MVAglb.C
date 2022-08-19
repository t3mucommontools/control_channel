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
    gErrorIgnoreLevel = kInfo;

    //open input rootfile
    TFile *fin = new TFile("/eos/user/f/fsimone/input_files/t3mminitree_xgb_2017_13may22_ptmu3_control.root", "READ");
    cout<<"opened input file t3mminitree_xgb_2017_13may22_ptmu3_control.root"<<endl;
    TTree *tin = (TTree*)fin->Get("OutputTree");
    TString var = "MVA_glbmu2";

    //IMPORTANT: event selection 
    TString common_cut = "bs_sv_d2Dsig>2.0 && mu3_pt > 1.2 && " 
                         "((mu1_pt>3.5 && mu1_eta<1.2) || (mu1_pt>2.0 && mu1_eta>=1.2 && mu1_eta<=2.4)) && "
                         "((mu2_pt>3.5 && mu2_eta<1.2) || (mu2_pt>2.0 && mu2_eta>=1.2 && mu2_eta<=2.4)) && "
                         "abs(phiMass-1.02)<0.045 && "
                         "l1double_DoubleMu0_fired";

    //first step: set optimal binning for MVA score to cut on
    TH1F *h_MVAmu2_prelim;
    TString binning_MVA = "(50, 0.1, 0.8)";
    tin->Draw("MVA_glbmu2>>h_MVAmu2_prelim"+binning_MVA, "("+common_cut+")");
    h_MVAmu2_prelim = (TH1F *)gDirectory->Get("h_MVAmu2_prelim");
    Int_t nbins = 8;
    Double_t edges[nbins+1];
    TString edges_s[nbins+1];
    Double_t proba[nbins+1];
    //compute x values corresponding to evenly distributed probabilities
    for (int i=0; i<nbins+1; i++) proba[i] = float(i)/float(nbins);
    h_MVAmu2_prelim->GetQuantiles(nbins+1, edges, proba);

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
    TH1D *SF_mva1 = new TH1D("SF_mva1", "SF vs MVA", ncut, 0.1, 0.8 ) ;
    TH1D *SF_mva2 = new TH1D("SF_mva2", "SF vs MVA", ncut, 0.1, 0.8 ) ;
    TH2D *SF_mva  = new TH2D("SF_mva",  "SF vs MVA-eta", ncut, 0.1, 0.8, 2, 0, 2.4) ;

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
    TString binning_mass = "(62, 1.72, 2.01)";
    Int_t n_bins = 62;
	
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
        TString invmass_all_data  = "puFactor*(tripletMass<2.01 && tripletMass>1.62 && "+common_cut+eta_cut+" && isMC==0";
        TString invmass_peak_MC   = "puFactor*(tripletMass<2.01 && tripletMass>1.93 && "+common_cut+eta_cut+" && isMC==1";

        //Fill histograms for each MVA cut
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
        
        x.setRange("R3",1.72,1.82); //background    
        x.setRange("R4",1.90,1.925); //background    
        x.setRange("R5",2.0,2.01); //background 
        
        // Fit a Crystal ball p.d.f to the data
        RooRealVar mass1("mass1","Central value of Gaussian",1.86774,1.75,2.01);
        RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0001,0.3);
        
        RooRealVar alpha1("alpha1","alpha value CB",1.9,0.1,20);
        RooRealVar n1("n1", "n1", 2, 0, 10);
        RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
        
        // Fit a Crystal ball p.d.f to the data
        RooRealVar mass2("mass2","Central value of Gaussian",1.96724,1.85,2.51);
        RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.3);
        
        RooRealVar alpha2("alpha2","alpha value CB",1.9,0.1,20);
        RooRealVar n2("n2", "n2", 2, 0, 20);
        RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
           
        // Fit an exponential to the background
        RooRealVar a("a", "a", -5, -30, 5.0);
        RooExponential bg_exp("bg_exp", "bg_exp", x, a);
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

        //for each bin: do plot and store integrals:
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
            TLatex* bin_label = new TLatex(0.20,0.64, bdt_cutlabels[i]+", "+eta_region);
            bin_label->SetNDC(kTRUE);
            bin_label->SetTextSize(0.04);
            bin_label->SetNDC(kTRUE);
            bin_label->Draw("same");
            
            //add CMS label to plot
            TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
            cmslabel->SetNDC(kTRUE);
            cmslabel->Draw("same");

            //save plot
            fout->WriteObject(c5,category+"_"+eta_region);
            c5->SaveAs("plots_MVAcut/dsphipi_fit_perMVA_new"+category+"_"+eta_region+".png");

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
        bin_label->SetNDC(kTRUE);
        bin_label->Draw("same");

        //add CMS label to plot
        TLatex* cmslabel = new TLatex(0.20,0.81, "#bf{CMS Preliminary}");
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");

        fout->WriteObject(c5,"full_"+eta_region);
        c5->SaveAs("plots_MVAcut/dsphipi_fit_perMVA_new_full_"+eta_region+".png");

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
    }
    fout->cd();
    fout->WriteObject(SF_mva1,"SF_mva1");
    fout->WriteObject(SF_mva2,"SF_mva2");
    fout->WriteObject(SF_mva,"SF_mva");
    fout->Close();
    
    return;
 
}
