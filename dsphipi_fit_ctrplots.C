#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>
//#include "../Control_common.h"

using namespace RooFit;

//N.B. inputpath_DsPhiPi = path of MC ttree "FinalTree_Control". 2mu+1trk mass is named "tripletMass"
//     inputpath_datarun_control[i] = same for data, where the index i runs on data taking periods (e.g. 2017B, 2017C etc ..)
//     nrun = number of data taking periods

//MC normalization factor:
// lumi*xsection_mc*BR/N_MC
// where: BR = Ds->PhiPi BR
//        N_MC = number of events generated in the MC sample
//        xsection_mc = cross section as given by XGenAnalyser (includes filter efficiencies)

void dsphipi_fit_ctrplots() 
{
    bool doAllRuns = false; //if true, full dataset will be used. if false, dsphipi fit and control plot will be produced separately for each run
    bool doPUrew = true;

    //open root files where to store output plots
    TFile *fout = new TFile("dsphipi_all_output.root", "RECREATE");
    fout->cd();

    TFile *f_mc = new TFile(inputpath_DsPhiPi, "READ");
    TTree *tmc = (TTree*)f_mc->Get("FinalTree_Control");

    int nrun = sizeof(inputpath_datarun_control)/sizeof(inputpath_datarun_control[0]);

    TChain *tdata_all = new TChain("FinalTree_Control");

    for(auto i=0; i<nrun; i++){
        TFile *f = new TFile(inputpath_datarun_control[i],"READ");
        tdata_all->Add(inputpath_datarun_control[i]);
        f->Close();
    }
    Double_t lumi = 0;
    for(int k=0; k<nrun; k++) lumi = lumi+ Lumi_data_control[k];
    if(doAllRuns) nrun = 1;
    TString run_lable = "";

    for(int i=0; i<nrun; i++){
        TFile *f = new TFile(inputpath_datarun_control[i],"READ");
        TTree *tdata = (TTree*)f->Get("FinalTree_Control");

        TString invmass_SB   = "";
        TString invmass_peak = "";
        if(doPUrew){
        invmass_SB   = "puFactor*(tripletMass<1.80 && tripletMass>1.70)";
        invmass_peak = "puFactor*(tripletMass<2.01 && tripletMass>1.93)";
        } else {
        invmass_SB   = "(tripletMass<1.80 && tripletMass>1.70)";
        invmass_peak = "(tripletMass<2.01 && tripletMass>1.93)";
        }
        TString binning_mass = "(72, 1.65, 2.01)";

        TH1F *h_tripletmass_mc;
        TH1F *h_tripletmass;
        TH1F *h_tripletmass_sign;
        TH1F *h_tripletmass_bkg;

        if(doAllRuns){
            run_lable = "2017";
            tdata_all->Draw("tripletMass>>h_tripletmass"+binning_mass);
            tdata_all->Draw("tripletMass>>h_tripletmass_bkg"+binning_mass, invmass_SB);
            tdata_all->Draw("tripletMass>>h_tripletmass_sign"+binning_mass, invmass_peak);
        }else{
            run_lable = run_name_control[i];
            lumi = Lumi_data_control[i];
            tdata->Draw("tripletMass>>h_tripletmass"+binning_mass);
            tdata->Draw("tripletMass>>h_tripletmass_bkg"+binning_mass, invmass_SB);
            tdata->Draw("tripletMass>>h_tripletmass_sign"+binning_mass, invmass_peak);
        }
        tmc->Draw("tripletMass>>h_tripletmass_mc"+binning_mass, invmass_peak);
    
        h_tripletmass     = (TH1F *)gDirectory->Get("h_tripletmass");
        h_tripletmass_bkg = (TH1F *)gDirectory->Get("h_tripletmass_bkg");
        h_tripletmass_sign = (TH1F *)gDirectory->Get("h_tripletmass_sign");
        h_tripletmass_mc = (TH1F *)gDirectory->Get("h_tripletmass_mc");


        //drawing triplet mass in MC for region mass selection and integral
        TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
        h_tripletmass_mc->Draw();
        auto f1  = new TF1("f1","gaus",1.93,2.01);
        h_tripletmass_mc->Fit("f1", "R");
        f1->Draw("same");
        
        Double_t n_mc_peak = f1->Integral(1.93, 2.01) / h_tripletmass_mc->Integral(h_tripletmass_mc->FindFixBin(1.93),h_tripletmass_mc->FindFixBin(2.01),"width") * h_tripletmass_mc->Integral(h_tripletmass_mc->FindFixBin(1.93),h_tripletmass_mc->FindFixBin(2.01));
        cout<<"n_mc_peak "<<n_mc_peak<<endl;
        c3->Update();
        if(i==0) fout->WriteObject(c3,"2mu1trk_invmass_mc");

        TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
        h_tripletmass->Draw();
        c1->Update();
        h_tripletmass_sign->Draw("same");
        h_tripletmass_bkg->Draw("same");
        c1->Update();

        // Declare observable x
        TCanvas *c5 = new TCanvas("c5","c5",150,10,990,660);
        c5->Update();
        RooRealVar x("x","2mu+1trk inv. mass (GeV)",1.65,2.01);

        // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
        RooDataHist dh("dh","dh",x,Import(*h_tripletmass)); 

        // Make plot of binned dataset showing Poisson error bars (RooFit default)
        RooPlot* frame = x.frame(Title(" "));
        dh.plotOn(frame);
        //dh.statOn(frame,Layout(0.55,0.99,0.8)) ;

        //set ranges
        x.setRange("R1",1.83,1.89); //first peak D+(1.87GeV)
        x.setRange("R2",1.93,2.01); //second peak Ds(1.97)
        x.setRange("R3",1.65,1.84); //background    
        x.setRange("R4",1.89,1.925); //background    
        x.setRange("R5",1.99,2.01); //background    
        x.setRange("R6",1.65,2.01); //full range    

        RooRealVar mass2("mass2","Central value of Gaussian",1.965,1.94,2.0);
        RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0.0001,0.06);
        RooRealVar mass1("mass1","Central value of Gaussian",1.87,1.84,1.89);
        RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0.0001,0.1);


        // Fit a Crystal ball p.d.f to the data
        RooRealVar alpha2("alpha2","alpha value CB",1,-20,20);
        RooRealVar n2("n2", "n2", 2, 0, 20);
        RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
        signal_CB2.fitTo(dh, Range("R2"));
        
        
        // Fit a Crystal ball p.d.f to the data
        RooRealVar alpha1("alpha1","alpha value CB",1,-20,20);
        RooRealVar n1("n1", "n1", 2, 0, 5);
        RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
        signal_CB1.fitTo(dh, Range("R1"));
       

        // Fit an exponential to the background
        Double_t a_high = 0.;
        //if(i==0) a_high = 1.; 
        RooRealVar a("a", "a", -5, -20, a_high);
        RooExponential bg_exp("bg_exp", "bg_exp", x, a);
        bg_exp.fitTo(dh, Range("R3,R4,R5"));
        //bg_exp.fitTo(dh, Range("R3,R4"));
        //bg_exp.fitTo(dh, Range("R3"));


        // Combine the models
        RooRealVar nsig2("nsig2","#signal2 events",80000,500,5000000);
        RooRealVar nsig1("nsig1","#signal1 events",40000,500,5000000);
        RooRealVar nbkg("nbkg","#background events",400000,500,10000000);
        RooAddPdf model("model","g+a",RooArgList(signal_CB1,signal_CB2,bg_exp),RooArgList(nsig1,nsig2,nbkg));
        RooFitResult * r = model.fitTo(dh, Save(true));
        r->Print();
        model.paramOn(frame,Layout(0.12, 0.4, 0.9));

        //plot
        model.plotOn(frame);
        model.plotOn(frame, Components(bg_exp), LineColor(kGreen), LineStyle(kDashed));
        model.plotOn(frame, Components(RooArgSet(signal_CB2, signal_CB1)), LineColor(kRed), LineStyle(kDashed) );
        frame->Draw();
        //Chi2
        cout<<"frame->chiSquare(3) "<<frame->chiSquare(3)<<endl;
        RooChi2Var chi2("chi2","chi2",model,dh);
        cout << "chi2.getVal() " << chi2.getVal() << endl ;

        //Compute integrals
        x.setRange("signal",1.93,2.01);
        x.setRange("sideband",1.7,1.8);

        //fraction of total events in 1.93,2.01 (n_signal_region_events/n_total_events)
        RooAbsReal* fsigregion_model = model.createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fs = fsigregion_model->getVal();
        Double_t fs_err = fsigregion_model->getPropagatedError(*r);
        //fraction of total events in 1.70,1.80 (n_sideband_region_events/n_total_events)
        RooAbsReal* fsidebandregion_model = model.createIntegral(x,NormSet(x),Range("sideband")); 

        //fraction of background events in 1.93,2.01 
        RooAbsReal* fsigregion_bkg = bg_exp.createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fb = fsigregion_bkg->getVal();
        Double_t fb_err = fsigregion_bkg->getPropagatedError(*r);
        //fraction of background events in 1.70, 1.80 
        RooAbsReal* fsidebandregion_bkg = bg_exp.createIntegral(x,NormSet(x),Range("sideband")); 


        Double_t nsigevents = fs * (nsig2.getVal()+nsig1.getVal()+nbkg.getVal()) - fb*nbkg.getVal(); 
        Double_t nsig_err = pow( pow(fs_err,2) * pow(nsig2.getVal()+nsig1.getVal()+nbkg.getVal(),2)  + ( pow(nsig2.getPropagatedError(*r),2)+pow(nsig1.getPropagatedError(*r),2)+pow(nbkg.getPropagatedError(*r),2)) * pow(fs,2) + pow(fb_err,2) * pow(nbkg.getVal(),2) + pow(nbkg.getPropagatedError(*r),2)*pow(fb,2) , 0.5);

        Double_t fsig = nsigevents/(fsigregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal()));

        //cout<<"fsigregion_model "<<  fsigregion_model->getVal()<<endl;
        //cout<<"fsigregion_bkg "<<  fsigregion_bkg->getVal()<<endl;

        //cout<<"fsidebandregion_model "<<  fsidebandregion_model->getVal()<<endl;
        //cout<<"fsidebandregion_bkg "<<  fsidebandregion_bkg->getVal()<<endl;

        cout<<"n background events in 1.70,1.80 "<< fsidebandregion_bkg->getVal()*nbkg.getVal() <<endl;
        cout<<"n total events in 1.70,1.80 "<< fsidebandregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal())<<endl;

        cout<<"+++++++++++++++++++++++++"<<endl;
        cout<<"Ds yield "+run_lable<<" lumi "<<lumi<<endl;
        cout<<"n signal events in 1.93,2.01 "<< nsigevents<<endl;
        cout<<"error on signal events in 1.93,2.01 "<< nsig_err<<endl;
        cout<<"fraction of signal events in 1.93,2.01 "<< fsig <<endl;
        cout<<"n background events in 1.93,2.01 "<< fsigregion_bkg->getVal()*nbkg.getVal() <<endl;
        Double_t ntotalevents = fsigregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal());
        cout<<"n total events in 1.93,2.01 "<< ntotalevents <<endl;

        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << lumi;
        std::string strLumi = stream.str();
        TLatex* cmslabel = new TLatex(0.15,0.81, "#bf{CMS Preliminary}");
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");
        TLatex* text = new TLatex(0.5,0.91, "\n\\text{data }"+run_lable+"\n\\text{ }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
        text->SetNDC(kTRUE);
        text->Draw("same");
        c5->Update();
        c5->SaveAs("plots_dec/DsPhiPi_invmass_"+run_lable+".png");

        TDirectory *dirRun = fout->mkdir(run_lable);
        dirRun->cd();
        dirRun->WriteObject(c5,"2mu1trk_invmass_data");

        //CONTROL PLOTS

        //variable names
        TString var[] = {
                        "Ptmu1","Ptmu2","Ptmu3","Etamu1","Etamu2","Etamu3",
                     //    "SVx", "SVy", "SVz",
                     //    "RefVx1", "RefVy1", "RefVz1",
                     //    "Pmu3","cLP","tKink","segmComp","fv_nC","d0sig","fv_dphi3D","fv_d3D","fv_d3Dsig",
                     //      "bs_sv_d3D","bs_sv_d3Dsig","pv_sv_dxy","pv_sv_dxy_sig",
                     //      "pv_sv_dxy/pv_sv_dxy_sig",
                     //    "mindca_iso","trkRel"
                         };

        int n = sizeof(var)/sizeof(var[0]);
        TH1F *hdata_bkg[n];
        TH1F *hdata_sgn[n];
        TH1F *hmc_sgn[n];
    
        TString binning;

        for(int k = 0; k<n; k++){
            TString varname = var[k];
            cout<<run_lable<<" "<<varname<<endl;    
            TString s = std::to_string(k);
   
            //depending on the variable, you might need to change the binning 
            if(varname=="Ptmu1" || varname=="Ptmu2" || varname=="Ptmu3") binning = "(60,0,30)";
            if(varname=="Etamu1" || varname=="Etamu2" || varname=="Etamu3") binning = "(60,0,3)";
            if(varname=="Pmu3") binning = "(100,0,50)";
            if(varname=="cLP") binning = "(100,0,50)";
            if(varname=="segmComp") binning = "(100,-0.1,1.1)";
            if(varname=="tKink") binning = "(100,0,250)";
            if(varname=="fv_nC") binning = "(104,-0.1,5.1)";
            if(varname=="fv_dphi3D") binning = "(50,0,0.25)";
            if(varname=="fv_d3Dsig") binning = "(100,0,200)";
            if(varname=="d0sig") binning = "(60,0,20)";
            if(varname=="mindca_iso") binning = "(100,0,1)";
            if(varname=="trkRel") binning = "(50,0,1)";
            if(varname=="nMatchesMu3") binning = "(20,0,10)";
            if(varname=="tripletMassReso") binning = "(80,0,0.02)";
            if(varname=="fv_dphi3D") binning = "(100,-0.1,2)";
            if(varname=="fv_d3Dsig") binning = "(80,-0.1,100)";
            if(varname=="fv_d3D") binning = "(40,0,8)";
            if(varname=="fv_d3D/fv_d3Dsig") binning = "(50,0,0.12)";
            if(varname=="bs_sv_d3Dsig") binning = "(80,-0.1,150)";
            if(varname=="bs_sv_d3D") binning = "(50,0,5)";
            if(varname=="bs_sv_d3D/bs_sv_d3Dsig") binning = "(50,0,0.06)";
            if(varname=="abs(dxy1/dxyErr1)" || varname=="abs(dxy2/dxyErr2)" || varname=="abs(dxy3/dxyErr3)" )  binning = "(50,0,40)";
            if(varname=="dxyErr1" || varname=="dxyErr2" || varname=="dxyErr3" )  binning = "(50,0,0.02)";
            if(varname.Contains("mu_pt")) binning = "(80,0,40)";
            if(varname.Contains("mu_eta")) binning = "(80,-2.4,2.4)";
            if(varname.Contains("BestTrackPt")) binning = "(80,0,40)";
            if(varname.Contains("BestTrackEta")) binning = "(80,-2.4,2.4)";
            if(varname.Contains("BestTrackPhi")) binning = "(80,0,3.2)";
            if(varname.Contains("BestTrack") && varname.Contains("Err")) binning = "(80,0,0.0025)";
            if(varname.Contains("BestTrackPtErr")) binning = "(80,0,0.5)";
            if(varname=="abs(RefVz1 - SVz)") binning = "(40,0,8)";
            if(varname=="SVz"|| varname=="RefVz1") binning = "(40,-15,15)";
            if(varname=="SVx" || varname=="SVy") binning = "(40,-1.5,1.5)";
            if(varname=="RefVx1" || varname=="RefVy1") binning = "(120,-0.1,0.1)";
            if(varname=="pv_sv_dxy") binning = "(40,0,3)";
            if(varname=="pv_sv_dxy_sig") binning = "(80,0,100)";
            if(varname=="pv_sv_dxy/pv_sv_dxy_sig") binning = "(80,0,0.04)";

            if(doAllRuns){
                tdata_all->Draw(varname+">>hdata_bkg"+s+binning, invmass_SB);
                tdata_all->Draw(varname+">>hdata_sgn"+s+binning, invmass_peak);
            }else{
                tdata->Draw(varname+">>hdata_bkg"+s+binning, invmass_SB);
                tdata->Draw(varname+">>hdata_sgn"+s+binning, invmass_peak);
            }
            tmc->Draw(varname+">>hmc_sgn"+s+binning, invmass_peak);
    
            hdata_bkg[k] = (TH1F *)gDirectory->Get("hdata_bkg"+s);
            hdata_sgn[k] = (TH1F *)gDirectory->Get("hdata_sgn"+s);
            hmc_sgn[k] = (TH1F *)gDirectory->Get("hmc_sgn"+s);

            TCanvas *c2 = new TCanvas("c2","c2",150,10,990,660);
            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);
            // Upper plot will be in pad1
            TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
            pad1->SetBottomMargin(0); // Upper and lower plot are joined
            pad1->SetGridx();      // Vertical grid
            pad1->Draw();          // Draw the upper pad: pad1
            pad1->cd();            // pad1 becomes the current pad

            //Normalizing Monte Carlo 
            Double_t wNorm = lumi*xsection_mc*BR/N_MC;
            cout<<"wNorm = lumi*xsection_mc*BR/N_MC = "<<wNorm<<endl;
            hmc_sgn[k]->Scale(wNorm);
    
            //scaling the SB distribution to number of background events in 1.93,2.01
            Double_t normSB = hdata_bkg[k]->GetEntries();
            hdata_bkg[k]->Scale(fsigregion_bkg->getVal()*nbkg.getVal()/normSB);
    
            cout<<"Entries in  hdata_sgn[k] before SB subtraction "<<hdata_sgn[k]->GetEntries()<<endl;
            hdata_sgn[k]->Add(hdata_bkg[k],-1); //subtract h2 from h1 : h1->Add(h2,-1)
    
            //rescaling
            hdata_sgn[k]->Scale( hmc_sgn[k]->Integral()/hdata_sgn[k]->Integral() );
 
            cout<<"Entries in  hdata_sgn[k] after SB subtraction "<<hdata_sgn[k]->GetEntries()<<endl;
            cout<<"Entries in  hmc_sgn[k] after rescaling "<<hmc_sgn[k]->GetEntries()<<endl;
            //plot makeup
            double Y_max = std::max(hmc_sgn[k]->GetMaximum(), hdata_sgn[k]->GetMaximum());
            Y_max = Y_max*1.05;
            hmc_sgn[k]->GetYaxis()->SetRangeUser(0, Y_max);

            hmc_sgn[k]->GetYaxis()->SetTitle("a.u.");
            hmc_sgn[k]->GetYaxis()->SetTitleSize(20);
            hmc_sgn[k]->GetYaxis()->SetTitleFont(43);
    
            hmc_sgn[k]->SetLineColor(kBlue);
            hmc_sgn[k]->SetLineWidth(3);
            hmc_sgn[k]->SetFillStyle(3004);
            hmc_sgn[k]->SetFillColor(kBlue);
            hdata_sgn[k]->SetLineColor(kRed);
            hdata_sgn[k]->SetLineWidth(3);
            hdata_sgn[k]->SetFillStyle(3005);
            hdata_sgn[k]->SetFillColor(kRed);

    
            hmc_sgn[k]->Draw("hist");
            hdata_sgn[k]->Draw("hist same");
    
            hmc_sgn[k]->SetStats(0);
    
            TLegend*leg = new TLegend(0.6,0.65,0.9,0.9);
            leg->AddEntry(hmc_sgn[k],"MC DsPhiPi","f");
            leg->AddEntry(hdata_sgn[k],"data "+run_lable+" (SB subtracted)","f");
            leg->Draw();

            // lower plot will be in pad2
            c2->cd();          // Go back to the main canvas before defining pad2
            TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
            pad2->SetTopMargin(0);
            pad2->SetBottomMargin(0.2);
            pad2->SetGridx(); // vertical grid
            pad2->Draw();
            pad2->cd();       // pad2 becomes the current pad    
            // Define the ratio plot
            TH1F *h_x_ratio = (TH1F*)hdata_sgn[k]->Clone("h_x_ratio");
            h_x_ratio->Sumw2(); 
            h_x_ratio->Divide(hmc_sgn[k]);
            h_x_ratio->SetStats(0);
            // Ratio plot settings
            gStyle->SetLineWidth(2);
            h_x_ratio->SetTitle(""); // Remove the ratio title
            h_x_ratio->GetYaxis()->SetTitle("ratio data/MC");
            h_x_ratio->GetYaxis()->SetRangeUser(-0.5,2);
            if(varname=="pv_sv_dxy/pv_sv_dxy_sig") h_x_ratio->GetXaxis()->SetTitle("Error on PV-SV distance on transverse plane (cm)");
            h_x_ratio->SetLineColor(kBlack);
            h_x_ratio->GetYaxis()->SetTitleSize(20);
            h_x_ratio->GetYaxis()->SetTitleFont(43);
            h_x_ratio->GetYaxis()->SetTitleOffset(1.25);
            h_x_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            h_x_ratio->GetYaxis()->SetLabelSize(15);
         
            // X axis ratio plot settings
            h_x_ratio->GetXaxis()->SetTitle(varname);
            h_x_ratio->GetXaxis()->SetTitleSize(20);
            h_x_ratio->GetXaxis()->SetTitleFont(43);
            h_x_ratio->GetXaxis()->SetTitleOffset(3);
            h_x_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            h_x_ratio->GetXaxis()->SetLabelSize(15);

            //Compute weighted average ratio
            Double_t mean = 0;
            Double_t std_dev = 0;
            for (int c=1; c<=(h_x_ratio->GetNbinsX()); c++)
            {
                if(h_x_ratio->GetBinContent(c) == 0 || h_x_ratio->GetBinError(c) == 0) continue;
                mean = mean + h_x_ratio->GetBinContent(c) / ( h_x_ratio->GetBinError(c) * h_x_ratio->GetBinError(c) );
                std_dev = std_dev + 1/( h_x_ratio->GetBinError(c) * h_x_ratio->GetBinError(c));
            }
            mean = mean/std_dev;
            std_dev = 1/std_dev;
            //Get mean value and error of ratio plot
            cout<<var[k]+run_lable+" Mean: "<<mean<<endl;
            cout<<var[k]+run_lable+" StdDev: "<<std_dev<<endl;
            //Draw line corresponding to mean value on ratio plot
            TLine l;
            h_x_ratio->Draw("ep");
            l.DrawLine(0,mean,h_x_ratio->GetXaxis()->GetXmax(),mean);            
            h_x_ratio->Draw("same");

            c2->cd();
            c2->Update();
            c2->SaveAs("plots_dec/"+varname+"_"+run_lable+".png");
            //c2->SaveAs("plots/"+varname+"_"+run_lable+".pdf");
            dirRun->WriteObject(c2,varname+run_lable);
            fout->cd();
        }    
    }    
    fout->Write();
    fout->Close();
}
