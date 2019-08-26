#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1F.h"
#include <cmath>

using namespace RooFit;

void dsphipi_fit() 
{
    TFile *f = new TFile("AnalysedTree_2loose_DsPhiPi_Data2017C_all.root","READ");
    TH1F *h_tripletmass;
    h_tripletmass = (TH1F*)f->Get("StepByStep/Triplet/Mass triplet_cut8");

    // Declare observable x
    TCanvas *c5 = new TCanvas("c5","c5",150,10,990,660);
    c5->Update();
    RooRealVar x("x","2mu+1trk inv. mass (GeV)",1.65,2.2);

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,Import(*h_tripletmass)); 

    // Make plot of binned dataset showing Poisson error bars (RooFit default)
    RooPlot* frame = x.frame(Title(" "));
    dh.plotOn(frame);
    //dh.statOn(frame,Layout(0.55,0.99,0.8)) ;

    //set ranges
    x.setRange("R1",1.84,1.89); //first peak D+(1.87GeV)
    x.setRange("R2",1.94,1.99); //second peak Ds(1.97)
    x.setRange("R3",1.65,1.84); //background    
    x.setRange("R4",1.89,1.92); //background    
    x.setRange("R5",1.99,2.02); //background    
    x.setRange("R6",1.65,2.02); //full range    

    RooRealVar mass2("mass2","Central value of Gaussian",1.958,1.94,1.99);
    RooRealVar sigma2("sigma2","Width of Gaussian",0.01,0,0.1);
    
    RooRealVar mass1("mass1","Central value of Gaussian",1.865,1.84,1.89);
    RooRealVar sigma1("sigma1","Width of Gaussian",0.01,0,0.1);


    // Fit a Crystal ball p.d.f to the data
    RooRealVar alpha2("alpha2","alpha value CB",1,-20,20);
    RooRealVar n2("n2", "n2", 2, 0, 5);
    RooCBShape signal_CB2("cb2", "The signal distribution", x, mass2, sigma2, alpha2, n2); 
    signal_CB2.fitTo(dh, Range("R2"));
    
    
    // Fit a Crystal ball p.d.f to the data
    RooRealVar alpha1("alpha1","alpha value CB",1,-20,20);
    RooRealVar n1("n1", "n1", 2, 0, 5);
    RooCBShape signal_CB1("cb1", "The signal distribution", x, mass1, sigma1, alpha1, n1); 
    signal_CB1.fitTo(dh, Range("R1"));
   

    // Fit an exponential to the background
    RooRealVar a("a", "a", -5, -10, 0.);
    RooExponential bg_exp("bg_exp", "bg_exp", x, a);
    bg_exp.fitTo(dh, Range("R3,R4,R5"));


    // Combine the models
    RooRealVar nsig2("nsig2","#signal2 events",80000,0.,500000);
    RooRealVar nsig1("nsig1","#signal1 events",40000,0.,500000);
    RooRealVar nbkg("nbkg","#background events",400000,0.,1000000);
    RooAddPdf model("model","g+a",RooArgList(signal_CB1,signal_CB2,bg_exp),RooArgList(nsig1,nsig2,nbkg));
    RooFitResult * r = model.fitTo(dh, Save(true));
    r->Print();
    model.paramOn(frame,Layout(0.12, 0.5, 0.95));

    //plot
    model.plotOn(frame);
    model.plotOn(frame, Components(bg_exp), LineColor(kGreen), LineStyle(kDashed));
    model.plotOn(frame, Components(RooArgSet(signal_CB2, signal_CB1)), LineColor(kRed), LineStyle(kDashed) );
    frame->Draw();

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

    cout<<"n background events in 1.70,1.80 "<< fsidebandregion_bkg->getVal()*nbkg.getVal() <<endl;
    cout<<"n total events in 1.70,1.80 "<< fsidebandregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal())<<endl;

    cout<<"n signal events in 1.93,2.01 "<< nsigevents<<endl;
    cout<<"error on signal events in 1.93,2.01 "<< nsig_err<<endl;
    cout<<"fraction of signal events in 1.93,2.01 "<< fsig <<endl;
    cout<<"n background events in 1.93,2.01 "<< fsigregion_bkg->getVal()*nbkg.getVal() <<endl;
    Double_t ntotalevents = fsigregion_model->getVal()*(nsig2.getVal()+nsig1.getVal()+nbkg.getVal());
    cout<<"n total events in 1.93,2.01 "<< ntotalevents <<endl;

}
