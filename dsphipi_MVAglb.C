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
#include "Control_common.h"

using namespace RooFit;

//gErrorIgnoreLevel = kFatal;

void dsphipi_MVAglb() 
{
    bool doPUrew = true;
    
    gErrorIgnoreLevel = kInfo;

    //open root files where to store all plots
    TFile *fout = new TFile("plots_MVAcut/dsphipi_MVAglb_2018_output.root", "RECREATE");
    fout->cd();
    
    TString bdt_cutnos[] = {"0.1","0.47", "0.5385", "0.566", "0.583", "0.611", "0.624", "0.638", "0.6545", "0.8"};

    //set here list of bdt cuts
    TString bdt_cutlist[] = {"MVA_glbmu2>="+bdt_cutnos[0]+" && MVA_glbmu2<"+bdt_cutnos[1],"MVA_glbmu2>="+bdt_cutnos[1]+" && MVA_glbmu2<"+bdt_cutnos[2], "MVA_glbmu2>="+bdt_cutnos[2]+" && MVA_glbmu2<"+bdt_cutnos[3], "MVA_glbmu2>="+bdt_cutnos[3]+" && MVA_glbmu2<"+bdt_cutnos[4], "MVA_glbmu2>="+bdt_cutnos[4]+" && MVA_glbmu2<"+bdt_cutnos[5], "MVA_glbmu2>="+bdt_cutnos[5]+" && MVA_glbmu2<"+bdt_cutnos[6], "MVA_glbmu2>="+bdt_cutnos[6]+" && MVA_glbmu2<"+bdt_cutnos[7], "MVA_glbmu2>="+bdt_cutnos[7]+" && MVA_glbmu2<"+bdt_cutnos[8], "MVA_glbmu2>="+bdt_cutnos[8]+" && MVA_glbmu2<"+bdt_cutnos[9]};
    int ncut = sizeof(bdt_cutlist)/sizeof(bdt_cutlist[0]);
    
    TString eta_cut = " && Etamu2<1.2";
    //TString eta_cut = " && Etamu2>=1.2";

    cout<<"Size of ncut: "<<ncut<<endl;
    
    TString filename_text = "dsphipi_yield_MVAglb_2018.txt";
    //open text file where to write Ds yields
    ofstream fout_yield(filename_text);
    
    Double_t Dsyield_data[ncut];
    Double_t Dsyield_data_err[ncut];
    Double_t Dsyield_MC[ncut];
    Double_t Dsyield_MC_err[ncut];
    Double_t Bkgyield_data[ncut];
    Double_t Bkgyield_data_err[ncut];
  
    TFile *fin = new TFile("MVAmu_control_2018_outputTree.root", "READ");
    cout<<"opened input file MVAmu_control_2018_outputTree.root"<<endl;
    TTree *tin = (TTree*)fin->Get("outputTree");
    
    TH1F *h_tripletmass[ncut];//all data range with mva cut
    TH1F *h_tripletmass_mc[ncut];//peak MC range with mva cut
    TH1F *h_tripletmass_full;
    TH1F *h_tripletmass_sign[ncut];//peak data range with mva cut
	
    Double_t events_SB = 0;
    TString binning_mass = "(62, 1.72, 2.01)";
    Int_t n_bins = 62;
	
    Double_t lumi_full = 3.0; //fb
    TString run_lable = "2018";

    TString common_cut = " && bs_sv_d2Dsig>2.0 && Ptmu3 > 1.2 && " 
                         "((Ptmu1>3.5 && Etamu1<1.2) || (Ptmu1>2.0 && Etamu1>=1.2 && Etamu1<=2.4)) && "
                         "((Ptmu2>3.5 && Etamu2<1.2) || (Ptmu2>2.0 && Etamu2>=1.2 && Etamu2<=2.4)) && "
                         "abs(phiMass-1.02)<0.045 && "
                         "!(l1double_DoubleMu4_fired && !l1double_DoubleMu0_fired)";

    TString invmass_SB   = "";
    TString invmass_peak = "";
    TString invmass_all_data = "";
    TString invmass_peak_MC = "";
    TString invmass_peak_data = "";
    
	
    if(doPUrew){
        invmass_all_data  = "puFactor*(tripletMass<2.02 && tripletMass>1.62"+common_cut+eta_cut+" && isMC==0";
        invmass_SB   = "puFactor*(tripletMass<1.80 && tripletMass>1.73"+common_cut+eta_cut+" && isMC==0";
        invmass_peak = "puFactor*(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==1";
        invmass_peak_MC = "puFactor*(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==1";
        invmass_peak_data = "puFactor*(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==0";
    } else {
        invmass_all_data  = "(tripletMass<2.02 && tripletMass>1.62"+common_cut+eta_cut+" && isMC==0";
        invmass_SB   = "(tripletMass<1.80 && tripletMass>1.73"+common_cut+eta_cut+" && isMC==0";
        invmass_peak = "(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==1";
        invmass_peak_MC = "(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==1";
        invmass_peak_data = "(tripletMass<2.01 && tripletMass>1.93"+common_cut+eta_cut+" && isMC==0";
    }

    //Fill histograms for each cut
    for(int i=0; i<ncut; i++){
        TString s = std::to_string(i-1);
        
        TString bdt_cut = bdt_cutlist[i];
        TString bdt_cut_label = bdt_cutlist[i];
        bdt_cut_label = bdt_cut_label.ReplaceAll(".", "p");
        bdt_cut_label = bdt_cut_label.ReplaceAll(">", "_");
        bdt_cut_label = bdt_cut_label.ReplaceAll("<", "_");
    		
        tin->Draw("tripletMass>>h_tripletmass["+s+"]"+binning_mass, invmass_all_data+bdt_cut);
        h_tripletmass[i]     = (TH1F *)gDirectory->Get("h_tripletmass["+s+"]");
        
        tin->Draw("tripletMass>>h_tripletmass_sign["+s+"]"+binning_mass, invmass_peak_data+bdt_cut);
        h_tripletmass_sign[i]     = (TH1F *)gDirectory->Get("h_tripletmass_sign["+s+"]");
        
        tin->Draw("tripletMass>>h_tripletmass_mc["+s+"]"+binning_mass, invmass_peak_MC+bdt_cut);
        h_tripletmass_mc[i]     = (TH1F *)gDirectory->Get("h_tripletmass_mc["+s+"]");
        
        cout<<"Events passing selections "<<bdt_cutlist[i]<<" = "<<h_tripletmass[i]->GetEntries()<<endl;
    }
	
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
        nbkg[i] = new RooRealVar("nbkg_"+category,"#bkacground events",int(entries/2),5,entries);
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
        //fraction of total events in 1.70,1.80 (n_sideband_region_events/n_total_events)
        RooAbsReal* fsidebandregion_model = model[i]->createIntegral(x,NormSet(x),Range("sideband")); 

        //fraction of background events in 1.93,2.01
        RooAbsReal* fsigregion_bkg = ((RooAbsPdf*)model[i]->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fb = fsigregion_bkg->getVal();
        Double_t fb_err = fsigregion_bkg->getPropagatedError(*r);
        //fraction of background events in 1.70, 1.80 
        RooAbsReal* fsidebandregion_bkg = ((RooAbsPdf*)model[i]->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("sideband")); 

        Dsyield_data[i] = fs * (nsig2[i]->getVal()+nsig1[i]->getVal()+nbkg[i]->getVal()) - fb*nbkg[i]->getVal();;
        Dsyield_data_err[i] = pow( pow(fs_err,2) * pow(nsig2[i]->getVal()+nsig1[i]->getVal()+nbkg[i]->getVal(),2)  + ( pow(nsig2[i]->getPropagatedError(*r),2)+pow(nsig1[i]->getPropagatedError(*r),2)+pow(nbkg[i]->getPropagatedError(*r),2)) * pow(fs,2) + pow(fb_err,2) * pow(nbkg[i]->getVal(),2) + pow(nbkg[i]->getPropagatedError(*r),2)*pow(fb,2) , 0.5);;

        Bkgyield_data[i] = fb*nbkg[i]->getVal();
        Bkgyield_data_err[i] = sqrt( pow(fb_err, 2.0) * pow(nbkg[i]->getVal(), 2.0) + pow(fb, 2.0) * pow(nbkg[i]->getPropagatedError(*r), 2.0));
    }

    //drawing triplet mass in MC for region mass selection and integral
    TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
    h_tripletmass_mc[0]->Draw();
    auto f1  = new TF1("f1","gaus",1.93,2.01);
    h_tripletmass_mc[0]->Fit("f1", "R");
    f1->Draw("same");
    c3->Update();
    fout->WriteObject(c3,"3glb_invmass_mc");
    Double_t n_mc_peak = f1->Integral(1.93, 2.01) / h_tripletmass_mc[0]->Integral(h_tripletmass_mc[0]->FindFixBin(1.93),h_tripletmass_mc[0]->FindFixBin(2.01),"width") * h_tripletmass_mc[0]->Integral(h_tripletmass_mc[0]->FindFixBin(1.93),h_tripletmass_mc[0]->FindFixBin(2.01));

    for(int i = 0; i<ncut; i++){
        //cout<<"n_mc_peak "<<n_mc_peak<<endl;
        //cout<<"scaled to lumi: "<<n_mc_peak*lumi[i]*xsection_mc*BR/N_MC<<endl;
        Dsyield_MC[i] = 0.90*h_tripletmass_mc[i]->GetEntries()*lumi_full*xsection_mc*BR/N_MC;
        //Dsyield_MC[i] = h_tripletmass_mc[i]->GetEntries();
        Dsyield_MC_err[i] = sqrt(0.90*h_tripletmass_mc[i]->GetEntries())*(lumi_full*xsection_mc*BR/N_MC); 
        //Dsyield_MC_err[i] = (f1->IntegralError(1.93, 2.01))*(lumi[i]*xsection_mc*BR/N_MC); 
        //Dsyield_MC[i] = h_tripletmass_mc[0]->GetEntries();

        fout_yield<<Dsyield_data[i]<<"\t"<<Dsyield_data_err[i]<<"\n";

        cout<<"\n"<<bdt_cutlist[i]<<" lumi="<<lumi_full<<endl;
        cout<<"=================\noverall data/MC scale factor:"<<endl;
        cout<<"data: "<<Dsyield_data[i]<<" +- "<<Dsyield_data_err[i]<<"\n";
        cout<<"MC: "<<Dsyield_MC[i]<<" +- "<<Dsyield_MC_err[i]<<endl;
        cout<<"scale factor: "<<Dsyield_data[i]/Dsyield_MC[i]<<" +- "<<sqrt(pow((Dsyield_data_err[i]/Dsyield_MC[i]),2.0) + pow((Dsyield_data[i]/(Dsyield_MC[i]*Dsyield_MC[i]))*Dsyield_MC_err[i],2.0))<<endl;
    }
return;
 
}
