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
    bool doPUrew = false;
    
    gErrorIgnoreLevel = kInfo;

    //open root files where to store all plots!!
    TFile *fout = new TFile("plots_BDTcut/dsphipi_all_output_"+TMVA_inputpath_control+"_BDT_cut.root", "RECREATE");
    fout->cd();
    
    TString bdt_cutnos[] = {"0.1","0.47", "0.5385", "0.566", "0.583", "0.611", "0.624", "0.638", "0.6545", "0.8"};

    //set here list of bdt cuts
    //TString bdt_cutlist[] = {"", "bdt>-0.02","bdt>0", "bdt>0.02", "bdt>0.03", "bdt>0.04", "bdt>0.06", "bdt>0.07", "bdt>0.075", "bdt>0.08", "bdt>0.10"};
    TString bdt_cutlist[] = {"", "MVA_glbmu2>="+bdt_cutnos[0]+" && MVA_glbmu2<"+bdt_cutnos[1],"MVA_glbmu2>="+bdt_cutnos[1]+" && MVA_glbmu2<"+bdt_cutnos[2], "MVA_glbmu2>="+bdt_cutnos[2]+" && MVA_glbmu2<"+bdt_cutnos[3], "MVA_glbmu2>="+bdt_cutnos[3]+" && MVA_glbmu2<"+bdt_cutnos[4], "MVA_glbmu2>="+bdt_cutnos[4]+" && MVA_glbmu2<"+bdt_cutnos[5], "MVA_glbmu2>="+bdt_cutnos[5]+" && MVA_glbmu2<"+bdt_cutnos[6], "MVA_glbmu2>="+bdt_cutnos[6]+" && MVA_glbmu2<"+bdt_cutnos[7], "MVA_glbmu2>="+bdt_cutnos[7]+" && MVA_glbmu2<"+bdt_cutnos[8], "MVA_glbmu2>="+bdt_cutnos[8]+" && MVA_glbmu2<"+bdt_cutnos[9]};
    int ncut = sizeof(bdt_cutlist)/sizeof(bdt_cutlist[0]);
    
    cout<<"Size of ncut: "<<ncut<<endl;
    
    TString filename_text = "dsphipi_yield_"+TMVA_inputpath_control+"perMVA.txt";
    //open text file where to write Ds yields
    ofstream fout_yield(filename_text);
    
    
    //New code start
    
    Double_t Dsyield_data[ncut-1];
    Double_t Dsyield_data_err[ncut-1];
    Double_t Dsyield_MC[ncut-1];
    Double_t Dsyield_MC_err[ncut-1];
    Double_t Bkgyield_data[ncut-1];
    Double_t Bkgyield_data_err[ncut-1];
    
    Double_t Dsyield_data_full;
    Double_t Dsyield_data_err_full;
    Double_t Dsyield_MC_full;
    Double_t Dsyield_MC_err_full;
    Double_t Bkgyield_data_full;
    Double_t Bkgyield_data_err_full;
	
	  //TFile *fin = new TFile("../MVA_control_2018_30july_dphi3DoutputTree.root", "READ");
	  TFile *fin = new TFile("MVAmu_control_2018_outputTree.root", "READ");
    cout<<"opened input file MVAmu_control_2018_outputTree.root"<<endl;
    TTree *tin = (TTree*)fin->Get("outputTree");
    
    TH1F *h_tripletmass[ncut-1];//all data range with mva cut
	  TH1F *h_tripletmass_mc[ncut-1];//peak MC range with mva cut
    TH1F *h_tripletmass_mc_full;//peak MC range no mva cut
    TH1F *h_tripletmass_full;
	  TH1F *h_tripletmass_sign[ncut-1];//peak data range with mva cut
	
	  int nrun = sizeof(inputpath_datarun_control)/sizeof(inputpath_datarun_control[0]);
	  TChain *tdata_all = new TChain("FinalTree_Control");
	  TChain *tdata_all_mu1 = new TChain("TreeMu1");
    TChain *tdata_all_mu2 = new TChain("TreeMu2");
	
	  for(auto i=0; i<nrun; i++){
        //TFile *f = new TFile(inputpath_datarun_control[i],"READ");
        //tdata_all->Add(inputpath_datarun_control[i]);
		  //tdata_all_mu1->Add(inputpath_datarun_control[i]); 
        //tdata_all_mu2->Add(inputpath_datarun_control[i]); 
        //f->Close();
    }
	
    Double_t events_SB = 0;
    TString binning_mass = "(62, 1.72, 2.01)";
    Int_t n_bins = 62;
	
	  Double_t lumi[nrun];
    Double_t lumi_full = 0.0;
    TString run_lable[ncut-1];
    TString run_lable_full;

    for(int k=0; k<ncut-1; k++){
        
        TString bdt_cut_label = bdt_cutlist[k+1];
        
        bdt_cut_label = bdt_cut_label.ReplaceAll(".", "p");
    		bdt_cut_label = bdt_cut_label.ReplaceAll(">", "_");
    		bdt_cut_label = bdt_cut_label.ReplaceAll("<", "_");
    }
    for(int k=0; k<nrun; k++){
        lumi[k] = Lumi_data_control[k];
        //run_lable[k] = run_name_control[k];
        //run_lable[k] = bdt_cut_label;
        lumi_full += lumi[k];
        run_lable_full = "2018";
    }
    
    //lumi_full = 59.7;
	
	  TString common_cut = " && bs_sv_d2Dsig>2.0 && Ptmu3 > 1.2 && " 
                         "((Ptmu1>3.5 && Etamu1<1.2) || (Ptmu1>2.0 && Etamu1>=1.2 && Etamu1<=2.4)) && "
                         "((Ptmu2>3.5 && Etamu2<1.2) || (Ptmu2>2.0 && Etamu2>=1.2 && Etamu2<=2.4)) && "
                         "abs(phiMass-1.02)<0.045 && "
                         "!(l1double_DoubleMu4_fired && !l1double_DoubleMu0_fired)";
                         //"(l1double_DoubleMu4_fired || l1double_DoubleMu0_fired)";

    TString invmass_SB   = "";
    TString invmass_peak = "";
    TString invmass_all_data = "";
	  TString invmass_peak_MC = "";
    TString invmass_peak_data = "";
    
    TString eta_cut = " && Etamu2<1.2";
    //TString eta_cut = " && Etamu2>=1.2";
	
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

    for(int i=0; i<ncut; i++){
        TString s = std::to_string(i-1);
        
        TString bdt_cut = bdt_cutlist[i];
        TString bdt_cut_label = bdt_cutlist[i];
    		
        if(bdt_cut==""){
            
    			if(bdt_cut!="") bdt_cut = +"&& "+bdt_cut+")";
    			else            bdt_cut = ")";
    
    			bdt_cut_label = bdt_cut_label.ReplaceAll(".", "p");
    			bdt_cut_label = bdt_cut_label.ReplaceAll(">", "_");
    			bdt_cut_label = bdt_cut_label.ReplaceAll("<", "_");
    			
    			//cout<<"tripletMass>>h_tripletmass"+binning_mass+", "+invmass_all_data+bdt_cut<<endl;
    			
    			//Fill histogram for full 2018
    			tin->Draw("tripletMass>>h_tripletmass_full"+binning_mass, invmass_all_data+bdt_cut);
    			h_tripletmass_full = (TH1F *)gDirectory->Get("h_tripletmass_full");
    			
    			//Fill histogram for MC
    			tin->Draw("tripletMass>>h_tripletmass_mc_full"+binning_mass, invmass_peak_MC+bdt_cut);
    			h_tripletmass_mc_full = (TH1F *)gDirectory->Get("h_tripletmass_mc_full");
    
    			cout<<"Events passing selections "<<bdt_cutlist[i]<<" = "<<h_tripletmass_full->GetEntries()<<endl;
    		
    		}
    		
    		else{
    		
    			if(bdt_cut!="") bdt_cut = +"&& "+bdt_cut+")";
    			else            bdt_cut = ")";
    
    			bdt_cut_label = bdt_cut_label.ReplaceAll(".", "p");
    			bdt_cut_label = bdt_cut_label.ReplaceAll(">", "_");
    			bdt_cut_label = bdt_cut_label.ReplaceAll("<", "_");
    			
    			//cout<<"tripletMass>>h_tripletmass"+binning_mass+", "+invmass_all_data+bdt_cut<<endl;
    
    			tin->Draw("tripletMass>>h_tripletmass["+s+"]"+binning_mass, invmass_all_data+bdt_cut);
    			h_tripletmass[i-1]     = (TH1F *)gDirectory->Get("h_tripletmass["+s+"]");
    			
    			tin->Draw("tripletMass>>h_tripletmass_sign["+s+"]"+binning_mass, invmass_peak_data+bdt_cut);
    			h_tripletmass_sign[i-1]     = (TH1F *)gDirectory->Get("h_tripletmass_sign["+s+"]");
    			
    			tin->Draw("tripletMass>>h_tripletmass_mc["+s+"]"+binning_mass, invmass_peak_MC+bdt_cut);
    			h_tripletmass_mc[i-1]     = (TH1F *)gDirectory->Get("h_tripletmass_mc["+s+"]");
    			
    			cout<<"Events passing selections "<<bdt_cutlist[i]<<" = "<<h_tripletmass[i-1]->GetEntries()<<endl;
    		
    		}
		
    }
	
    std::map<std::string, TH1 *> hmap;
    RooCategory c("c", "c");
    for(int i = 0; i<ncut-1; i++){
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
    RooRealVar *nsig2[ncut-1];
    RooRealVar *nsig1[ncut-1];
    RooRealVar *nbkg[ncut-1];
    RooAddPdf *model[ncut-1];

    // Associate model with the categories. In our case, model is the same for all categories
    for(int i = 0; i<ncut-1; i++){
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
    for(int i = 0; i<ncut-1; i++){
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
        TLatex* text_lumi = new TLatex(0.10,0.91, "\n\\text{data }"+run_lable[i]+"\n\\text{    }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
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
    h_tripletmass_mc_full->Draw();
    auto f1  = new TF1("f1","gaus",1.93,2.01);
    h_tripletmass_mc_full->Fit("f1", "R");
    f1->Draw("same");
    c3->Update();
    fout->WriteObject(c3,"3glb_invmass_mc");
    Double_t n_mc_peak = f1->Integral(1.93, 2.01) / h_tripletmass_mc_full->Integral(h_tripletmass_mc_full->FindFixBin(1.93),h_tripletmass_mc_full->FindFixBin(2.01),"width") * h_tripletmass_mc_full->Integral(h_tripletmass_mc_full->FindFixBin(1.93),h_tripletmass_mc_full->FindFixBin(2.01));

    for(int i = 0; i<ncut-1; i++){
        //cout<<"n_mc_peak "<<n_mc_peak<<endl;
        //cout<<"scaled to lumi: "<<n_mc_peak*lumi[i]*xsection_mc*BR/N_MC<<endl;
        Dsyield_MC[i] = 0.90*h_tripletmass_mc[i]->GetEntries()*lumi_full*xsection_mc*BR/N_MC;
        //Dsyield_MC[i] = h_tripletmass_mc[i]->GetEntries();
        Dsyield_MC_err[i] = sqrt(0.90*h_tripletmass_mc[i]->GetEntries())*(lumi_full*xsection_mc*BR/N_MC); 
        //Dsyield_MC_err[i] = (f1->IntegralError(1.93, 2.01))*(lumi[i]*xsection_mc*BR/N_MC); 
        //Dsyield_MC[i] = h_tripletmass_mc_full->GetEntries();

        fout_yield<<Dsyield_data[i]<<"\t"<<Dsyield_data_err[i]<<"\n";

        cout<<"\n"<<bdt_cutlist[i]<<" lumi="<<lumi_full<<endl;
        cout<<"=================\noverall data/MC scale factor:"<<endl;
        cout<<"data: "<<Dsyield_data[i]<<" +- "<<Dsyield_data_err[i]<<"\n";
        cout<<"MC: "<<Dsyield_MC[i]<<" +- "<<Dsyield_MC_err[i]<<endl;
        cout<<"scale factor: "<<Dsyield_data[i]/Dsyield_MC[i]<<" +- "<<sqrt(pow((Dsyield_data_err[i]/Dsyield_MC[i]),2.0) + pow((Dsyield_data[i]/(Dsyield_MC[i]*Dsyield_MC[i]))*Dsyield_MC_err[i],2.0))<<endl;
    }
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
    model_full->plotOn(frame_full, Components(RooArgSet(signal_CB2, signal_CB1)), LineColor(kRed), LineStyle(kDashed) );
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
    TLatex* text_lumi = new TLatex(0.10,0.91, "\n\\text{data }"+run_lable_full+"\n\\text{    }\n\\mathscr{L}="+strLumi+"\\text{fb}^{-1}");
    text_lumi->SetTextSize(0.05);
    text_lumi->SetNDC(kTRUE);
    text_lumi->Draw("same");

    Double_t Chi2 = chi2.getVal()/NDOF;
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

    fout->WriteObject(c5,"full");
    c5->SaveAs("dsphipi_fit_perMVA_new_full.png");

    //Integrals
        //fraction of total events in 1.93,2.01 (n_signal_region_events/n_total_events)
        RooAbsReal* fsigregion_model = model_full->createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fs = fsigregion_model->getVal();
        Double_t fs_err = fsigregion_model->getPropagatedError(*r_full);
        //fraction of total events in 1.70,1.80 (n_sideband_region_events/n_total_events)
        RooAbsReal* fsidebandregion_model = model_full->createIntegral(x,NormSet(x),Range("sideband")); 

        //fraction of background events in 1.93,2.01
        RooAbsReal* fsigregion_bkg = ((RooAbsPdf*)model_full->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("signal")); 
        Double_t fb = fsigregion_bkg->getVal();
        Double_t fb_err = fsigregion_bkg->getPropagatedError(*r);
        //fraction of background events in 1.70, 1.80 
        RooAbsReal* fsidebandregion_bkg = ((RooAbsPdf*)model_full->getComponents()->find("bg_exp"))->createIntegral(x,NormSet(x),Range("sideband")); 

        Dsyield_data_full = fs * (nsig2_full->getVal()+nsig1_full->getVal()+nbkg_full->getVal()) - fb*nbkg_full->getVal();;
        Dsyield_data_err_full = pow( pow(fs_err,2) * pow(nsig2_full->getVal()+nsig1_full->getVal()+nbkg_full->getVal(),2)  + ( pow(nsig2_full->getPropagatedError(*r),2)+pow(nsig1_full->getPropagatedError(*r),2)+pow(nbkg_full->getPropagatedError(*r),2)) * pow(fs,2) + pow(fb_err,2) * pow(nbkg_full->getVal(),2) + pow(nbkg_full->getPropagatedError(*r),2)*pow(fb,2) , 0.5);;

        Bkgyield_data_full = fb*nbkg_full->getVal();
        Bkgyield_data_err_full = sqrt( pow(fb_err, 2.0) * pow(nbkg_full->getVal(), 2.0) + pow(fb, 2.0) * pow(nbkg_full->getPropagatedError(*r), 2.0));
        //cout<<"n_mc_peak "<<n_mc_peak<<endl;
        //cout<<"scaled to lumi: "<<n_mc_peak*lumi_full*xsection_mc*BR/N_MC<<endl;
        Dsyield_MC_full = n_mc_peak*lumi_full*xsection_mc*BR/N_MC;
        Dsyield_MC_err_full = sqrt(n_mc_peak)*(lumi_full*xsection_mc*BR/N_MC); 
        //Dsyield_MC_err[i] = (f1->IntegralError(1.93, 2.01))*(lumi[i]*xsection_mc*BR/N_MC); 
        //Dsyield_MC[i] = h_tripletmass_mc_full->GetEntries();

        fout_yield<<Dsyield_data_full<<"\t"<<Dsyield_data_err_full<<"\n";

        cout<<"\n"<<run_lable_full<<" lumi="<<lumi_full<<endl;
        cout<<"=================\noverall data/MC scale factor:"<<endl;
        cout<<"data: "<<Dsyield_data_full<<" +- "<<Dsyield_data_err_full<<"\n";
        cout<<"MC: "<<Dsyield_MC_full<<" +- "<<Dsyield_MC_err_full<<endl;
        cout<<"scale factor: "<<Dsyield_data_full/Dsyield_MC_full<<" +- "<<sqrt(pow((Dsyield_data_err_full/Dsyield_MC_full),2.0) + pow((Dsyield_data_full/(Dsyield_MC_full*Dsyield_MC_full))*Dsyield_MC_err_full,2.0))<<endl;
    
return;
    //Drawing control plots using yields from full 2018
    //List of variables
        //variables plot
        TString var[] = {
//                         "MuonIDeval_Mu1.MuonID",
//                         "MuonIDeval_Mu2.MuonID",
                         "Ptmu1",
                         "Ptmu2","Ptmu3","abs(Etamu1)","abs(Etamu2)","abs(Etamu3)",
                //         "Pt_tripl","abs(Eta_tripl)","cLP>30?30:cLP",
                         "tKink>80?80:tKink",
                         "segmComp",
                         "fv_nC>25?25:fv_nC",
                         //"fv_dphi3D>0.15?0.15:fv_dphi3D","fv_d3Dsig>100?100:fv_d3Dsig",
                         "fv_dphi3D","fv_d3Dsig",
                         "fv_d3D",
                         "bs_sv_d2Dsig",
                       //  "pv_sv_dxy_sig",
                       //  "pv_sv_dxy",
                       //  "pv_sv_dxy_err",
                       //  "pv_sv_dxy/pv_sv_dxy_sig",
                 //        "d0sig",
                         "abs(dxy1)/dxyErr1",
                         "abs(dxy2)/dxyErr2",
                         "mindca_iso>0.5?0.5:mindca_iso", "trkRel>10?10:trkRel",
                         "TreeMu1.mu_sumPt03",
                         "TreeMu2.mu_sumPt03",
                         "TreeMu1.mu_sumPt03/TreeMu1.mu_pt",
                         "TreeMu2.mu_sumPt03/TreeMu2.mu_pt",
                         "TreeMu1.mu_nTracks03",
                         "TreeMu2.mu_nTracks03"
                        };
        TString var_names[] = {
  //                   "MVAmuID_glb_mu1",
  //                   "MVAmuID_glb_mu2",
                     "p_{T} (#mu_{1}) (GeV)",
                     "p_{T} (#mu_{2}) (GeV)",
                     "p_{T} (trk) (GeV)",
                     "|#eta| (#mu_{1})",
                     "|#eta| (#mu_{2})",
                     "|#eta| (trk)",
                //     "p_{T} (#mu#mu#pi)",
                //     "|#eta| (#mu#mu#pi)",
                //     "#chi^{2} of muon inner-outer tracks position matching (cLP_max)",
                     "muon track Kink (tKink_max)",
                     "muon segment compatibility (segmComp_min)",
                     "SV fit #chi^{2} (fv_nC)",
                     "pointing angle [SV-PV, p(#mu#mu#pi)] (fv_dphi3D)","significance of 3-dim SV-PV displacement (fv_d3Dsig)",
                     "SV-PV 3-dim displacement (fv_d3D)",
                     "significance of SV transverse displacement wrt the beamspot",
                //     "significance of SV-PV transverse displacement",
                //     "SV-PV transverse displacement",
                //     "uncertainty on the SV-PV transverse displacement",
                //     "uncertainty on the SV-PV transverse displacement",
                //     "d0sig",
                     "significance of #mu_{1} transverse displacement (d0sig #mu_{1})",
                     "significance of #mu_{2} transverse displacement (d0sig #mu_{1})",
                     "mindca_iso", "trkRel",
                     "summed p_{T} of the tracks in the #DeltaR<0.3 isolation cone (#mu_{1})",
                     "summed p_{T} of the tracks in the #DeltaR<0.3 isolation cone (#mu_{2})",
                     "summed p_{T} of the tracks in the #DeltaR<0.3 isolation cone / muon p_{T} (#mu_{1})",
                     "summed p_{T} of the tracks in the #DeltaR<0.3 isolation cone / muon p_{T} (#mu_{2})",
                     "number of charged tracks in the #DeltaR<0.3 isolation cone (#mu_{1})",
                     "number of charged tracks in the #DeltaR<0.3 isolation cone (#mu_{2})"
        };


    int n = sizeof(var)/sizeof(var[0]);
    TH1F *hdata_bkg[n];
    TH1F *hdata_bkg_plus[n];
    TH1F *hdata_bkg_minus[n];
    TH1F *hdata_sgn[n];
    TH1F *hdata_sgn_plus[n];
    TH1F *hdata_sgn_minus[n];
    TH1F *hmc_sgn[n];
    
    TString binning = "";

    for(int k = 0; k<n; k++){

        TString varname = var[k];
        TString varlable = var_names[k];
        cout<<run_lable[ncut-1-1]<<" "<<varname<<endl;
        TString s = std::to_string(k);

            if(varname=="Pt_tripl") binning = "(50,0,80)";
            if(varname=="Ptmu1" || varname=="Ptmu2" || varname=="Ptmu3") binning = "(60,0,30)";
            if(varname=="abs(Etamu1)" || varname=="abs(Etamu2)" || varname=="abs(Etamu3)" || varname=="abs(Eta_tripl)") binning = "(50,0.0,2.5)";
            if(varname.Contains("Pmu3")) binning = "(100,0,50)";
            if(varname.Contains("cLP")) binning = "(60,0,20)";
            if(varname.Contains("segmComp")) binning = "(100,-0.1,1.1)";
            if(varname.Contains("tKink")) binning = "(50,0,50)";
            if(varname.Contains("fv_nC")) binning = "(25,-0.1,5.1)";
            if(varname.Contains("d0sig")) binning = "(36,-0.1,18)";
            if(varname.Contains("mindca_iso")) binning = "(25,0,0.5)";
            if(varname.Contains("trkRel")) binning = "(30,0.05,10)";
            if(varname.Contains("nMatchesMu3")) binning = "(20,0,20)";
            if(varname.Contains("tripletMassReso")) binning = "(80,0,0.02)";
            if(varname.Contains("fv_dphi3D")) binning = "(80,-0.01,0.15)";
            if(varname.Contains("fv_d3Dsig")) binning = "(80,0.0,100)";
            if(varname=="fv_d3D")             binning = "(50,0,3.0)";
            if(varname.Contains("fv_d3D/fv_d3Dsig")) binning = "(50,0,0.12)";
            if(varname.Contains("bs_sv_d2D")) binning = "(80,-0.1,1.5)";
            if(varname.Contains("bs_sv_d2Dsig")) binning = "(80,-0.1,150)";
            if(varname.Contains("bs_sv_d2D/bs_sv_d2Dsig")) binning = "(50,0,0.06)";
            if(varname.Contains("abs(dxy1/dxyErr1)") || varname=="abs(dxy2/dxyErr2)" || varname=="abs(dxy3/dxyErr3)" )  binning = "(60,0,40)";
            if(varname.Contains("dxy1") || varname=="dxy2" || varname=="dxy3" )  binning = "(60,-0.1,2.5)";
            if(varname.Contains("dxyErr1") || varname=="dxyErr2" || varname=="dxyErr3" )  binning = "(50,0,0.02)";
            if(varname.Contains("dxyErr")) binning = "(61,-0.5,30)";
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
            if(varname=="pv_sv_dxy") binning = "(80,-0.1,2)";
            if(varname=="pv_sv_dxy_sig") binning = "(40,0,100)";
            if(varname=="pv_sv_dxy/pv_sv_dxy_sig" || varname=="pv_sv_dxy_err") binning = "(80,0,0.04)";
            if(varname.Contains("MuonIDeval")) binning = "(50,0.1,0.8)";
            if(varname.Contains("MuonID")) binning = "(50,0.1,0.8)";
            if(varname=="fv_nC>25?25:fv_nC") binning = "(50,-0.2,25.0)";
            if(varname.Contains("mu_sumPt03")) binning = "(80,-0.2,20.0)";
            if(varname.Contains("mu_nTracks03")) binning = "(20,0,20.0)";

        cout<<"binning "<<binning<<endl;
        tdata_all->Draw(varname+">>hdata_bkg"+s+binning, invmass_SB);
        tdata_all->Draw(varname+">>hdata_bkg_plus"+s+binning, invmass_SB);
        tdata_all->Draw(varname+">>hdata_bkg_minus"+s+binning, invmass_SB);
        tdata_all->Draw(varname+">>hdata_sgn"+s+binning, invmass_peak);
        tdata_all->Draw(varname+">>hdata_sgn_plus"+s+binning, invmass_peak);
        tdata_all->Draw(varname+">>hdata_sgn_minus"+s+binning, invmass_peak);

        tin->Draw(varname+">>hmc_sgn"+s+binning, invmass_peak);

        hdata_bkg[k] = (TH1F *)gDirectory->Get("hdata_bkg"+s);
        hdata_bkg_plus[k] = (TH1F *)gDirectory->Get("hdata_bkg_plus"+s);
        hdata_bkg_minus[k] = (TH1F *)gDirectory->Get("hdata_bkg_minus"+s);
        hdata_sgn[k] = (TH1F *)gDirectory->Get("hdata_sgn"+s);
        hdata_sgn_plus[k] = (TH1F *)gDirectory->Get("hdata_sgn_plus"+s);
        hdata_sgn_minus[k] = (TH1F *)gDirectory->Get("hdata_sgn_minus"+s);
        hmc_sgn[k] = (TH1F *)gDirectory->Get("hmc_sgn"+s);

        TCanvas *c2 = new TCanvas("c2","c2",150,10,800,800);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        // Upper plot will be in pad1
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->SetGridx();      // Vertical grid
        pad1->Draw();          // Draw the upper pad: pad1
        pad1->cd();            // pad1 becomes the current pad
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);

        Double_t normMC = hmc_sgn[k]->GetEntries();
        //Normalizing Monte Carlo 
        Double_t wNorm = lumi_full*xsection_mc*BR/N_MC;
        //Double_t wNorm = lumi_full*xsection_mc*BR/N_MC  *  n_mc_peak/hmc_sgn[k]->GetEntries();
        cout<<"wNorm = lumi_full*xsection_mc*BR/N_MC = "<<wNorm<<endl;
        hmc_sgn[k]->Scale(wNorm);

        //scaling the SB distribution to number of background events in 1.93,2.01
        Double_t normSB = hdata_bkg[k]->GetEntries();
        hdata_bkg[k]->Scale(Bkgyield_data_full/normSB);
        hdata_bkg_plus[k]->Scale( (Bkgyield_data_full/normSB) * 1.10);
        hdata_bkg_minus[k]->Scale( (Bkgyield_data_full/normSB) * 0.90);

        cout<<"Entries in  hdata_sgn[k] before SB subtraction "<<hdata_sgn[k]->GetEntries()<<endl;
        hdata_sgn[k]->Add(hdata_bkg[k],-1); //subtract h2 from h1 : h1->Add(h2,-1)
        hdata_sgn_plus[k]->Add(hdata_bkg_plus[k],-1); //subtract h2 from h1 : h1->Add(h2,-1)
        hdata_sgn_minus[k]->Add(hdata_bkg_minus[k],-1); //subtract h2 from h1 : h1->Add(h2,-1)

        //Rescaling to same integral
        hmc_sgn[k]->Scale( 1.0 / hmc_sgn[k]->Integral());
        hdata_sgn[k]->Scale( hmc_sgn[k]->Integral()/hdata_sgn[k]->Integral() );
        hdata_sgn_plus[k]->Scale( hmc_sgn[k]->Integral()/hdata_sgn_plus[k]->Integral() );
        hdata_sgn_minus[k]->Scale( hmc_sgn[k]->Integral()/hdata_sgn_minus[k]->Integral() );

        //plot makeup
        double Y_max = std::max(hmc_sgn[k]->GetMaximum(), hdata_sgn[k]->GetMaximum());
        Y_max = Y_max*1.4;
        hmc_sgn[k]->GetYaxis()->SetRangeUser(0, Y_max);

        hmc_sgn[k]->GetYaxis()->SetTitle("a.u.");
        hmc_sgn[k]->GetYaxis()->SetTitleSize(22);
        hmc_sgn[k]->GetYaxis()->SetTitleFont(43);
        hmc_sgn[k]->GetYaxis()->SetTitleOffset(1.7);

        hmc_sgn[k]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        hmc_sgn[k]->GetXaxis()->SetTitleOffset(1.5);

        hmc_sgn[k]->SetLineColor(kBlue);
        hmc_sgn[k]->SetLineWidth(3);
        hmc_sgn[k]->SetFillStyle(3004);
        hmc_sgn[k]->SetFillColor(kBlue);
        hdata_sgn[k]->SetLineColor(kRed);
        hdata_sgn[k]->SetLineWidth(3);
        hdata_sgn[k]->SetFillStyle(3005);
        hdata_sgn[k]->SetFillColor(kRed);

        hdata_sgn_plus[k]->SetLineColor(kBlack);
        hdata_sgn_plus[k]->SetFillStyle(3001);
        hdata_sgn_plus[k]->SetLineWidth(2);
        hdata_sgn_plus[k]->SetLineStyle(2);

        hdata_sgn_minus[k]->SetLineColor(kBlack);
        hdata_sgn_minus[k]->SetFillStyle(3001);
        hdata_sgn_minus[k]->SetLineWidth(2);
        hdata_sgn_minus[k]->SetLineStyle(3);

        hmc_sgn[k]->SetTitle(""); // Remove the title
        hmc_sgn[k]->Draw("hist");
        hdata_sgn[k]->Draw("hist same");
        hdata_sgn_plus[k]->Draw("hist same");
        hdata_sgn_minus[k]->Draw("hist same"); 

        hmc_sgn[k]->SetStats(0);

        Double_t x_low = 0.45, x_high = 0.89, y_low = 0.50, y_high = 0.89; //top right 
        if(varname.Contains("mu_timeAtIpInOutErr") || varname.Contains("Chi2") || varname.Contains("chi2") || 
          varname.Contains("Probability")||varname.Contains("globalDeltaEtaPhi") || varname.Contains("segmComp") ||
          varname.Contains("mu_pt") || varname.Contains("Numberofvalidpixelhits") || 
              varname.Contains("mu_eta") || varname.Contains("eval") ||
              varname.Contains("mu_sumPt03") || varname.Contains("mu_nTracks03"))
                x_low = 0.12, x_high = 0.57, y_low = 0.50, y_high = 0.89; //top left 
        TLegend*leg = new TLegend(x_low, y_low, x_high, y_high);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.040);
        //leg->SetHeader(varlable, "L");
        leg->AddEntry(hmc_sgn[k],"D_{s}#rightarrow#phi(#mu#mu)#pi MC","f");
        leg->AddEntry(hdata_sgn[k],"data "+run_lable[ncut-1-1]+" (SB subtracted)","f");
        leg->AddEntry(hdata_sgn_plus[k],"data "+run_lable[ncut-1-1]+" (SB +10\% subt.)","f");
        leg->AddEntry(hdata_sgn_minus[k],"data "+run_lable[ncut-1-1]+" (SB -10\% subt.)","f");
        leg->Draw();

        //extratext
        std::string strLumi = "37.9 \\text{ fb}^{-1} (13 \\text{ TeV})";
        //TLatex* cmslabel = new TLatex(0.1,0.91, "#bf{CMS Preliminary}");
        TLatex* cmslabel = new TLatex(0.1,0.91, "#bf{CMS} work in progress");
        cmslabel->SetTextFont(42);
        cmslabel->SetTextSize(0.045); 
        cmslabel->SetNDC(kTRUE);
        cmslabel->Draw("same");
        TLatex* text = new TLatex(0.65,0.91, strLumi.c_str());
        text->SetTextFont(42);
        text->SetTextSize(0.045); 
        text->SetNDC(kTRUE);
        text->Draw("same");

        //K-S consistency test
        Double_t KS = hdata_sgn[k]->KolmogorovTest(hmc_sgn[k]);
        std::stringstream stream_KS;
        stream_KS << std::fixed << std::setprecision(3) << KS;
        std::string strKS = stream_KS.str();
        TString KStstring = "K-S test = "+strKS;
        TLatex* text_KS = new TLatex(x_low+0.02,y_low-0.05, KStstring);
        text_KS->SetTextSize(0.04);
        text_KS->SetNDC(kTRUE);
        //text_KS->Draw("same");

        // lower plot will be in pad2
        c2->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.25);
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
        h_x_ratio->GetYaxis()->SetRangeUser(0.0,2);
        h_x_ratio->SetLineColor(kBlack);
        h_x_ratio->GetYaxis()->SetTitleSize(22);
        h_x_ratio->GetYaxis()->SetTitleFont(43);
        h_x_ratio->GetYaxis()->SetTitleOffset(1.7);
        h_x_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_x_ratio->GetYaxis()->SetLabelSize(15);

        // X axis ratio plot settings
        h_x_ratio->GetXaxis()->SetTitle(varlable);
        h_x_ratio->GetXaxis()->SetTitleSize(24);
        h_x_ratio->GetXaxis()->SetTitleFont(43);
        h_x_ratio->GetXaxis()->SetTitleOffset(3.0);
        h_x_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_x_ratio->GetXaxis()->SetLabelSize(15);
        //Compute weighted average ratio
        Double_t mean = 0;
        Double_t std_dev = 0;
        for (int c=1; c<=(h_x_ratio->GetNbinsX()); c++)
        {
            if(h_x_ratio->GetBinContent(c) == 0 || h_x_ratio->GetBinError(c) == 0) continue;
            //cout<<c<<" "<<h_x_ratio->GetBinContent(c)<<" +- "<<h_x_ratio->GetBinError(c)<<endl;
            mean = mean + h_x_ratio->GetBinContent(c) / ( h_x_ratio->GetBinError(c) * h_x_ratio->GetBinError(c) );
            std_dev = std_dev + 1/( h_x_ratio->GetBinError(c) * h_x_ratio->GetBinError(c));
        }
        mean = mean/std_dev;
        std_dev = 1/std_dev;
        //Get mean value and error of ratio plot
        cout<<var[k]+run_lable[ncut-1-1]+" Mean: "<<mean<<endl;
        cout<<var[k]+run_lable[ncut-1-1]+" StdDev: "<<std_dev<<endl;
        //Draw line corresponding to mean value on ratio plot
        TLine l;
        h_x_ratio->Draw("ep");
        l.DrawLine(h_x_ratio->GetXaxis()->GetXmin(),mean,h_x_ratio->GetXaxis()->GetXmax(),mean);
        h_x_ratio->Draw("same");

        c2->cd();
        c2->Update();
        varname = varname.ReplaceAll(".", "_");
        varname = varname.ReplaceAll(":", "_");
        varname = varname.ReplaceAll(">", "_");
        varname = varname.ReplaceAll("<", "_");
        varname = varname.ReplaceAll("?", "_");
        varname = varname.ReplaceAll("(", "_");
        varname = varname.ReplaceAll(")", "_");
        varname = varname.ReplaceAll("/", "_");
        h_x_ratio->SetName(varname+"_"+TMVA_inputpath_control+"_"+run_lable[ncut-1-1]);
        h_x_ratio->Write();
        c2->SaveAs("plots_AN/"+varname+"_"+TMVA_inputpath_control+"_"+run_lable[ncut-1-1]+".png");
        fout->cd();
        fout->WriteObject(c2,varname+run_lable[ncut-1-1]);

    }
 
}
