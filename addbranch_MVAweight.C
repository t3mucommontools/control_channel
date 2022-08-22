
#include "TH1F.h"
#include <cmath>
#include <string> 

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
        //cout<<"SF="<<val.value<<" err="<<val.error<<endl;
    }else{
        val.value = 1;
        val.error = 0;
    }
    return val;
}

void addbranch_MVAweight(TString inputfile, TString year) 
{
    //open input files //weights
    TString SF_filename = "/eos/user/f/fsimone/Tau23Mu_anatools/MVAweight_dsphipi/control_channel/"+year+"/plots_MVAcut/dsphipi_MVAglb_"+year+"_output.root";
    TFile *f2 = new TFile(SF_filename, "read");
    TH2F* SF_h = dynamic_cast<TH2F*> (f2->Get("SF_mva"));
    std::cout<<"Opened input file: "<<SF_filename<<" with input TH2 "<<SF_h<<std::endl;

    //open input files //minitree
    TString tree_filename = "/eos/user/f/fsimone/input_files/"+inputfile;
    TFile *f = new TFile(tree_filename, "read");
    TTree *t = (TTree*)f->Get("outputTree");
    std::cout<<"Opened input file: "<<tree_filename<<std::endl;

    //create output file //minitree
    TString tree_newfilename = "/eos/user/f/fsimone/input_files/minitree_"+year+"_MVAweights.root";
    TFile *fnew = new TFile(tree_newfilename, "recreate");
    TTree *tnew = t->CloneTree();
    std::cout<<"Created output file: "<<tree_newfilename<<std::endl;
    f->Close();
    fnew->cd();

    TBranch* br_MVA1 = (TBranch*)tnew->GetListOfBranches()->FindObject("MuonIDeval_Mu1.MuonID");
    TBranch* br_MVA2 = (TBranch*)tnew->GetListOfBranches()->FindObject("MuonIDeval_Mu2.MuonID");
    TBranch* br_MVA3 = (TBranch*)tnew->GetListOfBranches()->FindObject("MuonIDeval_Mu3.MuonID");
    if(!(br_MVA1 && br_MVA2 && br_MVA3)) { std::cout<<"ERROR: not existing branches"<<std::endl; return; }

    Double_t evt=0;
    Double_t isMC=0;
    Double_t mu1_eta=0;
    Double_t mu2_eta=0;
    Double_t mu3_eta=0;
    Double_t mu1_mvaglb=0;
    Double_t mu2_mvaglb=0;
    Double_t mu3_mvaglb=0;

    tnew->SetBranchAddress("evt",&evt);
    tnew->SetBranchAddress("isMC",&isMC);
    tnew->SetBranchAddress("Etamu1",&mu1_eta);
    tnew->SetBranchAddress("Etamu2",&mu2_eta);
    tnew->SetBranchAddress("Etamu3",&mu3_eta);
    tnew->SetBranchAddress("MuonIDeval_Mu1.MuonID",&mu1_mvaglb);
    tnew->SetBranchAddress("MuonIDeval_Mu2.MuonID",&mu2_mvaglb);
    tnew->SetBranchAddress("MuonIDeval_Mu3.MuonID",&mu3_mvaglb);

    //new branches
    Double_t mu1_MVAweight = 0;
    Double_t mu2_MVAweight = 0;
    Double_t mu3_MVAweight = 0;
    Double_t mu1_MVAweighterr = 0;
    Double_t mu2_MVAweighterr = 0;
    Double_t mu3_MVAweighterr = 0;

    auto mu1_MVAweight_b = tnew->Branch("mu1_MVAweight", &mu1_MVAweight, "mu1_MVAweight/D");
    auto mu2_MVAweight_b = tnew->Branch("mu2_MVAweight", &mu2_MVAweight, "mu2_MVAweight/D");
    auto mu3_MVAweight_b = tnew->Branch("mu3_MVAweight", &mu3_MVAweight, "mu3_MVAweight/D");
    auto mu1_MVAweighterr_b = tnew->Branch("mu1_MVAweighterr", &mu1_MVAweighterr, "mu1_MVAweighterr/D");
    auto mu2_MVAweighterr_b = tnew->Branch("mu2_MVAweighterr", &mu2_MVAweighterr, "mu2_MVAweighterr/D");
    auto mu3_MVAweighterr_b = tnew->Branch("mu3_MVAweighterr", &mu3_MVAweighterr, "mu3_MVAweighterr/D");

    for(int i=0; i<tnew->GetEntries(); i++){ 
         tnew->GetEntry(i);
         sfMuon SF1 = GetMuonSF(SF_h, mu1_mvaglb, mu1_eta);
         sfMuon SF2 = GetMuonSF(SF_h, mu2_mvaglb, mu2_eta);
         sfMuon SF3 = GetMuonSF(SF_h, mu3_mvaglb, mu3_eta);

         mu1_MVAweight = (isMC>0&&isMC<4) ? SF1.value : 1.0;
         mu2_MVAweight = (isMC>0&&isMC<4) ? SF2.value : 1.0;
         mu3_MVAweight = (isMC>0&&isMC<4) ? SF3.value : 1.0;
         mu1_MVAweighterr = (isMC>0&&isMC<4) ? SF1.error : 0.0;
         mu2_MVAweighterr = (isMC>0&&isMC<4) ? SF2.error : 0.0;
         mu3_MVAweighterr = (isMC>0&&isMC<4) ? SF3.error : 0.0;

         mu1_MVAweight_b->Fill();
         mu2_MVAweight_b->Fill();
         mu3_MVAweight_b->Fill();
         mu1_MVAweighterr_b->Fill();
         mu2_MVAweighterr_b->Fill();
         mu3_MVAweighterr_b->Fill();
    }
    f2->Close();
    tnew->Write("", TObject::kOverwrite);
    fnew->Close();
}
