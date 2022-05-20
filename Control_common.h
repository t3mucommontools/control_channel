#include <iostream>

using namespace std;
    TString workdir_control = "/lustrehome/fsimone/MVA_2018/MVA_control/";
 
//TMVA Training options
    TString TMVA_outputpath_control = "MVA_control_2018_6nov_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    TString method_control = "BDT";
    TString TMVA_weightfilename_control = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath_control = "MVA_control_2018_6nov_";  //name to load TMVA results for evaluation

//data rootfiles

    TString inputpath_datarun_control[] = {
           ////Enlarged phi -> mumu mass window 0.97 - 1.07
           ////added dphi2D
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1657/AnalysedTree_data_2018A_control_30july.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1658/AnalysedTree_data_2018B_control_30july.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1659/AnalysedTree_data_2018C_control_30july.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1700/AnalysedTree_data_2018D_control_30july.root",
           //Enlarged phi -> mumu mass window 0.95 - 1.09
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1024/AnalysedTree_data_2018A_control_16june.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1025/AnalysedTree_data_2018B_control_16june.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1026/AnalysedTree_data_2018C_control_16june.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1027/AnalysedTree_data_2018D_control_16june.root",
           ////UL //DoubleMu0 or DoubleMu4 - added TreeMu1 TreeMu2 - added dxy dxyErr
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210605_1119/AnalysedTree_data_2018A_control_5june.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210605_1120/AnalysedTree_data_2018B_control_5june.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210605_1121/AnalysedTree_data_2018C_control_5june.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210605_1122/AnalysedTree_data_2018D_control_5june.root",
           ////UL //DoubleMu0 or DoubleMu4 - added TreeMu1 TreeMu2
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210531_2142/AnalysedTree_data_2018A_control_31may.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210531_2143/AnalysedTree_data_2018B_control_31may.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210531_2144/AnalysedTree_data_2018C_control_31may.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210531_2145/AnalysedTree_data_2018D_control_31may.root",
           ////UL //DoubleMu0 or DoubleMu4 (not tirple)
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210222_1744/AnalysedTree_data_2018A_control_UL_22feb.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210222_1745/AnalysedTree_data_2018B_control_UL_22feb.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210222_1746/AnalysedTree_data_2018C_control_UL_22feb.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210222_1747/AnalysedTree_data_2018D_control_UL_22feb.root"
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210521_1035/AnalysedTree_data_2018D_control_UL_21may.root"
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210528_2221/AnalysedTree_data_2018D_control_28may.root"

           ////DsPhiPi 2018 Data - no cut on SV chi2 - isGlobal+isPF - prescaled HLT //only DoubleMu0
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2155/AnalysedTree_data_2018A_control_14oct.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2156/AnalysedTree_data_2018B_control_14oct.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2157/AnalysedTree_data_2018C_control_14oct.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2158/AnalysedTree_data_2018D_control_14oct.root"
           ////DsPhiPi 2018 Data - no cut on SV chi2 - isGlobal+isPF - prescaled HLT
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1033/AnalysedTree_data_2018A_control_5ott.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1033/AnalysedTree_data_2018B_control_5ott.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1558/AnalysedTree_data_2018C_control_5ott.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1033/AnalysedTree_data_2018D_control_5ott.root"
           };

    //UL //DoubleMu0 or DoubleMu4
    //Enlarged phi -> mumu mass window 0.97 - 1.07 //added dphi2D
    TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210730_1656/AnalysedTree_MC_2018DsPhiPi_control_30july.root";
    //Enlarged phi -> mumu mass window 0.95 - 1.09
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210616_1023/AnalysedTree_MC_2018DsPhiPi_control_16june.root";
    //added TreeMu1 TreeMu2  - added dxy dxyErr
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210605_1118/AnalysedTree_MC_2018DsPhiPi_control_5june.root";
    ////added TreeMu1 TreeMu2
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210531_2142/AnalysedTree_MC_2018DsPhiPi_control_31may.root";
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210528_1635/AnalysedTree_MC_2018DsPhiPi_control_28may.root";
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210521_1037/AnalysedTree_MC_2018DsPhiPi_control_UL_21may.root";
    //UL //DoubleMu0 or DoubleMu4 (not tirple)
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20210222_1743/AnalysedTree_MC_2018DsPhiPi_control_UL_22feb.root";
    ////DsPhiPi 2018 MC - no cut on SV chi2 - isGlobal+isPF - prescaled HLT //only DoubleMu0
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201014_2154/AnalysedTree_MC_2018DsPhiPi_control_14oct.root";
    ////DsPhiPi 2018 MC - no cut on SV chi2 - isGlobal+isPF
    //TString inputpath_DsPhiPi = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20201005_1032/AnalysedTree_MC_2018DsPhiPi_control_5ott.root";


//Coefficients for signal normalisation
// Norm = Lumi_data[i]*xsection_mc*BR/N_MC;
    TString run_name_control[] = {"2018A", "2018B", "2018C", "2018D"};
    Double_t Lumi_data_control[] = {0.7025, 0.3522, 0.3458, 1.58803}; //recorded lumi by HLT_DoubleMu3_Trk_*
    //Double_t Lumi_data_control[] = {13.98, 7.06, 6.90, 31.75};

    ////rereco
    //Double_t xsection_mc = 2.32e10; //Ds Production Cross section
    //int N_MC = 2050462;  //rereco Total number of events in MC sample

    //UL
    Double_t xsection_mc = 1.06e10; //Ds Production Cross section
    int N_MC = 2964234;  //rereco Total number of events in MC sample

    Double_t BR = 1.29e-5;  //Branching ratio Ds to Phi Pi


//TMVA settings
    //common preselections
    TCut common_cut_control = "bs_sv_d2Dsig>2.0 && Ptmu3 > 1.2 && " 
                              "((Ptmu1>3.5 && Etamu1<1.2) || (Ptmu1>2.0 && Etamu1>=1.2 && Etamu1<=2.4)) && "
                              "((Ptmu2>3.5 && Etamu2<1.2) || (Ptmu2>2.0 && Etamu2>=1.2 && Etamu2<=2.4)) && "
                              "abs(phiMass-1.02)<0.045 && "
                              "!(l1double_DoubleMu4_fired && !l1double_DoubleMu0_fired)";
    //input variables
    TString BDTinVar_control = "BDTinputVar_control_v4.txt";
    TString BDTspecVar_control = "../BDTspecVar_2018_noevt.txt";
//Utility to read variables from text file
void readVarName_control(std::vector<TString> &var_name, std::vector<TString> &var_def,  TString filename ){
     TString var[2];
     ifstream inputVar;
     inputVar.open(filename);
     while (!inputVar.fail() && !inputVar.eof()){
         inputVar >> var[0] >> var[1];
         if(! (var[0]=="" || var[1]=="" )){
             var_name.push_back(var[0]);
             var_def.push_back(var[1]);
         }
     }
     inputVar.close();
}
