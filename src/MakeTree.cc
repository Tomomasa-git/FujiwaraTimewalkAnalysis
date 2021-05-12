/* * * * * * * * * * * *
 * MakeTree.cc         *
 *	2021. 03. 27 (Sat) *
 *	T. FUJIWARA        *
 * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;

#include "TApplication.h"
#include "TApplicationImp.h"
#include "TCanvasImp.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"
#include "TObject.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TCurlyLine.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TArc.h"
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./Setting.h"

/////////////
typedef struct
{
  Int_t tdc, qdc; // [ch]
  Double_t t;     // [nsec]
  Double_t e;     // arbitrary
} TQ;
/////////////
const string Label[10] = {"Ref1-1", "Ref1-2", "Ref1-3",
                          "Ref2-1", "Ref2-2", "Ref2-3",
                          "ToF1L" , "ToF1R" ,
                          "ToF2L" , "ToF2R"};

//////////	PATH
const string PATH_param_p = "../param/pede/";
const string PATH_param_t = "../param/twlk/";

const string PATH_RootFile = "../root/";
const int NofCanv = 8;
const int NofDet  = 2;
const int CHassign[NofDet]={1, 4};
const int NofRun=99;
/////////////

class MakeTree{
	public:
		 MakeTree(int runnum=1068);
		~MakeTree();
		void SetRoot();
		void ImportPedestal();
		void Read();
		void Export();

		int CallData(int i);

	private:
		TFile* ifp;
		TTree* tree;
		
		TFile* ofp;
		TTree* tree_Ex;

		Setting* MySetting;

		int Run;
		int EffectiveRunNum;
		string Root_In;
		string Root_Ex;

		//Input from original rootfiles		
		int ENum;
		int EvSum;
		int EvNum_In;
		int tdc_In[32];
		int qdc_In[32];

		//Export new rootfiles
		int EvTot;
		int EvNum_Ex;
		int tdc_Ex[32];
		int qdc_Ex[32];
		TQ Ref1;	//ch0
		TQ Ref2;	//ch1
		TQ ToF1;	//ch6
		TQ ToF2;	//ch8
		int Clock;	//ch15		

		double TDCReso=0.035;	//V775 param. 1 ch. = 0.035 ns = 35 ps

		//Pedestal Param
		string Param_In;
		ifstream ifs;

		//Mean of Pedestal (ch/0.1 pC)
		double Pedestal_Val_Ref1;
		double Pedestal_Val_Ref2;
		double Pedestal_Val_ToF1;
		double Pedestal_Val_ToF2;

		//Width of Pedestal (ch/0.1 pC)
		double Pedestal_Wid_Ref1;
		double Pedestal_Wid_Ref2;
		double Pedestal_Wid_ToF1;
		double Pedestal_Wid_ToF2;
};

/////////////
/////////////
/////////////
//----------------------------------------------------------
MakeTree::MakeTree(int runnum){
	MySetting = new Setting();
	MySetting->Setting_Gene();
	
	Run=runnum;

	Root_Ex = Form("../root/wararoot/wararoot%04d.root",Run);
	ofp = new TFile(Root_Ex.c_str(), "recreate");
	tree_Ex = new TTree("tree", "tree_for_analysis_Reference-det.");
	tree_Ex->Branch( "evnum", &EvNum_Ex, "NumEvent/I"                );
	tree_Ex->Branch( "tdc"  , tdc_Ex   , "V775_1[32]/I"              );
	tree_Ex->Branch( "qdc"  , qdc_Ex   , "V792_1[32]/I"              );
	tree_Ex->Branch( "Ref1" , &Ref1    , "TDC/I:QDC/I:Time/D:Edep/D" );
	tree_Ex->Branch( "Ref2" , &Ref2    , "TDC/I:QDC/I:Time/D:Edep/D" );

//	In this setup, There is not ToF det.
//	tree_Ex->Branch( "ToF1" , &ToF1    , "TDC/I:QDC/I:Time/D:Edep/D" );
//	tree_Ex->Branch( "ToF2" , &ToF2    , "TDC/I:QDC/I:Time/D:Edep/D" );
	EvTot=0;
}

//----------------------------------------------------------
MakeTree::~MakeTree(){
	delete MySetting;
}

//----------------------------------------------------------
void MakeTree::SetRoot(){
	Root_In = PATH_RootFile+Form("run%04d.root", Run);
	ifp  = new TFile( Root_In.c_str(), "READONLY" );
	tree = (TTree*)ifp->Get("tree");
	tree -> SetBranchAddress( "tdc"  , tdc_In    );
	tree -> SetBranchAddress( "qdc"  , qdc_In    );
	tree -> SetBranchAddress( "evnum", &EvNum_In );
	ENum = tree->GetEntries();
	EvTot += ENum;
	cout<<Form("run%04d.root >>", Run)<<"  "<<ENum<<"  "<<EvTot<<endl;
}

//----------------------------------------------------------
void MakeTree::ImportPedestal(){
	Param_In = PATH_param_p+Form("Pedestal_run%04d_000.dat", Run);
	ifs.open(Param_In.c_str(), ios_base::in);
	ifs>>Pedestal_Val_Ref1>>Pedestal_Wid_Ref1;	cout<<Pedestal_Val_Ref1<<endl;
	ifs>>Pedestal_Val_Ref2>>Pedestal_Wid_Ref2;	cout<<Pedestal_Val_Ref2<<endl;
//	ifs>>Pedestal_Val_ToF1>>Pedestal_Wid_ToF1;	cout<<Pedestal_Val_ToF1<<endl;
//	ifs>>Pedestal_Val_ToF2>>Pedestal_Wid_ToF2;	cout<<Pedestal_Val_ToF2<<endl;
	cout<<"**********"<<endl;
}

//----------------------------------------------------------
void MakeTree::Read(){
	for(int i=0; i<EvTot; i++){
		tree_Ex->GetEntry(i);
		ifp->cd();
		tree->GetEntry(i);
		
		//For each detectors
		////Ref1
		Ref1.tdc = tdc_In[0];
		Ref1.qdc = qdc_In[0];
		Ref1.t   = TDCReso*tdc_In[0];
		Ref1.e   = 0.1*((double)qdc_In[0]-Pedestal_Val_Ref1);
	
		////Ref2
		Ref2.tdc = tdc_In[1];
		Ref2.qdc = qdc_In[1];
		Ref2.t   = TDCReso*tdc_In[1];
		Ref2.e   = 0.1*((double)qdc_In[1]-Pedestal_Val_Ref2);
	
		////ToF1
//		ToF1.tdc = tdc_In[2];
//		ToF1.qdc = qdc_In[2];
//		ToF1.t   = TDCReso*tdc_In[2];
//		ToF1.e   = 0.1*((double)qdc_In[2]-Pedestal_Val_ToF1);
	
		////ToF2
//		ToF2.tdc = tdc_In[3];
//		ToF2.qdc = qdc_In[3];
//		ToF2.t   = TDCReso*tdc_In[3];
//		ToF2.e   = 0.1*((double)qdc_In[3]-Pedestal_Val_ToF2);

		for(int seg=0; seg<32; seg++){
			tdc_Ex[seg] = tdc_In[seg];
			qdc_Ex[seg] = qdc_In[seg];
		}

		EvNum_Ex=i;
		
		tree_Ex->Fill();
	}
}

//----------------------------------------------------------
int MakeTree::CallData(int i){
	int Ret;
	int j=0;
	while(1){
		if(i<EvSum){
			Ret=j;
			break;
		}else;
		j++;
	}

	return Ret;
}

//----------------------------------------------------------
void MakeTree::Export(){
	ofp->cd();
	tree_Ex->Write();
	ofp->Close();
}

//----------------------------------------------------------

//----------------------------------------------------------

//----------------------------------------------------------

//----------------------------------------------------------

//----------------------------------------------------------

//----------------------------------------------------------
int main(int argc, char** argv){
	int Mode;
	if(argc>1){
		Mode=atoi(argv[1]);
	}else{
		Mode=1068;
	}

	TApplication* theApp;
	MakeTree* ana;

	theApp = new TApplication( "App", &argc, argv );
	ana = new MakeTree(Mode);
	ana->SetRoot();
	ana->ImportPedestal();
	ana->Read();
	ana->Export();

	gSystem->Exit(-1);
	theApp->Run();
	
	return 0;
}


