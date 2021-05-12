/* * * * * * * * * * * *
 * GetPede.cc          *
 * 2020. 08. 27 (Sat)  *
 * T. Fujiwara         *
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
#include "TRandom.h"
#include "TStyle.h"
#include "TObject.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
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
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

#include "./Setting.h"


#define NofMPPC   2
#define NofScinti 2 

/*-----------path------------*/
string const PATH_param_p = "../param/pede/";
string const PATH_param_q = "../param/qdc/";
string const PATH_param_t = "../param/twlk/";


TString const PATH_dat         = "../root/";
TString const PATH_fig_General = "../fig/";
TString const PATH_fig_l       = "../fig/leng/";
TString const PATH_fig_q       = "../fig/qdc/";
TString const PATH_fig_p       = "../fig/pede/";
TString const PATH_fig_g       = "../fig/gain/";


//TString const Label[NofMPPC][3] =  { {"Ref1.QDC", "Ref1.TDC", "Ref1"},
//                                     {"Ref2.QDC", "Ref2.TDC", "Ref2"},
//                                     {"ToF1.QDC", "ToF1.TDC", "ToF1"},
//                                     {"ToF2.QDC", "ToF2.TDC", "ToF2"}};
TString const Label[NofMPPC][3] =  { {"Ref1 QDC", "Ref1 TDC", "Ref1"},
                                     {"Ref2 QDC", "Ref2 TDC", "Ref2"}};

int const ChID[NofMPPC] = {0, 1};
int const ClockCh = 15;

/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
class GetPede
{
	public:
		 GetPede();
		~GetPede();
		void SetRoot(int RNum);
		void DefCanv();
		void DefObj();
		void FillHist();
		void DrawHist();
		void FitHist();
		void Export();

	private:
		Setting *Set;
		TFile *ifp;
		TTree *tree;
		TCanvas *ca;
		TH1D *h_Pede[NofMPPC];
		TF1 *f_Pede[NofMPPC];
		TFrame *Fr[NofMPPC];
		TLatex *Lat;

		int ItNum;
		int RunNum;
		string RunName_Full;
		string RunName_str;
		TString RunName_tst;
		TString figure;
		TString figure_open, figure_close;
		string param;
		
		int EvN;
		int TDC[32];
		int QDC[32];
		int Entry = 0;

		bool Trig_Clock;
	
		int NofIt;
		double par_tmp;
		double FitRange;
		double FitPar[3];
		double PedeCenter[NofMPPC][2];	//0: Value, 1: Error
		double PedeWidth[NofMPPC][2]; 	//0: Value, 1: Error
		
};

////////////////////////////////////////////////////////////////////////////////////////////////////
GetPede::GetPede(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	
	Set = new Setting();
	Set -> Setting_Gene(1);

	Trig_Clock=false;
	NofIt=5;
//	FitRange = 10.;
	FitRange = 5.;
}
//------------------------------------------------------------------------------------------------//
GetPede::~GetPede(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
}
//------------------------------------------------------------------------------------------------//
void GetPede::SetRoot(int RNum){
	cout<<"==========     GetPede::SetRoot     =========="<<endl;
	RunNum = RNum;
	RunName_Full = Form("run%d.root", RunNum);
	RunName_Full = PATH_dat+RunName_Full;
	RunName_str  = Form("run%d.root", RunNum);
	RunName_tst  = Form("run%d"     , RunNum);
	figure = Form("Pedestal_run%04d_%03d.pdf", RunNum, ItNum);
	figure = PATH_fig_p+figure;
	figure_open  = figure+"[";
	figure_close = figure+"]";
	param = Form("Pedestal_run%04d_%03d.dat", RunNum, ItNum);
	param = PATH_param_p+param;
	
	ifp = new TFile( RunName_Full.c_str(), "READONLY" );
	tree = (TTree*)ifp->Get("tree");
	tree -> SetBranchAddress( "evnum", &EvN );
	tree -> SetBranchAddress( "qdc", &QDC );
	tree -> SetBranchAddress( "tdc", &TDC );
	Entry = tree->GetEntries();
	cout<<"======================================"<<endl;
	cout<<"    FILE IMPORT: "<<RunName_str<<endl;
	cout<<"    Entry: "<<Entry<<endl;
	cout<<"======================================"<<endl;
}
//------------------------------------------------------------------------------------------------//
void GetPede::DefCanv(){
	cout<<"==========     GetPede::DefCanv     =========="<<endl;
	ca = new TCanvas( "ca", "ca", 1600, 800);
	ca -> Divide(2,2);
}
//------------------------------------------------------------------------------------------------//
void GetPede::DefObj(){
	cout<<"==========     GetPede::DefObj     =========="<<endl;
	for(int i=0; i<NofMPPC; i++){
		h_Pede[i] = new TH1D( Form("h_Pede[%d]",i), Form("h_Pede[%d]",i), 300, 0, 300 );
		Set->Setting_Hist1D( h_Pede[i], RunName_tst+": "+Label[i][2], Label[i][0]+" [ch/0.1 pC]", "Counts/ch", 602, 1, 62, kCyan-10, 1003 );
		
		f_Pede[i] = new TF1( Form("f_Pede[%d]",i), "gaus(0)", 0., 200. );
		Set -> Setting_Func( f_Pede[i], 6, 1 );
	}

	Lat = new TLatex();
	Set -> Setting_Latex( Lat, 62, 22, 6, 0.05 );
}
//------------------------------------------------------------------------------------------------//
void GetPede::FillHist(){
	cout<<"==========     GetPede::FillHist     =========="<<endl;
	for(int i=0; i<Entry; i++){
		tree -> GetEntry(i);
		if(TDC[ClockCh]>0){
			Trig_Clock = true;
		}else Trig_Clock = false;
		if(Trig_Clock){
			for(int j=0; j<NofMPPC; j++){h_Pede[j]->Fill(QDC[j]);}
		}else;
	}
}
//------------------------------------------------------------------------------------------------//
void GetPede::DrawHist(){
	cout<<"==========     GetPede::DrawHist     =========="<<endl;
	ca->cd(1);
	for(int i=0; i<NofMPPC; i++){
		ca->cd(i+1);
		gPad->SetLogy(1);
		h_Pede[i]->Draw("");
		gPad->Update();
	}
}
//------------------------------------------------------------------------------------------------//
void GetPede::FitHist(){
	cout<<"==========     GetPede::FitHist     =========="<<endl;
	for(int i=0; i<NofMPPC; i++){
		ca->cd(i+1);
		par_tmp = h_Pede[i]->GetBinCenter( h_Pede[i]->GetMaximumBin() );
		for(int j=0; j<NofIt; j++){
			h_Pede[i] -> Fit( f_Pede[i], "", "", par_tmp-FitRange, par_tmp+FitRange );
			f_Pede[i]->GetParameters(FitPar);
			par_tmp = FitPar[1];
		}
		PedeCenter[i][0] = f_Pede[i]->GetParameter(1);
		PedeCenter[i][1] = f_Pede[i]->GetParError(1);
		PedeWidth[i][0]  = f_Pede[i]->GetParameter(2);
		PedeWidth[i][1]  = f_Pede[i]->GetParError(2);
		Lat -> DrawLatexNDC( .75, .50, Form("Pedestal: %.1lf #pm %.1lf [ch]", PedeCenter[i][0], PedeWidth[i][0] ) );
		gPad->Update();
		Fr[i] = gPad->GetFrame();
		Fr[i] -> SetFillStyle(0);
		Fr[i] -> Draw("same");
	}
}
//------------------------------------------------------------------------------------------------//
void GetPede::Export(){
	cout<<"==========     GetPede::Export     =========="<<endl;

	ca->Print( figure_open, "pdf" );
	ca->Print( figure );
	ca->Print( figure_close, "pdf" );

	ofstream ofs( param, ios_base::trunc );
	for(int i=0; i<NofMPPC; i++){ofs<<PedeCenter[i][0]<<"   "<<PedeWidth[i][0]<<endl;}
}
//------------------------------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv){
	int Num;
	
	if(argc==2){
		Num=atoi(argv[1]);
	}else{
		cout<<"  !!!WARNING!!!   "<<endl;
		cout<<"  This program needs one argument to specify runnumber."<<endl;
		cout<<"  Please input runnumber(Only number) >>>"<<endl;
		cin>>Num;
	}
	cout<<"-+-+-+-+-+-+-+-"<<endl;
	cout<<Form("Target of Analysys >> run%d.root", Num)<<endl;
	cout<<"-+-+-+-+-+-+-+-"<<endl;


	TApplication *theApp;
	GetPede *ana;
	
	theApp = new TApplication( "App", &argc, argv );

	ana    = new GetPede();
	ana -> SetRoot(Num);
	ana -> DefCanv();
	ana -> DefObj();
	ana -> FillHist();
	ana -> DrawHist();
	ana -> FitHist();
	ana -> Export();
	delete ana;

	//====	END of Job	====//
	gSystem->Exit(1);
	theApp->Run();

	return 0;
}
