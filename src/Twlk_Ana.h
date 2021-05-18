/* * * * * * * * * * * *
 * Twlk_Ana.h          *
 *	2021. 03. 27 (Sat) *
 *	T. FUJIWARA        *
 * * * * * * * * * * * */
#ifndef Twlk_Ana_h
#define Twlk_Ana_h 1

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
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

//////////	PATH

/////////////
typedef struct
{
  Int_t tdc, qdc; // [ch]
  Double_t t;     // [nsec]
  Double_t e;     // arbitrary
} TQ;
/////////////
const string Label[4] = {"Ref1", "Ref2",
                         "ToF1", "ToF2"};

//////////	PATH
const string PATH_param_p = "../param/pede/";
const string PATH_param_t = "../param/twlk/";

const string PATH_RootFile = "../root/";
/////////////
const int NofCanv = 18;
const int NofComb = 3;
const int NofMPPC = 4;
const int NofDet  = 2;
const int CHassign[NofDet]={0, 1};
	//0: Ref12
	//1: Ref22


class Twlk_Ana{
	public:
		 Twlk_Ana(int rnum, int phctype_1=8, int phctype_2=8);
		~Twlk_Ana();
		void TwlkMan();
		void SetTwlkRef1(int type=0);
		void SetTwlkRef2(int type=0);
		void SetRoot();
		void ExtractEvent();
		void ImportPedestal();
		void DefineCanv();
		void DefineObj();
		void FillHist();
		void Draw();
		void Fit();
		void Fit_Ref();
		void SearchBest();
		void DecTOF();
		void Export();
		
	//	int GetTwlkMode(){return TwlkMode;}
		int GetTwlkMode(int id=0);
		double Twlk_Correction( double x, double* parameter, double Coeff, int Mode=0 );

		bool Trig_Ref1FToF1(int* t);
		bool Trig_FToF1Ref2(int* t);
		bool Trig_Ref1Ref2(int* t);

	private:
		int RunNum;

		TFile* ifp;
		TTree* tree;
		TChain* chain;
		TCanvas* Ca[NofCanv];
		TCanvas* Ca_dummy;
		TH1D* h_1d_qdc_full[NofDet];
		TH1D* h_1d_qdc_wcut[NofDet];
		TH1D* h_1d_qdc_qcut[NofDet];
		TH1D* h_1d_tdc_full[NofDet];
		TH1D* h_1d_tdc_wcut[NofDet];
		TH2D* h_2d_rawtq[NofDet];
		TH2D* h_2d_rawtq_fit[NofDet];
		TH2D* h_2d_dectq_full[NofDet];
		TH2D* h_2d_dectq_wcut[NofDet];
		TH2D* h_2d_dectq_tmp[NofDet];
		TH1D* h_1d_rawtq_slice[NofDet];
		TH1D* h_1d_dectq_slice_tmp[NofDet];
		TH1D* h_1d_rawtof;
		TH1D* h_1d_dectof_full;
		TH1D* h_1d_dectof_wcut;
		TH1D* h_1d_dectof_tmp;
		TH2D* h_2d_twlkfact;
		TH2D* h_2d_twlkfact_Chis;
		TH2D* h_2d_twlkfact_RMS;
		TH2D* h_2d_twlkfact_Skew;
		TH2D* h_2d_twlkfact_ChisConsFit;
		TGraph* gr_twlkfact;
		TGraph* gr_twlkfact_Chis;

		TFrame* fr_2d_rawtq[NofDet];
		TFrame* fr_2d_rawtq_fit[NofDet];
		TFrame* fr_2d_dectq_full[NofDet];
		TFrame* fr_2d_dectq_wcut[NofDet];
		TFrame* fr_2d_twlkfact;
		TFrame* fr_2d_twlkfact_Chis;

		TF1* f_twlk[NofDet];
		TF1* f_tof[3];
		TF1* f_tof_tmp;
		TF1* f_const_tmp[NofDet];

		TLine* Ln_RefQDCCut[NofDet];

		TLegend* Leg_qdc[NofDet];
		TLegend* Leg_tdc[NofDet];
		TLatex* Lat;
		TPaveText* Pt[NofDet];

		TFile *ofp_obj, *ofp_dat;
		TTree* tree_Ex;
		ofstream ofs_twlk;
		ofstream ofs_FitSlicesPoints[NofDet];

		Setting* MySetting;

		string RootFile;
		string Figure;
		string Figure_Open;
		string Figure_Close;
		string Param_Twlk;
		string Param_Pede;
		string RootFile_Ex;
		string RootFile_Ex_obj;
		string RootFile_Ex_dat;

		TQ Ref[2];
			//0: Ref1
			//1: Ref2
		TQ ToF1;
		TQ ToF2;
	
		int Ent;
		int EvNum;
		int tdc[32];
		int qdc[32];
		double decTime[NofDet];
		double rtof;
		double dtof;
		vector<int> EffectiveEvent;
		int EffectiveEvCounter=0;
		int MaxEventItr=0;

		ifstream ifs_pede;
		double pede_v[NofMPPC];
		double pede_w[NofMPPC];

		int TwlkTDCGate[NofDet][2] = { {1900, 2200}, {1900, 2200}};
		double TwlkTOFAxisLim[2] ={-7., 7.};

		int TwlkMode=0;
		int Ref1TwlkMode=0;
		int Ref2TwlkMode=0;
		bool TwlkInitialParSetFlag=false;
		bool TwlkInitialParSetFlag_Ref[NofDet]={false, false};
		double TwlkInitialFitPar[NofDet][9];

		double TwlkQDCPeak[NofDet];
		double TwlkQDCCutForFit[NofDet][2]={{30.,70.},{50.,160.}};

		TString TwlkFitSliceOpt[NofDet] = {"QNRG4", "QNRG4"};
		TString TwlkFuncType[NofDet];
		int TwlkNPar[NofDet];
//		double TwlkFitMin[NofDet] = {1.E-2, 1.E-2};
		double TwlkFitMin[NofDet] = {6., 6.};
		double TwlkFitMax[NofDet] = {250., 250.};
		double Twlk[NofDet][9];
		double TwlkFactor[NofDet];
		double TwlkFact_Tmp[NofDet];
		double TwlkFact_Best[NofDet];
		double TwlkSigVal_Tmp;
		double TwlkSigVal_Best=999.;
		double TwlkSigErr_Tmp;
		double TwlkSigErr_Best=999.;
		double TwlkChisNDF_Tmp;
		double TwlkChisNDF_Best=99999.;
		double TwlkRMS_Tmp;
		double TwlkRMS_Best = 9999.;
		double TwlkSkew_Tmp;
		double TwlkSkew_Best = 9999.;
		double TwlkChisNDF_ConsFit_Sep_Tmp[NofDet];
		double TwlkChisNDF_ConsFit_Sep_Best[NofDet]={9999., 9999.};
		double TwlkChisNDF_ConsFit_Ave_Tmp;
		double TwlkChisNDF_ConsFit_Ave_Best=9999.;


		double TwlkDelta=0.10;

		double TwlkRefQDCCutLim[NofDet]={20., 20.};
		
		int TwlkQDCNBin=205;
		double TwlkQDCMax=400.;

		//For search best factor
		double Par_tmp;
		int FitItrTmp=3;
		double FitWid_Tmp=2.0;

		int EvNum_Ex;
		int tdc_Ex[32];
		int qdc_Ex[32];
		TQ Ex_Ref1;
		TQ Ex_Ref2;
};
#endif
