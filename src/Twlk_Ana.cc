/* * * * * * * * * * * *
 * Twlk_Ana.cc         *
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
#include "./Twlk_Ana.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
Twlk_Ana::Twlk_Ana(int rnum, int phctype_1, int phctype_2){
	int ItNum=0;
	MySetting = new Setting();
	MySetting->Setting_Gene();

	gStyle->SetPadRightMargin(.100);
	gStyle->SetPadLeftMargin(.120);
	gStyle->SetStatX(.875);
	gStyle->SetStatY(.900);

	RunNum = rnum;

	Ref1TwlkMode=phctype_1;
	Ref2TwlkMode=phctype_2;
	
	if(RunNum==99999){
		ItNum=1;
	}else if(RunNum==1068){
		ItNum=7;
	}else if(RunNum==1069){
		ItNum=2;
	}else if(RunNum==1070){
		ItNum=3;
	}else if(RunNum==1071){
		ItNum=1;
	}else if(RunNum==4000){
		ItNum=3;
	}else;
	
	RootFile        = Form("../root/wararoot/wararoot%04d.root", RunNum);
	Figure          = Form("../fig/twlk/TwlkAnaRef_wararoot%04d_%03d.pdf", RunNum, ItNum);
	Figure_Open     = Figure+"[";
	Figure_Close    = Figure+"]";
	RootFile_Ex_obj = Form("../root/Export/ana_wararoot%04d_%03d_obj.root", RunNum, ItNum);
	RootFile_Ex_dat = Form("../root/Export/ana_wararoot%04d_%03d_dat.root", RunNum, ItNum);
	Param_Twlk      = Form("../param/twlk/ana_wararoot%04d_%03d.dat"      , RunNum, ItNum);
//	Param_Pede = Form("../param/pede/pedestal_run%04d_%03d.dat", RunNum, ItNum);
	
	TwlkDelta=0.025;
	TwlkQDCNBin=410;

	if(RunNum==1069){
		TwlkQDCNBin=205;
		TwlkFitMax[0] = 120.;
		TwlkFitMax[1] = 380.;
		TwlkFitMin[0] = 5.;
		TwlkFitMin[1] = 20.;
		TwlkFitSliceOpt[0] = "QNRG5";
		TwlkFitSliceOpt[1] = "QNRG4";

		TwlkRefQDCCutLim[1] = 40.;
	}else if(RunNum==1070){
		TwlkFitMax[0] = 40.;
		TwlkFitMin[0] = 3.;
		TwlkFitSliceOpt[0] = "QNRG4";
		TwlkRefQDCCutLim[1] = 40.;
	}else if(RunNum==1071){
		TwlkFitMax[0] = 60.;
		TwlkFitMin[0] = 1.;
	}else if(RunNum==4000){
		TwlkQDCNBin = 205;
		TwlkFitMax[0] = 150.;
		TwlkFitMin[0] = 5.;
		TwlkFitSliceOpt[0] = "QNRG4";
		TwlkFitSliceOpt[1] = "QNRG4";
	}else;
}

//--------------------------------------------------------------------------------------------------
Twlk_Ana::~Twlk_Ana(){
	delete MySetting;
}

//--------------------------------------------------------------------------------------------------
double Twlk_Ana::Twlk_Correction( double x, double* parameter, double Coeff, int Mode ){
	double ret=0.;
	switch(Mode){
		case 0:
		ret = Coeff*parameter[0]/sqrt(x);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;
		
		case 1:
		ret = Coeff*parameter[0]*exp(parameter[1]*x);
		break;

		case 2:
		ret = Coeff*parameter[0]*log(x);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;

		case 3:
		ret = Coeff*parameter[0]*pow(x,parameter[1]);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;

		case 4:
		ret = Coeff*(parameter[0]/sqrt(x)+parameter[1]*x);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;

		case 5:
		ret = Coeff*(parameter[0]*exp(parameter[1]*x)+parameter[2]*x);
		break;

		case 6:
		ret = Coeff*(parameter[0]*log(x)+parameter[1]*x);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;

		case 7:
		ret = Coeff*(parameter[0]*pow(x,parameter[1]) + parameter[2]*x);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;

		case 8:
		ret = Coeff*parameter[0]/sqrt(x-parameter[1]);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;
		
		default:
		ret = Coeff*parameter[0]/sqrt(x);
		if(x<1.E-4){
			ret=0.;
		}else;
		break;
	}
	return ret;
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::TwlkMan(){
	SetTwlkRef1(Ref1TwlkMode);
	SetTwlkRef2(Ref2TwlkMode);
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::SetTwlkRef1(int type){
	switch(type){
		case 0:
		TwlkFuncType[0] = "[0]/sqrt(x)+[1]";
		TwlkNPar[0]=2;
		break;
		
		case 1:
		TwlkFuncType[0] = "[0]*exp([1]*x)+[2]";
		TwlkNPar[0]=3;
		break;
		
		case 2:
		TwlkFuncType[0] = "[0]*log(x)+[1]";
		TwlkNPar[0]=2;
		break;

		case 3:
		TwlkFuncType[0] = "[0]*pow(x,[1])+[2]";
		TwlkNPar[0]=3;
		TwlkInitialParSetFlag_Ref[0] = true;
		TwlkInitialFitPar[0][0] = -2.0;
		TwlkInitialFitPar[0][1] = -0.5;
		TwlkInitialFitPar[0][2] = 0.5;
		break;
		
		case 4:
		TwlkFuncType[0] = "[0]/sqrt(x)+[1]*x+[2]";
		TwlkNPar[0]=3;
		break;
		
		case 5:
		TwlkFuncType[0] = "[0]*exp([1]*x)+[2]*x+[3]";
		TwlkNPar[0]=4;
		break;

		case 6:
		TwlkFuncType[0] = "[0]*log(x)+[1]*x+[2]";
		TwlkNPar[0]=3;
		break;

		case 7:
		TwlkFuncType[0] = "[0]*pow(x,[1])+[2]*x+[3]";
		TwlkNPar[0]=4;
		TwlkInitialParSetFlag_Ref[0] = true;
		TwlkInitialFitPar[0][0] = -2.0;
		TwlkInitialFitPar[0][1] = -0.5;
		TwlkInitialFitPar[0][2] = 0.5;
		TwlkInitialFitPar[0][3] = 0.1;
		break;
		
		case 8:
		TwlkFuncType[0] = "[0]/sqrt(x-[1])+[2]";
		TwlkNPar[0]=3;
		TwlkInitialParSetFlag_Ref[0] = true;
		TwlkInitialFitPar[0][0] = -2.0;
		TwlkInitialFitPar[0][1] = 0.;
		TwlkInitialFitPar[0][2] = 1.0;
		break;

		default:
		TwlkFuncType[0] = "[0]/sqrt(x)+[1]";
		TwlkNPar[0]=2;
		break;
	}
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::SetTwlkRef2(int type){
	switch(type){
		case 0:
		TwlkFuncType[1] = "[0]/sqrt(x)+[1]";
		TwlkNPar[1]=2;
		break;
		
		case 1:
		TwlkFuncType[1] = "[0]*exp([1]*x)+[2]";
		TwlkNPar[1]=3;
		break;
		
		case 2:
		TwlkFuncType[1] = "[0]*log(x)+[1]";
		TwlkNPar[1]=2;
		break;

		case 3:
		TwlkFuncType[1] = "[0]*pow(x,[1])+[2]";
		TwlkNPar[1]=3;
		TwlkInitialParSetFlag_Ref[1] = true;
		TwlkInitialFitPar[1][0] = -2.0;
		TwlkInitialFitPar[1][1] = -0.5;
		TwlkInitialFitPar[1][2] = 0.5;
		break;
		
		case 4:
		TwlkFuncType[1] = "[0]/sqrt(x)+[1]*x+[2]";
		TwlkNPar[1]=3;
		break;
		
		case 5:
		TwlkFuncType[1] = "[0]*exp([1]*x)+[2]*x+[3]";
		TwlkNPar[1]=4;
		break;

		case 6:
		TwlkFuncType[1] = "[0]*log(x)+[1]*x+[2]";
		TwlkNPar[1]=3;
		break;

		case 7:
		TwlkFuncType[1] = "[0]*pow(x,[1])+[2]*x+[3]";
		TwlkNPar[1]=4;
		TwlkInitialParSetFlag_Ref[1] = true;
		TwlkInitialFitPar[1][0] = -2.0;
		TwlkInitialFitPar[1][1] = -0.5;
		TwlkInitialFitPar[1][2] = 0.5;
		TwlkInitialFitPar[1][3] = 0.1;
		break;
		
		case 8:
		TwlkFuncType[1] = "[0]/sqrt(x-[1])+[2]";
		TwlkNPar[1]=3;
		TwlkInitialParSetFlag_Ref[1] = true;
		TwlkInitialFitPar[1][0] = -2.0;
		TwlkInitialFitPar[1][1] = 5.;
		TwlkInitialFitPar[1][2] = 0.5;
		break;

		default:
		TwlkFuncType[1] = "[0]/sqrt(x)+[1]";
		TwlkNPar[0]=2;
		break;
	}
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::SetRoot(){
	cout<<Form("====     run%04d     =====", RunNum)<<endl;
	ifp = new TFile(RootFile.c_str(), "READONLY");
	tree = (TTree*)ifp->Get("tree");
	tree->SetBranchAddress( "tdc"   , tdc    );
	tree->SetBranchAddress( "qdc"   , qdc    );

	tree->SetBranchAddress( "evnum" , &EvNum );

	for(int i=0; i<NofDet; i++){tree->SetBranchAddress( Form("Ref%d",i+1), &Ref[i] );}

	Ent = tree->GetEntries();	cout<<Ent<<endl;
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::ExtractEvent(){
	cout<<"Extract effective event..."<<endl;
	for(int i=0; i<Ent; i++){
		tree->GetEntry(i);
		if(
			Trig_Ref1Ref2(tdc)
		){
			EffectiveEvent.push_back(i);
			EffectiveEvCounter++;
		}else;
	}
	MaxEventItr=EffectiveEvent.size();
	cout<<"Total effective event >> "<<MaxEventItr<<endl;
	
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::ImportPedestal(){
	ifs_pede.open(Param_Pede.c_str(), ios_base::in);
	if(!ifs_pede){
		cout<<"There is no file: "<<Param_Pede<<endl;
		return;
	}else;
	for(int i=0; i<NofMPPC; i++){
		ifs_pede>>pede_v[i]>>pede_w[i];	cout<<pede_v[i]<<"     "<<pede_w[i]<<endl;
	}
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::DefineCanv(){
	for(int i=0; i<NofCanv-4; i++){Ca[i] = new TCanvas( Form("Ca[%d]",i), Form("Ca[%d]",i), 1602, 824 );}
	Ca[NofCanv-4] = new TCanvas( Form("Ca[%d]", NofCanv-4), Form("Ca[%d]", NofCanv-4), 1602, 1624 );
	Ca[NofCanv-3] = new TCanvas( Form("Ca[%d]", NofCanv-3), Form("Ca[%d]", NofCanv-3), 1602, 1624 );
	Ca[NofCanv-2] = new TCanvas( Form("Ca[%d]", NofCanv-2), Form("Ca[%d]", NofCanv-2), 1602, 1624 );
	Ca[NofCanv-1] = new TCanvas( Form("Ca[%d]", NofCanv-1), Form("Ca[%d]", NofCanv-1), 1602, 1624 );
	Ca[0]->Divide(2,2);
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::DefineObj(){
	double TwlkQDCDiv = (TwlkQDCMax+10.)/TwlkQDCNBin;

	for(int i=0; i<NofDet; i++){
		//1-D histgram for QDC, TDC
		h_1d_qdc_full[i] = new TH1D( Form("h_1d_qdc_full[%d]",i), Form("h_1d_qdc_full[%d]",i), TwlkQDCNBin, -10., TwlkQDCMax );
		h_1d_qdc_wcut[i] = new TH1D( Form("h_1d_qdc_wcut[%d]",i), Form("h_1d_qdc_wcut[%d]",i), TwlkQDCNBin, -10., TwlkQDCMax );
		h_1d_tdc_full[i] = new TH1D( Form("h_1d_tdc_full[%d]",i), Form("h_1d_tdc_full[%d]",i), 4100, 100., 4200. );
		h_1d_tdc_wcut[i] = new TH1D( Form("h_1d_tdc_wcut[%d]",i), Form("h_1d_tdc_wcut[%d]",i), 4100, 100., 4200. );
		MySetting->Setting_Hist1D( h_1d_qdc_full[i], Form("run%04d: ", RunNum)+Label[CHassign[i]]+" QDC", Label[CHassign[i]]+" QDC (pC)"       , Form("Counts/%.2lf pC", TwlkQDCDiv), 2  , 1, 42, 616, 1001 );
		MySetting->Setting_Hist1D( h_1d_qdc_wcut[i], Form("run%04d: ", RunNum)+Label[CHassign[i]]+" QDC", Label[CHassign[i]]+" QDC (pC)"       , Form("Counts/%.2lf pC", TwlkQDCDiv), 602, 1, 42, 422, 1001 );
		MySetting->Setting_Hist1D( h_1d_tdc_full[i], Form("run%04d: ", RunNum)+Label[CHassign[i]]+" TDC", Label[CHassign[i]]+" TDC (ch./35 ps)", "Counts/ch."                       , 2  , 1, 42, 616, 1001 );
		MySetting->Setting_Hist1D( h_1d_tdc_wcut[i], Form("run%04d: ", RunNum)+Label[CHassign[i]]+" TDC", Label[CHassign[i]]+" TDC (ch./35 ps)", "Counts/ch."                       , 602, 1, 42, 422, 1001 );
		h_1d_qdc_full[i]->SetFillColorAlpha(616, 0.60);
		h_1d_qdc_wcut[i]->SetFillColorAlpha(422, 0.80);
		h_1d_tdc_full[i]->SetFillColorAlpha(616, 0.60);
		h_1d_tdc_wcut[i]->SetFillColorAlpha(422, 0.80);

		//2-D histgram for Correlation between QDC and TOF
		h_2d_rawtq[i]      = new TH2D( Form("h_2d_rawtq[%d]"     , i), Form("h_2d_rawtq[%d]"     , i), TwlkQDCNBin, -10., TwlkQDCMax, 200, -7., 7. );
		h_2d_rawtq_fit[i]  = new TH2D( Form("h_2d_rawtq_fit[%d]" , i), Form("h_2d_rawtq_fit[%d]" , i), TwlkQDCNBin, -10., TwlkQDCMax, 200, -7., 7. );
		h_2d_dectq_full[i] = new TH2D( Form("h_2d_dectq_full[%d]", i), Form("h_2d_dectq_full[%d]", i), TwlkQDCNBin, -10., TwlkQDCMax, 200, -7., 7. );
		h_2d_dectq_wcut[i] = new TH2D( Form("h_2d_dectq_wcut[%d]", i), Form("h_2d_dectq_wcut[%d]", i), TwlkQDCNBin, -10., TwlkQDCMax, 200, -7., 7. );
		MySetting->Setting_Hist2D( h_2d_rawtq[i]     , Form("#it{run%04d: ", RunNum)+Label[CHassign[i]]+" QDC vs. TOF_{Ref1-Ref2}}", Label[CHassign[i]]+" QDC (pC)", "RawTOF_{Ref1-Ref2} (ns)", "", 1. );
		MySetting->Setting_Hist2D( h_2d_rawtq_fit[i] , Form("#it{run%04d: ", RunNum)+Label[CHassign[i]]+" QDC vs. TOF_{Ref1-Ref2}}", Label[CHassign[i]]+" QDC (pC)", "RawTOF_{Ref1-Ref2} (ns)", "", 1. );
		MySetting->Setting_Hist2D( h_2d_dectq_full[i], Form("#it{run%04d: ", RunNum)+Label[CHassign[i]]+" QDC vs. TOF_{Ref1-Ref2}}", Label[CHassign[i]]+" QDC (pC)", "DecTOF_{Ref1-Ref2} (ns)", "", 1. );
		MySetting->Setting_Hist2D( h_2d_dectq_wcut[i], Form("#it{run%04d: ", RunNum)+Label[CHassign[i]]+" QDC vs. TOF_{Ref1-Ref2}}", Label[CHassign[i]]+" QDC (pC)", "DecTOF_{Ref1-Ref2} (ns)", "", 1. );

		f_twlk[i] = new TF1( Form("f_twlk[%d]", i), TwlkFuncType[i], 0., 500. );
		MySetting->Setting_Func(f_twlk[i], 6);
		int Npar = f_twlk[i]->GetNpar();
		for(int j=0; j<Npar; j++){f_twlk[i]->SetParName(j, Form("#it{p_{%d}}",j));}
		if(TwlkInitialParSetFlag){
			for(int j=0; j<TwlkNPar[i]; j++){
				f_twlk[i]->SetParameter(j, TwlkInitialFitPar[i][j]);
			}
			if(TwlkMode==7){
				f_twlk[i]->SetParLimits(1, -9., 0.);
			}else;
		}else;

		Leg_qdc[i] = new TLegend( .55, .65, .85, .95, "Trigger Type: Effective Ent./Total Ent.");
		MySetting->Setting_Legend( Leg_qdc[i], 42, 22, 602, 0.040 );
	//	Leg_qdc[i]->SetTextSize(0);
		Leg_qdc[i]->SetBorderSize(0);
		Leg_qdc[i]->SetLineWidth(0);
		Leg_qdc[i]->SetFillStyle(0);

		Leg_tdc[i] = new TLegend( .55, .65, .85, .95, "Trigger Type: Effective Ent./Total Ent.");
		MySetting->Setting_Legend( Leg_tdc[i], 42, 22, 602, 0.040 );
	//	Leg_tdc[i]->SetTextSize(0);
		Leg_tdc[i]->SetBorderSize(0);
		Leg_tdc[i]->SetLineWidth(0);
		Leg_tdc[i]->SetFillStyle(0);

		Pt[i] = new TPaveText();
		Pt[i] -> SetX1NDC(.675);
		Pt[i] -> SetX2NDC(.875);
		Pt[i] -> SetY1NDC(.200);
		Pt[i] -> SetY2NDC(.425);
		Pt[i] -> SetTextFont(42);
		Pt[i] -> SetTextSize(0.0325);
		Pt[i] -> SetFillColorAlpha(422, 0.65);
		Pt[i] -> SetLineColor(602);
		Pt[i] -> SetLineWidth(1);
		Pt[i] -> SetLineStyle(1);

		Ln_RefQDCCut[i] = new TLine();
		MySetting->Setting_Line( Ln_RefQDCCut[i], 602, 1, 4 );
	}
	h_1d_rawtof      = new TH1D( "h_1d_rawtof"     , "h_1d_rawtof"     , 400, -7., 7.);
	h_1d_dectof_full = new TH1D( "h_1d_dectof_full", "h_1d_dectof_full", 200, -3.5, 3.5);
	h_1d_dectof_wcut = new TH1D( "h_1d_dectof_wcut", "h_1d_dectof_wcut", 200, -3.5, 3.5);
		
	h_2d_twlkfact             = new TH2D( "h_2d_twlkfact"            , "h_2d_twlkfact"            , 101, -0.0125, 2.5125, 101, -0.0125, 2.5125 );
	h_2d_twlkfact_Chis        = new TH2D( "h_2d_twlkfact_Chis"       , "h_2d_twlkfact_Chis"       , 101, -0.0125, 2.5125, 101, -0.0125, 2.5125 );
	h_2d_twlkfact_RMS         = new TH2D( "h_2d_twlkfact_RMS"        , "h_2d_twlkfact_RMS"        , 101, -0.0125, 2.5125, 101, -0.0125, 2.5125 );
	h_2d_twlkfact_Skew        = new TH2D( "h_2d_twlkfact_Skew"       , "h_2d_twlkfact_Skew"       , 101, -0.0125, 2.5125, 101, -0.0125, 2.5125 );
	h_2d_twlkfact_ChisConsFit = new TH2D( "h_2d_twlkfact_ChisConsFit", "h_2d_twlkfact_ChisConsFit", 101, -0.0125, 2.5125, 101, -0.0125, 2.5125 );

	MySetting->Setting_Hist1D( h_1d_rawtof     , Form("run%04d: RawTOF"            ,RunNum), "TOF_{Ref1-Ref2} (ns)", "Counts/0.035 ns", 843, 1, 42, 406, 1001 );
	MySetting->Setting_Hist1D( h_1d_dectof_full, Form("run%04d: DecTOF w/o QDC cut",RunNum), "TOF_{Ref1-Ref2} (ns)", "Counts/0.035 ns", 4  , 1, 42, 422, 1001 );
	MySetting->Setting_Hist1D( h_1d_dectof_wcut, Form("run%04d: DecTOF w/ QDC cut" ,RunNum), "TOF_{Ref1-Ref2} (ns)", "Counts/0.035 ns", 4  , 1, 42, 422, 1001 );
	h_1d_rawtof      -> GetXaxis() -> SetNdivisions(515);
	h_1d_dectof_full -> GetXaxis() -> SetNdivisions(515);
	h_1d_dectof_wcut -> GetXaxis() -> SetNdivisions(515);


	//=== h_2d_twlkfact ====
	MySetting->Setting_Hist2D( h_2d_twlkfact, "Factor Ref1 vs. Ref2 (#sigma_{Ref1-Ref2})", "Ref1", "Ref2", "#sigma_{Ref1-Ref2} (ns)", 0. );
	//Other Setting
	////X axis
	h_2d_twlkfact->GetXaxis()->SetLabelSize(0.025);
	////Y axis
	h_2d_twlkfact->GetYaxis()->SetTitleOffset( 1.30*h_2d_twlkfact->GetYaxis()->GetTitleOffset() );
	h_2d_twlkfact->GetYaxis()->SetLabelSize(0.025);
	////Z axis
	h_2d_twlkfact->GetZaxis()->SetLabelSize(0.020);
	h_2d_twlkfact->GetZaxis()->SetTitleSize(0.030);

	//=== h_2d_twlkfact_Chis ====
	MySetting->Setting_Hist2D( h_2d_twlkfact_Chis, "Factor Ref1 vs. Ref2 (#chi^{2}/NDF)", "Ref1", "Ref2", "#chi^{2}/NDF", 0. );
	//Other Setting
	////X axis
	h_2d_twlkfact_Chis->GetXaxis()->SetLabelSize(0.025);
	////Y axis
	h_2d_twlkfact_Chis->GetYaxis()->SetTitleOffset( 1.30*h_2d_twlkfact_Chis->GetYaxis()->GetTitleOffset() );
	h_2d_twlkfact_Chis->GetYaxis()->SetLabelSize(0.025);
	////Z axis
	h_2d_twlkfact_Chis->GetZaxis()->SetLabelSize(0.020);
	h_2d_twlkfact_Chis->GetZaxis()->SetTitleSize(0.030);

	//=== h_2d_twlkfact_RMS ====
	MySetting->Setting_Hist2D( h_2d_twlkfact_RMS, "Factor Ref1 vs. Ref2 (RMS)", "Ref1", "Ref2", "RMS_{Ref1-Ref2} (ns)", 0. );
	//Other Setting
	////X axis
	h_2d_twlkfact_RMS->GetXaxis()->SetLabelSize(0.025);
	////Y axis
	h_2d_twlkfact_RMS->GetYaxis()->SetTitleOffset( 1.30*h_2d_twlkfact_RMS->GetYaxis()->GetTitleOffset() );
	h_2d_twlkfact_RMS->GetYaxis()->SetLabelSize(0.025);
	////Z axis
	h_2d_twlkfact_RMS->GetZaxis()->SetLabelSize(0.020);
	h_2d_twlkfact_RMS->GetZaxis()->SetTitleSize(0.030);

	//=== h_2d_twlkfact_Skew ====
//	MySetting->Setting_Hist2D( h_2d_twlkfact_Skew, "Factor Ref1 vs. Ref2 (Skewness)", "Ref1", "Ref2", "Skewness", 0. );
	//Title
	h_2d_twlkfact_Skew->SetTitle("Factor Ref1 vs. Ref2 (Skewness)");
	h_2d_twlkfact_Skew->GetXaxis()->SetTitle("Ref1");
	h_2d_twlkfact_Skew->GetYaxis()->SetTitle("Ref2");
	h_2d_twlkfact_Skew->GetZaxis()->SetTitle("Skewness");
	//Other Setting
	////X axis
	h_2d_twlkfact_Skew->GetXaxis()->SetTitleFont(62);
	h_2d_twlkfact_Skew->GetXaxis()->SetTitleSize(0.05);
	h_2d_twlkfact_Skew->GetXaxis()->SetTitleOffset(0.80);
	h_2d_twlkfact_Skew->GetXaxis()->SetLabelFont(42);
	h_2d_twlkfact_Skew->GetXaxis()->SetLabelSize(0.025);
	h_2d_twlkfact_Skew->GetXaxis()->SetLabelOffset(0.01);
	h_2d_twlkfact_Skew->GetXaxis()->SetNdivisions(515);
	////Y axis
	h_2d_twlkfact_Skew->GetYaxis()->SetTitleFont(62);
	h_2d_twlkfact_Skew->GetYaxis()->SetTitleSize(0.05);
	h_2d_twlkfact_Skew->GetYaxis()->SetTitleOffset( 1.30*h_2d_twlkfact_Skew->GetYaxis()->GetTitleOffset() );
	h_2d_twlkfact_Skew->GetYaxis()->SetLabelFont(42);
	h_2d_twlkfact_Skew->GetYaxis()->SetLabelSize(0.025);
	h_2d_twlkfact_Skew->GetYaxis()->SetLabelOffset(0.01);
	h_2d_twlkfact_Skew->GetYaxis()->SetNdivisions(515);
	////Z axis
	h_2d_twlkfact_Skew->GetZaxis()->SetLabelSize(0.020);
	h_2d_twlkfact_Skew->GetZaxis()->SetLabelOffset(0.008);
	h_2d_twlkfact_Skew->GetZaxis()->SetTitleSize(0.030);
	h_2d_twlkfact_Skew->GetZaxis()->SetTitleFont(62);
	h_2d_twlkfact_Skew->GetZaxis()->SetLabelFont(42);

	MySetting->Setting_Hist2D( h_2d_twlkfact_ChisConsFit, "Factor Ref1 vs. Ref2 (#chi^{2}/NDF for Conts. fit)", "Ref1", "Ref2", "ln(#chi^{2}/NDF_{Const. fit})", 0. );
	h_2d_twlkfact_ChisConsFit->GetXaxis()->SetLabelSize(0.025);
	h_2d_twlkfact_ChisConsFit->GetYaxis()->SetLabelSize(0.025);
	h_2d_twlkfact_ChisConsFit->GetZaxis()->SetLabelSize(0.020);
	h_2d_twlkfact_ChisConsFit->GetZaxis()->SetTitleSize(0.030);
	h_2d_twlkfact_ChisConsFit->GetYaxis()->SetTitleOffset( 1.20*h_2d_twlkfact_Chis->GetYaxis()->GetTitleOffset() );


	gr_twlkfact = new TGraph();
	MySetting->Setting_Graph( gr_twlkfact,"Factor Ref1 vs. Ref2", "Ref1", "Ref2",6, 1, 42, 6, 29 );

	for(int i=0; i<3; i++){
		f_tof[i] = new TF1( Form("f_tof[%d]",i), "gaus(0)", -10., 10. );
		MySetting->Setting_Func(f_tof[i], 6);
		f_tof[i]->SetNpx(5E+4);
	}
	f_tof_tmp = new TF1( "f_tof_tmp", "gaus(0)", -1.E+2, 1.E+2);
	for(int i=0; i<NofDet; i++){f_const_tmp[i] = new TF1( Form("f_const_tmp[%d]",i), "[0]", -10., 400.);}

	Lat = new TLatex();
	MySetting->Setting_Latex( Lat, 42, 22, 602, 0.05 );
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::FillHist(){
	for(int i=0; i<Ent; i++){
		tree->GetEntry(i);
		h_1d_qdc_full[0]->Fill(Ref[0].e);
		h_1d_qdc_full[1]->Fill(Ref[1].e);
		h_1d_tdc_full[0]->Fill(Ref[0].tdc);
		h_1d_tdc_full[1]->Fill(Ref[1].tdc);
	}


	for(int i=0; i<MaxEventItr; i++){
		tree->GetEntry(EffectiveEvent[i]);
		rtof = Ref[0].t-Ref[1].t;
		if(
			   Ref[0].qdc<4095
			&& Ref[1].qdc<4095
		){
			h_1d_qdc_wcut[0]  -> Fill(Ref[0].e);
			h_1d_qdc_wcut[1]  -> Fill(Ref[1].e);
			h_1d_tdc_wcut[0]  -> Fill(Ref[0].tdc);
			h_1d_tdc_wcut[1]  -> Fill(Ref[1].tdc);
			h_1d_rawtof       -> Fill(rtof);
			h_2d_rawtq[0]     -> Fill(Ref[0].e, rtof);
			h_2d_rawtq[1]     -> Fill(Ref[1].e, rtof);
			h_2d_rawtq_fit[0] -> Fill(Ref[0].e, rtof);
			h_2d_rawtq_fit[1] -> Fill(Ref[1].e, rtof);
		}else;
	}
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::Draw(){
	for(int i=0; i<2; i++){
		Ca[0]->cd(i+1);
		gPad->SetLogy(1);
		h_1d_qdc_full[i]->SetStats(kFALSE);
		h_1d_qdc_wcut[i]->SetStats(kFALSE);
		Leg_qdc[i]->AddEntry( h_1d_qdc_full[i], Form("Full Event: %.0lf/%.0lf"        , h_1d_qdc_full[i]->GetEffectiveEntries(), h_1d_qdc_full[i]->GetEntries()), "fl");
		Leg_qdc[i]->AddEntry( h_1d_qdc_wcut[i], Form("Ref1 #otimes Ref2: %.0lf/%.0lf", h_1d_qdc_wcut[i]->GetEffectiveEntries(), h_1d_qdc_wcut[i]->GetEntries()), "fl");
		h_1d_qdc_full[i]->Draw("");
		h_1d_qdc_wcut[i]->Draw("same");
		gPad->Update();
		gPad->Modified();
		Leg_qdc[i]->Draw();
	
		Ca[0]->cd(i+3);
		gPad->SetLogy(1);
		h_1d_tdc_full[i]->SetStats(kFALSE);
		h_1d_tdc_wcut[i]->SetStats(kFALSE);
		Leg_tdc[i]->AddEntry( h_1d_qdc_full[i], Form("Full Event: %.0lf/%.0lf"        , h_1d_tdc_full[i]->GetEffectiveEntries(), h_1d_tdc_full[i]->GetEntries()), "fl");
		Leg_tdc[i]->AddEntry( h_1d_qdc_wcut[i], Form("Ref1 #otimes Ref2: %.0lf/%.0lf", h_1d_tdc_wcut[i]->GetEffectiveEntries(), h_1d_tdc_wcut[i]->GetEntries()), "fl");
		h_1d_tdc_full[i]->Draw("");
		h_1d_tdc_wcut[i]->Draw("same");
		gPad->Update();
		gPad->Modified();
		Leg_tdc[i]->Draw();
	}
	
	for(int i=0; i<NofDet; i++){
		Ca[i+1]->cd();
		gPad->SetLogz(1);
		h_2d_rawtq[i]->Draw("colz");
		gPad->Update();
		gPad->Modified();
		fr_2d_rawtq[i] = gPad->GetFrame();
		fr_2d_rawtq[i] -> SetFillStyle(0);
		fr_2d_rawtq[i] -> Draw();
	}

	Ca[3]->cd();
	gPad->SetRightMargin(0.075);
	h_1d_rawtof->Draw();
	Par_tmp = h_1d_rawtof->GetBinCenter( h_1d_rawtof->GetMaximumBin() );
	FitWid_Tmp=2.;
	for(int i=0; i<10; i++){
		h_1d_rawtof->Fit(f_tof[0], "", "", Par_tmp- FitWid_Tmp, Par_tmp+FitWid_Tmp);
		Par_tmp    = f_tof[0]->GetParameter(1);
		FitWid_Tmp = f_tof[0]->GetParameter(2);
		FitWid_Tmp *=5.;
	}
	Lat->DrawLatexNDC(.75, .50, Form("#sigma_{Ref1-Ref2} = %.1lf #pm %.1lf ps", 1000.*(f_tof[0]->GetParameter(2)), 1000.*(f_tof[0]->GetParError(2))));
	Lat->DrawLatexNDC(.75, .40, Form("RMS_{Ref1-Ref2} = %.1lf ps", 1000.*(h_1d_rawtof->GetRMS()) ) );

}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::Fit(){
	for(int i=0; i<NofDet; i++){
		Ca[i+4]->cd();
		Ca[i+4]->SetLogz(1);
		h_2d_rawtq_fit[i]->Draw("colz");
		
		TF1* f_FitSlices = new TF1("f_FitSlices", "gaus(0)", -7., 7.);
		h_2d_rawtq_fit[i]->FitSlicesY( f_FitSlices, 0, -1, 0, TwlkFitSliceOpt[i] );
		h_1d_rawtq_slice[i] = (TH1D*)gROOT->FindObject( Form("h_2d_rawtq_fit[%d]_1", i) );
		Ca[i+4]->cd();
		h_1d_rawtq_slice[i] -> SetLineColor(1);
		h_1d_rawtq_slice[i] -> Draw("samePES");
	
		for(int j=0; j<100; j++){h_1d_rawtq_slice[i]->Fit( f_twlk[i], "", "", TwlkFitMin[i], TwlkFitMax[i] );}
		int NPar = f_twlk[i]->GetNpar();
		for(int j=0; j<NPar; j++){Twlk[i][j]=f_twlk[i]->GetParameter(j);}
		gPad->Update();
		gPad->Modified();
		fr_2d_rawtq_fit[i] = gPad->GetFrame();
		fr_2d_rawtq_fit[i] -> SetFillStyle(0);
		fr_2d_rawtq_fit[i] -> Draw();
		gPad->Update();
		Pt[i] -> AddText( Form("#chi^{2}/NDF= %.1lf/%d", f_twlk[i]->GetChisquare(), f_twlk[i]->GetNDF()) );
		for(int j=0; j<NPar; j++){
			double fitVal = f_twlk[i]->GetParameter(j);
			double fitErr = f_twlk[i]->GetParError(j);
			Pt[i]->AddText( Form("#it{p}_{%d} = %.4lf #pm %.4lf", j, fitVal, fitErr));
		}
		Pt[i]->Draw();
	}

	for(int i=0; i<NofDet; i++){
		Ca[i+6]->cd();
		h_1d_rawtq_slice[i]->GetYaxis()->SetTitle("RawTOF_{Ref1-Ref2} (ns)");
		h_1d_rawtq_slice[i]->Draw();
	}
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::Fit_Ref(){
}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::SearchBest(){
	for(int i=0; i<101;i++){
		for(int j=0; j<101;j++){
			TwlkFactor[0] = 0. + TwlkDelta*i;
			TwlkFactor[1] = 0. + TwlkDelta*j;

			h_1d_dectof_tmp = new TH1D( "h_1d_dectof_tmp", "h_1d_dectof_tmp", 400, -7., 7. );
			
			for(int k=0; k<MaxEventItr; k++){
				tree->GetEntry(EffectiveEvent[k]);
				decTime[0] = Ref[0].t-Twlk_Correction( Ref[0].e, Twlk[0], TwlkFactor[0], Ref1TwlkMode );
				decTime[1] = Ref[1].t+Twlk_Correction( Ref[1].e, Twlk[1], TwlkFactor[1], Ref2TwlkMode );
				dtof = decTime[0]-decTime[1];
				h_1d_dectof_tmp->Fill(dtof);
				for(int seg=0; seg<NofDet; seg++){
					h_2d_dectq_tmp[seg]->Fill(Ref[seg].e, dtof);
				}
			}

			Ca_dummy = new TCanvas( "Ca_dummy", "Ca_dummy", 482, 264 );
			Par_tmp = h_1d_dectof_tmp->GetBinCenter( h_1d_dectof_tmp->GetMaximumBin() );
			for(int Itr=0; Itr<FitItrTmp; Itr++){
				h_1d_dectof_tmp->Fit( f_tof_tmp, "Q", "", Par_tmp-FitWid_Tmp, Par_tmp+FitWid_Tmp);
				Par_tmp    = f_tof_tmp->GetParameter(1);
				FitWid_Tmp = f_tof_tmp->GetParameter(2);
				FitWid_Tmp *= 5.;							//fit with 5sigma
			}
	
			//Sigma
			TwlkSigVal_Tmp = f_tof_tmp->GetParameter(2);
			TwlkSigErr_Tmp = f_tof_tmp->GetParError(2);

			//Chisquare/NDF
			TwlkChisNDF_Tmp = (f_tof_tmp->GetChisquare())/(f_tof_tmp->GetNDF());

			//RMS
			TwlkRMS_Tmp = h_1d_dectof_tmp->GetRMS();

			//Skewness
			TwlkSkew_Tmp = h_1d_dectof_tmp->GetSkewness();

			h_1d_dectof_tmp -> Delete();

			Ca_dummy -> Destructor();

			if( TwlkSigVal_Tmp<TwlkSigVal_Best ){
				TwlkSigVal_Best = TwlkSigVal_Tmp;
				TwlkSigErr_Best = TwlkSigErr_Tmp;
				for(int seg=0; seg<NofDet; seg++){TwlkFact_Best[seg] = TwlkFactor[seg];}
				cout<<"Factor Ref1 & Ref2 >> "<<TwlkFactor[0]<<", "<<TwlkFactor[1]<<": "<<"sigma_TOF = "<<TwlkSigVal_Best<<" ns"<<endl;
			}else;
			h_2d_twlkfact->Fill( TwlkFactor[0], TwlkFactor[1], TwlkSigVal_Tmp );

			if( TwlkChisNDF_Tmp<TwlkChisNDF_Best ){
				TwlkChisNDF_Best = TwlkChisNDF_Tmp;
			}else;
			h_2d_twlkfact_Chis->Fill( TwlkFactor[0], TwlkFactor[1], TwlkChisNDF_Tmp );

			if( TwlkRMS_Tmp<TwlkRMS_Best ){
				TwlkRMS_Best = TwlkRMS_Tmp;
			//	for(int seg=0; seg<NofDet; seg++){TwlkFact_Best[seg] = TwlkFactor[seg];}
			//	cout<<"Factor Ref1 & Ref2 >> "<<TwlkFactor[0]<<", "<<TwlkFactor[1]<<": "<<"RMS_TOF = "<<TwlkRMS_Best<<" ns"<<endl;
			}else;
			h_2d_twlkfact_RMS->Fill( TwlkFactor[0], TwlkFactor[1], TwlkRMS_Tmp );

			if( abs(TwlkSkew_Tmp)<abs(TwlkSkew_Best) ){
				TwlkSkew_Best = TwlkSkew_Tmp;
			}else;
			h_2d_twlkfact_Skew->Fill( TwlkFactor[0], TwlkFactor[1], TwlkSkew_Tmp );


		}
	}
	h_2d_twlkfact             -> SetMinimum(TwlkSigVal_Best);
	h_2d_twlkfact_Chis        -> SetMinimum(TwlkChisNDF_Best);
	h_2d_twlkfact_RMS         -> SetMinimum(TwlkRMS_Best);
	double MaximumTmp=0.;
	MaximumTmp=h_2d_twlkfact->GetMaximumBin();
	if(MaximumTmp>0.30){
		h_2d_twlkfact->SetMaximum(0.30);
	}else;
	MaximumTmp=h_2d_twlkfact_RMS->GetMaximumBin();
	if(MaximumTmp>0.8){
		h_2d_twlkfact_RMS->SetMaximum(0.8);
	}else;
	gr_twlkfact->SetPoint( 0, TwlkFact_Best[0], TwlkFact_Best[1]);
}
/*
*/

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::DecTOF(){

	//Fill
	for(int i=0; i<MaxEventItr; i++){
		tree->GetEntry( EffectiveEvent[i] );
		decTime[0] = Ref[0].t-Twlk_Correction( Ref[0].e, Twlk[0], TwlkFact_Best[0], Ref1TwlkMode );
		decTime[1] = Ref[1].t+Twlk_Correction( Ref[1].e, Twlk[1], TwlkFact_Best[1], Ref2TwlkMode );
		dtof = decTime[0]-decTime[1];
		h_2d_dectq_full[0] -> Fill(Ref[0].e, dtof);
		h_2d_dectq_full[1] -> Fill(Ref[1].e, dtof);
		h_1d_dectof_full   -> Fill(dtof);
		if(
			   Ref[0].e>TwlkRefQDCCutLim[0]
			&& Ref[1].e>TwlkRefQDCCutLim[1]
		){
			h_2d_dectq_wcut[0] -> Fill(Ref[0].e, dtof);
			h_2d_dectq_wcut[1] -> Fill(Ref[1].e, dtof);
			h_1d_dectof_wcut   -> Fill(dtof);
		}else;
	}

	//Draw 2D w/o QDC cut (Full event)
	for(int i=0; i<NofDet; i++){
		Ca[i+8]->cd();
		Ca[i+8]->SetLogz(1);
		h_2d_dectq_full[i]->Draw("colz");
		gPad->Update();
		gPad->Modified();
		Ln_RefQDCCut[i]->SetX1(TwlkRefQDCCutLim[i]);
		Ln_RefQDCCut[i]->SetX2(TwlkRefQDCCutLim[i]);
		Ln_RefQDCCut[i]->SetY1(gPad->GetUymin());
		Ln_RefQDCCut[i]->SetY2(gPad->GetUymax());
		Ln_RefQDCCut[i]->Draw();
		gPad->Update();
		gPad->Modified();
		fr_2d_dectq_full[i]=gPad->GetFrame();
		fr_2d_dectq_full[i]->SetFillStyle(0);
		fr_2d_dectq_full[i]->Draw();
	}

	Ca[10]->cd();
	gPad->SetRightMargin(.075);
	h_1d_dectof_full->Draw("");
	Par_tmp = h_1d_dectof_full->GetBinCenter( h_1d_dectof_full->GetMaximumBin() );
	FitWid_Tmp=2.;
	for(int i=0; i<10; i++){
		h_1d_dectof_full->Fit(f_tof[1], "", "", Par_tmp-FitWid_Tmp, Par_tmp+FitWid_Tmp);
		Par_tmp    = f_tof[1]->GetParameter(1);
		FitWid_Tmp = f_tof[1]->GetParameter(2);
		FitWid_Tmp *=5.;
	}
	Lat->DrawLatexNDC(.75, .50, Form("#sigma_{Ref1-Ref2} = %.1lf #pm %.1lf ps", 1000.*(f_tof[1]->GetParameter(2)), 1000.*(f_tof[1]->GetParError(2))));
	Lat->DrawLatexNDC(.75, .40, Form("RMS_{Ref1-Ref2} = %.1lf ps", 1000.*(h_1d_dectof_full->GetRMS()) ) );
	
	//
	for(int i=0; i<NofDet; i++){
		Ca[i+11]->cd();
		Ca[i+11]->SetLogz(1);
		h_2d_dectq_wcut[i]->Draw("colz");
		gPad->Update();
		gPad->Modified();
		Ln_RefQDCCut[i]->Draw();
		gPad->Update();
		gPad->Modified();
		fr_2d_dectq_wcut[i]=gPad->GetFrame();
		fr_2d_dectq_wcut[i]->SetFillStyle(0);
		fr_2d_dectq_wcut[i]->Draw();
	}

	Ca[13]->cd();
	gPad->SetRightMargin(.075);
	h_1d_dectof_wcut->Draw("");
	Par_tmp = h_1d_dectof_wcut->GetBinCenter( h_1d_dectof_wcut->GetMaximumBin() );
	FitWid_Tmp=2.;
	for(int i=0; i<10; i++){
		h_1d_dectof_wcut->Fit(f_tof[2], "", "", Par_tmp-FitWid_Tmp, Par_tmp+FitWid_Tmp);
		Par_tmp    = f_tof[2]->GetParameter(1);
		FitWid_Tmp = f_tof[2]->GetParameter(2);
		FitWid_Tmp *=5.;
	}
	Lat->DrawLatexNDC(.75, .50, Form("#sigma_{Ref1-Ref2} = %.1lf #pm %.1lf ps", 1000.*(f_tof[2]->GetParameter(2)), 1000.*(f_tof[2]->GetParError(2))));
	Lat->DrawLatexNDC(.75, .40, Form("RMS_{Ref1-Ref2} = %.1lf ps", 1000.*(h_1d_dectof_wcut->GetRMS()) ) );

	Ca[14]->cd();
	Ca[14]->SetGrid(0,0);
	gPad->SetRightMargin(.135);
	gPad->SetLeftMargin(.125);
	gPad->SetTopMargin(.105);
	gPad->SetBottomMargin(.155);
	h_2d_twlkfact->SetStats(kFALSE);
	h_2d_twlkfact->Draw("colz");
	gr_twlkfact->SetMarkerSize( 2.50*gr_twlkfact->GetMarkerSize() );
	gr_twlkfact->Draw("sameP");
	Lat->SetTextColor(6);
	Lat->DrawLatex( TwlkFact_Best[0], TwlkFact_Best[1]+0.20, Form("(%.3lf, %.3lf)", TwlkFact_Best[0], TwlkFact_Best[1]));

	Ca[15]->cd();
	Ca[15]->SetGrid(0,0);
	gPad->SetRightMargin(.135);
	gPad->SetLeftMargin(.125);
	gPad->SetTopMargin(.105);
	gPad->SetBottomMargin(.155);
	h_2d_twlkfact_Chis->SetStats(kFALSE);
	h_2d_twlkfact_Chis->Draw("colz");
	gr_twlkfact->Draw("sameP");
	Lat->SetTextColor(6);
	Lat->DrawLatex( TwlkFact_Best[0], TwlkFact_Best[1]+0.20, Form("(%.3lf, %.3lf)", TwlkFact_Best[0], TwlkFact_Best[1]));

	Ca[16]->cd();
	Ca[16]->SetGrid(0,0);
	gPad->SetRightMargin(.135);
	gPad->SetLeftMargin(.125);
	gPad->SetTopMargin(.105);
	gPad->SetBottomMargin(.155);
	h_2d_twlkfact_RMS->SetStats(kFALSE);
	h_2d_twlkfact_RMS->Draw("colz");
	gr_twlkfact->Draw("sameP");
	Lat->SetTextColor(6);
	Lat->DrawLatex( TwlkFact_Best[0], TwlkFact_Best[1]+0.20, Form("(%.3lf, %.3lf)", TwlkFact_Best[0], TwlkFact_Best[1]));

	Ca[17]->cd();
	Ca[17]->SetGrid(0,0);
	gPad->SetRightMargin(.135);
	gPad->SetLeftMargin(.125);
	gPad->SetTopMargin(.100);
	gPad->SetBottomMargin(.155);
	h_2d_twlkfact_Skew->SetStats(kFALSE);
	h_2d_twlkfact_Skew->Draw("colz");
	gr_twlkfact->Draw("sameP");
	Lat->SetTextColor(6);
	Lat->DrawLatex( TwlkFact_Best[0], TwlkFact_Best[1]+0.20, Form("(%.3lf, %.3lf)", TwlkFact_Best[0], TwlkFact_Best[1]));

}

//--------------------------------------------------------------------------------------------------
void Twlk_Ana::Export(){
	Ca[0]->Print(Figure_Open.c_str(), "pdf");
	for(int i=0; i<NofCanv; i++){Ca[i]->Print(Figure.c_str());}
	Ca[NofCanv-1]->Print(Figure_Close.c_str(), "pdf");

	ofp_obj = new TFile( RootFile_Ex_obj.c_str(), "recreate" );
	ofp_obj->cd();
	for(int i=0; i<NofCanv; i++){
		Ca[i]->Write( Form("Ca_%02d", i));
	}
	for(int i=0; i<NofDet; i++)h_1d_qdc_full[i]    -> Write( Form("h_1d_qdc_full_%02d"   , i));
	for(int i=0; i<NofDet; i++)h_1d_qdc_wcut[i]    -> Write( Form("h_1d_qdc_wcut_%02d"   , i));
	for(int i=0; i<NofDet; i++)h_1d_tdc_full[i]    -> Write( Form("h_1d_tdc_full_%02d"   , i));
	for(int i=0; i<NofDet; i++)h_1d_tdc_wcut[i]    -> Write( Form("h_1d_tdc_wcut_%02d"   , i));
	for(int i=0; i<NofDet; i++)h_2d_rawtq[i]       -> Write( Form("h_2d_rawtq_%02d"      , i));
	for(int i=0; i<NofDet; i++)h_2d_rawtq_fit[i]   -> Write( Form("h_2d_rawtq_fit_%02d"  , i));
	for(int i=0; i<NofDet; i++)h_2d_dectq_full[i]  -> Write( Form("h_2d_dectq_full_%02d" , i));
	for(int i=0; i<NofDet; i++)h_2d_dectq_wcut[i]  -> Write( Form("h_2d_dectq_wcut_%02d" , i));
	for(int i=0; i<NofDet; i++)f_twlk[i]           -> Write( Form("f_twlk_%02d"          , i));
	for(int i=0; i<NofDet; i++)h_1d_rawtq_slice[i] -> Write( Form("h_1d_rawtq_slice_%02d", i));
	h_1d_rawtof      -> Write("h_1d_rawtof");
	h_1d_dectof_full -> Write("h_1d_dectof_full");
	h_1d_dectof_wcut -> Write("h_1d_dectof_wcut");

	h_2d_twlkfact             -> Write("h_2d_twlkfact");
	h_2d_twlkfact_Chis        -> Write("h_2d_twlkfact_Chis");
	h_2d_twlkfact_RMS         -> Write("h_2d_twlkfact_RMS");
	h_2d_twlkfact_Skew        -> Write("h_2d_twlkfact_Skew");
	gr_twlkfact               -> Write("gr_twlkfact");

	ofp_obj->Close();

	ofp_dat = new TFile( RootFile_Ex_dat.c_str(), "recreate" );
	tree_Ex = new TTree( "tree", "tree" );

	tree_Ex -> Branch( "EvNum", &EvNum_Ex, "EventNumber/I"             );
	tree_Ex -> Branch( "tdc"  , tdc_Ex   , "V775_1[32]/I"              );
	tree_Ex -> Branch( "qdc"  , qdc_Ex   , "V792_1[32]/I"              );
	tree_Ex -> Branch( "Ref1" , &Ex_Ref1 , "TDC/I:QDC/I:Time/D:Edep/D" );
	tree_Ex -> Branch( "Ref2" , &Ex_Ref2 , "TDC/I:QDC/I:Time/D:Edep/D" );
	
	for(int i=0; i<Ent; i++){
		ofp_dat->cd();
		tree_Ex->GetEntry(i);
		
		ifp->cd();
		tree->GetEntry(i);

		////Ref1
		Ex_Ref1.tdc = Ref[0].tdc;
		Ex_Ref1.qdc = Ref[0].qdc;
		Ex_Ref1.t   = Ref[0].t-Twlk_Correction( Ref[0].e, Twlk[0], TwlkFactor[0], Ref1TwlkMode );
		Ex_Ref1.e   = Ref[0].e;
		////Ref2
		Ex_Ref2.tdc = Ref[1].tdc;
		Ex_Ref2.qdc = Ref[1].qdc;
		Ex_Ref2.t   = Ref[1].t+Twlk_Correction( Ref[1].e, Twlk[1], TwlkFactor[1], Ref2TwlkMode );
		Ex_Ref2.e   = Ref[1].e;
		for(int seg=0; seg<32; seg++){
			tdc_Ex[seg]=tdc[seg];
			qdc_Ex[seg]=qdc[seg];
		}
		EvNum_Ex=i;

		ofp_dat->cd();
		tree_Ex->Fill();
	}
	tree_Ex->Write();
	ofp_dat->Close();

	ofs_twlk.open( Param_Twlk.c_str(), ios_base::trunc );
	if(!ofs_twlk){
		return;
	}else;
	for(int i=0; i<NofDet; i++){
		ofs_twlk<<CHassign[i]<<"   ";
		ofs_twlk<<TwlkMode<<" ";
		for(int j=0; j<TwlkNPar[i]; j++){ofs_twlk<<Twlk[i][j]<<" ";}
		ofs_twlk<<TwlkFact_Best[i]<<endl;
	}
}

//--------------------------------------------------------------------------------------------------
/*
bool Twlk_Ana::Trig_Ref1FToF1(int* t){
h

	bool ret=true;

	if(
		   t[CHassign[0]]>0
		&& t[CHassign[2]]>0
		&& t[CHassign[3]]>0
	){
		ret=true;
	}else{
		ret=false;
	}

	return ret;
}

//--------------------------------------------------------------------------------------------------
bool Twlk_Ana::Trig_FToF1Ref2(int* t){
	bool ret=true;

	if(
		   t[CHassign[1]]>0
		&& t[CHassign[2]]>0
		&& t[CHassign[3]]>0
	){
		ret=true;
	}else{
		ret=false;
	}

	return ret;
}
*/
//--------------------------------------------------------------------------------------------------
bool Twlk_Ana::Trig_Ref1Ref2(int* t){
	bool ret=true;

	if(
		   t[CHassign[0]]>0
		&& t[CHassign[1]]>0
	){
		ret=true;
	}else{
		ret=false;
	}

	return ret;
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
int main(int argc, char** argv){
	int runnum=1068;
	int Mode_1=0;
	int Mode_2=0;

	if(argc==2){
		runnum=atoi(argv[1]);
	}else if(argc==3){
		runnum = atoi(argv[1]);
		Mode_1 = atoi(argv[2]);
	}else if(argc==4){
		runnum = atoi(argv[1]);
		Mode_1 = atoi(argv[2]);
		Mode_2 = atoi(argv[3]);
	}else{
		runnum = 1068;
		Mode_1 = 8;
		Mode_2 = 8;
	}


	TApplication* theApp;
	Twlk_Ana* ana;

	theApp = new TApplication( "App", &argc, argv );
//	ana = new Twlk_Ana(1, 2);
	ana = new Twlk_Ana(runnum, Mode_1, Mode_2);
	ana->TwlkMan();
	ana->SetRoot();
	ana->ExtractEvent();
	ana->DefineCanv();
	ana->DefineObj();
	ana->FillHist();
	ana->Draw();
	ana->Fit();
	ana->SearchBest();
	ana->DecTOF();
	ana->Export();

	delete ana;
	gSystem->Exit(-1);
	theApp->Run();

	return 0;
}
