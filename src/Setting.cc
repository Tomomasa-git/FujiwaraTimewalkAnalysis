//Setting.cc
////2020. 07. 08 (Wed)
////Tomomasa FUJIWARA
////	Reviced 2020. 11. 30 (Mon)

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
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./Setting.h"

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
Setting::Setting(){
	cout<<"Setting::Constructor Called"<<endl;
}
/////////////////////////////////////////////////////////////////////////////////////
Setting::~Setting(){
	cout<<"Setting::Deconstructor Called"<<endl;
}
/////////////////////////////////////////////////////////////////////////////////////
void Setting::Setting_Gene( int BatchFlag ){
	gROOT -> Reset();
	gROOT -> SetBatch(BatchFlag);

	//Pad and Grid
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetGridColor(kGray);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
//	gStyle->SetPadRightMargin(.05);
//	gStyle->SetPadLeftMargin(.135);
//	gStyle->SetPadTopMargin(.05);
	gStyle->SetPadTopMargin(TopMarg);
	gStyle->SetPadBottomMargin(BotMarg);
	gStyle->SetPadRightMargin(RigMarg);
	gStyle->SetPadLeftMargin(LefMarg);

	//StatBox
	gStyle->SetOptStat("ei");
	gStyle->SetStatW(0.15);
	gStyle->SetStatX(.875);
	gStyle->SetStatY(1.-TopMarg);
	gStyle->SetStatColor(StCol);
	gStyle->SetStatFont(StFont);
	gStyle->SetStatTextColor(603);
	gStyle->SetStatFormat("8.7g");
	gStyle->SetFitFormat("8.7g");
	gStyle->SetOptFit(1111);
	gStyle->SetStatBorderSize(3);

	//Title
	gStyle->SetTitleX(0.45);
	gStyle->SetTitleFont(62, "XYZ");
	gStyle->SetTitleFont(62, ""   );
	gStyle->SetTitleFontSize(0.040);
	//Label
	gStyle->SetLabelFont(62, "XYZ");

	//Color Gradient
	const int NCont = 99;
	const int NRGBs = 5;
	double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	double red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
//	gStyle -> SetNumberContours(255);
	

//	const Int_t NRGBs = 10; 
//	Double_t Stops[NRGBs] = { 0.000, 0.012, 0.025, 0.050, 0.200, 0.250, 0.700, 0.800, 0.900, 1.000 };
//	Double_t Red[NRGBs]   = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 1.000, 1.000, 1.000 }; 
//	Double_t Green[NRGBs] = { 0.000, 0.020, 0.025, 0.050, 0.600, 0.800, 0.850, 0.900, 1.000, 1.000 };
//	Double_t Blue[NRGBs]  = { 0.640, 0.640, 0.720, 0.768, 0.800, 0.800, 0.250, 0.200, 0.200, 0.500 }; 
//	TColor::CreateGradientColorTable( NRGBs, Stops, Red, Green, Blue, NCont );

	gStyle -> SetNumberContours(NCont);

	float const mycol[4][3] = { { .900, .975, 1.   },
                                { 1.  , .950, .985 },
                                { .925, 1.  , .950 },
                                { 1.  , 1.  , .925 } };
	for(int i=0; i<4; i++){gROOT->GetColor(11+i)->SetRGB( mycol[i][0], mycol[i][1], mycol[i][2]);}
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Hist1D( TH1D *hist, TString HTitle, TString XTitle, TString YTitle, int LCol, int LSty, int Font, int FCol, int FSty ){
	hist->SetTitle(HTitle);
	hist->SetLineColor(LCol);
	hist->SetLineWidth(1);
	hist->SetTitleSize(0.04,"");
//	hist->SetTitleFont(Font,"");
	hist->SetTitleFont(Font, "");
	hist->SetFillStyle(FSty);
	hist->SetFillColor(FCol);

	hist->GetXaxis()->SetTitle(XTitle);
//	hist->GetXaxis()->CenterTitle();
	hist->GetXaxis()->SetTitleFont(Font);
	hist->GetXaxis()->SetTitleOffset(0.90);
	hist->GetXaxis()->SetTitleSize(0.05);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelOffset(0.01);

	hist->GetYaxis()->SetTitle(YTitle);
//	hist->GetYaxis()->CenterTitle();
	hist->GetYaxis()->SetTitleFont(Font);
	hist->GetYaxis()->SetTitleOffset(0.90);
	hist->GetYaxis()->SetTitleSize(0.05);
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelOffset(0.01);
	((TGaxis*)hist->GetYaxis())->SetMaxDigits(6);

}
//---------------------------------------------------------------------------------//
void Setting::Setting_Hist2D( TH2D *hist, TString HTitle, TString XTitle, TString YTitle, TString ZTitle, double Min ){
	hist->SetTitle(HTitle);
	hist->SetLineWidth(1);
	hist->SetTitleSize(0.04,"");
	hist->SetTitleFont(62,"");
	hist->SetMinimum(Min);

	hist->GetXaxis()->SetTitle(XTitle);
//	hist->GetXaxis()->CenterTitle();
	hist->GetXaxis()->SetTitleFont(62);
	hist->GetXaxis()->SetTitleOffset(0.90);
	hist->GetXaxis()->SetTitleSize(0.05);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelOffset(0.01);
	hist->GetXaxis()->SetNdivisions(515);

	hist->GetYaxis()->SetTitle(YTitle);
//	hist->GetYaxis()->CenterTitle();
	hist->GetYaxis()->SetTitleFont(62);
	hist->GetYaxis()->SetTitleOffset(0.60);
	hist->GetYaxis()->SetTitleSize(0.05);
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelOffset(0.01);
	hist->GetYaxis()->SetNdivisions(515);
	((TGaxis*)hist->GetYaxis())->SetMaxDigits(6);

	hist->GetZaxis()->SetTitle(ZTitle);
//	hist->GetZaxis()->CenterTitle();
	hist->GetZaxis()->SetTitleFont(62);
	hist->GetZaxis()->SetLabelFont(42);
	hist->GetZaxis()->SetLabelSize(0.030);
	hist->GetZaxis()->SetLabelOffset(0.008);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Graph( TGraph *gr, TString GTitle, TString XTitle, TString YTitle, int LCol, int LSty, int Font, int MCol, int MSty, double MSiz ){
	gr->SetTitle(GTitle);
	gr->SetName(GTitle);

	gr->GetXaxis()->SetTitle(XTitle);
	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitleFont(Font);
	gr->GetXaxis()->SetTitleOffset(0.90);
	gr->GetXaxis()->SetTitleSize(0.06);
	gr->GetXaxis()->SetLabelFont(Font);
	gr->GetXaxis()->SetLabelOffset(0.01);

	gr->GetYaxis()->SetTitle(YTitle);
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleFont(Font);
	gr->GetYaxis()->SetTitleOffset(1.00);
	gr->GetYaxis()->SetTitleSize(0.06);
	gr->GetYaxis()->SetLabelFont(Font);
	gr->GetYaxis()->SetLabelOffset(0.01);
	((TGaxis*)gr->GetYaxis())->SetMaxDigits(5);

	gr->SetLineColor(LCol);
	gr->SetLineStyle(LSty);
	gr->SetLineWidth(1);
	gr->SetMarkerStyle(MSty);
	gr->SetMarkerColor(MCol);
	gr->SetMarkerSize(1.25);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_GError( TGraphErrors *gr, TString GTitle, TString XTitle, TString YTitle, int LCol, int LSty, int Font, int MCol, int MSty, double MSiz ){
	gr->SetTitle(GTitle);
	gr->SetName(GTitle);

	gr->GetXaxis()->SetTitle(XTitle);
//	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitleFont(Font);
	gr->GetXaxis()->SetTitleOffset(0.80);
	gr->GetXaxis()->SetTitleSize(0.06);
	gr->GetXaxis()->SetLabelFont(42);
	gr->GetXaxis()->SetLabelSize(0.025);
	gr->GetXaxis()->SetLabelOffset(0.01);

	gr->GetYaxis()->SetTitle(YTitle);
//	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleFont(Font);
	gr->GetYaxis()->SetTitleOffset(0.80);
	gr->GetYaxis()->SetTitleSize(0.06);
	gr->GetYaxis()->SetLabelFont(42);
	gr->GetYaxis()->SetLabelSize(0.025);
	gr->GetYaxis()->SetLabelOffset(0.01);
	((TGaxis*)gr->GetYaxis())->SetMaxDigits(6);

	gr->SetLineColor(LCol);
	gr->SetLineStyle(LSty);
	gr->SetLineWidth(1);
	gr->SetMarkerStyle(MSty);
	gr->SetMarkerColor(MCol);
	gr->SetMarkerSize(MSiz);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_GError( TGraphAsymmErrors *gr, TString GTitle, TString XTitle, TString YTitle, int LCol, int LSty, int Font, int MCol, int MSty, double MSiz ){
	gr->SetTitle(GTitle);
	gr->SetName(GTitle);

	gr->GetXaxis()->SetTitle(XTitle);
//	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitleFont(Font);
	gr->GetXaxis()->SetTitleOffset(0.80);
	gr->GetXaxis()->SetTitleSize(0.06);
	gr->GetXaxis()->SetLabelFont(Font);
	gr->GetXaxis()->SetLabelOffset(0.01);

	gr->GetYaxis()->SetTitle(YTitle);
//	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleFont(Font);
	gr->GetYaxis()->SetTitleOffset(0.60);
	gr->GetYaxis()->SetTitleSize(0.06);
	gr->GetYaxis()->SetLabelFont(Font);
	gr->GetYaxis()->SetLabelOffset(0.01);
	((TGaxis*)gr->GetYaxis())->SetMaxDigits(5);

	gr->SetLineColor(LCol);
	gr->SetLineStyle(LSty);
	gr->SetLineWidth(1);
	gr->SetMarkerStyle(MSty);
	gr->SetMarkerColor(MCol);
	gr->SetMarkerSize(MSiz);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_G2D( TGraph2D *gr, TString GTitle, TString XTitle, TString YTitle, TString ZTitle, int Npx, int Npy, int LCol, int LSty, int Font, int MCol, int MSty, double MSiz, double Min ){
	gr->SetTitle(GTitle);
	gr->SetName(GTitle);

	gr->SetLineColor(LCol);
	gr->SetLineStyle(LSty);
	gr->SetLineWidth(1);
	gr->SetMarkerStyle(MSty);
	gr->SetMarkerColor(MCol);
	gr->SetMarkerSize(MSiz);

	if(Npx>0){
		gr->SetNpx(Npx);
	}else;
	if(Npy>0){
		gr->SetNpy(Npy);
	}else;

	gr->SetMinimum(Min);

	gr->GetXaxis()->SetTitle(XTitle);
	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitleFont(Font);
	gr->GetXaxis()->SetTitleOffset(0.80);
	gr->GetXaxis()->SetTitleSize(0.05);
	gr->GetXaxis()->SetLabelFont(42);
	gr->GetXaxis()->SetLabelSize(0.030);
	gr->GetXaxis()->SetLabelOffset(0.01);

	gr->GetYaxis()->SetTitle(YTitle);
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleFont(Font);
	gr->GetYaxis()->SetTitleOffset(0.80);
	gr->GetYaxis()->SetTitleSize(0.05);
	gr->GetYaxis()->SetLabelFont(42);
	gr->GetYaxis()->SetLabelSize(0.030);
	gr->GetYaxis()->SetLabelOffset(0.01);
	((TGaxis*)gr->GetYaxis())->SetMaxDigits(6);

	gr->GetZaxis()->SetTitle(ZTitle);
	gr->GetZaxis()->SetTitleFont(Font);
	gr->GetZaxis()->CenterTitle();
	gr->GetZaxis()->SetTitleSize(0.05);
	gr->GetZaxis()->SetLabelFont(42);
	gr->GetZaxis()->SetLabelSize(0.030);
	gr->GetZaxis()->SetLabelOffset(0.01);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_G2D( TGraph2DErrors *gr, TString GTitle, TString XTitle, TString YTitle, TString ZTitle, int Npx, int Npy, int LCol, int LSty, int Font, int MCol, int MSty, double MSiz, double Min ){
	gr->SetTitle(GTitle);
	gr->SetName(GTitle);

	gr->SetLineColor(LCol);
	gr->SetLineStyle(LSty);
	gr->SetLineWidth(1);
	gr->SetMarkerStyle(MSty);
	gr->SetMarkerColor(MCol);
	gr->SetMarkerSize(MSiz);

	if(Npx>0){
		gr->SetNpx(Npx);
	}else;
	if(Npy>0){
		gr->SetNpy(Npy);
	}else;

	gr->SetMinimum(Min);

	gr->GetXaxis()->SetTitle(XTitle);
	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitleFont(Font);
	gr->GetXaxis()->SetTitleOffset(0.80);
	gr->GetXaxis()->SetTitleSize(0.05);
	gr->GetXaxis()->SetLabelFont(42);
	gr->GetXaxis()->SetLabelSize(0.030);
	gr->GetXaxis()->SetLabelOffset(0.01);

	gr->GetYaxis()->SetTitle(YTitle);
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleFont(Font);
	gr->GetYaxis()->SetTitleOffset(0.80);
	gr->GetYaxis()->SetTitleSize(0.05);
	gr->GetYaxis()->SetLabelFont(42);
	gr->GetYaxis()->SetLabelSize(0.030);
	gr->GetYaxis()->SetLabelOffset(0.01);
	((TGaxis*)gr->GetYaxis())->SetMaxDigits(6);

	gr->GetZaxis()->SetTitle(ZTitle);
	gr->GetZaxis()->SetTitleFont(Font);
	gr->GetZaxis()->CenterTitle();
	gr->GetZaxis()->SetTitleSize(0.05);
	gr->GetZaxis()->SetLabelFont(42);
	gr->GetZaxis()->SetLabelSize(0.030);
	gr->GetZaxis()->SetLabelOffset(0.01);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Func( TF1 *func, int LCol, int LSty ){
	func->SetLineColor(LCol);
	func->SetLineStyle(LSty);
	func->SetLineWidth(1);
	func->SetNpx(1E+4);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Legend( TLegend *leg, int Font, int Al, int Col, double Size ){
	leg->SetTextSize(Size);
	leg->SetTextFont(Font);
	leg->SetTextAlign(Al);
	leg->SetTextColor(Col);
	leg->SetBorderSize(4);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Latex( TLatex *lat , int Font, int Al, int Col, double Size ){
	lat->SetTextSize(Size);
	lat->SetTextFont(Font);
	lat->SetTextAlign(Al);
	lat->SetTextColor(Col);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Line( TLine *lin, int LCol, int Wid, int Sty ){
	lin->SetLineColor(LCol);
	lin->SetLineWidth(Wid);
	lin->SetLineStyle(Sty);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Line( TPolyLine3D *lin, int LCol, int Wid, int Sty ){
	lin->SetLineColor(LCol);
	lin->SetLineWidth(Wid);
	lin->SetLineStyle(Sty);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Box( TBox* box, int LCol, int Wid, int Sty, int FCol, int FSty ){
	//Line attribute
	box->SetLineColor(LCol);
	box->SetLineWidth(Wid);
	box->SetLineStyle(Sty);

	//Fill attribute
	box->SetFillColor(FCol);
	box->SetFillStyle(FSty);
}
//---------------------------------------------------------------------------------//
void Setting::Setting_Pave( TPaveText* Pt, int Font, int TCol, int TAli, double TSize, int LCol, int LWid, int LSty, int FCol, int FSty, int BSize ){
	Pt->SetTextFont(Font);
	Pt->SetTextColor(TCol);
	Pt->SetTextAlign(TAli);
	Pt->SetTextSize(TSize);
	Pt->SetLineColor(LCol);
	Pt->SetLineWidth(LWid);
	Pt->SetLineStyle(LSty);
	Pt->SetFillColor(FCol);
	Pt->SetFillStyle(FSty);
	Pt->SetBorderSize(BSize);
}

