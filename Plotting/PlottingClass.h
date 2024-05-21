//******************************************************************************
// Plotting Class for exporting histograms (1D,2D,Ratios) in the root framework
//  Author: Nicolas Strangmann (nicolas.strangmann@cern.ch) Created: 10.10.2020
//******************************************************************************

#ifndef DRAWN
#define DRAWN

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TObjString.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TColor.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TEllipse.h"
#include "TCurlyLine.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include "TDatabasePDG.h"
#include "../Analysis/Logging.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++ Plotting ++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  The Plotting-class itself is never used but all the other classes inherit functions and attributes from it

double StdTextSize = 65;
double TitleSizeScaling = 1.2;

class Plotting
{
public:
  Plotting(); // Empty constructor

  ~Plotting(); // Empty destructor

  //  Universal function for setting the legends relative position on the canvas
  void SetLegend(double x1 = 0.15, double x2 = 0.4, double y1 = 0.7, double y2 = 0.9, bool l1hollow = false, Double_t InLegendTextMod = 1.);
  void SetLegend2(double x1 = 0.15, double x2 = 0.4, double y1 = 0.7, double y2 = 0.9, bool l2hollow = false, Double_t InLegend2TextMod = 1.);

  //  Setting a border to 42 triggers the AutoSetAxisRanges funtion
  void SetAxisRange(double xlow = 42, double xup = 42, double ylow = 42, double yup = 42, double zlow = 42, double zup = 42);

  //  Adds a string to vector Latex that will be drawn on canvas in function Plot(). Position in relative coordinates. Use ; to split string into seperate lines.
  void NewLatex(const double PositX = 0.2, const double PositY = 0.2, TString text = "", const double TextSize = StdTextSize, const double dDist = 0.06, const int font = 43, const int color = kBlack, const int angle = 0);

  //  Adds a line to vector lines that will be drawn on canvas in function Plot(). The lines coordinates relate to the axis. When setting negative style a curly line will be drawn instead of the straight one.
  void NewLine(double x1 = 0, double y1 = 0, double x2 = 1, double y2 = 1, int style = 1, int size = 2, int color = kBlack, TString label = "", bool InSecondLegend = false);

  //  Set the relative empty space between hist and the edges aswell as the canvas dimensions in pixels
  void SetMargins(double low = 0.12, double left = 0.1, double up = 0.025, double right = 0.025, int cw = 2000, int ch = 1500, double s = 1. / 3.);
  //   void SetMargins(double low = 0.1, double left = 0.1, double up = 0.02, double right = 0.02, int cw = 800, int ch = 600, double s = 1. / 3.);

  void AddEMCalOutline();
  void AddPHOSOutline();

  std::vector<TObject *> AdditionalPlottingObjects; // Manually added plotting objects (Additional legends or so)

protected:
  TCanvas *Canvas = nullptr; //  The canvas that all classes plot on

  TH2D *hDummy = nullptr; //  Empty histogram with the correct axis ranges and labels plotted first

  TLegend *leg = nullptr, *leg2 = nullptr;

  //  Vectors of elements that are added by New... functions and are drawn in the Plot() functions

  std::vector<TObject *> PlottingObjects;

  std::vector<TString> DrawOption;  //  A histogram is plotted using the corresponding DrawOption ("p","h",..)
  std::vector<TString> LegendLabel; //  Strings corresponding to histograms are added to legend

  //  The following are standard settings that can be changes by calling Set.. functions before Plot()
  double AxisRange[3][2] = {{42, 42}, {42, 42}, {42, 42}}; // xlow,xuo,ylow,yup,zlow,zup (in PlottingRatio z=ratio)
  TString AxisLabel[3] = {"x", "y", "Ratio"};
  double AxisLabelOffset[2] = {1., 1.25};                  // Third component not needed for ratio trivial, for 2D z has no Label currently
  double LegendBorders[2][2] = {{0.2, 0.4}, {0.7, 0.9}};   //  xlow,xup,ylow,yup in relative units (0-1)
  double Legend2Borders[2][2] = {{0.7, 0.95}, {0.7, 0.9}}; //  xlow,xup,ylow,yup in relative units (0-1)
  bool L1hollow = true, L2hollow = true, Leg2ManSet = false;
  double CanvasMargins[2][2] = {{0.1, 0.025}, {0.12, 0.025}}; //  left,right,low,up in relative units
  //   int CanvasDimensions[2] = {800, 600};                 // Dimension given in pixels (Standard defined by the ALICE collaboration)
  int CanvasDimensions[2] = {2000, 1500}; // Dimension given in pixels (Same ratio as ALICE standard but twice as large for better png quality)
  double Split = 1. / 3.;
  int Ndigits = 3;
  Double_t LegendTextMod = 1.;
  Double_t Legend2TextMod = 1.;

  //  If no style and or color are set these 10 standard styles and colors are used one after the other
  int AutoStyle[10] = {20, 21, 34, 33, 27, 24, 28, 22, 23, 29};
  int AutoStyleLine[10] = {1, 7, 9, 2, 8, 1, 7, 9, 2, 8};
  int AutoColor[10] = {kBlue + 1, kRed + 1, kGreen + 2, kBlack, kOrange + 2, kCyan + 3, kTeal - 7, kPink + 2, kYellow + 3, kSpring + 4};
  int counter = 0; //  Count up after each added histogram and use e.g. AutoStyle[counter] -> Unique colors and styles

  //  Create the legend leg using LegendBorders that will be drawn in Plot()
  void InitializeLegend();
  void SetLeg2Position();

  //  Adjusts the x and y axis range depending on the histograms that will be drawn
  void AutoSetAxisRanges(bool logy);

  template <typename T>
  void SetStyle(T x, int style, double size, int color, TString opt, int Counter);

  //  Converts the given DrawOptions to good parametes for the legend reference symbols
  TString LegendDrawOption(TString DrawOpt);
};

Plotting::Plotting()
{
}

Plotting::~Plotting()
{
}

//  The following 3 functions simply copy the user given settings into attributes of the Plotting class
void Plotting::SetMargins(double low, double left, double up, double right, int cw, int ch, double s)
{
  CanvasMargins[0][0] = left;
  CanvasMargins[0][1] = right;
  CanvasMargins[1][0] = low;
  CanvasMargins[1][1] = up;
  CanvasDimensions[0] = cw;
  CanvasDimensions[1] = ch;
  Split = s;
} //  These parameters will be used when Plot() calls InitializeCanvas

void Plotting::SetAxisRange(double xlow, double xup, double ylow, double yup, double zlow, double zup)
{
  AxisRange[0][0] = xlow;
  AxisRange[0][1] = xup;
  AxisRange[1][0] = ylow;
  AxisRange[1][1] = yup;
  AxisRange[2][0] = zlow; //  Not used in Plotting 1D
  AxisRange[2][1] = zup;  // In PlottingRatio this gives the range of the ratio
} //  These parameters will be used when Plot() calls InitializeAxis

void Plotting::SetLegend(double x1, double x2, double y1, double y2, bool l1hollow, Double_t InLegendTextMod)
{
  LegendBorders[0][0] = x1;
  LegendBorders[0][1] = x2;
  LegendBorders[1][0] = y1;
  LegendBorders[1][1] = y2;
  L1hollow = l1hollow;
  LegendTextMod = InLegendTextMod;
} //  These parameters will be used when Plot() calls InitializeLegend

void Plotting::SetLegend2(double x1, double x2, double y1, double y2, bool l2hollow, Double_t InLegend2TextMod)
{
  Legend2Borders[0][0] = x1;
  Legend2Borders[0][1] = x2;
  Legend2Borders[1][0] = y1;
  Legend2Borders[1][1] = y2;
  L2hollow = l2hollow;
  Leg2ManSet = true;
  Legend2TextMod = InLegend2TextMod;
} //  These parameters will be used when Plot() calls InitializeLegend

void Plotting::SetLeg2Position()
{
  Legend2Borders[0][0] = LegendBorders[0][0] - 0.04;
  Legend2Borders[0][1] = LegendBorders[0][0];
  Legend2Borders[1][0] = LegendBorders[1][0];
  Legend2Borders[1][1] = LegendBorders[1][1];
  L2hollow = L1hollow;
} //  These parameters will be used when Plot() calls InitializeLegend

void Plotting::InitializeLegend()
{
  leg = new TLegend(LegendBorders[0][0], LegendBorders[1][0], LegendBorders[0][1], LegendBorders[1][1]);
  leg->SetHeader(""); //  Remove title of legend
  leg->SetTextFont(43);
  leg->SetTextSize(LegendTextMod * StdTextSize);
  leg->SetBorderSize(0);   //  Remove black rectangle around legend
  leg->SetFillStyle(1001); // Solid white background to make legend readable. Set to 0 to make it hollow
  leg->SetTextAlign(12);   // Align text left - center
  if (L1hollow)
    leg->SetFillStyle(0);

  bool bLegendTable = false;
  for (int iLegEntry = 0; iLegEntry < (int)LegendLabel.size(); iLegEntry++)
  {
    if (LegendLabel.at(iLegEntry).BeginsWith("42") && !Leg2ManSet)
      bLegendTable = true;
  }
  if (bLegendTable)
    SetLeg2Position();

  leg2 = new TLegend(Legend2Borders[0][0], Legend2Borders[1][0], Legend2Borders[0][1], Legend2Borders[1][1]);
  leg2->SetHeader(""); //  Remove title of legend
  leg2->SetTextFont(43);
  leg2->SetTextSize(Legend2TextMod * StdTextSize);
  leg2->SetBorderSize(0);   //  Remove black rectangle around legend
  leg2->SetFillStyle(1001); // Solid white background to make legend readable. Set to 0 to make it hollow
  leg2->SetTextAlign(12);   // Align text left - center
  if (L2hollow)
    leg2->SetFillStyle(0);
  if (bLegendTable)
    leg2->SetMargin(1);
  // leg->SetMargin(0.1 / (LegendBorders[0][1] - LegendBorders[0][0]));
}

void Plotting::AutoSetAxisRanges(bool logy)
{
  double Range[3][2] = {{1E12, -1E12}, {1E12, -1E12}, {1E12, -1E12}}; // Set the minima large enough and the maxima low enough, so that they will be overwritten in the first iteration of the following object loop

  for (int io = 0; io < (int)PlottingObjects.size(); io++)
  {
    TString ClassName = (TString)PlottingObjects.at(io)->ClassName();
    TH1F *h = nullptr;
    if (ClassName.Contains("TH1"))
    {
      h = (TH1F *)PlottingObjects.at(io);
    }
    else if (ClassName.Contains("TGraph"))
    {
      h = ((TGraph *)PlottingObjects.at(io))->GetHistogram();
    }
    else
      continue;
    Range[0][0] = h->GetBinLowEdge(1) < Range[0][0] ? h->GetBinLowEdge(1) : Range[0][0];
    Range[0][1] = h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX()) > Range[0][1] ? h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX()) : Range[0][1];
    Range[1][0] = (h->GetMinimum() < Range[1][0] && (!logy || h->GetMinimum() > 0)) ? h->GetMinimum() : Range[1][0];
    Range[1][1] = h->GetMaximum() > Range[1][1] ? h->GetMaximum() : Range[1][1];
  }

  Range[1][1] = logy ? 2 * Range[1][1] : Range[1][1] + (Range[1][1] - Range[1][0]) / 10;                                                                                    //  With log scales a factor 2 isn't too much
  Range[1][0] = logy ? 0.5 * Range[1][0] : (Range[1][1] - 2 * Range[1][0] > 0 ? (Range[1][0] > 0 ? 0 : 1.1 * Range[1][0]) : Range[1][0] - (Range[1][1] - Range[1][0]) / 8); //  max-2*min>0 ? min is small->go down to 0 : all bins rather full

  for (int iA = 0; iA < 2; iA++)
  {
    for (int iE = 0; iE < 2; iE++)
    {
      if (AxisRange[iA][iE] > 41.99 && AxisRange[iA][iE] < 42.01)
        AxisRange[iA][iE] = Range[iA][iE];
    }
  }
  if (AxisRange[1][0] == AxisRange[1][1])
    AxisRange[1][1] = AxisRange[1][0] + 1;
}

void Plotting::NewLatex(const double PositX, const double PositY, TString text, const double TextSize, const double dDist, const int font, const int color, const int angle)
{
  std::vector<TString> LatStr;             //  Each element corresponds to a line of the printed latex string
  TObjArray *textStr = text.Tokenize(";"); //  The semicolon seperates the string into different lines
  for (int i = 0; i < textStr->GetEntries(); i++)
  {
    TObjString *tempObj = (TObjString *)textStr->At(i);
    LatStr.push_back(tempObj->GetString());
  }

  //  Loop thru the latex lines and set the formatting
  for (int i = 0; i < (int)LatStr.size(); ++i)
  {
    TLatex *latex = new TLatex(PositX, PositY - i * dDist, LatStr[i]);
    latex->SetNDC();
    latex->SetTextFont(font);
    latex->SetTextColor(color);
    latex->SetTextSize(TextSize);
    latex->SetTextAngle(angle);
    if (PositX > 0.81)
      latex->SetTextAlign(32);
    PlottingObjects.push_back(latex);
    DrawOption.push_back("");
    LegendLabel.push_back("");
  }
}

void Plotting::NewLine(double x1, double y1, double x2, double y2, int style, int size, int color, TString label, bool InSecondLegend)
{

  //  Curly lines can be used to draw photons or similar
  if (style < 0)
  {
    TCurlyLine *line = new TCurlyLine(x1, y1, x2, y2);
    line->SetLineColor(color);
    line->SetLineWidth(size);
    line->SetWaveLength(-0.02 * style); //  Standard wavelength is 0.02 -> Style -1
    PlottingObjects.push_back(line);
  }
  else
  {
    TLine *line = new TLine(x1, y1, x2, y2);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(size);
    PlottingObjects.push_back(line);
  }
  LegendLabel.push_back(Form("%s%s", InSecondLegend ? "42" : "", label.Data()));
  DrawOption.push_back("l");
}

void Plotting::AddEMCalOutline()
{
  int lineWidth = 1;
  int lineColor = kGray + 2;
  int lineStyle = 1;
  NewLine(-0.7, 1.4, -0.7, 3.26, lineStyle, lineWidth, lineColor);
  NewLine(-0.7, 1.4, +0.7, 1.4, lineStyle, lineWidth, lineColor);
  NewLine(+0.7, 1.4, +0.7, 3.26, lineStyle, lineWidth, lineColor);
  NewLine(-0.7, 3.26, +0.7, 3.26, lineStyle, lineWidth, lineColor);

  NewLine(-0.7, 5.7, +0.7, 5.7, lineStyle, lineWidth, lineColor);     // rechts
  NewLine(+0.7, 4.54, +0.7, 5.7, lineStyle, lineWidth, lineColor);    // oben
  NewLine(-0.7, 4.54, -0.7, 5.7, lineStyle, lineWidth, lineColor);    // unten
  NewLine(+0.23, 4.54, +0.23, 5.58, lineStyle, lineWidth, lineColor); // "oben"
  NewLine(-0.23, 4.54, -0.23, 5.58, lineStyle, lineWidth, lineColor); // "unten"
  NewLine(-0.23, 5.58, +0.23, 5.58, lineStyle, lineWidth, lineColor); // "links"
  NewLine(-0.7, 4.54, -0.23, 4.54, lineStyle, lineWidth, lineColor);  // links oben
  NewLine(0.7, 4.54, +0.23, 4.54, lineStyle, lineWidth, lineColor);   // links unten
}

void Plotting::AddPHOSOutline()
{
  int lineWidth = 1;
  int lineColor = kGray + 2;
  int lineStyle = 2;

  NewLine(-0.12, 4.36, +0.12, 4.36, lineStyle, lineWidth, lineColor); // links
  NewLine(-0.12, 4.36, -0.12, 5.58, lineStyle, lineWidth, lineColor); // links
  NewLine(+0.12, 4.36, +0.12, 5.58, lineStyle, lineWidth, lineColor); // links
}

TString Plotting::LegendDrawOption(TString UserDrawOpt)
{
  TString LegendDrawOpt = "p"; //  If user gives no DrawOption use p as standard
  if ((int)*(UserDrawOpt.Data()))
    LegendDrawOpt = UserDrawOpt;
  if (UserDrawOpt.Contains("h"))
    LegendDrawOpt = "l"; //  Histogram style hists should have a line in the legend
  if (UserDrawOpt.Contains("E1"))
    LegendDrawOpt = "pE1"; //  Histogram style hists should have a line in the legend
  if (UserDrawOpt.Contains("E1") && UserDrawOpt.Contains("f"))
    LegendDrawOpt = "fpE1"; //  Histogram style hists should have a line in the legend
  if (UserDrawOpt.Contains("E1") && UserDrawOpt.Contains("z"))
    LegendDrawOpt = "fpE1"; //  Histogram style hists should have a line in the legend
  if (UserDrawOpt.Contains("E1") && UserDrawOpt.Contains("y"))
    LegendDrawOpt = "p"; //  Histogram style hists should have a line in the legend
  if (UserDrawOpt.Contains("E0") && UserDrawOpt.Contains("y"))
    LegendDrawOpt = "p"; //  Histogram style hists should have a line in the legend
  if (UserDrawOpt.Contains("EX0"))
    LegendDrawOpt = "pE1"; //  Histogram style hists should have a line in the legend
  return LegendDrawOpt;
}

template <typename T>
void Plotting::SetStyle(T x, int style, double size, int color, TString opt, int Counter)
{
  x->SetMarkerStyle((style == -1) ? AutoStyle[Counter] : style);
  x->SetLineStyle((opt.EndsWith("2") && !opt.EndsWith("E2")) ? 2 : ((opt.Contains("h") || opt.Contains("l") || opt.Contains("E2")) && style < 10 && style > 0) ? style
                                                                                                                                                               : (style == -1 && (opt.Contains("h") || opt.Contains("l") || opt.Contains("E2")) ? AutoStyleLine[Counter] : 1));
  if (opt.Contains("f"))
    x->SetFillColorAlpha((color == -1) ? AutoColor[Counter] : color, 0.7);
  x->SetMarkerColor((color == -1) ? AutoColor[Counter] : color);
  x->SetLineColor((color == -1) ? AutoColor[Counter] : color);
  x->SetMarkerSize(size);
  x->SetLineWidth((opt.EndsWith("2") && !opt.EndsWith("E2")) || opt.Contains("E1p") ? size - 1 : (opt.Contains("lp") ? size / 2. : size));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++ Plotting 1D ++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Plotting1D : public Plotting
{
public:
  Plotting1D(); // This empty constructor always has to be called when plotting 1D

  ~Plotting1D(); // Destructor

  //  After adding all histograms, functions and graph using the New.. functions create the actual plot
  void Plot(TString name = "dummy.pdf", bool logx = false, bool logy = false, bool TransCanvas = false);

  //  Add a histogram to the hists vector and put all its settings into different vectors
  void New(TH1 *h = nullptr, TString label = "", int style = -1, double size = 2, int color = -1, TString opt = "p", bool InSecondLegend = false);
  void New(TF1 *f = nullptr, TString label = "", int style = -1, double size = 2, int color = -1, TString opt = "l", bool InSecondLegend = false);
  void New(TGraphAsymmErrors *g = nullptr, TString label = "", int style = -1, double size = 2, int color = -1, TString opt = "p", bool InSecondLegend = false);

  void Add(TObject *AdditionalObject);

  //  Store the user wishes for labels and offsets in the AxisLabel and AxisLabelOffset attributes. They will later be used in InitializeAxis.
  void SetAxisLabel(TString labelx = "", TString labely = "", double offsetx = 1., double offsety = 1.25, int ndigits = 3);

private:
  //  Create the canvas using the standard dimensions and margins, if they were not set by SetMargins
  void InitializeCanvas(bool logx, bool logy, bool TransCanvas);

  //  Create the hdummy that will be plotted first and give it the Set xis ranges and labels
  void InitializeAxis(bool logx, bool logy, bool TransCanvas);
};

Plotting1D::Plotting1D()
{
}

Plotting1D::~Plotting1D()
{
}

void Plotting1D::Plot(TString name, bool logx, bool logy, bool TransCanvas)
{
  if (PlottingObjects.size() < 1)
    ERRORETURN("Nothing added for plotting.")

  InitializeCanvas(logx, logy, TransCanvas); //  Creating Canvas with margins
  InitializeAxis(logx, logy, TransCanvas);   //  Create the hDummy and set its axis label + ranges
  hDummy->Draw();                            //  Draw the just set axis (label) on the Canvas
  InitializeLegend();                        //  Create leg and set its dimensions + format

  //  Loop thru all elements of all vectors and plot them on top of the empty hDummy and add them to leg
  //----------------------------------------------------------------------------
  for (int i = 0; i < (int)PlottingObjects.size(); ++i)
  {
    PlottingObjects.at(i)->Draw(Form("same %s", ((TString)DrawOption.at(i)).Data()));
    if (LegendLabel.at(i).BeginsWith("42"))
    {
      LegendLabel.at(i).Replace(0, 2, "");
      leg2->AddEntry(PlottingObjects.at(i), LegendLabel.at(i) == "" ? " " : LegendLabel.at(i).Data(), LegendDrawOption(DrawOption.at(i)));
    }
    else if ((int)*(LegendLabel.at(i).Data()))
      leg->AddEntry(PlottingObjects.at(i), LegendLabel.at(i).Data(), LegendDrawOption(DrawOption.at(i)).Data());
  }
  //----------------------------------------------------------------------------
  //  Now that everything is drawn just add the legend and print it.
  if (leg->GetNRows() > 1)
  {
    if (leg->GetNRows() > 10)
      leg->SetTextSize(leg->GetTextSize() / 1.5);
    leg->Draw("same");
  }
  if (leg2->GetNRows() > 1)
  {
    if (leg2->GetNRows() > 10)
      leg2->SetTextSize(leg2->GetTextSize() / 1.5);
    leg2->Draw("same");
  }
  for (int i = 0; i < (int)AdditionalPlottingObjects.size(); i++)
    AdditionalPlottingObjects.at(i)->Draw("same");

  Canvas->SaveAs(name);
  delete hDummy;
  delete Canvas;
}

void Plotting1D::New(TH1 *h, TString label, int style, double size, int color, TString opt, bool InSecondLegend)
{
  if (!h)
    ERRORETURN("New was given a invalid (nullptr) TH1.");
  PlottingObjects.push_back(h);
  SetStyle<TH1 *>(h, style, size, color, opt, counter);
  DrawOption.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  LegendLabel.push_back(Form("%s%s", InSecondLegend ? "42" : "", label.Data()));
  if ((style == -1) || (color == -1))
    counter++;
}

void Plotting1D::New(TF1 *f, TString label, int style, double size, int color, TString opt, bool InSecondLegend)
{
  if (!f)
    ERRORETURN("New was given a invalid (nullptr) TF1.");
  PlottingObjects.push_back(f);
  SetStyle<TF1 *>(f, style, size, color, opt, counter);
  DrawOption.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  LegendLabel.push_back(Form("%s%s", InSecondLegend ? "42" : "", label.Data()));
}

void Plotting1D::New(TGraphAsymmErrors *g, TString label, int style, double size, int color, TString opt, bool InSecondLegend)
{
  if (!g)
    ERRORETURN("New was given a invalid (nullptr) TGraphAsymmError.");
  PlottingObjects.push_back(g);
  SetStyle<TGraphAsymmErrors *>(g, style, size, color, opt, counter);
  LegendLabel.push_back(Form("%s%s", InSecondLegend ? "42" : "", label.Data()));
  DrawOption.push_back(opt);
  if ((style == -1) || (color == -1))
    counter++;
}

void Plotting1D::InitializeCanvas(bool logx, bool logy, bool TransCanvas)
{
  Canvas = new TCanvas("Canvas", "Canvas", CanvasDimensions[0], CanvasDimensions[1]);
  Canvas->SetLeftMargin(CanvasMargins[0][0]);
  Canvas->SetRightMargin(CanvasMargins[0][1]);
  Canvas->SetBottomMargin(CanvasMargins[1][0]);
  Canvas->SetTopMargin(CanvasMargins[1][1]);
  if (TransCanvas)
    Canvas->SetFillStyle(0);

  //  Set ticks at regular intervals on every edge of the histogram (also right and top)
  gPad->SetTickx();
  gPad->SetTicky();

  Canvas->cd();
  Canvas->SetLogx(logx);
  Canvas->SetLogy(logy);
}

void Plotting1D::SetAxisLabel(TString labelx, TString labely, double offsetx, double offsety, int ndigits)
{
  AxisLabel[0] = labelx;
  AxisLabel[1] = labely;
  AxisLabelOffset[0] = offsetx;
  AxisLabelOffset[1] = offsety;
  Ndigits = ndigits;
}

void Plotting1D::InitializeAxis(bool logx, bool logy, bool TransCanvas)
{
  AutoSetAxisRanges(logy); //  If any AxisRanges are still set to 42 -> Autoset them

  hDummy = new TH2D("hDummy", "hDummy", 1000, AxisRange[0][0], AxisRange[0][1], 1000, AxisRange[1][0], AxisRange[1][1]);
  TAxis *Axis[3] = {hDummy->GetXaxis(), hDummy->GetYaxis(), hDummy->GetZaxis()};
  hDummy->SetTitle("");
  hDummy->SetStats(0);
  if (logx && AxisRange[0][0] > 0.05 && AxisRange[0][1] < 100)
    hDummy->SetLabelOffset(-0.01);
  if (logx)
    hDummy->SetLabelOffset(-0.005);
  hDummy->GetYaxis()->SetMaxDigits(Ndigits);

  for (int iA = 0; iA < 3; iA++)
  {
    Axis[iA]->SetRangeUser(AxisRange[iA][0], AxisRange[iA][1]);
    Axis[iA]->SetTitleSize(TitleSizeScaling * StdTextSize);
    Axis[iA]->SetLabelSize(StdTextSize);
    Axis[iA]->SetLabelFont(43);
    Axis[iA]->SetTitleFont(43);
    if (TransCanvas)
    {
      Axis[iA]->SetTitleColor(kWhite);
      Axis[iA]->SetLabelColor(kWhite);
    }
    if (iA == 2)
      continue;
    Axis[iA]->SetTitle(AxisLabel[iA]);
    Axis[iA]->SetTitleOffset(AxisLabelOffset[iA]);
  }
  // if (logx)
  //   hDummy->GetXaxis()->SetMoreLogLabels();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++ Plotting 2D ++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Plotting2D : public Plotting
{
public:
  Plotting2D(); // This empty constructor always has to be called when plotting 2D

  ~Plotting2D(); //  Destructor

  //  After adding the histograms using the NewHist function create the actual plot
  void Plot(TString name = "dummy.pdf", bool logx = false, bool logy = false, bool logz = false, bool TransCanvas = false, int numcontours = 100);

  //  The standard palette is kBird, but there are also other nice 2D plotting styles (https://root.cern.ch/doc/master/classTColor.html)
  void New(TH2F *h = nullptr, TString opt = "COLZ", int palette = kBird);
  void New(TH2D *h = nullptr, TString opt = "COLZ", int palette = kBird); //  Convert to TH2F and call that NewHist function

  void New(TGraphAsymmErrors *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p");

  //  Add a new function to the funcs vector that will be drawn when calling Plot()
  void New(TF1 *f = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "l");
  void New(TF2 *f = nullptr, TString opt = "l");
  // void NewFunc(TF2* f = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "l");

  void SetAxisLabel(TString labelx = "", TString labely = "", double offsetx = 1., double offsety = 1.);

protected:
  TH2F *hist = nullptr;
  TF2 *func2 = nullptr;
  TString DrawOptionf2 = "";

  void InitializeCanvas(bool logx, bool logy, bool logz, bool TransCanvas);

  void InitializeAxis(bool logx, bool TransCanvas);
};

Plotting2D::Plotting2D()
{
  CanvasMargins[0][1] *= 4.5; //  To leave room for the z axis
}

Plotting2D::~Plotting2D()
{
}

void Plotting2D::Plot(TString name, bool logx, bool logy, bool logz, bool TransCanvas, int numcontours)
{

  if (!hist)
    FATAL("No hist added for plotting.");

  InitializeCanvas(logx, logy, logz, TransCanvas); // Creating Canvas with margins
  InitializeAxis(logx, TransCanvas);
  InitializeLegend();
  gStyle->SetNumberContours(numcontours);

  hist->Draw(Form("same,%s", ((TString)DrawOption.at(0)).Data()));

  for (int i = 0; i < (int)PlottingObjects.size(); ++i)
  {
    PlottingObjects.at(i)->Draw(Form("same %s", ((TString)DrawOption.at(i + 1)).Data()));
    if ((int)*(LegendLabel.at(i).Data()))
      leg->AddEntry(PlottingObjects.at(i), LegendLabel.at(i).Data(), LegendDrawOption(DrawOption.at(i + 1)).Data());
  }

  if (func2)
    func2->Draw(Form("same,%s", ((TString)DrawOptionf2).Data()));

  if (leg->GetNRows() > 1)
    leg->Draw("same");

  Canvas->SaveAs(name);
  delete Canvas;
}

void Plotting2D::New(TH2F *h, TString opt, int palette)
{
  if (!h)
    FATAL("NewHist was given a nullptr.");
  hist = h;
  gStyle->SetPalette(palette);
  DrawOption.push_back((opt == "p") ? "p" : opt);
}

void Plotting2D::New(TH2D *h, TString opt, int palette)
{
  TH2F *hd = (TH2F *)h;
  New(hd, opt, palette);
}

void Plotting2D::New(TGraphAsymmErrors *g, TString label, int style, double size, int color, TString opt)
{
  if (!g)
    FATAL("New was given a invalid (nullptr) TGraphAsymmError.");
  PlottingObjects.push_back(g);
  SetStyle<TGraphAsymmErrors *>(g, style, size, color, opt, counter);
  DrawOption.push_back((opt == "p") ? "p" : opt);
  LegendLabel.push_back(label);
  if ((style == -1) || (color == -1))
    counter++;
}

void Plotting2D::New(TF1 *f, TString label, int style, double size, int color, TString opt)
{
  if (!f)
    FATAL("New was given a invalid (nullptr) TF1.");
  PlottingObjects.push_back(f);
  SetStyle<TF1 *>(f, style, size, color, opt, counter);
  LegendLabel.push_back(label);
  DrawOption.push_back(opt);
  // DrawOption.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  if ((style == -1) || (color == -1))
    counter++;
}

void Plotting2D::New(TF2 *f, TString opt)
{
  if (!f)
    FATAL("New was given a invalid (nullptr) TF2.");

  func2 = f;
  // gStyle->SetPalette(palette);
  DrawOptionf2 = opt;
}

void Plotting2D::InitializeCanvas(bool logx, bool logy, bool logz, bool TransCanvas)
{
  if (Canvas)
    delete Canvas; //  This should never happen, but better safe than sorry.

  Canvas = new TCanvas("Canvas", "Canvas", CanvasDimensions[0], CanvasDimensions[1]);
  Canvas->SetLeftMargin(CanvasMargins[0][0]);
  Canvas->SetRightMargin(CanvasMargins[0][1]);
  Canvas->SetBottomMargin(CanvasMargins[1][0]);
  Canvas->SetTopMargin(CanvasMargins[1][1]);

  //  Set ticks at regular intervals on every edge of the histogram (also right and top)
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetOptStat(0);

  if (TransCanvas)
    Canvas->SetFillStyle(0);

  Canvas->cd();
  Canvas->SetLogx(logx);
  Canvas->SetLogy(logy);
  Canvas->SetLogz(logz);
}

void Plotting2D::SetAxisLabel(TString labelx, TString labely, double offsetx, double offsety)
{
  AxisLabel[0] = labelx;
  AxisLabel[1] = labely;
  AxisLabelOffset[0] = offsetx;
  AxisLabelOffset[1] = offsety;
}

void Plotting2D::InitializeAxis(bool logx, bool TransCanvas)
{
  hist->SetStats(0);
  hist->SetTitle("");
  if (logx && AxisRange[0][0] > 0.05 && AxisRange[0][1] < 100)
    hist->SetLabelOffset(-0.01);
  hist->GetYaxis()->SetMaxDigits(3);
  TAxis *Axis[3] = {hist->GetXaxis(), hist->GetYaxis(), hist->GetZaxis()};
  for (int iA = 0; iA < 3; iA++)
  {
    Axis[iA]->SetTitleSize(TitleSizeScaling * StdTextSize);
    Axis[iA]->SetLabelSize(StdTextSize);
    Axis[iA]->SetLabelFont(43);
    Axis[iA]->SetTitleFont(43);

    if (TransCanvas)
    {
      Axis[iA]->SetTitleColor(kWhite);
      Axis[iA]->SetLabelColor(kWhite);
    }
    if (iA == 2 && AxisRange[2][0] == 42)
      continue;
    Axis[iA]->SetRangeUser(AxisRange[iA][0], AxisRange[iA][1]);
    if (iA == 2)
      continue;
    Axis[iA]->SetTitle(AxisLabel[iA]);
    Axis[iA]->SetTitleOffset(iA == 1 ? 1.1 * AxisLabelOffset[iA] : AxisLabelOffset[iA]);
  }
  if (logx)
    hist->GetXaxis()->SetMoreLogLabels();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++ Plotting Ratio +++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PlottingRatio : public Plotting
{
public:
  PlottingRatio(); // This empty constructor always has to be called when plotting Ratio

  ~PlottingRatio(); //  Destructor

  void Plot(TString name = "dummy.pdf", bool logx = false, bool logy = false, bool logz = false, bool TransCanvas = false);

  //  Add histograms to the upper pad
  void New(TH1 *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p", bool InSecondLegend = kFALSE);
  void New(TGraphAsymmErrors *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p");
  void New(TF1 *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "l");

  void NewRatioLine(double x1 = 2, double y1 = 1, double x2 = 20, double y2 = 1, int style = 1, int size = 1, int color = kBlack, TString label = "");

  //  Add histograms to the lower pad
  void NewRatio(TH1 *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p");
  void NewRatio(TGraphAsymmErrors *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p");
  void NewRatio(TF1 *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "l");

  //  To remove the label conflict where y and ratio axis meet, add a white box there. This function can move that box (e.g. when margins are changed) or set to red to visualize the pad.
  void SetWhite(double low = 0.32, double left = 0.04, double up = 0.35, double right = 0.095, bool red = false);

  //  The ratios have a seperate legend. It's position can be set analogously to SetLegend using this function
  void SetLegendR(double x1 = 0.7, double x2 = 0.95, double y1 = 0.15, double y2 = 0.25, bool bhollow = false);

  //  The ratio label can be set indiviually, but its offset is set via the y axis label offset
  void SetAxisLabel(TString labelx = "", TString labely = "", TString labelz = "", double offsetx = 1., double offsety = 1., int ndigits = 3);

protected:
  //  Add the pads that are drawn on the canvas in Plot() as attributes
  TPad *HistoPad = nullptr; //  The upper pad containing all hists added by NewHist and TopFuncs
  TPad *RatioPad = nullptr; //  The lower pad containing all ratios added by NewRatio and BotFuncs
  TPad *WhitePad = nullptr; //  A white filled rectangle hiding the y axis label conflict at their border

  double WhiteBorders[2][2] = {{0.05, 0.13}, {0.365, 0.39}}; //  Defines where the white pad is drawn to hide the label conflict
  bool wred = false;                                         //  Can be set to true to make WhitePad red (visible)

  TH2D *rDummy = nullptr; // Dummy for the ratio histogram analogously to hdummy
  TLegend *legR = nullptr;
  bool LRhollow = false;
  double RatioLegendBorders[2][2] = {{0.7, 0.95}, {0.15, 0.25}};

  //  Vectors containing the functions for the upper pad as well as the functions and ratios for the lower pad
  std::vector<TObject *> RatiObjects;
  std::vector<TString> LegendLabelR; //  Ratio
  std::vector<TString> DrawOptionR;

  //  The bottom pad has its own counter because it can and often should have similar colors as the top.
  //  Starts at one, because there is often no ratio to standard, so now the top histos have the same colors and styles as their bot ratios.
  int counteR = 0;

  //  Creates all three pads and the canvas that they are on. The canvas dimension attributes are NOT used, instead a standard size of 1000x1000 is used
  void InitializeCanvas(bool logx, bool logy, bool logz, bool TransCanvas);

  void InitializeAxis(bool logx, bool logy, bool TransCanvas);

  void InitializeLegendR(); //  Creates the legR and sets its coordinates according to RatioLegendBorders
};

PlottingRatio::PlottingRatio()
{
  double s = 4. / 3.;
  CanvasMargins[0][0] *= s;
  CanvasMargins[1][0] /= s;
  CanvasMargins[0][1] *= s;
  CanvasMargins[1][1] /= s;
}

PlottingRatio::~PlottingRatio()
{
}

void PlottingRatio::Plot(TString name, bool logx, bool logy, bool logz, bool TransCanvas)
{

  if (PlottingObjects.size() < 1)
    FATAL("No hists added for plotting.");
  if (RatiObjects.size() < 1)
    FATAL("No ratios added for plotting.");

  InitializeCanvas(logx, logy, logz, TransCanvas); // Creating Canvas with margins
  InitializeAxis(logx, logy, TransCanvas);
  hDummy->GetYaxis()->SetMaxDigits(Ndigits);
  hDummy->Draw();

  InitializeLegend();
  InitializeLegendR();
  //  Print all hists and top funcs on the HistoPad
  //----------------------------------------------------------------------------
  // for (int i = 0; i < (int)PlottingObjects.size(); ++i) {
  //   PlottingObjects.at(i)->Draw(Form("same %s", ((TString)DrawOption.at(i)).Data()));
  //   if ((int)*(LegendLabel.at(i).Data()))
  //     leg->AddEntry(PlottingObjects.at(i), LegendLabel.at(i).Data(), LegendDrawOption(DrawOption.at(i)).Data());
  // }

  for (int i = 0; i < (int)PlottingObjects.size(); ++i)
  {
    PlottingObjects.at(i)->Draw(Form("same %s", ((TString)DrawOption.at(i)).Data()));
    if (LegendLabel.at(i).BeginsWith("42"))
    {
      LegendLabel.at(i).Replace(0, 2, "");
      leg2->AddEntry(PlottingObjects.at(i), LegendLabel.at(i) == "" ? " " : LegendLabel.at(i).Data(), LegendDrawOption(DrawOption.at(i)));
    }
    else if ((int)*(LegendLabel.at(i).Data()))
      leg->AddEntry(PlottingObjects.at(i), LegendLabel.at(i).Data(), LegendDrawOption(DrawOption.at(i)).Data());
  }

  //----------------------------------------------------------------------------
  //  The upper pad is now filled. Create and cd to the lower pad now
  Canvas->cd();
  RatioPad->Draw();
  RatioPad->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  RatioPad->SetLogx(logx);
  RatioPad->SetLogy(logz);
  if (logx)
    rDummy->GetXaxis()->SetMoreLogLabels();
  rDummy->Draw();

  //  Print all ratios, bot funcs and lines on the HistoPad
  //----------------------------------------------------------------------------
  for (int i = 0; i < (int)RatiObjects.size(); ++i)
  {
    RatiObjects.at(i)->Draw(Form("same %s", ((TString)DrawOptionR.at(i)).Data()));
    if ((int)*(LegendLabelR.at(i).Data()))
      legR->AddEntry(RatiObjects.at(i), LegendLabelR.at(i).Data(), LegendDrawOption(DrawOptionR.at(i)).Data());
  }

  //----------------------------------------------------------------------------
  //  Both pads are now filled. Create the white rectangle hiding the axis label conflict now
  Canvas->cd();
  if (wred)
    WhitePad->SetFillColor(kRed);
  WhitePad->Draw(); //  To eliminate the label conflict on y axis where plots meet

  Canvas->cd(); //  cd back to canvas in order to draw legend and Latex over entire canvas in relative coordinates
  Canvas->Update();

  if (leg->GetNRows() > 1)
  {
    if (leg->GetNRows() > 10)
      leg->SetTextSize(leg->GetTextSize() / 1.5);
    leg->Draw("same");
  }
  if (leg2->GetNRows() > 1)
  {
    if (leg2->GetNRows() > 10)
      leg2->SetTextSize(leg2->GetTextSize() / 1.5);
    leg2->Draw("same");
  }
  if (legR->GetNRows() > 1)
    legR->Draw("same");

  Canvas->SaveAs(name);
  delete hDummy;
  delete rDummy;
  delete Canvas;
}

void PlottingRatio::New(TH1 *h, TString label, int style, double size, int color, TString opt, bool InSecondLegend)
{
  if (!h)
    FATAL("New was given a invalid (nullptr) TH1.");

  PlottingObjects.push_back(h);
  SetStyle<TH1 *>(h, style, size, color, opt, counter);
  DrawOption.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  LegendLabel.push_back(Form("%s%s", InSecondLegend ? "42" : "", label.Data()));
  if ((style == -1) || (color == -1))
    counter++;
}

void PlottingRatio::New(TGraphAsymmErrors *g, TString label, int style, double size, int color, TString opt)
{
  if (!g)
    FATAL("New was given a invalid (nullptr) TGraphAsymmErrors.");

  PlottingObjects.push_back(g);
  SetStyle<TGraphAsymmErrors *>(g, style, size, color, opt, counter);
  DrawOption.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  LegendLabel.push_back(label);
  if ((style == -1) || (color == -1))
    counter++;
}

void PlottingRatio::NewRatio(TH1 *h, TString label, int style, double size, int color, TString opt)
{
  if (!h)
    FATAL("NewRatio was given a nullptr.");

  RatiObjects.push_back(h);
  SetStyle<TH1 *>(h, style, size, color, opt, counteR);
  DrawOptionR.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  LegendLabelR.push_back(label);
  if ((style == -1) || (color == -1))
    counteR++;
}

void PlottingRatio::NewRatio(TGraphAsymmErrors *g, TString label, int style, double size, int color, TString opt)
{
  if (!g)
    FATAL("NewRatio was given a nullptr.");

  RatiObjects.push_back(g);
  SetStyle<TGraphAsymmErrors *>(g, style, size, color, opt, counteR);
  DrawOptionR.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  LegendLabelR.push_back(label);
  if ((style == -1) || (color == -1))
    counteR++;
}

void PlottingRatio::New(TF1 *f, TString label, int style, double size, int color, TString opt)
{
  if (!f)
    FATAL("New was given a nullptr.");

  PlottingObjects.push_back(f);
  LegendLabel.push_back(label);
  DrawOption.push_back(opt);
  SetStyle<TF1 *>(f, style, size, color, opt, counter);

  if ((style == -1) || (color == -1))
    counter++;
}

void PlottingRatio::NewRatioLine(double x1, double y1, double x2, double y2, int style, int size, int color, TString label)
{
  TLine *line = new TLine(x1, y1, x2, y2);
  LegendLabelR.push_back(label);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(size);
  DrawOptionR.push_back("l");
  RatiObjects.push_back(line);
}

void PlottingRatio::NewRatio(TF1 *f, TString label, int style, double size, int color, TString opt)
{
  if (!f)
    FATAL("NewBotFunc was given a nullptr.");

  RatiObjects.push_back(f);

  SetStyle<TF1 *>(f, style, size, color, opt, counteR);
  DrawOptionR.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);
  LegendLabelR.push_back(label);

  if ((style == -1) || (color == -1))
    counteR++;
}

void PlottingRatio::SetAxisLabel(TString labelx, TString labely, TString labelz, double offsetx, double offsety, int ndigits)
{
  AxisLabel[0] = labelx;
  AxisLabel[1] = labely;
  AxisLabel[2] = labelz;
  AxisLabelOffset[0] = offsetx;
  AxisLabelOffset[1] = offsety;
  Ndigits = ndigits;
}

void PlottingRatio::InitializeAxis(bool logx, bool logy, bool TransCanvas)
{
  AutoSetAxisRanges(logy);

  hDummy = new TH2D("hDummy", "hDummy", 1000, AxisRange[0][0], AxisRange[0][1], 1000, AxisRange[1][0], AxisRange[1][1]);
  hDummy->SetTitle("");
  hDummy->SetStats(0);

  double Range[2] = {1E12, -1E12}; // Set the minima large enough and the maxima low enough, so that they will be overwritten in the first iteration of the following object loop

  for (int io = 0; io < (int)RatiObjects.size(); io++)
  {
    TString ClassName = (TString)RatiObjects.at(io)->ClassName();
    TH1F *h = nullptr;
    if (ClassName.Contains("TH1"))
    {
      h = (TH1F *)RatiObjects.at(io);
    }
    else if (ClassName.Contains("TGraph"))
    {
      h = ((TGraph *)RatiObjects.at(io))->GetHistogram();
    }
    else
      continue;
    Range[0] = (h->GetMinimum() < Range[0] && (!logy || h->GetMinimum() > 0)) ? h->GetMinimum() : Range[0];
    Range[1] = h->GetMaximum() > Range[1] ? h->GetMaximum() : Range[1];
  }

  Range[1] = logy ? 2 * Range[1] : Range[1] + (Range[1] - Range[0]) / 10;                                                                        //  With log scales a factor 2 isn't too much
  Range[0] = logy ? 0.5 * Range[0] : (Range[1] - 2 * Range[0] > 0 ? (Range[0] > 0 ? 0 : 1.1 * Range[0]) : Range[0] - (Range[1] - Range[0]) / 8); //  max-2*min>0 ? min is small->go down to 0 : all bins rather full
  for (int iE = 0; iE < 2; iE++)
  {
    if (AxisRange[2][iE] > 41.99 && AxisRange[2][iE] < 42.01)
      AxisRange[2][iE] = Range[iE];
  }

  rDummy = new TH2D("rDummy", "rDummy", 1000, AxisRange[0][0], AxisRange[0][1], 1000, AxisRange[2][0], AxisRange[2][1]);
  rDummy->SetTitle("");
  rDummy->SetStats(0);
  if (logx && AxisRange[0][0] > 0.5 && AxisRange[0][1] < 100)
    rDummy->SetLabelOffset(-0.025);

  TAxis *Axis[3] = {rDummy->GetXaxis(), hDummy->GetYaxis(), rDummy->GetYaxis()};

  for (int iA = 0; iA < 3; iA++)
  {
    Axis[iA]->SetRangeUser(AxisRange[iA][0], AxisRange[iA][1]);
    Axis[iA]->SetTitleSize(TitleSizeScaling * StdTextSize);
    Axis[iA]->SetLabelSize(StdTextSize);
    Axis[iA]->SetLabelFont(43);
    Axis[iA]->SetTitleFont(43);

    if (TransCanvas)
    {
      Axis[iA]->SetTitleColor(kWhite);
      Axis[iA]->SetLabelColor(kWhite);
    }
    Axis[iA]->SetTitle(AxisLabel[iA]);
    if (iA == 2)
      Axis[iA]->SetTitleOffset(1 + (AxisLabelOffset[1] - 1) / (1.5 - Split));
    else
      Axis[iA]->SetTitleOffset(AxisLabelOffset[iA]);
  }

  rDummy->GetXaxis()->SetTickSize(0.03 / Split);
  hDummy->GetXaxis()->SetTickSize(0.03 / (1 - Split));
  // rDummy->GetYaxis()->SetNdivisions(8);

  if (logx)
    hDummy->GetXaxis()->SetMoreLogLabels();
}

void PlottingRatio::SetWhite(double low, double left, double up, double right, bool red)
{
  WhiteBorders[0][0] = left;
  WhiteBorders[0][1] = right;
  WhiteBorders[1][0] = low;
  WhiteBorders[1][1] = up;
  wred = red;
}

void PlottingRatio::InitializeCanvas(bool logx, bool logy, bool logz, bool TransCanvas)
{

  if (Canvas)
    delete Canvas; //  This should never happen, but better safe than sorry.

  Canvas = new TCanvas("Canvas", "Canvas", CanvasDimensions[1], CanvasDimensions[0]);

  double PadSplit = (Split - Split * CanvasMargins[1][1] - Split * CanvasMargins[1][0]) + CanvasMargins[1][0];

  HistoPad = new TPad("HistoPad", "HistoPad", 0.0, PadSplit, 1, 1);
  RatioPad = new TPad("RatioPad", "RatioPad", 0.0, 0.0, 1, PadSplit);
  WhitePad = new TPad("WhitePad", "WhitePad", WhiteBorders[0][0], WhiteBorders[1][0], WhiteBorders[0][1], WhiteBorders[1][1]);
  if (TransCanvas)
  {
    Canvas->SetFillStyle(0);
    HistoPad->SetFillStyle(0);
    RatioPad->SetFillStyle(0);
    WhitePad->SetFillStyle(0);
  }

  HistoPad->SetTopMargin(CanvasMargins[1][1] / (1 - PadSplit));
  HistoPad->SetRightMargin(CanvasMargins[0][1]);
  HistoPad->SetLeftMargin(CanvasMargins[0][0]);
  HistoPad->SetBottomMargin(0);
  RatioPad->SetTopMargin(0);
  RatioPad->SetRightMargin(CanvasMargins[0][1]);
  RatioPad->SetLeftMargin(CanvasMargins[0][0]);
  RatioPad->SetBottomMargin(CanvasMargins[1][0] / PadSplit);

  Canvas->cd();
  HistoPad->Draw();
  HistoPad->cd();

  gPad->SetTickx();
  gPad->SetTicky();

  HistoPad->SetLogy(logy);
  HistoPad->SetLogx(logx);
}

void PlottingRatio::SetLegendR(double x1, double x2, double y1, double y2, bool bhollow)
{
  RatioLegendBorders[0][0] = x1;
  RatioLegendBorders[0][1] = x2;
  RatioLegendBorders[1][0] = y1;
  RatioLegendBorders[1][1] = y2;
  LRhollow = bhollow;
}

void PlottingRatio::InitializeLegendR()
{
  legR = new TLegend(RatioLegendBorders[0][0], RatioLegendBorders[1][0], RatioLegendBorders[0][1], RatioLegendBorders[1][1]);
  legR->SetHeader("");
  legR->SetTextFont(43);
  legR->SetTextSize(0.6 * StdTextSize);
  legR->SetBorderSize(0);
  legR->SetFillStyle(LRhollow ? 0 : 1001);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++ Plotting Grid ++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PlottingGrid : public Plotting
{
public:
  PlottingGrid(); // This empty constructor always has to be called when plotting Ratio

  ~PlottingGrid(); //  Destructor

  void Plot(TString name = "dummy.pdf", bool logx = false, bool logy = false);

  //  Add histograms to the upper pad
  void New(TH1F *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p");
  void New(TH1D *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p");
  void New(TGraphAsymmErrors *g = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "p");

  //  Adds a line to vector lines that will be drawn on canvas in function Plot(). The lines coordinates relate to the axis. When setting negative style a curly line will be drawn instead of the straight one.
  void NewLine(double x1 = 0, double y1 = 0, double x2 = 1, double y2 = 1, int style = 1, int color = kBlack, double width = 1, TString label = "");

  //  The ratios have a seperate legend. It's position can be set analogously to SetLegend using this function
  void SetGridLegend(double x1 = 0.7, double x2 = 0.95, double y1 = 0.15, double y2 = 0.25);

  //  Add functions to the tfuncs and bfuncs vector, which are added to the top and bottom pad when drawing in Plot()
  void NewFunc(TF1 *h = nullptr, TString label = "", int style = -1, double size = 1, int color = -1, TString opt = "l");

  void NextPad(TString Title, bool bFade = false);

  void NewPadTex(const double PositX = 0.2, const double PositY = 0.2, TString text = "", const double TextSize = 0.035, const double dDist = 0.05, const int font = 42, const int color = kBlack);

  void SetAxisLabel(TString labelx = "", TString labely = "", double offsetx = 0.9, double offsety = 1.);

protected:
  int ColumnsRows[31][2] = {{1, 1}, {1, 1}, {2, 1}, {2, 2}, {2, 2}, {3, 2}, {3, 2}, {4, 2}, {4, 2}, {3, 3}, {4, 3}, {4, 3}, {4, 3}, {5, 3}, {5, 3}, {5, 3}, {5, 4}, {5, 4}, {5, 4}, {5, 4}, {5, 4}, {6, 4}, {6, 4}, {6, 4}, {6, 4}, {6, 5}, {6, 5}, {6, 5}, {6, 5}, {6, 5}, {6, 5}};

  int Columns = 0;
  int Rows = 0;

  std::vector<TString> LegendLabelF;
  std::vector<TString> DrawOptionF;
  std::vector<TString> LegendLabelG;
  std::vector<TString> DrawOptionG;
  std::vector<TString> LegendLabelL;

  //  Add the pads that are drawn on the canvas in Plot() as attributes
  std::vector<TPad *> pads;
  std::vector<TPad *> Fadepads;

  int PadCounter = 1;

  TLegend *legR = nullptr;
  double GridLegendBorders[2][2] = {{0.7, 0.95}, {0.15, 0.25}};

  //  Vectors containing the functions for the upper pad as well as the functions and ratios for the lower pad
  std::vector<std::vector<TH1F>> hists;
  std::vector<TH1F> histsinpad; // Contains histograms in one pad
  std::vector<std::vector<TGraphAsymmErrors>> graphs;
  std::vector<TGraphAsymmErrors> graphsinpad; // Contains histograms in one pad
  std::vector<std::vector<TF1>> funcs;
  std::vector<TF1> funcsinpad; // Contains functions in one pad
  std::vector<std::vector<TLine>> lines;
  std::vector<TLine> linesinpad; // Contains lines in one pad
  std::vector<std::vector<TLatex>> padtex;
  std::vector<TLatex> texinpad;

  double FadeScale = 2. / 3.;
  std::vector<bool> bFadePad;

  //  Creates all three pads and the canvas that they are on. The canvas dimension attributes are NOT used, instead a standard size of 1000x1000 is used
  void InitializeCanvas(bool logx, bool logy);

  int NPads;
};

PlottingGrid::PlottingGrid()
{
  AxisLabelOffset[0] = 0.9;
}

PlottingGrid::~PlottingGrid()
{
}

void PlottingGrid::NextPad(TString Title, bool bFade)
{
  histsinpad.at(0).SetTitle(Title);

  if (histsinpad.size() > 0)
  {
    double max = histsinpad.at(0).GetMaximum(); //  Estimated upper y border

    //  For logarithmic y axis the minimum has be be larger than 0
    double min = histsinpad.at(0).GetMinimum();

    //  Find the smalles and largest bin in all loaded histograms (larger than 0 for logy)
    for (int i = 0; i < (int)histsinpad.size(); i++)
    {
      max = histsinpad.at(i).GetMaximum() > max ? histsinpad.at(i).GetMaximum() : max;
      min = histsinpad.at(i).GetMinimum() < min ? histsinpad.at(i).GetMinimum() : min;
    }

    max = max + (max - min) / 10;                                                  //  Leave space between highest/lowest bin and the axis borders so you can see every bin
    min = (max - 2 * min > 0 ? (min > 0 ? 0 : 1.1 * min) : min - (max - min) / 8); //  max-2*min>0 ? min is small->go down to 0 : all bins rather full

    if (AxisRange[1][0] == 42 && AxisRange[1][1] == 42)
      histsinpad.at(0).GetYaxis()->SetRangeUser(min, max);
    if (AxisRange[1][0] == 42 && AxisRange[1][1] != 42)
      histsinpad.at(0).GetYaxis()->SetRangeUser(min, AxisRange[1][1]);
    if (AxisRange[1][0] != 42 && AxisRange[1][1] == 42)
      histsinpad.at(0).GetYaxis()->SetRangeUser(AxisRange[1][0], max);
    if (AxisRange[1][0] != 42 && AxisRange[1][1] != 42)
      histsinpad.at(0).GetYaxis()->SetRangeUser(AxisRange[1][0], AxisRange[1][1]);
  }

  bFadePad.push_back(bFade);

  padtex.push_back(texinpad);
  texinpad.clear();
  hists.push_back(histsinpad);
  histsinpad.clear();
  graphs.push_back(graphsinpad);
  graphsinpad.clear();
  funcs.push_back(funcsinpad);
  funcsinpad.clear();
  lines.push_back(linesinpad);
  linesinpad.clear();
  PadCounter++;
  counter = 0;
}

void PlottingGrid::NewLine(double x1, double y1, double x2, double y2, int style, int color, double width, TString label)
{

  TLine *line = new TLine(x1, y1, x2, y2);
  if (PadCounter == 1)
    LegendLabelL.push_back(label);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(width);
  linesinpad.push_back(*line);
}

void PlottingGrid::SetAxisLabel(TString labelx, TString labely, double offsetx, double offsety)
{
  AxisLabel[0] = labelx;
  AxisLabel[1] = labely;
  AxisLabelOffset[0] = offsetx;
  AxisLabelOffset[1] = offsety;
}

void PlottingGrid::Plot(TString name, bool logx, bool logy)
{

  if (hists.size() < 1)
    FATAL("No hists added for plotting.");

  NPads = PadCounter; // Add one to the existing ones for the legend

  InitializeCanvas(logx, logy); // Creating Canvas with margins
  gStyle->SetEndErrorSize(0);

  hists.insert(hists.begin() + Columns - 1, histsinpad);
  graphs.insert(graphs.begin() + Columns - 1, graphsinpad);
  funcs.insert(funcs.begin() + Columns - 1, funcsinpad);
  lines.insert(lines.begin() + Columns - 1, linesinpad);
  padtex.insert(padtex.begin() + Columns - 1, texinpad);
  bFadePad.insert(bFadePad.begin() + Columns - 1, false);

  leg = new TLegend(LegendBorders[0][0], LegendBorders[1][0], LegendBorders[0][1], LegendBorders[1][1]);
  leg->SetHeader(""); //  Remove title of legend
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->SetBorderSize(0); //  Remove black rectangle around legend
  leg->SetFillStyle(0);  // Solid white background to make legend readable. Set to 0 to make it hollow

  for (int iPad = 0; iPad < (int)pads.size(); ++iPad)
  {

    Canvas->cd();
    pads.at(iPad)->Draw();
    pads.at(iPad)->cd();

    if (iPad == Columns - 1)
    {

      for (int i = 0; i < (int)PlottingObjects.size(); ++i)
        PlottingObjects.at(i)->Draw("same");

      for (int i = 0; i < (int)hists.at(0).size(); ++i)
      {
        if ((int)*(LegendLabel.at(i).Data()))
        {
          leg->AddEntry(&hists.at(0).at(i), LegendLabel.at(i).Data(), LegendDrawOption(DrawOption.at(i))); //  Dont add anything to the legend if LegendLabel is empty
        }
      }
      for (int i = 0; i < (int)graphs.at(0).size(); ++i)
      {
        if ((int)*(LegendLabelG.at(i).Data()))
          leg->AddEntry(&graphs.at(0).at(i), LegendLabelG.at(i).Data(), LegendDrawOption(DrawOptionG.at(i))); //  Dont add anything to the legend if LegendLabel is empty
      }
      for (int i = 0; i < (int)funcs.at(0).size(); ++i)
      {
        if ((int)*(LegendLabelF.at(i).Data()))
          leg->AddEntry(&funcs.at(0).at(i), LegendLabelF.at(i).Data(), LegendDrawOption(DrawOptionF.at(i))); //  Dont add anything to the legend if LegendLabel is empty
      }

      for (int i = 0; i < (int)lines.at(0).size(); ++i)
      {
        if ((int)*(LegendLabelL.at(i).Data()))
          leg->AddEntry(&lines.at(0).at(i), LegendLabelL.at(i).Data(), "l");
      }

      leg->Draw();

      continue;
    }

    gPad->SetTickx();
    gPad->SetTicky();
    for (int ihist = 0; ihist < (int)hists.at(iPad).size(); ++ihist)
      hists.at(iPad).at(ihist).Draw(Form("same %s", ((TString)DrawOption.at(ihist)).Data()));
    for (int igraph = 0; igraph < (int)graphs.at(iPad).size(); ++igraph)
      graphs.at(iPad).at(igraph).Draw(Form("same %s", ((TString)DrawOptionG.at(igraph)).Data()));
    for (int ifunc = 0; ifunc < (int)funcs.at(iPad).size(); ++ifunc)
      funcs.at(iPad).at(ifunc).Draw(Form("same %s", ((TString)DrawOptionF.at(ifunc)).Data()));
    for (int iline = 0; iline < (int)lines.at(iPad).size(); ++iline)
      lines.at(iPad).at(iline).Draw("same");
    for (int itex = 0; itex < (int)padtex.at(iPad).size(); ++itex)
      padtex.at(iPad).at(itex).Draw("same");

    if (bFadePad.at(iPad))
    {
      Fadepads.at(iPad)->SetFillColorAlpha(kWhite, FadeScale);
      Fadepads.at(iPad)->Draw();
    }
  }

  Canvas->SaveAs(name);
  delete Canvas;
}

void PlottingGrid::New(TH1F *h, TString label, int style, double size, int color, TString opt)
{

  if (!h)
    FATAL("NewHist was given a nullptr.");

  h->SetStats(0);
  h->GetYaxis()->SetTitleSize(TitleSizeScaling * h->GetYaxis()->GetLabelSize());
  h->GetXaxis()->SetTitleSize(TitleSizeScaling * h->GetXaxis()->GetLabelSize());
  h->GetXaxis()->SetTitle(AxisLabel[0]);
  h->GetYaxis()->SetTitle(AxisLabel[1]);
  h->GetXaxis()->SetTitleOffset(AxisLabelOffset[0]);
  h->GetYaxis()->SetTitleOffset(AxisLabelOffset[1]);
  h->GetYaxis()->SetMaxDigits(3);

  if (PadCounter == 1)
    DrawOption.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);

  if (PadCounter == 1)
    LegendLabel.push_back(label);

  SetStyle<TH1 *>(h, style, size, color, opt, counter);
  h->SetLineWidth(1);
  histsinpad.push_back(*h);
  counter++;
}

void PlottingGrid::New(TGraphAsymmErrors *g, TString label, int style, double size, int color, TString opt)
{

  if (!g)
    FATAL("NewHist was given a nullptr.");

  g->SetStats(0);
  g->GetYaxis()->SetTitleSize(TitleSizeScaling * g->GetYaxis()->GetLabelSize());
  g->GetXaxis()->SetTitleSize(TitleSizeScaling * g->GetXaxis()->GetLabelSize());
  g->GetXaxis()->SetTitle(AxisLabel[0]);
  g->GetYaxis()->SetTitle(AxisLabel[1]);
  g->GetXaxis()->SetTitleOffset(AxisLabelOffset[0]);
  g->GetYaxis()->SetTitleOffset(AxisLabelOffset[1]);
  g->GetYaxis()->SetMaxDigits(3);

  if (PadCounter == 1)
    DrawOptionG.push_back((opt == "l" || opt == "c") ? opt + " hist" : opt);

  if (PadCounter == 1)
    LegendLabelG.push_back(label);

  SetStyle<TGraphAsymmErrors *>(g, style, size, color, opt, counter);
  g->SetLineWidth(1);
  counter++;
  graphsinpad.push_back(*g);
}

void PlottingGrid::New(TH1D *h, TString label, int style, double size, int color, TString opt)
{
  TH1F *hd = (TH1F *)h;
  New(hd, label, style, size, color, opt);
}

void PlottingGrid::NewPadTex(const double PositX, const double PositY, TString text, const double TextSize, const double dDist, const int font, const int color)
{
  std::vector<TString> LatStr;             //  Each element corresponds to a line of the printed latex string
  TObjArray *textStr = text.Tokenize(";"); //  The semicolon seperates the string into different lines
  for (int i = 0; i < textStr->GetEntries(); i++)
  {
    TObjString *tempObj = (TObjString *)textStr->At(i);
    LatStr.push_back(tempObj->GetString());
  }

  //  Loop thru the latex lines and set the formatting
  for (int i = 0; i < (int)LatStr.size(); ++i)
  {
    TLatex *textmp = new TLatex(PositX, PositY - i * dDist, LatStr[i]);
    texinpad.push_back(*textmp);
    texinpad.at(texinpad.size() - 1).SetNDC();
    texinpad.at(texinpad.size() - 1).SetTextFont(font);
    texinpad.at(texinpad.size() - 1).SetTextColor(color);
    texinpad.at(texinpad.size() - 1).SetTextSize(TextSize);
  }
}

void PlottingGrid::NewFunc(TF1 *f, TString label, int style, double size, int color, TString opt)
{

  if (!f)
    FATAL("NewTopFunc was given a nullptr.");

  if (PadCounter == 1)
    LegendLabelF.push_back(label);
  if (PadCounter == 1)
    DrawOptionF.push_back(opt);

  f->SetMarkerStyle((style == -1) ? AutoStyle[counter] : style);
  f->SetLineStyle((style == -1) ? AutoStyleLine[counter] : style);
  f->SetMarkerColor((color == -1) ? AutoColor[counter] : color);
  f->SetLineColor((color == -1) ? AutoColor[counter] : color);
  f->SetMarkerSize(size);
  f->SetLineWidth(size);
  funcsinpad.push_back(*f);

  counter++;
}

void PlottingGrid::InitializeCanvas(bool logx, bool logy)
{

  if (Canvas)
    delete Canvas; //  This should never happen, but better safe than sorry.

  Columns = ColumnsRows[NPads][0];
  Rows = ColumnsRows[NPads][1];

  Canvas = new TCanvas("Canvas", "Canvas", 1600, 900);

  int Column = 0;     // Left to right
  int Row = Rows - 1; // Up to down
  for (int iPad = 0; iPad < NPads; iPad++)
  {
    TPad *Pad = new TPad(Form("Pad_%d", iPad), Form("Pad_%d", iPad), ((double)Column) / ((double)Columns), ((double)Row) / ((double)Rows), ((double)Column + 1) / ((double)Columns), ((double)Row + 1) / ((double)Rows));
    Pad->SetTopMargin(0.09);
    Pad->SetRightMargin(0.01);
    Pad->SetLeftMargin(0.09);
    Pad->SetBottomMargin(0.092);
    Pad->SetLogx(logx);
    Pad->SetLogy(logy);

    pads.push_back(Pad);

    TPad *FadePad = new TPad(Form("FadePad_%d", iPad), Form("FadePad_%d", iPad), 0, 0, 1, 1);
    Fadepads.push_back(FadePad);

    Column++;
    if (Column == Columns)
    {
      Column = 0;
      Row--;
      if (Row < 0)
        break;
    }
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++ Plotting Paint +++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PlottingPaint : public Plotting
{

public:
  //  Using this class one can draw ellipses, angles, lines and curly lines
  void Plot(TString name = "You_forgot_the_name_..._dummy.pdf");

  void SetCanvas(double cw = 1200, double ch = 1200);

  //  Draw an incomplete ellipse, emulating an angle for example to draw a particle decay angle
  void NewAngle(double x = 0.5, double y = 0.5, double r1 = 0.3, double r2 = 0.2,
                double phimin = 60, double phimax = 120, double theta = 0);

protected:
  void InitializeCanvas();
};

void PlottingPaint::Plot(TString name)
{

  InitializeCanvas(); //  Creating Canvas with margins

  for (int i = 0; i < (int)PlottingObjects.size(); ++i)
    PlottingObjects.at(i)->Draw("same");

  Canvas->SaveAs(name);
  delete Canvas;
}

void PlottingPaint::NewAngle(double x, double y, double r1, double r2, double phimin, double phimax, double theta)
{

  TEllipse *angle = new TEllipse(x, y, r1, r2, phimin, phimax, theta);
  angle->SetNoEdges();
  PlottingObjects.push_back(angle);
}

void PlottingPaint::InitializeCanvas()
{

  if (Canvas)
    delete Canvas; //  This should never happen, but better safe than sorry.

  Canvas = new TCanvas("Canvas", "Canvas", CanvasDimensions[0], CanvasDimensions[1]);
  Canvas->cd();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++ Extra functions ++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TString ReturnTimehms(int s)
{
  int h = s / 3600;
  s = s % 3600;
  int m = s / 60;
  s = s % 60;
  return Form("%d hours, %d minutes and %d seconds", h, m, s);
}

void PrintProgress(int i, int N, int Steps = 100)
{

  int barWidth = 100;

  static int Progress = 0;

  static clock_t t = clock();

  if (((double)i) / ((double)N) * Steps > Progress)
  {

    cout.flush();
    cout << " [" << Progress * 100 / Steps << "%]"
         << "[";
    int pos = barWidth * (((double)Progress) / Steps);
    for (int i = 0; i < barWidth; ++i)
    {
      if (i < pos)
        cout << "|";
      else
        cout << " ";
    }
    cout << "] - " << ReturnTimehms(((clock() - t) / CLOCKS_PER_SEC) * ((double)(N - i) / (double)i)) << " left.  "
         << "\r";

    Progress++;
  }
  if (i == N - 1)
  {
    cout.flush();
    cout << "[" << Progress << "%]"
         << "[";
    int pos = barWidth * (Progress / 100.);
    for (int i = 0; i < barWidth; ++i)
    {
      cout << "|";
    }
    cout << "] - "
         << "Processing finished in " << ReturnTimehms((clock() - t) / CLOCKS_PER_SEC) << endl;
  }
}

void PrintProgressNumber(int i, int N, int Steps = 100)
{
  static int Progress = 0;
  if (((double)i) / ((double)N) * Steps > Progress)
  {
    cout << "[" << Form("%.1f", Progress * 100. / Steps) << "%]" << endl;
    Progress++;
  }
}

#endif
