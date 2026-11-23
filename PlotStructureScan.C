// Figures of paper Section 15 from the output of RunStructureScan.C.
// Plotting code taken from the original scripts (PlotSummary.C,
// PlotTemplateGradient.C, PlotTemplateGradientX.C); only file names changed.
// Usage: root -l -b -q 'PlotStructureScan.C("inDir")'  -> PNGs with the paper's names in inDir
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TH2D.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

void PlotSummary(const char* filename = "CentralFraction_vs_Rin.txt", const char* yTitle = "Central fraction [%]", const char* outName = nullptr)
{
    std::ifstream f(filename);
    if (!f.is_open()) { printf("Cannot open %s\n", filename); return; }

    std::string paramLine, headerLine;
    std::getline(f, paramLine);
    std::getline(f, headerLine);

    // Config names = header tokens after the first ("Rin[Ohm]")
    std::vector<std::string> configs;
    {
        std::istringstream hs(headerLine);
        std::string tok;
        std::getline(hs, tok, '\t');   // discard "Rin[Ohm]"
        while (std::getline(hs, tok, '\t')) configs.push_back(tok);
    }
    int nConfig = (int)configs.size();
    if (nConfig == 0) { printf("No configuration columns found in %s\n", filename); return; }

    std::vector<TGraph*> graphs(nConfig);
    for (int i = 0; i < nConfig; i++) graphs[i] = new TGraph();

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        std::string tok;
        std::vector<double> v;
        while (std::getline(ls, tok, '\t'))
            if (!tok.empty()) { char* e; double d = strtod(tok.c_str(), &e); if (e != tok.c_str()) v.push_back(d); }
        if ((int)v.size() < 1) continue;
        double rin = v[0];
        for (int i = 0; i < nConfig && i + 1 < (int)v.size(); i++)
            graphs[i]->SetPoint(graphs[i]->GetN(), rin, v[i + 1]);
    }
    f.close();

    gStyle->SetOptStat(0);

    const int   kCol[5] = { kBlack, kRed+1, kAzure+1, kGreen+2, kMagenta+1 };
    const int   kMrk[5] = { 20, 21, 22, 23, 33 };

    TMultiGraph *mg = new TMultiGraph(Form("mgCF_%s", filename), "");
    TLegend *leg = new TLegend(0.65, 0.15, 0.90, 0.40);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);

    for (int i = 0; i < nConfig; i++) {
        TGraph *g = graphs[i];
        g->SetMarkerStyle(kMrk[i % 5]); g->SetMarkerSize(1.3);
        g->SetMarkerColor(kCol[i % 5]); g->SetLineColor(kCol[i % 5]);
        g->SetLineWidth(2); g->SetLineStyle(2);
        mg->Add(g, "PL");
        leg->AddEntry(g, configs[i].c_str(), "p");
    }

    TCanvas *c1 = new TCanvas(Form("c1_%s", filename), "c1 - Central-4 fraction vs Rin", 800, 600);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.12);

    mg->SetTitle(TString::Format(";Rin [#Omega];%s", yTitle));
    mg->Draw("A");
    mg->GetYaxis()->SetRangeUser(0, 105);
    mg->GetXaxis()->SetTitleSize(0.045); mg->GetYaxis()->SetTitleSize(0.045);
    leg->Draw();

    // Split the parameter line in two at the nearest space to its midpoint,
    // so it fits within the canvas width instead of running off the edge.
    std::string line1 = paramLine, line2;
    if (paramLine.size() > 60) {
        size_t mid = paramLine.size() / 2;
        size_t sp = paramLine.rfind(' ', mid);
        if (sp == std::string::npos) sp = paramLine.find(' ', mid);
        if (sp != std::string::npos) {
            line1 = paramLine.substr(0, sp);
            line2 = paramLine.substr(sp + 1);
        }
    }
    TLatex ltx; ltx.SetNDC(); ltx.SetTextSize(0.024);
    ltx.DrawLatex(0.13, 0.965, line1.c_str());
    if (!line2.empty()) ltx.DrawLatex(0.13, 0.935, line2.c_str());

    TString outPng = filename;
    if (outPng.EndsWith(".txt")) outPng.Resize(outPng.Length() - 4);
    outPng += ".png";
    if (outName) outPng = outName;
    c1->Update();
    c1->SaveAs(outPng.Data());
    printf("Plot saved to %s\n", outPng.Data());
}

// Shared by PlotUnbalance() and PlotTemplateFit(): reads a "sigma/nu vs Rin"
// file and plots both quantities side by side in one canvas.
//
// File format:
//   line 1: space-separated "key=value" fixed simulation parameters
//   line 2: header, tab-separated: Rin[Ohm]  <cfg1>_sigma  <cfg1>_nu  <cfg2>_sigma  <cfg2>_nu  ...
//   line 3+: data rows, tab-separated: Rin  sigma1 nu1  sigma2 nu2  ...
void PlotSigmaNu(const char* filename, const char* canvasTitle, const char* outPng,
                  const char* suffix1 = "_sigma", const char* suffix2 = "_nu",
                  const char* yTitle1 = "#sigma_{t-Student} [#mum]", const char* yTitle2 = "#nu_{t-Student}",
                  double yMax1 = -1, double yMax2 = -1, const char* xTitle = "Rin [#Omega]",
                  double legX1 = 0.18, double legY1 = 0.684)
{
    std::ifstream f(filename);
    if (!f.is_open()) { printf("Cannot open %s\n", filename); return; }

    std::string paramLine, headerLine;
    std::getline(f, paramLine);
    std::getline(f, headerLine);

    // Config names = header tokens after the first, taking every other
    // token (sigma column) and stripping the "_sigma" suffix.
    std::vector<std::string> configs;
    {
        std::istringstream hs(headerLine);
        std::string tok;
        std::getline(hs, tok, '\t');   // discard "Rin[Ohm]"
        int col = 0;
        while (std::getline(hs, tok, '\t')) {
            if (col % 2 == 0) {
                size_t p = tok.rfind(suffix1);
                configs.push_back(p != std::string::npos ? tok.substr(0, p) : tok);
            }
            col++;
        }
    }
    int nConfig = (int)configs.size();
    if (nConfig == 0) { printf("No configuration columns found in %s\n", filename); return; }

    std::vector<TGraph*> gSigma(nConfig), gNu(nConfig);
    for (int i = 0; i < nConfig; i++) { gSigma[i] = new TGraph(); gNu[i] = new TGraph(); }

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        std::string tok;
        std::vector<double> v;
        while (std::getline(ls, tok, '\t'))
            if (!tok.empty()) { char* e; double d = strtod(tok.c_str(), &e); if (e != tok.c_str()) v.push_back(d); }
        if ((int)v.size() < 1) continue;
        double rin = v[0];
        for (int i = 0; i < nConfig; i++) {
            int iSigma = 1 + 2*i, iNu = 2 + 2*i;
            if (iSigma < (int)v.size()) gSigma[i]->SetPoint(gSigma[i]->GetN(), rin, v[iSigma]);
            if (iNu    < (int)v.size()) gNu[i]->SetPoint(gNu[i]->GetN(), rin, v[iNu]);
        }
    }
    f.close();

    gStyle->SetOptStat(0);

    const int kCol[5] = { kBlack, kRed+1, kAzure+1, kGreen+2, kMagenta+1 };
    const int kMrk[5] = { 20, 21, 22, 23, 33 };

    auto buildMG = [&](std::vector<TGraph*> &graphs, const char* name) -> TMultiGraph* {
        TMultiGraph *mg = new TMultiGraph(name, "");
        for (int i = 0; i < nConfig; i++) {
            TGraph *g = graphs[i];
            g->SetMarkerStyle(kMrk[i % 5]); g->SetMarkerSize(1.3);
            g->SetMarkerColor(kCol[i % 5]); g->SetLineColor(kCol[i % 5]);
            g->SetLineWidth(2); g->SetLineStyle(2);
            mg->Add(g, "PL");
        }
        return mg;
    };

    TMultiGraph *mgSigma = buildMG(gSigma, Form("mgS_%s", outPng));
    TMultiGraph *mgNu    = buildMG(gNu, Form("mgN_%s", outPng));

    // Common y-range for both pads, computed directly from the point values
    // (not from a drawn histogram) so margins/titles only need to be set once,
    // before Draw -- setting pad margins after Draw() leaves the already-built
    // frame geometry stale and clips the axis labels.
    auto maxOf = [](std::vector<TGraph*> &graphs) {
        double m = 0;
        for (auto *g : graphs) {
            double *y = g->GetY();
            for (int i = 0; i < g->GetN(); i++) m = std::max(m, y[i]);
        }
        return m;
    };
    double autoMax = std::max(maxOf(gSigma), maxOf(gNu)) * 1.35;
    double useMax1 = (yMax1 > 0) ? yMax1 : autoMax;
    double useMax2 = (yMax2 > 0) ? yMax2 : autoMax;

    TCanvas *c2 = new TCanvas(Form("c2_%s", outPng), canvasTitle, 1300, 600);
    c2->Divide(2, 1);

    c2->cd(1);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.12);
    mgSigma->SetTitle(TString::Format(";%s;%s", xTitle, yTitle1));
    mgSigma->Draw("A");
    mgSigma->GetYaxis()->SetTitleOffset(1.3);
    mgSigma->GetYaxis()->CenterTitle(true);
    mgSigma->SetMinimum(0);
    mgSigma->SetMaximum(useMax1);
    mgSigma->GetXaxis()->SetTitleSize(0.045); mgSigma->GetYaxis()->SetTitleSize(0.045);
    TLegend *legS = new TLegend(legX1, legY1, legX1+0.39, legY1+0.14);
    legS->SetBorderSize(0); legS->SetFillStyle(0); legS->SetTextSize(0.032);
    for (int i = 0; i < nConfig; i++) legS->AddEntry(gSigma[i], configs[i].c_str(), "p");
    legS->Draw();

    c2->cd(2);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.12);
    mgNu->SetTitle(TString::Format(";%s;%s", xTitle, yTitle2));
    mgNu->Draw("A");
    mgNu->GetYaxis()->SetTitleOffset(1.3);
    mgNu->GetYaxis()->CenterTitle(true);
    mgNu->SetMinimum(0);
    mgNu->SetMaximum(useMax2);
    mgNu->GetXaxis()->SetTitleSize(0.045); mgNu->GetYaxis()->SetTitleSize(0.045);
    TLegend *legN = new TLegend(legX1, legY1, legX1+0.39, legY1+0.14);
    legN->SetBorderSize(0); legN->SetFillStyle(0); legN->SetTextSize(0.032);
    for (int i = 0; i < nConfig; i++) legN->AddEntry(gNu[i], configs[i].c_str(), "p");
    legN->Draw();

    // Parameter annotation split in two lines, drawn on the canvas itself
    // (not a sub-pad) so it spans the full width above both plots.
    c2->cd();
    std::string line1 = paramLine, line2;
    if (paramLine.size() > 100) {
        size_t mid = paramLine.size() / 2;
        size_t sp = paramLine.rfind(' ', mid);
        if (sp == std::string::npos) sp = paramLine.find(' ', mid);
        if (sp != std::string::npos) { line1 = paramLine.substr(0, sp); line2 = paramLine.substr(sp + 1); }
    }
    TLatex ltx; ltx.SetNDC(); ltx.SetTextSize(0.018);
    ltx.DrawLatex(0.02, 0.975, line1.c_str());
    if (!line2.empty()) ltx.DrawLatex(0.02, 0.955, line2.c_str());

    c2->Update();
    c2->SaveAs(outPng);
    printf("Plot saved to %s\n", outPng);
}


static bool LoadTemplate(const char* path, int &nx, int &ny, double &pp, TH2D* &hVal)
{
    std::ifstream f(path);
    if (!f.is_open()) { printf("Cannot open %s\n", path); return false; }
    f >> nx >> ny >> pp;
    hVal = new TH2D(Form("val_%s", path), "", nx, 0, pp, ny, 0, pp);
    int bx, by; double v[4];
    while (f >> bx >> by >> v[0] >> v[1] >> v[2] >> v[3])
        hVal->SetBinContent(bx, by, v[0]);   // v[0] = pad7 = top-left electrode
    f.close();
    return true;
}

// Gradient magnitude via central differences (one-sided at the edges),
// expressed in "fraction per mm" for readability.
static TH2D* BuildGradientMap(TH2D* hVal, int nx, int ny, double pp, const char* name)
{
    double binw_mm = (pp / nx) / 1000.0;   // bin width in mm
    TH2D *hGrad = new TH2D(name, "", nx, 0, pp, ny, 0, pp);
    for (int bx = 1; bx <= nx; bx++) {
        for (int by = 1; by <= ny; by++) {
            double vxm = hVal->GetBinContent(std::max(1, bx-1), by);
            double vxp = hVal->GetBinContent(std::min(nx, bx+1), by);
            double vym = hVal->GetBinContent(bx, std::max(1, by-1));
            double vyp = hVal->GetBinContent(bx, std::min(ny, by+1));
            int nxsteps = std::min(nx, bx+1) - std::max(1, bx-1);
            int nysteps = std::min(ny, by+1) - std::max(1, by-1);
            double gx = (nxsteps > 0) ? (vxp - vxm) / (nxsteps * binw_mm) : 0;
            double gy = (nysteps > 0) ? (vyp - vym) / (nysteps * binw_mm) : 0;
            hGrad->SetBinContent(bx, by, std::sqrt(gx*gx + gy*gy));
        }
    }
    return hGrad;
}


static bool LoadLeftSum(const char* path, int &nx, int &ny, double &pp, TH2D* &hSum)
{
    std::ifstream f(path);
    if (!f.is_open()) { printf("Cannot open %s\n", path); return false; }
    f >> nx >> ny >> pp;
    hSum = new TH2D(Form("sum_%s", path), "", nx, 0, pp, ny, 0, pp);
    int bx, by; double v[4];
    while (f >> bx >> by >> v[0] >> v[1] >> v[2] >> v[3])
        hSum->SetBinContent(bx, by, v[0] + v[2]);   // pad7 (top-left) + pad6 (bottom-left)
    f.close();
    return true;
}

// d/dx only (central differences, one-sided at the edges), in "per mm".
static TH2D* BuildGradX(TH2D* hSum, int nx, int ny, double pp, const char* name)
{
    double binw_mm = (pp / nx) / 1000.0;
    TH2D *hGx = new TH2D(name, "", nx, 0, pp, ny, 0, pp);
    for (int bx = 1; bx <= nx; bx++) {
        for (int by = 1; by <= ny; by++) {
            double vxm = hSum->GetBinContent(std::max(1, bx-1), by);
            double vxp = hSum->GetBinContent(std::min(nx, bx+1), by);
            int nsteps = std::min(nx, bx+1) - std::max(1, bx-1);
            double gx = (nsteps > 0) ? (vxp - vxm) / (nsteps * binw_mm) : 0;
            hGx->SetBinContent(bx, by, gx);
        }
    }
    return hGx;
}


void PlotStructureScan(const char *d = "sec15_out") {
    gStyle->SetOptStat(0);
    PlotSummary(Form("%s/CentralFraction_vs_Rin.txt", d), "Central fraction [%]", Form("%s/centralfraction_vs_rin.png", d));
    PlotSummary(Form("%s/SideFraction_vs_Rin.txt", d), "Central fraction [%]", Form("%s/sidefraction_vs_rin.png", d));
    PlotSigmaNu(Form("%s/UnbalanceFit_vs_Rin.txt", d), "Unbalance reconstruction", Form("%s/unbalance_comparison.png", d),
                "_sigma", "_nu", "#sigma_{t-Student} [#mum]", "#nu_{t-Student}", 80, 30, "Rin [#Omega]", 0.234, 0.5472);
    PlotSigmaNu(Form("%s/TemplateFit_vs_Rin.txt", d), "Template reconstruction", Form("%s/templatefit_vs_rin.png", d),
                "_sigma", "_nu", "#sigma_{t-Student} [#mum]", "#nu_{t-Student}", 80, 30);
    PlotSigmaNu(Form("%s/InterpadResistance_vs_Rin.txt", d), "Inter-pad resistance", Form("%s/interpadresistance_vs_rin.png", d),
                "_Rside", "_Rdiag", "R_{side} [#Omega]", "R_{diag} [#Omega]");

    const char *lab[3] = {"trench", "no_structure", "RStrip1000"};
    TH2D *hVal[3], *hSum[3], *hGx[3];
    for (int i = 0; i < 3; i++) {
        int nx, ny; double pp;
        TString path = Form("%s/template_%s_Rin50.txt", d, lab[i]);
        if (!LoadTemplate(path.Data(), nx, ny, pp, hVal[i])) return;
        if (!LoadLeftSum(path.Data(), nx, ny, pp, hSum[i])) return;
        hGx[i] = BuildGradX(hSum[i], nx, ny, pp, Form("gx_%d", i));
    }
    TCanvas *c1 = new TCanvas("cT1", "pad7 fraction", 1500, 500); c1->Divide(3, 1);
    for (int i = 0; i < 3; i++) {
        c1->cd(i+1); gPad->SetRightMargin(0.15);
        hVal[i]->SetTitle(TString::Format("%s: pad7 fraction;x [#mum];y [#mum]", lab[i]));
        hVal[i]->SetMinimum(0); hVal[i]->SetMaximum(1.0); hVal[i]->Draw("COLZ");
    }
    c1->SaveAs(Form("%s/template_pad7_comparison.png", d));
    TCanvas *c2 = new TCanvas("cT2", "left sum", 1500, 500); c2->Divide(3, 1);
    for (int i = 0; i < 3; i++) {
        c2->cd(i+1); gPad->SetRightMargin(0.15);
        hSum[i]->SetTitle(TString::Format("%s: frac_{7}+frac_{6} (left);x [#mum];y [#mum]", lab[i]));
        hSum[i]->SetMinimum(0); hSum[i]->SetMaximum(1.0); hSum[i]->Draw("COLZ");
    }
    c2->SaveAs(Form("%s/leftsum_comparison.png", d));
    TCanvas *c3 = new TCanvas("cT3", "d(left sum)/dx", 1500, 500); c3->Divide(3, 1);
    for (int i = 0; i < 3; i++) {
        c3->cd(i+1); gPad->SetRightMargin(0.15);
        hGx[i]->SetTitle(TString::Format("%s: d(frac_{7}+frac_{6})/dx [1/mm];x [#mum];y [#mum]", lab[i]));
        hGx[i]->SetMinimum(-4); hGx[i]->SetMaximum(0); hGx[i]->Draw("COLZ");
        double s = 0; int nb = hGx[i]->GetNbinsX()*hGx[i]->GetNbinsY();
        for (int bx = 1; bx <= hGx[i]->GetNbinsX(); bx++) for (int by = 1; by <= hGx[i]->GetNbinsY(); by++) s += hGx[i]->GetBinContent(bx, by);
        printf("%-14s mean(dSum/dx)=%.4f /mm\n", lab[i], s/nb);
    }
    c3->SaveAs(Form("%s/leftsum_gradient_comparison.png", d));
}
