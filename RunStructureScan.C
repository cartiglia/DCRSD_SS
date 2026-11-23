// Comparison of the inter-electrode structures vs R_in (paper Section 15).
//
// Five configurations: trench, no structure, resistive strip with total pad-to-pad
// resistance 500, 1000, 2000 Ohm (strip width p/20 = 20 squares: 25, 50, 100 Ohm/sq). Pitch 500 um, Rsheet = 1000 Ohm/sq,
// pad radius 30 um, R_in = 0, 25, 50, 100 Ohm, default grid (h = p/100).
// For each configuration and R_in:
//   - central fraction: charge fraction on the 4 central electrodes for a hit at
//     the pixel centre and near the pixel side (local (10, 250) um);
//   - unbalance reconstruction: sigma and nu of a Student's t fit to the x residual;
//   - charge-sharing templates (50 x 50 bins) and template (chi2) reconstruction:
//     sigma and nu of a Student's t fit to the x residual;
// and for each configuration the inter-electrode resistance (through the
// electrode disks, others floating), R_eff = R_geometric + 2 R_in.
// Reconstruction: amplitude 50 a.u., noise 1 a.u., no Landau, 5000 events
// uniformly in the central pixel (fixed seed). Hits on an electrode, or fully
// inside the trench (no gain), are skipped.
//
// Output (in outDir): CentralFraction_vs_Rin.txt, SideFraction_vs_Rin.txt,
// UnbalanceFit_vs_Rin.txt, TemplateFit_vs_Rin.txt, InterpadResistance_vs_Rin.txt,
// templates template_<config>_Rin<R>.txt (format: "nx ny PP", then
// "bx by f7 f11 f6 f10"), read by PlotStructureScan.C.
// Usage: root -l -q 'RunStructureScan.C("outDir")'   (the GUI must be created: no -b)
#include "DCRSD_SurfaceSimulator.c"
#include <TH1D.h>
#include <TF1.h>
#include <TRandom3.h>

namespace sscan {

const double kPP = 500, kRsheet = 1000, kAmp = 50, kNoise = 1;
const int    kNEvents = 5000, kNx = 50;
const int    pIDs[4] = {7, 11, 6, 10};

struct Cfg { const char *label; int structure; double rstrip; };  // rstrip = total pad-to-pad strip resistance [Ohm]

bool OnPad(RSDStripPanel *p, double x, double y) {
    for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++)
        if (TMath::Hypot(x - r*kPP, y - c*kPP) <= 30.0) return true;
    return false;
}

void FitStudent(TH1D &h, double &sigma, double &nu, bool templFit) {
    double xlo = h.GetXaxis()->GetXmin(), xhi = h.GetXaxis()->GetXmax();
    if (!templFit) {   // as in the unbalance scan
        TF1 f("fSt_u", "[0]*TMath::Gamma(([3]+1.0)/2.0)/(TMath::Sqrt([3]*TMath::Pi())*TMath::Gamma([3]/2.0)*[2])"
                       "*pow(1.0+1.0/[3]*pow((x-[1])/[2],2),-([3]+1.0)/2.0)", xlo, xhi);
        double rms = h.GetRMS();
        f.SetParameters(h.Integral("width"), h.GetBinCenter(h.GetMaximumBin()), rms, 2.0);
        f.SetParLimits(2, 0.1, rms*5); f.SetParLimits(3, 1.0, 100.0);
        h.Fit(&f, "RQ0");
        sigma = f.GetParameter(2); nu = f.GetParameter(3);
    } else {           // as in the template scan
        h.Fit("gaus", "Q0");
        TF1 *g = h.GetFunction("gaus");
        double sig0 = g ? g->GetParameter(2) : kPP/20.0;
        TF1 f("fSt_t", "[0]/[2] * TMath::Student((x-[1])/[2], [3])", xlo, xhi);
        f.SetParameters(h.Integral("width"), 0.0, sig0, 3.0);
        f.SetParLimits(2, 1e-3, kPP); f.SetParLimits(3, 0.5, 200.0);
        h.Fit(&f, "RQ0");
        sigma = f.GetParameter(2); nu = f.GetParameter(3);
    }
}

// Unbalance reconstruction (Eqs. unbalance_x): x from the left/right sums.
void Unbalance(RSDStripPanel *p, double &sigma, double &nu) {
    TRandom3 rand(4357);
    TH1D h("hUnb", "", 60, -150, 150);
    for (int ev = 0; ev < kNEvents; ev++) {
        double x = kPP + rand.Rndm()*kPP, y = kPP + rand.Rndm()*kPP;
        if (OnPad(p, x, y)) continue;
        std::vector<double> cur;
        if (p->StaticCurrents(x, y, cur) <= 0) continue;
        double A6 = cur[6]*kAmp, A7 = cur[7]*kAmp, A10 = cur[10]*kAmp, A11 = cur[11]*kAmp;
        double s4nl = A6+A7+A10+A11;
        double xr_nl = (s4nl > 1e-9) ? kPP*1.5 + ((A10+A11)-(A6+A7))/s4nl*kPP/2.0 : kPP*1.5;
        double a6 = A6 + rand.Gaus(0, kNoise), a7 = A7 + rand.Gaus(0, kNoise);
        double a10 = A10 + rand.Gaus(0, kNoise), a11 = A11 + rand.Gaus(0, kNoise);
        double s4 = a6+a7+a10+a11;
        double xr = (s4 > 1e-9) ? kPP*1.5 + ((a10+a11)-(a6+a7))/s4*kPP/2.0 : xr_nl;
        h.Fill(xr - x);
    }
    FitStudent(h, sigma, nu, false);
}

void BuildTemplates(RSDStripPanel *p, TH2D *t[4]) {
    double step = kPP/kNx;
    for (int q = 0; q < 4; q++) { t[q] = new TH2D(Form("tmpl_%d_%d", q, gRandom->Integer(1000000)), "", kNx, 0, kPP, kNx, 0, kPP); t[q]->SetDirectory(0); }
    for (int bx = 0; bx < kNx; bx++) for (int by = 0; by < kNx; by++) {
        double x = kPP + (bx + 0.5)*step, y = kPP + (by + 0.5)*step;
        std::vector<double> cur;
        p->StaticCurrents(x, y, cur);
        double s4 = 0; for (int q = 0; q < 4; q++) s4 += cur[pIDs[q]];
        for (int q = 0; q < 4; q++) t[q]->SetBinContent(bx+1, by+1, (s4 > 0) ? cur[pIDs[q]]/s4 : 0);
    }
}

void SaveTemplates(TH2D *t[4], const TString &path) {
    std::ofstream f(path.Data());
    f << kNx << " " << kNx << " " << kPP << "\n";
    for (int bx = 1; bx <= kNx; bx++) for (int by = 1; by <= kNx; by++) {
        f << bx << " " << by;
        for (int q = 0; q < 4; q++) f << " " << t[q]->GetBinContent(bx, by);
        f << "\n";
    }
}

// Template (chi2) reconstruction, as RSDStripPanel::DoTemplateReconstruction.
void TemplateReco(RSDStripPanel *p, TH2D *t[4], double &sigma, double &nu) {
    TRandom3 rand(9871);
    TH1D h("hTpl", "", 100, -200., 200.);
    for (int ev = 0; ev < kNEvents; ev++) {
        double x = kPP + rand.Rndm()*kPP, y = kPP + rand.Rndm()*kPP;
        if (OnPad(p, x, y)) continue;
        std::vector<double> cur;
        if (p->StaticCurrents(x, y, cur) <= 0) continue;
        double A[4], sum = 0;
        for (int q = 0; q < 4; q++) { A[q] = kAmp*cur[pIDs[q]] + rand.Gaus(0, kNoise); sum += A[q]; }
        for (int q = 0; q < 4; q++) A[q] /= sum;
        auto chi2 = [&](int bx, int by) {
            double c = 0;
            for (int q = 0; q < 4; q++) { double d = A[q] - t[q]->GetBinContent(bx, by); c += d*d; }
            return c;
        };
        double c0 = 1e9; int bxb = 1, byb = 1;
        for (int bx = 1; bx <= kNx; bx++) for (int by = 1; by <= kNx; by++) {
            double c = chi2(bx, by);
            if (c < c0) { c0 = c; bxb = bx; byb = by; }
        }
        double cm = (bxb > 1) ? chi2(bxb-1, byb) : c0, cp = (bxb < kNx) ? chi2(bxb+1, byb) : c0;
        double den = cm - 2*c0 + cp, dx = (den > 0) ? 0.5*(cm - cp)/den : 0;
        dx = std::max(-0.5, std::min(0.5, dx));
        double xr = t[0]->GetXaxis()->GetBinCenter(bxb) + dx*t[0]->GetXaxis()->GetBinWidth(bxb);
        h.Fill(xr - (x - kPP));
    }
    FitStudent(h, sigma, nu, true);
}

} // namespace sscan

void RunStructureScan(const char *outDir = "sec15_out") {
    using namespace sscan;
    gSystem->mkdir(outDir, kTRUE);
    RSDStripPanel *p = new RSDStripPanel();
    Cfg cfgs[5] = { {"trench", 1, 0}, {"no_structure", 0, 0},
                    {"RStrip500", 2, 500}, {"RStrip1000", 2, 1000}, {"RStrip2000", 2, 2000} };
    const double rins[4] = {0, 25, 50, 100};
    double cen[5][4], side[5][4], uS[5][4], uN[5][4], tS[5][4], tN[5][4], rS[5], rD[5];

    for (int ic = 0; ic < 5; ic++) {
        for (int ir = 0; ir < 4; ir++) {
            p->ConfigureForCheck(kPP, cfgs[ic].structure, kRsheet, 50, rins[ir], cfgs[ic].rstrip);
            std::vector<double> f;
            p->StaticHit(1.5*kPP, 1.5*kPP, f);       cen[ic][ir]  = 100*(f[6]+f[7]+f[10]+f[11]);
            p->StaticHit(kPP + 10, kPP + 250, f);     side[ic][ir] = 100*(f[6]+f[7]+f[10]+f[11]);
            Unbalance(p, uS[ic][ir], uN[ic][ir]);
            TH2D *t[4]; BuildTemplates(p, t);
            SaveTemplates(t, TString::Format("%s/template_%s_Rin%.0f.txt", outDir, cfgs[ic].label, rins[ir]));
            TemplateReco(p, t, tS[ic][ir], tN[ic][ir]);
            for (int q = 0; q < 4; q++) delete t[q];
            printf("SCAN %-13s Rin=%3.0f  central=%.2f%%  side=%.2f%%  unbalance sigma=%.1f nu=%.2f  template sigma=%.1f nu=%.2f\n",
                   cfgs[ic].label, rins[ir], cen[ic][ir], side[ic][ir], uS[ic][ir], uN[ic][ir], tS[ic][ir], tN[ic][ir]);
            fflush(stdout);
        }
        p->DoResistanceMatrix();
        rS[ic] = p->LastRside(); rD[ic] = p->LastRdiag();
        printf("SCAN %-13s R_side=%.1f R_diag=%.1f Ohm (geometric)\n", cfgs[ic].label, rS[ic], rD[ic]);
        fflush(stdout);
    }

    TString par = "Pitch=500um PadRadius=30um Rsheet=1000Ohm/sq StripTotalRes=500/1000/2000Ohm(StripSheetRes=25/50/100Ohm/sq,20sq) TrenchWidth=5um_no_gain";
    auto header = [&](std::ofstream &o, const char *extra, bool pair, const char *s1, const char *s2) {
        o << par << " " << extra << "\nRin[Ohm]";
        for (int ic = 0; ic < 5; ic++) {
            if (pair) o << "\t" << cfgs[ic].label << s1 << "\t" << cfgs[ic].label << s2;
            else      o << "\t" << cfgs[ic].label;
        }
        o << "\n";
    };
    {
        std::ofstream o(TString::Format("%s/CentralFraction_vs_Rin.txt", outDir).Data());
        header(o, "EvalPoint=pixel_center", false, "", "");
        for (int ir = 0; ir < 4; ir++) { o << rins[ir]; for (int ic = 0; ic < 5; ic++) o << "\t" << cen[ic][ir]; o << "\n"; }
    }
    {
        std::ofstream o(TString::Format("%s/SideFraction_vs_Rin.txt", outDir).Data());
        header(o, "EvalPoint=pixel_side_local(10,250)", false, "", "");
        for (int ir = 0; ir < 4; ir++) { o << rins[ir]; for (int ic = 0; ic < 5; ic++) o << "\t" << side[ic][ir]; o << "\n"; }
    }
    {
        std::ofstream o(TString::Format("%s/UnbalanceFit_vs_Rin.txt", outDir).Data());
        header(o, "Amplitude=50 Noise=1 NEvents=5000", true, "_sigma", "_nu");
        for (int ir = 0; ir < 4; ir++) { o << rins[ir]; for (int ic = 0; ic < 5; ic++) o << "\t" << uS[ic][ir] << "\t" << uN[ic][ir]; o << "\n"; }
    }
    {
        std::ofstream o(TString::Format("%s/TemplateFit_vs_Rin.txt", outDir).Data());
        header(o, "Amplitude=50 Noise=1 NEvents=5000 TemplateBins=50", true, "_sigma", "_nu");
        for (int ir = 0; ir < 4; ir++) { o << rins[ir]; for (int ic = 0; ic < 5; ic++) o << "\t" << tS[ic][ir] << "\t" << tN[ic][ir]; o << "\n"; }
    }
    {
        std::ofstream o(TString::Format("%s/InterpadResistance_vs_Rin.txt", outDir).Data());
        header(o, "Method=through_electrode_disks_others_floating R_eff=R_geometric+2*Rin", true, "_Rside", "_Rdiag");
        for (int ir = 0; ir < 4; ir++) { o << rins[ir]; for (int ic = 0; ic < 5; ic++) o << "\t" << rS[ic] + 2*rins[ir] << "\t" << rD[ic] + 2*rins[ir]; o << "\n"; }
    }
    printf(">>> Done: results in %s/\n", outDir);
}
