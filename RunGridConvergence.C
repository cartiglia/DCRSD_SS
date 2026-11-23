// Grid-convergence study (referee comment B1).
// Rsheet = 1000 Ohm/sq, R_in = 50 Ohm, pad radius 30 um, h = p/ndiv.
//
// part = 0: sharing, all structures (no structure, trench, 500 Ohm strip),
//           ndiv = 50, 100 (paper default), 200, compared with ndiv = 200:
//           - f_left = (P6+P7)/(P6+P7+P10+P11) along x at mid-height and its mean
//             gradient (sets the position resolution);
//           - central fraction (P6+P7+P10+P11)/total, averaged over a 5x5 grid of
//             hits in the central pixel (signal leakage).
// part = 1: inter-electrode resistances R_side, R_diag (Section 12), all
//           structures, ndiv = 50, 100, 200.
// part = 2: trench of fixed physical width (one node pitch of the default grid,
//           p/100): sharing with ndiv = 100 and 300 (3 rows), and the effect of
//           the missing gain layer under the trench (ndiv = 100).
// Usage: root -l -q 'RunGridConvergence.C(part)'   (the GUI must be created: no -b)
#include "DCRSD_SurfaceSimulator.c"

struct GridResult { std::vector<double> fl; double grad, cen; };

GridResult Sharing(RSDStripPanel *p, double pitch) {
    GridResult r;
    const int nx = 10;
    std::vector<double> xs(nx);
    r.fl.assign(nx, 0.0);
    for (int i = 0; i < nx; i++) {
        xs[i] = 0.05 + 0.1*i;
        std::vector<double> f;
        p->StaticHit(pitch*(1 + xs[i]), 1.5*pitch, f);
        r.fl[i] = (f[6]+f[7])/(f[6]+f[7]+f[10]+f[11]);
    }
    double sg = 0;
    for (int i = 0; i < nx; i++) {
        int a = std::max(0, i-1), b = std::min(nx-1, i+1);
        sg += std::abs((r.fl[b]-r.fl[a])/((xs[b]-xs[a])*pitch*1e-3));
    }
    r.grad = sg/nx;
    double sc = 0;
    for (int ix = 0; ix < 5; ix++) for (int iy = 0; iy < 5; iy++) {
        std::vector<double> f;
        p->StaticHit(pitch*(1.1+0.2*ix), pitch*(1.1+0.2*iy), f);
        sc += f[6]+f[7]+f[10]+f[11];
    }
    r.cen = sc/25.;
    return r;
}

void PrintSharing(const char *tag, double pitch, const char *sname, const char *lab,
                  const GridResult &r, const GridResult &ref) {
    double dmax = 0;
    for (size_t i = 0; i < ref.fl.size(); i++) dmax = std::max(dmax, std::abs(r.fl[i]-ref.fl[i]));
    printf("%s pitch=%.0f struct=\"%s\" %s  grad=%.4f (%.4f)  max|df_left|=%.4f  central=%.4f (%+.4f)\n",
           tag, pitch, sname, lab, r.grad, r.grad/ref.grad, dmax, r.cen, r.cen-ref.cen);
}

void RunGridConvergence(int part = 0, double rsheet = 1000, double rin = 50) {
    RSDStripPanel *p = new RSDStripPanel();
    const char *sname[3] = {"no structure", "trench", "strip 500 Ohm"};
    if (part == 0 || part == 1) {
        const int nd[3] = {50, 100, 200};
        for (double pitch : {300., 500., 1000.})
            for (int st : {0, 1, 2}) {
                std::vector<GridResult> res; double rs[3], rd[3];
                for (int g = 0; g < 3; g++) {
                    p->SetGridDivisions(nd[g]);
                    p->ConfigureForCheck(pitch, st, rsheet, 50, rin);
                    if (part == 0) res.push_back(Sharing(p, pitch));
                    else { p->DoResistanceMatrix(); rs[g] = p->LastRside(); rd[g] = p->LastRdiag(); }
                }
                for (int g = 0; g < 3; g++) {
                    TString lab = TString::Format("ndiv=%d h=%.1f", nd[g], pitch/nd[g]);
                    if (part == 0) PrintSharing("CONV", pitch, sname[st], lab.Data(), res[g], res[2]);
                    else printf("RES pitch=%.0f struct=\"%s\" %s  Rside=%.1f (%.4f)  Rdiag=%.1f (%.4f)  ratio=%.3f\n",
                                pitch, sname[st], lab.Data(), rs[g], rs[g]/rs[2], rd[g], rd[g]/rd[2], rd[g]/rs[g]);
                }
                fflush(stdout);
            }
    }
    if (part == 2) {
        for (double pitch : {300., 500., 1000.}) {
            p->SetGridDivisions(300);
            p->ConfigureForCheck(pitch, 1, rsheet, 50, rin);
            GridResult ref = Sharing(p, pitch);
            p->SetGridDivisions(100);
            p->ConfigureForCheck(pitch, 1, rsheet, 50, rin);
            GridResult r100 = Sharing(p, pitch);
            p->SetTrenchNoGain(false);
            GridResult r100g = Sharing(p, pitch);
            p->SetTrenchNoGain(true);
            PrintSharing("TRENCH", pitch, "trench", "ndiv=300 (ref)          ", ref, ref);
            PrintSharing("TRENCH", pitch, "trench", "ndiv=100                 ", r100, ref);
            PrintSharing("TRENCH", pitch, "trench", "ndiv=100 gain under trench", r100g, r100);
            fflush(stdout);
        }
    }
}
