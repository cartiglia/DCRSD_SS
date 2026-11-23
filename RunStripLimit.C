// Limit check: resistive strip with increasing total pad-to-pad resistance must
// tend to the no-structure case. Pitch 500 um, Rsheet 1000 Ohm/sq, pad radius 30 um.
// For each R_strip (total, Ohm): central fraction at the pixel centre for
// R_in = 0 and 50 Ohm, and the geometric inter-electrode resistance.
// Usage: root -l -q RunStripLimit.C < /dev/null   (GUI must be created: no -b)
#include "DCRSD_SurfaceSimulator.c"

void RunStripLimit() {
    RSDStripPanel *p = new RSDStripPanel();
    const double rs[] = {500, 2000, 10000, 20000, 40000, 100000, 400000, 2000000};
    const int n = sizeof(rs)/sizeof(rs[0]);
    std::vector<double> f;
    for (int k = -1; k < n; k++) {
        int st = (k < 0) ? 0 : 2;
        double r = (k < 0) ? 500 : rs[k];
        double cen[2];
        for (int ir = 0; ir < 2; ir++) {
            p->ConfigureForCheck(500, st, 1000, 50, ir ? 50 : 0, r);
            p->StaticHit(750, 750, f);
            cen[ir] = 100*(f[6]+f[7]+f[10]+f[11]);
        }
        p->DoResistanceMatrix();
        printf("LIMIT %-14s central(Rin0)=%.2f%% central(Rin50)=%.2f%% R_side=%.1f R_diag=%.1f\n",
               st ? TString::Format("Rstrip=%.0f", r).Data() : "no_structure",
               cen[0], cen[1], p->LastRside(), p->LastRdiag());
        fflush(stdout);
    }
    gApplication->Terminate(0);
}
