// Domain-size convergence study (referee item B2).
//
// The simulator models a fixed array of 3x3 pixels (4x4 electrodes), with the outer
// electrodes centred on the domain boundary and no current across the boundary
// (Neumann). This macro re-implements the same static network model for an array
// of NP x NP pixels ((NP+1)^2 electrodes, NP odd) and evaluates, for the central
// pixel, the quantities used in the paper:
//   - fraction of the signal collected by the four electrodes of the central pixel,
//     for a hit at the pixel centre and near a pixel side (local (10, 250) um);
//   - fraction collected by the two left electrodes, normalised to the four
//     electrodes of the pixel, along the horizontal mid-line, and its gradient
//     d(f_left)/dx at the pixel centre and 25 um from the pixel edge;
//   - inter-electrode resistance (through the electrode disks, others floating)
//     between central side and diagonal neighbours.
// Model (identical to DCRSD_SurfaceSimulator.c): square grid h = p/nd; electrode =
// disk of radius a, nodes shorted (1e6 S), only the centre node drains via R_in
// (R_in = 0: all electrode nodes at V = 0); resistive strip = band |d| <= 0.025 p
// along the electrode rows/columns (boundary ones included), strip sheet
// resistance R_strip,total/20 in parallel to the surface; trench = node row on the
// inner pixel boundaries, edges across it cut except within a of an electrode,
// no gain (no injected charge) on trench nodes; bilinear charge injection.
// With NP = 3 the results coincide with those of the simulator.
//
// Usage: root -l -b -q 'RunDomainSize.C(nd)'   (no GUI needed)
#include <TMatrixDSparse.h>
#include <TDecompSparse.h>
#include <TVectorD.h>
#include <TMath.h>
#include <vector>
#include <cstdio>

namespace dom {

class Solver : public TDecompSparse {
public: void MaxPivoting() { SetThresholdPivoting(kThresholdPivotingMax); }
};

struct Geo {
    int NP, nd, n;          // pixels per side, nodes per pitch, nodes per side
    double pp, h, a, rsh, rstrip, rin;
    int structure;          // 0 none, 1 trench, 2 strip
    std::vector<int> type;  // 0 surface, k>0 electrode id
    int Pad(int r, int c) const { return r*(NP+1) + c + 1; }
    int Node(int i, int j) const { return i*n + j; }
    bool TrenchLine(int k) const { return k % nd == 0 && k > 0 && k < n-1; }   // width p/100 at nd = 100: one row
    bool NearPad(double along) const {
        for (int c = 0; c <= NP; c++) if (std::abs(along - c*pp) <= a) return true;
        return false;
    }
    bool TrenchNode(int i, int j) const {
        if (structure != 1) return false;
        return (TrenchLine(j) && !NearPad(i*h)) || (TrenchLine(i) && !NearPad(j*h));
    }
    double G(int i, int j, int di, int dj) const {
        if (structure == 1) {
            if (TrenchLine(j) != TrenchLine(j+dj) && !NearPad(i*h)) return 0;
            if (TrenchLine(i) != TrenchLine(i+di) && !NearPad(j*h)) return 0;
            return 1.0/rsh;
        }
        if (structure == 2) {
            double xm = (i + 0.5*di)*h, ym = (j + 0.5*dj)*h;
            double d = (dj == 0) ? ym : xm;   // horizontal edge: horizontal strips at y = k p
            for (int k = 0; k <= NP; k++)
                if (std::abs(d - k*pp) <= pp*0.025) return 1.0/rsh + 1.0/rstrip;
        }
        return 1.0/rsh;
    }
    void Build(int NP_, int nd_, double pp_, double a_, double rsh_, int st, double rstripTot, double rin_) {
        NP = NP_; nd = nd_; pp = pp_; a = a_; rsh = rsh_; structure = st;
        rstrip = rstripTot/20.0; rin = rin_; n = NP*nd + 1; h = pp/nd;
        type.assign(n*n, 0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            for (int r = 0; r <= NP; r++) for (int c = 0; c <= NP; c++)
                if (TMath::Hypot(i*h - r*pp, j*h - c*pp) <= a) type[Node(i,j)] = Pad(r,c);
    }
    int CentreNode(int pad) const {
        int r = (pad-1)/(NP+1), c = (pad-1)%(NP+1);
        return Node(r*nd, c*nd);
    }
};

const int di[4] = {1,-1,0,0}, dj[4] = {0,0,1,-1};
const double kGshort = 1e6;

// free = true: matrix for the resistance (no R_in, no Dirichlet), weak reference at node 0.
void Matrix(const Geo &g, bool free, TMatrixDSparse &M) {
    int N = g.n*g.n;
    std::vector<int> ir, ic; std::vector<double> v;
    ir.reserve(5*N); ic.reserve(5*N); v.reserve(5*N);
    bool useRin = g.rin > 0;
    for (int i = 0; i < g.n; i++) for (int j = 0; j < g.n; j++) {
        int idx = g.Node(i,j), pk = g.type[idx];
        if (!free && pk > 0 && !useRin) { ir.push_back(idx); ic.push_back(idx); v.push_back(1.0); continue; }
        double sum = 0;
        if (!free && pk > 0 && idx == g.CentreNode(pk)) sum += 1.0/g.rin;
        if (free && idx == 0) sum += 1e-9;
        for (int k = 0; k < 4; k++) {
            int ni = i+di[k], nj = j+dj[k];
            if (ni < 0 || ni >= g.n || nj < 0 || nj >= g.n) continue;
            int nidx = g.Node(ni,nj), npk = g.type[nidx];
            double c;
            if (pk > 0 && npk > 0) { if (pk != npk) continue; c = kGshort; }
            else c = g.G(i, j, di[k], dj[k]);
            sum += c;
            bool dir = (!free && !useRin && npk > 0);   // neighbour at V = 0
            if (!dir) { ir.push_back(idx); ic.push_back(nidx); v.push_back(-c); }
        }
        ir.push_back(idx); ic.push_back(idx); v.push_back(sum);
    }
    M.ResizeTo(N, N);
    M.SetMatrixArray(v.size(), ir.data(), ic.data(), v.data());
}

// Currents per electrode for a unit charge at (x,y); returns the injected weight.
double Hit(const Geo &g, Solver &s, double x, double y, std::vector<double> &cur) {
    int N = g.n*g.n, npad = (g.NP+1)*(g.NP+1);
    cur.assign(npad+1, 0.0);
    TVectorD B(N); B.Zero();
    double gx = x/g.h, gy = y/g.h;
    int i0 = std::max(0, std::min(g.n-2, (int)std::floor(gx)));
    int j0 = std::max(0, std::min(g.n-2, (int)std::floor(gy)));
    double fx = gx - i0, fy = gy - j0, win = 0;
    int ii[4] = {i0, i0+1, i0, i0+1}, jj[4] = {j0, j0, j0+1, j0+1};
    double w[4] = {(1-fx)*(1-fy), fx*(1-fy), (1-fx)*fy, fx*fy};
    for (int k = 0; k < 4; k++) {
        if (g.TrenchNode(ii[k], jj[k])) continue;
        int idx = g.Node(ii[k], jj[k]); win += w[k];
        if (g.rin <= 0 && g.type[idx] > 0) cur[g.type[idx]] += w[k]; else B(idx) += w[k];
    }
    if (win <= 0) return 0;
    bool ok; TVectorD V = s.Solve(B, ok);
    if (g.rin > 0) {
        for (int p = 1; p <= npad; p++) cur[p] += V(g.CentreNode(p))/g.rin;
    } else {
        for (int i = 0; i < g.n; i++) for (int j = 0; j < g.n; j++) {
            int idx = g.Node(i,j); if (g.type[idx]) continue;
            for (int k = 0; k < 4; k++) {
                int ni = i+di[k], nj = j+dj[k];
                if (ni < 0 || ni >= g.n || nj < 0 || nj >= g.n) continue;
                int npk = g.type[g.Node(ni,nj)];
                if (npk) cur[npk] += V(idx)*g.G(i, j, di[k], dj[k]);
            }
        }
    }
    return win;
}

void Run(int NP, int nd, int st, double rstripTot, double rin, const char *lab) {
    const double pp = 500, a = 30, rsh = 1000;
    Geo g; g.Build(NP, nd, pp, a, rsh, st, rstripTot, rin);
    int c0 = (NP-1)/2;                                   // central pixel: electrodes c0, c0+1
    double X0 = c0*pp, Y0 = c0*pp;                       // lower-left corner of the central pixel
    int pLL = g.Pad(c0,c0), pLU = g.Pad(c0,c0+1), pRL = g.Pad(c0+1,c0), pRU = g.Pad(c0+1,c0+1);

    TMatrixDSparse M; Matrix(g, false, M);
    Solver s; s.SetMatrix(M); s.MaxPivoting(); s.Decompose();
    std::vector<double> cur;
    auto frac = [&](double lx, double ly, double &fcen, double &fleft) {
        Hit(g, s, X0+lx, Y0+ly, cur);
        double tot = 0; for (size_t k = 1; k < cur.size(); k++) tot += cur[k];
        double c4 = cur[pLL]+cur[pLU]+cur[pRL]+cur[pRU];
        fcen = 100*c4/tot; fleft = (cur[pLL]+cur[pLU])/c4;
    };
    double fc, fs, fl, dummy, fl1, fl2;
    frac(250, 250, fc, fl);
    frac(10, 250, fs, dummy);
    frac(240, 250, dummy, fl1); frac(260, 250, dummy, fl2);
    double gc = (fl2 - fl1)/0.020;                       // per mm
    frac(15, 250, dummy, fl1); frac(35, 250, dummy, fl2);
    double ge = (fl2 - fl1)/0.020;
    double fl100; frac(100, 250, dummy, fl100);
    double fl100c; frac(100, 100, dummy, fl100c);

    TMatrixDSparse Mf; Matrix(g, true, Mf);
    Solver sf; sf.SetMatrix(Mf); sf.MaxPivoting(); sf.Decompose();
    auto R = [&](int pa, int pb) {
        TVectorD B(g.n*g.n); B.Zero(); B(g.CentreNode(pa)) = 1; B(g.CentreNode(pb)) = -1;
        bool ok; TVectorD V = sf.Solve(B, ok); return V(g.CentreNode(pa)) - V(g.CentreNode(pb));
    };
    double Rs = R(pRU, pLU), Rd = R(pRU, pLL);
    printf("DOM %-13s NP=%d nd=%d Rin=%3.0f  central=%.2f%% side=%.2f%%  fleft(250,250)=%.4f fleft(100,250)=%.4f fleft(100,100)=%.4f  "
           "grad_centre=%.3f grad_edge=%.3f /mm  Rside=%.1f Rdiag=%.1f\n",
           lab, NP, nd, rin, fc, fs, fl, fl100, fl100c, gc, ge, Rs, Rd);
    fflush(stdout);
}

} // namespace dom

void RunDomainSize(int nd = 100, int npMax = 7) {
    struct C { const char *lab; int st; double rs; } cf[4] =
        { {"no_structure", 0, 0}, {"trench", 1, 0}, {"strip500", 2, 500}, {"strip2000", 2, 2000} };
    for (int NP = 3; NP <= npMax; NP += 2)
        for (auto &c : cf)
            for (double rin : {0., 50.})
                dom::Run(NP, nd, c.st, c.rs, rin, c.lab);
}
