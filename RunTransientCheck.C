// Batch check of the quasi-static approximation (paper Section 8).
// For pitch 300/500/1000 um, trench and no structure: static vs transient sharing
// (signal area and peak amplitude). Pulse and thickness from the GUI defaults
// (trapezoid 400/100/500 ps, d = 50 um); Rsheet = 1000 Ohm/sq, R_in = 50 Ohm.
// Usage: root -l -q RunTransientCheck.C
#include "DCRSD_SurfaceSimulator.c"

void RunTransientCheck(double rsheet = 1000, double thick = 50, double rin = 50) {
    RSDStripPanel *p = new RSDStripPanel();
    for (double pitch : {300., 500., 1000.})
        for (int structure : {0, 1}) {
            p->ConfigureForCheck(pitch, structure, rsheet, thick, rin);
            p->DoTransientSharing();
        }
}
