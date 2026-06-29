#include <TMatrixDSparse.h>
#include <TVectorD.h>
#include <TDecompSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGFrame.h>
#include <TGNumberEntry.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TApplication.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TRandom3.h> 
#include <TF1.h> 
#include <TEllipse.h>
#include <TLine.h>
#include <TMarker.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <queue> 

// Physical Constants
const double eps0 = 8.854e-12; // F/m
const double epsSi = 11.7;
// pitch_m removed: use Pitch()*1e-6 instead

class RSDStripPanel : public TGMainFrame {
private:
    TRootEmbeddedCanvas *fEcanvas; 
    
    // Inputs
    TGNumberEntry       *fNumX, *fNumY, *fNumRsheet, *fNumThick, *fNumDur;
    TGNumberEntry       *fNumRStripPad, *fNumRStripMid, *fNumStripGap, *fNumStripGapPad, *fNumPixelPitch;
    TGNumberEntry       *fNumAmplitude;
    TGNumberEntry       *fNumNoise, *fNumAmpOffset, *fNumConstTerm, *fNumScanStep, *fNumThreshold, *fNumRecoEvents;
    TGNumberEntry       *fNumRin;
    TGNumberEntry       *fNumPadRadius;
    TGTextEntry         *fCfgFile;
    TGCheckButton       *fChkTrench    = nullptr;
    TGCheckButton       *fChkNoStruct = nullptr;
    TGCheckButton       *fChkLandau = nullptr;
    TGLabel             *fLblStripR;
    TGTextButton        *fStatusLabel;
    TGTextButton        *fBtnReco   = nullptr;
    TGTextButton        *fBtnExport = nullptr;

    const int n = 301;
    double PP()        { return fNumPixelPitch->GetNumber(); }   // pixel Pitch() [um]
    double Pitch()     { return PP() / 100.0; }                  // grid Pitch() [um/node]
    double Dim()       { return 3.0 * PP(); }                    // total detector size [um]
    double PadRadius()  { return fNumPadRadius->GetNumber(); }   // pad radius [um]
    double Amplitude()  { return fNumAmplitude->GetNumber(); }   // signal amplitude [a.u.]
    const int nTot = n * n;

    // --- CACHE ---
    TDecompSparse  *fSolverStatic = nullptr;
    TMatrixDSparse *fMatrixStatic = nullptr; 
    double fLastRsheetStat = -1, fLastRStripPStat = -1, fLastRStripMStat = -1, fLastGapStat = -1, fLastGapPadStat = -1, fLastPPStat = -1, fLastRinStat = -1, fLastPadRadiusStat = -1;
    bool   fLastTrenchStat = false;

    TDecompSparse  *fSolverTrans = nullptr;
    TMatrixDSparse *fMatrixTrans = nullptr; 
    double fLastRsheetTrans = -1, fLastRStripPTrans = -1, fLastRStripMTrans = -1, fLastThickTrans = -1, fLastGapTrans = -1, fLastPPTrans = -1, fLastRinTrans = -1;
    bool   fLastTrenchTrans = false;
    double fLastDtTrans = -1, fLastCnodeTrans = -1;

    std::vector<int> fTypeMap;
    std::vector<int> fPadOwnerMap;

    TH2D *fTemplates[4] = {};
    bool fTemplatesReady = false;
    TCanvas *cRmat = nullptr;
    std::vector<TGTextButton*> fDepButtons; // enabled only after first Run

public:
 RSDStripPanel() : TGMainFrame(gClient->GetRoot(), 1650, 900), fTypeMap(301*301, 0), fPadOwnerMap(301*301, 0), fTemplatesReady(false) {
        SetCleanup(kDeepCleanup);
        gStyle->SetPalette(kBird);
        gStyle->SetOptStat(0);

        // Build type maps with default pitch (fNumPixelPitch not yet initialized)
        BuildTypeMaps(500.0);

        // --- GUI CONSTRUCTION ---
        TGHorizontalFrame *hFrame = new TGHorizontalFrame(this);
        TGHorizontalFrame *controlFrame = new TGHorizontalFrame(hFrame, 1020, 900);
        TGVerticalFrame *vCol1 = new TGVerticalFrame(controlFrame, 340, 900);
        TGVerticalFrame *vCol2 = new TGVerticalFrame(controlFrame, 340, 900);
        TGVerticalFrame *vCol3 = new TGVerticalFrame(controlFrame, 320, 900);
        
        TGLayoutHints *lblL = new TGLayoutHints(kLHintsLeft | kLHintsTop, 20, 5, 10, 0);
        TGLayoutHints *entL = new TGLayoutHints(kLHintsLeft | kLHintsTop, 20, 5, 2, 5);
        TGLayoutHints *btnL = new TGLayoutHints(kLHintsLeft | kLHintsTop, 20, 5, 10, 0);
        TGLayoutHints *headL = new TGLayoutHints(kLHintsCenterX | kLHintsTop, 0, 0, 15, 5);
        auto BoldHead = [&](TGCompositeFrame *f, const char *txt) {
            TGLabel *l = new TGLabel(f, txt);
            l->SetTextFont("-*-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
            f->AddFrame(l, headL);
        };

        BoldHead(vCol1, "--- Injection Position ---");
        vCol1->AddFrame(new TGLabel(vCol1, "Position X [um]:"), lblL);
        fNumX = new TGNumberEntry(vCol1, 750, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 1500);
        vCol1->AddFrame(fNumX, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "Position Y [um]:"), lblL);
        fNumY = new TGNumberEntry(vCol1, 750, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 1500);
        vCol1->AddFrame(fNumY, entL);

        BoldHead(vCol1, "--- Geometry ---");
        vCol1->AddFrame(new TGLabel(vCol1, "Pixel Pitch [um]:"), lblL);
        fNumPixelPitch = new TGNumberEntry(vCol1, 500, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 50, 1500);
        vCol1->AddFrame(fNumPixelPitch, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "Thickness [um]:"), lblL);
        fNumThick = new TGNumberEntry(vCol1, 50, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 10, 1000);
        vCol1->AddFrame(fNumThick, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "Pad Radius [um]:"), lblL);
        fNumPadRadius = new TGNumberEntry(vCol1, 60, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 1, 500);
        vCol1->AddFrame(fNumPadRadius, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "# bins per pixel:"), lblL);
        fNumScanStep = new TGNumberEntry(vCol1, 100, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 5, 200);
        vCol1->AddFrame(fNumScanStep, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "Surface Resistivity [Ohm/sq]:"), lblL);
        fNumRsheet = new TGNumberEntry(vCol1, 1000, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 1, 100000);
        vCol1->AddFrame(fNumRsheet, entL);

        BoldHead(vCol1, "--- Analysis Parameters ---");
        vCol1->AddFrame(new TGLabel(vCol1, "Signal Amplitude [a.u.]:"), lblL);
        fNumAmplitude = new TGNumberEntry(vCol1, 100, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 0.1, 1e6);
        vCol1->AddFrame(fNumAmplitude, entL);

        fChkLandau = new TGCheckButton(vCol1, "Landau distribution");
        vCol1->AddFrame(fChkLandau, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "Noise RMS [a.u.]:"), lblL);
        fNumNoise = new TGNumberEntry(vCol1, 1.0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0.0, 10.0);
        vCol1->AddFrame(fNumNoise, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "Amplitude offset [a.u.]:"), lblL);
        fNumAmpOffset = new TGNumberEntry(vCol1, 0.0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0.0, 50.0);
        vCol1->AddFrame(fNumAmpOffset, entL);

        vCol1->AddFrame(new TGLabel(vCol1, "Constant term [um] (ref smear):"), lblL);
        fNumConstTerm = new TGNumberEntry(vCol1, 0.0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0.0, 500.0);
        vCol1->AddFrame(fNumConstTerm, entL);

        fNumPixelPitch->Connect("ValueSet(Long_t)", "RSDStripPanel", this, "UpdateStripRLabel()");

        BoldHead(vCol2, "--- Inter pixel ---");
        fChkTrench = new TGCheckButton(vCol2, "Trench between pixels");
        vCol2->AddFrame(fChkTrench, new TGLayoutHints(kLHintsLeft | kLHintsTop, 20, 5, 10, 5));
        fChkTrench->Connect("Toggled(Bool_t)", "RSDStripPanel", this, "DoTrenchToggle()");

        fChkNoStruct = new TGCheckButton(vCol2, "No structure (bulk only)");
        vCol2->AddFrame(fChkNoStruct, new TGLayoutHints(kLHintsLeft | kLHintsTop, 20, 5, 4, 5));
        fChkNoStruct->Connect("Toggled(Bool_t)", "RSDStripPanel", this, "DoNoStructToggle()");

        fLblStripR = new TGLabel(vCol2, "Total R_strip: 9999.9 kOhm");
        vCol2->AddFrame(fLblStripR, entL);

        vCol2->AddFrame(new TGLabel(vCol2, "Strip R @ Pad (Edge) [Ohm/sq]:"), lblL);
        fNumRStripPad = new TGNumberEntry(vCol2, 25, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 1, 100000);
        vCol2->AddFrame(fNumRStripPad, entL);

        vCol2->AddFrame(new TGLabel(vCol2, "Strip R @ Mid (Center) [Ohm/sq]:"), lblL);
        fNumRStripMid = new TGNumberEntry(vCol2, 25, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 1, 100000);
        vCol2->AddFrame(fNumRStripMid, entL);

        vCol2->AddFrame(new TGLabel(vCol2, "Strip Gap at midpoint [um]:"), entL);
        fNumStripGap = new TGNumberEntry(vCol2, 0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 3000);
        vCol2->AddFrame(fNumStripGap, entL);

        vCol2->AddFrame(new TGLabel(vCol2, "Strip Gap at pad [um]:"), lblL);
        fNumStripGapPad = new TGNumberEntry(vCol2, 0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 490);
        vCol2->AddFrame(fNumStripGapPad, entL);

        fNumRStripPad->Connect("ValueSet(Long_t)", "RSDStripPanel", this, "UpdateStripRLabel()");
        fNumRStripMid->Connect("ValueSet(Long_t)", "RSDStripPanel", this, "UpdateStripRLabel()");
        UpdateStripRLabel();

        BoldHead(vCol2, "--- Electronics ---");
        vCol2->AddFrame(new TGLabel(vCol2, "Input Impedance R_in [Ohm]:"), lblL);
        fNumRin = new TGNumberEntry(vCol2, 50, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 1000000);
        vCol2->AddFrame(fNumRin, entL);

        vCol2->AddFrame(new TGLabel(vCol2, "Pulse Duration [ns]:"), lblL);
        fNumDur = new TGNumberEntry(vCol2, 2.0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 0.1, 100.0);
        vCol2->AddFrame(fNumDur, entL);

        BoldHead(vCol2, "--- Threshold scan ---");
        vCol2->AddFrame(new TGLabel(vCol2, "Threshold for scans [%]:"), lblL);
        fNumThreshold = new TGNumberEntry(vCol2, 5.0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 0.1, 50.0);
        vCol2->AddFrame(fNumThreshold, entL);

        TGTextButton *btnScan6Limit = new TGTextButton(vCol2, "Pad 6 < threshold");
        btnScan6Limit->Connect("Clicked()", "RSDStripPanel", this, "DoScanPad6Limit()");
        vCol2->AddFrame(btnScan6Limit, btnL); btnScan6Limit->Resize(250, 30);

        TGTextButton *btnScan67Limit = new TGTextButton(vCol2, "Pad 6 & 7 < threshold");
        btnScan67Limit->Connect("Clicked()", "RSDStripPanel", this, "DoScan67Limit()");
        vCol2->AddFrame(btnScan67Limit, btnL); btnScan67Limit->Resize(250, 30);

        BoldHead(vCol2, "--- Diagnostics ---");
        TGTextButton *btnResMatrix = new TGTextButton(vCol2, "Inter-pad resistance");
        btnResMatrix->Connect("Clicked()", "RSDStripPanel", this, "DoResistanceMatrix()");
        vCol2->AddFrame(btnResMatrix, btnL);
        btnResMatrix->Resize(250, 30);
        btnResMatrix->SetBackgroundColor(0xfff0cc);

        fStatusLabel = new TGTextButton(vCol3, "  Ready  ");
        fStatusLabel->SetBackgroundColor(0xd0f0d0);
        fStatusLabel->SetForegroundColor(0x000000);
        vCol3->AddFrame(fStatusLabel, new TGLayoutHints(kLHintsCenterX | kLHintsTop, 5, 5, 8, 4));

        TGTextButton *btnGrid = new TGTextButton(vCol3, "Compute resistance matrix (mesh = 1% pitch)");
        btnGrid->Connect("Clicked()", "RSDStripPanel", this, "DoCalculateGrid()");
        vCol3->AddFrame(btnGrid, btnL); btnGrid->Resize(250, 40);
        btnGrid->SetBackgroundColor(0x44bb44);
        btnGrid->SetForegroundColor(0xffffff);

        BoldHead(vCol3, "--- Single Event ---");
        TGTextButton *btnSim = new TGTextButton(vCol3, "Inject charge and compute sharing");
        btnSim->Connect("Clicked()", "RSDStripPanel", this, "DoInjectCharge()");
        vCol3->AddFrame(btnSim, btnL); btnSim->Resize(250, 40);
        btnSim->SetBackgroundColor(0xffe4b0);

        TGTextButton *btnTrans = new TGTextButton(vCol3, "Inject charge and compute signals (4 central pads)");
        btnTrans->Connect("Clicked()", "RSDStripPanel", this, "DoTransient()");
        vCol3->AddFrame(btnTrans, btnL); btnTrans->Resize(250, 40);
        btnTrans->SetBackgroundColor(0xffe4b0);

        BoldHead(vCol3, "--- 2D Scans ---");
        TGTextButton *btnScanC = new TGTextButton(vCol3, "Fraction of charge in the 4 central pads");
        btnScanC->Connect("Clicked()", "RSDStripPanel", this, "DoMapScanCentral()");
        vCol3->AddFrame(btnScanC, btnL); btnScanC->Resize(250, 30);
        btnScanC->SetBackgroundColor(0xc0e8ff);

        TGTextButton *btnScan6 = new TGTextButton(vCol3, "Fraction of charge in pad 6");
        btnScan6->Connect("Clicked()", "RSDStripPanel", this, "DoMapScanPad6()");
        vCol3->AddFrame(btnScan6, btnL); btnScan6->Resize(250, 30);
        btnScan6->SetBackgroundColor(0xc0e8ff);

        TGTextButton *btnScanMax1 = new TGTextButton(vCol3, "Fraction of charge in any central pads");
        btnScanMax1->Connect("Clicked()", "RSDStripPanel", this, "DoMapScanMax1()");
        vCol3->AddFrame(btnScanMax1, btnL); btnScanMax1->Resize(250, 30);
        btnScanMax1->SetBackgroundColor(0xc0e8ff);

        TGTextButton *btnScanMax2 = new TGTextButton(vCol3, "Fraction of charge in the 2 highest central pads");
        btnScanMax2->Connect("Clicked()", "RSDStripPanel", this, "DoMapScanMax2()");
        vCol3->AddFrame(btnScanMax2, btnL); btnScanMax2->Resize(250, 30);
        btnScanMax2->SetBackgroundColor(0xc0e8ff);

        BoldHead(vCol3, "--- Position Reconstruction ---");
        TGTextButton *btnRes = new TGTextButton(vCol3, "Compute Unbalance reconstruction");
        btnRes->Connect("Clicked()", "RSDStripPanel", this, "DoResolutionAnalysis()");
        vCol3->AddFrame(btnRes, btnL);
        btnRes->Resize(250, 45);
        btnRes->SetBackgroundColor(0xddddff);

        TGTextButton *btnSharing = new TGTextButton(vCol3, "Compute sharing templates");
        btnSharing->Connect("Clicked()", "RSDStripPanel", this, "DoSharingTemplates()");
        vCol3->AddFrame(btnSharing, btnL);
        btnSharing->Resize(250, 30);
        btnSharing->SetBackgroundColor(0xccffcc);

        fBtnExport = new TGTextButton(vCol3, "Export sharing template to file");
        fBtnExport->Connect("Clicked()", "RSDStripPanel", this, "DoExportTemplate()");
        vCol3->AddFrame(fBtnExport, btnL);
        fBtnExport->Resize(250, 30);
        fBtnExport->SetBackgroundColor(0xccffcc);

        fBtnReco = new TGTextButton(vCol3, "Compute sharing templates resolution");
        TGTextButton *btnReco = fBtnReco;
        btnReco->Connect("Clicked()", "RSDStripPanel", this, "DoTemplateReconstruction()");
        vCol3->AddFrame(btnReco, btnL);
        btnReco->Resize(250, 30);
        btnReco->SetBackgroundColor(0xccffcc);

        vCol3->AddFrame(new TGLabel(vCol3, "N events for Template Reco:"), lblL);
        fNumRecoEvents = new TGNumberEntry(vCol3, 5000, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 100, 100000);
        vCol3->AddFrame(fNumRecoEvents, entL);

        BoldHead(vCol3, "--- Cache ---");
        TGTextButton *btnLoadCache = new TGTextButton(vCol3, "Load params from cache file");
        btnLoadCache->Connect("Clicked()", "RSDStripPanel", this, "DoLoadFromCache()");
        vCol3->AddFrame(btnLoadCache, btnL); btnLoadCache->Resize(250, 30);

        BoldHead(vCol3, "--- Config File ---");
        {
          TGHorizontalFrame *hCfgName = new TGHorizontalFrame(vCol3);
          fCfgFile = new TGTextEntry(hCfgName, "ConfigFiles/DCRSD_config.txt");
          fCfgFile->Resize(210, 22);
          hCfgName->AddFrame(fCfgFile, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));
          TGTextButton *btnBrowse = new TGTextButton(hCfgName, "...");
          btnBrowse->Connect("Clicked()", "RSDStripPanel", this, "DoBrowseConfig()");
          btnBrowse->Resize(30, 22);
          hCfgName->AddFrame(btnBrowse, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
          vCol3->AddFrame(hCfgName, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0));
        }
        TGHorizontalFrame *hCfg = new TGHorizontalFrame(vCol3);
        TGTextButton *btnSaveCfg = new TGTextButton(hCfg, "Save Config");
        btnSaveCfg->Connect("Clicked()", "RSDStripPanel", this, "DoSaveConfig()");
        btnSaveCfg->Resize(122, 30);
        hCfg->AddFrame(btnSaveCfg, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
        TGTextButton *btnLoadCfg = new TGTextButton(hCfg, "Load Config");
        btnLoadCfg->Connect("Clicked()", "RSDStripPanel", this, "DoLoadConfig()");
        btnLoadCfg->Resize(122, 30);
        hCfg->AddFrame(btnLoadCfg, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
        vCol3->AddFrame(hCfg, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));

        // All buttons except "Calculate grid" are disabled until DoCalculateGrid() completes
        fDepButtons = { btnSim, btnTrans, btnScanC, btnScan6, btnScanMax1, btnScanMax2,
                        btnScan6Limit, btnScan67Limit, btnRes, btnSharing, btnResMatrix };
        for (auto btn : fDepButtons) btn->SetEnabled(kFALSE);
        fBtnReco->SetEnabled(kFALSE);
        fBtnExport->SetEnabled(kFALSE);  // enabled only after DoSharingTemplates()


        controlFrame->AddFrame(vCol1, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 5, 5, 5, 5));
        controlFrame->AddFrame(vCol2, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 5, 5, 5, 5));
        controlFrame->AddFrame(vCol3, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 5, 5, 5, 5));

        TGTextButton *btnQuit = new TGTextButton(vCol3, "Quit");
        btnQuit->Connect("Clicked()", "TApplication", gApplication, "Terminate()");
        btnQuit->SetBackgroundColor(0xcc2222);
        btnQuit->SetForegroundColor(0xffffff);
        vCol3->AddFrame(btnQuit, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 10, 5));
        btnQuit->Resize(120, 35);

        hFrame->AddFrame(controlFrame, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 10, 10, 10, 10));
        fEcanvas = new TRootEmbeddedCanvas("Ecanvas", hFrame, 600, 600);
        hFrame->AddFrame(fEcanvas, new TGLayoutHints(kLHintsCenterY, 10, 10, 10, 10));
        AddFrame(hFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

        SetWindowName("DCRSD Surface Simulator - v1.0");
        MapSubwindows(); Resize(GetDefaultSize()); MapWindow();
        fEcanvas->GetCanvas()->Divide(2,2);
    }

    ~RSDStripPanel() {
        if(fSolverStatic) delete fSolverStatic; if(fMatrixStatic) delete fMatrixStatic;
        if(fSolverTrans) delete fSolverTrans; if(fMatrixTrans) delete fMatrixTrans;
    }

    void BuildTypeMaps(double pp, double pad_r = 60.0) {
        std::fill(fTypeMap.begin(), fTypeMap.end(), 0);
        std::fill(fPadOwnerMap.begin(), fPadOwnerMap.end(), 0);
        double p = pp / 100.0; // grid pitch
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double x = i * p, y = j * p;
                int idx = i * n + j;
                bool is_s = false;
                double min_dp = 1e9;
                int closest_pad = -1;
                for (int r = 0; r < 4; r++) {
                    double ds = std::min(std::abs(x - r*pp), std::abs(y - r*pp));
                    if (ds <= pp*0.025) is_s = true;
                    for (int c = 0; c < 4; c++) {
                        double dp = TMath::Hypot(x - r*pp, y - c*pp);
                        if (dp <= pad_r) fTypeMap[idx] = (r * 4 + c) + 1;
                        if (dp < min_dp) { min_dp = dp; closest_pad = (r*4+c)+1; }
                    }
                }
                if (is_s) fPadOwnerMap[idx] = closest_pad;
            }
        }
    }

    bool TrenchEnabled() { return fChkTrench && fChkTrench->IsOn(); }

    // Center node of pad k (1..16): the single representative that drains to GND via R_in
    int PadCenterNode(int pad_k) {
        int r = (pad_k-1)/4, c = (pad_k-1)%4;
        int pi = std::max(0, std::min(n-1, (int)TMath::Nint(r * PP() / Pitch())));
        int pj = std::max(0, std::min(n-1, (int)TMath::Nint(c * PP() / Pitch())));
        return pi*n + pj;
    }

    void SetStatus(const char* msg, bool busy) {
        ULong_t col = busy ? 0xffcc00 : 0xd0f0d0;
        fStatusLabel->SetText(new TGHotString(msg));
        fStatusLabel->ChangeBackground(col);         // aggiorna anche la X11 window background
        fStatusLabel->SetForegroundColor(0x000000);
        gClient->NeedRedraw(fStatusLabel, kTRUE);
        gClient->NeedRedraw((TGWindow*)fStatusLabel->GetParent(), kTRUE);
        if (busy) gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kWatch));
        else       gVirtualX->SetCursor(GetId(), kNone);
        gSystem->ProcessEvents();
        gVirtualX->Update(1);
        gSystem->Sleep(busy ? 80 : 50);
        gSystem->ProcessEvents();
        gVirtualX->Update(1);
    }

    void DoTrenchToggle() {
        bool trench = TrenchEnabled();
        fNumRStripPad->SetState(!trench);
        fNumRStripMid->SetState(!trench);
        fNumStripGap->SetState(!trench);
        fNumStripGapPad->SetState(!trench);
        UpdateStripRLabel();
    }

    void DoNoStructToggle() {
        bool ns = fChkNoStruct->IsOn();
        if (ns)
            fNumStripGap->SetNumber(PP());
        else
            fNumStripGap->SetNumber(0);
        fChkTrench->SetEnabled(!ns);
        bool disableStrips = ns || TrenchEnabled();
        fNumRStripPad->SetState(!disableStrips);
        fNumRStripMid->SetState(!disableStrips);
        fNumStripGap->SetState(!disableStrips);
        fNumStripGapPad->SetState(!disableStrips);
        UpdateStripRLabel();
    }

    void UpdateStripRLabel() {
        if (TrenchEnabled()) {
            fLblStripR->SetText("Total R_strip: N/A (trench active)");
            return;
        }
        double R_sP   = fNumRStripPad->GetNumber();
        double R_sM   = fNumRStripMid->GetNumber();
        double pp          = PP();
        // Strip width = 2 * (PP*0.025) from GetEdgeConductance threshold
        double strip_width = 2.0 * 0.025 * pp;   // = 0.05 * pp  (always 20 sq pad-to-pad)
        double N_sq        = pp / strip_width;    // = 20, independent of PP
        // Resistance: N_sq/2 squares per half, average resistivity (R_sP+R_sM)/2
        double R_tot       = (N_sq / 2.0) * (R_sP + R_sM);
        TString txt;
        if (R_tot >= 1e6)
            txt = TString::Format("Total R_strip: %.2f MOhm", R_tot / 1e6);
        else if (R_tot >= 1e3)
            txt = TString::Format("Total R_strip: %.1f kOhm", R_tot / 1e3);
        else
            txt = TString::Format("Total R_strip: %.0f Ohm", R_tot);
        fLblStripR->SetText(txt.Data());
        fLblStripR->Resize(fLblStripR->GetDefaultWidth(), fLblStripR->GetDefaultHeight());
        ((TGCompositeFrame*)fLblStripR->GetParent())->Layout();
    }

    // Returns the conductance of the edge between node (i,j) and node (i+di, j+dj).
    //
    // Key rules:
    //   0. Trench mode: edges crossing y=500,1000 (vertical) or x=500,1000 (horizontal)
    //      have zero conductance, except within pad_r of an electrode (trench stops at pad).
    //   1. Near a pad (within pad_radius + half_Pitch()): always R_bulk.
    //      The metal pad provides a perfect contact — strip stops here.
    //   2. Edges ALONG a strip direction: g = 1/R_bulk + 1/R_strip (parallel conductors).
    //   3. Edges PERPENDICULAR to a strip: g = 1/R_bulk (bulk only).
    //   4. Bulk: g = 1/R_bulk.
    double GetEdgeConductance(int i, int j, int di, int dj) {
        // Rule 0: Trench along electrode lines (x,y = 500, 1000 µm).
        // Cuts cross-trench bulk edges, but leaves gaps at the electrode pads.
        if (TrenchEnabled()) {
            const double pad_r = PadRadius();
            if (di == 0) {
                // Vertical edge: check if it crosses y = 500 or y = 1000
                for (int k = 1; k <= 2; k++) {
                    if ((j < k*100) != (j + dj < k*100)) {
                        double xmid = i * Pitch();
                        bool near_pad = false;
                        for (int c = 0; c < 4 && !near_pad; c++)
                            near_pad = std::abs(xmid - c*PP()) <= pad_r;
                        if (!near_pad) return 0.0;
                    }
                }
            } else {
                // Horizontal or diagonal edge: check both x and y boundary crossings
                for (int k = 1; k <= 2; k++) {
                    if ((i < k*100) != (i + di < k*100)) {
                        double ymid = j * Pitch();
                        bool near_pad = false;
                        for (int c = 0; c < 4 && !near_pad; c++)
                            near_pad = std::abs(ymid - c*PP()) <= pad_r;
                        if (!near_pad) return 0.0;
                    }
                    if (dj != 0 && (j < k*100) != (j + dj < k*100)) {
                        double xmid = i * Pitch();
                        bool near_pad = false;
                        for (int c = 0; c < 4 && !near_pad; c++)
                            near_pad = std::abs(xmid - c*PP()) <= pad_r;
                        if (!near_pad) return 0.0;
                    }
                }
            }
        }
        double R_bulk  = fNumRsheet->GetNumber();
        double R_sP    = fNumRStripPad->GetNumber();
        double R_sM    = fNumRStripMid->GetNumber();
        double xmid = (i + 0.5*di) * Pitch();
        double ymid = (j + 0.5*dj) * Pitch();

        // With trench active, strips have no effect (only bulk)
        if (TrenchEnabled()) return 1.0/R_bulk;

        // dj==0  →  horizontal edge (x-direction): check horizontal strips at y = k*500
        // di==0  →  vertical edge   (y-direction): check vertical   strips at x = k*500
        double gap_half     = fNumStripGap->GetNumber() / 2.0;
        double gap_pad      = fNumStripGapPad->GetNumber();   // exclusion radius around each electrode
        if (dj == 0) {
            for (int k = 0; k < 4; k++) {
                if (std::abs(ymid - k*PP()) <= PP()*0.025) {
                    // Gap at midpoint between electrodes
                    if (gap_half > 0)
                        for (int c = 0; c < 3; c++)
                            if (std::abs(xmid - (c*PP() + PP()/2.0)) < gap_half) return 1.0/R_bulk;
                    // Gap at pad: strip ends before reaching the electrode
                    double min_d = 1e9;
                    for (int c = 0; c < 4; c++) {
                        double d = std::abs(xmid - c*PP());
                        if (d < min_d) min_d = d;
                    }
                    if (gap_pad > 0 && min_d < gap_pad) return 1.0/R_bulk;
                    double d = std::min(min_d, PP()/2.0);
                    double r_strip = R_sP + (R_sM - R_sP) * (2.0*d/PP());
                    return 1.0/R_bulk + 1.0/r_strip;
                }
            }
        } else {
            for (int k = 0; k < 4; k++) {
                if (std::abs(xmid - k*PP()) <= PP()*0.025) {
                    // Gap at midpoint between electrodes
                    if (gap_half > 0)
                        for (int c = 0; c < 3; c++)
                            if (std::abs(ymid - (c*PP() + PP()/2.0)) < gap_half) return 1.0/R_bulk;
                    // Gap at pad: strip ends before reaching the electrode
                    double min_d = 1e9;
                    for (int c = 0; c < 4; c++) {
                        double d = std::abs(ymid - c*PP());
                        if (d < min_d) min_d = d;
                    }
                    if (gap_pad > 0 && min_d < gap_pad) return 1.0/R_bulk;
                    double d = std::min(min_d, PP()/2.0);
                    double r_strip = R_sP + (R_sM - R_sP) * (2.0*d/PP());
                    return 1.0/R_bulk + 1.0/r_strip;
                }
            }
        }

        return 1.0/R_bulk;
    }

    double GetLocalResistance(double x, double y) {
        double R_bulk = fNumRsheet->GetNumber();
        double R_sP = fNumRStripPad->GetNumber();
        double R_sM = fNumRStripMid->GetNumber();
        double gap_half = fNumStripGap->GetNumber() / 2.0;
        double gap_pad  = fNumStripGapPad->GetNumber();
        double min_dist_pad_center = 1e9;
        bool is_s = false;

        for(int k=0; k<4; k++) {
            bool on_hstrip = std::abs(y - k*PP()) <= PP()*0.025;
            bool on_vstrip = std::abs(x - k*PP()) <= PP()*0.025;

            // Check midpoint gap
            if (on_hstrip && gap_half > 0)
                for(int c=0; c<3; c++)
                    if(std::abs(x - (c*PP()+PP()/2.0)) < gap_half) { on_hstrip = false; break; }
            if (on_vstrip && gap_half > 0)
                for(int c=0; c<3; c++)
                    if(std::abs(y - (c*PP()+PP()/2.0)) < gap_half) { on_vstrip = false; break; }

            // Check pad gap: strip ends before reaching each electrode
            if (on_hstrip && gap_pad > 0)
                for(int c=0; c<4; c++)
                    if(std::abs(x - c*PP()) < gap_pad) { on_hstrip = false; break; }
            if (on_vstrip && gap_pad > 0)
                for(int c=0; c<4; c++)
                    if(std::abs(y - c*PP()) < gap_pad) { on_vstrip = false; break; }

            if (on_hstrip || on_vstrip) is_s = true;

            for(int c=0; c<4; c++) {
                double dp = TMath::Hypot(x - k*PP(), y - c*PP());
                if(dp - PadRadius() <= 0) is_s = true;
                if(dp < min_dist_pad_center) min_dist_pad_center = dp;
            }
        }

        if (is_s) {
            double d = (min_dist_pad_center > PP()/2.0) ? PP()/2.0 : min_dist_pad_center;
            return R_sP + (R_sM - R_sP) * (2.0*d/PP());
        }
        return R_bulk;
    }

    void CalculatePropagation(double injX, double injY, std::vector<double>& out_delays, std::vector<TGraph*>& out_paths) {
        double thick_um = fNumThick->GetNumber();
        
        int targets[4] = {6, 7, 10, 11};
        int targetCenters[4] = {100*n + 100, 100*n + 200, 200*n + 100, 200*n + 200};

        for (int idx = 0; idx < 4; idx++) {
            int targetID = targets[idx];
            int targetCenter = targetCenters[idx];

            int injI = TMath::Nint(injX/Pitch());
            int injJ = TMath::Nint(injY/Pitch());
            if(injI < 0) injI = 0; if(injI >= n) injI = n - 1;
            if(injJ < 0) injJ = 0; if(injJ >= n) injJ = n - 1;
            int startNode = injI * n + injJ;

            std::vector<double> dist(nTot, 1e15); 
            std::vector<int> parent(nTot, -1);
            dist[startNode] = 0.0;

            using pdi = std::pair<double, int>;
            std::priority_queue<pdi, std::vector<pdi>, std::greater<pdi>> pq;
            pq.push({0.0, startNode});

            int di[8] = {1, -1, 0, 0, 1, 1, -1, -1};
            int dj[8] = {0, 0, 1, -1, 1, -1, 1, -1};
            double ds_mm[8] = { 0.005, 0.005, 0.005, 0.005, 0.005*1.414213, 0.005*1.414213, 0.005*1.414213, 0.005*1.414213 };

            while(!pq.empty()) {
                double d = pq.top().first;
                int u = pq.top().second;
                pq.pop();

                if (d > dist[u]) continue;
                if (u == targetCenter) break; 

                int ui = u / n;
                int uj = u % n;
                double r_u = GetLocalResistance(ui * Pitch(), uj * Pitch());

                for (int k = 0; k < 8; k++) {
                    int vi = ui + di[k];
                    int vj = uj + dj[k];
                    if (vi >= 0 && vi < n && vj >= 0 && vj < n) {
                        int v = vi * n + vj;
                        
                        int pad_at_v = fTypeMap[v];
                        int strip_owner = fPadOwnerMap[v];
                        
                        bool is_foreign_pad = (pad_at_v > 0 && pad_at_v != targetID);
                        bool is_foreign_strip = (strip_owner > 0 && strip_owner != targetID);
                        
                        if ((is_foreign_pad || is_foreign_strip) && v != startNode) {
                            continue; 
                        }
                        
                        double r_v = GetLocalResistance(vi * Pitch(), vj * Pitch());
                        double r_avg = (r_u + r_v) / 2.0; 
                        
                        double local_v_factor = 0.56 * std::sqrt((r_avg / PP()*1.5) * (50.0 / thick_um));
                        double cost = ds_mm[k] * local_v_factor;
                        
                        if (dist[u] + cost < dist[v]) {
                            dist[v] = dist[u] + cost;
                            parent[v] = u;
                            pq.push({dist[v], v});
                        }
                    }
                }
            }

            out_delays.push_back(dist[targetCenter]);
            
            TGraph* gr = new TGraph();
            int curr = targetCenter;
            int pt = 0;
            while(curr != -1 && curr != startNode) {
                double cx = (curr / n) * Pitch();
                double cy = (curr % n) * Pitch();
                gr->SetPoint(pt++, cx, cy);
                curr = parent[curr];
            }
            if (curr == startNode) {
                double cx = (startNode / n) * Pitch();
                double cy = (startNode % n) * Pitch();
                gr->SetPoint(pt++, cx, cy);
            }
            gr->SetLineWidth(3);
            out_paths.push_back(gr);
        }
    }

    // Collects pad currents from solution vector V.
    // R_in=0: V[pad]=0, current from edge flows; R_in>0: I[pad]=V[pad]/R_in
    void CollectCurrents(const TVectorD &V, std::vector<double> &cur, double &iTot) {
        std::fill(cur.begin(), cur.end(), 0.0); iTot = 0;
        double R_in = fNumRin->GetNumber();
        if (R_in > 0) {
            // Only the representative (center) node of each pad drains via R_in
            for (int k=1; k<=16; k++) {
                double I = V(PadCenterNode(k)) / R_in;
                cur[k] += I; iTot += I;
            }
        } else {
            // R_in=0: V[pad]=0, collect from mesh-adjacent surface nodes
            const int di[]={1,-1,0,0}, dj[]={0,0,1,-1};
            for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
                int idx = i*n+j;
                if (fTypeMap[idx]==0) for (int k=0; k<4; k++) {
                    int ni=i+di[k], nj=j+dj[k];
                    if (ni>=0 && ni<n && nj>=0 && nj<n && fTypeMap[ni*n+nj]>0) {
                        double val = V(idx)*GetEdgeConductance(i,j,di[k],dj[k]);
                        cur[fTypeMap[ni*n+nj]] += val; iTot += val;
                    }
                }
            }
        }
    }

    TDecompSparse* GetSolver(bool transient, double dt = 1e-11, double C_node = 1e-15) {
        double R_sh = fNumRsheet->GetNumber(); double R_sP = fNumRStripPad->GetNumber();
        double R_sM = fNumRStripMid->GetNumber(); double th = fNumThick->GetNumber();

        if (!transient) {
            bool trench = TrenchEnabled();
            double gap = fNumStripGap->GetNumber(); double gap_pad = fNumStripGapPad->GetNumber(); double pp = PP();
            double pad_r = PadRadius();
            if (pp != fLastPPStat || pad_r != fLastPadRadiusStat) BuildTypeMaps(pp, pad_r);
            double rin = fNumRin->GetNumber();
            if (fSolverStatic && fMatrixStatic && R_sh==fLastRsheetStat && R_sP==fLastRStripPStat && R_sM==fLastRStripMStat && gap==fLastGapStat && gap_pad==fLastGapPadStat && pp==fLastPPStat && trench==fLastTrenchStat && rin==fLastRinStat && pad_r==fLastPadRadiusStat) return fSolverStatic;
            if (fSolverStatic) delete fSolverStatic; if (fMatrixStatic) delete fMatrixStatic;
            fMatrixStatic = nullptr;
            // --- On-disk cache ---
            gSystem->mkdir("rsd_cache", kTRUE);
            gSystem->mkdir(TString::Format("rsd_cache/PP%.0f", pp).Data(), kTRUE);
            TString cacheDir;
            if (trench)
                cacheDir = TString::Format("rsd_cache/PP%.0f/trench", pp);
            else
                cacheDir = TString::Format("rsd_cache/PP%.0f/RsP%.0f_RsM%.0f", pp, R_sP, R_sM);
            gSystem->mkdir(cacheDir.Data(), kTRUE);
            TString cacheKey;
            if (trench)
                cacheKey = TString::Format(
                    "PP%.0f_Rsh%.0f_tr1_Rin%.2f_Rpad%.0f",
                    pp, R_sh, rin, pad_r);
            else
                cacheKey = TString::Format(
                    "PP%.0f_Rsh%.0f_RsP%.0f_RsM%.0f_gap%.1f_gapP%.1f_tr0_Rin%.2f_Rpad%.0f",
                    pp, R_sh, R_sP, R_sM, gap, gap_pad, rin, pad_r);
            TString cachePath = cacheDir + "/matrix_" + cacheKey + ".root";
            if (!gSystem->AccessPathName(cachePath.Data())) {  // file exists
                TFile *fcache = TFile::Open(cachePath.Data(), "READ");
                if (fcache && !fcache->IsZombie()) {
                    TMatrixDSparse *mtmp = (TMatrixDSparse*)fcache->Get("KCL");
                    if (mtmp && mtmp->GetNrows() == nTot) {
                        fMatrixStatic = new TMatrixDSparse(*mtmp);  // clone from cache
                        std::cout << ">>> Matrix loaded from cache: " << cachePath << std::endl;
                    }
                    fcache->Close(); delete fcache;
                }
            }
            if (!fMatrixStatic) {  // not in cache: compute and save
                fMatrixStatic = new TMatrixDSparse(nTot, nTot);
                BuildMatrix(*fMatrixStatic, false, 0.0);
                TFile *fcache = TFile::Open(cachePath.Data(), "RECREATE");
                if (fcache && !fcache->IsZombie()) {
                    fcache->cd();
                    fMatrixStatic->Write("KCL");
                    fcache->Close(); delete fcache;
                    std::cout << ">>> Matrix saved to cache: " << cachePath << std::endl;
                }
            }
            fSolverStatic = new TDecompSparse(); fSolverStatic->SetMatrix(*fMatrixStatic); fSolverStatic->Decompose();
            fLastRsheetStat = R_sh; fLastRStripPStat = R_sP; fLastRStripMStat = R_sM; fLastGapStat = gap; fLastGapPadStat = gap_pad; fLastPPStat = pp; fLastTrenchStat = trench; fLastRinStat = rin; fLastPadRadiusStat = pad_r;
            return fSolverStatic;
        } else {
            bool trench = TrenchEnabled();
            double gap = fNumStripGap->GetNumber(); double pp = PP();
        double rin = fNumRin->GetNumber();
            if (fSolverTrans && fMatrixTrans && R_sh==fLastRsheetTrans && R_sP==fLastRStripPTrans && R_sM==fLastRStripMTrans && th==fLastThickTrans && gap==fLastGapTrans && pp==fLastPPTrans && dt==fLastDtTrans && C_node==fLastCnodeTrans && trench==fLastTrenchTrans && rin==fLastRinTrans) return fSolverTrans;
            if (fSolverTrans) delete fSolverTrans; if (fMatrixTrans) delete fMatrixTrans;
            fMatrixTrans = new TMatrixDSparse(nTot, nTot); BuildMatrix(*fMatrixTrans, true, C_node/dt);
            fSolverTrans = new TDecompSparse(); fSolverTrans->SetMatrix(*fMatrixTrans); fSolverTrans->Decompose();
            fLastRsheetTrans = R_sh; fLastRStripPTrans = R_sP; fLastRStripMTrans = R_sM; fLastThickTrans = th; fLastGapTrans = gap; fLastPPTrans = pp; fLastDtTrans = dt; fLastCnodeTrans = C_node; fLastTrenchTrans = trench; fLastRinTrans = rin;
            return fSolverTrans;
        }
    }

    void BuildMatrix(TMatrixDSparse &M, bool transient, double diagAdd) {
        double R_in = fNumRin->GetNumber();
        bool use_Rin = (R_in > 0);
        // Short-circuit conductance within metal pad (enforces equipotential)
        const double G_short = 1e6;
        const int di[] = {1,-1,0,0}, dj[] = {0,0,1,-1};
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                int idx = i*n + j;
                int pad_k = fTypeMap[idx];

                // Pad node, R_in=0: force V=0 (perfect sink)
                if (pad_k > 0 && !use_Rin) { M(idx, idx) = 1.0 + diagAdd; continue; }

                // Pad node, R_in>0: equipotential metal pad model
                if (pad_k > 0 && use_Rin) {
                    double sum_g = diagAdd;
                    // Only the center node of this pad drains to GND via R_in
                    if (idx == PadCenterNode(pad_k)) sum_g += 1.0/R_in;
                    for (int k2=0; k2<4; k2++) {
                        int ni=i+di[k2], nj=j+dj[k2];
                        if (ni>=0&&ni<n&&nj>=0&&nj<n) {
                            int npad = fTypeMap[ni*n+nj];
                            if (npad == 0) {
                                // Surface neighbor: current flows in from resistive layer
                                double cond = GetEdgeConductance(i, j, di[k2], dj[k2]);
                                sum_g += cond; M(idx, ni*n+nj) = -cond;
                            } else if (npad == pad_k) {
                                // Same-pad neighbor: short circuit (equipotential metal)
                                sum_g += G_short; M(idx, ni*n+nj) = -G_short;
                            }
                            // Different pad: no connection (separate metal contacts)
                        }
                    }
                    M(idx, idx) = sum_g;
                    continue;
                }

                // Surface node: regular KCL
                double sum_g = diagAdd;
                for (int k2=0; k2<4; k2++) {
                    int ni=i+di[k2], nj=j+dj[k2];
                    if (ni>=0&&ni<n&&nj>=0&&nj<n) {
                        double cond = GetEdgeConductance(i, j, di[k2], dj[k2]);
                        sum_g += cond;
                        if (fTypeMap[ni*n+nj] == 0 || use_Rin) M(idx, ni*n+nj) = -cond;
                    }
                }
                M(idx, idx) = sum_g;
            }
        }
    }


    void DoCalculateGrid() {
        SetStatus("  Computing...  ", true);
        GetSolver(false);  // build and cache the static matrix/solver
        for (auto btn : fDepButtons) btn->SetEnabled(kTRUE);
        fBtnReco->SetEnabled(kFALSE);  // requires new template computation
        fTemplatesReady = false;

        // Draw local resistivity map in pad (1,1) already divided as after injection
        TCanvas *c = fEcanvas->GetCanvas(); c->Clear(); c->Divide(2,2);
        delete gROOT->FindObject("hResGrid");
        TH2D *hResGrid = new TH2D("hResGrid", "Local Resistivity [Ohm/sq]", n, 0, Dim(), n, 0, Dim());
        for (int i=0; i<n; i++)
            for (int j=0; j<n; j++) {
                int idx = i*n+j;
                hResGrid->SetBinContent(i+1, j+1, (fTypeMap[idx]>0) ? 0 : GetLocalResistance(i*Pitch(), j*Pitch()));
            }
        c->cd(1); gPad->SetLogz();
        hResGrid->Draw("COLZ");
        struct PadColor { int r, c, color; };
        PadColor padColors[] = {
            {0,0,kGray},{0,1,kGray},{0,2,kGray},{0,3,kGray},
            {1,0,kGray},{1,1,kRed},{1,2,kBlue},{1,3,kGray},
            {2,0,kGray},{2,1,kGreen+2},{2,2,kMagenta},{2,3,kGray},
            {3,0,kGray},{3,1,kGray},{3,2,kGray},{3,3,kGray}
        };
        for (auto &p : padColors) {
            TEllipse *el = new TEllipse(p.r*PP(), p.c*PP(), PadRadius(), PadRadius());
            el->SetFillColor(p.color); el->SetLineColor(p.color); el->Draw("SAME");
        }
        c->Modified();
        c->Update();
        gSystem->ProcessEvents();
        gVirtualX->Update(1);
        SetStatus("  Ready  ", false);
    }

void DoInjectCharge() {
    SetStatus("  Computing...  ", true);
  TDecompSparse *solver = GetSolver(false);
    double injX = fNumX->GetNumber(), injY = fNumY->GetNumber();
    TVectorD B(nTot); B.Zero();
    { int injI = std::max(0, std::min(n-1, TMath::Nint(injX/Pitch())));
      int injJ = std::max(0, std::min(n-1, TMath::Nint(injY/Pitch())));
      B(injI*n + injJ) = 1.0; }
    bool ok; TVectorD V = solver->Solve(B, ok);
    
    TCanvas *c = fEcanvas->GetCanvas(); c->Clear(); c->Divide(2,2);
    delete gROOT->FindObject("hRes"); delete gROOT->FindObject("hPot");
    delete gROOT->FindObject("hI"); delete gROOT->FindObject("hI6");
    delete gROOT->FindObject("hI7"); delete gROOT->FindObject("hI10"); delete gROOT->FindObject("hI11");
    c->cd(1); gPad->SetLogz(); TH2D *hRes = new TH2D("hRes", "Local Resistivity [Ohm/sq]", n, 0, Dim(), n, 0, Dim());
    c->cd(2); TH2D *hPot = new TH2D("hPot", "Static Potential Map", n, 0, Dim(), n, 0, Dim());
    
    std::vector<double> I_out(17, 0.0);
    double iTotal = 0;
    CollectCurrents(V, I_out, iTotal);

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            int idx = i*n+j;
            double r_local = GetLocalResistance(i*Pitch(), j*Pitch());
            hRes->SetBinContent(i+1, j+1, (fTypeMap[idx]>0) ? 0 : r_local);
            hPot->SetBinContent(i+1, j+1, V(idx));
        }
    }
    
    if (iTotal < 0.9)
        std::cout << ">>> WARNING: Total collected current = " << iTotal << " (expected 1.0)" << std::endl;
    double fracCentral = (iTotal > 0) ? (I_out[6]+I_out[7]+I_out[10]+I_out[11])/iTotal : 0;
    printf(">>> Central 4 pads fraction: %.1f %%\n", fracCentral * 100.0);
    
    // Pad colors: 6=red, 7=blue, 10=green+2, 11=magenta, others=gray
    struct PadColor { int r, c, color; };
    PadColor padColors[] = {
        {0,0,kGray},{0,1,kGray},{0,2,kGray},{0,3,kGray},
        {1,0,kGray},{1,1,kRed},{1,2,kBlue},{1,3,kGray},
        {2,0,kGray},{2,1,kGreen+2},{2,2,kMagenta},{2,3,kGray},
        {3,0,kGray},{3,1,kGray},{3,2,kGray},{3,3,kGray}
    };
    c->cd(1); hRes->Draw("COLZ");
    for (auto &p : padColors) {
        TEllipse *el = new TEllipse(p.r*PP(), p.c*PP(), PP()*0.05, PP()*0.05);
        el->SetFillColor(p.color); el->SetLineColor(p.color);
        el->Draw("SAME");
    }
    c->cd(2); hPot->Draw("COLZ");
    for (auto &p : padColors) {
        TEllipse *el = new TEllipse(p.r*PP(), p.c*PP(), PP()*0.05, PP()*0.05);
        el->SetFillColor(p.color); el->SetLineColor(p.color);
        el->Draw("SAME");
    }
    c->cd(3);
    TH1D *hI    = new TH1D("hI",    "Charge Fractions;Pad ID;Fraction", 16, 0.5, 16.5);
    TH1D *hI6   = new TH1D("hI6",   "", 16, 0.5, 16.5);
    TH1D *hI7   = new TH1D("hI7",   "", 16, 0.5, 16.5);
    TH1D *hI10  = new TH1D("hI10",  "", 16, 0.5, 16.5);
    TH1D *hI11  = new TH1D("hI11",  "", 16, 0.5, 16.5);
    for(int k=1; k<=16; k++) {
        double frac = (iTotal>0) ? I_out[k]/iTotal : 0;
        hI->SetBinContent(k, frac);
        if (k==6)  hI6->SetBinContent(k, frac);
        if (k==7)  hI7->SetBinContent(k, frac);
        if (k==10) hI10->SetBinContent(k, frac);
        if (k==11) hI11->SetBinContent(k, frac);
    }
    hI->SetFillColor(kGray); hI->Draw("BAR TEXT0");
    hI6->SetFillColor(kRed);      hI6->Draw("BAR SAME");
    hI7->SetFillColor(kBlue);     hI7->Draw("BAR SAME");
    hI10->SetFillColor(kGreen+2); hI10->Draw("BAR SAME");
    hI11->SetFillColor(kMagenta); hI11->Draw("BAR SAME");
    c->Update();
    SetStatus("  Ready  ", false);
}



    void DoTransient() {
        SetStatus("  Computing...  ", true);
        double thick_um = fNumThick->GetNumber();
        double thick = thick_um * 1e-6;
        double duration = fNumDur->GetNumber() * 1e-9; 
        double C_node = eps0 * epsSi * (Pitch()*1e-6) * (Pitch()*1e-6) / thick;
        double dt = 20e-12; int nSteps = (int)(duration / dt);
        
        TDecompSparse *solver = GetSolver(true, dt, C_node);
        TVectorD V(nTot); V.Zero(); TVectorD B(nTot); 
        double injX = fNumX->GetNumber(); double injY = fNumY->GetNumber();
        int injI = std::max(0, std::min(n-1, TMath::Nint(injX/Pitch())));
        int injJ = std::max(0, std::min(n-1, TMath::Nint(injY/Pitch())));
        int injIdx = injI*n + injJ;

        TGraph *gP6  = new TGraph(); 
        TGraph *gP7  = new TGraph(); 
        TGraph *gP10 = new TGraph(); 
        TGraph *gP11 = new TGraph();
        TGraph *gInput = new TGraph(); 

        gP6->SetLineColor(kRed); gP6->SetLineWidth(2);
        gP7->SetLineColor(kBlue); gP7->SetLineWidth(2);
        gP10->SetLineColor(kGreen+2); gP10->SetLineWidth(2);
        gP11->SetLineColor(kMagenta); gP11->SetLineWidth(2);
        gInput->SetLineColor(kBlack); gInput->SetLineStyle(2); gInput->SetLineWidth(2);

        std::vector<double> delays;
        std::vector<TGraph*> paths; 
        CalculatePropagation(injX, injY, delays, paths);
        for(auto p : paths) delete p; 
        
        double delay6  = delays[0];
        double delay7  = delays[1];
        double delay10 = delays[2];
        double delay11 = delays[3];

        double maxSig = 0;

        for (int t=0; t<nSteps; t++) {
            for(int i=0; i<nTot; i++) B(i) = (C_node/dt) * V(i);
            double inputVal = 0; if (t > 5 && t <= 10) { inputVal = 10.0; B(injIdx) += inputVal; }
            bool ok; V = solver->Solve(B, ok); 

            double i6=0, i7=0, i10=0, i11=0;
            const int di[]={1,-1,0,0}, dj[]={0,0,1,-1};
            for (int i=0; i<n; i++) { for (int j=0; j<n; j++) {
                int idx = i*n+j;
                if (fTypeMap[idx] == 0) {
                    for(int k=0; k<4; k++){
                        int ni=i+di[k], nj=j+dj[k];
                        if(ni>=0 && ni<n && nj>=0 && nj<n && fTypeMap[ni*n+nj]>0) {
                            double cond = GetEdgeConductance(i, j, di[k], dj[k]);
                            double cur  = V(idx) * cond;
                            int id = fTypeMap[ni*n+nj];
                            if (id==6) i6+=cur; else if (id==7) i7+=cur; else if (id==10) i10+=cur; else if (id==11) i11+=cur;
                        }
                    }
                }
            }}
            
            double time_ns = t * dt * 1e9;
            
            if (delay6 < 1e14) gP6->SetPoint(t, time_ns + delay6, i6); else gP6->SetPoint(t, time_ns, i6);
            if (delay7 < 1e14) gP7->SetPoint(t, time_ns + delay7, i7); else gP7->SetPoint(t, time_ns, i7);
            if (delay10 < 1e14) gP10->SetPoint(t, time_ns + delay10, i10); else gP10->SetPoint(t, time_ns, i10);
            if (delay11 < 1e14) gP11->SetPoint(t, time_ns + delay11, i11); else gP11->SetPoint(t, time_ns, i11);
            
            maxSig = std::max({maxSig, i6, i7, i10, i11});
            gInput->SetPoint(t, time_ns, inputVal);
            if (t % 100 == 0) gSystem->ProcessEvents(); 
        }

        double maxIn = TMath::MaxElement(gInput->GetN(), gInput->GetY());
        if (maxIn > 0 && maxSig > 0) { double scaleFactor = maxSig / maxIn; for (int k=0; k<gInput->GetN(); k++) gInput->GetY()[k] *= scaleFactor; }
        
        TCanvas *cScope = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cScope");
        if (!cScope) cScope = new TCanvas("cScope", "c1 - RSD Oscilloscope", 480, 360);
        cScope->cd(); cScope->Clear(); gPad->SetGrid();
        
        double valid_delay6 = (delay6 < 1e14) ? delay6 : 0;
        double valid_delay7 = (delay7 < 1e14) ? delay7 : 0;
        double valid_delay10 = (delay10 < 1e14) ? delay10 : 0;
        double valid_delay11 = (delay11 < 1e14) ? delay11 : 0;

        double max_t = (nSteps * dt * 1e9) + std::max({valid_delay6, valid_delay7, valid_delay10, valid_delay11});
        TH1F *hFrame = cScope->DrawFrame(0, 0, max_t, (maxSig>0)?maxSig*1.1:1.0);
        hFrame->SetTitle("Transient with Fermat Diffraction Delay;Time [ns];Current [a.u.]");

        gP6->Draw("L SAME"); gP7->Draw("L SAME"); gP10->Draw("L SAME"); gP11->Draw("L SAME"); gInput->Draw("L SAME");
        
        TLegend *leg = new TLegend(0.75, 0.6, 0.95, 0.9);
        leg->AddEntry(gInput, "Input Pulse (T0)", "l"); 
        leg->AddEntry(gP6, "Pad 6", "l"); 
        leg->AddEntry(gP7, "Pad 7", "l");
        leg->AddEntry(gP10, "Pad 10", "l"); 
        leg->AddEntry(gP11, "Pad 11", "l");
        leg->Draw(); cScope->Update(); cScope->RaiseWindow();
        SetStatus("  Ready  ", false);
    }

    void DoMapScanCentral() { RunScan(0); }
    void DoMapScanPad6()    { RunScan(1); }
    void DoScanPad6Limit()  { RunScan(2); }
    void DoScan67Limit()    { RunScan(3); }
    void DoMapScanMax1()    { RunScan(4); }
    void DoMapScanMax2()    { RunScan(5); }

    void RunScan(int mode) {
        SetStatus("  Computing...  ", true);
        TDecompSparse *solver = GetSolver(false);
        int    nx      = (int)fNumScanStep->GetNumber();
        double scanMin = (mode == 0) ? PP()      : 0.0;
        double scanMax = (mode == 0) ? 2.0*PP()  : Dim();
        const double step = PP() / nx;
        int nBins = (int)((scanMax - scanMin) / step) + 1;
        TString title;
        if (mode==0) title = "Central 4 Pads Fraction";
        else if (mode==1) title = "Pad 6 Fraction";
        else if (mode==2) title = "Limit: Pad 6 < 5%";
        else if (mode==3) title = "Blind Zone: Pad 6 & 7 < 5%";
        else if (mode==4) title = "Max 1 Pad Fraction";
        else if (mode==5) title = "Max 2 Pads Fraction";

        const char* cNames[] = { "cScan0","cScan1","cScan2","cScan3","cScan4","cScan5" };
        const char* cTitles[]= {
            "c2a - Central 4 pads fraction",
            "c2b - Pad 6 fraction",
            "c2c - Pad 6 < threshold",
            "c2d - Pad 6 & 7 < threshold",
            "c2e - Max 1 pad fraction",
            "c2f - Max 2 pads fraction"
        };
        TCanvas *cScan = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(cNames[mode]);
        int cW = (mode==0) ? 840 : 600;
        int cH = (mode==0) ? 840 : (mode==1 || mode==4 || mode==5) ? 600 : 420;
        if (mode == 0 && cScan) { delete cScan; cScan = nullptr; }  // force correct size for 2x2
        if (!cScan) cScan = new TCanvas(cNames[mode], cTitles[mode], cW, cH);
        cScan->cd(); cScan->Clear();
        if (mode == 0) cScan->Divide(2,2);

        double histMin = scanMin;
        double histMax = scanMax;

        TString hScanName = TString::Format("hScan%d", mode);
        delete gROOT->FindObject(hScanName);
        TH2D *hScan = new TH2D(hScanName, title + ";X [um];Y [um]",
                               nBins, histMin-step/2., histMax+step/2.,
                               nBins, histMin-step/2., histMax+step/2.);
        hScan->SetMinimum(0.7); hScan->SetMaximum(1.0);
        hScan->SetStats(0);

        TString hLeakName = TString::Format("hLeak%d", mode);
        delete gROOT->FindObject(hLeakName);
        TH2D *hLeak = nullptr;
        if (mode == 0) {
            hLeak = new TH2D(hLeakName, "Signal Leakage to outer pixels;X [um];Y [um]",
                             nBins, histMin-step/2., histMax+step/2.,
                             nBins, histMin-step/2., histMax+step/2.);
            hLeak->SetMinimum(0.0); hLeak->SetMaximum(0.3);
            hLeak->SetStats(0);
        }

        std::vector<double> sumReadX(nBins, 0), sum2ReadX(nBins, 0);
        std::vector<double> sumLostX(nBins, 0), sum2LostX(nBins, 0);
        std::vector<int>    cntProjX(nBins, 0);

        if (mode >= 2 && mode <= 3) gStyle->SetPalette(kTemperatureMap); else gStyle->SetPalette(kBird);

        for (int ix = 0; ix < nBins; ix++) {
            double x = scanMin + ix * step;
            for (int iy = 0; iy < nBins; iy++) {
                double y = scanMin + iy * step;
                TVectorD B_s(nTot); B_s.Zero(); int idx = TMath::Nint(x/Pitch())*n + TMath::Nint(y/Pitch()); if (idx >= nTot) idx = nTot - 1; B_s(idx) = 1.0;
                bool ok; TVectorD V_s = solver->Solve(B_s, ok);

                std::vector<double> i_out(17, 0.0); double iTotal = 0;
                CollectCurrents(V_s, i_out, iTotal);
                
                double val = 0; double f6 = (iTotal>0)?i_out[6]/iTotal:0; double f7 = (iTotal>0)?i_out[7]/iTotal:0;
                std::vector<double> central = {i_out[6], i_out[7], i_out[10], i_out[11]};
                std::sort(central.begin(), central.end(), std::greater<double>());

                if (mode == 0) val = (i_out[6]+i_out[7]+i_out[10]+i_out[11])/iTotal;
                else if (mode == 1) val = f6;
                else if (mode == 2) val = (f6 < fNumThreshold->GetNumber()/100.0) ? 1.0 : 0.0;
                else if (mode == 3) val = (f6 < fNumThreshold->GetNumber()/100.0 && f7 < fNumThreshold->GetNumber()/100.0) ? 1.0 : 0.0;
                else if (mode == 4) val = (iTotal>0) ? central[0]/iTotal : 0;
                else if (mode == 5) val = (iTotal>0) ? (central[0]+central[1])/iTotal : 0;

                hScan->Fill(x, y, val);
                if (hLeak) hLeak->Fill(x, y, iTotal > 0 ? 1.0 - val : 0.0);
                if (mode == 0 && x >= PP() && x <= 2*PP() && y >= PP()+PadRadius() && y <= 2*PP()-PadRadius()) {
                    sumReadX[ix] += val;
                    sum2ReadX[ix] += val*val;
                    double lval = iTotal > 0 ? 1.0 - val : 0.0;
                    sumLostX[ix] += lval;
                    sum2LostX[ix] += lval*lval;
                    cntProjX[ix]++;
                }
            }
            if (ix % 2 == 0) {
                cScan->cd(mode==0 ? 1 : 0); hScan->Draw("COLZ");
                if (hLeak) { cScan->cd(2); hLeak->Draw("COLZ"); }
                cScan->Update(); gSystem->ProcessEvents();
            }
        }
        cScan->cd(mode==0 ? 1 : 0); hScan->Draw("COLZ");
        if (hLeak) {
            cScan->cd(2); hLeak->Draw("COLZ");
            double sumLeak = 0; int nFilled = 0;
            for (int bx=1; bx<=hLeak->GetNbinsX(); bx++)
                for (int by=1; by<=hLeak->GetNbinsY(); by++) {
                    sumLeak += hLeak->GetBinContent(bx, by); nFilled++;
                }
            if (nFilled > 0)
                printf(">>> Mean signal leakage = %.2f %%  (average over central pixel surface)\n",
                       100.0 * sumLeak / nFilled);

            // Panels 3 & 4: X projection of fraction read and fraction lost (avg over Y)
            delete gROOT->FindObject("hProjReadX"); delete gROOT->FindObject("hProjLostX");
            TH1D *hProjRead = new TH1D("hProjReadX",
                "Fraction read vs X (avg over Y);X [#mum];Mean fraction",
                nBins, histMin-step/2., histMax+step/2.);
            TH1D *hProjLost = new TH1D("hProjLostX",
                "Fraction lost vs X (avg over Y);X [#mum];Mean fraction",
                nBins, histMin-step/2., histMax+step/2.);
            hProjRead->SetStats(0); hProjLost->SetStats(0);
            for (int ix = 0; ix < nBins; ix++) {
                int cnt = cntProjX[ix];
                if (cnt > 0) {
                    double mR = sumReadX[ix]/cnt;
                    double sR = sqrt(std::max(0.0, sum2ReadX[ix]/cnt - mR*mR));
                    hProjRead->SetBinContent(ix+1, mR);
                    hProjRead->SetBinError(ix+1, sR/sqrt((double)cnt));
                    double mL = sumLostX[ix]/cnt;
                    double sL = sqrt(std::max(0.0, sum2LostX[ix]/cnt - mL*mL));
                    hProjLost->SetBinContent(ix+1, mL);
                    hProjLost->SetBinError(ix+1, sL/sqrt((double)cnt));
                }
            }
            cScan->cd(3); gPad->SetGrid();
            hProjRead->SetMinimum(0); hProjRead->SetMaximum(1);
            hProjRead->GetXaxis()->SetRangeUser(PP(), 2*PP());
            hProjRead->SetMarkerStyle(20); hProjRead->SetMarkerSize(0.8);
            hProjRead->SetMarkerColor(kTeal-5); hProjRead->SetLineColor(kTeal-5);
            hProjRead->Draw("E1");

            cScan->cd(4); gPad->SetGrid();
            hProjLost->SetMinimum(0); hProjLost->SetMaximum(1);
            hProjLost->GetXaxis()->SetRangeUser(PP(), 2*PP());
            hProjLost->SetMarkerStyle(20); hProjLost->SetMarkerSize(0.8);
            hProjLost->SetMarkerColor(kOrange-3); hProjLost->SetLineColor(kOrange-3);
            hProjLost->Draw("E1");
        }
        cScan->Update(); cScan->RaiseWindow();
        SetStatus("  Ready  ", false);
    }

    void DoResolutionAnalysis() {
        SetStatus("  Computing...  ", true);
        TDecompSparse *solver = GetSolver(false);
        TRandom3 rand(0);

        double noise_rms = fNumNoise->GetNumber();
        int    nEvents   = (int)fNumRecoEvents->GetNumber();

        TCanvas *cRes = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cRes");
        if (!cRes) cRes = new TCanvas("cRes", "c3 - Resolution Analysis", 800, 800);
        cRes->cd(); cRes->Clear(); cRes->Divide(2,2);

        TGraph *gRec  = new TGraph(); gRec->SetMarkerStyle(20);  gRec->SetMarkerSize(0.6);  gRec->SetMarkerColor(kRed);

        delete gROOT->FindObject("hErrXMap");   delete gROOT->FindObject("hRes1D");
        delete gROOT->FindObject("gScatterX");  delete gROOT->FindObject("hSumE2Unb");
        delete gROOT->FindObject("hCntUnb");
        const int nMapBins = 30;
        TH2D *hSumE2 = new TH2D("hSumE2Unb", "", nMapBins, 0, PP(), nMapBins, 0, PP());
        TH2D *hCnt   = new TH2D("hCntUnb",   "", nMapBins, 0, PP(), nMapBins, 0, PP());
        TH2D *hErrX  = new TH2D("hErrXMap",
                                 Form("X Resolution Map (Noise=%.3f);X [#mum];Y [#mum];RMS #DeltaX [#mum]", noise_rms),
                                 nMapBins, 0, PP(), nMapBins, 0, PP());
        hErrX->SetStats(0);

        TH1D *hRes1D = new TH1D("hRes1D", "Residuals (Unbalance);X_{rec}-X_{true} [#mum];Entries", 60, -150, 150);

        TGraph *gScatterX = new TGraph(); gScatterX->SetName("gScatterX");
        gScatterX->SetTitle("X linearity (no noise);X_{true} [#mum];X_{reco} [#mum]");
        gScatterX->SetMarkerStyle(7); gScatterX->SetMarkerColor(kBlue+1);

        int count = 0, nSkipped = 0;
        std::cout << "--- Unbalance Resolution (" << nEvents << " random events, Noise=" << noise_rms << ") ---" << std::endl;

        for (int ev = 0; ev < nEvents; ev++) {
            double x = PP() + rand.Rndm() * PP();
            double y = PP() + rand.Rndm() * PP();

            int ii  = std::max(0, std::min(n-1, TMath::Nint(x/Pitch())));
            int jj  = std::max(0, std::min(n-1, TMath::Nint(y/Pitch())));
            int idx = ii*n + jj;

            // Skip events landing on an electrode pad (no sharing in the real sensor)
            if (fTypeMap[idx] > 0) { nSkipped++; continue; }

            TVectorD B_s(nTot); B_s.Zero(); B_s(idx) = 1.0;
            bool ok; TVectorD V_s = solver->Solve(B_s, ok);
            if (!ok) { nSkipped++; continue; }

            std::vector<double> cur4(17, 0.0); double iTot4 = 0;
            CollectCurrents(V_s, cur4, iTot4);
            { double ao = fNumAmpOffset->GetNumber() / Amplitude();
              cur4[6] += ao; cur4[7] += ao; cur4[10] += ao; cur4[11] += ao; }
            double amp = Amplitude();
            double A6 = cur4[6]*amp, A7 = cur4[7]*amp, A10 = cur4[10]*amp, A11 = cur4[11]*amp;

            // No-noise linearity
            double s4nl = A6+A7+A10+A11;
            double xr_nl = (s4nl > 1e-9) ? PP()*1.5 + ((A10+A11)-(A6+A7))/s4nl * PP()/2.0 : PP()*1.5;
            double yr_nl = (s4nl > 1e-9) ? PP()*1.5 + ((A7+A11)-(A6+A10))/s4nl * PP()/2.0 : PP()*1.5;
            gScatterX->SetPoint(gScatterX->GetN(), x - PP(), xr_nl - PP());

            // Noisy reconstruction
            double L   = fChkLandau->IsOn() ? rand.Landau(1.0, 0.3) : 1.0;
            double a6  = L*A6  + rand.Gaus(0, noise_rms);
            double a7  = L*A7  + rand.Gaus(0, noise_rms);
            double a10 = L*A10 + rand.Gaus(0, noise_rms);
            double a11 = L*A11 + rand.Gaus(0, noise_rms);
            double s4  = a6+a7+a10+a11;
            double xr  = (s4 > 1e-9) ? PP()*1.5 + ((a10+a11)-(a6+a7))/s4 * PP()/2.0 : xr_nl;
            double yr  = (s4 > 1e-9) ? PP()*1.5 + ((a7+a11)-(a6+a10))/s4 * PP()/2.0 : yr_nl;

            double errX = xr - x;
            hRes1D->Fill(errX);
            hSumE2->Fill(x - PP(), y - PP(), errX*errX);
            hCnt  ->Fill(x - PP(), y - PP());
            gRec ->SetPoint(count, xr - PP(), yr - PP());
            count++;

            if (ev % 1000 == 0) gSystem->ProcessEvents();
        }
        std::cout << "  Accepted: " << count << "  Skipped (on pad): " << nSkipped << std::endl;

        // Build RMS resolution map from accumulated sum-of-squares
        for (int bx = 1; bx <= nMapBins; bx++)
            for (int by = 1; by <= nMapBins; by++) {
                double cnt = hCnt->GetBinContent(bx, by);
                if (cnt > 0)
                    hErrX->SetBinContent(bx, by, std::sqrt(hSumE2->GetBinContent(bx, by) / cnt));
            }

        cRes->cd(1); gPad->SetGrid();
        TH1F *hr = cRes->cd(1)->DrawFrame(0, 0, PP(), PP());
        hr->SetTitle(Form("Distortion Grid (Noise=%.3f);X [#mum];Y [#mum]", noise_rms));
        gRec->Draw("P SAME");

        cRes->cd(2); hErrX->Draw("COLZ");

        cRes->cd(3);
        hRes1D->SetFillColor(kCyan-9);
        hRes1D->SetStats(1);
        hRes1D->Draw();
        {
            double xlo = hRes1D->GetXaxis()->GetXmin();
            double xhi = hRes1D->GetXaxis()->GetXmax();
            delete gROOT->FindObject("tStudRes");
            TF1 *tStud = new TF1("tStudRes",
                "[0]*TMath::Gamma(([3]+1.0)/2.0)"
                "/(TMath::Sqrt([3]*TMath::Pi())*TMath::Gamma([3]/2.0)*[2])"
                "*pow(1.0+1.0/[3]*pow((x-[1])/[2],2),-([3]+1.0)/2.0)",
                xlo, xhi);
            double peak = hRes1D->GetBinCenter(hRes1D->GetMaximumBin());
            double rms  = hRes1D->GetRMS();
            tStud->SetParameters(hRes1D->Integral("width"), peak, rms, 2.0);
            tStud->SetParLimits(2, 0.1, rms*5);   // sigma > 0
            tStud->SetParLimits(3, 1.0, 100.0);   // nu >= 1
            tStud->SetLineColor(kRed); tStud->SetLineWidth(2);
            hRes1D->Fit(tStud, "R");
            double sigma = tStud->GetParameter(2);
            double nu    = tStud->GetParameter(3);
            TLatex ltx; ltx.SetNDC(); ltx.SetTextSize(0.04);
            ltx.DrawLatex(0.15, 0.85, Form("#sigma = %.1f #mum", sigma));
            ltx.DrawLatex(0.15, 0.79, Form("#nu = %.1f", nu));
        }

        cRes->cd(4); gPad->SetGrid();
        TH1F *hrLin = cRes->cd(4)->DrawFrame(0, 0, PP(), PP());
        hrLin->SetTitle("X linearity (no noise);X_{true} [#mum];X_{reco} [#mum]");
        gScatterX->Draw("P SAME");
        TLine *lDiagC3 = new TLine(0, 0, PP(), PP());
        lDiagC3->SetLineColor(kRed); lDiagC3->SetLineStyle(2); lDiagC3->Draw();

        cRes->Update(); cRes->RaiseWindow();
        SetStatus("  Ready  ", false);
    }



void DoSharingTemplates() {
    SetStatus("  Computing...  ", true);
    TDecompSparse *solver = GetSolver(false);
    int    nx   = (int)fNumScanStep->GetNumber();
    int    ny   = nx;
    double step = PP() / nx;

    double min_ax = 0;
    double max_ax = PP();
    
    int pIDs[4] = {7, 11, 6, 10};

    for(int i=0; i<4; i++) {
        if(fTemplatesReady && fTemplates[i]) delete fTemplates[i];
        fTemplates[i] = new TH2D(Form("fTemp_%d", pIDs[i]),
                         Form("Sharing Template Pad %d (offset=%.3f);X [um];Y [um]", pIDs[i], fNumAmpOffset->GetNumber()),
                         nx, min_ax, max_ax, ny, min_ax, max_ax);
        fTemplates[i]->SetStats(0);
    }

    
    TCanvas *cT = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cT");
    if (!cT) cT = new TCanvas("cT", "c4 - Sharing Templates (Central Pixel)", 600, 540);
    cT->cd(); cT->Clear(); cT->Divide(2,2);
    
    std::cout << ">>> Computing Sharing Templates in the central pixel (500-1000 um)..." << std::endl;

    // Loop at bin centres of the central pixel [PP, 2*PP]
    for (double x = PP()+step/2.0; x <= PP()*2.0-step/2.0+1e-6; x += step) {
        for (double y = PP()+step/2.0; y <= PP()*2.0-step/2.0+1e-6; y += step) {
            
            TVectorD B(nTot); B.Zero();
            { int ii = std::max(0, std::min(n-1, TMath::Nint(x/Pitch())));
              int jj = std::max(0, std::min(n-1, TMath::Nint(y/Pitch())));
              B(ii*n + jj) = 1.0; }
            
            bool ok; 
            TVectorD V = solver->Solve(B, ok);
            
            std::vector<double> cur(17, 0.0);
            double iTot = 0;
            CollectCurrents(V, cur, iTot);

            // normalize to 1 over the 4 central pads
            double ampOffset = fNumAmpOffset->GetNumber() / Amplitude();
            for(int p=0; p<4; p++) cur[pIDs[p]] += ampOffset;
            double iTot4 = 0;
            for(int p=0; p<4; p++) {
	      iTot4 += cur[pIDs[p]];
            }

            for(int p=0; p<4; p++) {
	      fTemplates[p]->Fill(x - PP(), y - PP(), (iTot4 > 0) ? cur[pIDs[p]]/iTot4 : 0);
            }
        }
        gSystem->ProcessEvents(); // prevent GUI from freezing
        std::cout << "Progress: " << (int)((x - PP()) / PP() * 100.0) << "%" << std::endl;
    }
    
    fTemplatesReady = true;
    fBtnReco->SetEnabled(kTRUE);
    fBtnExport->SetEnabled(kTRUE);

    for(int i=0; i<4; i++) {
        cT->cd(i+1);
        fTemplates[i]->SetMinimum(0.0);
        fTemplates[i]->SetMaximum(1.0);
        fTemplates[i]->Draw("COLZ");
    }
    cT->Update();
    std::cout << ">>> Sharing Templates complete." << std::endl;
    SetStatus("  Ready  ", false);
}

 void DoLoadFromCache() {
    TGFileInfo fi;
    const char *types[] = { "ROOT files", "*.root", "All files", "*", nullptr, nullptr };
    fi.fFileTypes = types;
    fi.fIniDir    = StrDup("rsd_cache");
    TString origDir = gSystem->WorkingDirectory();
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    gSystem->ChangeDirectory(origDir.Data());
    if (!fi.fFilename) return;

    std::string path(fi.fFilename);
    auto slash = path.rfind('/');
    std::string fn = (slash == std::string::npos) ? path : path.substr(slash + 1);
    auto dot = fn.rfind(".root");
    if (dot != std::string::npos) fn = fn.substr(0, dot);

    // Extract number after key; character after key must start with digit, '-', or '.'
    auto extract = [&](const std::string &key, double &val) -> bool {
        size_t pos = fn.find(key);
        if (pos == std::string::npos) return false;
        pos += key.size();
        if (pos >= fn.size() || (!std::isdigit((unsigned char)fn[pos]) && fn[pos] != '-' && fn[pos] != '.')) return false;
        auto end = fn.find('_', pos);
        std::string tok = (end == std::string::npos) ? fn.substr(pos) : fn.substr(pos, end - pos);
        try { val = std::stod(tok); return true; } catch (...) { return false; }
    };

    // _gap must not be followed by 'P' (to distinguish from _gapP)
    auto extractGap = [&](double &val) -> bool {
        size_t pos = 0;
        while ((pos = fn.find("_gap", pos)) != std::string::npos) {
            size_t a = pos + 4;
            if (a < fn.size() && fn[a] != 'P') {
                auto end = fn.find('_', a);
                std::string tok = (end == std::string::npos) ? fn.substr(a) : fn.substr(a, end - a);
                try { val = std::stod(tok); return true; } catch (...) {}
            }
            pos++;
        }
        return false;
    };

    double pp = 500, rsh = 1000, rsp = 25, rsm = 25, gap = 0, gapp = 0, tr = 0, rin = 0, rpad = 60;
    if (!extract("_PP", pp) || !extract("_Rsh", rsh) || !extract("_tr", tr)
        || !extract("_Rin", rin) || !extract("_Rpad", rpad)) {
        std::cout << "ERROR: cannot parse cache filename: " << fn << std::endl;
        return;
    }

    bool isTrench = (tr > 0.5);
    if (!isTrench) {
        extract("_RsP", rsp);
        extract("_RsM", rsm);
        extractGap(gap);
        extract("_gapP", gapp);
    }

    fNumPixelPitch->SetNumber((int)pp);
    fNumRsheet->SetNumber(rsh);
    fNumRin->SetNumber(rin);
    fNumPadRadius->SetNumber(rpad);
    fChkTrench->SetState(isTrench ? kButtonDown : kButtonUp);
    if (!isTrench) {
        fNumRStripPad->SetNumber(rsp);
        fNumRStripMid->SetNumber(rsm);
        fNumStripGap->SetNumber(gap);
        fNumStripGapPad->SetNumber(gapp);
    }
    DoTrenchToggle();
    std::cout << ">>> Cache params loaded: PP=" << (int)pp << " Rsh=" << rsh
              << (isTrench ? " Trench" : Form(" RsP=%.0f RsM=%.0f gap=%.1f gapP=%.1f", rsp, rsm, gap, gapp))
              << " Rin=" << rin << " Rpad=" << rpad << std::endl;
    DoCalculateGrid();
}

 void DoBrowseConfig() {
    TGFileInfo fi;
    const char *types[] = { "Config files", "*.txt", "All files", "*", nullptr, nullptr };
    fi.fFileTypes = types;
    fi.fIniDir    = StrDup("ConfigFiles");
    TString origDir = gSystem->WorkingDirectory();
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    gSystem->ChangeDirectory(origDir.Data());
    if (fi.fFilename) fCfgFile->SetText(fi.fFilename);
}

 void DoSaveConfig() {
    gSystem->mkdir("ConfigFiles", kTRUE);
    std::ofstream f(fCfgFile->GetText());
    if (!f) { std::cout << "ERROR: cannot write " << fCfgFile->GetText() << std::endl; return; }
    f << "PP="          << fNumPixelPitch->GetNumber()  << "\n";
    f << "Rsheet="      << fNumRsheet->GetNumber()      << "\n";
    f << "RStripPad="   << fNumRStripPad->GetNumber()   << "\n";
    f << "RStripMid="   << fNumRStripMid->GetNumber()   << "\n";
    f << "StripGap="    << fNumStripGap->GetNumber()    << "\n";
    f << "StripGapPad=" << fNumStripGapPad->GetNumber() << "\n";
    f << "Rin="         << fNumRin->GetNumber()         << "\n";
    f << "PadRadius="   << fNumPadRadius->GetNumber()   << "\n";
    f << "Trench="      << (fChkTrench->IsOn() ? 1 : 0) << "\n";
    f << "Amplitude="   << fNumAmplitude->GetNumber()   << "\n";
    f << "Noise="       << fNumNoise->GetNumber()       << "\n";
    f << "AmpOffset="   << fNumAmpOffset->GetNumber()   << "\n";
    f << "ConstTerm="   << fNumConstTerm->GetNumber()   << "\n";
    f << "ScanStep="    << fNumScanStep->GetNumber()    << "\n";
    f << "Threshold="   << fNumThreshold->GetNumber()   << "\n";
    f << "Thick="       << fNumThick->GetNumber()       << "\n";
    f << "Dur="         << fNumDur->GetNumber()         << "\n";
    f << "RecoEvents="  << fNumRecoEvents->GetNumber()  << "\n";
    f << "Landau="      << (fChkLandau->IsOn() ? 1 : 0) << "\n";
    f.close();
    std::cout << ">>> Config saved to " << fCfgFile->GetText() << std::endl;
}

 void DoLoadConfig() {
    std::ifstream f(fCfgFile->GetText());
    if (!f) { std::cout << "ERROR: cannot open " << fCfgFile->GetText() << std::endl; return; }
    std::string key; double val;
    while (f >> key) {
        auto eq = key.find('=');
        if (eq == std::string::npos) continue;
        std::string name = key.substr(0, eq);
        val = std::stod(key.substr(eq + 1));
        if      (name == "PP")          fNumPixelPitch->SetNumber((int)val);
        else if (name == "Rsheet")      fNumRsheet->SetNumber(val);
        else if (name == "RStripPad")   fNumRStripPad->SetNumber(val);
        else if (name == "RStripMid")   fNumRStripMid->SetNumber(val);
        else if (name == "StripGap")    fNumStripGap->SetNumber(val);
        else if (name == "StripGapPad") fNumStripGapPad->SetNumber(val);
        else if (name == "Rin")         fNumRin->SetNumber(val);
        else if (name == "PadRadius")   fNumPadRadius->SetNumber(val);
        else if (name == "Trench")      fChkTrench->SetState(val ? kButtonDown : kButtonUp);
        else if (name == "Amplitude")   fNumAmplitude->SetNumber(val);
        else if (name == "Noise")       fNumNoise->SetNumber(val);
        else if (name == "AmpOffset")   fNumAmpOffset->SetNumber(val);
        else if (name == "ConstTerm")   fNumConstTerm->SetNumber(val);
        else if (name == "ScanStep")    fNumScanStep->SetNumber(val);
        else if (name == "Threshold")   fNumThreshold->SetNumber(val);
        else if (name == "Thick")       fNumThick->SetNumber(val);
        else if (name == "Dur")         fNumDur->SetNumber(val);
        else if (name == "RecoEvents")  fNumRecoEvents->SetNumber((int)val);
        else if (name == "Landau")      fChkLandau->SetState(val ? kButtonDown : kButtonUp);
    }
    f.close();
    DoTrenchToggle();  // sync enabled/disabled state of strip fields
    std::cout << ">>> Config loaded from " << fCfgFile->GetText() << std::endl;
}

 void DoExportTemplate() {
    if (!fTemplatesReady) {
        std::cout << "ERROR: Run Sharing Templates first!" << std::endl;
        return;
    }
    int nbinx = fTemplates[0]->GetNbinsX();
    int nbiny = fTemplates[0]->GetNbinsY();

    // Template axes are in [0, PP] (shifted from absolute [PP, 2PP])
    int lowlim = 0;
    int uplim  = (int)PP();

    // Electrode order: pIDs={7,11,6,10} → [0]=E7(top-left), [1]=E11(top-right),
    //                                      [2]=E6(bottom-left), [3]=E10(bottom-right)
    // Output order matching ChRSD[pix][0..3] (ch0,ch1,ch2,ch3): fT[0],fT[1],fT[3],fT[2]
    int outIdx[4] = {0, 1, 3, 2};

    TString key;
    if (fChkTrench->IsOn())
        key = TString::Format("PP%.0f_Rsh%.0f_tr1_Rin%.2f_Rpad%.0f",
            PP(), fNumRsheet->GetNumber(), fNumRin->GetNumber(), PadRadius());
    else
        key = TString::Format("PP%.0f_Rsh%.0f_RsP%.0f_RsM%.0f_gap%.1f_gapP%.1f_tr0_Rin%.2f_Rpad%.0f",
            PP(), fNumRsheet->GetNumber(), fNumRStripPad->GetNumber(),
            fNumRStripMid->GetNumber(), fNumStripGap->GetNumber(),
            fNumStripGapPad->GetNumber(), fNumRin->GetNumber(), PadRadius());
    double ampOff = fNumAmpOffset->GetNumber();
    TString fname = TString("DCRSD_Template_") + key
                    + TString::Format("_off%.3f", ampOff) + ".txt";

    std::ofstream f(fname.Data());
    f << "DCRSD Sharing Template export " << key << "\n";
    f << lowlim << " " << uplim << " " << nbinx << " " << nbiny << "\n";
    for (int bx = 1; bx <= nbinx; bx++)
        for (int by = 1; by <= nbiny; by++) {
            f << bx << " " << by;
            for (int p = 0; p < 4; p++)
                f << " " << fTemplates[outIdx[p]]->GetBinContent(bx, by) * 100.0;
            for (int p = 0; p < 4; p++)
                f << " " << 0.0;  // delays not computed by DCRSD
            f << "\n";
        }
    f.close();
    std::cout << ">>> Template written to " << fname.Data()
              << " (" << nbinx << "x" << nbiny << " bins)" << std::endl;
}

 void DoTemplateReconstruction() {
    SetStatus("  Computing...  ", true);
    if (!fTemplatesReady) {
        std::cout << "ERROR: Generate Sharing Templates first!" << std::endl;
        return;
    }

    TDecompSparse *solver = GetSolver(false);
    TRandom3 rand(0);
    double noise = fNumNoise->GetNumber();
    
    int nEvents = (int)fNumRecoEvents->GetNumber();

    TCanvas *cReco = new TCanvas("cReco", "c5 - Template Reconstruction Analysis", 1260, 840);
    cReco->Divide(3,2);

    delete gROOT->FindObject("hResMapReco"); delete gROOT->FindObject("hErrXReco");
    delete gROOT->FindObject("hSumE2Reco"); delete gROOT->FindObject("hCntReco");
    delete gROOT->FindObject("hProjXReco"); delete gROOT->FindObject("gLinStrip");
    delete gROOT->FindObject("hRecMapReco");
    TH2D *hSumE2 = new TH2D("hSumE2Reco", "", 30, 0, PP(), 30, 0, PP());
    TH2D *hCnt   = new TH2D("hCntReco",   "", 30, 0, PP(), 30, 0, PP());
    TH2D *hResMap = new TH2D("hResMapReco",
                              Form("X Resolution Map (Noise=%.1f);X True [#mum];Y True [#mum];RMS #DeltaX [#mum]", noise),
                              30, 0, PP(), 30, 0, PP());

    // global X residual histogram
    TH1D *hErrX = new TH1D("hErrXReco", Form("X Resolution (Noise=%.1f);#Delta X [#mum];Events", noise), 100, -200., 200.);

    TGraph *gLin = new TGraph(); gLin->SetName("gLinStrip");
    gLin->SetTitle("X linearity;X_{true} [#mum];X_{reco} [#mum]");
    gLin->SetMarkerStyle(7); gLin->SetMarkerColor(kBlue+1);

    TH2D *hRecMap = new TH2D("hRecMapReco",
        Form("Reconstructed Position Density (Noise=%.1f);X_{reco} [#mum];Y_{reco} [#mum];Entries", noise),
        30, 0, PP(), 30, 0, PP());
    hRecMap->SetStats(0);

    std::cout << ">>> Running " << nEvents << " random events..." << std::endl;

    int nSkippedReco = 0;
    for (int ev = 0; ev < nEvents; ev++) {
        double x_true = PP() + rand.Rndm() * PP();
        double y_true = PP() + rand.Rndm() * PP();

        int ii = std::max(0, std::min(n-1, TMath::Nint(x_true/Pitch())));
        int jj = std::max(0, std::min(n-1, TMath::Nint(y_true/Pitch())));
        if (fTypeMap[ii*n + jj] > 0) { nSkippedReco++; continue; }

        TVectorD B(nTot); B.Zero();
        B(ii*n + jj) = 1.0;
        bool ok; TVectorD V = solver->Solve(B, ok);

        std::vector<double> cur(17, 0.0); double iTot = 0;
        CollectCurrents(V, cur, iTot);

        int pIDs[4] = {7, 11, 6, 10};
        double ampOffset = fNumAmpOffset->GetNumber() / Amplitude();
        for(int p=0; p<4; p++) cur[pIDs[p]] += ampOffset;
        double iTot4 = 0;
        for(int p=0; p<4; p++) iTot4 += cur[pIDs[p]];
        double L = fChkLandau->IsOn() ? rand.Landau(1.0, 0.3) : 1.0;
        double amp = Amplitude();
        double A_noisy[4];
        double sum_noisy = 0;
        for(int p=0; p<4; p++) {
            double sig = (iTot4 > 0) ? cur[pIDs[p]]/iTot4 : 0;
            A_noisy[p] = amp * sig + rand.Gaus(0, noise) / L;
            sum_noisy += A_noisy[p];
        }

        for(int p=0; p<4; p++) A_noisy[p] /= sum_noisy;

        // template matching: find minimum chi2
        double min_chi2 = 1e9;
        int bx_best = 1, by_best = 1;

        int nBinsX = fTemplates[0]->GetNbinsX();
        int nBinsY = fTemplates[0]->GetNbinsY();

        for(int bx = 1; bx <= nBinsX; bx++) {
            for(int by = 1; by <= nBinsY; by++) {
                double chi2 = 0;
                for(int p=0; p<4; p++) {
                    double t_val = fTemplates[p]->GetBinContent(bx, by);
                    chi2 += pow(A_noisy[p] - t_val, 2);
                }
                if(chi2 < min_chi2) {
                    min_chi2 = chi2;
                    bx_best = bx; by_best = by;
                }
            }
        }

        // parabolic interpolation around chi2 minimum; removes scan-direction discretisation bias
        auto chi2_at = [&](int bx, int by) {
            double c = 0;
            for(int p=0; p<4; p++) {
                double t = fTemplates[p]->GetBinContent(bx, by);
                c += pow(A_noisy[p] - t, 2);
            }
            return c;
        };
        double bw_x = fTemplates[0]->GetXaxis()->GetBinWidth(bx_best);
        double bw_y = fTemplates[0]->GetYaxis()->GetBinWidth(by_best);
        double c0   = min_chi2;
        double cx_m = (bx_best > 1)      ? chi2_at(bx_best-1, by_best) : c0;
        double cx_p = (bx_best < nBinsX) ? chi2_at(bx_best+1, by_best) : c0;
        double cy_m = (by_best > 1)      ? chi2_at(bx_best, by_best-1) : c0;
        double cy_p = (by_best < nBinsY) ? chi2_at(bx_best, by_best+1) : c0;
        double dx = 0, dy = 0;
        double denom_x = cx_m - 2*c0 + cx_p;
        double denom_y = cy_m - 2*c0 + cy_p;
        if(denom_x > 0) dx = 0.5*(cx_m - cx_p) / denom_x;
        if(denom_y > 0) dy = 0.5*(cy_m - cy_p) / denom_y;
        // Clamp to ±0.5 bin
        dx = std::max(-0.5, std::min(0.5, dx));
        dy = std::max(-0.5, std::min(0.5, dy));
        double x_rec = fTemplates[0]->GetXaxis()->GetBinCenter(bx_best) + dx * bw_x;
        double y_rec = fTemplates[0]->GetYaxis()->GetBinCenter(by_best) + dy * bw_y;

        // x_rec/y_rec in [0,PP] (shifted templates); x_true/y_true in [PP,2PP]
        // Constant term: smear the reference position to simulate telescope misalignment.
        // Total resolution in quadrature: σ_tot² = σ_intrinsic² + σ_const²
        double sigma_ct = fNumConstTerm->GetNumber();
        double ref_smear_x = (sigma_ct > 0) ? rand.Gaus(0, sigma_ct) : 0.0;
        double ref_smear_y = (sigma_ct > 0) ? rand.Gaus(0, sigma_ct) : 0.0;
        double errX = x_rec - (x_true - PP()) - ref_smear_x;
        double errY = y_rec - (y_true - PP()) - ref_smear_y;
        hSumE2->Fill(x_true - PP(), y_true - PP(), errX*errX);
        hCnt->Fill(x_true - PP(), y_true - PP());
        hErrX->Fill(errX);
        gLin->SetPoint(gLin->GetN(), x_true - PP(), x_rec);
        hRecMap->Fill(x_rec, y_rec);

        if (ev % 500 == 0) std::cout << "Processed " << ev << " / " << nEvents << " events..." << std::endl;
    }
    std::cout << "  Accepted: " << (nEvents - nSkippedReco) << "  Skipped (on pad): " << nSkippedReco << std::endl;

    // per-bin RMS: sqrt(sum(errX²)/N)
    for (int bx = 1; bx <= hResMap->GetNbinsX(); bx++)
        for (int by = 1; by <= hResMap->GetNbinsY(); by++) {
            double cnt = hCnt->GetBinContent(bx, by);
            if (cnt > 0)
                hResMap->SetBinContent(bx, by, sqrt(hSumE2->GetBinContent(bx, by) / cnt));
        }
    hResMap->SetStats(0);
    hResMap->SetMinimum(0);
    hResMap->SetMaximum(30);


    // X projection: RMS(errX) vs X_true, integrated over Y, with error bars = RMS/sqrt(2N)
    TH1D *hProjX = new TH1D("hProjXReco",
        Form("X Resolution vs X_{true} (Noise=%.1f);X_{true} [#mum];RMS #DeltaX [#mum]", noise),
        hSumE2->GetNbinsX(), 0, PP());
    hProjX->SetStats(0);
    for (int bx = 1; bx <= hSumE2->GetNbinsX(); bx++) {
        double sumE2_x = 0, cnt_x = 0;
        for (int by = 1; by <= hSumE2->GetNbinsY(); by++) {
            sumE2_x += hSumE2->GetBinContent(bx, by);
            cnt_x   += hCnt->GetBinContent(bx, by);
        }
        if (cnt_x > 0) {
            double rms = sqrt(sumE2_x / cnt_x);
            hProjX->SetBinContent(bx, rms);
            hProjX->SetBinError(bx, rms / sqrt(2.0 * cnt_x));
        }
    }

    // --- DISEGNO (3x2: P1=distortion grid, P2=ResMap, P3=ProjX, P4=gauss, P5=t-Student, P6=linearita') ---

    // Panel 1: 2D density map of reconstructed positions
    cReco->cd(1);
    gPad->SetRightMargin(0.15);
    hRecMap->Draw("COLZ");
    hRecMap->GetXaxis()->SetNdivisions(5, kFALSE);

    // Panel 2: 2D RMS resolution map
    cReco->cd(2);
    gPad->SetRightMargin(0.15);
    hResMap->Draw("COLZ");
    hResMap->GetXaxis()->SetNdivisions(5, kFALSE);

    // Panel 3: RMS(errX) vs X_true with error bars = RMS/sqrt(2N)
    cReco->cd(3); gPad->SetGrid();
    hProjX->SetMinimum(0);
    hProjX->SetMarkerStyle(20); hProjX->SetMarkerSize(0.8);
    hProjX->SetMarkerColor(kTeal-5); hProjX->SetLineColor(kTeal-5);
    hProjX->SetFillColor(0);
    hProjX->Draw("E1");
    hProjX->GetXaxis()->SetNdivisions(5, kFALSE);

    // Panel 4: 1D residuals with gaussian fit
    cReco->cd(4);
    hErrX->SetFillStyle(0);
    hErrX->SetStats(1);
    hErrX->Draw();
    hErrX->Fit("gaus", "Q");
    TF1 *fitG = hErrX->GetFunction("gaus");
    if(fitG) {
        TLatex txt; txt.SetTextSize(0.05);
        txt.DrawLatex(-PP()/6.25, hErrX->GetMaximum()*0.85, Form("#sigma = %.1f #mum", fitG->GetParameter(2)));
    }

    // Panel 5: t-Student fit on residuals clone
    cReco->cd(5);
    delete gROOT->FindObject("hErrXStudReco");
    TH1D *hErrXS = (TH1D*)hErrX->Clone("hErrXStudReco");
    hErrXS->SetTitle(Form("X Resolution - t-Student (Noise=%.1f);#Delta X [#mum];Events", noise));
    hErrXS->SetFillStyle(0);
    hErrXS->SetStats(1);
    hErrXS->Draw();
    {
        double xLo = hErrXS->GetXaxis()->GetXmin();
        double xHi = hErrXS->GetXaxis()->GetXmax();
        double sig0  = fitG ? fitG->GetParameter(2) : PP()/20.0;
        double norm0 = hErrXS->Integral("width");
        TF1 *fitT = new TF1("fStudReco",
            "[0]/[2] * TMath::Student((x-[1])/[2], [3])",
            xLo, xHi);
        fitT->SetParNames("Norm", "Mean", "Scale", "#nu");
        fitT->SetParameters(norm0, 0.0, sig0, 3.0);
        fitT->SetParLimits(2, 1e-3, PP());
        fitT->SetParLimits(3, 0.5, 200.0);
        fitT->SetLineColor(kRed+1); fitT->SetLineWidth(2);
        hErrXS->Fit(fitT, "RQ");
        fitT->Draw("SAME");
        double nu_fit   = fitT->GetParameter(3);
        double sig_fit  = fitT->GetParameter(2);
        double enu_fit  = fitT->GetParError(3);
        double esig_fit = fitT->GetParError(2);
        printf("  t-Student fit:  sigma_scale = %.2f +/- %.2f um,  nu = %.2f +/- %.2f\n",
               sig_fit, esig_fit, nu_fit, enu_fit);
        if (nu_fit > 2.0)
            printf("  sigma_RMS = %.2f um\n", sig_fit * sqrt(nu_fit / (nu_fit - 2.0)));
    }

    // Panel 6: X linearity
    cReco->cd(6); gPad->SetGrid();
    TH1F *hfLin5 = cReco->cd(6)->DrawFrame(0, 0, PP(), PP());
    hfLin5->SetTitle("X linearity;X_{true} [#mum];X_{reco} [#mum]");
    gLin->Draw("P SAME");
    TLine *lDiagC5 = new TLine(0, 0, PP(), PP());
    lDiagC5->SetLineColor(kRed); lDiagC5->SetLineStyle(2); lDiagC5->Draw();

    cReco->Update();
    std::cout << ">>> Analysis complete." << std::endl;
    SetStatus("  Ready  ", false);
}

void DoResistanceMatrix() {
    SetStatus("  Computing...  ", true);
    std::cout << "\n>>> Computing inter-electrode resistance (2-terminal, all others floating)..." << std::endl;

    // Pure 2-terminal measurement: inject 1A into P11, ground only the sink pad.
    // All other 14 pads floating. R = V(P11) - V(sink).
    // Pairs: P11(2,2) <-> P7(1,2) [side], P11(2,2) <-> P6(1,1) [diagonal]

    const int di[] = {1,-1,0,0}, dj[] = {0,0,1,-1};

    auto padNode = [&](int r, int c) {
        int pi = std::max(0, std::min(n-1, TMath::Nint(r * PP() / Pitch())));
        int pj = std::max(0, std::min(n-1, TMath::Nint(c * PP() / Pitch())));
        return pi*n + pj;
    };

    int nodeP11 = padNode(2, 2);
    int nodeP7  = padNode(1, 2);  // side neighbor
    int nodeP6  = padNode(1, 1);  // diagonal neighbor

    // Free KCL matrix: no special pad treatment, all nodes follow pure KCL
    TMatrixDSparse Mfree(nTot, nTot);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i*n + j;
            double sum_g = 0;
            for (int kk = 0; kk < 4; kk++) {
                int ni2 = i+di[kk], nj2 = j+dj[kk];
                if (ni2>=0 && ni2<n && nj2>=0 && nj2<n) {
                    double cond = GetEdgeConductance(i, j, di[kk], dj[kk]);
                    sum_g += cond;
                    Mfree(idx, ni2*n+nj2) = -cond;
                }
            }
            Mfree(idx, idx) = sum_g;
        }
    }

    // Kelvin measurement: inject +1A into P11, extract -1A at sink.
    // System is non-singular thanks to a weak voltage reference (1 nS) at a corner node.
    auto addWeakRef = [&](TMatrixDSparse &M) { M(0, 0) += 1e-9; };

    double R_side = 0, R_diag = 0;

    // --- R_side: P11 ↔ P7 ---
    {
        TMatrixDSparse M(Mfree); addWeakRef(M);
        TDecompSparse sol; sol.SetMatrix(M);
        if (!sol.Decompose()) { std::cerr << ">>> ERROR R_side" << std::endl; SetStatus("  Ready  ", false); return; }
        TVectorD B(nTot); B.Zero(); B(nodeP11) = +1.0; B(nodeP7) = -1.0;
        bool ok; TVectorD V = sol.Solve(B, ok);
        R_side = V(nodeP11) - V(nodeP7);
        printf("  R_side  P11(2,2) <-> P7(1,2)  = %.1f Ohm\n", R_side);
    }

    // --- R_diag: P11 ↔ P6 ---
    {
        TMatrixDSparse M(Mfree); addWeakRef(M);
        TDecompSparse sol; sol.SetMatrix(M);
        if (!sol.Decompose()) { std::cerr << ">>> ERROR R_diag" << std::endl; SetStatus("  Ready  ", false); return; }
        TVectorD B(nTot); B.Zero(); B(nodeP11) = +1.0; B(nodeP6) = -1.0;
        bool ok; TVectorD V = sol.Solve(B, ok);
        R_diag = V(nodeP11) - V(nodeP6);
        printf("  R_diag  P11(2,2) <-> P6(1,1)  = %.1f Ohm\n", R_diag);
    }

    // Bar chart
    if (!cRmat) cRmat = new TCanvas("cRmat", "c6 - Inter-electrode Resistance - Central Pixel", 360, 300);
    cRmat->cd(); cRmat->Clear();
    delete gROOT->FindObject("hRvec");
    TH1D *hRvec = new TH1D("hRvec",
        Form("Inter-electrode Resistance [#Omega]  (PP=%.0f #mum, R_{sh}=%.0f, R_{sP}=%.0f, R_{sM}=%.0f);Pair;R [#Omega]",
             PP(), fNumRsheet->GetNumber(), fNumRStripPad->GetNumber(), fNumRStripMid->GetNumber()),
        2, 0, 2);
    hRvec->GetXaxis()->SetBinLabel(1, "R_{side}  P11#leftrightarrowP7");
    hRvec->GetXaxis()->SetBinLabel(2, "R_{diag}  P11#leftrightarrowP6");
    hRvec->SetBinContent(1, R_side);
    hRvec->SetBinContent(2, R_diag);
    hRvec->SetFillColor(kAzure+1);
    hRvec->SetMinimum(0);
    hRvec->GetXaxis()->SetLabelSize(0.05);
    gStyle->SetPaintTextFormat(".0f");
    hRvec->Draw("BAR TEXT0");
    cRmat->Update();
    std::cout << ">>> Inter-electrode resistance computed." << std::endl;
    SetStatus("  Ready  ", false);
}




};

void DCRSD_SurfaceSimulator() {
    new RSDStripPanel();
}
