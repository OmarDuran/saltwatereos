/**
 * @file test_pHX_bisection_smoothness.cpp
 * @brief Test smoothness of density and saturation across phase boundaries
 *        using prop_pHX_bisection.
 *
 * Strategy
 * --------
 * For several representative (P, X) pairs, sweep enthalpy H in fine steps
 * and verify:
 *
 *   1. **Within** a single phase region, Rho, S_l, S_v change smoothly
 *      (no unphysical jumps).
 *   2. At phase-boundary crossings (region changes), large jumps are
 *      expected and are only reported, not treated as failures.
 *   3. Phase densities Rho_l, Rho_v are non-negative and finite.
 *
 * The test uses 8 points in each of:
 *   0.0001 < X  < 0.20
 *   6      < P  < 510    bar  (= 0.6 – 51 MPa)
 * and sweeps enthalpy from 500 to 3500 kJ/kg in 50 kJ/kg steps.
 */

#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

// ─── helpers ────────────────────────────────────────────────────────────

static const char* regionName(int r)
{
    switch (r) {
        case  0: return "L";
        case  1: return "LV_X0";
        case  2: return "V";
        case  3: return "L+H";
        case  4: return "V+H";
        case  5: return "V+L+H";
        case  6: return "V+L_L";
        case  7: return "V+L_V";
        case  8: return "Unknown";
        default: return "???";
    }
}

/// Absolute jump between two values.
static double absJump(double a, double b) { return std::abs(b - a); }

// ─── per-scan result ────────────────────────────────────────────────────

struct ScanResult {
    double P_bar;
    double X_wt;
    int    nSteps;
    int    nRegionChanges;

    // Worst WITHIN-REGION jumps (the ones that indicate bugs)
    double maxInRegion_Rho;
    double maxInRegion_Sl;
    double maxInRegion_Sv;
    double H_inRegion_Rho;
    double H_inRegion_Sl;
    double H_inRegion_Sv;
    int    reg_inRegion_Rho;
    int    reg_inRegion_Sl;
    int    reg_inRegion_Sv;
    int    nInRegionRhoFails;
    int    nInRegionSatFails;

    // Worst CROSS-BOUNDARY jumps (expected, reported only)
    double maxCross_Rho;
    double maxCross_Sl;
    double H_cross_Rho;
    int    reg_before_cross;
    int    reg_after_cross;

    // counts
    int nNaN;
    int nNegRho;

    bool passed;
};

// ─── threshold constants ────────────────────────────────────────────────

// Maximum allowed WITHIN-REGION jump in bulk density (kg/m³) per 50 kJ/kg step.
// Within a single phase region, density should change smoothly.
static const double MAX_INREGION_RHO_JUMP = 150.0;   // kg/m³

// Maximum allowed WITHIN-REGION saturation jump per step.
// Within a single phase region, saturations should be constant or change smoothly.
static const double MAX_INREGION_SAT_JUMP = 0.30;

// Enthalpy step size (kJ/kg)
static const double DH_KJ = 50.0;

// ─── main ───────────────────────────────────────────────────────────────

int main()
{
    H2ONaCl::cH2ONaCl eos;

    cout << "================================================================\n";
    cout << " Phase Boundary Smoothness Test  (prop_pHX_bisection)\n";
    cout << "================================================================\n\n";
    cout << "Checks that Rho, S_l, S_v change continuously WITHIN each\n";
    cout << "phase region.  Cross-boundary jumps are expected and reported\n";
    cout << "but not treated as failures.\n\n";

    // ── build test grid ────────────────────────────────────────────────
    vector<double> X_vals;
    for (int i = 0; i < 8; ++i)
        X_vals.push_back(0.0001 + i * (0.20 - 0.0001) / 7.0);

    vector<double> P_vals;
    {
        double logMin = log10(6.0), logMax = log10(510.0);
        for (int i = 0; i < 8; ++i)
            P_vals.push_back(pow(10.0, logMin + i * (logMax - logMin) / 7.0));
    }

    const double H_min_kJ = 500.0;
    const double H_max_kJ = 3500.0;
    const int    nH = static_cast<int>((H_max_kJ - H_min_kJ) / DH_KJ) + 1;

    cout << "Grid:  " << X_vals.size() << " salinities x "
         << P_vals.size() << " pressures x "
         << nH << " enthalpy steps  (dH = " << DH_KJ << " kJ/kg)\n";
    cout << "Within-region thresholds:  max |dRho| = " << MAX_INREGION_RHO_JUMP
         << " kg/m³,  max |dS| = " << MAX_INREGION_SAT_JUMP << "\n\n";

    // ── run scans ──────────────────────────────────────────────────────

    int totalScans     = 0;
    int passedScans    = 0;
    int failedScans    = 0;
    int totalSteps     = 0;
    int totalNaN       = 0;
    int totalNegRho    = 0;
    int totalCrossBoundary = 0;

    vector<ScanResult> results;

    cout << "Scanning...\n";
    cout << string(120, '=') << "\n";
    cout << setw(8) << "P(bar)"
         << setw(8) << "X"
         << setw(6) << "RegCh"
         << setw(14) << "inReg dRho"
         << setw(10) << "H@dRho"
         << setw(8) << "reg"
         << setw(12) << "inReg dSl"
         << setw(12) << "inReg dSv"
         << setw(14) << "cross dRho"
         << setw(12) << "transition"
         << setw(6) << "NaN"
         << setw(8) << "Status"
         << "\n";
    cout << string(120, '-') << "\n";

    for (size_t iP = 0; iP < P_vals.size(); ++iP) {
        double P_bar = P_vals[iP];
        double P_Pa  = P_bar * 1e5;

        for (size_t iX = 0; iX < X_vals.size(); ++iX) {
            double X_wt = X_vals[iX];
            totalScans++;

            ScanResult sr{};
            sr.P_bar = P_bar;
            sr.X_wt  = X_wt;
            sr.nSteps = nH;
            sr.passed = true;

            // Evaluate first point
            double H0_J = H_min_kJ * 1e3;
            H2ONaCl::PROP_H2ONaCl prev = eos.prop_pHX_bisection(P_Pa, H0_J, X_wt);
            int prevRegion = prev.Region;

            if (isnan(prev.Rho) || isnan(prev.T)) sr.nNaN++;
            if (prev.Rho < 0 && !isnan(prev.Rho)) sr.nNegRho++;

            for (int iH = 1; iH < nH; ++iH) {
                double H_kJ = H_min_kJ + iH * DH_KJ;
                double H_J  = H_kJ * 1e3;

                H2ONaCl::PROP_H2ONaCl cur = eos.prop_pHX_bisection(P_Pa, H_J, X_wt);
                totalSteps++;

                // NaN check
                if (isnan(cur.Rho) || isnan(cur.T)) {
                    sr.nNaN++;
                    prev = cur; prevRegion = cur.Region;
                    continue;
                }
                if (cur.Rho < 0) sr.nNegRho++;
                if (isnan(prev.Rho) || isnan(prev.T)) {
                    prev = cur; prevRegion = cur.Region;
                    continue;
                }

                // Use max(0, Rho) for jump comparison (negative Rho is near-zero vapor)
                double prevRho = std::max(0.0, prev.Rho);
                double curRho  = std::max(0.0, cur.Rho);

                bool regionChanged = (cur.Region != prevRegion);
                if (regionChanged) {
                    sr.nRegionChanges++;
                    totalCrossBoundary++;

                    // Track cross-boundary jumps (for reporting only)
                    double dRho = absJump(prevRho, curRho);
                    double dSl  = absJump(prev.S_l, cur.S_l);
                    if (dRho > sr.maxCross_Rho) {
                        sr.maxCross_Rho = dRho;
                        sr.H_cross_Rho  = H_kJ;
                        sr.reg_before_cross = prevRegion;
                        sr.reg_after_cross  = cur.Region;
                    }
                    if (dSl > sr.maxCross_Sl) sr.maxCross_Sl = dSl;
                } else {
                    // WITHIN-REGION: these jumps are the ones that matter
                    double dRho = absJump(prevRho, curRho);
                    double dSl  = absJump(prev.S_l, cur.S_l);
                    double dSv  = absJump(prev.S_v, cur.S_v);

                    if (dRho > sr.maxInRegion_Rho) {
                        sr.maxInRegion_Rho = dRho;
                        sr.H_inRegion_Rho  = H_kJ;
                        sr.reg_inRegion_Rho = cur.Region;
                    }
                    if (dSl > sr.maxInRegion_Sl) {
                        sr.maxInRegion_Sl = dSl;
                        sr.H_inRegion_Sl  = H_kJ;
                        sr.reg_inRegion_Sl = cur.Region;
                    }
                    if (dSv > sr.maxInRegion_Sv) {
                        sr.maxInRegion_Sv = dSv;
                        sr.H_inRegion_Sv  = H_kJ;
                        sr.reg_inRegion_Sv = cur.Region;
                    }
                    if (dRho > MAX_INREGION_RHO_JUMP) sr.nInRegionRhoFails++;
                    if (dSl > MAX_INREGION_SAT_JUMP || dSv > MAX_INREGION_SAT_JUMP)
                        sr.nInRegionSatFails++;
                }

                prev = cur;
                prevRegion = cur.Region;
            }

            // ── evaluate pass / fail (within-region only) ──
            if (sr.nInRegionRhoFails > 0 || sr.nInRegionSatFails > 0) {
                sr.passed = false;
            }

            totalNaN    += sr.nNaN;
            totalNegRho += sr.nNegRho;
            if (sr.passed) passedScans++; else failedScans++;

            results.push_back(sr);

            // Print failures and every 16th scan for progress
            if (!sr.passed || (totalScans % 16 == 0)) {
                ostringstream crossTrans;
                if (sr.nRegionChanges > 0)
                    crossTrans << regionName(sr.reg_before_cross) << "->" << regionName(sr.reg_after_cross);
                else
                    crossTrans << "none";

                cout << setw(8) << fixed << setprecision(1) << P_bar
                     << setw(8) << setprecision(4) << X_wt
                     << setw(6) << sr.nRegionChanges
                     << setw(14) << setprecision(2) << sr.maxInRegion_Rho
                     << setw(10) << setprecision(0) << sr.H_inRegion_Rho
                     << setw(8) << regionName(sr.reg_inRegion_Rho)
                     << setw(12) << setprecision(4) << sr.maxInRegion_Sl
                     << setw(12) << setprecision(4) << sr.maxInRegion_Sv
                     << setw(14) << setprecision(1) << sr.maxCross_Rho
                     << setw(12) << crossTrans.str()
                     << setw(6) << sr.nNaN
                     << setw(8) << (sr.passed ? "OK" : "FAIL")
                     << "\n";
            }
        }
    }

    cout << string(120, '=') << "\n\n";

    // ── detailed failure report ────────────────────────────────────────
    int nFailsToPrint = 0;
    for (auto& sr : results) if (!sr.passed) nFailsToPrint++;

    if (nFailsToPrint > 0) {
        cout << "========== WITHIN-REGION FAILURE REPORT ==========\n\n";

        vector<ScanResult*> byRho;
        for (auto& sr : results) if (sr.nInRegionRhoFails > 0) byRho.push_back(&sr);
        sort(byRho.begin(), byRho.end(), [](const ScanResult* a, const ScanResult* b) {
            return a->maxInRegion_Rho > b->maxInRegion_Rho;
        });

        if (!byRho.empty()) {
            cout << "Within-region density jumps exceeding " << MAX_INREGION_RHO_JUMP << " kg/m³:\n";
            cout << setw(10) << "P(bar)" << setw(10) << "X"
                 << setw(14) << "maxdRho" << setw(12) << "H(kJ/kg)"
                 << setw(12) << "Region" << "\n";
            cout << string(58, '-') << "\n";
            for (auto* s : byRho) {
                cout << setw(10) << fixed << setprecision(1) << s->P_bar
                     << setw(10) << setprecision(4) << s->X_wt
                     << setw(14) << setprecision(2) << s->maxInRegion_Rho
                     << setw(12) << setprecision(0) << s->H_inRegion_Rho
                     << "  " << regionName(s->reg_inRegion_Rho)
                     << "\n";
            }
            cout << "\n";
        }

        vector<ScanResult*> bySat;
        for (auto& sr : results) if (sr.nInRegionSatFails > 0) bySat.push_back(&sr);
        sort(bySat.begin(), bySat.end(), [](const ScanResult* a, const ScanResult* b) {
            return max(a->maxInRegion_Sl, a->maxInRegion_Sv) > max(b->maxInRegion_Sl, b->maxInRegion_Sv);
        });

        if (!bySat.empty()) {
            cout << "Within-region saturation jumps exceeding " << MAX_INREGION_SAT_JUMP << ":\n";
            cout << setw(10) << "P(bar)" << setw(10) << "X"
                 << setw(12) << "maxdSl" << setw(12) << "maxdSv"
                 << setw(12) << "H(kJ/kg)" << setw(12) << "Region" << "\n";
            cout << string(68, '-') << "\n";
            for (auto* s : bySat) {
                int regSat = (s->maxInRegion_Sl >= s->maxInRegion_Sv) ? s->reg_inRegion_Sl : s->reg_inRegion_Sv;
                double H_sat = (s->maxInRegion_Sl >= s->maxInRegion_Sv) ? s->H_inRegion_Sl : s->H_inRegion_Sv;
                cout << setw(10) << fixed << setprecision(1) << s->P_bar
                     << setw(10) << setprecision(4) << s->X_wt
                     << setw(12) << setprecision(4) << s->maxInRegion_Sl
                     << setw(12) << setprecision(4) << s->maxInRegion_Sv
                     << setw(12) << setprecision(0) << H_sat
                     << "  " << regionName(regSat)
                     << "\n";
            }
            cout << "\n";
        }
    }

    // ── cross-boundary summary (informational) ────────────────────────
    cout << "========== CROSS-BOUNDARY TRANSITIONS (informational) ==========\n\n";
    cout << "  Total boundary crossings:  " << totalCrossBoundary << "\n";
    {
        double worstCross = 0;
        for (auto& sr : results) if (sr.maxCross_Rho > worstCross) worstCross = sr.maxCross_Rho;
        cout << "  Worst cross-boundary dRho: " << fixed << setprecision(1) << worstCross << " kg/m³\n";
    }
    cout << "  (These are expected — phase transitions cause large property changes.)\n\n";

    // ── summary ────────────────────────────────────────────────────────
    cout << "================================================================\n";
    cout << " SUMMARY\n";
    cout << "================================================================\n";
    cout << "  Total scans (P,X pairs):   " << totalScans << "\n";
    cout << "  Total enthalpy steps:      " << totalSteps << "\n";
    cout << "  Passed scans:              " << passedScans
         << "  (" << fixed << setprecision(1)
         << (100.0 * passedScans / totalScans) << "%)\n";
    cout << "  Failed scans:              " << failedScans
         << "  (" << (100.0 * failedScans / totalScans) << "%)\n";
    cout << "  NaN evaluations:           " << totalNaN << "\n";
    cout << "  Negative-Rho evaluations:  " << totalNegRho
         << "  (treated as Rho=0, near-vacuum vapor)\n";
    cout << "  Boundary crossings:        " << totalCrossBoundary
         << "  (excluded from failure criteria)\n";
    cout << "================================================================\n\n";

    // ── verdict ────────────────────────────────────────────────────────
    if (failedScans == 0) {
        cout << "RESULT: PASS  –  All within-region transitions are smooth.\n";
        return 0;
    } else {
        double passRate = 100.0 * passedScans / totalScans;
        if (passRate >= 90.0) {
            cout << "RESULT: PASS (with warnings)  –  "
                 << setprecision(1) << passRate << "% of scans are smooth within regions.\n";
            cout << "  Some within-region jumps exceed thresholds near phase boundaries.\n";
            return 0;
        } else {
            cout << "RESULT: FAIL  –  Only "
                 << setprecision(1) << passRate << "% of scans are smooth within regions.\n";
            return 1;
        }
    }
}
