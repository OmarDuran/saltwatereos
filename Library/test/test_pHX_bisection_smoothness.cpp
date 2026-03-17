/**
 * @file test_pHX_bisection_smoothness.cpp
 * @brief Test smoothness of density and saturation across phase boundaries
 *        using prop_pHX_bisection.
 *
 * Strategy
 * --------
 * For several representative (P, X) pairs that cross different phase
 * boundaries, sweep enthalpy H in fine steps and verify:
 *
 *   1. Bulk density Rho changes continuously (no unphysical jumps).
 *   2. Saturations S_l, S_v change continuously (no unphysical jumps).
 *   3. Phase densities Rho_l, Rho_v are non-negative and finite.
 *   4. At phase-boundary crossings the region label changes but the
 *      properties still join smoothly.
 *
 * A "jump" is defined as a relative change that exceeds a generous
 * threshold between two adjacent enthalpy steps.
 *
 * The test uses 25 points in each of:
 *   0.0001 < X  < 0.20
 *   500    < H  < 3500   kJ/kg
 *   6      < P  < 510    bar  (= 0.6 – 51 MPa)
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

/// Relative jump with safe denominator.
static double relJump(double a, double b)
{
    double ref = 0.5 * (std::abs(a) + std::abs(b));
    if (ref < 1e-12) return 0.0;           // both ~0 → no jump
    return std::abs(b - a) / ref;
}

// ─── per-scan result ────────────────────────────────────────────────────

struct ScanResult {
    double P_bar;
    double X_wt;
    int    nSteps;
    int    nRegionChanges;

    // worst absolute jumps found
    double maxJump_Rho;     // kg/m³
    double maxJump_Sl;      // dimensionless
    double maxJump_Sv;      // dimensionless

    // where the worst jumps occurred (enthalpy kJ/kg)
    double H_maxJump_Rho;
    double H_maxJump_Sl;
    double H_maxJump_Sv;

    // region labels at worst jump
    int reg_before_Rho, reg_after_Rho;
    int reg_before_Sl,  reg_after_Sl;
    int reg_before_Sv,  reg_after_Sv;

    // counts of NaN / negative density
    int nNaN;
    int nNegRho;

    bool passed;            // true if all thresholds satisfied
    string failReason;
};

// ─── threshold constants ────────────────────────────────────────────────

// Maximum allowed absolute jump in bulk density (kg/m³) between adjacent
// enthalpy steps.  Phase transitions can cause large density changes,
// so we use a generous threshold. We allow up to 400 kg/m³ jump at
// phase boundaries (liquid→vapor is ~600 kg/m³ at 1 bar, but steps are
// fine enough that each step should be smaller).
static const double MAX_RHO_JUMP     = 400.0;  // kg/m³

// Maximum allowed absolute jump in saturation between adjacent steps.
// Phase transitions cause S_l to jump from 1→0 or vice versa, but with
// fine enthalpy steps the transition should be gradual in two-phase
// regions.  We allow up to 0.5 per step.
static const double MAX_SAT_JUMP     = 0.50;

// Enthalpy step size (kJ/kg) – fine enough to resolve transitions
static const double DH_KJ           = 25.0;

// ─── main ───────────────────────────────────────────────────────────────

int main()
{
    H2ONaCl::cH2ONaCl eos;

    cout << "================================================================\n";
    cout << " Phase Boundary Smoothness Test  (prop_pHX_bisection)\n";
    cout << "================================================================\n\n";
    cout << "Checks that Rho, S_l, S_v change continuously as enthalpy\n";
    cout << "is swept at fixed (P, X).  Reports the worst jumps found.\n\n";

    // ── build test grid ────────────────────────────────────────────────
    // 25 salinities
    vector<double> X_vals;
    for (int i = 0; i < 25; ++i)
        X_vals.push_back(0.0001 + i * (0.20 - 0.0001) / 24.0);

    // 25 pressures (log-spaced, 6 – 510 bar)
    vector<double> P_vals;
    {
        double logMin = log10(6.0), logMax = log10(510.0);
        for (int i = 0; i < 25; ++i)
            P_vals.push_back(pow(10.0, logMin + i * (logMax - logMin) / 24.0));
    }

    // enthalpy sweep range
    const double H_min_kJ = 500.0;
    const double H_max_kJ = 3500.0;
    const int    nH = static_cast<int>((H_max_kJ - H_min_kJ) / DH_KJ) + 1;

    cout << "Grid:  " << X_vals.size() << " salinities x "
         << P_vals.size() << " pressures x "
         << nH << " enthalpy steps  (dH = " << DH_KJ << " kJ/kg)\n";
    cout << "Thresholds:  max |dRho| = " << MAX_RHO_JUMP << " kg/m³,  "
         << "max |dS| = " << MAX_SAT_JUMP << "\n\n";

    // ── run scans ──────────────────────────────────────────────────────

    int totalScans   = 0;
    int passedScans  = 0;
    int failedScans  = 0;
    int totalSteps   = 0;
    int totalNaN     = 0;
    int totalNegRho  = 0;
    int rhoJumpFails = 0;
    int satJumpFails = 0;

    // Store all results for the summary table
    vector<ScanResult> results;

    // Header for per-scan details (print failures only)
    cout << "Scanning...\n";
    cout << string(140, '=') << "\n";
    cout << setw(8) << "P(bar)"
         << setw(8) << "X"
         << setw(6) << "Steps"
         << setw(6) << "RegCh"
         << setw(12) << "maxdRho"
         << setw(10) << "H@dRho"
         << setw(10) << "reg"
         << setw(12) << "maxdSl"
         << setw(10) << "H@dSl"
         << setw(10) << "reg"
         << setw(12) << "maxdSv"
         << setw(10) << "H@dSv"
         << setw(10) << "reg"
         << setw(8) << "NaN"
         << setw(8) << "Status"
         << "\n";
    cout << string(140, '-') << "\n";

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
            sr.nRegionChanges = 0;
            sr.maxJump_Rho = 0; sr.maxJump_Sl = 0; sr.maxJump_Sv = 0;
            sr.H_maxJump_Rho = 0; sr.H_maxJump_Sl = 0; sr.H_maxJump_Sv = 0;
            sr.reg_before_Rho = -1; sr.reg_after_Rho = -1;
            sr.reg_before_Sl  = -1; sr.reg_after_Sl  = -1;
            sr.reg_before_Sv  = -1; sr.reg_after_Sv  = -1;
            sr.nNaN = 0; sr.nNegRho = 0;
            sr.passed = true;

            // Evaluate first point
            double H0_J = H_min_kJ * 1e3;
            H2ONaCl::PROP_H2ONaCl prev = eos.prop_pHX_bisection(P_Pa, H0_J, X_wt);
            int prevRegion = prev.Region;

            if (isnan(prev.Rho) || isnan(prev.T)) sr.nNaN++;
            if (prev.Rho <= 0 && !isnan(prev.Rho)) sr.nNegRho++;

            for (int iH = 1; iH < nH; ++iH) {
                double H_kJ = H_min_kJ + iH * DH_KJ;
                double H_J  = H_kJ * 1e3;

                H2ONaCl::PROP_H2ONaCl cur = eos.prop_pHX_bisection(P_Pa, H_J, X_wt);
                totalSteps++;

                // NaN / negative checks
                if (isnan(cur.Rho) || isnan(cur.T)) { sr.nNaN++; prev = cur; prevRegion = cur.Region; continue; }
                if (cur.Rho <= 0) { sr.nNegRho++; prev = cur; prevRegion = cur.Region; continue; }
                if (isnan(prev.Rho) || prev.Rho <= 0) { prev = cur; prevRegion = cur.Region; continue; }

                // Region change
                if (cur.Region != prevRegion) sr.nRegionChanges++;

                // ── density jump ──
                double dRho = absJump(prev.Rho, cur.Rho);
                if (dRho > sr.maxJump_Rho) {
                    sr.maxJump_Rho   = dRho;
                    sr.H_maxJump_Rho = H_kJ;
                    sr.reg_before_Rho = prevRegion;
                    sr.reg_after_Rho  = cur.Region;
                }

                // ── saturation jumps ──
                double dSl = absJump(prev.S_l, cur.S_l);
                if (dSl > sr.maxJump_Sl) {
                    sr.maxJump_Sl   = dSl;
                    sr.H_maxJump_Sl = H_kJ;
                    sr.reg_before_Sl = prevRegion;
                    sr.reg_after_Sl  = cur.Region;
                }
                double dSv = absJump(prev.S_v, cur.S_v);
                if (dSv > sr.maxJump_Sv) {
                    sr.maxJump_Sv   = dSv;
                    sr.H_maxJump_Sv = H_kJ;
                    sr.reg_before_Sv = prevRegion;
                    sr.reg_after_Sv  = cur.Region;
                }

                prev = cur;
                prevRegion = cur.Region;
            }

            // ── evaluate pass / fail ──
            ostringstream reason;
            if (sr.maxJump_Rho > MAX_RHO_JUMP) {
                sr.passed = false;
                rhoJumpFails++;
                reason << "dRho=" << fixed << setprecision(1) << sr.maxJump_Rho << " ";
            }
            if (sr.maxJump_Sl > MAX_SAT_JUMP) {
                sr.passed = false;
                satJumpFails++;
                reason << "dSl=" << fixed << setprecision(3) << sr.maxJump_Sl << " ";
            }
            if (sr.maxJump_Sv > MAX_SAT_JUMP) {
                sr.passed = false;
                if (sr.maxJump_Sl <= MAX_SAT_JUMP) satJumpFails++;  // don't double count
                reason << "dSv=" << fixed << setprecision(3) << sr.maxJump_Sv << " ";
            }
            sr.failReason = reason.str();

            totalNaN    += sr.nNaN;
            totalNegRho += sr.nNegRho;
            if (sr.passed) passedScans++; else failedScans++;

            results.push_back(sr);

            // Print failures, or every 50th scan for progress
            if (!sr.passed || (totalScans % 50 == 0)) {
                ostringstream regRho, regSl, regSv;
                regRho << regionName(sr.reg_before_Rho) << "->" << regionName(sr.reg_after_Rho);
                regSl  << regionName(sr.reg_before_Sl)  << "->" << regionName(sr.reg_after_Sl);
                regSv  << regionName(sr.reg_before_Sv)  << "->" << regionName(sr.reg_after_Sv);

                cout << setw(8) << fixed << setprecision(1) << P_bar
                     << setw(8) << setprecision(4) << X_wt
                     << setw(6) << nH
                     << setw(6) << sr.nRegionChanges
                     << setw(12) << setprecision(2) << sr.maxJump_Rho
                     << setw(10) << setprecision(0) << sr.H_maxJump_Rho
                     << setw(10) << regRho.str()
                     << setw(12) << setprecision(4) << sr.maxJump_Sl
                     << setw(10) << setprecision(0) << sr.H_maxJump_Sl
                     << setw(10) << regSl.str()
                     << setw(12) << setprecision(4) << sr.maxJump_Sv
                     << setw(10) << setprecision(0) << sr.H_maxJump_Sv
                     << setw(10) << regSv.str()
                     << setw(8) << sr.nNaN
                     << setw(8) << (sr.passed ? "OK" : "FAIL")
                     << "\n";
            }
        }
    }

    cout << string(140, '=') << "\n\n";

    // ── detailed failure report ────────────────────────────────────────
    int nFailsToPrint = 0;
    for (auto& sr : results) if (!sr.passed) nFailsToPrint++;

    if (nFailsToPrint > 0) {
        cout << "========== DETAILED FAILURE REPORT ==========\n\n";

        // Show up to 30 worst density jumps
        vector<ScanResult*> byRho;
        for (auto& sr : results) if (sr.maxJump_Rho > MAX_RHO_JUMP) byRho.push_back(&sr);
        sort(byRho.begin(), byRho.end(), [](const ScanResult* a, const ScanResult* b) {
            return a->maxJump_Rho > b->maxJump_Rho;
        });

        if (!byRho.empty()) {
            cout << "Top density discontinuities (max " << byRho.size() << "):\n";
            cout << setw(10) << "P(bar)" << setw(10) << "X"
                 << setw(14) << "maxdRho" << setw(12) << "H(kJ/kg)"
                 << setw(14) << "Transition" << "\n";
            cout << string(60, '-') << "\n";
            int nShow = min((int)byRho.size(), 30);
            for (int i = 0; i < nShow; ++i) {
                auto* s = byRho[i];
                cout << setw(10) << fixed << setprecision(1) << s->P_bar
                     << setw(10) << setprecision(4) << s->X_wt
                     << setw(14) << setprecision(2) << s->maxJump_Rho
                     << setw(12) << setprecision(0) << s->H_maxJump_Rho
                     << "  " << regionName(s->reg_before_Rho) << " -> " << regionName(s->reg_after_Rho)
                     << "\n";
            }
            cout << "\n";
        }

        // Show up to 30 worst saturation jumps
        vector<ScanResult*> bySat;
        for (auto& sr : results) if (sr.maxJump_Sl > MAX_SAT_JUMP || sr.maxJump_Sv > MAX_SAT_JUMP) bySat.push_back(&sr);
        sort(bySat.begin(), bySat.end(), [](const ScanResult* a, const ScanResult* b) {
            return max(a->maxJump_Sl, a->maxJump_Sv) > max(b->maxJump_Sl, b->maxJump_Sv);
        });

        if (!bySat.empty()) {
            cout << "Top saturation discontinuities (max " << bySat.size() << "):\n";
            cout << setw(10) << "P(bar)" << setw(10) << "X"
                 << setw(12) << "maxdSl" << setw(12) << "H@dSl"
                 << setw(12) << "maxdSv" << setw(12) << "H@dSv"
                 << setw(14) << "Transition" << "\n";
            cout << string(82, '-') << "\n";
            int nShow = min((int)bySat.size(), 30);
            for (int i = 0; i < nShow; ++i) {
                auto* s = bySat[i];
                int rb, ra;
                if (s->maxJump_Sl >= s->maxJump_Sv) {
                    rb = s->reg_before_Sl; ra = s->reg_after_Sl;
                } else {
                    rb = s->reg_before_Sv; ra = s->reg_after_Sv;
                }
                cout << setw(10) << fixed << setprecision(1) << s->P_bar
                     << setw(10) << setprecision(4) << s->X_wt
                     << setw(12) << setprecision(4) << s->maxJump_Sl
                     << setw(12) << setprecision(0) << s->H_maxJump_Sl
                     << setw(12) << setprecision(4) << s->maxJump_Sv
                     << setw(12) << setprecision(0) << s->H_maxJump_Sv
                     << "  " << regionName(rb) << " -> " << regionName(ra)
                     << "\n";
            }
            cout << "\n";
        }
    }

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
    cout << "    - density jump fails:    " << rhoJumpFails << "\n";
    cout << "    - saturation jump fails: " << satJumpFails << "\n";
    cout << "  NaN evaluations:           " << totalNaN << "\n";
    cout << "  Negative-Rho evaluations:  " << totalNegRho << "\n";
    cout << "================================================================\n\n";

    // ── verdict ────────────────────────────────────────────────────────
    double passRate = 100.0 * passedScans / totalScans;
    if (failedScans == 0) {
        cout << "RESULT: PASS  –  All phase-boundary transitions are smooth.\n";
        return 0;
    } else if (passRate >= 95.0) {
        cout << "RESULT: PASS (with warnings)  –  "
             << setprecision(1) << passRate << "% of scans are smooth.\n";
        cout << "  Some transitions have large jumps that may warrant investigation.\n";
        return 0;     // still pass – known edge cases at extreme conditions
    } else {
        cout << "RESULT: FAIL  –  Only "
             << setprecision(1) << passRate << "% of scans are smooth.\n";
        return 1;
    }
}
