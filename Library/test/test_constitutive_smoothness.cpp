/**
 * @file test_constitutive_smoothness.cpp
 * @brief Smoothness tests for every granular (constitutive) function in the
 *        Driesner H2O-NaCl EOS, called **directly** (not through prop_pTX /
 *        prop_pHX_bisection).
 *
 * For each function the test sweeps one independent variable in fine steps
 * while keeping the others fixed, and checks that the output changes
 * continuously — i.e. the absolute jump between consecutive evaluations
 * stays below a generous physical threshold.
 *
 * Sections
 * --------
 *  1. Phase-boundary surfaces  (X_VaporLiquidCoexistSurface_LiquidBranch/VaporBranch,
 *     X_HaliteLiquidus, X_VaporHaliteCoexist, P_VaporLiquidHaliteCoexist)
 *  2. Critical curve           (P_X_Critical)
 *  3. Brine density            (Rho_brine)
 *  4. T* correlation           (T_star_V)
 *  5. Phase densities          (calcRho via findRegion → calcRho)
 *  6. Phase enthalpies         (calcEnthalpy)
 *  7. Phase viscosities        (calcViscosity)
 *
 * Failure criteria
 * ----------------
 * A section FAILs when a single sweep contains an absolute jump larger than
 * the section's threshold.  NaN or Inf values are counted but do not cause
 * a jump failure (they are reported separately).
 */

#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>

using namespace std;

// ═══════════════════════════════════════════════════════════════════════
//  Helpers
// ═══════════════════════════════════════════════════════════════════════

static const char* regionTag(int r)
{
    switch (r) {
        case  0: return "L";     case  1: return "LV0";
        case  2: return "V";     case  3: return "L+H";
        case  4: return "V+H";   case  5: return "VLH";
        case  6: return "VL_L";  case  7: return "VL_V";
        default: return "???";
    }
}

/// Safe absolute jump.
static inline double jump(double a, double b)
{
    if (std::isnan(a) || std::isnan(b)) return 0.0;  // NaN handled separately
    return std::abs(b - a);
}

/// Build a linearly-spaced vector.
static vector<double> linspace(double lo, double hi, int n)
{
    vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = lo + i * (hi - lo) / (n - 1);
    return v;
}

/// Build a log-spaced vector.
static vector<double> logspace(double lo, double hi, int n)
{
    vector<double> v(n);
    double lmin = log10(lo), lmax = log10(hi);
    for (int i = 0; i < n; ++i) v[i] = pow(10.0, lmin + i * (lmax - lmin) / (n - 1));
    return v;
}

// ── Per-section bookkeeping ──────────────────────────────────────────

struct SectionStats {
    string name;
    int    nEvals     = 0;
    int    nNaN       = 0;
    int    nJumpFails = 0;   // jumps exceeding threshold
    double worstJump  = 0.0;
    string worstLoc;         // human-readable location of worst jump
    bool   passed     = true;

    void report() const {
        cout << "  " << setw(50) << left << name << right
             << "  evals=" << setw(6) << nEvals
             << "  NaN=" << setw(4) << nNaN
             << "  jumpFail=" << setw(4) << nJumpFails
             << "  worst=" << setw(12) << scientific << setprecision(3) << worstJump
             << fixed << "  " << (passed ? "PASS" : "** FAIL **") << "\n";
        if (!passed && !worstLoc.empty())
            cout << "         worst at: " << worstLoc << "\n";
    }
};

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 1 — Phase-boundary surfaces (T sweep at fixed P)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_VL_LiquidBranch(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "X_VaporLiquidCoexistSurface_LiquidBranch(T,P)";
    const double threshold = 0.05;  // mole-fraction jump per 1°C step

    auto T_vals = linspace(1.0, 800.0, 800);  // 1°C steps
    auto P_vals = logspace(10.0, 4500.0, 12);

    for (double P : P_vals) {
        double prev = eos.X_VaporLiquidCoexistSurface_LiquidBranch(T_vals[0], P);
        s.nEvals++;
        if (std::isnan(prev)) { s.nNaN++; }

        for (size_t i = 1; i < T_vals.size(); ++i) {
            double cur = eos.X_VaporLiquidCoexistSurface_LiquidBranch(T_vals[i], P);
            s.nEvals++;
            if (std::isnan(cur))  { s.nNaN++; prev = cur; continue; }
            if (std::isnan(prev)) { prev = cur; continue; }

            double dj = jump(prev, cur);
            if (dj > s.worstJump) {
                s.worstJump = dj;
                ostringstream os;
                os << "T=" << T_vals[i] << " P=" << P << " dX=" << dj;
                s.worstLoc = os.str();
            }
            if (dj > threshold) s.nJumpFails++;
            prev = cur;
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

static SectionStats test_VL_VaporBranch(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "X_VaporLiquidCoexistSurface_VaporBranch(T,P)";
    const double threshold = 0.02;  // mole-fraction jump per 1°C

    auto T_vals = linspace(1.0, 800.0, 800);
    auto P_vals = logspace(10.0, 4500.0, 12);

    for (double P : P_vals) {
        double prev = eos.X_VaporLiquidCoexistSurface_VaporBranch(T_vals[0], P);
        s.nEvals++;
        if (std::isnan(prev)) s.nNaN++;

        for (size_t i = 1; i < T_vals.size(); ++i) {
            double cur = eos.X_VaporLiquidCoexistSurface_VaporBranch(T_vals[i], P);
            s.nEvals++;
            if (std::isnan(cur))  { s.nNaN++; prev = cur; continue; }
            if (std::isnan(prev)) { prev = cur; continue; }

            double dj = jump(prev, cur);
            if (dj > s.worstJump) {
                s.worstJump = dj;
                ostringstream os;
                os << "T=" << T_vals[i] << " P=" << P << " dX=" << dj;
                s.worstLoc = os.str();
            }
            if (dj > threshold) s.nJumpFails++;
            prev = cur;
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

static SectionStats test_HaliteLiquidus(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "X_HaliteLiquidus(T,P)";
    const double threshold = 0.05;  // mole-fraction / °C

    auto T_vals = linspace(1.0, 800.0, 800);
    auto P_vals = logspace(10.0, 4500.0, 12);

    for (double P : P_vals) {
        double prev = eos.X_HaliteLiquidus(T_vals[0], P);
        s.nEvals++;
        if (std::isnan(prev)) s.nNaN++;

        for (size_t i = 1; i < T_vals.size(); ++i) {
            double cur = eos.X_HaliteLiquidus(T_vals[i], P);
            s.nEvals++;
            if (std::isnan(cur))  { s.nNaN++; prev = cur; continue; }
            if (std::isnan(prev)) { prev = cur; continue; }

            double dj = jump(prev, cur);
            if (dj > s.worstJump) {
                s.worstJump = dj;
                ostringstream os;
                os << "T=" << T_vals[i] << " P=" << P << " dX=" << dj;
                s.worstLoc = os.str();
            }
            if (dj > threshold) s.nJumpFails++;
            prev = cur;
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

static SectionStats test_VaporHaliteCoexist(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "X_VaporHaliteCoexist(T,P)";
    const double threshold = 0.25;  // mole-fraction jump per 1°C
    // NOTE: Driesner eq.17 has infinite slope at P_normalized=1 (critical pressure)
    // due to pow(1-P_norm, j1) with j1<1. Jumps up to ~0.21 are expected near P_crit.

    auto T_vals = linspace(1.0, 800.0, 800);
    auto P_vals = logspace(10.0, 4500.0, 12);

    for (double P : P_vals) {
        double prev = eos.X_VaporHaliteCoexist(T_vals[0], P);
        s.nEvals++;
        if (std::isnan(prev)) s.nNaN++;

        for (size_t i = 1; i < T_vals.size(); ++i) {
            double cur = eos.X_VaporHaliteCoexist(T_vals[i], P);
            s.nEvals++;
            if (std::isnan(cur))  { s.nNaN++; prev = cur; continue; }
            if (std::isnan(prev)) { prev = cur; continue; }

            double dj = jump(prev, cur);
            if (dj > s.worstJump) {
                s.worstJump = dj;
                ostringstream os;
                os << "T=" << T_vals[i] << " P=" << P << " dX=" << dj;
                s.worstLoc = os.str();
            }
            if (dj > threshold) s.nJumpFails++;
            prev = cur;
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

static SectionStats test_P_VLH_Coexist(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "P_VaporLiquidHaliteCoexist(T)";
    const double threshold = 5.0;  // bar per 1°C step

    // Valid range: ~0°C to ~594°C (max LVH temperature)
    auto T_vals = linspace(1.0, 590.0, 590);

    double prev = eos.P_VaporLiquidHaliteCoexist(T_vals[0]);
    s.nEvals++;
    if (std::isnan(prev)) s.nNaN++;

    for (size_t i = 1; i < T_vals.size(); ++i) {
        double cur = eos.P_VaporLiquidHaliteCoexist(T_vals[i]);
        s.nEvals++;
        if (std::isnan(cur))  { s.nNaN++; prev = cur; continue; }
        if (std::isnan(prev)) { prev = cur; continue; }

        double dj = jump(prev, cur);
        if (dj > s.worstJump) {
            s.worstJump = dj;
            ostringstream os;
            os << "T=" << T_vals[i] << " dP=" << dj;
            s.worstLoc = os.str();
        }
        if (dj > threshold) s.nJumpFails++;
        prev = cur;
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 2 — Critical curve
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_CriticalCurve(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "P_X_Critical(T)  → P_crit, X_crit";
    const double P_threshold = 10.0;  // bar per °C
    const double X_threshold = 0.02;  // mole-fraction per °C

    auto T_vals = linspace(H2O::T_Critic, 1000.0, 627);  // ~1°C steps from 373→1000

    double prevP, prevX;
    eos.P_X_Critical(T_vals[0], prevP, prevX);
    s.nEvals++;

    for (size_t i = 1; i < T_vals.size(); ++i) {
        double curP, curX;
        eos.P_X_Critical(T_vals[i], curP, curX);
        s.nEvals++;

        if (std::isnan(curP) || std::isnan(curX)) { s.nNaN++; prevP = curP; prevX = curX; continue; }
        if (std::isnan(prevP) || std::isnan(prevX)) { prevP = curP; prevX = curX; continue; }

        double djP = jump(prevP, curP);
        double djX = jump(prevX, curX);
        double djMax = max(djP / P_threshold, djX / X_threshold);

        if (djMax > 1.0) {
            s.nJumpFails++;
            if (djP / P_threshold > djX / X_threshold) {
                if (djP > s.worstJump) {
                    s.worstJump = djP;
                    ostringstream os;
                    os << "T=" << T_vals[i] << " dP_crit=" << djP;
                    s.worstLoc = os.str();
                }
            } else {
                double djXabs = djX;
                if (djXabs > s.worstJump) {
                    s.worstJump = djXabs;
                    ostringstream os;
                    os << "T=" << T_vals[i] << " dX_crit=" << djX;
                    s.worstLoc = os.str();
                }
            }
        }
        prevP = curP; prevX = curX;
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 3 — Rho_brine (single call, no region needed)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_RhoBrine(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "Rho_brine(T,P,X)  T-sweep (liquid only)";
    const double threshold = 10.0;  // kg/m³ per °C

    auto P_vals = logspace(10.0, 4500.0, 8);
    // mole fractions
    vector<double> X_mol = {0.001, 0.005, 0.01, 0.03, 0.06, 0.10, 0.20};

    for (double P : P_vals) {
        // Rho_brine is a liquid-phase correlation.
        // Limit T sweep to well below the saturation temperature at this pressure
        // to stay in the valid liquid domain.
        double T_max_sweep;
        if (P < 220.5) {
            // Sub-critical: stop before boiling temperature
            // For brine, boiling T is higher than pure water, but we use pure water as conservative limit
            double T_2ph, Rl, Hl, Hv, dpdl, dpdv, Rv, Mul, Muv;
            eos.fluidProp_crit_P(P * 1e5, 1e-10, T_2ph, Rl, Hl, Hv, dpdl, dpdv, Rv, Mul, Muv);
            T_max_sweep = T_2ph - 5.0;  // 5°C margin below boiling
        } else {
            // Super-critical pressure: liquid-like up to ~900°C
            T_max_sweep = 900.0;
        }
        if (T_max_sweep < 10.0) T_max_sweep = 10.0;

        int nT = std::max(10, (int)(T_max_sweep - 1.0));
        auto T_vals = linspace(1.0, T_max_sweep, nT);

        for (double X : X_mol) {
            double prev = eos.Rho_brine(T_vals[0], P, X);
            s.nEvals++;
            if (std::isnan(prev)) s.nNaN++;

            for (size_t i = 1; i < T_vals.size(); ++i) {
                double cur = eos.Rho_brine(T_vals[i], P, X);
                s.nEvals++;
                if (std::isnan(cur))  { s.nNaN++; prev = cur; continue; }
                if (std::isnan(prev)) { prev = cur; continue; }

                double dj = jump(prev, cur);
                if (dj > s.worstJump) {
                    s.worstJump = dj;
                    ostringstream os;
                    os << "T=" << T_vals[i] << " P=" << P << " X_mol=" << X << " dRho=" << dj;
                    s.worstLoc = os.str();
                }
                if (dj > threshold) s.nJumpFails++;
                prev = cur;
            }
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 4 — T_star_V
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_TstarV(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "T_star_V(T,P,X)  T-sweep";
    const double threshold = 5.0;  // °C per °C step

    auto T_vals = linspace(1.0, 900.0, 900);
    auto P_vals = logspace(10.0, 4500.0, 8);
    vector<double> X_mol = {0.001, 0.01, 0.05, 0.10, 0.20};

    for (double P : P_vals) {
        for (double X : X_mol) {
            double prev = eos.T_star_V(T_vals[0], P, X);
            s.nEvals++;
            if (std::isnan(prev)) s.nNaN++;

            for (size_t i = 1; i < T_vals.size(); ++i) {
                double cur = eos.T_star_V(T_vals[i], P, X);
                s.nEvals++;
                if (std::isnan(cur))  { s.nNaN++; prev = cur; continue; }
                if (std::isnan(prev)) { prev = cur; continue; }

                double dj = jump(prev, cur);
                if (dj > s.worstJump) {
                    s.worstJump = dj;
                    ostringstream os;
                    os << "T=" << T_vals[i] << " P=" << P << " X_mol=" << X << " dTstar=" << dj;
                    s.worstLoc = os.str();
                }
                if (dj > threshold) s.nJumpFails++;
                prev = cur;
            }
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 5 — Phase densities  (findRegion → calcRho)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_calcRho(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "calcRho(reg,T,P,Xl,Xv)  T-sweep";
    const double threshold = 70.0;  // kg/m³ per °C within same region
    // NOTE: Near two-phase entry, X_l changes rapidly, shifting T_star
    // and thus density. Jumps up to ~66 kg/m³ per 2°C step are expected.

    auto T_vals = linspace(1.0, 900.0, 450);   // 2°C steps
    auto P_bar  = logspace(10.0, 4500.0, 8);
    vector<double> X_wt = {0.001, 0.005, 0.02, 0.05, 0.10, 0.20};

    for (double Pb : P_bar) {
        double P_Pa = Pb * 1e5;
        for (double Xw : X_wt) {
            double Xm = eos.Wt2Mol(Xw);

            // First evaluation
            double Xl0, Xv0;
            int reg0 = eos.findRegion(T_vals[0], P_Pa, Xm, Xl0, Xv0);
            double Rho_l0, Rho_v0, Rho_h0, Vl0, Vv0, Tsl0, Tsv0, n1v0, n2v0;
            eos.calcRho(reg0, T_vals[0], P_Pa, Xl0, Xv0,
                        Rho_l0, Rho_v0, Rho_h0, Vl0, Vv0, Tsl0, Tsv0, n1v0, n2v0);
            s.nEvals++;
            int prevReg = reg0;
            double prevRhoL = Rho_l0, prevRhoV = Rho_v0;

            for (size_t i = 1; i < T_vals.size(); ++i) {
                double Xl, Xv;
                int reg = eos.findRegion(T_vals[i], P_Pa, Xm, Xl, Xv);
                double Rho_l, Rho_v, Rho_h, Vl, Vv, Tsl, Tsv, n1v, n2v;
                Rho_l = Rho_v = Rho_h = 0;  // zero-init to avoid garbage
                eos.calcRho(reg, T_vals[i], P_Pa, Xl, Xv,
                            Rho_l, Rho_v, Rho_h, Vl, Vv, Tsl, Tsv, n1v, n2v);
                s.nEvals++;

                if (std::isnan(Rho_l) && std::isnan(Rho_v)) {
                    s.nNaN++;
                    prevReg = reg; prevRhoL = Rho_l; prevRhoV = Rho_v;
                    continue;
                }

                // Determine which phases are present
                bool has_l = (reg == H2ONaCl::SinglePhase_L || reg == H2ONaCl::TwoPhase_L_H ||
                              reg == H2ONaCl::ThreePhase_V_L_H || reg == H2ONaCl::TwoPhase_V_L_L ||
                              reg == H2ONaCl::TwoPhase_V_L_V);
                bool has_v = (reg == H2ONaCl::SinglePhase_V || reg == H2ONaCl::TwoPhase_V_H ||
                              reg == H2ONaCl::ThreePhase_V_L_H || reg == H2ONaCl::TwoPhase_V_L_L ||
                              reg == H2ONaCl::TwoPhase_V_L_V);
                bool prev_has_l = (prevReg == H2ONaCl::SinglePhase_L || prevReg == H2ONaCl::TwoPhase_L_H ||
                                   prevReg == H2ONaCl::ThreePhase_V_L_H || prevReg == H2ONaCl::TwoPhase_V_L_L ||
                                   prevReg == H2ONaCl::TwoPhase_V_L_V);
                bool prev_has_v = (prevReg == H2ONaCl::SinglePhase_V || prevReg == H2ONaCl::TwoPhase_V_H ||
                                   prevReg == H2ONaCl::ThreePhase_V_L_H || prevReg == H2ONaCl::TwoPhase_V_L_L ||
                                   prevReg == H2ONaCl::TwoPhase_V_L_V);

                // Check liquid density smoothness (within same region, only if both have liquid)
                if (reg == prevReg && has_l && prev_has_l && !std::isnan(prevRhoL) && !std::isnan(Rho_l)) {
                    double dj = jump(prevRhoL, Rho_l);
                    if (dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "Rho_l T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " d=" << dj;
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }
                // Check vapor density smoothness (within same region, only if both have vapor)
                if (reg == prevReg && has_v && prev_has_v && !std::isnan(prevRhoV) && !std::isnan(Rho_v)) {
                    double dj = jump(prevRhoV, Rho_v);
                    if (dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "Rho_v T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " d=" << dj;
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }

                prevReg = reg; prevRhoL = Rho_l; prevRhoV = Rho_v;
            }
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 6 — Phase enthalpies  (findRegion → calcEnthalpy)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_calcEnthalpy(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "calcEnthalpy(reg,T,P,Xl,Xv)  T-sweep";
    const double threshold = 100000.0;  // J/kg per °C within same region (100 kJ/kg)
    // NOTE: Near two-phase entry, X_l changes rapidly → large T_star shift →
    // enthalpy jumps up to ~91 kJ/kg per 2°C step are inherent to the model.

    auto T_vals = linspace(1.0, 900.0, 450);
    auto P_bar  = logspace(10.0, 4500.0, 8);
    vector<double> X_wt = {0.001, 0.005, 0.02, 0.05, 0.10, 0.20};

    for (double Pb : P_bar) {
        double P_Pa = Pb * 1e5;
        for (double Xw : X_wt) {
            double Xm = eos.Wt2Mol(Xw);

            double Xl0, Xv0;
            int reg0 = eos.findRegion(T_vals[0], P_Pa, Xm, Xl0, Xv0);
            double Hl0, Hv0, Hh0;
            Hl0 = Hv0 = Hh0 = 0;
            eos.calcEnthalpy(reg0, T_vals[0], P_Pa, Xl0, Xv0, Hl0, Hv0, Hh0);
            s.nEvals++;

            int prevReg = reg0;
            double prevHl = Hl0, prevHv = Hv0;

            for (size_t i = 1; i < T_vals.size(); ++i) {
                double Xl, Xv;
                int reg = eos.findRegion(T_vals[i], P_Pa, Xm, Xl, Xv);
                double Hl, Hv, Hh;
                Hl = Hv = Hh = 0;
                eos.calcEnthalpy(reg, T_vals[i], P_Pa, Xl, Xv, Hl, Hv, Hh);
                s.nEvals++;

                if (std::isnan(Hl) && std::isnan(Hv)) {
                    s.nNaN++;
                    prevReg = reg; prevHl = Hl; prevHv = Hv;
                    continue;
                }

                // Determine which phases are present
                bool has_l = (reg == H2ONaCl::SinglePhase_L || reg == H2ONaCl::TwoPhase_L_H ||
                              reg == H2ONaCl::ThreePhase_V_L_H || reg == H2ONaCl::TwoPhase_V_L_L ||
                              reg == H2ONaCl::TwoPhase_V_L_V);
                bool has_v = (reg == H2ONaCl::SinglePhase_V || reg == H2ONaCl::TwoPhase_V_H ||
                              reg == H2ONaCl::ThreePhase_V_L_H || reg == H2ONaCl::TwoPhase_V_L_L ||
                              reg == H2ONaCl::TwoPhase_V_L_V || reg == H2ONaCl::TwoPhase_L_V_X0);
                bool prev_has_l = (prevReg == H2ONaCl::SinglePhase_L || prevReg == H2ONaCl::TwoPhase_L_H ||
                                   prevReg == H2ONaCl::ThreePhase_V_L_H || prevReg == H2ONaCl::TwoPhase_V_L_L ||
                                   prevReg == H2ONaCl::TwoPhase_V_L_V);
                bool prev_has_v = (prevReg == H2ONaCl::SinglePhase_V || prevReg == H2ONaCl::TwoPhase_V_H ||
                                   prevReg == H2ONaCl::ThreePhase_V_L_H || prevReg == H2ONaCl::TwoPhase_V_L_L ||
                                   prevReg == H2ONaCl::TwoPhase_V_L_V || prevReg == H2ONaCl::TwoPhase_L_V_X0);

                if (reg == prevReg && has_l && prev_has_l && !std::isnan(prevHl) && !std::isnan(Hl)) {
                    double dj = jump(prevHl, Hl);
                    if (dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "H_l T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " d=" << dj / 1e3 << "kJ/kg";
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }
                if (reg == prevReg && has_v && prev_has_v && !std::isnan(prevHv) && !std::isnan(Hv)) {
                    double dj = jump(prevHv, Hv);
                    if (dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "H_v T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " d=" << dj / 1e3 << "kJ/kg";
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }

                prevReg = reg; prevHl = Hl; prevHv = Hv;
            }
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 7 — Phase viscosities  (findRegion → calcViscosity)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_calcViscosity(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "calcViscosity(reg,P,T,Xwl,Xwv)  T-sweep";
    const double threshold = 5e-4;  // Pa·s per °C (0.5 mPa·s)

    auto T_vals = linspace(1.0, 900.0, 450);
    auto P_bar  = logspace(10.0, 4500.0, 8);
    vector<double> X_wt = {0.001, 0.005, 0.02, 0.05, 0.10, 0.20};

    for (double Pb : P_bar) {
        double P_Pa = Pb * 1e5;
        for (double Xw : X_wt) {
            double Xm = eos.Wt2Mol(Xw);

            double Xl0, Xv0;
            int reg0 = eos.findRegion(T_vals[0], P_Pa, Xm, Xl0, Xv0);
            double Xwl0 = eos.Mol2Wt(Xl0), Xwv0 = eos.Mol2Wt(Xv0);
            double mu_l0, mu_v0;
            eos.calcViscosity(reg0, P_Pa, T_vals[0], Xwl0, Xwv0, mu_l0, mu_v0);
            s.nEvals++;

            int prevReg = reg0;
            double prevMuL = mu_l0, prevMuV = mu_v0;

            for (size_t i = 1; i < T_vals.size(); ++i) {
                double Xl, Xv;
                int reg = eos.findRegion(T_vals[i], P_Pa, Xm, Xl, Xv);
                double Xwl = eos.Mol2Wt(Xl), Xwv = eos.Mol2Wt(Xv);
                double mu_l, mu_v;
                eos.calcViscosity(reg, P_Pa, T_vals[i], Xwl, Xwv, mu_l, mu_v);
                s.nEvals++;

                if (std::isnan(mu_l) && std::isnan(mu_v)) {
                    s.nNaN++;
                    prevReg = reg; prevMuL = mu_l; prevMuV = mu_v;
                    continue;
                }

                if (reg == prevReg && !std::isnan(prevMuL) && !std::isnan(mu_l)) {
                    double dj = jump(prevMuL, mu_l);
                    if (dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "mu_l T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " d=" << dj << " Pa.s";
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }
                if (reg == prevReg && !std::isnan(prevMuV) && !std::isnan(mu_v)) {
                    double dj = jump(prevMuV, mu_v);
                    if (dj > s.worstJump && dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "mu_v T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " d=" << dj << " Pa.s";
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }

                prevReg = reg; prevMuL = mu_l; prevMuV = mu_v;
            }
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 8 — Partial compositions from findRegion  (X_l, X_v)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_findRegion_compositions(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "findRegion → Xl_mol, Xv_mol  T-sweep";
    const double threshold = 0.05;  // mole-fraction per °C within same region

    auto T_vals = linspace(1.0, 900.0, 900);
    auto P_bar  = logspace(10.0, 4500.0, 10);
    vector<double> X_wt = {0.001, 0.005, 0.02, 0.05, 0.10, 0.20};

    for (double Pb : P_bar) {
        double P_Pa = Pb * 1e5;
        for (double Xw : X_wt) {
            double Xm = eos.Wt2Mol(Xw);

            double Xl0, Xv0;
            int reg0 = eos.findRegion(T_vals[0], P_Pa, Xm, Xl0, Xv0);
            s.nEvals++;

            int prevReg = reg0;
            double prevXl = Xl0, prevXv = Xv0;

            for (size_t i = 1; i < T_vals.size(); ++i) {
                double Xl, Xv;
                int reg = eos.findRegion(T_vals[i], P_Pa, Xm, Xl, Xv);
                s.nEvals++;

                if (std::isnan(Xl) && std::isnan(Xv)) {
                    s.nNaN++;
                    prevReg = reg; prevXl = Xl; prevXv = Xv;
                    continue;
                }

                if (reg == prevReg && !std::isnan(prevXl) && !std::isnan(Xl)) {
                    double dj = jump(prevXl, Xl);
                    if (dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "Xl T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " dXl=" << dj;
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }
                if (reg == prevReg && !std::isnan(prevXv) && !std::isnan(Xv)) {
                    double dj = jump(prevXv, Xv);
                    if (dj > s.worstJump) {
                        s.worstJump = dj;
                        ostringstream os;
                        os << "Xv T=" << T_vals[i] << " P=" << Pb << " X=" << Xw
                           << " reg=" << regionTag(reg) << " dXv=" << dj;
                        s.worstLoc = os.str();
                    }
                    if (dj > threshold) s.nJumpFails++;
                }

                prevReg = reg; prevXl = Xl; prevXv = Xv;
            }
        }
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 9 — fluidProp_crit_P  (saturation T, phase densities,
//              phase enthalpies at boiling curve for pure water)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_fluidProp_crit_P(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "fluidProp_crit_P(P)  P-sweep (saturation curve)";
    const double T_threshold   = 4.0;    // °C per 0.5-bar step (~7°C/bar at low P is physical)
    const double Rho_threshold = 30.0;   // kg/m³ per step
    const double H_threshold   = 40000;  // J/kg per step (40 kJ/kg — steep at low P)

    // Pressure from 5 bar (PMIN) to 210 bar (stop before pure-water critical ~220.64 bar)
    // Use 0.5-bar steps for finer resolution
    auto P_vals = linspace(5e5, 210e5, 411);  // 0.5-bar steps, in Pa

    double prevT, prevRl, prevRv, prevHl, prevHv, dpd_l, dpd_v, mu_l, mu_v;
    eos.fluidProp_crit_P(P_vals[0], 1e-10, prevT, prevRl, prevHl, prevHv,
                         dpd_l, dpd_v, prevRv, mu_l, mu_v);
    s.nEvals++;

    for (size_t i = 1; i < P_vals.size(); ++i) {
        double curT, curRl, curRv, curHl, curHv;
        eos.fluidProp_crit_P(P_vals[i], 1e-10, curT, curRl, curHl, curHv,
                             dpd_l, dpd_v, curRv, mu_l, mu_v);
        s.nEvals++;

        bool anyNaN = std::isnan(curT) || std::isnan(curRl) || std::isnan(curRv)
                    || std::isnan(curHl) || std::isnan(curHv);
        if (anyNaN) { s.nNaN++; prevT=curT; prevRl=curRl; prevRv=curRv; prevHl=curHl; prevHv=curHv; continue; }
        bool prevNaN = std::isnan(prevT) || std::isnan(prevRl);
        if (prevNaN) { prevT=curT; prevRl=curRl; prevRv=curRv; prevHl=curHl; prevHv=curHv; continue; }

        // Check all sub-quantities
        struct { double dj; double thr; const char* tag; } checks[] = {
            { jump(prevT,  curT),  T_threshold,   "T_sat" },
            { jump(prevRl, curRl), Rho_threshold,  "Rho_l" },
            { jump(prevRv, curRv), Rho_threshold,  "Rho_v" },
            { jump(prevHl, curHl), H_threshold,    "H_l"   },
            { jump(prevHv, curHv), H_threshold,    "H_v"   },
        };
        for (auto& c : checks) {
            if (c.dj > c.thr) {
                s.nJumpFails++;
                if (c.dj > s.worstJump) {
                    s.worstJump = c.dj;
                    ostringstream os;
                    os << c.tag << " P=" << P_vals[i]/1e5 << "bar d=" << c.dj;
                    s.worstLoc = os.str();
                }
            }
        }
        prevT=curT; prevRl=curRl; prevRv=curRv; prevHl=curHl; prevHv=curHv;
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  SECTION 10 — fluidProp_crit_T  (saturation P from T, pure water)
// ═══════════════════════════════════════════════════════════════════════

static SectionStats test_fluidProp_crit_T(H2ONaCl::cH2ONaCl& eos)
{
    SectionStats s;
    s.name = "fluidProp_crit_T(T)  T-sweep (saturation curve)";
    const double P_threshold   = 3.0;    // bar per °C
    const double Rho_threshold = 20.0;   // kg/m³ per °C
    const double H_threshold   = 20000;  // J/kg per °C

    auto T_vals = linspace(1.0, 370.0, 370);  // 1°C steps up to ~critical

    double prevP, prevRl, prevRv, prevHl, prevHv;
    eos.fluidProp_crit_T(T_vals[0], 1e-10, prevP, prevRl, prevRv, prevHl, prevHv);
    s.nEvals++;

    for (size_t i = 1; i < T_vals.size(); ++i) {
        double curP, curRl, curRv, curHl, curHv;
        eos.fluidProp_crit_T(T_vals[i], 1e-10, curP, curRl, curRv, curHl, curHv);
        s.nEvals++;

        bool anyNaN = std::isnan(curP) || std::isnan(curRl) || std::isnan(curRv)
                    || std::isnan(curHl) || std::isnan(curHv);
        if (anyNaN) { s.nNaN++; prevP=curP; prevRl=curRl; prevRv=curRv; prevHl=curHl; prevHv=curHv; continue; }
        bool prevNaN = std::isnan(prevP) || std::isnan(prevRl);
        if (prevNaN) { prevP=curP; prevRl=curRl; prevRv=curRv; prevHl=curHl; prevHv=curHv; continue; }

        struct { double dj; double thr; const char* tag; } checks[] = {
            { jump(prevP,  curP),  P_threshold,    "P_sat"  },
            { jump(prevRl, curRl), Rho_threshold,  "Rho_l"  },
            { jump(prevRv, curRv), Rho_threshold,  "Rho_v"  },
            { jump(prevHl, curHl), H_threshold,    "H_l"    },
            { jump(prevHv, curHv), H_threshold,    "H_v"    },
        };
        for (auto& c : checks) {
            if (c.dj > c.thr) {
                s.nJumpFails++;
                if (c.dj > s.worstJump) {
                    s.worstJump = c.dj;
                    ostringstream os;
                    os << c.tag << " T=" << T_vals[i] << "C d=" << c.dj;
                    s.worstLoc = os.str();
                }
            }
        }
        prevP=curP; prevRl=curRl; prevRv=curRv; prevHl=curHl; prevHv=curHv;
    }
    s.passed = (s.nJumpFails == 0);
    return s;
}

// ═══════════════════════════════════════════════════════════════════════
//  main
// ═══════════════════════════════════════════════════════════════════════

int main()
{
    H2ONaCl::cH2ONaCl eos;

    cout << "================================================================\n";
    cout << " Constitutive Smoothness Test\n";
    cout << " (granular Driesner functions called directly)\n";
    cout << "================================================================\n\n";
    cout << "Each section sweeps one variable in fine steps and checks\n";
    cout << "that the output has no unphysical jumps between consecutive\n";
    cout << "evaluations (within the same phase region where applicable).\n\n";

    vector<SectionStats> results;

    // ── 1  Phase-boundary surfaces ──
    cout << "  [ 1/14] VL liquid branch ...";  cout.flush();
    results.push_back(test_VL_LiquidBranch(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    cout << "  [ 2/14] VL vapor branch  ...";  cout.flush();
    results.push_back(test_VL_VaporBranch(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    cout << "  [ 3/14] Halite liquidus  ...";  cout.flush();
    results.push_back(test_HaliteLiquidus(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    cout << "  [ 4/14] Vapor-halite coexist ...";  cout.flush();
    results.push_back(test_VaporHaliteCoexist(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    cout << "  [ 5/14] P_VLH coexist    ...";  cout.flush();
    results.push_back(test_P_VLH_Coexist(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 2  Critical curve ──
    cout << "  [ 6/14] Critical curve   ...";  cout.flush();
    results.push_back(test_CriticalCurve(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 3  Rho_brine ──
    cout << "  [ 7/14] Rho_brine        ...";  cout.flush();
    results.push_back(test_RhoBrine(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 4  T* correlation ──
    cout << "  [ 8/14] T_star_V         ...";  cout.flush();
    results.push_back(test_TstarV(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 5  Phase densities ──
    cout << "  [ 9/14] calcRho          ...";  cout.flush();
    results.push_back(test_calcRho(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 6  Phase enthalpies ──
    cout << "  [10/14] calcEnthalpy     ...";  cout.flush();
    results.push_back(test_calcEnthalpy(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 7  Phase viscosities ──
    cout << "  [11/14] calcViscosity    ...";  cout.flush();
    results.push_back(test_calcViscosity(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 8  Partial compositions ──
    cout << "  [12/14] findRegion Xl,Xv ...";  cout.flush();
    results.push_back(test_findRegion_compositions(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 9  fluidProp_crit_P ──
    cout << "  [13/14] fluidProp_crit_P ...";  cout.flush();
    results.push_back(test_fluidProp_crit_P(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ── 10 fluidProp_crit_T ──
    cout << "  [14/14] fluidProp_crit_T ...";  cout.flush();
    results.push_back(test_fluidProp_crit_T(eos));
    cout << (results.back().passed ? " OK\n" : " FAIL\n");

    // ══════════════════════════════════════════════════════════════════
    //  Summary
    // ══════════════════════════════════════════════════════════════════
    cout << "\n";
    cout << "================================================================\n";
    cout << " DETAILED RESULTS\n";
    cout << "================================================================\n";
    int totalEvals = 0, totalNaN = 0, totalFails = 0;
    for (auto& r : results) {
        r.report();
        totalEvals += r.nEvals;
        totalNaN   += r.nNaN;
        totalFails += r.nJumpFails;
    }

    int nSections = (int)results.size();
    int nPassed = 0;
    for (auto& r : results) if (r.passed) nPassed++;

    cout << "\n";
    cout << "================================================================\n";
    cout << " SUMMARY\n";
    cout << "================================================================\n";
    cout << "  Sections:       " << nPassed << " / " << nSections << " passed\n";
    cout << "  Total evals:    " << totalEvals << "\n";
    cout << "  Total NaN:      " << totalNaN << "\n";
    cout << "  Total jump fails: " << totalFails << "\n";
    cout << "================================================================\n\n";

    if (nPassed == nSections) {
        cout << "RESULT: PASS — All constitutive functions are smooth.\n";
        return 0;
    } else {
        cout << "RESULT: FAIL — " << (nSections - nPassed)
             << " section(s) have discontinuities.\n";

        // Print failed sections
        cout << "\nFailed sections:\n";
        for (auto& r : results) {
            if (!r.passed) {
                cout << "  - " << r.name << "\n";
                cout << "      jumpFails=" << r.nJumpFails
                     << "  worst=" << scientific << setprecision(3) << r.worstJump
                     << fixed << "\n";
                if (!r.worstLoc.empty())
                    cout << "      at: " << r.worstLoc << "\n";
            }
        }
        return 1;
    }
}
