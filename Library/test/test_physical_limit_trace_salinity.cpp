#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

struct TestCondition {
    double P_bar;
    double T_C;
    string description;
    string region_type;
};

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Physical Limit Test: Very Low Salinity → Pure Water\n";
    cout << "====================================================\n\n";
    
    cout << "Objective: Verify that X=1e-5 reproduces pure water properties\n";
    cout << "Method: Compare H2ONaCl(X=1e-5) vs H2ONaCl(X=0) using IAPWS\n";
    cout << "Expected: Properties should match within O(1e-5) relative error\n\n";
    
    cout << "NOTE: Pure water (X=0) uses IAPWS directly, which may identify\n";
    cout << "      regions differently than the Driesner & Heinrich correlations.\n";
    cout << "      Tests are only performed where region identification matches.\n";
    cout << "      The key test is that thermodynamic properties (rho, H) converge.\n\n";
    
    const double X_trace = 1.0e-5;  // Trace salinity
    const double X_pure = 0.0;       // Pure water
    const double rel_tol = 1.0e-4;   // 0.01% tolerance (conservative)
    
    // Define test conditions covering different regions
    // NOTE: Avoid saturation conditions as phase boundaries shift slightly with salinity
    // Focus on single-phase regions far from boundaries
    vector<TestCondition> conditions = {
        // Subcooled liquid (well below saturation)
        {50.0,  50.0,  "Subcooled liquid, low T", "liquid"},
        {100.0, 100.0, "Subcooled liquid, moderate T", "liquid"},
        {200.0, 250.0, "Subcooled liquid, high P", "liquid"},
        {300.0, 300.0, "Subcooled liquid, very high P", "liquid"},
        
        // Superheated vapor (well above saturation)
        {10.0,  400.0, "Superheated vapor, low P", "vapor"},
        {20.0,  500.0, "Superheated vapor, moderate P", "vapor"},
        {50.0,  600.0, "Superheated vapor, high T", "vapor"},
        
        // Supercritical region (above critical point)
        {250.0, 400.0, "Supercritical, near critical", "supercritical"},
        {300.0, 500.0, "Supercritical, moderate", "supercritical"},
        {400.0, 600.0, "Supercritical, high P/T", "supercritical"},
        {500.0, 700.0, "Supercritical, very high P/T", "supercritical"}
    };
    
    int total_tests = 0;
    int failed_tests = 0;
    int skipped_tests = 0;
    double max_rho_error = 0.0;
    double max_h_error = 0.0;
    double max_mu_error = 0.0;
    
    cout << "Test Results:\n";
    cout << "=============\n\n";
    
    for(size_t i = 0; i < conditions.size(); i++) {
        TestCondition tc = conditions[i];
        double P_Pa = tc.P_bar * 1e5;
        double T_K = tc.T_C + 273.15;
        
        cout << "Test " << (i+1) << ": " << tc.description << " (" << tc.region_type << ")\n";
        cout << "  P = " << tc.P_bar << " bar, T = " << tc.T_C << " C\n";
        
        // Get properties for pure water (reference)
        H2ONaCl::PROP_H2ONaCl prop_pure = eos.prop_pTX(P_Pa, T_K, X_pure, true);
        
        // Get properties for trace salinity
        H2ONaCl::PROP_H2ONaCl prop_trace = eos.prop_pTX(P_Pa, T_K, X_trace, true);
        
        if(isnan(prop_pure.Rho) || isnan(prop_trace.Rho)) {
            cout << "  SKIPPED: Invalid properties returned\n\n";
            continue;
        }
        
        // Check if regions match (important for physical limit test)
        if(prop_pure.Region != prop_trace.Region) {
            cout << "  SKIPPED: Different regions identified (Pure: " << prop_pure.Region
                 << ", Trace: " << prop_trace.Region << ")\n";
            cout << "           This indicates the conditions are near a phase boundary.\n\n";
            skipped_tests++;
            continue;
        }
        
        cout << "  Region: " << prop_pure.Region << " (" << eos.getPhaseRegionName(prop_pure.Region) << ")\n";
        
        cout << "  Pure water (X=0):      Rho=" << fixed << setprecision(3) << prop_pure.Rho
             << " kg/m3, H=" << setprecision(1) << prop_pure.H/1e3 << " kJ/kg";
        if(prop_pure.Mu_l > 0) cout << ", Mu_l=" << scientific << setprecision(3) << prop_pure.Mu_l << " Pa·s";
        if(prop_pure.Mu_v > 0) cout << ", Mu_v=" << setprecision(3) << prop_pure.Mu_v << " Pa·s";
        cout << "\n";
        
        cout << "  Trace salinity (X=1e-5): Rho=" << fixed << setprecision(3) << prop_trace.Rho
             << " kg/m3, H=" << setprecision(1) << prop_trace.H/1e3 << " kJ/kg";
        if(prop_trace.Mu_l > 0) cout << ", Mu_l=" << scientific << setprecision(3) << prop_trace.Mu_l << " Pa·s";
        if(prop_trace.Mu_v > 0) cout << ", Mu_v=" << setprecision(3) << prop_trace.Mu_v << " Pa·s";
        cout << "\n";
        
        // Calculate relative errors
        double rho_err = fabs(prop_trace.Rho - prop_pure.Rho) / prop_pure.Rho;
        double h_err = fabs(prop_trace.H - prop_pure.H) / fabs(prop_pure.H);
        
        cout << "  Relative errors: Rho=" << scientific << setprecision(3) << rho_err
             << " (" << setprecision(2) << rho_err*100 << "%)";
        cout << ", H=" << setprecision(3) << h_err << " (" << setprecision(2) << h_err*100 << "%)\n";
        
        // Check viscosities if available
        bool mu_l_valid = (prop_pure.Mu_l > 0 && prop_trace.Mu_l > 0);
        bool mu_v_valid = (prop_pure.Mu_v > 0 && prop_trace.Mu_v > 0);
        double mu_l_err = 0.0, mu_v_err = 0.0;
        
        if(mu_l_valid) {
            mu_l_err = fabs(prop_trace.Mu_l - prop_pure.Mu_l) / prop_pure.Mu_l;
            cout << "                   Mu_l=" << setprecision(3) << mu_l_err
                 << " (" << setprecision(2) << mu_l_err*100 << "%)";
        }
        if(mu_v_valid) {
            mu_v_err = fabs(prop_trace.Mu_v - prop_pure.Mu_v) / prop_pure.Mu_v;
            if(mu_l_valid) cout << ", ";
            cout << "Mu_v=" << setprecision(3) << mu_v_err
                 << " (" << setprecision(2) << mu_v_err*100 << "%)";
        }
        if(mu_l_valid || mu_v_valid) cout << "\n";
        
        total_tests++;
        bool test_passed = true;
        
        // Check if errors are within tolerance
        if(rho_err > rel_tol) {
            cout << "  FAILED: Density error exceeds tolerance!\n";
            failed_tests++;
            test_passed = false;
        }
        if(h_err > rel_tol) {
            cout << "  FAILED: Enthalpy error exceeds tolerance!\n";
            if(test_passed) failed_tests++;
            test_passed = false;
        }
        if(mu_l_valid && mu_l_err > rel_tol) {
            cout << "  WARNING: Liquid viscosity error exceeds tolerance\n";
        }
        if(mu_v_valid && mu_v_err > rel_tol) {
            cout << "  WARNING: Vapor viscosity error exceeds tolerance\n";
        }
        
        if(test_passed) {
            cout << "  PASSED\n";
        }
        
        max_rho_error = max(max_rho_error, rho_err);
        max_h_error = max(max_h_error, h_err);
        if(mu_l_valid) max_mu_error = max(max_mu_error, mu_l_err);
        if(mu_v_valid) max_mu_error = max(max_mu_error, mu_v_err);
        
        cout << "\n";
    }
    
    // Now test prop_pHX_bisection
    cout << "====================================================\n";
    cout << "Testing prop_pHX_bisection with trace salinity\n";
    cout << "====================================================\n\n";
    
    int pHX_tests = 0;
    int pHX_failed = 0;
    
    // Expanded test grid covering P-H Cartesian space more comprehensively
    // Pressure range: 10 - 800 bar
    // Temperature range: 25 - 800 C (to cover wide enthalpy range)
    vector<TestCondition> pHX_conditions = {
        // Low pressure region (10-50 bar) - Subcooled liquid to superheated vapor
        {10.0,  50.0,  "Low P - Subcooled liquid", "liquid"},
        {10.0,  180.0, "Low P - Near saturation liquid", "liquid"},
        {10.0,  300.0, "Low P - Superheated vapor", "vapor"},
        {10.0,  500.0, "Low P - High superheat vapor", "vapor"},
        {10.0,  700.0, "Low P - Very high superheat vapor", "vapor"},
        
        {25.0,  50.0,  "25 bar - Subcooled liquid", "liquid"},
        {25.0,  150.0, "25 bar - Moderate subcooled liquid", "liquid"},
        {25.0,  350.0, "25 bar - Superheated vapor", "vapor"},
        {25.0,  600.0, "25 bar - High superheat vapor", "vapor"},
        
        {50.0,  50.0,  "50 bar - Cold subcooled liquid", "liquid"},
        {50.0,  150.0, "50 bar - Warm subcooled liquid", "liquid"},
        {50.0,  250.0, "50 bar - Hot subcooled liquid", "liquid"},
        {50.0,  400.0, "50 bar - Superheated vapor", "vapor"},
        {50.0,  650.0, "50 bar - High superheat vapor", "vapor"},
        
        // Moderate pressure region (75-150 bar)
        {75.0,  100.0, "75 bar - Subcooled liquid", "liquid"},
        {75.0,  200.0, "75 bar - Warm liquid", "liquid"},
        {75.0,  450.0, "75 bar - Superheated vapor", "vapor"},
        {75.0,  700.0, "75 bar - High superheat vapor", "vapor"},
        
        {100.0, 50.0,  "100 bar - Cold liquid", "liquid"},
        {100.0, 100.0, "100 bar - Subcooled liquid", "liquid"},
        {100.0, 200.0, "100 bar - Warm liquid", "liquid"},
        {100.0, 300.0, "100 bar - Hot liquid/near critical", "liquid"},
        {100.0, 500.0, "100 bar - Superheated vapor", "vapor"},
        {100.0, 750.0, "100 bar - Very high superheat", "vapor"},
        
        {150.0, 100.0, "150 bar - Subcooled liquid", "liquid"},
        {150.0, 250.0, "150 bar - Hot liquid", "liquid"},
        {150.0, 450.0, "150 bar - Vapor/supercritical", "vapor"},
        {150.0, 700.0, "150 bar - High temperature", "vapor"},
        
        // Near-critical and supercritical region (200-300 bar)
        {200.0, 100.0, "200 bar - Cold liquid", "liquid"},
        {200.0, 200.0, "200 bar - Warm liquid", "liquid"},
        {200.0, 350.0, "200 bar - Near critical", "supercritical"},
        {200.0, 500.0, "200 bar - Supercritical", "supercritical"},
        {200.0, 700.0, "200 bar - High T supercritical", "supercritical"},
        
        {250.0, 150.0, "250 bar - Subcooled liquid", "liquid"},
        {250.0, 300.0, "250 bar - Hot liquid/critical", "supercritical"},
        {250.0, 400.0, "250 bar - Supercritical", "supercritical"},
        {250.0, 550.0, "250 bar - High T supercritical", "supercritical"},
        {250.0, 750.0, "250 bar - Very high T", "supercritical"},
        
        {300.0, 200.0, "300 bar - Warm liquid", "liquid"},
        {300.0, 350.0, "300 bar - Supercritical", "supercritical"},
        {300.0, 500.0, "300 bar - Moderate supercritical", "supercritical"},
        {300.0, 700.0, "300 bar - High T supercritical", "supercritical"},
        
        // High pressure region (400-600 bar)
        {400.0, 150.0, "400 bar - Cold compressed liquid", "liquid"},
        {400.0, 300.0, "400 bar - Hot liquid", "liquid"},
        {400.0, 450.0, "400 bar - Supercritical", "supercritical"},
        {400.0, 600.0, "400 bar - High supercritical", "supercritical"},
        {400.0, 800.0, "400 bar - Very high T", "supercritical"},
        
        {500.0, 200.0, "500 bar - Warm compressed liquid", "liquid"},
        {500.0, 350.0, "500 bar - Hot liquid/supercritical", "supercritical"},
        {500.0, 550.0, "500 bar - Supercritical", "supercritical"},
        {500.0, 750.0, "500 bar - High T supercritical", "supercritical"},
        
        {600.0, 250.0, "600 bar - Compressed liquid", "liquid"},
        {600.0, 400.0, "600 bar - Hot supercritical", "supercritical"},
        {600.0, 600.0, "600 bar - High supercritical", "supercritical"},
        {600.0, 800.0, "600 bar - Very high T", "supercritical"},
        
        // Very high pressure region (700-800 bar)
        {700.0, 300.0, "700 bar - Hot compressed liquid", "liquid"},
        {700.0, 500.0, "700 bar - Supercritical", "supercritical"},
        {700.0, 700.0, "700 bar - High T supercritical", "supercritical"},
        
        {800.0, 350.0, "800 bar - Very high P liquid", "liquid"},
        {800.0, 550.0, "800 bar - Very high P supercritical", "supercritical"},
        {800.0, 750.0, "800 bar - Extreme conditions", "supercritical"}
    };
    
    for(size_t i = 0; i < pHX_conditions.size(); i++) {
        TestCondition tc = pHX_conditions[i];
        double P_Pa = tc.P_bar * 1e5;
        double T_K = tc.T_C + 273.15;
        
        cout << "pHX Test " << (i+1) << ": " << tc.description << "\n";
        cout << "  P = " << tc.P_bar << " bar, T = " << tc.T_C << " C\n";
        
        // Get enthalpy at this condition for pure water
        H2ONaCl::PROP_H2ONaCl prop_pure_pTX = eos.prop_pTX(P_Pa, T_K, X_pure, false);
        double H_target = prop_pure_pTX.H;
        
        if(isnan(H_target) || H_target == 0.0) {
            cout << "  SKIPPED: Invalid H from pure water\n\n";
            continue;
        }
        
        // Use prop_pHX_bisection with trace salinity
        H2ONaCl::PROP_H2ONaCl prop_trace_pHX = eos.prop_pHX_bisection(P_Pa, H_target, X_trace);
        
        // Use prop_pHX_bisection with pure water
        H2ONaCl::PROP_H2ONaCl prop_pure_pHX = eos.prop_pHX_bisection(P_Pa, H_target, X_pure);
        
        cout << "  H_target = " << fixed << setprecision(2) << H_target/1e3 << " kJ/kg\n";
        cout << "  Pure water:      T=" << setprecision(2) << prop_pure_pHX.T << " C, Rho="
             << setprecision(3) << prop_pure_pHX.Rho << " kg/m3\n";
        cout << "  Trace salinity:  T=" << setprecision(2) << prop_trace_pHX.T << " C, Rho="
             << setprecision(3) << prop_trace_pHX.Rho << " kg/m3\n";
        
        double T_err = fabs(prop_trace_pHX.T - prop_pure_pHX.T) / (prop_pure_pHX.T + 273.15);
        double rho_err = fabs(prop_trace_pHX.Rho - prop_pure_pHX.Rho) / prop_pure_pHX.Rho;
        
        cout << "  Relative errors: T=" << scientific << setprecision(3) << T_err
             << " (" << setprecision(2) << T_err*100 << "%)";
        cout << ", Rho=" << setprecision(3) << rho_err
             << " (" << setprecision(2) << rho_err*100 << "%)\n";
        
        pHX_tests++;
        
        if(T_err > rel_tol || rho_err > rel_tol) {
            cout << "  FAILED: Errors exceed tolerance!\n";
            pHX_failed++;
        } else {
            cout << "  PASSED\n";
        }
        
        cout << "\n";
    }
    
    // NEW: Viscosity Continuity Test
    cout << "====================================================\n";
    cout << "Testing Viscosity Continuity Across Phase Regions\n";
    cout << "====================================================\n\n";
    cout << "Checking for discontinuous jumps in Mu_v (vapor viscosity)\n";
    cout << "across two-phase and single-phase regions.\n\n";
    
    int visc_tests = 0;
    int visc_jumps_detected = 0;
    double max_mu_v_jump = 0.0;
    
    // Test 1: Temperature scan at moderate pressure
    cout << "Viscosity Test 1: Temperature scan at P=100 bar, X=0.05\n";
    cout << "---------------------------------------------------------------\n";
    cout << setw(10) << "T (°C)"
         << setw(12) << "Region"
         << setw(15) << "Mu_v (Pa·s)"
         << setw(15) << "ΔMu_v (%)"
         << setw(12) << "S_l"
         << setw(12) << "S_v" << "\n";
    cout << string(80, '-') << "\n";
    
    double P_visc1 = 100.0e5;  // 100 bar
    double X_visc1 = 0.05;     // 5% salinity
    double prev_Mu_v = 0.0;
    int prev_region = -1;
    
    for(double T = 200.0; T <= 400.0; T += 5.0) {
        double T_K = T + 273.15;
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pTX(P_visc1, T_K, X_visc1, true);
        
        double delta_Mu_v = 0.0;
        if(prev_Mu_v > 1e-12 && prop.Mu_v > 1e-12) {
            delta_Mu_v = 100.0 * (prop.Mu_v - prev_Mu_v) / prev_Mu_v;
            
            // Flag significant jumps (> 15%) especially at region transitions
            if(std::abs(delta_Mu_v) > 15.0) {
                visc_jumps_detected++;
                max_mu_v_jump = max(max_mu_v_jump, std::abs(delta_Mu_v));
            }
        }
        
        cout << setw(10) << fixed << setprecision(1) << T
             << setw(12) << prop.Region
             << setw(15) << scientific << setprecision(4) << prop.Mu_v
             << setw(15) << fixed << setprecision(2) << delta_Mu_v
             << setw(12) << setprecision(4) << prop.S_l
             << setw(12) << setprecision(4) << prop.S_v;
        
        if(std::abs(delta_Mu_v) > 15.0 && prev_region >= 0) {
            cout << "  *** JUMP at region " << prev_region << "→" << prop.Region << " ***";
        }
        cout << "\n";
        
        prev_Mu_v = prop.Mu_v;
        prev_region = prop.Region;
        visc_tests++;
    }
    cout << "\n";
    
    // Test 2: Salinity scan at two-phase conditions
    cout << "Viscosity Test 2: Salinity scan at P=50 bar, T=350 °C\n";
    cout << "---------------------------------------------------------------\n";
    cout << setw(10) << "X_wt"
         << setw(12) << "Region"
         << setw(15) << "Mu_v (Pa·s)"
         << setw(15) << "ΔMu_v (%)"
         << setw(12) << "S_l"
         << setw(12) << "S_v" << "\n";
    cout << string(80, '-') << "\n";
    
    double P_visc2 = 50.0e5;    // 50 bar
    double T_visc2 = 350.0 + 273.15;  // 350°C
    prev_Mu_v = 0.0;
    prev_region = -1;
    
    for(double X = 0.0; X <= 0.2501; X += 0.01) {
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pTX(P_visc2, T_visc2, X, true);
        
        double delta_Mu_v = 0.0;
        if(prev_Mu_v > 1e-12 && prop.Mu_v > 1e-12) {
            delta_Mu_v = 100.0 * (prop.Mu_v - prev_Mu_v) / prev_Mu_v;
            
            if(std::abs(delta_Mu_v) > 15.0) {
                visc_jumps_detected++;
                max_mu_v_jump = max(max_mu_v_jump, std::abs(delta_Mu_v));
            }
        }
        
        cout << setw(10) << fixed << setprecision(4) << X
             << setw(12) << prop.Region
             << setw(15) << scientific << setprecision(4) << prop.Mu_v
             << setw(15) << fixed << setprecision(2) << delta_Mu_v
             << setw(12) << setprecision(4) << prop.S_l
             << setw(12) << setprecision(4) << prop.S_v;
        
        if(std::abs(delta_Mu_v) > 15.0 && prev_region >= 0) {
            cout << "  *** JUMP at region " << prev_region << "→" << prop.Region << " ***";
        }
        cout << "\n";
        
        prev_Mu_v = prop.Mu_v;
        prev_region = prop.Region;
        visc_tests++;
    }
    cout << "\n";
    
    // Test 3: Near critical point scan (trace salinity)
    cout << "Viscosity Test 3: Near critical point at P=220 bar, X=0.0001\n";
    cout << "---------------------------------------------------------------\n";
    cout << setw(10) << "T (°C)"
         << setw(12) << "Region"
         << setw(15) << "Mu_v (Pa·s)"
         << setw(15) << "ΔMu_v (%)" << "\n";
    cout << string(60, '-') << "\n";
    
    double P_visc3 = 220.0e5;   // 220 bar (near critical)
    double X_visc3 = 0.0001;    // Trace salinity
    prev_Mu_v = 0.0;
    prev_region = -1;
    
    for(double T = 360.0; T <= 390.0; T += 1.0) {
        double T_K = T + 273.15;
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pTX(P_visc3, T_K, X_visc3, true);
        
        double delta_Mu_v = 0.0;
        if(prev_Mu_v > 1e-12 && prop.Mu_v > 1e-12) {
            delta_Mu_v = 100.0 * (prop.Mu_v - prev_Mu_v) / prev_Mu_v;
            
            if(std::abs(delta_Mu_v) > 15.0) {
                visc_jumps_detected++;
                max_mu_v_jump = max(max_mu_v_jump, std::abs(delta_Mu_v));
            }
        }
        
        cout << setw(10) << fixed << setprecision(1) << T
             << setw(12) << prop.Region
             << setw(15) << scientific << setprecision(4) << prop.Mu_v
             << setw(15) << fixed << setprecision(2) << delta_Mu_v;
        
        if(std::abs(delta_Mu_v) > 15.0 && prev_region >= 0) {
            cout << "  *** JUMP at region " << prev_region << "→" << prop.Region << " ***";
        }
        cout << "\n";
        
        prev_Mu_v = prop.Mu_v;
        prev_region = prop.Region;
        visc_tests++;
    }
    cout << "\n";
    
    cout << "Viscosity continuity test summary:\n";
    cout << "  Total viscosity checks: " << visc_tests << "\n";
    cout << "  Discontinuous jumps detected (>15%): " << visc_jumps_detected << "\n";
    if(visc_jumps_detected > 0) {
        cout << "  Maximum viscosity jump: " << fixed << setprecision(1) << max_mu_v_jump << "%\n";
        cout << "\n  ⚠ WARNING: Vapor viscosity has discontinuities!\n";
        cout << "  This indicates jumps in Mu_v at phase boundaries.\n";
    } else {
        cout << "\n  ✓ SUCCESS: No significant viscosity discontinuities detected\n";
    }
    cout << "\n";
    
    // NEW: Boundary Smoothness Test - Steam/Two-Phase Transition
    cout << "====================================================\n";
    cout << "Boundary Smoothness Test: Steam → Two-Phase Transition\n";
    cout << "====================================================\n\n";
    cout << "Objective: Verify smoothness of the boundary between pure steam\n";
    cout << "           and two-phase regions at trace salinity (X=1e-4)\n";
    cout << "Method: Temperature scan at constant pressure, comparing with IAPWS\n";
    cout << "Expected: Properties should transition smoothly across phase boundary\n\n";
    
    const double X_boundary = 1.0e-4;  // Trace salinity for boundary test
    int boundary_tests = 0;
    int boundary_issues = 0;
    
    // Test at a single representative pressure
    double P_bar = 100.0;
    double P_Pa = P_bar * 1e5;
    
    cout << "\nBoundary Test at P = " << P_bar << " bar\n";
    cout << "----------------------------------------\n";
    cout << setw(10) << "T (°C)"
         << setw(15) << "Region (EOS)"
         << setw(15) << "Region (IAPWS)"
         << setw(15) << "Rho (kg/m³)"
         << setw(15) << "Rho_IAPWS"
         << setw(15) << "ΔRho (%)"
         << setw(15) << "H (kJ/kg)"
         << setw(15) << "H_IAPWS"
         << setw(15) << "ΔH (%)" << "\n";
    cout << string(130, '-') << "\n";
    
    // Temperature scan from well below to well above saturation
    double T_start = 250.0;
    double T_end = 400.0;
    double T_step = 5.0;
    
    // Store data for derivative analysis
    vector<double> temps;
    vector<int> regions_eos;
    vector<int> regions_iapws;
    vector<double> rhos_eos;
    vector<double> rhos_iapws;
    vector<double> hs_eos;
    vector<double> hs_iapws;
    
    for(double T = T_start; T <= T_end; T += T_step) {
        double T_K = T + 273.15;
        
        // Get properties with trace salinity
        H2ONaCl::PROP_H2ONaCl prop_trace = eos.prop_pTX(P_Pa, T_K, X_boundary, true);
        
        // Get properties with pure water (IAPWS reference)
        H2ONaCl::PROP_H2ONaCl prop_iapws = eos.prop_pTX(P_Pa, T_K, X_pure, true);
        
        if(isnan(prop_trace.Rho) || isnan(prop_iapws.Rho)) {
            continue;
        }
        
        temps.push_back(T);
        regions_eos.push_back(prop_trace.Region);
        regions_iapws.push_back(prop_iapws.Region);
        rhos_eos.push_back(prop_trace.Rho);
        rhos_iapws.push_back(prop_iapws.Rho);
        hs_eos.push_back(prop_trace.H);
        hs_iapws.push_back(prop_iapws.H);
        
        double rho_err = fabs(prop_trace.Rho - prop_iapws.Rho) / prop_iapws.Rho * 100.0;
        double h_err = fabs(prop_trace.H - prop_iapws.H) / fabs(prop_iapws.H) * 100.0;
        
        cout << setw(10) << fixed << setprecision(1) << T
             << setw(15) << prop_trace.Region
             << setw(15) << prop_iapws.Region
             << setw(15) << setprecision(3) << prop_trace.Rho
             << setw(15) << setprecision(3) << prop_iapws.Rho
             << setw(15) << setprecision(4) << rho_err
             << setw(15) << setprecision(2) << prop_trace.H/1e3
             << setw(15) << setprecision(2) << prop_iapws.H/1e3
             << setw(15) << setprecision(4) << h_err;
        
        // Flag transitions
        if(temps.size() > 1 && regions_eos[regions_eos.size()-2] != prop_trace.Region) {
            cout << "  ← Transition";
        }
        
        // Check for excessive errors near boundaries
        if(rho_err > 1.0 || h_err > 1.0) {
            cout << "  ⚠";
            boundary_issues++;
        }
        
        cout << "\n";
        boundary_tests++;
    }
    
    // Analyze smoothness: check for jumps in derivatives
    cout << "\nSmoothness Analysis (derivatives):\n";
    cout << setw(10) << "T (°C)"
         << setw(15) << "dRho/dT"
         << setw(15) << "dRho/dT_IAPWS"
         << setw(15) << "d²Rho/dT²"
         << setw(15) << "dH/dT"
         << setw(15) << "dH/dT_IAPWS" << "\n";
    cout << string(90, '-') << "\n";
    
    for(size_t i = 1; i < temps.size() - 1; i++) {
        // First derivative (central difference)
        double dT = temps[i+1] - temps[i-1];
        double drho_dT_eos = (rhos_eos[i+1] - rhos_eos[i-1]) / dT;
        double drho_dT_iapws = (rhos_iapws[i+1] - rhos_iapws[i-1]) / dT;
        double dh_dT_eos = (hs_eos[i+1] - hs_eos[i-1]) / dT;
        double dh_dT_iapws = (hs_iapws[i+1] - hs_iapws[i-1]) / dT;
        
        // Second derivative for rho
        double d2rho_dT2 = 0.0;
        if(i > 1 && i < temps.size() - 2) {
            double dT_sq = T_step * T_step;
            d2rho_dT2 = (rhos_eos[i+1] - 2.0*rhos_eos[i] + rhos_eos[i-1]) / dT_sq;
        }
        
        cout << setw(10) << fixed << setprecision(1) << temps[i]
             << setw(15) << setprecision(4) << drho_dT_eos
             << setw(15) << setprecision(4) << drho_dT_iapws
             << setw(15) << setprecision(4) << d2rho_dT2
             << setw(15) << setprecision(2) << dh_dT_eos
             << setw(15) << setprecision(2) << dh_dT_iapws;
        
        // Flag large second derivatives (indicating non-smoothness)
        if(fabs(d2rho_dT2) > 5.0) {
            cout << "  ⚠ Large curvature";
        }
        
        cout << "\n";
    }
    
    cout << "\nBoundary smoothness test summary:\n";
    cout << "  Total boundary points checked: " << boundary_tests << "\n";
    cout << "  Points with large errors (>1%): " << boundary_issues << "\n";
    if(boundary_issues > 0) {
        cout << "\n  ⚠ WARNING: Boundary region shows discontinuities or large errors!\n";
        cout << "  The steam/two-phase transition may not be smooth at trace salinity.\n";
    } else {
        cout << "\n  ✓ SUCCESS: Boundary transition appears smooth\n";
    }
    cout << "\n";
    
    // Summary
    cout << "====================================================\n";
    cout << "Summary:\n";
    cout << "  prop_pTX tests: " << total_tests << " executed, "
         << (total_tests - failed_tests) << " passed, " << failed_tests << " failed, "
         << skipped_tests << " skipped (region mismatch)\n";
    cout << "  prop_pHX_bisection tests: " << pHX_tests << " executed, "
         << (pHX_tests - pHX_failed) << " passed, " << pHX_failed << " failed\n";
    cout << "  Maximum density error: " << scientific << setprecision(3) << max_rho_error
         << " (" << setprecision(2) << max_rho_error*100 << "%)\n";
    cout << "  Maximum enthalpy error: " << setprecision(3) << max_h_error
         << " (" << setprecision(2) << max_h_error*100 << "%)\n";
    cout << "  Maximum viscosity error: " << setprecision(3) << max_mu_error
         << " (" << setprecision(2) << max_mu_error*100 << "%)\n";
    
    cout << "\n====================================================\n";
    
    // Return status
    int total_failed = failed_tests + pHX_failed;
    if(total_failed > 0 || boundary_issues > 0) {
        cout << "TEST FAILED: " << total_failed << " test(s) failed, " 
             << boundary_issues << " boundary issue(s) detected\n";
        return 1;
    } else {
        cout << "ALL TESTS PASSED\n";
        return 0;
    }
}
