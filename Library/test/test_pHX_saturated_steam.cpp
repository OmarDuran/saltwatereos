#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    const double X_test = 1.0e-4;  // 0.0001 = 0.01 wt%
    const double X_pure = 0.0;
    
    cout << "========================================================================\n";
    cout << "Comprehensive Test: prop_pHX_bisection on Saturated Steam Curve\n";
    cout << "========================================================================\n\n";
    cout << "Testing at X = 0.0001 (0.01 wt%)\n";
    cout << "Objective: Verify S_v = 1.0 along entire saturation curve\n";
    cout << "Method: Use IAPWS H_v at saturation, call prop_pHX_bisection\n";
    cout << "Expected: Region=2 (vapor), S_v=1.0, S_l=0.0 at all points\n\n";
    
    // Test along saturation curve with dense temperature sampling
    vector<double> test_temps;
    
    // Low temperature range (100-200°C): every 10°C
    for(double T = 100.0; T <= 200.0; T += 10.0) {
        test_temps.push_back(T);
    }
    
    // Mid temperature range (200-300°C): every 5°C
    for(double T = 205.0; T <= 300.0; T += 5.0) {
        test_temps.push_back(T);
    }
    
    // High temperature range (300-370°C): every 2°C (near critical)
    for(double T = 302.0; T <= 370.0; T += 2.0) {
        test_temps.push_back(T);
    }
    
    int total_tests = 0;
    int failures = 0;
    int region_errors = 0;
    int sv_errors = 0;
    int sl_errors = 0;
    int consistency_errors = 0;
    
    double max_rho_diff = 0.0;
    double max_T_diff = 0.0;
    
    cout << "Testing " << test_temps.size() << " points along saturation curve:\n\n";
    cout << string(140, '=') << "\n";
    cout << setw(8) << "T(C)"
         << setw(10) << "P_sat(bar)"
         << setw(12) << "H_v(kJ/kg)"
         << setw(10) << "Region"
         << setw(10) << "S_v"
         << setw(10) << "S_l"
         << setw(12) << "Rho(kg/m³)"
         << setw(12) << "Rho_IAPWS"
         << setw(10) << "ΔRho(%)"
         << setw(10) << "T_conv(C)"
         << setw(10) << "ΔT(°C)"
         << setw(15) << "Status" << "\n";
    cout << string(140, '=') << "\n";
    
    for(size_t i = 0; i < test_temps.size(); i++) {
        double T_C = test_temps[i];
        double T_K = T_C + 273.15;
        
        // Get saturation pressure
        double P_sat_bar = eos.m_water.P_Boiling(T_C);
        double P_sat_Pa = P_sat_bar * 1e5;
        
        // Get IAPWS reference properties at saturation (pure water)
        H2ONaCl::PROP_H2ONaCl prop_iapws = eos.prop_pTX(P_sat_Pa, T_K, X_pure, false);
        
        if(isnan(prop_iapws.H_v) || prop_iapws.H_v <= 0) {
            continue;
        }
        
        double H_v_sat = prop_iapws.H_v;  // Vapor enthalpy at saturation
        
        // Test prop_pHX_bisection with saturated vapor enthalpy
        H2ONaCl::PROP_H2ONaCl prop_pHX;
        try {
            prop_pHX = eos.prop_pHX_bisection(P_sat_Pa, H_v_sat, X_test);
        } catch(...) {
            cout << setw(8) << fixed << setprecision(1) << T_C
                 << setw(10) << setprecision(2) << P_sat_bar
                 << setw(12) << setprecision(1) << H_v_sat/1e3
                 << "  EXCEPTION in prop_pHX\n";
            failures++;
            continue;
        }
        
        if(isnan(prop_pHX.Rho) || isnan(prop_pHX.T)) {
            cout << setw(8) << fixed << setprecision(1) << T_C
                 << setw(10) << setprecision(2) << P_sat_bar
                 << setw(12) << setprecision(1) << H_v_sat/1e3
                 << "  NaN returned\n";
            failures++;
            continue;
        }
        
        // Calculate errors
        double rho_diff = fabs(prop_pHX.Rho - prop_iapws.Rho_v) / prop_iapws.Rho_v * 100.0;
        double T_diff = fabs(prop_pHX.T - T_C);  // Temperature difference in °C (prop_pHX.T is in Celsius)
        
        max_rho_diff = max(max_rho_diff, rho_diff);
        max_T_diff = max(max_T_diff, T_diff);
        
        // Check for issues
        string status = "";
        bool has_error = false;
        
        // Check region
        if(prop_pHX.Region != 2) {
            status = "WRONG_REGION";
            region_errors++;
            has_error = true;
        }
        
        // Check S_v
        if(fabs(prop_pHX.S_v - 1.0) > 0.01) {
            if(!status.empty()) status += ",";
            status += "BAD_S_v";
            sv_errors++;
            has_error = true;
        }
        
        // Check S_l
        if(fabs(prop_pHX.S_l) > 0.01) {
            if(!status.empty()) status += ",";
            status += "BAD_S_l";
            sl_errors++;
            has_error = true;
        }
        
        // Check density accuracy
        if(rho_diff > 5.0) {
            if(!status.empty()) status += ",";
            status += "HIGH_RHO_ERR";
            has_error = true;
        }
        
        // Check temperature convergence
        if(T_diff > 1.0) {  // More than 1K difference
            if(!status.empty()) status += ",";
            status += "HIGH_T_ERR";
            has_error = true;
        }
        
        if(has_error) {
            failures++;
        }
        
        if(status.empty()) {
            status = "OK";
        }
        
        total_tests++;
        
        // Print every 5th point or if there's an issue
        if(i % 5 == 0 || has_error) {
            cout << setw(8) << fixed << setprecision(1) << T_C
                 << setw(10) << setprecision(2) << P_sat_bar
                 << setw(12) << setprecision(1) << H_v_sat/1e3
                 << setw(10) << prop_pHX.Region
                 << setw(10) << setprecision(4) << prop_pHX.S_v
                 << setw(10) << setprecision(4) << prop_pHX.S_l
                 << setw(12) << setprecision(3) << prop_pHX.Rho
                 << setw(12) << setprecision(3) << prop_iapws.Rho_v
                 << setw(10) << setprecision(4) << rho_diff
                 << setw(10) << setprecision(2) << prop_pHX.T
                 << setw(10) << setprecision(4) << T_diff
                 << setw(15) << status << "\n";
        }
    }
    
    cout << string(140, '=') << "\n\n";
    
    // Comparison test: prop_pTX vs prop_pHX consistency
    cout << "========================================================================\n";
    cout << "Consistency Test: prop_pTX vs prop_pHX_bisection\n";
    cout << "========================================================================\n\n";
    
    cout << "Testing consistency at selected temperatures:\n\n";
    cout << string(120, '=') << "\n";
    cout << setw(8) << "T(C)"
         << setw(10) << "P(bar)"
         << setw(10) << "Reg_pTX"
         << setw(10) << "Reg_pHX"
         << setw(10) << "S_v_pTX"
         << setw(10) << "S_v_pHX"
         << setw(12) << "Rho_pTX"
         << setw(12) << "Rho_pHX"
         << setw(10) << "ΔRho(%)"
         << setw(15) << "Status" << "\n";
    cout << string(120, '=') << "\n";
    
    vector<double> check_temps = {120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0};
    
    for(size_t i = 0; i < check_temps.size(); i++) {
        double T_C = check_temps[i];
        double T_K = T_C + 273.15;
        double P_sat_bar = eos.m_water.P_Boiling(T_C);
        double P_sat_Pa = P_sat_bar * 1e5;
        
        // Get IAPWS reference
        H2ONaCl::PROP_H2ONaCl prop_iapws = eos.prop_pTX(P_sat_Pa, T_K, X_pure, false);
        double H_v_sat = prop_iapws.H_v;
        
        // Test with prop_pTX
        H2ONaCl::PROP_H2ONaCl prop_pTX = eos.prop_pTX(P_sat_Pa, T_K, X_test, false);
        
        // Test with prop_pHX
        H2ONaCl::PROP_H2ONaCl prop_pHX = eos.prop_pHX(P_sat_Pa, H_v_sat, X_test);
        
        double rho_diff = fabs(prop_pTX.Rho - prop_pHX.Rho) / prop_pTX.Rho * 100.0;
        
        string status = "OK";
        if(prop_pTX.Region != prop_pHX.Region) {
            status = "REGION_MISMATCH";
            consistency_errors++;
        } else if(fabs(prop_pTX.S_v - prop_pHX.S_v) > 0.01) {
            status = "S_v_MISMATCH";
            consistency_errors++;
        } else if(rho_diff > 1.0) {
            status = "RHO_MISMATCH";
            consistency_errors++;
        }
        
        cout << setw(8) << fixed << setprecision(1) << T_C
             << setw(10) << setprecision(2) << P_sat_bar
             << setw(10) << prop_pTX.Region
             << setw(10) << prop_pHX.Region
             << setw(10) << setprecision(4) << prop_pTX.S_v
             << setw(10) << setprecision(4) << prop_pHX.S_v
             << setw(12) << setprecision(3) << prop_pTX.Rho
             << setw(12) << setprecision(3) << prop_pHX.Rho
             << setw(10) << setprecision(4) << rho_diff
             << setw(15) << status << "\n";
    }
    
    cout << string(120, '=') << "\n\n";
    
    // Summary
    cout << "========================================================================\n";
    cout << "TEST SUMMARY\n";
    cout << "========================================================================\n\n";
    
    cout << "Saturation curve test:\n";
    cout << "  Total points tested: " << total_tests << "\n";
    cout << "  Points with issues: " << failures << "\n";
    cout << "    - Region errors: " << region_errors << "\n";
    cout << "    - S_v errors (S_v ≠ 1.0): " << sv_errors << "\n";
    cout << "    - S_l errors (S_l ≠ 0.0): " << sl_errors << "\n";
    cout << "  Max density error: " << fixed << setprecision(4) << max_rho_diff << "%\n";
    cout << "  Max temperature error: " << setprecision(4) << max_T_diff << " K\n\n";
    
    cout << "Consistency test:\n";
    cout << "  prop_pTX vs prop_pHX mismatches: " << consistency_errors << "\n\n";
    
    int total_failures = failures + consistency_errors;
    
    if(total_failures == 0) {
        cout << "✓✓✓ ALL TESTS PASSED ✓✓✓\n";
        cout << "prop_pHX_bisection correctly identifies saturated steam (S_v=1.0)\n";
        cout << "along the entire saturation curve for X=0.0001\n";
        return 0;
    } else {
        cout << "⚠ TEST FAILED: " << total_failures << " issue(s) detected\n";
        cout << "prop_pHX_bisection has problems with saturated steam identification\n";
        return 1;
    }
}
