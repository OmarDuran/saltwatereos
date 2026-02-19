#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    const double X_test = 1.0e-4;  // 0.01 wt%
    const double X_pure = 0.0;
    
    cout << "Boundary Smoothness Test: Two-Phase to Single-Phase Vapor\n";
    cout << "==========================================================\n\n";
    cout << "Testing boundary at X=0.0001 (0.01 wt%)\n";
    cout << "Scanning pressure around saturation at multiple temperatures\n\n";
    
    vector<double> test_temps = {150.0, 200.0, 250.0, 300.0, 350.0};
    
    for(size_t t_idx = 0; t_idx < test_temps.size(); t_idx++) {
        double T_C = test_temps[t_idx];
        double T_K = T_C + 273.15;
        
        // Get saturation pressure
        double P_sat_bar = eos.m_water.P_Boiling(T_C);
        double P_sat_Pa = P_sat_bar * 1e5;
        
        cout << "\nT = " << T_C << " C, P_sat = " << fixed << setprecision(2) << P_sat_bar << " bar\n";
        cout << string(100, '=') << "\n";
        cout << setw(12) << "P/P_sat"
             << setw(10) << "P(bar)"
             << setw(10) << "Region"
             << setw(12) << "Rho(kg/m³)"
             << setw(10) << "S_v"
             << setw(10) << "S_l"
             << setw(15) << "Status" << "\n";
        cout << string(100, '-') << "\n";
        
        // Scan from below saturation (0.95*P_sat) to above (1.05*P_sat)
        for(double P_ratio = 0.90; P_ratio <= 1.10; P_ratio += 0.01) {
            double P_test_bar = P_sat_bar * P_ratio;
            double P_test_Pa = P_test_bar * 1e5;
            
            H2ONaCl::PROP_H2ONaCl prop = eos.prop_pTX(P_test_Pa, T_K, X_test, false);
            
            if(isnan(prop.Rho)) continue;
            
            string status = "";
            
            // Check for issues
            if(P_ratio < 0.99 && prop.Region != 2) {
                status = "Expected Region 2 (vapor)";
            }
            if(P_ratio > 1.01 && prop.Region != 0) {
                status = "Expected Region 0 (liquid)";
            }
            
            // Check saturation consistency
            if(prop.Region == 2) {
                if(fabs(prop.S_v - 1.0) > 0.01) status = "BAD S_v";
                if(fabs(prop.S_l) > 0.01) status = "BAD S_l";
            }
            if(prop.Region == 0) {
                if(fabs(prop.S_l - 1.0) > 0.01) status = "BAD S_l";
                if(fabs(prop.S_v) > 0.01) status = "BAD S_v";
            }
            
            // Highlight region transitions
            static int prev_region = -1;
            if(prev_region >= 0 && prev_region != prop.Region) {
                status = "→ TRANSITION";
            }
            prev_region = prop.Region;
            
            // Print every point near saturation
            if(fabs(P_ratio - 1.0) < 0.05 || !status.empty()) {
                cout << setw(12) << fixed << setprecision(4) << P_ratio
                     << setw(10) << setprecision(2) << P_test_bar
                     << setw(10) << prop.Region
                     << setw(12) << setprecision(3) << prop.Rho
                     << setw(10) << setprecision(4) << prop.S_v
                     << setw(10) << setprecision(4) << prop.S_l
                     << setw(15) << status << "\n";
            }
        }
        cout << "\n";
    }
    
    // Now test at exact saturation with IAPWS comparison
    cout << "\n\nDetailed Saturation Point Analysis\n";
    cout << "===================================\n\n";
    
    for(double T_C = 100.0; T_C <= 370.0; T_C += 20.0) {
        double T_K = T_C + 273.15;
        double P_sat_bar = eos.m_water.P_Boiling(T_C);
        double P_sat_Pa = P_sat_bar * 1e5;
        
        // Get IAPWS properties
        H2ONaCl::PROP_H2ONaCl prop_iapws = eos.prop_pTX(P_sat_Pa, T_K, X_pure, false);
        
        // Get X=0.0001 properties
        H2ONaCl::PROP_H2ONaCl prop_test = eos.prop_pTX(P_sat_Pa, T_K, X_test, false);
        
        cout << "T=" << setw(6) << fixed << setprecision(1) << T_C << "C, P=" << setw(6) << setprecision(2) << P_sat_bar << " bar: ";
        cout << "IAPWS Region=" << prop_iapws.Region << ", X=0.0001 Region=" << prop_test.Region;
        cout << ", S_v=" << setprecision(4) << prop_test.S_v << ", S_l=" << prop_test.S_l;
        
        if(prop_test.Region != 2 || fabs(prop_test.S_v - 1.0) > 0.01) {
            cout << " ⚠ ISSUE";
        }
        cout << "\n";
    }
    
    return 0;
}
