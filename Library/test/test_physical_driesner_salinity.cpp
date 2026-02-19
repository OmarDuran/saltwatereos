#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Driesner Correlation Integrity Test: prop_pHX\n";
    cout << "==============================================\n\n";
    
    cout << "Objective: Verify integrity of Driesner correlations for X > 0.0001\n";
    cout << "Method: Scan P-H-X space and ensure no NaN or invalid properties\n";
    cout << "Testing: prop_pHX(P, H, X_wt)\n\n";
    
    // Define test grid
    // Salinity: 10 points from 0.0001 to 0.2 (0.01% to 20%)
    vector<double> X_values;
    for(int i = 0; i < 10; i++) {
        double X = 0.0001 + i * (0.2 - 0.0001) / 9.0;
        X_values.push_back(X);
    }
    
    // Pressure: Range from 6 to 510 MPa (60 to 5100 bar)
    // Use logarithmic spacing for better coverage
    vector<double> P_values_bar;
    double P_min_bar = 60.0;   // 6 MPa
    double P_max_bar = 5100.0; // 510 MPa
    for(int i = 0; i < 10; i++) {
        double log_P = log10(P_min_bar) + i * (log10(P_max_bar) - log10(P_min_bar)) / 9.0;
        P_values_bar.push_back(pow(10.0, log_P));
    }
    
    // Enthalpy: 10 points from 500 to 3500 kJ/kg
    vector<double> H_values_kJ;
    for(int i = 0; i < 10; i++) {
        double H_kJ = 500.0 + i * (3500.0 - 500.0) / 9.0;
        H_values_kJ.push_back(H_kJ);
    }
    
    cout << "Test Grid:\n";
    cout << "  Salinity (X):  " << X_values.size() << " points from " 
         << fixed << setprecision(4) << X_values[0] << " to " << X_values.back() << "\n";
    cout << "  Pressure (P):  " << P_values_bar.size() << " points from " 
         << setprecision(1) << P_values_bar[0] << " to " << P_values_bar.back() << " bar\n";
    cout << "  Enthalpy (H):  " << H_values_kJ.size() << " points from " 
         << setprecision(0) << H_values_kJ[0] << " to " << H_values_kJ.back() << " kJ/kg\n";
    cout << "  Total tests:   " << (X_values.size() * P_values_bar.size() * H_values_kJ.size()) << "\n\n";
    
    int total_tests = 0;
    int failed_tests = 0;
    int nan_tests = 0;
    int invalid_region_tests = 0;
    int converged_tests = 0;
    
    cout << "Running tests...\n";
    cout << string(100, '=') << "\n";
    cout << setw(6) << "Test#"
         << setw(8) << "X"
         << setw(10) << "P(bar)"
         << setw(10) << "H(kJ/kg)"
         << setw(10) << "T(C)"
         << setw(12) << "Rho(kg/m3)"
         << setw(8) << "Region"
         << setw(12) << "Status"
         << "\n";
    cout << string(100, '-') << "\n";
    
    int test_num = 0;
    int print_interval = 100; // Print every 100th test or failures
    
    for(size_t iX = 0; iX < X_values.size(); iX++) {
        double X_wt = X_values[iX];
        
        for(size_t iP = 0; iP < P_values_bar.size(); iP++) {
            double P_bar = P_values_bar[iP];
            double P_Pa = P_bar * 1e5;
            
            for(size_t iH = 0; iH < H_values_kJ.size(); iH++) {
                double H_kJ = H_values_kJ[iH];
                double H_J = H_kJ * 1e3;
                
                test_num++;
                total_tests++;
                
                // Call prop_pHX
                H2ONaCl::PROP_H2ONaCl prop;
                bool test_failed = false;
                string status = "OK";
                
                try {
                    prop = eos.prop_pHX(P_Pa, H_J, X_wt);
                    
                    // Check for NaN in critical properties
                    if(isnan(prop.Rho) || isnan(prop.T) || isnan(prop.H)) {
                        nan_tests++;
                        failed_tests++;
                        test_failed = true;
                        status = "NaN";
                    }
                    // Check for invalid density (negative or zero)
                    else if(prop.Rho <= 0.0) {
                        failed_tests++;
                        test_failed = true;
                        status = "RHO<=0";
                    }
                    // Check for invalid temperature (outside valid range)
                    else if(prop.T < -50.0 || prop.T > 1100.0) {
                        failed_tests++;
                        test_failed = true;
                        status = "T_RANGE";
                    }
                    // Check for invalid region
                    else if(prop.Region < 0 || prop.Region > 10) {
                        invalid_region_tests++;
                        failed_tests++;
                        test_failed = true;
                        status = "BAD_REG";
                    }
                    else {
                        converged_tests++;
                    }
                    
                } catch(...) {
                    failed_tests++;
                    test_failed = true;
                    status = "EXCEPTION";
                    prop.T = 0.0;
                    prop.Rho = 0.0;
                    prop.Region = -1;
                }
                
                // Print results for failed tests or at intervals
                if(test_failed || (test_num % print_interval == 0)) {
                    cout << setw(6) << test_num
                         << setw(8) << fixed << setprecision(4) << X_wt
                         << setw(10) << setprecision(1) << P_bar
                         << setw(10) << setprecision(0) << H_kJ
                         << setw(10) << setprecision(2) << prop.T
                         << setw(12) << setprecision(3) << prop.Rho
                         << setw(8) << prop.Region
                         << setw(12) << status
                         << "\n";
                }
            }
        }
    }
    
    cout << string(100, '=') << "\n\n";
    
    // Summary statistics
    cout << "Test Summary:\n";
    cout << "=============\n";
    cout << "  Total tests:           " << total_tests << "\n";
    cout << "  Converged successfully: " << converged_tests 
         << " (" << fixed << setprecision(2) << (100.0 * converged_tests / total_tests) << "%)\n";
    cout << "  Failed tests:          " << failed_tests 
         << " (" << setprecision(2) << (100.0 * failed_tests / total_tests) << "%)\n";
    cout << "    - NaN results:       " << nan_tests << "\n";
    cout << "    - Invalid region:    " << invalid_region_tests << "\n";
    cout << "    - Other failures:    " << (failed_tests - nan_tests - invalid_region_tests) << "\n";
    cout << "\n";
    
    // Additional diagnostics by salinity
    cout << "Convergence by Salinity:\n";
    cout << string(60, '-') << "\n";
    cout << setw(10) << "X"
         << setw(15) << "Tests"
         << setw(15) << "Success"
         << setw(15) << "Failed"
         << setw(15) << "Success %"
         << "\n";
    cout << string(60, '-') << "\n";
    
    for(size_t iX = 0; iX < X_values.size(); iX++) {
        double X_wt = X_values[iX];
        int X_total = 0;
        int X_success = 0;
        int X_failed = 0;
        
        for(size_t iP = 0; iP < P_values_bar.size(); iP++) {
            double P_bar = P_values_bar[iP];
            double P_Pa = P_bar * 1e5;
            
            for(size_t iH = 0; iH < H_values_kJ.size(); iH++) {
                double H_kJ = H_values_kJ[iH];
                double H_J = H_kJ * 1e3;
                
                X_total++;
                
                try {
                    H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX(P_Pa, H_J, X_wt);
                    
                    if(!isnan(prop.Rho) && !isnan(prop.T) && prop.Rho > 0.0 && 
                       prop.T >= -50.0 && prop.T <= 1100.0 && 
                       prop.Region >= 0 && prop.Region <= 10) {
                        X_success++;
                    } else {
                        X_failed++;
                    }
                } catch(...) {
                    X_failed++;
                }
            }
        }
        
        double success_pct = 100.0 * X_success / X_total;
        cout << setw(10) << fixed << setprecision(4) << X_wt
             << setw(15) << X_total
             << setw(15) << X_success
             << setw(15) << X_failed
             << setw(15) << setprecision(2) << success_pct
             << "\n";
    }
    
    cout << string(60, '=') << "\n\n";
    
    // Overall result
    if(failed_tests == 0) {
        cout << "✓✓✓ SUCCESS ✓✓✓\n";
        cout << "All tests passed! Driesner correlations are valid across tested P-H-X space.\n";
        return 0;
    } else if(failed_tests < total_tests * 0.05) {
        cout << "⚠ MOSTLY PASSED ⚠\n";
        cout << "Most tests passed (" << setprecision(2) << (100.0 * converged_tests / total_tests) 
             << "%), but some failures detected.\n";
        cout << "Review failed tests to determine if they are in physically valid regions.\n";
        return 1;
    } else {
        cout << "✗✗✗ FAILED ✗✗✗\n";
        cout << "Significant failures detected (" << setprecision(2) << (100.0 * failed_tests / total_tests) << "%).\n";
        cout << "Driesner correlations may have issues in the tested P-H-X space.\n";
        return 1;
    }
}
