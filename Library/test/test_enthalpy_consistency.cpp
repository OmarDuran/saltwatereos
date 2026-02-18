#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

struct TestCase {
    double P_bar;
    double T_C;
    double X_wt;
    string description;
};

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Testing Enthalpy Consistency in prop_pHX_bisection\n";
    cout << "===================================================\n\n";
    
    cout << "Objective: Ensure specified H matches computed H within tolerance\n";
    cout << "Tolerance: 1e-4 relative error (0.01%)\n\n";
    
    // Define test cases covering different regions and conditions
    vector<TestCase> testCases = {
        {100.0, 300.0, 0.05, "Single phase liquid, moderate salinity"},
        {100.0, 600.0, 0.01, "High temperature, low salinity"},
        {250.0, 400.0, 0.0,  "Pure water above critical P"},
        {250.0, 400.0, 0.01, "Saline water above critical P"},
        {50.0,  200.0, 0.1,  "Low pressure, moderate salinity"},
        {200.0, 500.0, 0.05, "High pressure, high temperature"},
        {10.0,  100.0, 0.0,  "Low pressure pure water"},
        {100.0, 350.0, 0.15, "Moderate P/T, high salinity"},
        {150.0, 450.0, 0.02, "Near critical conditions"},
        {300.0, 700.0, 0.05, "Very high P/T"}
    };
    
    int total_tests = 0;
    int failed_tests = 0;
    double max_error = 0.0;
    double sum_error = 0.0;
    
    const double rel_tol = 1e-4; // 0.01% relative error
    
    cout << "Test Results:\n";
    cout << "-------------\n\n";
    
    for(size_t i = 0; i < testCases.size(); i++) {
        TestCase tc = testCases[i];
        
        cout << "Test " << (i+1) << ": " << tc.description << "\n";
        cout << "  Input: P=" << tc.P_bar << " bar, T=" << tc.T_C << " C, X=" << tc.X_wt << "\n";
        
        // Step 1: Get properties at (P, T, X) to obtain H
        H2ONaCl::PROP_H2ONaCl prop_pTX = eos.prop_pTX(tc.P_bar * 1e5, tc.T_C + 273.15, tc.X_wt, false);
        
        double H_target = prop_pTX.H;
        
        if(isnan(H_target) || H_target == 0.0) {
            cout << "  SKIPPED: Invalid H from prop_pTX\n\n";
            continue;
        }
        
        cout << "  H_target from pTX: " << fixed << setprecision(2) << H_target/1e3 << " kJ/kg\n";
        
        // Step 2: Use prop_pHX_bisection with the target H
        H2ONaCl::PROP_H2ONaCl prop_pHX = eos.prop_pHX_bisection(tc.P_bar * 1e5, H_target, tc.X_wt);
        
        double H_computed = prop_pHX.H;
        
        cout << "  H_computed from pHX: " << setprecision(2) << H_computed/1e3 << " kJ/kg\n";
        
        // Step 3: Check consistency
        double abs_error = fabs(H_computed - H_target);
        double rel_error = abs_error / fabs(H_target);
        
        cout << "  Absolute error: " << scientific << setprecision(3) << abs_error << " J/kg\n";
        cout << "  Relative error: " << setprecision(3) << rel_error * 100 << " %\n";
        
        total_tests++;
        
        if(rel_error > rel_tol) {
            cout << "  FAILED: Relative error exceeds tolerance!\n";
            failed_tests++;
        } else {
            cout << "  PASSED\n";
        }
        
        max_error = max(max_error, rel_error);
        sum_error += rel_error;
        
        cout << "\n";
    }
    
    // Summary statistics
    cout << "===================================================\n";
    cout << "Summary:\n";
    cout << "  Total tests: " << total_tests << "\n";
    cout << "  Passed: " << (total_tests - failed_tests) << "\n";
    cout << "  Failed: " << failed_tests << "\n";
    cout << "  Maximum relative error: " << scientific << setprecision(3) << max_error * 100 << " %\n";
    cout << "  Average relative error: " << setprecision(3) << (sum_error / total_tests) * 100 << " %\n";
    cout << "===================================================\n\n";
    
    if(failed_tests > 0) {
        cout << "FAILED: " << failed_tests << " test(s) exceeded tolerance!\n";
        return 1;
    } else {
        cout << "SUCCESS: All enthalpy values are consistent within tolerance!\n";
        return 0;
    }
}
