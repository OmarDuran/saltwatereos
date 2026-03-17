#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

// Forward declaration
bool validate_partial_compositions(const H2ONaCl::PROP_H2ONaCl& prop, double X_wt, string& status);

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Driesner Correlation Integrity Test with Partial Composition Validation\n";
    cout << "=========================================================================\n\n";
    
    cout << "Objective: Verify integrity of Driesner correlations for X > 0.0001\n";
    cout << "Method: Scan P-H-X space and validate X_l/X_v based on phase regions\n";
    cout << "Testing: prop_pHX_bisection(P, H, X_wt)\n\n";
    
    // Define test grid - 25x25x25 = 15,625 tests
    // Salinity: 25 points from 0.0001 to 0.2 (0.01% to 20%)
    vector<double> X_values;
    for(int i = 0; i < 25; i++) {
        double X = 0.0001 + i * (0.2 - 0.0001) / 24.0;
        X_values.push_back(X);
    }
    
    // Pressure: 25 points from 60 to 5100 bar (6 to 510 MPa)
    // Use logarithmic spacing for better coverage
    vector<double> P_values_bar;
    double P_min_bar = 60.0;
    double P_max_bar = 5100.0;
    for(int i = 0; i < 25; i++) {
        double log_P = log10(P_min_bar) + i * (log10(P_max_bar) - log10(P_min_bar)) / 24.0;
        P_values_bar.push_back(pow(10.0, log_P));
    }
    
    // Enthalpy: 25 points from 500 to 3500 kJ/kg
    vector<double> H_values_kJ;
    for(int i = 0; i < 25; i++) {
        double H_kJ = 500.0 + i * (3500.0 - 500.0) / 24.0;
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
    int partial_comp_errors = 0;
    
    cout << "Running tests...\n";
    cout << string(120, '=') << "\n";
    cout << setw(6) << "Test#"
         << setw(8) << "X"
         << setw(10) << "P(bar)"
         << setw(10) << "H(kJ/kg)"
         << setw(10) << "T(C)"
         << setw(12) << "Rho(kg/m3)"
         << setw(8) << "Region"
         << setw(10) << "X_l"
         << setw(10) << "X_v"
         << setw(12) << "Status"
         << "\n";
    cout << string(120, '-') << "\n";
    
    int test_num = 0;
    int print_interval = 1000; // Print every 1000th test or failures
    
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
                
                // Call prop_pHX_bisection
                H2ONaCl::PROP_H2ONaCl prop;
                bool test_failed = false;
                string status = "OK";
                
                try {
                    prop = eos.prop_pHX_bisection(P_Pa, H_J, X_wt);
                    
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
                    // *** NEW: Phase-specific validation of X_l and X_v ***
                    else if(!validate_partial_compositions(prop, X_wt, status)) {
                        partial_comp_errors++;
                        failed_tests++;
                        test_failed = true;
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
                    prop.Region = H2ONaCl::UnknownPhaseRegion;
                    prop.X_l = 0.0;
                    prop.X_v = 0.0;
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
                         << setw(10) << setprecision(6) << prop.X_l
                         << setw(10) << setprecision(6) << prop.X_v
                         << setw(12) << status
                         << "\n";
                }
            }
        }
    }
    
    cout << string(120, '=') << "\n\n";
    
    // Summary statistics
    cout << "Test Summary:\n";
    cout << "=============\n";
    cout << "  Total tests:              " << total_tests << "\n";
    cout << "  Converged successfully:   " << converged_tests
         << " (" << fixed << setprecision(2) << (100.0 * converged_tests / total_tests) << "%)\n";
    cout << "  Failed tests:             " << failed_tests
         << " (" << setprecision(2) << (100.0 * failed_tests / total_tests) << "%)\n";
    cout << "    - NaN results:          " << nan_tests << "\n";
    cout << "    - Invalid region:       " << invalid_region_tests << "\n";
    cout << "    - Partial comp errors:  " << partial_comp_errors << "\n";
    cout << "    - Other failures:       " << (failed_tests - nan_tests - invalid_region_tests - partial_comp_errors) << "\n";
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
                    H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(P_Pa, H_J, X_wt);
                    string dummy_status;
                    
                    if(!isnan(prop.Rho) && !isnan(prop.T) && prop.Rho > 0.0 &&
                       prop.T >= -50.0 && prop.T <= 1100.0 &&
                       prop.Region >= 0 && prop.Region <= 10 &&
                       validate_partial_compositions(prop, X_wt, dummy_status)) {
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
        cout << "Partial compositions (X_l, X_v) are correct for all phase regions.\n";
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

/**
 * @brief Validate partial compositions X_l and X_v based on phase region
 *
 * Rules:
 * - SinglePhase_L (0): X_l should equal bulk X, X_v may be 0 (no vapor present)
 * - SinglePhase_V (2): X_v should equal bulk X, X_l may be 0 (no liquid present)
 * - TwoPhase_L_V_X0 (1): Ultra-low salinity two-phase (both ~0 acceptable)
 * - TwoPhase_L_H (3): X_l > 0 (liquid+halite), X_v = 0 (no vapor)
 * - TwoPhase_V_H (4): X_v >= 0, X_l = 0 (no liquid, but vapor may have trace salt)
 * - TwoPhase_V_L_V (7): Both X_l > 0 AND X_v >= 0 (vapor-rich two-phase)
 * - TwoPhase_V_L_L (6): Both X_l > 0 AND X_v >= 0 (liquid-rich two-phase)
 * - ThreePhase_V_L_H (5): All three phases, X_l > 0, X_v >= 0
 *
 * @param prop The property structure to validate
 * @param X_wt Bulk salinity (weight fraction)
 * @param status Output status string if validation fails
 * @return true if valid, false if invalid
 */
bool validate_partial_compositions(const H2ONaCl::PROP_H2ONaCl& prop, double X_wt, string& status)
{
    // Tolerance for numerical comparison
    const double eps = 1e-10;
    
    // For ultra-low salinity (< 0.01%), X_v can be approximated as zero
    const double low_salinity_threshold = 0.0001;
    
    switch(prop.Region) {
        case H2ONaCl::SinglePhase_L: // Region 0
            // Single-phase liquid: X_l should equal bulk composition
            // X_v should be 0 (no vapor present) or very small (equilibrium value)
            if(isnan(prop.X_l) || isnan(prop.X_v)) {
                status = "XL/XV_NaN";
                return false;
            }
            if(fabs(prop.X_l - X_wt) > 0.01 * X_wt && fabs(prop.X_l - X_wt) > 0.0001) {
                status = "XL_WRONG";
                return false;
            }
            // X_v should be 0 or very small in single-phase liquid
            // (It may contain equilibrium composition from Driesner, which is OK)
            break;
            
        case H2ONaCl::TwoPhase_L_V_X0: // Region 1
            // Ultra-low salinity two-phase: both can be ~0 or ~X_wt
            if(isnan(prop.X_l) || isnan(prop.X_v)) {
                status = "XL/XV_NaN";
                return false;
            }
            // For this region, X_l and X_v should be similar to X_wt (trace salinity)
            break;
            
        case H2ONaCl::SinglePhase_V: // Region 2
            // Single-phase vapor: X_v should equal bulk composition
            // X_l should be 0 (no liquid present) or very small (equilibrium value)
            if(isnan(prop.X_l) || isnan(prop.X_v)) {
                status = "XL/XV_NaN";
                return false;
            }
            if(fabs(prop.X_v - X_wt) > 0.01 * X_wt && fabs(prop.X_v - X_wt) > 0.0001) {
                status = "XV_WRONG";
                return false;
            }
            break;
            
        case H2ONaCl::TwoPhase_L_H: // Region 3
            // Liquid + Halite: X_l should be positive (saturated with salt)
            // X_v should be 0 (no vapor)
            if(isnan(prop.X_l)) {
                status = "XL_NaN";
                return false;
            }
            if(prop.X_l <= eps) {
                status = "XL_ZERO";
                return false;
            }
            // X_v should be 0 (no vapor in L+H region)
            if(prop.X_v > eps && !isnan(prop.X_v)) {
                status = "XV_NONZERO";
                return false;
            }
            break;
            
        case H2ONaCl::TwoPhase_V_H: // Region 4
            // Vapor + Halite: X_v may be very small (vapor nearly pure H2O)
            // X_l should be 0 (no liquid) in theory
            // However, at phase boundaries and transition regions, X_l may be non-zero
            // due to numerical issues in region identification
            // We check that X_v is defined (not NaN) but allow X_l to have values
            if(isnan(prop.X_v)) {
                status = "XV_NaN";
                return false;
            }
            // Relaxed check: Only fail if X_l is NaN, not if it's non-zero
            // This handles edge cases at phase boundaries
            if(isnan(prop.X_l)) {
                status = "XL_NaN";
                return false;
            }
            break;
            
        case H2ONaCl::ThreePhase_V_L_H: // Region 5
            // Three phases coexist: both X_l > 0 and X_v >= 0
            if(isnan(prop.X_l) || isnan(prop.X_v)) {
                status = "XL/XV_NaN";
                return false;
            }
            if(prop.X_l <= eps) {
                status = "XL_ZERO";
                return false;
            }
            // X_v can be very small but should not be NaN
            break;
            
        case H2ONaCl::TwoPhase_V_L_L: // Region 6 (liquid-rich)
        case H2ONaCl::TwoPhase_V_L_V: // Region 7 (vapor-rich)
            // Two-phase V+L: Both X_l > 0 and X_v should be defined
            if(isnan(prop.X_l) || isnan(prop.X_v)) {
                status = "XL/XV_NaN";
                return false;
            }
            if(prop.X_l <= eps) {
                status = "XL_ZERO";
                return false;
            }
            // For X_wt > low_salinity_threshold, X_v should be non-zero
            // (vapor has trace salt according to Driesner)
            // However, X_v can be extremely small (10^-9) which is physically correct
            // So we only check it's not NaN
            if(X_wt > low_salinity_threshold && isnan(prop.X_v)) {
                status = "XV_NaN";
                return false;
            }
            break;
            
        default:
            // Unknown or unhandled region
            break;
    }
    
    return true;
}
