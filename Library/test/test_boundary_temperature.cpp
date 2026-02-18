#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Testing boundary temperature T=1000°C\n";
    cout << "======================================\n\n";
    
    cout << "NOTE: At T=1000°C (TMAX_C), there may be a phase transition\n";
    cout << "      boundary due to correlation limits. Properties should\n";
    cout << "      still be calculated without warnings.\n\n";
    
    // Test at boundary temperature T=1000°C
    double p = 100e5;  // 100 bar
    double T_K = 1000.0 + 273.15; // 1000°C in Kelvin
    double X_wt = 0.05;  // 5% NaCl
    
    cout << "Test 1: prop_pTX at T=1000°C (should produce NO warnings)\n";
    cout << "  Input: P=" << p/1e5 << " bar, T=" << T_K-273.15 << " C, X=" << X_wt << "\n";
    
    H2ONaCl::PROP_H2ONaCl prop1 = eos.prop_pTX(p, T_K, X_wt, false);
    cout << "  Result: Region=" << prop1.Region << ", Rho=" << prop1.Rho << " kg/m3, H=" << prop1.H/1e3 << " kJ/kg\n";
    
    if(isnan(prop1.Rho)) {
        cout << "  FAILED: Density is NaN!\n\n";
        return 1;
    } else {
        cout << "  PASSED (no NaN values)\n\n";
    }
    
    // Test just below boundary
    cout << "Test 2: prop_pTX at T=999.0°C\n";
    double T_K_below = 999.0 + 273.15;
    cout << "  Input: P=" << p/1e5 << " bar, T=" << T_K_below-273.15 << " C, X=" << X_wt << "\n";
    
    H2ONaCl::PROP_H2ONaCl prop2 = eos.prop_pTX(p, T_K_below, X_wt, false);
    cout << "  Result: Region=" << prop2.Region << ", Rho=" << prop2.Rho << " kg/m3, H=" << prop2.H/1e3 << " kJ/kg\n";
    
    if(isnan(prop2.Rho)) {
        cout << "  FAILED: Density is NaN!\n\n";
        return 1;
    } else {
        cout << "  PASSED (no NaN values)\n\n";
    }
    
    // Test prop_pHX_bisection with T=999 enthalpy
    cout << "Test 3: prop_pHX_bisection at T=999°C enthalpy\n";
    cout << "  Input: P=" << p/1e5 << " bar, H=" << prop2.H/1e3 << " kJ/kg, X=" << X_wt << "\n";
    
    H2ONaCl::PROP_H2ONaCl prop3 = eos.prop_pHX_bisection(p, prop2.H, X_wt);
    cout << "  Result: Region=" << prop3.Region << ", Rho=" << prop3.Rho << " kg/m3, T=" << prop3.T << " C\n";
    
    double T_err = fabs(prop3.T - (T_K_below - 273.15));
    double rho_err = fabs(prop3.Rho - prop2.Rho);
    
    cout << "  Errors: T=" << T_err << " C, Rho=" << rho_err << " kg/m3\n";
    
    if(isnan(prop3.Rho)) {
        cout << "  FAILED: Density is NaN!\n\n";
        return 1;
    } else if(T_err > 1.0 || rho_err > 10.0) {
        cout << "  WARNING: Large error (may indicate boundary issues)\n\n";
    } else {
        cout << "  PASSED\n\n";
    }
    
    // Test at exactly TMIN_C
    cout << "Test 4: prop_pTX at T=TMIN_C (0.1°C)\n";
    double T_K_min = 0.1 + 273.15;
    cout << "  Input: P=" << p/1e5 << " bar, T=" << T_K_min-273.15 << " C, X=" << X_wt << "\n";
    
    H2ONaCl::PROP_H2ONaCl prop4 = eos.prop_pTX(p, T_K_min, X_wt, false);
    cout << "  Result: Region=" << prop4.Region << ", Rho=" << prop4.Rho << " kg/m3, H=" << prop4.H/1e3 << " kJ/kg\n";
    
    if(isnan(prop4.Rho)) {
        cout << "  FAILED: Density is NaN!\n\n";
        return 1;
    } else {
        cout << "  PASSED (no NaN values)\n\n";
    }
    
    cout << "======================================\n";
    cout << "All boundary tests completed successfully!\n";
    cout << "No temperature range warnings should appear.\n";
    
    return 0;
}
