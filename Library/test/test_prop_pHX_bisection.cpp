#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Testing prop_pHX_bisection implementation\n";
    cout << "==========================================\n\n";
    
    // Test case 1: Above critical point (pure water)
    cout << "Test 1: Above critical point (pure water)\n";
    double p1 = 25e6;  // 250 bar (above critical)
    double T1 = 400 + 273.15; // 400°C
    double X1 = 0.0;
    
    // Get properties at pTX to get the enthalpy
    H2ONaCl::PROP_H2ONaCl prop_pTX1 = eos.prop_pTX(p1, T1, X1);
    cout << "  pTX: P=" << p1/1e5 << " bar, T=" << T1-273.15 << " C, X=" << X1 << endl;
    cout << "    Region: " << prop_pTX1.Region << ", Rho=" << prop_pTX1.Rho << " kg/m3, H=" << prop_pTX1.H/1e3 << " kJ/kg\n";
    
    // Now call prop_pHX_bisection with the same P, H, X
    H2ONaCl::PROP_H2ONaCl prop_pHX1 = eos.prop_pHX_bisection(p1, prop_pTX1.H, X1);
    cout << "  pHX: P=" << p1/1e5 << " bar, H=" << prop_pTX1.H/1e3 << " kJ/kg, X=" << X1 << endl;
    cout << "    Region: " << prop_pHX1.Region << ", Rho=" << prop_pHX1.Rho << " kg/m3, T=" << prop_pHX1.T << " C\n";
    
    double rho_err1 = fabs(prop_pHX1.Rho - prop_pTX1.Rho);
    double T_err1 = fabs(prop_pHX1.T - (T1-273.15));
    cout << "  Errors: Rho=" << rho_err1 << " kg/m3, T=" << T_err1 << " C\n";
    
    if(isnan(prop_pHX1.Rho)) {
        cout << "  FAILED: Density is NaN!\n\n";
        return 1;
    } else if(rho_err1 > 1.0) {
        cout << "  WARNING: Large density error!\n\n";
    } else {
        cout << "  PASSED\n\n";
    }
    
    // Test case 2: Above critical point with small salinity
    cout << "Test 2: Above critical point with small salinity\n";
    double p2 = 25e6;  // 250 bar
    double T2 = 400 + 273.15; // 400°C (same as Test 1)
    double X2 = 0.01;  // 1% NaCl (small salinity)
    
    H2ONaCl::PROP_H2ONaCl prop_pTX2 = eos.prop_pTX(p2, T2, X2);
    cout << "  pTX: P=" << p2/1e5 << " bar, T=" << T2-273.15 << " C, X=" << X2 << endl;
    cout << "    Region: " << prop_pTX2.Region << ", Rho=" << prop_pTX2.Rho << " kg/m3, H=" << prop_pTX2.H/1e3 << " kJ/kg\n";
    
    H2ONaCl::PROP_H2ONaCl prop_pHX2 = eos.prop_pHX_bisection(p2, prop_pTX2.H, X2);
    cout << "  pHX: P=" << p2/1e5 << " bar, H=" << prop_pTX2.H/1e3 << " kJ/kg, X=" << X2 << endl;
    cout << "    Region: " << prop_pHX2.Region << ", Rho=" << prop_pHX2.Rho << " kg/m3, T=" << prop_pHX2.T << " C\n";
    
    double rho_err2 = fabs(prop_pHX2.Rho - prop_pTX2.Rho);
    double T_err2 = fabs(prop_pHX2.T - (T2-273.15));
    cout << "  Errors: Rho=" << rho_err2 << " kg/m3, T=" << T_err2 << " C\n";
    
    if(isnan(prop_pHX2.Rho)) {
        cout << "  FAILED: Density is NaN!\n\n";
        return 1;
    } else if(rho_err2 > 1.0) {
        cout << "  WARNING: Large density error!\n\n";
    } else {
        cout << "  PASSED\n\n";
    }
    
    // Test case 3: Two-phase region (if applicable)
    cout << "Test 3: Two-phase region\n";
    double p3 = 10e6;  // 100 bar
    double T3 = 300 + 273.15; // 300°C
    double X3 = 0.05;  // 5% NaCl
    
    H2ONaCl::PROP_H2ONaCl prop_pTX3 = eos.prop_pTX(p3, T3, X3);
    cout << "  pTX: P=" << p3/1e5 << " bar, T=" << T3-273.15 << " C, X=" << X3 << endl;
    cout << "    Region: " << prop_pTX3.Region << ", Rho=" << prop_pTX3.Rho << " kg/m3, H=" << prop_pTX3.H/1e3 << " kJ/kg\n";
    
    H2ONaCl::PROP_H2ONaCl prop_pHX3 = eos.prop_pHX_bisection(p3, prop_pTX3.H, X3);
    cout << "  pHX: P=" << p3/1e5 << " bar, H=" << prop_pTX3.H/1e3 << " kJ/kg, X=" << X3 << endl;
    cout << "    Region: " << prop_pHX3.Region << ", Rho=" << prop_pHX3.Rho << " kg/m3, T=" << prop_pHX3.T << " C\n";
    
    double rho_err3 = fabs(prop_pHX3.Rho - prop_pTX3.Rho);
    double T_err3 = fabs(prop_pHX3.T - (T3-273.15));
    cout << "  Errors: Rho=" << rho_err3 << " kg/m3, T=" << T_err3 << " C\n";
    
    if(isnan(prop_pHX3.Rho)) {
        cout << "  FAILED: Density is NaN!\n\n";
        return 1;
    } else if(rho_err3 > 1.0) {
        cout << "  WARNING: Large density error!\n\n";
    } else {
        cout << "  PASSED\n\n";
    }
    
    cout << "==========================================\n";
    cout << "All tests completed!\n";
    
    return 0;
}
