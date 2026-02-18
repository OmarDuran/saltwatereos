#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Testing for negative liquid saturation at high enthalpy\n";
    cout << "=======================================================\n\n";
    
    cout << "Expected: At high enthalpy (superheated vapor), S_l should be 0 or very small positive\n";
    cout << "Problem: Getting negative S_l values\n\n";
    
    int failed_count = 0;
    int total_tests = 0;
    
    // Test 1: High enthalpy, low salinity, moderate pressure
    cout << "Test 1: High enthalpy, low salinity, moderate pressure\n";
    double p1 = 100e5;  // 100 bar
    double T1 = 600.0 + 273.15; // 600°C - superheated
    double X1 = 0.01;  // 1% NaCl
    
    H2ONaCl::PROP_H2ONaCl prop1 = eos.prop_pTX(p1, T1, X1, false);
    cout << "  pTX: P=" << p1/1e5 << " bar, T=" << T1-273.15 << " C, X=" << X1 << "\n";
    cout << "  Result: Region=" << prop1.Region << ", S_l=" << prop1.S_l << ", S_v=" << prop1.S_v << "\n";
    cout << "  H=" << prop1.H/1e3 << " kJ/kg, Rho=" << prop1.Rho << " kg/m3\n";
    
    total_tests++;
    if(prop1.S_l < -1e-10) {
        cout << "  FAILED: S_l is negative (" << prop1.S_l << ")\n\n";
        failed_count++;
    } else if(prop1.S_l > 0.5) {
        cout << "  WARNING: S_l unexpectedly high for superheated vapor\n\n";
    } else {
        cout << "  PASSED\n\n";
    }
    
    // Test 2: Use prop_pHX_bisection with high enthalpy
    cout << "Test 2: prop_pHX_bisection with high enthalpy from Test 1\n";
    double H_high = prop1.H * 1.2; // 20% higher than calculated
    cout << "  Input: P=" << p1/1e5 << " bar, H=" << H_high/1e3 << " kJ/kg, X=" << X1 << "\n";
    
    H2ONaCl::PROP_H2ONaCl prop2 = eos.prop_pHX_bisection(p1, H_high, X1);
    cout << "  Result: Region=" << prop2.Region << ", S_l=" << prop2.S_l << ", S_v=" << prop2.S_v << "\n";
    cout << "  T=" << prop2.T << " C, Rho=" << prop2.Rho << " kg/m3\n";
    
    total_tests++;
    if(prop2.S_l < -1e-10) {
        cout << "  FAILED: S_l is negative (" << prop2.S_l << ")\n\n";
        failed_count++;
    } else if(prop2.S_l > 0.5) {
        cout << "  WARNING: S_l unexpectedly high for high enthalpy\n\n";
    } else {
        cout << "  PASSED\n\n";
    }
    
    // Test 3: Very high enthalpy, pure water
    cout << "Test 3: Very high enthalpy, pure water\n";
    double p3 = 50e5;  // 50 bar
    double H3 = 3500e3; // 3500 kJ/kg - very high
    double X3 = 0.0;    // Pure water
    
    cout << "  Input: P=" << p3/1e5 << " bar, H=" << H3/1e3 << " kJ/kg, X=" << X3 << "\n";
    
    H2ONaCl::PROP_H2ONaCl prop3 = eos.prop_pHX_bisection(p3, H3, X3);
    cout << "  Result: Region=" << prop3.Region << ", S_l=" << prop3.S_l << ", S_v=" << prop3.S_v << "\n";
    cout << "  T=" << prop3.T << " C, Rho=" << prop3.Rho << " kg/m3\n";
    
    total_tests++;
    if(prop3.S_l < -1e-10) {
        cout << "  FAILED: S_l is negative (" << prop3.S_l << ")\n\n";
        failed_count++;
    } else if(prop3.S_l > 0.5) {
        cout << "  WARNING: S_l unexpectedly high for superheated vapor\n\n";
    } else {
        cout << "  PASSED\n\n";
    }
    
    // Test 4: Range of high enthalpies with low salinity
    cout << "Test 4: Range of high enthalpies, low salinity\n";
    cout << "  P=100 bar, X=0.01\n";
    cout << "  H(kJ/kg)\tS_l\tS_v\tT(C)\tRegion\n";
    cout << "  --------\t---\t---\t----\t------\n";
    
    for(double H = 2000e3; H <= 4000e3; H += 500e3) {
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(100e5, H, 0.01);
        cout << "  " << fixed << setprecision(0) << H/1e3 
             << "\t\t" << setprecision(4) << prop.S_l 
             << "\t" << prop.S_v 
             << "\t" << setprecision(1) << prop.T
             << "\t" << prop.Region << "\n";
        
        total_tests++;
        if(prop.S_l < -1e-10) {
            cout << "       ^^^ FAILED: Negative S_l!\n";
            failed_count++;
        }
    }
    cout << "\n";
    
    // Test 5: Check two-phase region behavior
    cout << "Test 5: Two-phase region (should have 0 < S_l < 1)\n";
    double p5 = 100e5;
    double T5 = 300.0 + 273.15; // 300°C
    double X5 = 0.05;
    
    H2ONaCl::PROP_H2ONaCl prop5 = eos.prop_pTX(p5, T5, X5, false);
    cout << "  pTX: P=" << p5/1e5 << " bar, T=" << T5-273.15 << " C, X=" << X5 << "\n";
    cout << "  Result: Region=" << prop5.Region << ", S_l=" << prop5.S_l << ", S_v=" << prop5.S_v << ", S_h=" << prop5.S_h << "\n";
    
    total_tests++;
    if(prop5.Region >= 5 && prop5.Region <= 8) { // Two-phase regions
        if(prop5.S_l < -1e-10 || prop5.S_l > 1.0001) {
            cout << "  FAILED: S_l out of range [0,1] for two-phase region\n\n";
            failed_count++;
        } else {
            cout << "  PASSED\n\n";
        }
    } else {
        cout << "  Note: Not a two-phase region, skipping range check\n\n";
    }
    
    // Summary
    cout << "=======================================================\n";
    cout << "Summary: " << (total_tests - failed_count) << "/" << total_tests << " tests passed\n";
    
    if(failed_count > 0) {
        cout << "FAILED: Found " << failed_count << " cases with negative or invalid saturation!\n";
        return 1;
    } else {
        cout << "SUCCESS: All saturation values are valid!\n";
        return 0;
    }
}
