#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "==========================================================\n";
    cout << "Testing Vapor Viscosity at Supercritical Conditions\n";
    cout << "==========================================================\n\n";
    
    const double X_wt = 0.0001;
    
    cout << "Testing at X = " << X_wt << "\n";
    cout << "Scanning high pressure (P > 221 bar) and high enthalpy (H > 2950 kJ/kg)\n\n";
    
    cout << setw(10) << "P (bar)"
         << setw(12) << "H (kJ/kg)"
         << setw(12) << "T (°C)"
         << setw(12) << "Region"
         << setw(15) << "Mu_v (Pa·s)"
         << setw(15) << "Mu_l (Pa·s)"
         << setw(10) << "S_v" << "\n";
    cout << string(95, '=') << "\n";
    
    // Test at pressures above critical (221 bar)
    vector<double> pressures = {221, 225, 230, 240, 250, 260, 280, 300};
    vector<double> enthalpies = {2900, 2950, 3000, 3100, 3200, 3300, 3400, 3500};
    
    int zero_count = 0;
    int hardcoded_count = 0;
    
    for(auto P_bar : pressures) {
        for(auto H_kJ : enthalpies) {
            double P_Pa = P_bar * 1e5;
            double H_J = H_kJ * 1e3;
            
            H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(P_Pa, H_J, X_wt);
            
            cout << setw(10) << fixed << setprecision(0) << P_bar
                 << setw(12) << setprecision(0) << H_kJ
                 << setw(12) << setprecision(1) << prop.T
                 << setw(12) << prop.Region
                 << setw(15) << scientific << setprecision(4) << prop.Mu_v
                 << setw(15) << setprecision(4) << prop.Mu_l
                 << setw(10) << fixed << setprecision(4) << prop.S_v;
            
            if(prop.Mu_v == 0.0) {
                cout << "  *** ZERO ***";
                zero_count++;
            } else if(abs(prop.Mu_v - 1.0e-4) < 1e-10) {
                cout << "  *** HARDCODED 1e-4 ***";
                hardcoded_count++;
            }
            cout << "\n";
        }
    }
    
    cout << "\n" << string(95, '=') << "\n";
    cout << "Summary:\n";
    cout << "  Tests with Mu_v = 0: " << zero_count << "\n";
    cout << "  Tests with Mu_v = 1e-4 (hardcoded): " << hardcoded_count << "\n";
    
    if(zero_count > 0 || hardcoded_count > 0) {
        cout << "\n⚠ ISSUE DETECTED: Vapor viscosity has incorrect values!\n";
        cout << "This is likely due to fallback logic returning wrong defaults.\n";
        return 1;
    }
    
    cout << "\n✓ All vapor viscosity values are reasonable.\n";
    return 0;
}
