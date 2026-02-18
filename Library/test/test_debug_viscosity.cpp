#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "==========================================================\n";
    cout << "Debugging: What happens at P=98.8 vs P=99.0 bar?\n";
    cout << "==========================================================\n\n";
    
    const double X_wt = 0.0001;
    const double H_target = 1940.0e3;  // J/kg
    
    double P1_Pa = 98.8e5;
    double P2_Pa = 99.0e5;
    
    cout << "Getting properties with prop_pHX_bisection:\n\n";
    
    H2ONaCl::PROP_H2ONaCl prop1 = eos.prop_pHX_bisection(P1_Pa, H_target, X_wt);
    H2ONaCl::PROP_H2ONaCl prop2 = eos.prop_pHX_bisection(P2_Pa, H_target, X_wt);
    
    cout << "P = 98.8 bar:\n";
    cout << "  Region: " << prop1.Region << "\n";
    cout << "  T = " << fixed << setprecision(4) << prop1.T << " °C\n";
    cout << "  T_K = " << setprecision(4) << (prop1.T + 273.15) << " K\n";
    cout << "  Mu_v = " << scientific << setprecision(6) << prop1.Mu_v << " Pa·s\n";
    cout << "  Mu_l = " << setprecision(6) << prop1.Mu_l << " Pa·s\n";
    cout << "  X_l = " << fixed << setprecision(6) << prop1.X_l << "\n";
    cout << "  X_v = " << setprecision(6) << prop1.X_v << "\n\n";
    
    cout << "P = 99.0 bar:\n";
    cout << "  Region: " << prop2.Region << "\n";
    cout << "  T = " << fixed << setprecision(4) << prop2.T << " °C\n";
    cout << "  T_K = " << setprecision(4) << (prop2.T + 273.15) << " K\n";
    cout << "  Mu_v = " << scientific << setprecision(6) << prop2.Mu_v << " Pa·s\n";
    cout << "  Mu_l = " << setprecision(6) << prop2.Mu_l << " Pa·s\n";
    cout << "  X_l = " << fixed << setprecision(6) << prop2.X_l << "\n";
    cout << "  X_v = " << setprecision(6) << prop2.X_v << "\n\n";
    
    // Now try calling prop_pTX directly with these temperatures
    cout << "Testing prop_pTX directly at the same conditions:\n\n";
    
    double T1_K = prop1.T + 273.15;
    double T2_K = prop2.T + 273.15;
    
    H2ONaCl::PROP_H2ONaCl direct1 = eos.prop_pTX(P1_Pa, T1_K, X_wt, true);
    H2ONaCl::PROP_H2ONaCl direct2 = eos.prop_pTX(P2_Pa, T2_K, X_wt, true);
    
    cout << "prop_pTX at P=98.8 bar, T=" << fixed << setprecision(2) << prop1.T << "°C:\n";
    cout << "  Mu_v = " << scientific << setprecision(6) << direct1.Mu_v << " Pa·s\n";
    cout << "  Mu_l = " << setprecision(6) << direct1.Mu_l << " Pa·s\n\n";
    
    cout << "prop_pTX at P=99.0 bar, T=" << fixed << setprecision(2) << prop2.T << "°C:\n";
    cout << "  Mu_v = " << scientific << setprecision(6) << direct2.Mu_v << " Pa·s\n";
    cout << "  Mu_l = " << setprecision(6) << direct2.Mu_l << " Pa·s\n\n";
    
    // Check if it's an issue with calcViscosity_ph being called
    cout << "Ratio Mu_v (99.0/98.8): " << fixed << setprecision(2) 
         << (prop2.Mu_v / prop1.Mu_v) << "x\n";
    
    cout << "\nThe viscosity is calculated differently in prop_pHX_bisection!\n";
    cout << "Let me check if calcViscosity_ph is being used instead of calcViscosity...\n";
    
    return 0;
}
