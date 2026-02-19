#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "==========================================================\n";
    cout << "Direct Test of water_mu_pT at Problematic Conditions\n";
    cout << "==========================================================\n\n";
    
    // Test water_mu_pT directly at the exact conditions where jump occurs
    
    cout << "Testing pure water viscosity calculation:\n\n";
    cout << setw(12) << "P (bar)"
         << setw(12) << "T (°C)"
         << setw(15) << "T (K)"
         << setw(18) << "Mu from pT"
         << setw(15) << "T_sat (K)" << "\n";
    cout << string(80, '=') << "\n";
    
    // The problematic range is around P=99 bar, T=310°C
    for(double P_bar = 98.0; P_bar <= 100.0; P_bar += 0.2) {
        double P_Pa = P_bar * 1e5;
        double T_C = 309.0 + (P_bar - 98.0) * 0.5;  // Approximate T from our earlier tests
        double T_K = T_C + 273.15;
        
        // Use public method mu_pTX with X=0 for pure water
        double mu = eos.mu_pTX(P_Pa, T_K, 0.0);
        
        // Also need to know what T_sat is at this pressure
        // We can get this by checking pure water saturation
        H2ONaCl::PROP_H2ONaCl sat_prop = eos.prop_pTX(P_Pa, T_K - 5.0, 0.0, false);
        
        cout << setw(12) << fixed << setprecision(1) << P_bar
             << setw(12) << setprecision(2) << T_C
             << setw(15) << setprecision(3) << T_K
             << setw(18) << scientific << setprecision(6) << mu
             << setw(15) << fixed << setprecision(3) << (sat_prop.T + 273.15) << "\n";
    }
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "\nDetailed analysis at P = 98.8 bar vs P = 99.0 bar:\n\n";
    
    double P1 = 98.8e5;
    double P2 = 99.0e5;
    double T1 = 310.144 + 273.15;  // From our earlier diagnostic
    double T2 = 310.293 + 273.15;
    
    // Use public method mu_pTX with X=0 for pure water
    double mu1 = eos.mu_pTX(P1, T1, 0.0);
    double mu2 = eos.mu_pTX(P2, T2, 0.0);
    
    cout << "P = 98.8 bar, T = " << fixed << setprecision(3) << (T1 - 273.15) << " °C:\n";
    cout << "  mu_pTX(X=0) returned: " << scientific << setprecision(6) << mu1 << " Pa·s\n\n";
    
    cout << "P = 99.0 bar, T = " << fixed << setprecision(3) << (T2 - 273.15) << " °C:\n";
    cout << "  mu_pTX(X=0) returned: " << scientific << setprecision(6) << mu2 << " Pa·s\n\n";
    
    cout << "Ratio: " << fixed << setprecision(2) << (mu2 / mu1) << "x\n";
    cout << "Change: " << setprecision(1) << ((mu2/mu1 - 1.0) * 100.0) << "%\n\n";
    
    if(mu2 / mu1 > 2.0) {
        cout << "⚠ CONFIRMED: Pure water viscosity has a discontinuity!\n";
        cout << "This is the root cause of the vapor viscosity jump.\n\n";
        cout << "The issue is in the viscosity calculation for pure water.\n";
        cout << "At P≈99 bar, T≈310°C, the function is likely:\n";
        cout << "1. Near saturation line and triggering fallback logic\n";
        cout << "2. The fallback is selecting wrong phase (liquid vs vapor)\n";
        cout << "3. Or the dt tolerance (0.5 K) is still too small for this condition\n";
    }
    
    return 0;
}
