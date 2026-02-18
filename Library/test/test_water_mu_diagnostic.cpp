#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "==========================================================\n";
    cout << "Diagnostic: Detailed Properties at Viscosity Jump\n";
    cout << "==========================================================\n\n";
    
    const double X_wt = 0.0001;
    const double H_target = 1940.0e3;  // J/kg
    
    cout << "Examining detailed properties around the jump (98-100 bar):\n\n";
    cout << setw(10) << "P (bar)"
         << setw(12) << "T (°C)"
         << setw(15) << "Mu_v (Pa·s)"
         << setw(15) << "Mu_l (Pa·s)"
         << setw(12) << "Region"
         << setw(10) << "S_l"
         << setw(10) << "S_v"
         << setw(12) << "Rho_v"
         << setw(12) << "Rho_l" << "\n";
    cout << string(120, '=') << "\n";
    
    // Test around the jump with very fine resolution
    double prev_Mu_v = 0.0;
    for(double P_bar = 98.0; P_bar <= 100.0; P_bar += 0.1) {
        double P_Pa = P_bar * 1e5;
        
        // Get full properties
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(P_Pa, H_target, X_wt);
        
        double delta_Mu_v = 0.0;
        if(prev_Mu_v > 1e-12 && prop.Mu_v > 1e-12) {
            delta_Mu_v = 100.0 * (prop.Mu_v - prev_Mu_v) / prev_Mu_v;
        }
        
        cout << setw(10) << fixed << setprecision(1) << P_bar
             << setw(12) << setprecision(3) << prop.T
             << setw(15) << scientific << setprecision(4) << prop.Mu_v
             << setw(15) << setprecision(4) << prop.Mu_l
             << setw(12) << prop.Region
             << setw(10) << fixed << setprecision(4) << prop.S_l
             << setw(10) << setprecision(4) << prop.S_v
             << setw(12) << setprecision(2) << prop.Rho_v
             << setw(12) << setprecision(2) << prop.Rho_l;
        
        if(std::abs(delta_Mu_v) > 50.0) {
            cout << "  *** JUMP: " << fixed << setprecision(1) << delta_Mu_v << "% ***";
        }
        cout << "\n";
        
        prev_Mu_v = prop.Mu_v;
    }
    
    cout << "\n" << string(120, '=') << "\n\n";
    
    // Specific comparison at the exact jump
    cout << "Detailed comparison at the jump:\n\n";
    
    double P1 = 98.8e5;
    double P2 = 98.9e5;
    
    H2ONaCl::PROP_H2ONaCl prop1 = eos.prop_pHX_bisection(P1, H_target, X_wt);
    H2ONaCl::PROP_H2ONaCl prop2 = eos.prop_pHX_bisection(P2, H_target, X_wt);
    
    cout << "At P = 98.8 bar:\n";
    cout << "  T = " << fixed << setprecision(3) << prop1.T << " °C\n";
    cout << "  Mu_v = " << scientific << setprecision(6) << prop1.Mu_v << " Pa·s\n";
    cout << "  Mu_l = " << setprecision(6) << prop1.Mu_l << " Pa·s\n";
    cout << "  Region: " << prop1.Region << "\n";
    cout << "  S_l = " << fixed << setprecision(4) << prop1.S_l << ", S_v = " << prop1.S_v << "\n";
    cout << "  Rho_v = " << setprecision(2) << prop1.Rho_v << " kg/m³, Rho_l = " << prop1.Rho_l << " kg/m³\n\n";
    
    cout << "At P = 98.9 bar:\n";
    cout << "  T = " << fixed << setprecision(3) << prop2.T << " °C\n";
    cout << "  Mu_v = " << scientific << setprecision(6) << prop2.Mu_v << " Pa·s\n";
    cout << "  Mu_l = " << setprecision(6) << prop2.Mu_l << " Pa·s\n";
    cout << "  Region: " << prop2.Region << "\n";
    cout << "  S_l = " << fixed << setprecision(4) << prop2.S_l << ", S_v = " << prop2.S_v << "\n";
    cout << "  Rho_v = " << setprecision(2) << prop2.Rho_v << " kg/m³, Rho_l = " << prop2.Rho_l << " kg/m³\n\n";
    
    double ratio = prop2.Mu_v / prop1.Mu_v;
    cout << "Mu_v ratio (98.9/98.8): " << fixed << setprecision(2) << ratio << "x\n";
    cout << "Change: " << setprecision(1) << ((ratio - 1.0) * 100.0) << "%\n\n";
    
    cout << "Analysis:\n";
    cout << "  Temperature change: " << setprecision(3) << (prop2.T - prop1.T) << " °C (smooth)\n";
    cout << "  Liquid saturation change: " << setprecision(4) << (prop2.S_l - prop1.S_l) << " (smooth)\n";
    cout << "  The jump appears to be in the vapor viscosity calculation method itself.\n";
    
    return 0;
}
