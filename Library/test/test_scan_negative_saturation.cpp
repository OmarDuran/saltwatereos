#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Scanning for negative liquid saturation\n";
    cout << "========================================\n\n";
    
    int neg_count = 0;
    int total_count = 0;
    
    // Scan different conditions to find negative saturations
    cout << "Scanning P-H-X space for negative S_l:\n\n";
    
    double pressures[] = {10e5, 50e5, 100e5, 200e5, 300e5};
    double salinities[] = {0.0, 0.001, 0.01, 0.05, 0.1};
    
    for(int ip = 0; ip < 5; ip++) {
        double p = pressures[ip];
        
        for(int ix = 0; ix < 5; ix++) {
            double X = salinities[ix];
            
            // Get enthalpy at high temperature
            H2ONaCl::PROP_H2ONaCl prop_high_T = eos.prop_pTX(p, 800.0 + 273.15, X, false);
            double H_base = prop_high_T.H;
            
            // Test range of enthalpies above the base
            for(double H_mult = 1.0; H_mult <= 1.5; H_mult += 0.1) {
                double H = H_base * H_mult;
                
                H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(p, H, X);
                
                total_count++;
                
                if(prop.S_l < -1e-10) {
                    cout << "NEGATIVE S_l found!\n";
                    cout << "  P=" << p/1e5 << " bar, H=" << H/1e3 << " kJ/kg, X=" << X << "\n";
                    cout << "  Region=" << prop.Region << ", T=" << prop.T << " C\n";
                    cout << "  S_l=" << prop.S_l << ", S_v=" << prop.S_v << ", S_h=" << prop.S_h << "\n";
                    cout << "  Rho=" << prop.Rho << " kg/m3\n";
                    cout << "  Rho_l=" << prop.Rho_l << ", Rho_v=" << prop.Rho_v << ", Rho_h=" << prop.Rho_h << "\n";
                    cout << "  X_l=" << prop.X_l << ", X_v=" << prop.X_v << "\n\n";
                    neg_count++;
                }
            }
        }
    }
    
    cout << "========================================\n";
    cout << "Scanned " << total_count << " conditions\n";
    cout << "Found " << neg_count << " cases with negative S_l\n\n";
    
    if(neg_count > 0) {
        cout << "ISSUE CONFIRMED: Negative liquid saturations detected!\n";
        return 1;
    } else {
        cout << "No negative saturations found in this scan.\n";
        cout << "Please provide specific P, H, X values where you see the issue.\n";
        return 0;
    }
}
