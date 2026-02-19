#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "==========================================================\n";
    cout << "Vapor Viscosity Investigation Along Constant P-H Line\n";
    cout << "==========================================================\n\n";
    
    cout << "Investigating Mu_v discontinuity reported at:\n";
    cout << "  X_NaCl = 0.0001\n";
    cout << "  H = 1940 kJ/kg\n";
    cout << "  P = 96 bar → Mu_v = 2.009e-5 Pa·s\n";
    cout << "  P = 103 bar → Mu_v = 8.08e-5 Pa·s (jump of ~4x!)\n\n";
    
    const double X_wt = 0.0001;
    const double H_target = 1940.0e3;  // J/kg
    
    cout << "Scanning pressure from 50 to 150 bar with fine resolution\n";
    cout << "to identify where the viscosity jump occurs.\n\n";
    
    cout << setw(10) << "P (bar)"
         << setw(12) << "Region"
         << setw(12) << "T (°C)"
         << setw(15) << "Mu_v (Pa·s)"
         << setw(15) << "ΔMu_v (%)"
         << setw(12) << "Rho (kg/m³)"
         << setw(10) << "S_l"
         << setw(10) << "S_v" << "\n";
    cout << string(110, '=') << "\n";
    
    double prev_Mu_v = 0.0;
    double prev_P = 0.0;
    int jump_count = 0;
    double max_jump = 0.0;
    double P_at_max_jump = 0.0;
    
    // Coarse scan first (every 2 bar)
    cout << "\nCoarse scan (2 bar steps):\n";
    cout << string(110, '-') << "\n";
    
    for(double P_bar = 50.0; P_bar <= 150.0; P_bar += 2.0) {
        double P_Pa = P_bar * 1e5;
        
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(P_Pa, H_target, X_wt);
        
        double delta_Mu_v = 0.0;
        if(prev_Mu_v > 1e-12 && prop.Mu_v > 1e-12) {
            delta_Mu_v = 100.0 * (prop.Mu_v - prev_Mu_v) / prev_Mu_v;
            
            if(std::abs(delta_Mu_v) > 50.0) {
                jump_count++;
                if(std::abs(delta_Mu_v) > max_jump) {
                    max_jump = std::abs(delta_Mu_v);
                    P_at_max_jump = P_bar;
                }
            }
        }
        
        cout << setw(10) << fixed << setprecision(1) << P_bar
             << setw(12) << prop.Region
             << setw(12) << setprecision(2) << prop.T
             << setw(15) << scientific << setprecision(4) << prop.Mu_v
             << setw(15) << fixed << setprecision(1) << delta_Mu_v
             << setw(12) << setprecision(2) << prop.Rho
             << setw(10) << setprecision(4) << prop.S_l
             << setw(10) << setprecision(4) << prop.S_v;
        
        if(std::abs(delta_Mu_v) > 50.0) {
            cout << "  *** LARGE JUMP ***";
        }
        cout << "\n";
        
        prev_Mu_v = prop.Mu_v;
        prev_P = P_bar;
    }
    
    cout << "\n" << string(110, '=') << "\n";
    cout << "Coarse scan results:\n";
    cout << "  Jumps detected (>50%): " << jump_count << "\n";
    if(jump_count > 0) {
        cout << "  Maximum jump: " << fixed << setprecision(1) << max_jump << "% at P = " << P_at_max_jump << " bar\n";
    }
    cout << "\n";
    
    // Fine scan around the reported problem area (95-105 bar)
    cout << "\nFine scan around reported problem (95-105 bar, 0.5 bar steps):\n";
    cout << string(110, '-') << "\n";
    
    prev_Mu_v = 0.0;
    int fine_jumps = 0;
    
    for(double P_bar = 95.0; P_bar <= 105.0; P_bar += 0.5) {
        double P_Pa = P_bar * 1e5;
        
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(P_Pa, H_target, X_wt);
        
        double delta_Mu_v = 0.0;
        if(prev_Mu_v > 1e-12 && prop.Mu_v > 1e-12) {
            delta_Mu_v = 100.0 * (prop.Mu_v - prev_Mu_v) / prev_Mu_v;
            
            if(std::abs(delta_Mu_v) > 20.0) {
                fine_jumps++;
            }
        }
        
        cout << setw(10) << fixed << setprecision(1) << P_bar
             << setw(12) << prop.Region
             << setw(12) << setprecision(2) << prop.T
             << setw(15) << scientific << setprecision(4) << prop.Mu_v
             << setw(15) << fixed << setprecision(1) << delta_Mu_v
             << setw(12) << setprecision(2) << prop.Rho
             << setw(10) << setprecision(4) << prop.S_l
             << setw(10) << setprecision(4) << prop.S_v;
        
        if(std::abs(delta_Mu_v) > 20.0) {
            cout << "  *** JUMP ***";
        }
        cout << "\n";
        
        prev_Mu_v = prop.Mu_v;
    }
    
    cout << "\n" << string(110, '=') << "\n";
    cout << "Fine scan results:\n";
    cout << "  Jumps detected (>20%): " << fine_jumps << "\n";
    cout << "\n";
    
    // Ultra-fine scan if jumps detected
    if(fine_jumps > 0) {
        cout << "\nUltra-fine scan (99.5-100.5 bar, 0.1 bar steps):\n";
        cout << string(110, '-') << "\n";
        
        prev_Mu_v = 0.0;
        int prev_region = -1;
        
        for(double P_bar = 99.5; P_bar <= 100.5; P_bar += 0.1) {
            double P_Pa = P_bar * 1e5;
            
            H2ONaCl::PROP_H2ONaCl prop = eos.prop_pHX_bisection(P_Pa, H_target, X_wt);
            
            double delta_Mu_v = 0.0;
            if(prev_Mu_v > 1e-12 && prop.Mu_v > 1e-12) {
                delta_Mu_v = 100.0 * (prop.Mu_v - prev_Mu_v) / prev_Mu_v;
            }
            
            cout << setw(10) << fixed << setprecision(2) << P_bar
                 << setw(12) << prop.Region
                 << setw(12) << setprecision(3) << prop.T
                 << setw(15) << scientific << setprecision(6) << prop.Mu_v
                 << setw(15) << fixed << setprecision(2) << delta_Mu_v
                 << setw(12) << setprecision(3) << prop.Rho
                 << setw(10) << setprecision(4) << prop.S_l
                 << setw(10) << setprecision(4) << prop.S_v;
            
            if(std::abs(delta_Mu_v) > 10.0 || prop.Region != prev_region) {
                cout << "  *** ";
                if(prop.Region != prev_region && prev_region >= 0) {
                    cout << "Region change " << prev_region << "→" << prop.Region;
                }
                if(std::abs(delta_Mu_v) > 10.0) {
                    cout << " JUMP";
                }
                cout << " ***";
            }
            cout << "\n";
            
            prev_Mu_v = prop.Mu_v;
            prev_region = prop.Region;
        }
    }
    
    cout << "\n" << string(110, '=') << "\n";
    cout << "Investigation complete.\n\n";
    
    // Specific check at reported conditions
    cout << "Verification of reported conditions:\n";
    cout << string(110, '-') << "\n";
    
    double P1 = 96.0e5;
    double P2 = 103.0e5;
    
    H2ONaCl::PROP_H2ONaCl prop1 = eos.prop_pHX_bisection(P1, H_target, X_wt);
    H2ONaCl::PROP_H2ONaCl prop2 = eos.prop_pHX_bisection(P2, H_target, X_wt);
    
    cout << "At P = 96 bar:\n";
    cout << "  Region: " << prop1.Region << " (" << eos.getPhaseRegionName(prop1.Region) << ")\n";
    cout << "  T = " << fixed << setprecision(2) << prop1.T << " °C\n";
    cout << "  Mu_v = " << scientific << setprecision(4) << prop1.Mu_v << " Pa·s\n";
    cout << "  Rho = " << fixed << setprecision(2) << prop1.Rho << " kg/m³\n";
    cout << "  S_l = " << setprecision(4) << prop1.S_l << ", S_v = " << prop1.S_v << "\n\n";
    
    cout << "At P = 103 bar:\n";
    cout << "  Region: " << prop2.Region << " (" << eos.getPhaseRegionName(prop2.Region) << ")\n";
    cout << "  T = " << fixed << setprecision(2) << prop2.T << " °C\n";
    cout << "  Mu_v = " << scientific << setprecision(4) << prop2.Mu_v << " Pa·s\n";
    cout << "  Rho = " << fixed << setprecision(2) << prop2.Rho << " kg/m³\n";
    cout << "  S_l = " << setprecision(4) << prop2.S_l << ", S_v = " << prop2.S_v << "\n\n";
    
    double ratio = prop2.Mu_v / prop1.Mu_v;
    cout << "Viscosity ratio (P=103/P=96): " << fixed << setprecision(2) << ratio << "x\n";
    cout << "Change: " << setprecision(1) << ((ratio - 1.0) * 100.0) << "%\n";
    
    if(ratio > 2.0) {
        cout << "\n⚠ CONFIRMED: Large viscosity jump detected!\n";
        cout << "   This is likely due to a phase transition or region boundary crossing.\n";
        return 1;
    } else {
        cout << "\n✓ No significant discontinuity found at these conditions.\n";
        return 0;
    }
}
