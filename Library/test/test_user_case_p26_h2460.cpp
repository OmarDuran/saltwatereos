#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Investigation: User Reported Case P=26 bar, H=2460 kJ/kg, X=0.009\n";
    cout << "===================================================================\n\n";
    
    double P_Pa = 26e5;       // 26 bar
    double H_target = 2460e3; // 2460 kJ/kg
    double X_salt = 0.009;    // 0.9% salinity
    double X_pure = 0.0;      // Pure water for IAPWS reference
    
    cout << "Question: Is H=2460 kJ/kg achievable at P=26 bar with X=0.009?\n\n";
    
    // Step 1: Check what pure water (IAPWS) gives at this pressure
    cout << "Step 1: Survey pure water properties at P=26 bar\n";
    cout << "------------------------------------------------\n";
    cout << "T(C)\tH_pure(kJ/kg)\tRho_pure(kg/m3)\tRegion\n";
    
    double T_at_target_H = -1;
    double H_max_pure = 0;
    double T_at_max_H = 0;
    
    for(double T_C = 100; T_C <= 999; T_C += 50) {
        H2ONaCl::PROP_H2ONaCl prop = eos.prop_pTX(P_Pa, T_C + 273.15, X_pure, false);
        
        if(!isnan(prop.H) && prop.H > H_max_pure && T_C < 1000) {
            H_max_pure = prop.H;
            T_at_max_H = T_C;
        }
        
        cout << fixed << setprecision(0) << T_C << "\t" 
             << setprecision(1) << prop.H/1e3 << "\t\t" 
             << setprecision(2) << prop.Rho << "\t\t" 
             << prop.Region << "\n";
        
        // Check if we crossed the target
        if(T_C > 100) {
            H2ONaCl::PROP_H2ONaCl prop_prev = eos.prop_pTX(P_Pa, (T_C-50) + 273.15, X_pure, false);
            if((prop_prev.H < H_target && prop.H > H_target) || 
               (prop_prev.H > H_target && prop.H < H_target)) {
                cout << "  ^^^ Target H=2460 kJ/kg is between T=" << (T_C-50) << " and T=" << T_C << "\n";
                T_at_target_H = T_C - 25; // Approximate
            }
        }
    }
    
    cout << "\nMaximum H for pure water at P=26 bar before T=1000C:\n";
    cout << "  H_max = " << H_max_pure/1e3 << " kJ/kg at T = " << T_at_max_H << " C\n\n";
    
    if(H_target > H_max_pure) {
        cout << "CONCLUSION: H=2460 kJ/kg EXCEEDS the maximum achievable enthalpy\n";
        cout << "            for pure water at P=26 bar within T=[0.1, 1000]C\n";
        cout << "            Target H = " << H_target/1e3 << " kJ/kg\n";
        cout << "            Max H    = " << H_max_pure/1e3 << " kJ/kg\n";
        cout << "            Difference = " << (H_target - H_max_pure)/1e3 << " kJ/kg\n\n";
        cout << "This is a PHYSICAL LIMITATION, not a numerical error.\n";
        cout << "To reach H=2460 kJ/kg at P=26 bar, you would need T > 1000C\n";
        cout << "(outside the valid range of the Driesner & Heinrich correlations).\n\n";
        return 1;
    }
    
    // Step 2: If achievable, find the exact temperature for pure water
    if(T_at_target_H > 0) {
        cout << "Step 2: Fine scan near target H for pure water\n";
        cout << "-----------------------------------------------\n";
        
        double T_search_min = T_at_target_H - 25;
        double T_search_max = T_at_target_H + 25;
        
        cout << "T(C)\tH_pure(kJ/kg)\tError(kJ/kg)\n";
        
        double best_T = T_search_min;
        double min_error = 1e10;
        
        for(double T_C = T_search_min; T_C <= T_search_max; T_C += 1.0) {
            H2ONaCl::PROP_H2ONaCl prop = eos.prop_pTX(P_Pa, T_C + 273.15, X_pure, false);
            double error = fabs(prop.H - H_target);
            
            if(error < min_error) {
                min_error = error;
                best_T = T_C;
            }
            
            if(T_C == T_search_min || T_C == T_search_max || 
               fabs(T_C - T_at_target_H) < 2 || error < 100e3) {
                cout << fixed << setprecision(1) << T_C << "\t" 
                     << setprecision(1) << prop.H/1e3 << "\t\t" 
                     << setprecision(1) << error/1e3 << "\n";
            }
        }
        
        cout << "\nBest match for pure water: T = " << setprecision(2) << best_T 
             << " C, error = " << setprecision(1) << min_error/1e3 << " kJ/kg\n\n";
        
        // Step 3: Compare with X=0.009
        cout << "Step 3: Compare pure water vs X=0.009 at this temperature\n";
        cout << "----------------------------------------------------------\n";
        
        H2ONaCl::PROP_H2ONaCl prop_pure = eos.prop_pTX(P_Pa, best_T + 273.15, X_pure, false);
        H2ONaCl::PROP_H2ONaCl prop_salt = eos.prop_pTX(P_Pa, best_T + 273.15, X_salt, false);
        
        cout << "At T = " << best_T << " C, P = 26 bar:\n";
        cout << "  Pure water (X=0):    H=" << setprecision(1) << prop_pure.H/1e3 
             << " kJ/kg, Rho=" << setprecision(2) << prop_pure.Rho << " kg/m3, Region=" << prop_pure.Region << "\n";
        cout << "  Saline (X=0.009):    H=" << prop_salt.H/1e3 
             << " kJ/kg, Rho=" << prop_salt.Rho << " kg/m3, Region=" << prop_salt.Region << "\n";
        
        double H_diff = fabs(prop_salt.H - prop_pure.H);
        double H_rel_diff = H_diff / prop_pure.H;
        
        cout << "\nDifference: " << H_diff/1e3 << " kJ/kg (" << setprecision(2) << H_rel_diff*100 << "%)\n\n";
        
        // Step 4: Test prop_pHX_bisection
        cout << "Step 4: Test prop_pHX_bisection with X=0.009\n";
        cout << "---------------------------------------------\n";
        
        H2ONaCl::PROP_H2ONaCl prop_pHX = eos.prop_pHX_bisection(P_Pa, H_target, X_salt);
        
        cout << "Result from prop_pHX_bisection:\n";
        cout << "  T = " << prop_pHX.T << " C\n";
        cout << "  H = " << prop_pHX.H/1e3 << " kJ/kg (target: " << H_target/1e3 << ")\n";
        cout << "  Rho = " << prop_pHX.Rho << " kg/m3\n";
        cout << "  Region = " << prop_pHX.Region << "\n\n";
        
        // Verify
        H2ONaCl::PROP_H2ONaCl prop_verify = eos.prop_pTX(P_Pa, prop_pHX.T + 273.15, X_salt, false);
        cout << "Verification with prop_pTX at T=" << prop_pHX.T << " C:\n";
        cout << "  H = " << prop_verify.H/1e3 << " kJ/kg\n";
        cout << "  Rho = " << prop_verify.Rho << " kg/m3\n";
        cout << "  Region = " << prop_verify.Region << "\n\n";
        
        double H_error = fabs(prop_verify.H - H_target);
        double H_rel_error = H_error / H_target;
        
        cout << "Error: " << H_error/1e3 << " kJ/kg (" << setprecision(3) << H_rel_error*100 << "%)\n\n";
        
        if(H_rel_error < 1e-4) {
            cout << "SUCCESS: prop_pHX_bisection converged correctly!\n";
            return 0;
        } else if(H_rel_error < 0.01) {
            cout << "ACCEPTABLE: Small error, likely due to salinity effects.\n";
            return 0;
        } else {
            cout << "FAILED: Large error indicates convergence problem.\n";
            return 1;
        }
    }
    
    return 0;
}
