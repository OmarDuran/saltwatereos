#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    const double X_test = 1.0e-4;
    const double X_pure = 0.0;
    
    cout << "prop_pHX_bisection Boundary Test at X=0.0001\n";
    cout << "=============================================\n\n";
    
    int total_tests = 0;
    int failures = 0;
    
    // Test along saturation curve using prop_pHX
    for(double T_C = 150.0; T_C <= 350.0; T_C += 25.0) {
        double T_K = T_C + 273.15;
        double P_sat_bar = eos.m_water.P_Boiling(T_C);
        double P_sat_Pa = P_sat_bar * 1e5;
        
        // Get IAPWS saturation enthalpy (vapor side)
        H2ONaCl::PROP_H2ONaCl prop_iapws = eos.prop_pTX(P_sat_Pa, T_K, X_pure, false);
        double H_v_sat = prop_iapws.H_v;
        
        // Test prop_pTX at saturation
        H2ONaCl::PROP_H2ONaCl prop_pTX = eos.prop_pTX(P_sat_Pa, T_K, X_test, false);
        
        // Test prop_pHX at saturation
        H2ONaCl::PROP_H2ONaCl prop_pHX = eos.prop_pHX(P_sat_Pa, H_v_sat, X_test);
        
        cout << "T=" << fixed << setprecision(1) << T_C << "C, P=" << setprecision(2) << P_sat_bar << " bar:\n";
        cout << "  prop_pTX: Region=" << prop_pTX.Region << ", S_v=" << setprecision(4) << prop_pTX.S_v 
             << ", S_l=" << prop_pTX.S_l << ", Rho=" << setprecision(3) << prop_pTX.Rho << " kg/m³\n";
        cout << "  prop_pHX: Region=" << prop_pHX.Region << ", S_v=" << setprecision(4) << prop_pHX.S_v 
             << ", S_l=" << prop_pHX.S_l << ", Rho=" << setprecision(3) << prop_pHX.Rho << " kg/m³\n";
        
        total_tests++;
        
        // Check consistency
        bool region_match = (prop_pTX.Region == prop_pHX.Region);
        bool sv_match = fabs(prop_pTX.S_v - prop_pHX.S_v) < 0.01;
        bool sl_match = fabs(prop_pTX.S_l - prop_pHX.S_l) < 0.01;
        
        if(!region_match || !sv_match || !sl_match) {
            cout << "  ⚠ MISMATCH!\n";
            failures++;
        } else if(prop_pHX.Region != 2 || fabs(prop_pHX.S_v - 1.0) > 0.01) {
            cout << "  ⚠ INCORRECT: Should be Region=2, S_v=1.0\n";
            failures++;
        } else {
            cout << "  ✓ OK\n";
        }
        cout << "\n";
    }
    
    cout << "Summary: " << total_tests << " tests, " << failures << " failures\n";
    
    if(failures == 0) {
        cout << "\n✓✓✓ SUCCESS: prop_pHX_bisection is consistent with prop_pTX!\n";
        return 0;
    } else {
        cout << "\n⚠ FAILED: prop_pHX_bisection has inconsistencies\n";
        return 1;
    }
}
