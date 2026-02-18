#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "Testing findRegion for above-critical conditions\n";
    cout << "=================================================\n\n";
    
    // Test case: Above critical point with small salinity
    double p = 25e6;  // 250 bar (above H2O critical pressure ~220.6 bar)
    double T = 400.0; // 400°C
    double X_wt = 0.01;  // 1% NaCl
    double X_mol = eos.Wt2Mol(X_wt);
    
    cout << "Input: P=" << p/1e5 << " bar, T=" << T << " C, X_wt=" << X_wt << "\n\n";
    
    H2ONaCl::PhaseRegion region = eos.findPhaseRegion_pTX(p, T + 273.15, X_wt);
    
    cout << "Result from findPhaseRegion_pTX:\n";
    cout << "  Region: " << region << " (" << eos.getPhaseRegionName(region) << ")\n\n";
    
    // Now check what prop_pTX returns
    H2ONaCl::PROP_H2ONaCl prop = eos.prop_pTX(p, T + 273.15, X_wt, false);
    cout << "prop_pTX result:\n";
    cout << "  Region: " << prop.Region << " (" << eos.getPhaseRegionName(prop.Region) << ")\n";
    cout << "  Rho: " << prop.Rho << " kg/m3\n";
    cout << "  Rho_l: " << prop.Rho_l << " kg/m3\n";
    cout << "  Rho_v: " << prop.Rho_v << " kg/m3\n";
    cout << "  H: " << prop.H/1e3 << " kJ/kg\n";
    cout << "  S_l: " << prop.S_l << ", S_v: " << prop.S_v << "\n";
    
    return 0;
}
