#include "H2ONaCl.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

// IAPWS-95 Critical Point Constants for Pure Water
const double P_CRIT_BAR = 220.64;      // bar (22.064 MPa)
const double T_CRIT_C = 373.946;       // °C
const double T_CRIT_K = 647.096;       // K
const double RHO_CRIT = 322.0;         // kg/m³
const double H_CRIT_KJ_KG = 2087.546;  // kJ/kg (approximate)

struct CriticalPointTest {
    double P_bar;
    double H_kJ_kg;
    string description;
    int expected_region_pure;  // Expected region for pure water (IAPWS)
    int expected_region_trace; // Expected region for trace salinity (Driesner)
};

string region_name(int region) {
    switch(region) {
        case 0: return "Single phase(Liquid)";
        case 1: return "Single phase(Vapor)";
        case 2: return "Vapor-Liquid Equilibrium";
        case 3: return "Liquid-Halite Equilibrium";
        case 4: return "Vapor-Halite Equilibrium";
        case 5: return "ThreePhase_V_L_H";
        default: return "Unknown/Other";
    }
}

int main()
{
    H2ONaCl::cH2ONaCl eos;
    
    cout << "==========================================================\n";
    cout << "Critical Point Region Identification Test\n";
    cout << "==========================================================\n\n";
    
    cout << "Objective: Verify region identification algorithm around the\n";
    cout << "           critical point of water matches IAPWS behavior\n\n";
    
    cout << "Critical Point of Pure Water (IAPWS-95):\n";
    cout << "  Tc = " << T_CRIT_C << " °C (" << T_CRIT_K << " K)\n";
    cout << "  Pc = " << P_CRIT_BAR << " bar (" << P_CRIT_BAR/10.0 << " MPa)\n";
    cout << "  ρc ≈ " << RHO_CRIT << " kg/m³\n";
    cout << "  Hc ≈ " << H_CRIT_KJ_KG << " kJ/kg\n\n";
    
    const double X_trace = 1.0e-5;  // Trace salinity
    const double X_pure = 0.0;       // Pure water
    
    // Create a grid around the critical point in P-H space
    // Pressure range: Pc ± 30 bar (190-250 bar)
    // Enthalpy range: Hc ± 200 kJ/kg (~1880-2280 kJ/kg)
    
    vector<double> P_offsets = {-30, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 30};  // bar
    vector<double> H_offsets = {-200, -150, -100, -50, -20, -10, 0, 10, 20, 50, 100, 150, 200};  // kJ/kg
    
    int total_tests = 0;
    int region_match = 0;
    int region_mismatch = 0;
    int both_liquid = 0;
    int both_vapor = 0;
    int both_twophase = 0;
    int critical_crossing = 0;
    
    cout << "Scanning P-H space around critical point:\n";
    cout << "  P range: " << (P_CRIT_BAR - 30) << " to " << (P_CRIT_BAR + 30) << " bar\n";
    cout << "  H range: " << (H_CRIT_KJ_KG - 200) << " to " << (H_CRIT_KJ_KG + 200) << " kJ/kg\n";
    cout << "  Grid points: " << P_offsets.size() << " × " << H_offsets.size()
         << " = " << (P_offsets.size() * H_offsets.size()) << " tests\n\n";
    
    cout << "Legend:\n";
    cout << "  Region 0: Single phase liquid\n";
    cout << "  Region 1: Single phase vapor\n";
    cout << "  Region 2: Two-phase (V-L)\n";
    cout << "  X: Region mismatch between pure and trace salinity\n";
    cout << "  ✓: Region match\n\n";
    
    cout << "Results Grid (P vs H):\n";
    cout << string(80, '=') << "\n";
    
    // Print header with H values
    cout << setw(8) << "P(bar)";
    for(double dH : H_offsets) {
        cout << setw(6) << "H" << showpos << setw(4) << (int)dH << noshowpos;
    }
    cout << "\n" << string(80, '-') << "\n";
    
    // Scan the grid
    for(double dP : P_offsets) {
        double P_bar = P_CRIT_BAR + dP;
        double P_Pa = P_bar * 1e5;
        
        cout << setw(8) << fixed << setprecision(1) << P_bar;
        
        for(double dH : H_offsets) {
            double H_kJ_kg = H_CRIT_KJ_KG + dH;
            double H_J_kg = H_kJ_kg * 1e3;
            
            total_tests++;
            
            // Test pure water (X=0)
            H2ONaCl::PROP_H2ONaCl prop_pure;
            try {
                prop_pure = eos.prop_pHX(P_Pa, H_J_kg, X_pure);
            } catch(...) {
                cout << setw(6) << "ERR";
                continue;
            }
            
            // Test trace salinity (X=1e-5)
            H2ONaCl::PROP_H2ONaCl prop_trace;
            try {
                prop_trace = eos.prop_pHX(P_Pa, H_J_kg, X_trace);
            } catch(...) {
                cout << setw(6) << "ERR";
                continue;
            }
            
            // Compare regions
            int region_pure = prop_pure.Region;
            int region_trace = prop_trace.Region;
            
            // Categorize the result
            string symbol;
            if(region_pure == region_trace) {
                region_match++;
                if(region_pure == 0) {
                    both_liquid++;
                    symbol = "L";  // Both liquid
                } else if(region_pure == 1) {
                    both_vapor++;
                    symbol = "V";  // Both vapor
                } else if(region_pure == 2) {
                    both_twophase++;
                    symbol = "2";  // Both two-phase
                } else {
                    symbol = "✓";  // Other matching
                }
            } else {
                region_mismatch++;
                // Check if crossing critical region
                if((region_pure == 0 && region_trace == 1) ||
                   (region_pure == 1 && region_trace == 0) ||
                   (region_pure == 2 || region_trace == 2)) {
                    critical_crossing++;
                    symbol = "*";  // Critical region crossing
                } else {
                    symbol = "X";  // Other mismatch
                }
            }
            
            cout << setw(6) << symbol;
        }
        cout << "\n";
    }
    
    cout << string(80, '=') << "\n\n";
    
    // Detailed statistics
    cout << "Detailed Results:\n";
    cout << string(60, '-') << "\n";
    cout << "Total tests performed:        " << total_tests << "\n";
    cout << "Region matches:               " << region_match
         << " (" << fixed << setprecision(1) << (100.0*region_match/total_tests) << "%)\n";
    cout << "  - Both liquid:              " << both_liquid << "\n";
    cout << "  - Both vapor:               " << both_vapor << "\n";
    cout << "  - Both two-phase:           " << both_twophase << "\n";
    cout << "Region mismatches:            " << region_mismatch
         << " (" << setprecision(1) << (100.0*region_mismatch/total_tests) << "%)\n";
    cout << "  - Critical region crossing: " << critical_crossing << "\n";
    cout << "  - Other mismatches:         " << (region_mismatch - critical_crossing) << "\n";
    cout << string(60, '-') << "\n\n";
    
    // Detailed point-by-point analysis for key locations
    cout << "\n" << string(80, '=') << "\n";
    cout << "Detailed Analysis at Key Points:\n";
    cout << string(80, '=') << "\n\n";
    
    vector<CriticalPointTest> key_tests = {
        // Exactly at critical point
        {P_CRIT_BAR, H_CRIT_KJ_KG, "At critical point (Pc, Hc)", -1, -1},
        
        // Just below critical pressure (subcritical)
        {P_CRIT_BAR - 5, H_CRIT_KJ_KG - 100, "Below Pc, liquid side", 0, 0},
        {P_CRIT_BAR - 5, H_CRIT_KJ_KG, "Below Pc, two-phase", 2, 2},
        {P_CRIT_BAR - 5, H_CRIT_KJ_KG + 100, "Below Pc, vapor side", 1, 1},
        
        // Just above critical pressure (supercritical)
        {P_CRIT_BAR + 5, H_CRIT_KJ_KG - 100, "Above Pc, liquid-like", 0, 0},
        {P_CRIT_BAR + 5, H_CRIT_KJ_KG, "Above Pc, near-critical", 0, 0},
        {P_CRIT_BAR + 5, H_CRIT_KJ_KG + 100, "Above Pc, vapor-like", 0, 0},
        
        // Well below critical pressure
        {P_CRIT_BAR - 20, H_CRIT_KJ_KG - 150, "P << Pc, liquid", 0, 0},
        {P_CRIT_BAR - 20, H_CRIT_KJ_KG + 150, "P << Pc, vapor", 1, 1},
        
        // Well above critical pressure
        {P_CRIT_BAR + 30, H_CRIT_KJ_KG - 150, "P >> Pc, dense fluid", 0, 0},
        {P_CRIT_BAR + 30, H_CRIT_KJ_KG + 150, "P >> Pc, supercritical", 0, 0},
    };
    
    int detailed_tests = 0;
    int detailed_match = 0;
    
    for(const auto& test : key_tests) {
        double P_Pa = test.P_bar * 1e5;
        double H_J_kg = test.H_kJ_kg * 1e3;
        
        cout << "Test: " << test.description << "\n";
        cout << "  P = " << setprecision(2) << test.P_bar << " bar, "
             << "H = " << setprecision(1) << test.H_kJ_kg << " kJ/kg\n";
        
        // Pure water
        H2ONaCl::PROP_H2ONaCl prop_pure = eos.prop_pHX(P_Pa, H_J_kg, X_pure);
        cout << "  Pure water (X=0):   Region " << prop_pure.Region
             << " (" << region_name(prop_pure.Region) << ")\n";
        cout << "    T = " << setprecision(2) << prop_pure.T << " °C, "
             << "ρ = " << setprecision(1) << prop_pure.Rho << " kg/m³\n";
        
        // Trace salinity
        H2ONaCl::PROP_H2ONaCl prop_trace = eos.prop_pHX(P_Pa, H_J_kg, X_trace);
        cout << "  Trace salt (X=1e-5): Region " << prop_trace.Region
             << " (" << region_name(prop_trace.Region) << ")\n";
        cout << "    T = " << setprecision(2) << prop_trace.T << " °C, "
             << "ρ = " << setprecision(1) << prop_trace.Rho << " kg/m³\n";
        
        detailed_tests++;
        if(prop_pure.Region == prop_trace.Region) {
            detailed_match++;
            cout << "  Result: ✓ MATCH\n";
        } else {
            cout << "  Result: ✗ MISMATCH\n";
        }
        
        // Calculate property differences
        double T_diff = abs(prop_trace.T - prop_pure.T);
        double rho_diff_rel = abs(prop_trace.Rho - prop_pure.Rho) / prop_pure.Rho;
        cout << "  ΔT = " << setprecision(3) << T_diff << " °C, "
             << "Δρ/ρ = " << setprecision(2) << (rho_diff_rel * 100) << "%\n";
        
        cout << "\n";
    }
    
    cout << string(80, '=') << "\n";
    cout << "Summary of Detailed Tests:\n";
    cout << "  Total: " << detailed_tests << "\n";
    cout << "  Matches: " << detailed_match << " ("
         << setprecision(1) << (100.0*detailed_match/detailed_tests) << "%)\n";
    cout << "  Mismatches: " << (detailed_tests - detailed_match) << " ("
         << setprecision(1) << (100.0*(detailed_tests-detailed_match)/detailed_tests) << "%)\n";
    cout << string(80, '=') << "\n\n";
    
    // Final assessment
    double match_rate = 100.0 * region_match / total_tests;
    
    if(match_rate >= 95.0) {
        cout << "✓ EXCELLENT: Region identification matches IAPWS in "
             << setprecision(1) << match_rate << "% of tests\n";
        cout << "  The algorithm correctly handles the critical region.\n";
        return 0;
    } else if(match_rate >= 85.0) {
        cout << "✓ GOOD: Region identification matches IAPWS in "
             << setprecision(1) << match_rate << "% of tests\n";
        cout << "  Minor discrepancies near critical point are acceptable.\n";
        return 0;
    } else if(match_rate >= 70.0) {
        cout << "⚠ ACCEPTABLE: Region identification matches IAPWS in "
             << setprecision(1) << match_rate << "% of tests\n";
        cout << "  Significant discrepancies observed near critical point.\n";
        return 0;
    } else {
        cout << "✗ POOR: Region identification matches IAPWS in only "
             << setprecision(1) << match_rate << "% of tests\n";
        cout << "  Algorithm needs improvement for critical region handling.\n";
        return 1;
    }
}
