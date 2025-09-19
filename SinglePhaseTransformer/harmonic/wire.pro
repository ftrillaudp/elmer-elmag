/*
Author: F. Trillaud <ftrillaudp@gmail.com>
Date: June 2025
NB: Based on inductor example by C. Geuzaine and model library: Lib_Magnetodynamcis2D_av_Cir.pro
*/

Include "wire.par";
Include "BH.pro"; // nonlinear BH caracteristic of magnetic material


Flag_NL = 0; // if 0 = linear (frequency and time domain), 1 == non-linear (time domain))
Flag_freq = 0; // frequency domain, linear case
If (Flag_freq)
	Flag_NL = 0;
EndIf


Group
{
/// Inductors
For i In {1:nbw}
	Omega_Cp~{i} = Region[{wirepID~{i}}];
	//~ Printf("Wire_p", wirepID~{i});
	Omega_Cp += Region[{wirepID~{i}}];
	//Gamma_Cp += Region[{edgewirepID~{i}}];
	Omega_Cn~{i} = Region[{wirenID~{i}}];
	//~ Printf("Wire_n", wirenID~{i});
	Omega_Cn += Region[{wirenID~{i}}];
	Omega_Cn~{i} = Region[{wirenID~{i}}];
	//Gamma_Cn += Region[{edgewirenID~{i}}];
EndFor
Omega_C = Region[{Omega_Cp, Omega_Cn}];
//Gamma_C = Region[{Gamma_Cp, Gamma_Cn}];

/// For the core:
Core = Region[{coreID}];
Omega_C += Region[{Core}];

/// Massive (M), stranded (S) codnuctors
Omega_C_M = Region[{ Omega_C }];
/// Massive (M), stranded (S) conductors
Omega_C_M = Region[{Omega_C}];
Omega_C_S = Region[{ }];

/// Air domain
Shell = Region[{shellID}];
Air = Region[{airID, shellID}];

/// Non ferromagnetic domain
Omega_nf = Region[{Air, Omega_Cp, Omega_Cn}];

/// Ferromagnetic domain
Omega_f = Region[{Core}];

/// Non conductive domain
Omega_CC = Region[{Air}];

/// Linear or non linear computation. If frequency domain, linear model by default
If (Flag_NL)
	Omega_nl = Region[{Core}];
Else
	Omega_nl = Region[{ }];
EndIf

/// Full domain
Omega = Region[{Omega_C, Omega_CC, Omega_f}];

/// Boundary
Gamma_dirichlet = Region[{boundaryID}];

/// Component regions for circuit model
Resistance_Cir  = Region[{ }]; // all resistances
Inductance_Cir  = Region[{ }] ; // all inductances
Capacitance_Cir = Region[{ }] ; // all capacitances
Diode_Cir = Region[{ }] ; // all diodes
SourceV_Cir = Region[{ }]; // all voltage sources
SourceI_Cir = Region[{ }]; // all current sources

/// Circuit
///  _____ PS ______ R _________________________
/// |											|
/// |											|
/// |___ Wire1 ___ Wire2 ____ ... ____ WireN____|

/// Power source (PS)
V_ps = Region[10001]; // arbitrary region number (not linked to the mesh)
SourceV_Cir += Region[{V_ps}];
/// Resitance (R)
R_1 = Region[10002];
Resistance_Cir += Region[{R_1}];
}


Function
{
mu0 = 4.0*Pi*1e-7;
nu[Omega_nf] = 1./mu0;
/// Definition of ferromagnetic material
If (Flag_NL)
	nu[Omega_f] = nu_M19[$1];
	dhdb_NL[Omega_f] = dhdb_M19_NL[$1];
	dhdb[Omega_f] = dhdb_M19[$1];
Else
	mu_r = 1; //1500.0; // M19
	nu[Omega_f] = 1./(mu_r*mu0);
EndIf

// To be defined separately for each coil portion, to fix the convention of
// positive current (1: along Oz, -1: along -Oz)
vDir[Omega_Cp] = 1;
vDir[Omega_Cn] = -1;
vDir[Core] = 1;

// Number of strands in wire
Ns[] = 1; // number of strands
Sc[] = SurfaceArea[]{wirepID_1}; // surface of the bundle of strands, surface of wire

freq = 60.; // frequency
w = 2*Pi*freq; // pulsation
Je[Omega_C_S] = Vector[0, 0, vDir[]]; // Direction of engineering current density

/// Time data
t_ini = 0;
nbp = 1; // Number of periods
t_fin = nbp /freq;
nbs = 100; // Number of time steps
dt = t_fin / nbs; // Time step
deg2rad = Pi/180; // degrees to radiants

/// For a correct definition of the voltage (depth)
CoefGeos[Omega_C] = vDir[] * corethickness;

/// Input RMS voltage
val_V_ps = 10.;
phase_V_ps = 0. * deg2rad;

/// Impedance value of the circuit
Resistance[R_1] = 1.0;
Inductance[] = 0;
Capacitance[] = 0;

/// Definiton of the electrical conductivity for the wires (Cu at 300 K)
sigma[Omega_Cp] = 5.79e7;
sigma[Omega_Cn] = 5.79e7;
sigma[Core] = 1e3;
}


Constraint
{
	{ Name Current_2D; Case { } }
	{ Name Voltage_2D; Case { { Region Core; Value 0.0;} } }
	{ Name Current_Cir; Case { } }
	{ Name Voltage_Cir;
		Case {
				{ 	
					Region V_ps; 
					Value val_V_ps; TimeFunction F_Sin_wt_p[]{w, phase_V_ps}; // Sine function
				}
		}
	}
}

Constraint
{
	{ Name ElectricalCircuit; Type Network;
		Case Circuit_1 {
						{ Region V_ps; Branch {1,2}; }
						{ Region R_1; Branch {2,3}; }
						{ Region Omega_Cp_1; Branch{3, 4}; }
						{ Region Omega_Cn_1; Branch{4, 1}; }
		}
	}
}


Include "integration.pro";
Include "jacobian.pro";
Include "avformulation.pro";
Include "postprocessing.pro";
