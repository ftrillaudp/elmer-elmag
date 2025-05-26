/*
Author: F. Trillaud <ftrillaudp@gmail.com>
Date: 04/17/2025
NB: Based on inductor example by C. Geuzaine
*/

Include "transformer.par";
Include "BH.pro"; // nonlinear BH caracteristic of magnetic material

Flag_NL = 1; // if 0 = air (frequency and time domain), 1 == iron (time domain))
Flag_freq = 0; // frequency domain, linear case
If (Flag_freq)
	Flag_NL = 0;
EndIf


Group
{
/// Inductors
Coil_HVr = Region[{coilID_HVr}];
Coil_HVl = Region[{coilID_HVl}];
Coils_HV = Region[{Coil_HVr, Coil_HVl}];
Coil_LVr = Region[{coilID_LVr}];
Coil_LVl = Region[{coilID_LVl}];
Coils_LV = Region[{Coil_LVr, Coil_LVl}];
Coils = Region[{Coils_HV, Coils_LV}];

/// For the core:
Core = Region[{coreID}];

Shell = Region[{shellID}];
Air = Region[{airID, shellID}];
Omega_nf = Region[{Air, Coils}]; // non ferromagnetic materials
Omega_f = Region[{Core}]; // Ferromagnetic material
Omega_CC = Region[{Air, Core}]; // Non conducting region
Omega_C = Region[{Coils}];

If (Flag_NL)
	Omega_nl = Region[{Core}]; // Non linear material
Else
	Omega_nl = Region[{ }];
EndIf

Omega = Region[{Omega_C, Omega_CC}];

/// Boundary:
Gamma_dirichlet = Region[{boundaryID}];
Gamma_treeCotreeGauge = Region[{boundaryID}];

// Empty Groups to be filled
Resistance_Cir  = Region[{}]; // all resistances
Inductance_Cir  = Region[{}] ; // all inductances
Capacitance_Cir = Region[{}] ; // all capacitances
SourceV_Cir = Region[{}]; // all voltage sources
SourceI_Cir = Region[{}]; // all current sources

// Primary side
V_HV = Region[10001]; // arbitrary region number (not linked to the mesh)
SourceV_Cir += Region[{V_HV}];
R_HV = Region[10002]; // arbitrary region number (not linked to the mesh)
Resistance_Cir += Region[{R_HV}];

// Secondary side
R_s = Region[10101]; // arbitrary region number (not linked to the mesh)
Resistance_Cir += Region[{R_s}];
L_l = Region[10102]; // arbitrary region number (not linked to the mesh)
Inductance_Cir += Region[{L_l}];

Impedances_Cir = Region[ {Resistance_Cir, Inductance_Cir} ];
// all circuit sources
Sources_Cir = Region[ {SourceV_Cir, SourceI_Cir} ];
// all circuit elements
Omega_Cir = Region[ {Impedances_Cir, Sources_Cir} ];
}


Function
{
mu0 = 4.0*Pi*1e-7;
nu[Omega_nf] = 1./mu0;
If (Flag_NL)
	nu[Core] = nu_M19[$1];
	dhdb_NL[Core] = dhdb_M19_NL[$1];
	dhdb[Core] = dhdb_M19[$1];
Else
	mu_r = 1500.0;
	nu[Core] = 1./(mu_r*mu0);
EndIf

// To be defined separately for each coil portion, to fix the convention of
// positive current (1: along Oz, -1: along -Oz)
vDir[Coil_HVr] = 1;
vDir[Coil_HVl] = -1;
vDir[Coil_LVr] = 1;
vDir[Coil_LVl] = -1;

N_coils[Coils_HV] = 840;
N_coils[Coils_LV] = 7;

Ae_coils[Coils_HV] = Ae_HV;
Ae_coils[Coils_LV] = Ae_LV; 

freq = 60.;
w = 2*Pi*freq; // pulsation
Je[Coils] = ((N_coils[] / Ae_coils[]) * Vector[0, 0, vDir[]]); // Engineering current density: F_Sin_wt_p[]{w,phase};

t_ini = 0;
nbp = 1; // Number of periods
t_fin = nbp /freq;
nbs = 100; // Number of time steps
dt = t_fin / nbs;

Nb_max_iter = 50;
relaxation_factor = 0.8;
stop_criterion = 1e-8;

NL_tol_abs = 1e-8;
NL_tol_rel = 1e-8;
NL_iter_max = 100;

sigma[Coils] = 5.79e7;

// For a correct definition of the voltage
CoefGeos[Coils] = vDir[] * corethickness;

deg2rad = Pi/180; // degrees to radiants
// Input RMS voltage (half of the voltage because of symmetry; half coils are
// defined!)
val_V_HV = 100.;
phase_V_HV = 0. * deg2rad;

// High value for an open-circuit test; Low value for a short-circuit test;
// any value in-between for any charge
Resistance[R_s] = 0.001;
Inductance[L_l] = 0.0001;

// End-winding primary winding resistance for more realistic primary coil
// model
Resistance[R_HV] = 3.7;
}

Include "integration.pro";
Include "jacobian.pro";


Constraint
{
	{ Name A_constraint; Type Assign;
		Case { { Region Gamma_dirichlet; Value 0.0; } }
	}
	{ Name Current_2D; Case { } }
	{ Name Voltage_2D; Case { } }
	{ Name Current_Cir; Case { } }
	{ Name Voltage_Cir;
		Case {
				{	Region V_HV; Value val_V_HV; 
					TimeFunction F_Sin_wt_p[]{w, phase_V_HV};
				}
		}
	}
	{ Name ElectricalCircuit; Type Network;
		Case Circuit_1 {
						{ Region V_HV; Branch {1,2}; }
						{ Region R_HV; Branch {2,3}; }
						{ Region Coil_HVr; Branch {3,4}; }
						{ Region Coil_HVl; Branch {4,1} ; }
		}
		Case Circuit_2 {
						{ Region R_s; Branch {1,2}; }
						{ Region L_l; Branch {2,3}; }
						{ Region Coil_LVr; Branch {3,4} ; }
						{ Region Coil_LVl; Branch {4,1} ; }
		}
	}
}


FunctionSpace {
	{ Name A_FS; Type Form1P;
		BasisFunction {
			{ Name s_n; NameOfCoef a_n; Function BF_PerpendicularEdge;
			Support Omega; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef a_n; EntityType NodesOf; NameOfConstraint A_constraint; }
		}
	}
	// Current in stranded coil (2D)
	{ Name J_FS; Type Vector;
		BasisFunction {
			{ Name sr; NameOfCoef ir; Function BF_RegionZ; Support Coils; Entity Coils; }
		}
		GlobalQuantity {
			{ Name Is; Type AliasOf; NameOfCoef ir; }
			{ Name Us; Type AssociatedWith; NameOfCoef ir; }
		}
		Constraint {
			{ NameOfCoef Us; EntityType Region; NameOfConstraint Voltage_2D; }
			{ NameOfCoef Is; EntityType Region; NameOfConstraint Current_2D; }
		}
	}
	// UZ and IZ for impedances
	{ Name Impedances_FS; Type Scalar;
		BasisFunction {
			{ Name sr; NameOfCoef ir; Function BF_Region; Support Omega_Cir; Entity Omega_Cir; }
		}
		GlobalQuantity {
			{ Name Iz; Type AliasOf; NameOfCoef ir; }
			{ Name Uz; Type AssociatedWith; NameOfCoef ir; }
		}
		Constraint {
			{ NameOfCoef Uz; EntityType Region; NameOfConstraint Voltage_Cir; }
			{ NameOfCoef Iz; EntityType Region; NameOfConstraint Current_Cir; }
		}
	}
}

// Dynamic Formulation (eddy currents)
Formulation {
	{ Name AV_F; Type FemEquation;
		Quantity {
			{ Name a; Type Local; NameOfSpace A_FS; }

			{ Name ir; Type Local; NameOfSpace J_FS; }
			{ Name Us; Type Global; NameOfSpace J_FS [Us]; }
			{ Name Is; Type Global; NameOfSpace J_FS [Is]; }

			{ Name Uz; Type Global; NameOfSpace Impedances_FS [Uz]; }
			{ Name Iz; Type Global; NameOfSpace Impedances_FS [Iz]; }
		}
		Equation {
			Integral { [ nu[] * Dof{d a} , {d a} ]; In Omega_nf; Jacobian JVol; Integration Integ; }
			If (Flag_NL)
				Integral { [ nu[{d a}] * {d a} , {d a} ]; In Omega_nl; Jacobian JVol; Integration Integ; }
				Integral { [ dhdb[{d a}] * Dof{d a} , {d a} ]; In Omega_nl; Jacobian JVol; Integration Integ; }
				Integral { [ - dhdb[{d a}] * {d a} , {d a} ]; In Omega_nl; Jacobian JVol; Integration Integ; }
			Else
				Integral { [ nu[] * Dof{d a} , {d a} ]; In Omega_f; Jacobian JVol; Integration Integ; }
			EndIf
			
			// js[0] should be of the form: N_coils[]/Ae_coils[] * Vector[0,0,1]
			Integral { [ - (Je[] * Vector[0, 0, 1]) * Dof{ir} , {a} ]; In Coils; Jacobian JVol; Integration Integ; }
			Integral { DtDof [ N_coils[]/Ae_coils[] * Dof{a} , {ir} ]; In Coils; Jacobian JVol; Integration Integ; }
			Integral { 
				[ N_coils[] / Ae_coils[]  / sigma[] * ((Je[] * Vector[0, 0, 1]) * Dof{ir}) , {ir} ]; In Coils;
				Jacobian JVol; Integration Integ;
			}
			
			GlobalTerm { [ Dof{Us} / CoefGeos[] , {Is} ]; In Coils; }

			GlobalTerm { NeverDt[ Dof{Uz} , {Iz} ]; In Resistance_Cir; }
			GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ]; In Resistance_Cir; }

			GlobalTerm { [ Dof{Uz} , {Iz} ]; In Inductance_Cir; }
			GlobalTerm { DtDof [ Inductance[] * Dof{Iz} , {Iz} ]; In Inductance_Cir; }

			GlobalTerm { [ 0. * Dof{Iz} , {Iz} ]; In Sources_Cir; }

			GlobalEquation {
				Type Network; NameOfConstraint ElectricalCircuit;
				{ Node {Is}; Loop {Us}; Equation {Us}; In Coils; }
				{ Node {Iz}; Loop {Uz}; Equation {Uz}; In Omega_Cir; }
			}
		}
	}
}


Resolution {
	{ Name resolution;
		System {
				{ Name AV_S; NameOfFormulation AV_F;
				If(Flag_freq)
					Type ComplexValue; Frequency freq;
				EndIf
				}
		}
		Operation {
			CreateDirectory["resu"]; // create directory to store result files
			If (Flag_freq)
				Generate[AV_S]; Solve[AV_S]; SaveSolution[AV_S];
			Else
				InitSolution[AV_S]; // provide initial condition
				TimeLoopTheta[t_ini, t_fin, dt, 1.]{
				// Euler implicit (1) -- Crank-Nicolson (0.5)
					Print[{$Time}, Format "Time %03g"]; 
					If (NbrRegions[Omega_nl])
						Generate[AV_S]; GetResidual[AV_S, $res0];
						Evaluate[ $res = $res0, $iter = 0 ];
						Print[{$iter, $res, $res / $res0},
						Format "Residual %03g: abs %14.12e rel %14.12e"];
						While[$res > NL_tol_abs && $res / $res0 > NL_tol_rel && $res / $res0 <= 1 && $iter < NL_iter_max]{
							Solve[AV_S]; Generate[AV_S]; GetResidual[AV_S, $res];
							Evaluate[ $iter = $iter + 1 ];
							Print[{$iter, $res, $res / $res0},
							Format "Residual %03g: abs %14.12e rel %14.12e"];
						}
					Else
						Generate[AV_S]; Solve[AV_S];
					EndIf
					SaveSolution[AV_S];
				}
			EndIf
		}
	}
}


PostProcessing
{
	{
		Name postProcessing;
		NameOfFormulation AV_F;
		NameOfSystem AV_S;
		Quantity
		{
			{ Name A; Value { Local { [Norm[{a}]]; In Omega; Jacobian JVol; } } }
			{ Name B; Value { Local { [Norm[{d a}]]; In Omega; Jacobian JVol; } } }
			{ Name J; Value { Term { [ (Je[] * Vector[0, 0, 1]) * {ir} ]; In Coils; Jacobian JVol; } } }
			{ Name vecA; Value { Local{ [{a}]; In Omega; Jacobian JVol; } } }
			{ Name vecB; Value { Local{ [{d a}]; In Omega; Jacobian JVol; } } }
			{ Name Flux; Value { Integral { [ CoefGeos[] * N_coils[] / Ae_coils[] * CompZ[{a}] ]; 
				In Coils; Jacobian JVol; Integration Integ; } } }
		}
	}
}

PostOperation
{
	{
		Name postOperation;
		NameOfPostProcessing postProcessing;
		Operation
		{
			Print[A, OnElementsOf Omega, File "resu/magneticVectorPotential.pos", Name "|A| (T-m)"];
			Print[vecA, OnElementsOf Omega, File "resu/magneticVectorPotentialVector.pos", Name "A (T)"];
			// Add "Smoothing" to be able to see the isolines, does not provide time dependent results
			Print[B, OnElementsOf Omega, File "resu/magneticFluxDensity.pos", Name "|B| (T)"];
			Print[vecB, OnElementsOf Omega, File "resu/magneticFluxDensityVector.pos", Name "B (T)"];
			Print[J, OnElementsOf Coils, File "resu/currentDensities.pos", Name "Je (A-m^-2)"];
		}
	}
}
