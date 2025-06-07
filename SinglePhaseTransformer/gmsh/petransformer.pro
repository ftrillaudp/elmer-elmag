/*
Author: F. Trillaud <ftrillaudp@gmail.com>
Date: June 2025
NB: Based on inductor example by C. Geuzaine
*/

Include "petransformer.par";
Include "BH.pro"; // nonlinear BH caracteristic of magnetic material

Flag_NL = 1; // if 0 = linear (frequency and time domain), 1 == non-linear (time domain))
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
	Gamma_Cp += Region[{edgewirepID~{i}}];
	Omega_Cn~{i} = Region[{wirenID~{i}}];
	//~ Printf("Wire_n", wirenID~{i});
	Omega_Cn += Region[{wirenID~{i}}];
	Omega_Cn~{i} = Region[{wirenID~{i}}];
	Gamma_Cn += Region[{edgewirenID~{i}}];
EndFor
Omega_C = Region[{Omega_Cp, Omega_Cn}];
Gamma_C = Region[{Gamma_Cp, Gamma_Cn}];

/// For the core:
Core = Region[{coreID}];
Omega_C += Region[{Core}];

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
Omega = Region[{Omega_C, Omega_CC}];

/// Boundary
Gamma_dirichlet = Region[{boundaryID}];

/// Component regions for circuit model
Resistance_Cir  = Region[{}]; // all resistances
Inductance_Cir  = Region[{}] ; // all inductances
Capacitance_Cir = Region[{}] ; // all capacitances
SourceV_Cir = Region[{}]; // all voltage sources
SourceI_Cir = Region[{}]; // all current sources

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
/// Circuit impedance
Impedances_Cir = Region[ {Resistance_Cir, Inductance_Cir} ];
/// all circuit sources
Sources_Cir = Region[ {SourceV_Cir, SourceI_Cir} ];
/// all circuit elements
Omega_Cir = Region[ {Impedances_Cir, Sources_Cir} ];
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
	mu_r = 1500.0; // M19
	nu[Omega_f] = 1./(mu_r*mu0);
EndIf

// To be defined separately for each coil portion, to fix the convention of
// positive current (1: along Oz, -1: along -Oz)
vDir[Omega_Cp] = 1;
vDir[Omega_Cn] = -1;
vDir[Core] = 1;

freq = 60.; // frequency
w = 2*Pi*freq; // pulsation
Je[Omega_C] = Vector[0, 0, vDir[]]; // Direction of engineering current density

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
/// Resistance value of the circuit
Resistance[R_1] = 3.7;

/// Definiton of the electrical conductivity for the wires (Cu at 300 K)
sigma[Omega_Cp] = 5.79e7;
sigma[Omega_Cn] = 5.79e7;
sigma[Core] = 1e3;

/// Parameters for non linear iterative scheme
Nb_max_iter = 50;
relaxation_factor = 0.8;
stop_criterion = 1e-8;
NL_tol_abs = 1e-8;
NL_tol_rel = 1e-8;
NL_iter_max = 100;
}


Include "integration.pro";

Include "jacobian.pro";


Constraint
{
	{ Name A_constraint; Type Assign;
		Case { { Region Gamma_dirichlet; Value 0.0; } }
	}
	{ Name Current_2D; Case { } }
	{ Name Voltage_2D; Case { { Region Core; Value 0; } } }
	{ Name Current_Cir; Case { } }
	{ Name Voltage_Cir;
		Case {
				{	Region V_ps; Value val_V_ps; 
					TimeFunction F_Sin_wt_p[]{w, phase_V_ps}; // Sine function
				}
		}
	}
	{ Name ElectricalCircuit; Type Network;
		Case Circuit_1 {
						{ Region V_ps; Branch {1,2}; }
						{ Region R_1; Branch {2,3}; }
						For i In {1:nbw}
							k1 = 3+i-1; k2 = 3+i;
							//~ Printf("k1", k1); Printf("k2", k2);
							{ Region Omega_Cp~{i}; Branch{k1, k2}; }
							k3 = 3+nbw+i-1; k4 = 3+nbw+i;
							If (k4 < 3+2*nbw)
								{ Region Omega_Cn~{i}; Branch{k3, k4}; }
								//~ Printf("k3", k3); Printf("k4", k4);
							Else
								{ Region Omega_Cn~{i}; Branch{k3, 1}; }
								//~ Printf("k3", k3); Printf("k4", 1);
							EndIf
						EndFor
						//~ { Region Omega_Cp_1; Branch{3, 4}; }
						//~ { Region Omega_Cp_2; Branch{4, 5}; }
						//~ { Region Omega_Cp_3; Branch{5, 6}; }
						//~ { Region Omega_Cp_4; Branch{6, 7}; }
						//~ { Region Omega_Cp_5; Branch{7, 8}; }
						//~ { Region Omega_Cp_6; Branch{8, 9}; }
						//~ { Region Omega_Cp_7; Branch{9, 10}; }
						//~ { Region Omega_Cp_8; Branch{10, 11}; }
						//~ { Region Omega_Cp_9; Branch{11, 12}; }
						//~ { Region Omega_Cp_10; Branch{12, 13}; }
						//~ { Region Omega_Cp_11; Branch{13, 14}; }
						//~ { Region Omega_Cp_12; Branch{14, 15}; }
						//~ { Region Omega_Cp_13; Branch{15, 16}; }
						//~ { Region Omega_Cp_14; Branch{16, 17}; }
						//~ { Region Omega_Cp_15; Branch{17, 18}; }
						//~ { Region Omega_Cp_16; Branch{18, 19}; }
						//~ { Region Omega_Cp_17; Branch{19, 20}; }
						
						//~ { Region Omega_Cn_1; Branch{20, 21}; }
						//~ { Region Omega_Cn_2; Branch{21, 22}; }
						//~ { Region Omega_Cn_3; Branch{22, 23}; }
						//~ { Region Omega_Cn_4; Branch{23, 24}; }
						//~ { Region Omega_Cn_5; Branch{24, 25}; }
						//~ { Region Omega_Cn_6; Branch{25, 26}; }
						//~ { Region Omega_Cn_7; Branch{26, 27}; }
						//~ { Region Omega_Cn_8; Branch{27, 28}; }
						//~ { Region Omega_Cn_9; Branch{28, 29}; }
						//~ { Region Omega_Cn_10; Branch{29, 30}; }
						//~ { Region Omega_Cn_11; Branch{30, 31}; }
						//~ { Region Omega_Cn_12; Branch{31, 32}; }
						//~ { Region Omega_Cn_13; Branch{32, 33}; }
						//~ { Region Omega_Cn_14; Branch{33, 34}; }
						//~ { Region Omega_Cn_15; Branch{34, 35}; }
						//~ { Region Omega_Cn_16; Branch{35, 36}; }
						//~ { Region Omega_Cn_17; Branch{36, 1}; }
		}
	}
}


FunctionSpace {
	{ Name A_FS; Type Form1P;
		BasisFunction {
			{ Name s_n; NameOfCoef a_n; Function BF_PerpendicularEdge;
			Support Omega; Entity NodesOf[All]; 
			}
		}
		Constraint {
			{ NameOfCoef a_n; EntityType NodesOf; NameOfConstraint A_constraint; }
		}
	}
	/// Current in stranded coil (2D)
	{ Name J_FS; Type Vector;
		BasisFunction {
			{ Name sr; NameOfCoef ir; Function BF_RegionZ; Support Omega_C; Entity Omega_C; }
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
	/// UZ and IZ for impedances
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


/// Dynamic Formulation
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
				Integral { [ + dhdb[{d a}] * Dof{d a} , {d a} ]; In Omega_nl; Jacobian JVol; Integration Integ; }
				Integral { [ - dhdb[{d a}] * {d a} , {d a} ]; In Omega_nl; Jacobian JVol; Integration Integ; }
			Else
				Integral { [ nu[] * Dof{d a} , {d a} ]; In Omega_f; Jacobian JVol; Integration Integ; }
			EndIf
			
			Integral { DtDof [ sigma[] * Dof{a} , {a} ]; In Omega_C; Jacobian JVol; Integration Integ; }
			
			/// js[0] should be of the form: Je * Vector[0,0,1]
			Integral { [ - (Je[] * Vector[0, 0, 1]) * Dof{ir} , {a} ]; In Omega_C; Jacobian JVol; Integration Integ; }
			Integral { DtDof [ Dof{a} , {ir} ]; In Omega_C; Jacobian JVol; Integration Integ; }
			Integral { 
				[ 1  / sigma[] * ((Je[] * Vector[0, 0, 1]) * Dof{ir}) , {ir} ]; In Omega_C;
				Jacobian JVol; Integration Integ;
			}
			
			/// Inductor
			GlobalTerm { [ Dof{Us} / CoefGeos[] , {Is} ]; In Omega_C; }
			/// Resistances
			GlobalTerm { NeverDt[ Dof{Uz} , {Iz} ]; In Resistance_Cir; }
			GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ]; In Resistance_Cir; }
			/// Inductances
			//GlobalTerm { [ Dof{Uz} , {Iz} ]; In Inductance_Cir; }
			//GlobalTerm { DtDof [ Inductance[] * Dof{Iz} , {Iz} ]; In Inductance_Cir; }
			/// Current sources
			GlobalTerm { [ 0. * Dof{Iz} , {Iz} ]; In Sources_Cir; }

			GlobalEquation {
				Type Network; NameOfConstraint ElectricalCircuit;
				{ Node {Is}; Loop {Us}; Equation {Us}; In Omega_C; }
				{ Node {Iz}; Loop {Uz}; Equation {Uz}; In Omega_Cir; }
			}
		}
	}
}


/// Solver
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


/// Derivation of the results
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
			{ Name J; Value { Term { [ (Je[] * Vector[0, 0, 1]) * {ir} ]; In Omega_C; Jacobian JVol; } } }
			{ Name vecA; Value { Local{ [{a}]; In Omega; Jacobian JVol; } } }
			{ Name vecB; Value { Local{ [{d a}]; In Omega; Jacobian JVol; } } }
			{ Name Flux; Value { Integral { [ CoefGeos[] * CompZ[{a}] ]; 
				In Omega_C; Jacobian JVol; Integration Integ; } } }
		}
	}
}


/// Saving results
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
			Print[J, OnElementsOf Omega_C, File "resu/currentDensities.pos", Name "Je (A-m^-2)"];
		}
	}
}
