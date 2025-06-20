
Group
{
	DefineGroup[
		Omega_C, // Conductor domain
		Gamma_C, // Surface of conductors
		Core, // Core domain
		Air, // "air" domain (include Shell domain)
		Shell, // shell transform domain
		Omega_C_M, // Massive conductors
		Omega_C_S, // Stranded conductors
		Omega_CC, // Non conducting domain
		Omega_nf, // Non ferromagnetic domain
		Omega_f, // Ferromagnetic domain
		Omega_nl, // Non linear domain
		Omega, // linear domain
		Gamma_dirichlet, // Boudary condition
		Resistance_Cir, // all resistances
		Inductance_Cir, // all inductances
		Capacitance_Cir, // all capacitances
		Diode_Cir, // all diodes
		SourceV_Cir, // all voltage sources
		SourceI_Cir, // all current sources
		Impedances_Cir, // Impedances
		Sources_Cir, // Sources
		Omega_Cir // All circuit domains
	];	
}


Group
{
	/// Circuit impedance
	Impedances_Cir = Region[ {Resistance_Cir, Inductance_Cir, Capacitance_Cir, Diode_Cir} ];
	/// all circuit sources
	Sources_Cir = Region[ {SourceV_Cir, SourceI_Cir} ];
	/// all circuit elements
	Omega_Cir = Region[ {Impedances_Cir, Sources_Cir} ];
}


Function
{
	/// Parameters for non linear iterative scheme
	Nb_max_iter = 50;
	relaxation_factor = 0.8;
	stop_criterion = 1e-8;
	NL_tol_abs = 1e-8;
	NL_tol_rel = 1e-8;
	NL_iter_max = 100;
}


Constraint
{
	{ Name A_constraint; Type Assign; Case { { Region Gamma_dirichlet; Value 0.0; } } }
}

/// A formulation
FunctionSpace {
	{ Name A_FS; Type Form1P;
		BasisFunction {
			{ Name s_n; NameOfCoef a_n; Function BF_PerpendicularEdge; Support Omega; Entity NodesOf[All]; }
			{ Name s_e; NameOfCoef a_e; Function BF_PerpendicularEdge_2E; Support Omega; Entity EdgesOf[All]; }
		}
		Constraint {
			{ NameOfCoef a_n; EntityType NodesOf; NameOfConstraint A_constraint; }
			{ NameOfCoef a_e; EntityType EdgesOf; NameOfConstraint A_constraint; }
		}
	}
}

/// Massive conductors: e = -Grad(V)
FunctionSpace {

	{ Name U_FS; Type Form1P;
		BasisFunction {
			{ Name sr; NameOfCoef ur; Function BF_RegionZ; Support Omega_C_M; Entity Omega_C_M; }
		}
		GlobalQuantity {
			{ Name U; Type AliasOf; NameOfCoef ur; }
			{ Name I; Type AssociatedWith; NameOfCoef ur; }
		}
		Constraint {
			{ NameOfCoef U; EntityType Region; NameOfConstraint Voltage_2D; }
			{ NameOfCoef I;	EntityType Region; NameOfConstraint Current_2D; }
		}
	}
}

/// Stranded conductors
FunctionSpace {

	{ Name J_FS; Type Vector;
		BasisFunction {
			{ Name sr; NameOfCoef ir; Function BF_RegionZ; Support Omega_C_S; Entity Omega_C_S; }
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
}

/// UZ and IZ for impedances
FunctionSpace {
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

			{ Name ur; Type Local; NameOfSpace U_FS; }
			{ Name I; Type Global; NameOfSpace U_FS [I]; }
			{ Name U; Type Global; NameOfSpace U_FS [U]; }
			
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
			
			/// Massive conductors: e = -Dt[{a}]-{ur}, with {ur} = Grad(v) constant in each region of Massive conductors
			Integral { DtDof [ sigma[] * Dof{a} , {a} ]; In Omega_C_M; Jacobian JVol; Integration Integ; }
			Integral { [ sigma[] * Dof{ur} / CoefGeos[] , {a} ]; In Omega_C_M; Jacobian JVol; Integration Integ; }
			Integral { DtDof [ sigma[] * Dof{a} , {ur} ]; In Omega_C_M; Jacobian JVol; Integration Integ; }
			Integral { [ sigma[] * Dof{ur} / CoefGeos[] , {ur} ]; In Omega_C_M; Jacobian JVol; Integration Integ; }
			
			/// Stranded conductors: Je * Vector[0,0,1]
			Integral { [ - (Je[] * Vector[0, 0, 1]) * Dof{ir} , {a} ]; In Omega_C_S; Jacobian JVol; Integration Integ; }
			Integral { DtDof [ (Ns[] / Sc[]) * Dof{a} , {ir} ]; In Omega_C_S; Jacobian JVol; Integration Integ; }
			Integral { 
				[ (Ns[] / Sc[])  / sigma[] * ((Je[] * Vector[0, 0, 1]) * Dof{ir}) , {ir} ]; In Omega_C_S;
				Jacobian JVol; Integration Integ;
			}
			
			/// Coupling FEM-circuits
			/// Massive conductors
			GlobalTerm { [ Dof{I} *(CoefGeos[]/Fabs[CoefGeos[]]) , {U} ]; In Omega_C_M; }
			/// Stranded conductors
			GlobalTerm { [ Dof{Us} / CoefGeos[] , {Is} ]; In Omega_C_S; }
			/// Resistances
			GlobalTerm { NeverDt[ Dof{Uz} , {Iz} ]; In Resistance_Cir; }
			GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ]; In Resistance_Cir; }
			/// Inductances
			GlobalTerm { [ Dof{Uz} , {Iz} ]; In Inductance_Cir; }
			GlobalTerm { DtDof [ Inductance[] * Dof{Iz} , {Iz} ]; In Inductance_Cir; }
			/// Capacitances
			GlobalTerm { NeverDt[ Dof{Iz} , {Iz} ]; In Capacitance_Cir; }
			GlobalTerm { DtDof [ Capacitance[] * Dof{Uz} , {Iz} ]; In Capacitance_Cir; }
			/// Diodes
			GlobalTerm { NeverDt[ Dof{Uz} , {Iz} ]; In Diode_Cir; }
			GlobalTerm { NeverDt[ Resistance[{Uz}] * Dof{Iz} , {Iz} ]; In Diode_Cir; }
			/// Current sources
			GlobalTerm { [ 0. * Dof{Iz} , {Iz} ]; In Sources_Cir; }

			/// Circuit equations
			GlobalEquation {
				Type Network; NameOfConstraint ElectricalCircuit;
				{ Node {I};  Loop {U};  Equation {I};  In Omega_C_M; }
				{ Node {Is}; Loop {Us}; Equation {Us}; In Omega_C_S; }
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
