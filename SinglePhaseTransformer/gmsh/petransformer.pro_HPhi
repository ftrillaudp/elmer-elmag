/*
Author: F. Trillaud <ftrillaudp@gmail.com>
Date: May 2025
*/


Include "petransformer.par";


Group {
/// Regions ///
Core = Region[{coreID}];
Air = Region[{airID}];
Shell = Region[{shellID}];
Dirichlet = Region[{boundaryID}];

For i In {1:nbw}
  Omega_Cp~{i} = Region[{wirepID~{i}}];
//  Printf("Wire_p", wirepID~{i});
  Omega_Cp += Region[{wirepID~{i}}];
  Gamma_Cp += Region[{edgewirepID~{i}}];
  Omega_Cn~{i} = Region[{wirenID~{i}}];
//  Printf("Wire_n", wirenID~{i});
  Omega_Cn += Region[{wirenID~{i}}];
  Omega_Cn~{i} = Region[{wirenID~{i}}];
  Gamma_Cn += Region[{edgewirenID~{i}}];
EndFor
Omega_C = Region[{Omega_Cp, Omega_Cn}];
Gamma_C = Region[{Gamma_Cp, Gamma_Cn}];

/// Cuts for positive wires
For i In {1:nbw}
    Cutp~{i} = Region[{(nbID+1+2*nbw-i)}];
    Cutsp += Region[(nbID+1+2*nbw-i)];
//    Printf("Cuts_p", nbID+1+2*nbw-i);
EndFor


/// Cuts for negative wires
For i In {1:nbw}
    Cutn~{i} = Region[{(nbID+1+nbw-i)}];
    Cutsn += Region[(nbID+1+nbw-i)];
//    Printf("Cuts_n", nbID+1+nbw-i);
EndFor
Cuts = Region[{Cutsp, Cutsn}];

/// Domains (set of regions) ///
Omega_CC = Region[{Air, Core, Shell}]; // non conductive domain
Omega = Region[{Omega_CC, Omega_C}]; // Full domain
Omega_nf = Region[{Air, Shell, Omega_C}];
Omega_f = Region[{Core}];

/// Electrical circuitry (dummy regions):
powerSupplyID = 10000;
resistanceID = 20000;
PowerSupply = Region[{powerSupplyID}];
Resistance = Region[{resistanceID}];
Network = Region[{PowerSupply, Resistance}];
Current_Cir = Region[{}];
/*
Elementary electrical circuitry. The resistance of the coil is Rcl.
 ________ PS __________________
|														   |
|														  Rcl
|__cut_n ... cut_2 -- cut_1 ___|
*/
}


Function {
eco_pos = 1;
visualization = 0;
flag_restart = 0;
/// Physical constants:
mu0 = 4*Pi*1e-7;
mur = 1500;

mu[Omega_nf] = mu0;
mu[Omega_f] = mur*mu0;

rho[Omega_C] = 1.86e-6;

freq = 60;
w = 2*Pi*freq; // pulsation
t_ini = 0;
nbp = 1; // Number of periods
t_fin = nbp /freq;
nbs = 100; // Number of time steps
dt = t_fin / nbs;
/// Voltage source:
V0 = 1;
phase_V = 0.;
timeFunction[] = 1; //F_Sin_wt_p[]{w, phase_V};

//Electrical circuitry:
Rcl = 0.1; // resistance of current leads leads
}

Printf("Votage input", V0);

Include "integration.pro";
Include "jacobian.pro";


Constraint {
 { Name H_constraint; Type Assign;
	Case { { Region Dirichlet; Value 0.0; } }
 }
 { Name voltageConstraint; Case {} }
 { Name currentConstraint; Case {} }
 { Name voltageSource; Case {} }
 {
  Name currentSource;
  Case { { Region PowerSupply; Value V0; TimeFunction F_Sin_wt_p[]{w, phase_V}; } }
 }
 {
  Name networkConstraint; Type Network; 
  Case Circuit
  {
    { Region PowerSupply; Branch {0, 1}; }
    { Region Resistance; Branch {1, 2}; }
    // Branching in series the wires. A trick is used to take care of the bug on the cut direction (inverse the sense of the branch).
    //For i In {1:nbw}
      //k1 = 1+i; k2 = 2+i;
      //{ Region Cutp~{i}; Branch{k1, k2}; }
      { Region Cutp_1; Branch{2, 3}; }
      { Region Cutp_2; Branch{3, 4}; }
      { Region Cutp_3; Branch{4, 5}; }
      { Region Cutp_4; Branch{5, 6}; }
      { Region Cutp_5; Branch{6, 7}; }
      { Region Cutp_6; Branch{7, 8}; }
      { Region Cutp_7; Branch{8, 9}; }
      { Region Cutp_8; Branch{9, 10}; }
      { Region Cutp_9; Branch{10, 11}; }
      { Region Cutp_10; Branch{11, 12}; }
      { Region Cutp_11; Branch{12, 13}; }
      { Region Cutp_12; Branch{13, 14}; }
      { Region Cutp_13; Branch{14, 15}; }
      { Region Cutp_14; Branch{15, 16}; }
      { Region Cutp_15; Branch{16, 17}; }
      { Region Cutp_16; Branch{17, 18}; }
      { Region Cutp_17; Branch{18, 19}; }
    //EndFor
   // For i In {1:nbw-1}
    //  k1 = nbw+1+i; k2 = nbw+2+i;
     // { Region Cutn~{i}; Branch{k2, k1}; }
      { Region Cutn_1; Branch{19, 20}; }
      { Region Cutn_2; Branch{20, 21}; }
      { Region Cutn_3; Branch{21, 22}; }
      { Region Cutn_4; Branch{22, 23}; }
      { Region Cutn_5; Branch{23, 24}; }
      { Region Cutn_6; Branch{24, 25}; }
      { Region Cutn_7; Branch{25, 26}; }
      { Region Cutn_8; Branch{26, 27}; }
      { Region Cutn_9; Branch{27, 28}; }
      { Region Cutn_10; Branch{28, 29}; }
      { Region Cutn_11; Branch{29, 30}; }
      { Region Cutn_12; Branch{30, 31}; }
      { Region Cutn_13; Branch{31, 32}; }
      { Region Cutn_14; Branch{32, 33}; }
      { Region Cutn_15; Branch{33, 34}; }
      { Region Cutn_16; Branch{34, 35}; }
      { Region Cutn_17; Branch{35, 0}; }
    //EndFor
    //{ Region Cutn~{nbw}; Branch{0, k2}; }
  }
 }
}


FunctionSpace {
 {
  Name H_FunctionSpace; Type Form1;
  BasisFunction
  {
   {
    Name wh_Edge; NameOfCoef ch_Edge;
    Function BF_Edge; Support Omega_C; Entity EdgesOf[All, Not Gamma_C];
   }
   {
    Name wh_Node; NameOfCoef ch_Node;
    Function BF_GradNode; Support Omega; Entity NodesOf[Omega_CC];
   }
   {
    Name wV_Edge; NameOfCoef cV_Edge;
    Function BF_GroupOfEdges; Support Omega; Entity GroupsOfEdgesOf[Cuts];
   }
  }
  GlobalQuantity
  {
   { Name I_ct; Type AliasOf; NameOfCoef cV_Edge; }
   { Name V_ct; Type AssociatedWith; NameOfCoef cV_Edge; }
  }
  Constraint
  {
    { NameOfCoef ch_Node; EntityType NodesOf; NameOfConstraint H_constraint; }
    { NameOfCoef I_ct; EntityType GroupsOfEdgesOf; NameOfConstraint currentConstraint; }
    { NameOfCoef V_ct; EntityType GroupsOfEdgesOf; NameOfConstraint voltageConstraint; }
  }
 }
 {
  Name network_FunctionSpace; Type Scalar;
  BasisFunction
  {
   {
    Name w_nodal; NameOfCoef c_nodal;
    Function BF_Region; Support Network; Entity Network;
    }
   }
    GlobalQuantity
    {
      { Name I_nt; Type AliasOf; NameOfCoef c_nodal; }
      { Name V_nt; Type AssociatedWith; NameOfCoef c_nodal; }
    }
    Constraint
    {
      { NameOfCoef I_nt; EntityType Region; NameOfConstraint currentSource; }
      { NameOfCoef V_nt; EntityType Region; NameOfConstraint voltageSource; }
    }
  }
}


Formulation {
 {
  Name H_Formulation; Type FemEquation;
  Quantity
  {
   { Name H; Type Local; NameOfSpace H_FunctionSpace; }
   { Name I_ct; Type Global; NameOfSpace H_FunctionSpace[I_ct]; }
   { Name V_ct; Type Global; NameOfSpace H_FunctionSpace[V_ct]; }
    { Name I_nt; Type Global; NameOfSpace network_FunctionSpace[I_nt]; }
    { Name V_nt; Type Global; NameOfSpace network_FunctionSpace[V_nt]; }
  }
  Equation
  {
    Integral
    {
      [ rho[]*{d H}, {d H} ]; In Omega_C;
      Jacobian JVol; Integration Integ;
    }
    Integral
    {
      DtDof[ mu[]*Dof{H}, {H} ]; In Omega;
      Jacobian JVol; Integration Integ;
    }
   // Cuts
   GlobalTerm { [Dof{V_ct},  {I_ct}]; In Cuts; }
   // U = R * I in the resistance on the electrical circuit
   GlobalTerm { NeverDt [ Dof{V_nt}, {I_nt} ]; In Resistance; }
   GlobalTerm { NeverDt [ Rcl*Dof{I_nt}, {I_nt} ]; In Resistance; }
   //~ GlobalTerm { [ 0. * Dof{I_nt} , {I_nt} ]; In Current_Cir; }
   // Connections between the DoFs of the FEM and the network models
   GlobalEquation
   {
    Type Network;
    NameOfConstraint networkConstraint;
    { Node {I_ct}; Loop {V_ct}; Equation {V_ct}; In Cuts; }
    { Node {I_nt}; Loop {V_nt}; Equation {V_nt}; In Network; }
   }
  }
 }
}


Resolution {
	{
 		Name resolution;
		System {{Name H_System; NameOfFormulation H_Formulation;}}
		Operation
		{
      CreateDirectory["resu"]; // create directory to store result files
      InitSolution[H_System]; // provide initial condition
      TimeLoopTheta[t_ini, t_fin, dt, 0.5]{
      // Euler implicit (1) -- Crank-Nicolson (0.5)
        Generate[H_System]; Solve[H_System]; SaveSolution[H_System];
        //~ Test[GetNumberRunTime[visualization]{"Output/PostProcessing/Visualization"}]{PostOperation[onFlightVisualization];}
      }
		}
  }
}


PostProcessing
{
    {
        Name postProcessing;
        NameOfFormulation H_Formulation;
        NameOfSystem H_System;
        Quantity
        {
            {Name currentDensity; Value{Local{[{d H}]; In Omega_C; Jacobian JVol;}}}
            {Name magneticFluxDensityM; Value{Local{[ mu[]*Norm[{H}] ]; In Omega; Jacobian JVol;}}}
            {Name magneticFluxDensityV; Value{Local{[ mu[]*{H} ]; In Omega; Jacobian JVol;}}}
            For i In{1:nbw}
              {Name In_ct~{i}; Value{Term{[ {I_ct} ]; In Cutn~{i};}}}
              {Name Vn_ct~{i}; Value{Term{[ {V_ct} ]; In Cutn~{i};}}}
              {Name Zn_ct~{i}; Value{Term{[ {V_ct}/{I_ct} ]; In Cutn~{i};}}}
              {Name current~{i}; Value{Integral{[ CompZ[{d H}] ]; In Omega_Cn~{i}; Jacobian JVol; Integration Integ;}}}
            EndFor
            {Name totalCurrent; Value{Integral{[ CompZ[{d H}] ]; In Omega_C; Jacobian JVol; Integration Integ;}}}
            {Name losses; Value{Integral{[ rho[]*SquNorm[{d H}] ]; In Omega_C; Jacobian JVol; Integration Integ;}}}
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
          If (eco_pos == 1)
            Print[losses[Omega_C], OnGlobal, Format TimeTable, File "resu/losses.dat", Name "Q (W)" ];
            Print[totalCurrent[Omega_C], OnGlobal, Format TimeTable, File "resu/totalCurrent.dat", Name "Total current (A-turns)"]; 
            Print[currentDensity, OnElementsOf Omega_C, File "resu/currentDensity.pos", Name "J (A-m⁻²)"];
            Print[magneticFluxDensityM, OnElementsOf Omega, File "resu/magneticFluxDensityM.pos", Name "|B| (T)"];
            Print[magneticFluxDensityV, OnElementsOf Omega, File "resu/magneticFluxDensityV.pos", Name "B (T)"];
          Else
            Print[magneticFluxDensityM, OnElementsOf Omega, File "resu/magneticFluxDensityM.pos", Name "|B| (T)"];
            Print[magneticFluxDensityV, OnElementsOf Omega, File "resu/magneticFluxDensityV.pos", Name "B (T)"];
            Print[magneticFluxDensityX, OnElementsOf Omega, File "resu/magneticFluxDensityX.pos", Name "|BX| (T)"];
            Print[magneticFluxDensityY, OnElementsOf Omega, File "resu/magneticFluxDensityY.pos", Name "|BY| (T)"];
            Print[currentDensity, OnElementsOf Omega_C, File "resu/currentDensity.pos", Name "J (A-m⁻²)"];
            Print[magneticFluxDensityA[Omega_C], OnGlobal, Format TimeTable, File "resu/magneticFluxDensityA.dat", Name "B_{av} (T)" ];
            Print[totalCurrent[Omega_C], OnGlobal, Format TimeTable, File "resu/totalCurrent.dat", Name "Total current (A-turns)"];
            Print[losses[Omega_C], OnGlobal, Format TimeTable, File "resu/losses.dat", Name "Q (W)" ];
            For i In{1:nbw}
              Print[I_ct~{i}, OnRegion Cutn~{i}, Format TimeTable, File "resu/Icut_", AppendExpressionToFileName i, AppendExpressionFormat "%g.dat"];
              Print[V_ct~{i}, OnRegion Cutn~{i}, Format TimeTable, File "resu/Vcut_", AppendExpressionToFileName i, AppendExpressionFormat "%g.dat"];
              Print[Z_ct~{i}, OnRegion Cutn~{i}, Format TimeTable, File "resu/Zcut_", AppendExpressionToFileName i, AppendExpressionFormat "%g.dat"];
              Print[current~{i}[Omega_Cn~{i}], OnGlobal, Format TimeTable, File "resu/current_", AppendExpressionToFileName i, AppendExpressionFormat "%g.dat"];
            EndFor
          EndIf
        }
    }
    {
        Name onFlightVisualization;
        NameOfPostProcessing postProcessing;
				LastTimeStepOnly visualization;
        Operation
        {
            If (eco_pos == 1)
              Print[losses[Omega_C], OnGlobal, Format Table, LastTimeStepOnly, File > "resu/losses.dat", SendToServer "Output/PostProcessing/1Losses (W)", Color "Red" ];
              Print[currentDensity, OnElementsOf Omega_C, File "resu/currentDensity.pos", Name "J (A-m⁻²)"];
              Print[magneticFluxDensityM, OnElementsOf Omega, File "resu/magneticFluxDensityM.pos", Name "|B| (T)"];
              Print[magneticFluxDensityV, OnElementsOf Omega, File "resu/magneticFluxDensityV.pos", Name "B (T)"];
            Else
              Print[currentDensity, OnElementsOf Omega_C, File "resu/currentDensity.pos", Name "J (A-m⁻²)"];
              Print[magneticFluxDensityM, OnElementsOf Omega, File "resu/magneticFluxDensityM.pos", Name "|B| (T)"];
              Print[totalCurrent[Omega_C], OnGlobal, Format TimeTable, LastTimeStepOnly, File "resu/lastTimeStep-totalCurrent.dat", SendToServer "Output/PostProcessing/0Total current (A-turns)", Color "Yellow"];
              Print[losses[Omega_C], OnGlobal, Format Table, LastTimeStepOnly, File > "resu/losses.dat", SendToServer "Output/PostProcessing/1Losses (W)", Color "Red" ];
            EndIf
        }
    }
}

If (flag_restart == 0)
  DefineConstant[ C_ = {"-solve -bin -v 3 -v2 -cpu", Name "GetDP/9ComputeCommand", Visible 1} ];
Else
  DefineConstant[ C_ = {"-solve -bin -v 3 -v2 -cpu -restart", Name "GetDP/9ComputeCommand", Visible 1} ];
EndIf
