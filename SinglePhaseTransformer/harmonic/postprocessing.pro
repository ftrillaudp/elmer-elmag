
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
			{ Name J; 
				Value { 
						Term { [ (Je[] * Vector[0, 0, 1]) * {ir} ]; In Omega_C_S; Jacobian JVol; }
						Term { [ -sigma[] * (Dt[{a}]+{ur}/CoefGeos[]) ]; In Omega_C_M; Jacobian JVol; }
				}
			}
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
			Print[B, OnElementsOf Omega_f, File "resu/magneticFluxDensityCore.pos", Name "|B| (T)"];
			Print[vecB, OnElementsOf Omega_f, File "resu/magneticFluxDensityVectorCore.pos", Name "B (T)"];
			Print[J, OnElementsOf Omega_C, File "resu/currentDensities.pos", Name "Je (A-m^-2)"];
		}
	}
}
