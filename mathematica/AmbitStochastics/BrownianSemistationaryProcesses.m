(* Mathematica Package *)

BeginPackage["AmbitStochastics`BrownianSemistationaryProcesses`", {"AmbitStochastics`Utilities`"}]
(* Exported symbols added here with SymbolName::usage *)  


PrepareBSSBrownianCoefficients::usage      = "FIXME";
PrepareBSSDeterministicCoefficients::usage = "FIXME";
SimulateBSSDeterministicIntegral::usage    = "FIXME";
SimulateBSSBrownianIntegral::usage         = "FIXME";


Begin["`Private`"] (* Begin Private Context *) 

Clear[PrepareBSSDeterministicCoefficients]
PrepareBSSDeterministicCoefficients[kernel_, delta_, n_Integer, order : (0 | 1), tbl : (Table | ParallelTable) : Table] /; delta > 0 && n > 0 :=
    Switch[order,
        0, 
        {
            tbl[delta * NIntegrate[kernel[(j + u) * delta], {u, 0, 1}], {j, 0, n - 1}]
        },
        1,
        {
            tbl[delta * NIntegrate[kernel[(j + u) * delta] * u,       {u, 0, 1}], {j, 0, n - 1}],
            tbl[delta * NIntegrate[kernel[(j + u) * delta] * (1 - u), {u, 0, 1}], {j, 0, n - 1}]
        }
    ];
    
PrepareBSSDeterministicCoefficients /: Parallelize[PrepareBSSDeterministicCoefficients[kernel_, delta_, n_, order_]] :=
    PrepareBSSDeterministicCoefficients[kernel, delta, n, order, ParallelTable];

Clear[PrepareBSSBrownianCoefficients]
PrepareBSSBrownianCoefficients[kernel_, delta_, n_Integer, order : (0 | 1), tbl : (Table | ParallelTable) : Table] /; delta > 0 && n > 0 :=
    Switch[order,
        0, 
        {
    	   Sqrt @ tbl[delta * NIntegrate[kernel[(j + u) * delta] ^ 2, {u, 0, 1}], {j, 0, n - 1}]
        },
    	1,
    	{
    		Sqrt @ tbl[delta * NIntegrate[kernel[(j + u) * delta] ^ 2 * u,       {u, 0, 1}], {j, 0, n - 1}],
            Sqrt @ tbl[delta * NIntegrate[kernel[(j + u) * delta] ^ 2 * (1 - u), {u, 0, 1}], {j, 0, n - 1}]
    	}
    ];

PrepareBSSBrownianCoefficients /: Parallelize[PrepareBSSBrownianCoefficients[kernel_, delta_, n_, order_]] :=
    PrepareBSSBrownianCoefficients[kernel, delta, n, order, ParallelTable];

Clear[SimulateBSSDeterministicIntegral]
SimulateBSSDeterministicIntegral[{G0_ ? NumberVectorQ}, a_] :=
    ListConvolve[G0, a];

SimulateBSSDeterministicIntegral[{G0_ ? NumberVectorQ, G1_ ? NumberVectorQ}, a_] :=
    ListConvolve[G0, Most[a]] + ListConvolve[G1, Rest[a]];

Clear[SimulateBSSBrownianIntegral]
SimulateBSSBrownianIntegral[{G0_ ? NumberVectorQ}, sigma_] :=
    ListConvolve[G0, sigma * RandomVariate[NormalDistribution[], {Length[sigma]}]];

SimulateBSSBrownianIntegral[{G0_ ? NumberVectorQ, G1_ ? NumberVectorQ}, sigma_] :=
    ListConvolve[G0, Most[sigma] * RandomVariate[NormalDistribution[], {Length[sigma] - 1}]] + 
    ListConvolve[G1, Rest[sigma] * RandomVariate[NormalDistribution[], {Length[sigma] - 1}]];


End[] (* End Private Context *)

EndPackage[]