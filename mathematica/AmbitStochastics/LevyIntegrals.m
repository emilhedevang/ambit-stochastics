(* Mathematica Package *)

BeginPackage["AmbitStochastics`LevyIntegrals`"]
(* Exported symbols added here with SymbolName::usage *)  

LevyIntegral::usage              = "FIXME"

Begin["`Private`"] (* Begin Private Context *) 


Clear[LevyIntegral]
LevyIntegral[kernel_ ? VectorQ, seed_, delta_, m_] :=
    ListConvolve[kernel, RandomVariate[seed[delta], m + Length[kernel] - 1]];

LevyIntegral[kernel_Function | kernel_Symbol, seed_, delta_, n_, m_] :=
    ListConvolve[Map[kernel, delta * (Range[0, n-1] + 0.5)], RandomVariate[seed[delta], m + n - 1]];


End[] (* End Private Context *)

EndPackage[]