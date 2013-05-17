(* Mathematica Package *)

BeginPackage["AmbitStochastics`LevyBases`"]
(* Exported symbols added here with SymbolName::usage *)  

GaussianSeed::usage              = "FIXME";
NormalInverseGaussianSeed::usage = "FIXME";
GeneralisedHyperbolicSeed::usage = "FIXME";

Begin["`Private`"] (* Begin Private Context *) 
Clear[GaussianSeed]
GaussianSeed[mean_, variance_][volume_] :=
    NormalDistribution[mean * volume, Sqrt[variance * volume]];

Clear[NormalInverseGaussianSeed]
NormalInverseGaussianSeed[alpha_, beta_, mu_, delta_][volume_] :=
    GeneralisedHyperbolicSeed[-1/2, alpha, beta, mu, delta][volume];

Clear[GeneralisedHyperbolicSeed]
GeneralisedHyperbolicSeed[lambda_, alpha_, beta_, mu_, delta_][volume_] :=
    HyperbolicDistribution[lambda, alpha, beta, delta * volume, mu * volume];

End[] (* End Private Context *)

EndPackage[]