(* Mathematica Package *)

BeginPackage["AmbitStochastics`Gamma`"]
(* Exported symbols added here with SymbolName::usage *)
  
GammaKernel::usage                         = "FIXME";
GammaKernelCorrelationFunction::usage      = "FIXME";
GammaKernelVariance::usage                 = "FIXME";
GammaKernelCovarianceFunction::usage       = "FIXME";
GammaKernelSpectralDensityFunction::usage  = "FIXME";

Begin["`Private`"] (* Begin Private Context *) 

Clear[GammaKernel]
GammaKernel[t_, nu_, lambda_] := 
    Boole[t > 0] * lambda ^ nu / Gamma[nu] * t ^ (nu - 1) * Exp[-lambda * t];

Clear[GammaKernelCorrelationFunction]
GammaKernelCorrelationFunction[t_, nu_, lambda_] /; t != 0 :=
    2 ^ (3/2 - nu) * (lambda * Abs[t]) ^ (nu - 1/2) * BesselK[nu - 1/2, lambda * Abs[t]] / Gamma[nu - 1/2]
GammaKernelCorrelationFunction[t_, nu_, lambda_] /; t == 0 :=
    1

Clear[GammaKernelVariance]
GammaKernelVariance[nu_, lambda_] := 
    2 ^ (1 - 2 nu) * lambda * Gamma[2 nu - 1] / Gamma[nu] ^ 2;

Clear[GammaKernelCovarianceFunction]
GammaKernelCovarianceFunction[t_, nu_, lambda_] :=
    GammaKernelVariance[nu, lambda] * GammaKernelCorrelationFunction[t, nu, lambda];

Clear[GammaKernelSpectralDensityFunction]
GammaKernelSpectralDensityFunction[omega_, nu_, lambda_] :=
    (1 + (2 * Pi * Abs[omega] / lambda) ^ 2) ^ (-nu);

End[] (* End Private Context *)

EndPackage[]