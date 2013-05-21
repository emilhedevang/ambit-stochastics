(* Mathematica Package *)

BeginPackage["AmbitStochastics`Gamma`", {"AmbitStochastics`Utilities`"}]
(* Exported symbols added here with SymbolName::usage *)
  
GammaKernel::usage                           = "FIXME";
GammaVariance::usage                         = "FIXME";
GammaCorrelationFunction::usage              = "FIXME";
GammaComplementaryCorrelationFunction::usage = "FIXME";
GammaCovarianceFunction::usage               = "FIXME";
GammaComplementaryCovarianceFunction         = "FIXME";
GammaSpectralDensityFunction::usage          = "FIXME";

GammaGammaKernel::usage                           = "FIXME";
GammaGammaVariance::usage                         = "FIXME";
GammaGammaCorrelationFunction::usage              = "FIXME";
GammaGammaComplementaryCorrelationFunction::usage = "FIXME";
GammaGammaCovarianceFunction::usage               = "FIXME";
GammaGammaComplementaryCovarianceFunction::usage  = "FIXME";
GammaGammaSpectralDensityFunction::usage          = "FIXME";

Begin["`Private`"] (* Begin Private Context *) 

(* Gamma *)

Clear[GammaKernel]
GammaKernel[t_, nu_, lambda_] :=
    Boole[t > 0] * lambda ^ nu / Gamma[nu] * t ^ (nu - 1) * Exp[-lambda * t];

Clear[GammaCorrelationFunction]
GammaCorrelationFunction[t_, nu_, lambda_] :=
    If[ t == 0,
        1,
        2 ^ (3/2 - nu) * (lambda * Abs[t]) ^ (nu - 1/2) * BesselK[nu - 1/2, lambda * Abs[t]] / Gamma[nu - 1/2]
    ];

Clear[GammaComplementaryCorrelationFunction]
GammaComplementaryCorrelationFunction[t_, nu_, lambda_] :=
    If[ t == 0,
        0,
        1 - 2 ^ (3/2 - nu) * (lambda * Abs[t]) ^ (nu - 1/2) * BesselK[nu - 1/2, lambda * Abs[t]] / Gamma[nu - 1/2]
    ];

    
Clear[GammaVariance]
GammaVariance[nu_, lambda_] :=
    2 ^ (1 - 2 nu) * lambda * Gamma[2 nu - 1] / Gamma[nu] ^ 2;

Clear[GammaCovarianceFunction]
GammaCovarianceFunction[t_, nu_, lambda_] :=
    GammaVariance[nu, lambda] * GammaCorrelationFunction[t, nu, lambda];

Clear[GammaComplementaryCovarianceFunction]
GammaComplementaryCovarianceFunction[t_, nu_, lambda_] :=
    GammaVariance[nu, lambda] * GammaComplementaryCorrelationFunction[t, nu, lambda];

Clear[GammaSpectralDensityFunction]
GammaSpectralDensityFunction[omega_, nu_, lambda_] :=
    (1 + (2 * Pi * Abs[omega] / lambda) ^ 2) ^ (-nu);


(* Gamma Gamma *)

Clear[GammaGammaKernel]
GammaGammaKernel[t_, nu1_, nu2_, lambda1_, lambda2_] :=
    Boole[t > 0] * lambda1 ^ nu1 * lambda2 ^ nu2 / Gamma[nu1 + nu2] * t ^ (nu1 + nu2 - 1) * Exp[-lambda1 * t] * Hypergeometric1F1[nu2, nu1 + nu2, -(lambda2 - lambda1) * t];
    
Clear[GammaGammaVariance]
GammaGammaVariance[nu1_, nu2_, lambda1_, lambda2_] :=
    lambda1 / Pi * NIntegrate[(1+u^2)^-nu1*(1+(u/(lambda2/lambda1))^2)^-nu2, {u, 0, Infinity}];

Clear[GammaGammaSpectralDensityFunction]
GammaGammaSpectralDensityFunction[omega_, nu1_, nu2_, lambda1_, lambda2_] :=
    (1 + (2 * Pi * Abs[omega] / lambda1) ^ 2) ^ (-nu1) * (1 + (2 * Pi * Abs[omega] / lambda2) ^ 2) ^ (-nu2);

Clear[GammaGammaCorrelationFunction]
GammaGammaCorrelationFunction[t_, nu1_, nu2_, lambda1_, lambda2_] :=
    GammaGammaCovarianceFunction[t, nu1, nu2, lambda1, lambda2] / GammaGammaVariance[nu1, nu2, lambda1, lambda2];

Clear[GammaGammaComplementaryCorrelationFunction]
GammaGammaComplementaryCorrelationFunction[t_, nu1_, nu2_, lambda1_, lambda2_] :=
    GammaGammaComplementaryCovarianceFunction[t, nu1, nu2, lambda1, lambda2] / GammaGammaVariance[nu1, nu2, lambda1, lambda2];

Clear[GammaGammaCovFunInt]
SetAttributes[GammaGammaCovFunInt, Listable]
GammaGammaCovFunInt[t_, nu_, mu_, kappa_] :=
    If[ Abs[t] > 1,
        NIntegrate[(1 + u^2)^-nu (1 + (u / kappa)^2)^-mu, {u, 0, Infinity}] - 
        NIntegrate[Null, {u, 0, Infinity},
            Method -> {
                "LevinRule",
                "AdditiveTerm" -> (1 + u^2)^-nu (1 + (u / kappa)^2)^-mu,
                "Amplitude" -> {-(1 + u^2)^-nu (1 + (u / kappa)^2)^-mu, 0},
                "Kernel" -> {Cos[t u], -Sin[t u]},
                "DifferentialMatrix" -> {{0, t}, {-t, 0}}
            },
            MaxRecursion -> 20
        ],
        NIntegrate[Null, {u, 0, Infinity},
            Method -> {
                "LevinRule",
                "AdditiveTerm" -> 0,
                "Amplitude" -> {(1 + u^2)^-nu (1 + (u / kappa)^2)^-mu, 0},
                "Kernel" -> {Cos[t u], -Sin[t u]},
                "DifferentialMatrix" -> {{0, t}, {-t, 0}}
            },
            MaxRecursion -> 20
        ]
    ];

Clear[GammaGammaCovarianceFunction]
GammaGammaCovarianceFunction[t_, nu1_, nu2_, lambda1_, lambda2_] :=
    lambda1 / Pi * GammaGammaCovFunInt[lambda1 * t, nu1, nu2, lambda2 / lambda1];

Clear[GammaGammaComplCovFunInt]
SetAttributes[GammaGammaComplCovFunInt, Listable]
GammaGammaComplCovFunInt[t_, nu_, mu_, kappa_] :=
    If[ Abs[t] > 1,
        NIntegrate[Null, {u, 0, Infinity},
            Method -> {
                "LevinRule",
                "AdditiveTerm" -> (1 + u^2)^-nu (1 + (u / kappa)^2)^-mu,
                "Amplitude" -> {-(1 + u^2)^-nu (1 + (u / kappa)^2)^-mu, 0},
                "Kernel" -> {Cos[t u], -Sin[t u]},
                "DifferentialMatrix" -> {{0, t}, {-t, 0}}
            },
            MaxRecursion -> 20
        ],
        NIntegrate[(1 + u^2)^-nu (1 + (u / kappa)^2)^-mu, {u, 0, Infinity}]-
        NIntegrate[Null, {u, 0, Infinity},
            Method -> {
                "LevinRule",
                "AdditiveTerm" -> 0,
                "Amplitude" -> {(1 + u^2)^-nu (1 + (u / kappa)^2)^-mu, 0},
                "Kernel" -> {Cos[t u], -Sin[t u]},
                "DifferentialMatrix" -> {{0, t}, {-t, 0}}
            },
            MaxRecursion -> 20
        ]
    ];

GammaGammaComplementaryCovarianceFunction[t_, nu1_, nu2_, lambda1_, lambda2_] :=
    lambda1 / Pi * GammaGammaComplCovFunInt[lambda1 * t, nu1, nu2, lambda2 / lambda1];
 
End[] (* End Private Context *)

EndPackage[]