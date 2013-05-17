(* Mathematica Package *)

BeginPackage["AmbitStochastics`GeneralisedHyperbolics`"]
(* Exported symbols added here with SymbolName::usage *)  

GHDistribution::usage = "GHDistribution[{lambda, alpha, beta, delta, mu}, parametrization] gives a convenient interface to four different parameterizations of the univariate GH distribution. The parametrization argument is optional and can be set to \"Native\" (default), \"ScaleInvariant\", \"ShapeTriangle\", and \"Unbounded\"."
ConvertGHParameters::usage = "ConvertGHParameters[parameters, paramtrization1 -> parametrization2] converts between the parametrizations \"Native\", \"ScaleInvariant\", \"ShapeTriangle\", and \"Unbounded\"."

GHMomentEstimate::usage = "GHMomentEstimate[{x, w}, lambda, parametrization] attempts to find a moment estimate of the GH parameters"
SymmetricGHMomentEstimate::usage = "SymmetricGHMomentEstimate[{x, w}, lambda, parametrization] attempts to find moment a estimate of the GH parameters under the assumption of symmetry."
GHMomentEstimateGivenShape::usage = "GHMomentEstimateGivenShape[{x, w}, lambda, steepness, asymmetry, parametrization] finds a moment estimate of the GH parameters given fixed steepness and asymmetry. Other possibilities are GHMomentEstimateGivenShape[{x, w}, lambda, shape, parametrization] where shape either {xi, chi} or a list of such pairs."

GHShapeTriangleGrid::usage = "GHShapeTriangleGrid[n] returns a uniform n-by-n grid in the shape triangle. GHShapeTriangleGrid[n, {xi0, xi1}, {chi0, chi1}] retains only the shapes with xi0 <= steepness <= xi1 and chi0 <= asymmetry <= chi1."

GHLogLikelihood::usage = "GHLogLikelihood[{x, w}, parameters, parametrization] calculates the normalized log likelihood where parameters is either a single list of the five GH parameters, or a list of such lists."
GHProfileMLE::usage = "GHProfileMLE[{x, w}, lambda, paramterization, truncation, start, maxMGFOrder] attemps to perform MLE of the four GH parameters when the index lambda is fixed. The default parametrization is \"Native\", the truncation can be either \"Untruncated\" (default) or \"Truncated\", and start can either be Automatic (default) or a set of initial points for the Nelder-Mead algorithm. maxMGFOrder specifies the extend to which the moment generating function should exist."

GHLeastSquaresLogPDF::usage = "GHLeastSquaresLogPDF[{x, pdf}, lambda, start, parametrization] attempts to find a least squares fit of the log pdf using start as a starting point."

Begin["`Private`"] (* Begin Private Context *) 

Clear[GHDistribution]
GHDistribution[parameters_, parametrization_ : "Native"] := 
    HyperbolicDistribution @@ ConvertGHParameters[parameters, parametrization -> "Native"];

(* Convert between the four different parametrisations *)

Clear[ConvertGHParameters]
ConvertGHParameters[{lambda_, alpha_, beta_, delta_, mu_}, "Native" -> "ScaleInvariant"] :=
    {
        lambda, 
        alpha * delta, 
        beta * delta, 
        delta, 
        mu
    };
        
ConvertGHParameters[{lambda_, alpha_, beta_, delta_, mu_}, "Native" -> "ShapeTriangle"] :=
    {
        lambda, 
        1 / Sqrt[1 + delta * Sqrt[alpha ^ 2 - beta ^ 2]], 
        (beta / alpha) / Sqrt[1 + delta * Sqrt[alpha ^ 2 - beta ^ 2]],
        delta,
        mu
    };

ConvertGHParameters[{lambda_, alpha_, beta_, delta_, mu_}, "Native" -> "Unbounded"] :=
    {
        lambda, 
        Log[alpha * delta], 
        Tan[Pi / 2 * beta / alpha], 
        Log[delta], 
        mu
    };

ConvertGHParameters[{lambda_, alphaBar_, betaBar_, delta_, mu_}, "ScaleInvariant" -> "Native"] :=
    {
        lambda, 
        alphaBar / delta, 
        betaBar / delta, 
        delta, 
        mu
    };

ConvertGHParameters[{lambda_, alphaBar_, betaBar_, delta_, mu_}, "ScaleInvariant" -> "ShapeTriangle"] :=
    {
        lambda, 
        1 / Sqrt[1 + Sqrt[alphaBar ^ 2 - betaBar ^ 2]], 
        betaBar / (alphaBar Sqrt[1 + Sqrt[alphaBar ^ 2 - betaBar ^ 2]]),
        delta, 
        mu
    };
 
ConvertGHParameters[{lambda_, alphaBar_, betaBar_, delta_, mu_}, "ScaleInvariant" -> "Unbounded"] :=
    {
        lambda, 
        Log[alphaBar], 
        Tan[(betaBar Pi)/(2 alphaBar)],
        Log[delta], 
        mu
    };

ConvertGHParameters[{lambda_, xi_, chi_, delta_, mu_}, "ShapeTriangle" -> "Native"] :=
    {
        lambda,
        (1 - xi ^ 2) / (delta xi Sqrt[xi ^ 2 - chi ^ 2]),
        chi (1 - xi ^ 2) / (delta xi ^ 2 Sqrt[xi ^ 2 - chi ^ 2]),
        delta,
        mu
    };
    
ConvertGHParameters[{lambda_, xi_, chi_, delta_, mu_}, "ShapeTriangle" -> "ScaleInvariant"] :=
    {
        lambda, 
        (1 - xi ^ 2) / (xi Sqrt[xi ^ 2 - chi ^ 2]), 
        chi (1 - xi ^ 2) / (xi ^ 2 Sqrt[xi ^ 2 - chi ^ 2]),
        delta, 
        mu
    };
 
ConvertGHParameters[{lambda_, xi_, chi_, delta_, mu_}, "ShapeTriangle" -> "Unbounded"] :=
    {
        lambda, 
        Log[(1 - xi ^ 2) / (xi Sqrt[xi ^ 2 - chi ^ 2])], 
        Tan[(chi Pi) / (2 xi)], 
        Log[delta], 
        mu
    };

ConvertGHParameters[{lambda_, a_, b_, d_, mu_}, "Unbounded" -> "Native"] :=
    {
        lambda,
        Exp[a - d],
        2 Exp[a - d] ArcTan[b] / Pi,
        Exp[d],
        mu
    };

ConvertGHParameters[{lambda_, a_, b_, d_, mu_}, "Unbounded" -> "ScaleInvariant"] :=
    {
        lambda,
        Exp[a],
        2 Exp[a] ArcTan[b] / Pi,
        Exp[d],
        mu
    };
    
ConvertGHParameters[{lambda_, a_, b_, d_, mu_}, "Unbounded" -> "ShapeTriangle"] := 
    With[{xi = 1 / Sqrt[1 + Exp[a] Sqrt[1 - 4 ArcTan[b] ^ 2 / Pi ^ 2]]},
        {
            lambda,
            xi,
            xi 2 ArcTan[b] / Pi,
            Exp[d],
            mu
        }
    ];

ConvertGHParameters[pp_, HoldPattern[p_ -> p_]] := pp;


(* Maximum likelihood estimation *)

Untangle[n_, x_]:= 
    Which[
        MatchQ[Dimensions[x], {n, _}], x,
        MatchQ[Dimensions[x], {_, n}], Transpose[x],
        True,                          Return[$Failed]
    ];

Clear[GHLogLikelihood]
GHLogLikelihood[xw_, parameters_ ? VectorQ, parametrization_ : "Native"] := 
    Module[ {x, w},
        {x, w} = Untangle[2, xw];
        w /= Total[w];
        GHLogLikelihoodScaleInvariant[x, w, Sequence @@ ConvertGHParameters[parameters, parametrization -> "ScaleInvariant"], 0.0]
    ];
    
GHLogLikelihood[xw_, parameters_ ? (ArrayQ[#, 2] &), parametrization_ : "Native"] :=
    Module[ {x, w},
        {x, w} = Untangle[2, xw];
        w /= Total[w];
        GHLogLikelihoodScaleInvariant[x, w, Sequence @@ ConvertGHParameters[#, parametrization -> "ScaleInvariant"], 0.0] & /@ parameters
    ];

(* Assumes that Total[w] == 1 *)
Clear[GHLogLikelihoodScaleInvariant]
GHLogLikelihoodScaleInvariant =
    Compile[ {{x, _Real, 1}, {w, _Real, 1}, lambda, alpha, beta, delta, mu, maxMGFOrder},
        Module[ {gamma, z, q, llh},
            llh = -$MaxMachineNumber;
            If[ alpha - maxMGFOrder * delta > beta > - alpha && alpha > maxMGFOrder * delta / 2 && delta > 0,
                (* Valid parameters *)
                gamma = Sqrt[alpha ^ 2 - beta ^ 2];
                z     = (x - mu) / delta;
                q     = Sqrt[1 + z^2];
                llh   = (lambda Log[gamma] - (lambda - 0.5) Log[alpha] - Log[delta] - 0.9189385332046727 (* Log[Sqrt[2 Pi]] *) - Log @ BesselK[lambda, gamma]) 
                      + Dot[w, beta z + Log @ BesselK[lambda - 0.5, alpha q] + (lambda - 0.5) Log[q]];
            ];
            llh
        ]
    ];

(* Assumes that Total[w] == 1 *)
(* Assumed uniformly spaced data *)
Clear[TruncatedGHLogLikelihoodScaleInvariant]
TruncatedGHLogLikelihoodScaleInvariant =
    Compile[ {{x, _Real, 1}, {w, _Real, 1}, lambda, alpha, beta, delta, mu, dx, maxMGFOrder},
        Module[ {gamma, z, q, logk, logkk, preFactor, pdf, llh},
            llh = -$MaxMachineNumber;
            If[ alpha - maxMGFOrder > beta > - alpha && alpha > maxMGFOrder / 2 && delta > 0,
                (* Valid parameters *)
                gamma = Sqrt[alpha ^ 2 - beta ^ 2];
                z     = (x - mu) / delta;
                q     = Sqrt[1 + z^2];
                logk  = Log @ BesselK[lambda, gamma];
                logkk = Log @ BesselK[lambda - 0.5, alpha q];
                If[ Log[$MinMachineNumber] <= Min[logk, logkk] && Max[logk, logkk] <= Log[$MaxMachineNumber],
                    (* Valid values of BesselK *) 
                    preFactor = gamma ^ lambda / (alpha ^ (lambda - 0.5) delta 2.5066282746310002 (* Sqrt[2 Pi] *) Exp[logk]); 
                    pdf       = preFactor Exp[beta z + logkk] q ^ (lambda - 0.5);
                    If[ Min[pdf] > 0.0,
                        (* Valid pdf values *) 
                        llh = Dot[w, Log[pdf]] - Log[dx Total[pdf]];
                    ];
                ];
            ];
            llh
        ]
    ];

Clear[GHProfileMLE]
GHProfileMLE[
    xw_ ? (ArrayDepth[#] == 2 &), 
    lambda_ ? NumberQ, 
    parametrization_ : "Native", 
    truncation : ("Untruncated" | "Truncated") : "Untruncated", 
    start_ : Automatic,
    maxMGFOrder_ : 0.0
] :=
    Module[ {x, w, startScInv, res, obj, dx, shape, me, llh, best},
        {x, w} = Untangle[2, xw]; 
        w /= Total[w];
        If[ start === Automatic
            ,
            shape = GHShapeTriangleGrid[20];
            me    = GHMomentEstimateGivenShape[{x, w}, lambda, shape, "ScaleInvariant"];
            llh   = GHLogLikelihood[{x, w}, me, "ScaleInvariant"];
            best  = Ordering[llh, -20];
            startScInv = me[[best, 2 ;; -1]];
            ,
            startScInv = Rest[ConvertGHParameters[If[Length[#] == 4, Prepend[#, lambda], #], parametrization -> "ScaleInvariant"]] & /@ If[ArrayDepth[start] == 1, {start}, start];
        ];
        If[ truncation === "Truncated"
            ,
            dx = Mean[Differences[Sort[x]]];
            obj[a_ ? NumberQ, b_ ? NumberQ, d_ ? NumberQ, m_ ? NumberQ] := obj[a, b, d, m] = TruncatedGHLogLikelihoodScaleInvariant[x, w, lambda, a, b, d, m, dx, maxMGFOrder]
            ,
            obj[a_ ? NumberQ, b_ ? NumberQ, d_ ? NumberQ, m_ ? NumberQ] := obj[a, b, d, m] = GHLogLikelihoodScaleInvariant[x, w, lambda, a, b, d, m, maxMGFOrder];
        ];
        res = NMaximize[
            {obj[alpha, beta, delta, mu], alpha > Abs[beta] >= 0 && delta > 0},
            {alpha, beta, delta, mu},
            Method -> {"NelderMead", "PostProcess" -> False, "InitialPoints" -> startScInv},
            MaxIterations -> 1000
        ];
        {
            res[[1]], 
            ConvertGHParameters[{lambda, alpha, beta, delta, mu} /. res[[2]], "ScaleInvariant" -> parametrization]
        }
    ];

(* Simple grid in the shape cone *)

Clear[GHShapeTriangleGrid]
GHShapeTriangleGrid[n_Integer?Positive, steepnessRange_ : {0, 1}, asymmetryRange_ : {-1, 1}] :=
    Select[
        With[{d = 1 / (n + 1)}, N @ Flatten[Table[{xi, chi}, {xi, d, 1 - d, d}, {chi, -xi + d, xi - d, d}], 1]],
        steepnessRange[[1]] <= #[[1]] <= steepnessRange[[2]] && asymmetryRange[[1]] <= #[[2]] <= asymmetryRange[[2]] &
    ];

    
Clear[GHMomentEstimateGivenShape]
GHMomentEstimateGivenShape[xw_ ? (ArrayDepth[#] == 2 &), lambda_, xi_ ? NumberQ, chi_ ? NumberQ, param_ : "Native"] :=
    GHMomentEstimateGivenShape[xw, lambda, {xi, chi}, param];
    
GHMomentEstimateGivenShape[xw_ ? (ArrayDepth[#] == 2 &), lambda_, shape_ ? (ArrayQ[#, 1 | 2, NumberQ] &), param_ : "Native"] := 
    Module[
        {x, w, mean, var, z, r, k0, k1, k2, kappa1, kappa2, delta, mu, xi, chi},
        {x, w} = Untangle[2, xw]; 
        {xi, chi} = If[ArrayDepth[shape] == 1, shape, Untangle[2, shape]];
        w /= Total[w];
        mean = Dot[x, w];
        var  = Dot[(x - mean) ^ 2, w];
        z    = xi ^ -2 - 1;
        r    = (chi / xi) ^ 2 * (1 - xi ^ 2) / (xi ^ 2 - chi ^ 2);
        k0   = BesselK[lambda, z];
        k1   = BesselK[lambda + 1, z];
        k2   = BesselK[lambda + 2, z];
        kappa1 = k1 / k0;
        kappa2 = k2 / k0;
        delta = Sqrt[var z / (kappa1 - r (kappa1 ^ 2 - kappa2))];
        mu    = mean - delta kappa1 chi / Sqrt[xi ^ 2 - chi ^ 2];
        If[ ArrayDepth[shape] == 1,
            ConvertGHParameters[{lambda, xi, chi, delta, mu}, "ShapeTriangle" -> param],
            MapThread[ConvertGHParameters[{lambda, ##}, "ShapeTriangle" -> param] &, {xi, chi, delta, mu}] 
        ]
    ];


(* Moment based parameters estimates *)

Clear[GHMomentEstimate]
GHMomentEstimate[xw_ ? (ArrayDepth[#] == 2 &), lambda_, parametrization_ : "Native"] := 
    Module[ {x, w, k, mean, cm2, cm3, cm4, emean, ecm2, ecm3, ecm4, res, obj},
        {x, w} = Untangle[2, xw]; 
        w /= Total[w];
        emean = Dot[x, w];
        ecm2  = Dot[(x - emean)^2, w];
        ecm3  = Dot[(x - emean)^3, w];
        ecm4  = Dot[(x - emean)^4, w];
        k[gamma_] := BesselK[lambda + 1, gamma]/BesselK[lambda, gamma];
        (* (Central) moments in the scale invariant parametrization *)
        mean[beta_, gamma_, delta_, mu_] := 
            mu + (beta delta k[gamma])/gamma;
        cm2[beta_,  gamma_, delta_] :=
            (beta^2 delta^2)/gamma^2 + (delta^2 (2 beta^2 + gamma^2 + 2 beta^2 lambda) k[gamma])/gamma^3 - (beta^2 delta^2 k[gamma]^2)/gamma^2;
        cm3[beta_,  gamma_, delta_] :=
            (beta*delta^3*(4*beta^2 + 3*gamma^2 + 2*beta^2*lambda))/gamma^4 - (2*beta*delta^3*(-4*beta^2 - 3*gamma^2 + beta^2*gamma^2 - 6*beta^2*lambda - 3*gamma^2*lambda - 2*beta^2*lambda^2)*k[gamma])/gamma^5 - (3*beta*delta^3*(2*beta^2 + gamma^2 + 2*beta^2*lambda)*k[gamma]^2)/gamma^4 + (2*beta^3*delta^3*k[gamma]^3)/gamma^3;
        cm4[beta_,  gamma_, delta_] :=
            (delta^4*(24*beta^4 + 24*beta^2*gamma^2 + beta^4*gamma^2 + 3*gamma^4 + 20*beta^4*lambda + 12*beta^2*gamma^2*lambda + 4*beta^4*lambda^2))/gamma^6 - (2*delta^4*(-24*beta^4 - 24*beta^2*gamma^2 + 4*beta^4*gamma^2 - 3*gamma^4 + 3*beta^2*gamma^4 - 44*beta^4*lambda - 36*beta^2*gamma^2*lambda + 2*beta^4*gamma^2*lambda - 3*gamma^4*lambda - 24*beta^4*lambda^2 - 12*beta^2*gamma^2*lambda^2 - 4*beta^4*lambda^3)*k[gamma])/gamma^7 + (2*beta^2*delta^4*(-16*beta^2 - 12*gamma^2 + beta^2*gamma^2 - 24*beta^2*lambda - 12*gamma^2*lambda - 8*beta^2*lambda^2)*k[gamma]^2)/gamma^6 + (6*beta^2*delta^4*(2*beta^2 + gamma^2 + 2*beta^2*lambda)*k[gamma]^3)/gamma^5 - (3*beta^4*delta^4*k[gamma]^4)/gamma^4;
        obj[b_ ? NumberQ, g_ ? NumberQ, d_ ? NumberQ, m_ ? NumberQ] := obj[b, g, d, m] = 
            If[ gamma <= 0 || delta <= 0, 
                $MaxMachineNumber, 
               (emean - mean[beta, gamma, delta, mu])^2
                + (ecm2 - cm2[beta, gamma, delta])^2
                + (ecm3 - cm3[beta, gamma, delta])^2
                + (ecm4 - cm4[beta, gamma, delta])^2
            ];
        res = NMinimize[{obj[beta, gamma, delta, mu], gamma > 0 && delta > 0}, {beta, gamma, delta, mu}];
        ConvertGHParameters[{lambda, Sqrt[gamma ^ 2 + beta ^ 2], beta, delta, mu} /. res[[2]], "ScaleInvariant" -> parametrization]
    ];
    

Clear[SymmetricGHMomentEstimate]
SymmetricGHMomentEstimate[xw_ ? (ArrayDepth[#] == 2 &), lambda_, parametrization_ : "Native"] := 
    Module[ {x, w, mean, var, kurt, minKurt, maxKurt, a, aa, dd, mm, res, obj},
        {x, w} = Untangle[2, xw]; 
        w /= Total[w];
        mean = Dot[x, w]; 
        var  = Dot[(x - mean) ^ 2, w]; 
        kurt = Dot[(x - mean) ^ 4, w] / var ^ 2; 
        minKurt = 3.0; 
        maxKurt = Which[
            lambda < -2, 3*(1 + lambda)/(2 + lambda),
            lambda >  0, 3 + 3/lambda,
            True       , Infinity
        ];
        obj[aaa_ ? NumberQ] := obj[aaa] = (3*BesselK[lambda, E^aaa]*BesselK[2 + lambda, E^aaa])/BesselK[1 + lambda, E^aaa]^2 - kurt;
        res = If[
            minKurt < kurt < maxKurt
            , 
            mm = mean; 
            aa = FindRoot[obj[a], {a, 0.}, PrecisionGoal -> 10][[1,2]]; 
            dd = (1/2)*(aa + Log[(var*BesselK[lambda, E^aa])/BesselK[1 + lambda, E^aa]]);
            {lambda, aa, 0., dd, mm}
            , 
            {lambda, 0., 0., 0., 0.}
        ];
        ConvertGHParameters[res, "Unbounded" -> parametrization]
    ];

(* Least squares fit to log pdf *)

Clear[GHLeastSquaresLogPDF]
GHLeastSquaresLogPDF[xpdf_ ? (ArrayDepth[#] == 2 &), lambda_, start_, parametrization_ : "Native"] :=
    Module[ {x, pdf, startScInv, obj, res},
        {x, pdf} = Untangle[2, xpdf];
        {x, pdf} = Transpose @ Select[Transpose @ {x, pdf}, Last[#] > 0 &];
        startScInv = Rest[ConvertGHParameters[If[Length[#] == 4, Prepend[#, lambda], #], parametrization -> "ScaleInvariant"]] & /@ If[ArrayDepth[start] == 1, {start}, start];
        obj[a_ ? NumberQ, b_ ? NumberQ, d_ ? NumberQ, m_ ? NumberQ] := obj[a, b, d, m] = GHLogPDFSquaredErrorScaleInvariant[x, pdf, lambda, a, b, d, m];
        res = NMinimize[
            {obj[alpha, beta, delta, mu], alpha > Abs[beta] >= 0 && delta > 0},
            {alpha, beta, delta, mu},
            Method -> {"NelderMead", "PostProcess" -> False, "InitialPoints" -> startScInv},
            MaxIterations -> 1000
        ];
        {
            res[[1]], 
            ConvertGHParameters[{lambda, alpha, beta, delta, mu} /. res[[2]], "ScaleInvariant" -> parametrization]
        }
    ];

Clear[GHLogPDFSquaredErrorScaleInvariant]
GHLogPDFSquaredErrorScaleInvariant =
    Compile[ {{x, _Real, 1}, {pdf, _Real, 1}, lambda, alpha, beta, delta, mu},
        Module[ {gamma, z, q, logpdf, logk, logkk, err},
            err = $MaxMachineNumber;
            If[ alpha > Abs[beta] >= 0 && delta > 0,
                (* Valid parameters *)
                gamma = Sqrt[alpha ^ 2 - beta ^ 2];
                z     = (x - mu) / delta;
                q     = Sqrt[1 + z^2];
                logk  = Log @ BesselK[lambda, gamma];
                logkk = Log @ BesselK[lambda - 0.5, alpha q];
                If[ Log[$MinMachineNumber] <= Min[logk, logkk] && Max[logk, logkk] <= Log[$MaxMachineNumber],
                    (* Valid values of BesselK *) 
                    logpdf = (lambda Log[gamma] - (lambda - 0.5) Log[alpha] - Log[delta] - 0.9189385332046727 (* Log[Sqrt[2 Pi]] *) - logk) + beta z + logkk + (lambda - 0.5) Log[q];
                    err = Total[(logpdf - Log[pdf]) ^ 2];
                ];
            ];
            err
        ]
    ];


End[] (* End Private Context *)

EndPackage[]