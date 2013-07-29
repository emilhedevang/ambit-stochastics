(* Mathematica Package *)

BeginPackage["AmbitStochastics`SpatialRandomVectorFields`"]
(* Exported symbols added here with SymbolName::usage *)  

OddKernelBoxIntegrals::usage       = "OddKernelBoxIntegrals[f, n, delta]."
OddKernelBoxIntegrals2::usage       = "OddKernelBoxIntegrals2[f, n, delta]."
OddKernelBoxCovarianceRoots::usage = "OddKernelBoxCovarianceRoots[W11, W12, cov]."
OddKernelBoxCovarianceRoots2::usage = "OddKernelBoxCovarianceRoots2[W11, W12, cov]."
SimulateSpatialRandomVectorField3D::usage = "SimulateSpatialRandomVectorField3D[boxCovRoot, {n1, n2, n3}]."

Begin["`Private`"] (* Begin Private Context *) 

$CompilationTarget = "WVM";

Clear[BoxIntegralApproximation]
BoxIntegralApproximation[f_, beta_, delta_] :=
    delta^3 * Dot[f[delta * (# + beta)] & /@ BoxIntegralApproximationOrbits, BoxIntegralApproximationWeights];

Clear[BoxIntegralApproximator]
BoxIntegralApproximator[f_, delta_] :=
    Compile[{{beta, _Real, 1}},
        delta^3 * Dot[f[delta * (# + beta)] & /@ BoxIntegralApproximationOrbits, BoxIntegralApproximationWeights],
        CompilationTarget -> $CompilationTarget
    ];


Clear[OddKernelBoxIntegrals]
OddKernelBoxIntegrals[f_, n_, delta_, tbl : (Table | ParallelTable) : Table] :=
    Module[ {W11, W12},
        W11 = Developer`ToPackedArray @ N @ Chop[#, $MinMachineNumber] & @ tbl[BoxIntegralApproximation[f[Sqrt @ Total[# * #]] * #[[1]] * #[[1]] &, {b1, b2, b3}, delta], {b1, 0, n}, {b2, 0, n}, {b3, 0, n}];
        W12 = Developer`ToPackedArray @ N @ Chop[#, $MinMachineNumber] & @ tbl[BoxIntegralApproximation[f[Sqrt @ Total[# * #]] * #[[1]] * #[[2]] &, {b1, b2, b3}, delta], {b1, 0, n}, {b2, 0, n}, {b3, 0, n}];
        {W11, W12}
    ];

Clear[OddKernelBoxIntegrals2]
OddKernelBoxIntegrals2[f_, n_, delta_, tbl : (Table | ParallelTable) : Table] :=
    Module[ {W11, W12},
        W11 = Developer`ToPackedArray @ N @ Chop[#, $MinMachineNumber] & @ 
            tbl[If[b3 <= b2, BoxIntegralApproximation[f[Sqrt @ Total[# * #]] * #[[1]] * #[[1]] &, {b1, b2, b3}, delta], 0.0], {b1, 0, n}, {b2, 0, n}, {b3, 0, n}];
        W12 = Developer`ToPackedArray @ N @ Chop[#, $MinMachineNumber] & @ 
            tbl[If[b2 <= b1, BoxIntegralApproximation[f[Sqrt @ Total[# * #]] * #[[1]] * #[[2]] &, {b1, b2, b3}, delta], 0.0], {b1, 0, n}, {b2, 0, n}, {b3, 0, n}];
        {W11, W12}
    ];

OddKernelBoxIntegrals /: Parallelize[OddKernelBoxIntegrals[f_, n_, delta_]] := 
    OddKernelBoxIntegrals[f, n, delta, ParallelTable];

Clear[OddKernelBoxCovariance]
OddKernelBoxCovariance = 
    Compile[{w11, w12, w13, w22, w23, w33, cov11, cov12, cov13, cov22, cov23, cov33},
        {
            {
                cov33 w22 - 2 cov23 w23 + cov22 w33, 
                -cov33 w12 + cov23 w13 + cov13 w23 - cov12 w33, 
                cov23 w12 - cov22 w13 - cov13 w22 + cov12 w23
            },
            {
                -cov33 w12 + cov23 w13 + cov13 w23 - cov12 w33, 
                cov33 w11 - 2 cov13 w13 + cov11 w33, 
                -cov23 w11 + cov13 w12 + cov12 w13 - cov11 w23
            }, 
            {
                cov23 w12 - cov22 w13 - cov13 w22 + cov12 w23, 
                -cov23 w11 + cov13 w12 + cov12 w13 - cov11 w23, 
                cov22 w11 - 2 cov12 w12 + cov11 w22
            }
        },
        CompilationTarget -> $CompilationTarget
    ];


Clear[OddKernelBoxCovarianceRoots]
OddKernelBoxCovarianceRoots[W11_, W12_, cov_, tbl : (Table | ParallelTable) : Table] :=
    With[ {n = Length[W11] - 1},
        Developer`ToPackedArray @ tbl[
            Re @ MatrixFunction[
                Sqrt,
                OddKernelBoxCovariance[
                    W11[[Abs[b1]+1,Abs[b2]+1,Abs[b3]+1]],
                    Sign[b1] Sign[b2] W12[[Abs[b1]+1,Abs[b2]+1,Abs[b3]+1]],
                    Sign[b1] Sign[b3] W12[[Abs[b1]+1,Abs[b3]+1,Abs[b2]+1]],             
                    W11[[Abs[b2]+1,Abs[b1]+1,Abs[b3]+1]],
                    Sign[b2] Sign[b3] W12[[Abs[b2]+1,Abs[b3]+1,Abs[b1]+1]],
                    W11[[Abs[b3]+1,Abs[b2]+1,Abs[b1]+1]],
                    cov[[1,1]], cov[[1,2]], cov[[1,3]], cov[[2,2]], cov[[2,3]], cov[[3,3]]
                ]
            ],
            {b1, -n, n}, {b2, -n, n}, {b3, -n, n}
        ]
    ];

Clear[OddKernelBoxCovarianceRoots2]
OddKernelBoxCovarianceRoots2[W11_, W12_, cov_, tbl : (Table | ParallelTable) : Table] :=
    With[ {n = Length[W11] - 1},
        Developer`ToPackedArray @ tbl[
            Re @ MatrixFunction[
                Sqrt,
                OddKernelBoxCovariance[
                    W11[[Abs[b1]+1,Max@Abs[{b2,b3}]+1,Min@Abs[{b2,b3}]+1]],
                    Sign[b1] Sign[b2] W12[[Max@Abs[{b1,b2}]+1,Min@Abs[{b1,b2}]+1,Abs[b3]+1]],
                    Sign[b1] Sign[b3] W12[[Max@Abs[{b1,b3}]+1,Min@Abs[{b1,b3}]+1,Abs[b2]+1]],             
                    W11[[Abs[b2]+1,Max@Abs[{b1,b3}]+1,Min@Abs[{b1,b3}]+1]],
                    Sign[b2] Sign[b3] W12[[Max@Abs[{b2,b3}]+1,Min@Abs[{b2,b3}]+1,Abs[b1]+1]],
                    W11[[Abs[b3]+1,Max@Abs[{b1,b2}]+1,Min@Abs[{b1,b2}]+1]],
                    cov[[1,1]], cov[[1,2]], cov[[1,3]], cov[[2,2]], cov[[2,3]], cov[[3,3]]
                ]
            ],
            {b1, -n, n}, {b2, -n, n}, {b3, -n, n}
        ]
    ];


SimulateSpatialRandomVectorField3D[boxCovRoot_, {n1_, n2_, n3_}, mode : ("Linear" | "Cyclic") : "Linear"] :=
    Module[ {dim, k, x11, x12, x13, x21, x22, x23, x31, x32, x33, z1, z2, z3},
        If[ mode == "Linear",
            dim =  {n1, n2, n3} + Dimensions[boxCovRoot, 3] - 1;
            k = {-1, 1};,
            dim =  {n1, n2, n3};
            k = {1, 1};
        ];
        PrintTemporary["z1"];
        z1  = RandomVariate[NormalDistribution[], dim];
        PrintTemporary["z2"];
        z2  = RandomVariate[NormalDistribution[], dim];
        PrintTemporary["z3"];
        z3  = RandomVariate[NormalDistribution[], dim];
        PrintTemporary["x11"];
        x11 = ListConvolve[boxCovRoot[[All, All, All, 1, 1]], z1, k];
        PrintTemporary["x12"];
        x12 = ListConvolve[boxCovRoot[[All, All, All, 1, 2]], z2, k];
        PrintTemporary["x13"];
        x13 = ListConvolve[boxCovRoot[[All, All, All, 1, 3]], z3, k];
        PrintTemporary["x21"];
        x21 = ListConvolve[boxCovRoot[[All, All, All, 2, 1]], z1, k];
        PrintTemporary["x22"];
        x22 = ListConvolve[boxCovRoot[[All, All, All, 2, 2]], z2, k];
        PrintTemporary["x23"];
        x23 = ListConvolve[boxCovRoot[[All, All, All, 2, 3]], z3, k];
        PrintTemporary["x31"];
        x31 = ListConvolve[boxCovRoot[[All, All, All, 3, 1]], z1, k];
        PrintTemporary["x32"];
        x32 = ListConvolve[boxCovRoot[[All, All, All, 3, 2]], z2, k];
        PrintTemporary["x33"];
        x33 = ListConvolve[boxCovRoot[[All, All, All, 3, 3]], z3, k];
        {
            x11 + x12 + x13, 
            x21 + x22 + x23,
            x31 + x32 + x33
        }
    ];    





BoxIntegralApproximationOrbits = 
    Developer`ToPackedArray @ {{0., 0., 0.}, {-0.47795365790226946, 0., 0.}, {0., -0.47795365790226946, 0.}, {0., 0., -0.47795365790226946}, {0., 0., 0.47795365790226946}, {0., 0.47795365790226946, 0.}, 
    {0.47795365790226946, 0., 0.}, {-0.20302858736911977, 0., 0.}, {0., -0.20302858736911977, 0.}, {0., 0., -0.20302858736911977}, {0., 0., 0.20302858736911977}, {0., 0.20302858736911977, 0.}, 
    {0.20302858736911977, 0., 0.}, {-0.4476273546261776, 0., 0.}, {0., -0.4476273546261776, 0.}, {0., 0., -0.4476273546261776}, {0., 0., 0.4476273546261776}, {0., 0.4476273546261776, 0.}, 
    {0.4476273546261776, 0., 0.}, {-0.125, 0., 0.}, {0., -0.125, 0.}, {0., 0., -0.125}, {0., 0., 0.125}, {0., 0.125, 0.}, {0.125, 0., 0.}, {-0.47795365790226946, -0.47795365790226946, 0.}, 
    {-0.47795365790226946, 0., -0.47795365790226946}, {-0.47795365790226946, 0., 0.47795365790226946}, {-0.47795365790226946, 0.47795365790226946, 0.}, 
    {0., -0.47795365790226946, -0.47795365790226946}, {0., -0.47795365790226946, 0.47795365790226946}, {0., 0.47795365790226946, -0.47795365790226946}, {0., 0.47795365790226946, 0.47795365790226946}, 
    {0.47795365790226946, -0.47795365790226946, 0.}, {0.47795365790226946, 0., -0.47795365790226946}, {0.47795365790226946, 0., 0.47795365790226946}, {0.47795365790226946, 0.47795365790226946, 0.}, 
    {-0.47795365790226946, -0.20302858736911977, 0.}, {-0.47795365790226946, 0., -0.20302858736911977}, {-0.47795365790226946, 0., 0.20302858736911977}, 
    {-0.47795365790226946, 0.20302858736911977, 0.}, {-0.20302858736911977, -0.47795365790226946, 0.}, {-0.20302858736911977, 0., -0.47795365790226946}, 
    {-0.20302858736911977, 0., 0.47795365790226946}, {-0.20302858736911977, 0.47795365790226946, 0.}, {0., -0.47795365790226946, -0.20302858736911977}, {0., -0.47795365790226946, 0.20302858736911977}, 
    {0., -0.20302858736911977, -0.47795365790226946}, {0., -0.20302858736911977, 0.47795365790226946}, {0., 0.20302858736911977, -0.47795365790226946}, {0., 0.20302858736911977, 0.47795365790226946}, 
    {0., 0.47795365790226946, -0.20302858736911977}, {0., 0.47795365790226946, 0.20302858736911977}, {0.20302858736911977, -0.47795365790226946, 0.}, {0.20302858736911977, 0., -0.47795365790226946}, 
    {0.20302858736911977, 0., 0.47795365790226946}, {0.20302858736911977, 0.47795365790226946, 0.}, {0.47795365790226946, -0.20302858736911977, 0.}, {0.47795365790226946, 0., -0.20302858736911977}, 
    {0.47795365790226946, 0., 0.20302858736911977}, {0.47795365790226946, 0.20302858736911977, 0.}, {-0.47795365790226946, -0.47795365790226946, -0.47795365790226946}, 
    {-0.47795365790226946, -0.47795365790226946, 0.47795365790226946}, {-0.47795365790226946, 0.47795365790226946, -0.47795365790226946}, 
    {-0.47795365790226946, 0.47795365790226946, 0.47795365790226946}, {0.47795365790226946, -0.47795365790226946, -0.47795365790226946}, 
    {0.47795365790226946, -0.47795365790226946, 0.47795365790226946}, {0.47795365790226946, 0.47795365790226946, -0.47795365790226946}, {0.47795365790226946, 0.47795365790226946, 0.47795365790226946}, 
    {-0.34303789878087815, -0.34303789878087815, -0.34303789878087815}, {-0.34303789878087815, -0.34303789878087815, 0.34303789878087815}, 
    {-0.34303789878087815, 0.34303789878087815, -0.34303789878087815}, {-0.34303789878087815, 0.34303789878087815, 0.34303789878087815}, 
    {0.34303789878087815, -0.34303789878087815, -0.34303789878087815}, {0.34303789878087815, -0.34303789878087815, 0.34303789878087815}, 
    {0.34303789878087815, 0.34303789878087815, -0.34303789878087815}, {0.34303789878087815, 0.34303789878087815, 0.34303789878087815}};
 
BoxIntegralApproximationWeights = 
    Developer`ToPackedArray @ {-0.2028842312372801, -0.08037737787038624, -0.08037737787038624, -0.08037737787038624, -0.08037737787038624, -0.08037737787038624, -0.08037737787038624, 0.0788999673604565, 0.0788999673604565, 
    0.0788999673604565, 0.0788999673604565, 0.0788999673604565, 0.0788999673604565, 0.057693384490972686, 0.057693384490972686, 0.057693384490972686, 0.057693384490972686, 0.057693384490972686, 
    0.057693384490972686, 0., 0., 0., 0., 0., 0., 0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 
    0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 0.004907147921572262, 0.022543144647178933, 0.022543144647178933, 
    0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 
    0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 
    0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.022543144647178933, 0.001770878225839135, 0.001770878225839135, 
    0.001770878225839135, 0.001770878225839135, 0.001770878225839135, 0.001770878225839135, 0.001770878225839135, 0.001770878225839135, 0.03143751436914348, 0.03143751436914348, 0.03143751436914348, 
    0.03143751436914348, 0.03143751436914348, 0.03143751436914348, 0.03143751436914348, 0.03143751436914348};
 


End[] (* End Private Context *)

EndPackage[]