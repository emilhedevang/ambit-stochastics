(* Mathematica Package *)

BeginPackage["AmbitStochastics`TrawlProcesses`", {"AmbitStochastics`Utilities`"}]
(* Exported symbols added here with SymbolName::usage *)  

CascadeTrawlKernel::usage   = "FIXME";
PrepareTrawlVolumes::usage  = "PrepareTrawlVolumes[kernel, upperBound, n] precomputes coefficients needed in the trawl simulation. The kernel is assumed to be decreasing and have support contained in [0, upperBound], and the discretisation mesh size is upperBound / n.";
SimulateTrawlProcess::usage = "SimulateTrawlProcess[volumes, seed, n] returns a simulation of length n of the trawl process corresponding to the precomputed volumes; seed specifies the Levy seed.";

Begin["`Private`"] (* Begin Private Context *) 

Clear[CascadeTrawlKernel]
CascadeTrawlKernel[t_, T_, L_, theta_] /; T > 0 && L > 1 && theta > 0 :=
    Boole[0 < t < T] * 1 / L * ((T ^ theta - t ^ theta) / ((T / L) ^ theta + t ^ theta)) ^ (1 / theta);
    
Clear[PrepareTrawlVolumes]
PrepareTrawlVolumes[kernel_, upperBound_, n_, tbl : (Table | ParallelTable) : Table] :=
    Append[Most @ # - Rest @ #, Last @ #] & @ tbl[2.0 * NIntegrate[kernel[s], {s, (i - 1) * upperBound / n, i * upperBound / n}], {i, n}];

PrepareTrawlVolumes /: Parallelize[PrepareTrawlVolumes[kernel_, upperBound_, n_]] :=
    PrepareTrawlVolumes[kernel, upperBound, n, ParallelTable];

Clear[MovingTotalC]
MovingTotalC =
    Compile[{{x, _Real, 1}, {n, _Integer, 0}},
        With[ {y = Prepend[Accumulate[x], 0.0]},
            Drop[y, n] - Drop[y, -n]
        ],
        CompilationTarget -> "C"
    ];

Clear[SimulateTrawlProcess]
SimulateTrawlProcess[volumes_, seed_, n_, foldGen : (FoldGenerate | ParallelFoldGenerate) : FoldGenerate] :=
    foldGen[
        Plus,
        MovingTotalC[RandomVariate[seed[#[[1]]], n + #[[2]] - 1], #[[2]]] &,
        ConstantArray[0.0, n],
        Transpose @ {volumes, Range @ Length @ volumes}
    ];
    
SimulateTrawlProcess /: Parallelize[SimulateTrawlProcess[kernel_, seed_, n_]] :=
    SimulateTrawlProcess[kernel, seed, n, ParallelFoldGenerate];


End[] (* End Private Context *)

EndPackage[]