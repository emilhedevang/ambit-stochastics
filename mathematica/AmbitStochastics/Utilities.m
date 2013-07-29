(* Mathematica Package *)

BeginPackage["AmbitStochastics`Utilities`"]
(* Exported symbols added here with SymbolName::usage *)  

ListNames::usage = "ListNames[]"

LogLogSlope::usage = "FIXME";

StructureFunction::usage  = "StructureFunction[x, lags, order]";
CorrelatorFunction::usage = "CorrelatorFunction[x, lags, {p, q}] calculates E[x(l)^p*x(0)^q] / (E[x^p]*E[x^q]) where l ranges over lags.";

NumberVectorQ::usage = "NumberVectorQ is equivalent to VectorQ[#, NumberQ] &.";

FoldGenerate::usage         = "FoldGenerate[op, gen, zero, xs] is equivalent to Fold[op[#1, gen[#2]] &, zero, xs], where zero should satisfy op[zero, x] = op[x, zero] = x for all x.";
ParallelFoldGenerate::usage = "ParallelFoldGenerate is a parallel version of FoldGenerate.";

LogScale::usage        = "LogScale[a, b, n] gives a list of n logarithmically spaced values from a to b.";
IntegerLogScale::usage = "IntegerLogScale[a, b, n] gives approximately n logarithmically spaced values from a to b."
LogTake::usage         = "LogTake[x, n] gives approximately n logarithmically spaced elements of x.";

Begin["`Private`"] (* Begin Private Context *) 

Clear[ListNames]
ListNames[pattern_ : "Global`*"] :=
    TableForm @ SortBy[Select[{#, ByteCount[Symbol[#]]} & /@ Names["Global`*"], #[[2]] > 0 &], -#[[2]] &];

Clear[LogLogSlope]
LogLogSlope[xy_] :=
    Transpose @ {xy[[1 ;; -2, 1]], Differences[Log[xy[[All, 2]]]]/ Differences[Log[xy[[All, 1]]]]};

Clear[StructureFunction]
StructureFunction[x_, lags_, order_, transform : (Identity | Abs) : Identity, tbl : (Table | ParallelTable) : Table] :=
    tbl[Mean[transform[(Drop[x, l] - Drop[x, -l])] ^ order], {l, lags}];
    
StructureFunction /: Parallelize[StructureFunction[x_, lags_, order_, transform : (Identity | Abs) : Identity]] :=
    StructureFunction[x, lags, order, transform, ParallelTable];

Clear[CorrelatorFunction]
CorrelatorFunction[x_, lags_, order : {p_, q_} : {1, 1}, tbl : (Table | ParallelTable) : Table] :=
	Module[{xp, xq, meanxp, meanxq},
		xp = x ^ p; meanxp = Mean[xp];
		xq = x ^ q; meanxq = Mean[xq];
		tbl[Mean[Drop[xp, l] * Drop[xq, -l]], {l, lags}] / (meanxp * meanxq)
	];

CorrelatorFunction /: Parallelize[CorrelatorFunction[x_, lags_, order : {p_, q_} : {1, 1}]] :=
    CorrelatorFunction[x, lags, order, ParallelTable];
    
    
Clear[NumberVectorQ]
NumberVectorQ[x_] :=
    VectorQ[x, NumberQ];

Clear[FoldGenerate]
FoldGenerate[op_, gen_, zero_, xs_] :=
    Fold[op[#1, gen[#2]] &, zero, xs];

FoldGenerate /: Parallelize[FoldGenerate[op_, gen_, zero_, xs_]] :=
    ParallelFoldGenerate[op, gen, zero, xs];

Clear[ParallelFoldGenerate]
ParallelFoldGenerate[op_, gen_, zero_, xs_] :=
    ParallelCombine[FoldGenerate[op, gen, zero, #] &, xs, Fold[op, zero, {##}] &];

Clear[LogScale]
LogScale[a_, b_, n_Integer] /; 0 < a < b && n > 1 :=
    a * (b / a) ^ (Range[0, n - 1] / (n - 1));
    
Clear[IntegerLogScale]
IntegerLogScale[a_Integer, b_Integer, n_Integer] /; 0 < a < b && n > 1 :=
    Union[Clip[Floor[LogScale[a, b, n]], {1, b}]];

Clear[LogTake]
LogTake[x_, n_Integer] /; n > 1 :=
    Part[x, IntegerLogScale[1, Length[x], n]];

End[] (* End Private Context *)

EndPackage[]