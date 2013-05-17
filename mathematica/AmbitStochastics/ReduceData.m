(* Mathematica Package *)

BeginPackage["AmbitStochastics`ReduceData`"]
(* Exported symbols added here with SymbolName::usage *)  

AdjoinUniformWeightsAndSort::usage = "AdjoinUniformWeightsAndSort[x] will sort and adjoin uniform weights to the data points such that the sum of the weights is 1."

ReduceDataBy::usage = "ReduceBy[{x, w}, n] will reduce an ordered data set by clustering n points at a time. The sum of the new weights the same as the sum of the old weights."

ReduceDataByUntil::usage = "ReduceByUntil[{x, w}, n, m] will iteratively reduce a data set by n until the length of the dataset is as small as possible while not being smaller than m."

ReduceDataTo::usage = "ReduceTo "

ReducedDataToPDF::usage = "ReducedToPDF[{x, w}] converts the weights in w to the corresponding pdf values. This requires that the data in x is evenly spaced."

Begin["`Private`"] (* Begin Private Context *) 

Clear[AdjoinUniformWeightsAndSort]
AdjoinUniformWeightsAndSort[x_] := 
    {Sort[x], ConstantArray[1.0 / Length[x], Length[x]]};

Clear[ReduceDataBy]
ReduceDataBy[{x_, w_}, n_] := (* assumes ordered x *)
    Module[ {px, pw, tx, tw},
        px = Partition[x, n, n, {1, 1}, {}];
        pw = Partition[w, n, n, {1, 1}, {}];
        tw = Total[pw, {2}];
        tx = MapThread[#1 . #2 / #3 &, {px, pw, tw}];
        {tx, tw}
    ];
    
Clear[ReduceDataByUntil]
ReduceDataByUntil[{x_, w_}, n_, m_] := 
    Most @ NestWhileList[ReduceDataBy[#, n] &, {x, w}, Length[#] >= m &];

Clear[ReduceDataTo]
ReduceDataTo = 
    Compile[ {{x, _Real, 1}, {n, _Integer, 0}},
        Module[ {min, max, dx, idx, count},
            min = Min[x];
            max = Max[x];
            dx  = (max - min) / n;
            idx = Floor[(x - min) / dx];
            count = Table[0.0, {n + 1}];
            Do[count[[i + 1]]++, {i, idx}];
            count[[-2]] += count[[-1]];
            {
                (Range[0, n - 1] + 0.5) * dx + min,
                Most[count]
            }
        ],
        CompilationTarget -> "C"
    ];

Clear[ReducedDataToPDF]
ReducedDataToPDF[{x_, w_}] := 
    {x, w / (Total[w] Mean[Differences[Sort[x]]])};



End[] (* End Private Context *)

EndPackage[]