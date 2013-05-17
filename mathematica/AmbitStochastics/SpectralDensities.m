(* Mathematica Package *)

BeginPackage["AmbitStochastics`SpectralDensities`"]
(* Exported symbols added here with SymbolName::usage *)  

SpectralDensity::usage = "FIXME";

Begin["`Private`"] (* Begin Private Context *) 

SpectralDensity[data_, samplingResolution_, segmentLength_, segmentOverlap_ : 0.5, windowFunction_ : HannWindow] :=
    Module[{p},
        p = PeriodogramArray[data - Mean[data], segmentLength, Floor[segmentLength * (1 - segmentOverlap)], windowFunction];
        p = samplingResolution * p[[2 ;; Floor[Length[p]/2 + 1]]];
        Transpose @ {Range[Length[p]]/(2 * samplingResolution * Length[p]), p} 
    ];

End[] (* End Private Context *)

EndPackage[]