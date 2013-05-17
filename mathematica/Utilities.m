(* ::Package:: *)

BeginPackage["AmbitStochastics`Utilities`"];

LogScale::usage = "LogScale[a, b, n] gives n points on a logarithmic scale from a to b.";

IntegerLogScale::usage = "IntegerLogScale[a, b, n] gives at most n positive integers on a log scale from a to b.";

LogTake::usage = "LogTake[x, n] takes at most n elements indexed by an integer log scale from 1 to Length[x].";

Begin["`Private`"]; (* Begin Private Context *)

(* Logarithmic scales *)
Clear[LogScale]
LogScale[a_, b_, n_] := 
    Module[ {alpha, beta},
        beta  = (b/a) ^ (1 / (n - 1));
        alpha = a / beta;
        alpha * beta ^ Range[n]
    ];
    
Clear[IntegerLogScale]
IntegerLogScale[a_, b_, n_] := 
    Union[Clip[Floor[LogScale[a, b, n]], {1, b}]];

Clear[LogTake]
LogTake[x_, n_] := 
  Part[x, IntegerLogScale[1, Length[x], n]];

End[]; (* End Private Context *)

EndPackage[];
