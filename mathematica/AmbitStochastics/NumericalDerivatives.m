(* Mathematica Package *)

BeginPackage["AmbitStochastics`NumericalDerivatives`"]
(* Exported symbols added here with SymbolName::usage *)  

NumericalPartialDerivative::usage = "NumericalPartialDerivative[tbl, idx, delta, \"Cyclic\"|\"Linear\"]"
EnergyDissipationCyclic::usage = "EnergyDissipationCyclic[vx, vy, vz, delta]"

Begin["`Private`"] (* Begin Private Context *) 

d8kernel = 
    {-0.0035714285714285713`, 0.0380952380952381`, -0.2`, 0.8`, 0.`, -0.8`, 0.2`, -0.0380952380952381`, 0.0035714285714285713`};

Clear[NumericalPartialDerivative]

NumericalPartialDerivative[tbl_ ? VectorQ, 1, delta_, mode : ("Cyclic" | "Linear") : "Linear"] :=
    Switch[mode,
        "Cyclic", delta^-1 * ListConvolve[d8kernel, tbl, {5, 5}],
        "Linear", delta^-1 * ListConvolve[d8kernel, tbl, {-1, 1}]
    ];

NumericalPartialDerivative[tbl_ ? (ArrayQ[#] && ArrayDepth[#] > 1 &), idx_, delta_, mode : ("Cyclic" | "Linear") : "Linear"] :=
    Module[ {n, perm, transp, dim, res},
	    n      = ArrayDepth[tbl];
	    perm   = If[ idx == n, Cycles[{}], Cycles[{{n, idx}}] ];
	    transp = Permute[Range[n], perm];
	    dim    = Most @ Permute[Dimensions[tbl], perm];
	    res    = Transpose[Developer`ToPackedArray @ tbl, transp];
	    Switch[mode,
	    	"Cyclic", Array[ res[[##]] = delta^-1 * ListConvolve[d8kernel, res[[##]], {5, 5}]; &, dim ];,
	    	"Linear", res = Developer`ToPackedArray @ Array[delta^-1 * ListConvolve[d8kernel, res[[##]], {-1, 1}]; &, dim ];
	    ];
	    Transpose[res, transp]
    ];

Clear[EnergyDissipationCyclic]

EnergyDissipationCyclic[vxyz_ ? (ArrayQ[#, 4] &), args___] :=
    EnergyDissipationCyclic[vxyz[[1]], vxyz[[2]], vxyz[[3]], args];
    
EnergyDissipationCyclic[vx_ ? (ArrayQ[#, 3] &), vy_ ? (ArrayQ[#, 3] &), vz_ ? (ArrayQ[#, 3] &), delta_, debug : (True | False) : False] := 
	Module[{endis, ctr = 0},
        If[Not[Dimensions[vx] == Dimensions[vy] == Dimensions[vz]], Return[$Failed]];
        endis = ConstantArray[0.0, Dimensions[vx]];
        If[debug, PrintTemporary[++ctr]];
        endis += 4 * NumericalPartialDerivative[vx, 1, delta, "Cyclic"] ^ 2;
        If[debug, PrintTemporary[++ctr]];
        endis += 4 * NumericalPartialDerivative[vy, 2, delta, "Cyclic"] ^ 2;
        If[debug, PrintTemporary[++ctr]];
        endis += 4 * NumericalPartialDerivative[vz, 3, delta, "Cyclic"] ^ 2;
        If[debug, PrintTemporary[++ctr]];
        endis += 2 * (NumericalPartialDerivative[vx, 2, delta, "Cyclic"] + NumericalPartialDerivative[vy, 1, delta, "Cyclic"]) ^ 2;
        If[debug, PrintTemporary[++ctr]];
        endis += 2 * (NumericalPartialDerivative[vx, 3, delta, "Cyclic"] + NumericalPartialDerivative[vz, 1, delta, "Cyclic"]) ^ 2;
        If[debug, PrintTemporary[++ctr]];
        endis += 2 * (NumericalPartialDerivative[vy, 3, delta, "Cyclic"] + NumericalPartialDerivative[vz, 2, delta, "Cyclic"]) ^ 2;
        endis
	];

End[] (* End Private Context *)

EndPackage[]