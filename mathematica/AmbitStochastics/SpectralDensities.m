(* Mathematica Package *)

BeginPackage["AmbitStochastics`SpectralDensities`"]
(* Exported symbols added here with SymbolName::usage *)  

SpectralDensity::usage = "FIXME";

SimplePeriodogram::usage = 
"SimplePeriodogram[x, delta] calculates the SimplePeriodogram where delta is the resolution.";

WOSASpectralDensityEstimate::usage = 
"WOSASpectralDensityEstimate[x, delta, taper, blocksize] estimates the spectral density using \
Welsh's Overlapping Segment Averages using the given taper and blocksize.";

RectangularTaper::usage = 
"RectangularTaper gives the rectangular taper, i.e., no tapering will be performed \
which is a bad idea.";

HanningTaper::usage = 
"HanningTaper gives the Hanning taper which is a good taper.";

ApproximateDPSS0Taper::usage = 
"ApproximateDPSS0Taper gives a taper using an approximation of the zeroth discrete \
prolate spheriodal sequence.";

AdjoinFrequencies::usage = 
"AdjoinFrequencies converts the output of the spectral density estimators into \
something suitable for plotting.";

DirectSpectralDensityEstimate::usage =
"DirectSpectralDensityEstimate[x, delta, taper]";


Begin["`Private`"] (* Begin Private Context *) 

SpectralDensity[data_, samplingResolution_, segmentLength_, segmentOverlap_ : 0.5, windowFunction_ : HannWindow] :=
    Module[{p},
        p = PeriodogramArray[data - Mean[data], segmentLength, Floor[segmentLength * (1 - segmentOverlap)], windowFunction];
        p = samplingResolution * p[[2 ;; Floor[Length[p]/2 + 1]]];
        Transpose @ {Range[Length[p]]/(2 * samplingResolution * Length[p]), p} 
    ];



Clear[AdjoinFrequencies]
AdjoinFrequencies[{s_, {omegaMin_, omegaMax_, deltaomega_}}] := 
  Transpose[{Range[deltaomega, omegaMax, deltaomega], Rest[s]}];

Clear[InnerSpectralDensityEstimate]
InnerSpectralDensityEstimate[x_, deltat_] := 
    Module[ {nn, kk, deltaomega, omegaMin, omegaMax, sdf},
        nn = Length[x];
        kk = Floor[nn/2];
        deltaomega = 1/(nn deltat); 
        omegaMin = 0; 
        omegaMax = kk * deltaomega;
        sdf = deltat Abs[Take[Fourier[x, FourierParameters -> {1, -1}], kk + 1]]^2;
        {sdf, {omegaMin, omegaMax, deltaomega}}
    ];

Clear[SimplePeriodogram]
SimplePeriodogram[x_, deltat_] := 
  DirectSpectralDensityEstimate[x, deltat, RectangularTaper];

Clear[DirectSpectralDensityEstimate]
DirectSpectralDensityEstimate[x_, deltat_, taper_] := 
    Module[ {h},
        h = If[ Head[taper] === List, 
            taper, 
            taper[Length[x]]
        ];
        InnerSpectralDensityEstimate[h * (x - Mean[x]), deltat]
    ];

Clear[WOSASpectralDensityEstimate]
WOSASpectralDensityEstimate[x_, deltat_, taper_, blockSize_, overlap_: 0.5] := 
    Module[ {h, blocks, sdfs},
    h = If[ Head[taper] === List, 
        taper, 
        taper[blockSize]
    ];
    blocks = Map[# - Mean[#] &, Partition[x, blockSize, Floor[(1 - overlap) * blockSize]]];
    sdfs = Map[InnerSpectralDensityEstimate[h * #, deltat] &, blocks];
    {Mean[First /@ sdfs], Last @ First @ sdfs}
    ];


(* Tapers *)

Clear[RectangularTaper]
RectangularTaper[n_] := 
    ConstantArray[N[n ^ (-1 / 2)], n];

Clear[HanningTaper]
HanningTaper[n_] := 
    Sqrt[2 / (3 * (n + 1))] * (1 - Cos[(2 * Pi) / (n + 1) * N @ Range[1, n]]);

Clear[ApproximateDPSS0Taper]
ApproximateDPSS0Taper[n_, nw_: 4] := 
    Module[ {w, t, b},
        w = nw / n;
        t = N @ Range[0, n - 1];
        b = BesselI[0, Pi * w * (n - 1) * Sqrt[1 - ((2 * t + 1 - n) / n) ^ 2]] / N @ BesselI[0, Pi * w * (n - 1)];
        b / Sqrt[Total[b ^ 2]]
    ];



End[] (* End Private Context *)

EndPackage[]

