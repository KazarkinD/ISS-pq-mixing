(* ::Package:: *)

(* ============================================
   iss_functions.m
   Custom functions for calculations within the framework of 
   Inverse seesaw ISS(p,q)
   
   version 1.0: 9 August 2026
   ============================================
*)

BeginPackage["ISS`"];

(**********************************************************************)
(*                         Main functions                             *)
(**********************************************************************)

logSpace[min_, max_, n_] := 10^Subdivide[Log10[min], Log10[max], n-1]

logUniformRandomComplex[{zmin_, zmax_}, dims_] := 
	Module[
	{xmin, xmax, ymin, ymax, rmatrix},
	
	xmin = Re[zmin];
	ymin = Im[zmin];
	xmax = Re[zmax];
	ymax = Im[zmax];
	
	rmatrix = Exp[RandomReal[{Log[xmin],Log[xmax]}, dims]] + I*Exp[RandomReal[{Log[ymin], Log[ymax]}, dims]]
]

(**********************************************************************)

massMatrixISSpq::usage = "Inverse seesaw mass matrix M = {{0, mD, 0}, {mD^T, 0, MR}, {0, MR^T, mu}"
massMatrixISSpq[p_:issP, q_:issQ, scales_:issScales, YD_, YR_, Ymu_]:=
If[	SymmetricMatrixQ[Ymu] && Dimensions[YD] === {3,p} && Dimensions[YR] === {p,q} && Dimensions[Ymu] === {q,q} && VectorQ[scales] && Length[scales] == 3,
	ArrayFlatten[
			{
				{ConstantArray[0, {3, 3}],       YD*scales[[1]],        ConstantArray[0, {3, q}] },
				{Transpose[YD*scales[[1]]],   ConstantArray[0, {p, p}],     YR*scales[[2]]       },
				{ConstantArray[0, {q, 3}],   Transpose[YR*scales[[2]]],     Ymu*scales[[3]]      }
			}
		],
	Print["massMatrixISSpq:ABORT: Yukawa input matrices or scales have incorrect form"];
	Abort[]
]

(**********************************************************************)

takagiDecomposition::usage = "Evaluate {{s_1, s_2, ...}, U} singular values s_i and unitary matrix U for input symmetric matrix M: M = UDU^T, D = diag(s_1,...)"
takagiDecomposition[matrix_?MatrixQ, prec_: 50]:= 
If[ SymmetricMatrixQ[matrix],
	Module[
		{U, S, P, V, M, Umix, ord, phases},
		
		(*SVD*)
		M = SetPrecision[matrix, prec];
		{U, S, V} = SingularValueDecomposition[M, Tolerance->0];
		
		(* Reordering with increasing of s_i *)
		ord = Ordering[Diagonal[S]]; 
		{U, S, V} = {U[[All, ord]], S[[ord, ord]], V[[All, ord]]};
			
		(* Rephasing gives real singular values *)
		phases = Arg[Diagonal[S]];
		P = DiagonalMatrix[Exp[-I*phases/2]];
		Umix = U . P;
		Return[{Abs[Diagonal[S]], Umix}]
	],
	Print["takagiDecomposition:ABORT: Matrix ", MatrixForm[matrix] ," is NOT symmetric \n"];
	Abort[]
]

(**********************************************************************)

getUUmixing::usage = "Numerical calculated sum over flavours of absolute squared values of mixing matrix |U|^2 with fixed HNL-generation for obvious mass matrix M"
getUUmixing[matrix_, prec_:50, Ihnl_]:= 
If[SymmetricMatrixQ[matrix],
	Module[
	{U, UUmix},
	
	U = takagiDecomposition[matrix, prec][[2]];

	UUmix = Re[Sum[Conjugate[U[[a,Ihnl+3]]]*U[[a,Ihnl+3]] ,{a, 1, 3}]];
	Return[UUmix]
],
	Print["getUUmixing:ABORT: Mass matrix is not symmetric"];
	Abort[]
]

(**********************************************************************)

getUUalphaMixing::usage = "Numerical calculated Sum over HNL-generations of absolute squared values of mixing matrix |U|^2 with fixed flavour alpha for obvious mass matrix M"
getUUalphaMixing[matrix_, prec_:50, alpha_]:= 
If[SymmetricMatrixQ[matrix],
	Module[
	{U, UUalpha, d},

	d = Dimensions[matrix][[1]];
	U = takagiDecomposition[matrix, prec][[2]];
	
	UUalpha = Re[Sum[Conjugate[U[[alpha,Ihnl+3]]]*U[[alpha,Ihnl+3]] ,{Ihnl, 1, d-3}]];
	Return[UUalpha]
	],
	Print["getUUalphaMixing:ABORT: Mass matrix is not symmetric"];
	Abort[]
]


(**********************************************************************)
(*       Semi-analytical methods for evaluation of mixing matrix      *)
(**********************************************************************)

getMassMatrixBlocks::usage = "Extract blocks {mD, MR, mu} from ISS(p,q) full mass matrix"
getMassMatrixBlocks[p_:issP, q_:issQ, massMatrix_]:= Module[
	{md, mr, mu},

	md = massMatrix[[1 ;; 3, 4 ;; 3+p]]; 
	mr = massMatrix[[4 ;; 3+p, 4+p ;; 3+p+q ]]; 
	mu = massMatrix[[4+p ;; 3+p+q, 4+p ;; 3+p+q]];

	Return[{md, mr, mu}]
]

getThetaMatrix::usage = "getThetaMatrix[P, Q, M] evaluates (analytical) anti-Hermitian matrix \[Theta] for ISS(P,Q), where M is mass matrix"
getThetaMatrix[p_:issP, q_:issQ, massMatrix_?MatrixQ]:= 
If[p + q + 3 == Length[massMatrix],
	Module[
		{theta1, theta2, md, mr, mu},

		{md, mr, mu} = getMassMatrixBlocks[p,q,massMatrix];
		
		theta1 = md . Inverse[mr . Inverse[mu] . Transpose[mr]]; 
		theta2 = theta1 . mr . Inverse[mu];
		
		Return[{theta1, theta2}]
	],
	Print["getThetaMatrix:ABORT: p and q are not consistent with the size of input mass matrix \n"];
	Abort[]
]
(**********************************************************************)

(* 
getThetaCasasIbarra::usage =
getThetaCasasIbarra[]:= 
 *)

(**********************************************************************)
	
getMixMatrix::usage = "getMixMatrix[{\[Theta]1,\[Theta]2}, \!\(\*SubscriptBox[\(U\), \(a\)]\), \!\(\*SubscriptBox[\(U\), \(S\)]\)] evaluates mixing matrix U = W.V (semi-analytical) using W = {{I, \[Theta]},{\!\(\*SuperscriptBox[\(\[Theta]\), \(+\)]\),I}} (first order accuracy), \[Theta] = {\[Theta]1, \[Theta]2}, U = {{I, \[Theta]},{\!\(\*SuperscriptBox[\(\[Theta]\), \(+\)]\),I}}.{{\!\(\*SubscriptBox[\(U\), \(a\)]\),0},{0,\!\(\*SubscriptBox[\(U\), \(s\)]\)}}"	
getMixMatrix[tildeTheta_, Uactive_?MatrixQ, Usterile_?MatrixQ]:= 
Module[
	{w, ww, vv, U, tTheta},

	tTheta = Join[tildeTheta[[1]], tildeTheta[[2]], 2];
	w = ArrayFlatten[{{0, -tTheta },{ConjugateTranspose[tTheta], 0}}];
	ww = IdentityMatrix[Length[w]] + w; (* only first term of W = exp(w) are taken*)
	(* ww = MatrixExp[w];*) (* Num. evaluation of full matrix W = exp(w)*)
	vv = ArrayFlatten[{{Uactive, 0},{0, Usterile}}];
	U = ww . vv;
	Return[U]
]

(**********************************************************************)

Utheta[p_?IntegerQ, q_?IntegerQ, massMatrix_?MatrixQ ] := 
If[ p + q + 3 == Length[massMatrix],
	Module[
		{md, mr, mu, mnu, u},

		x = massMatrix[[4 ;; 3+p+q, 4 ;; 3+p+q]]; (* sterile block *)
		{md, mr, mu} = getMassMatrixBlocks[p, q, massMatrix];
		mnu = md . Inverse[mr . Inverse[mu]. Transpose[mr]] . Transpose[md]; (* active neutrino mass operator *)
		mnu = 0.5*(mnu+Transpose[mnu]);
		
		u = getMixMatrix[
				getThetaMatrix[p, q, massMatrix], (* theta1 and theta2*)
				takagiDecomposition[mnu][[2]], (* U_nu from num. result*)
				takagiDecomposition[x][[2]] (* U_s from num. result*)
				];
		Return[u]
	],
	Print["Utheta:ABORT: p and q are not consistent with the size of input mass matrix \n"];
	Abort[]	
]

(**********************************************************************)
(***********       Casas-Ibarra parametrization             ***********)
(**********************************************************************)

yukawaDiracCasasIbarra::usage = "Yukawa matrix for active neutrino in ISS(p,q) using PMNS-data and Delta m^2-data (Casas-Ibarra parametrization)"
yukawaDiracCasasIbarra[hierarchy_, mlight_:0, p_:issP, q_:issQ, scales_:issScales, Omega_, YR_, Ymu_]:= 
Module[
	{UB, B, sngvlsB, Y},

	B = (scales[[2]]*scales[[2]])/scales[[3]]*YR.Inverse[Ymu].Transpose[YR]; (* B is (x|mu) without(!) inverse *)
	UB = takagiDecomposition[B][[2]]; (* B = UB.D.UB^T*)
	sngvlsB = takagiDecomposition[B][[1]];

	Switch[
		hierarchy, 	
		"NH",
			If[ p>2, (*for NH*)
				Y = 1/scales[[1]]*UpmnsNH.DiagonalMatrix[{Sqrt[mlight], Sqrt[mNH[[2]]], Sqrt[mNH[[3]]]}].Omega.DiagonalMatrix[Sqrt[sngvlsB]].Transpose[UB],

				If[ p == 2,
					(*p=2*) Y=1/scales[[1]]*UpmnsNH.DiagonalMatrix[{Sqrt[mlight], Sqrt[mNH[[2]]], Sqrt[mNH[[3]]]}].ArrayFlatten[{{ConstantArray[0, {1, 2}]}, {Omega}}].DiagonalMatrix[Sqrt[sngvlsB]].Transpose[UB],
					(*p<2*) Print["yukawaDiracCasasIbarra:ABORT: Less than 1 massive active neutrino"];
					Abort[]
				]
			],
		"IH",
			If[ p>2, (*for IH*)
				Y = 1/scales[[1]]*UpmnsIH.DiagonalMatrix[{Sqrt[mlight], Sqrt[mIH[[2]]], Sqrt[mlight]}].Omega.DiagonalMatrix[Sqrt[sngvlsB]].Transpose[UB],

				If[ p == 2,
					(*p=2*) Y=1/scales[[1]]*UpmnsIH.DiagonalMatrix[{Sqrt[mIH[[1]]], Sqrt[mIH[[2]]], Sqrt[mlight]}].ArrayFlatten[{{Omega}, {ConstantArray[0, {1, 2}]}}].DiagonalMatrix[Sqrt[sngvlsB]].Transpose[UB],
					(*p<2*) Print["yukawaDiracCasasIbarra:ABORT: Less than 1 massive active neutrino"];
					Abort[]
				]
			],
		_,
			Print["yukawaDiracCasasIbarra:ABORT: Wrong hierarchy input form, use NH or IH"];
			Abort[]
		];
	Return[Y]
]

(**********************************************************************)

massMatrixISSpqCasasIbarra::usage = "Inverse seesaw mass matrix with Casas-Ibarra parametrization"
massMatrixISSpqCasasIbarra[hierarchy_, mlight_:0, p_:issP, q_:issQ, scales_:issScales, Omega_, YR_, Ymu_]:= 
massMatrixISSpq[
	p,
	q, 
	scales,
	yukawaDiracCasasIbarra[
		hierarchy,
		mlight,
		p,
		q,
		scales,
		Omega,
		YR,
		Ymu
	], 
	YR, 
	Ymu
]

(**********************************************************************)
(***********       General analysis for random Yukawa       ***********)
(**********************************************************************)

randomYDirac::usage = "Yukawa matrix for Dirac-type mass term, non-symmetric & non-squared in general"
randomYDirac[Nrows_, Ncolumns_, YukawaRegime_:"uniform", zmin_:zMIN, zmax_:zMAX]:= 
If[YukawaRegime == "log-uniform", 
	logUniformRandomComplex[{zmin,zmax}, {Nrows, Ncolumns}],
	If[YukawaRegime == "uniform", 
		RandomComplex[{zmin,zmax}, {Nrows, Ncolumns}], 
		Print["Incorrect value of 'YukawaRegime'. Use 'log-uniform' or 'uniform'."]; Return[$Failed]
	]
];

(**********************************************************************)

randomYMajorana::usage = "Yukawa matrix for Majorana-type mass term, ONLY symmetric & squared matrix"
randomYMajorana[Nrows_, YukawaRegime_:"uniform", zmin_:zMIN, zmax_:zMAX]:= 
Module[
	{lower, mu},

	lower = Switch[
		YukawaRegime,
		"log-uniform", logUniformRandomComplex[{zmin, zmax}, {Nrows, Nrows}],
		"uniform", RandomComplex[{zmin, zmax}, {Nrows, Nrows}],
		_, Print["Incorrect value of 'YukawaRegime'."]; Return[$Failed]
	];
	mu = lower + Transpose[lower] - DiagonalMatrix[Diagonal[lower]]
]

(**********************************************************************)

randomMassMatrixISSpq::usage = "Random mass matrix for ISS(p,q) with fixed scales <mD>, <MR> and <mu> and Casas-Ibarra-parametrization"
randomMassMatrixISSpq[hierarchy_, Omega_, mlight_:0, p_:issP, q_:issQ, scales_:issScales, YukawaRegime_:"uniform", zmin_:zMIN, zmax_:zMAX]:=
Module[
	{rYD, rYR, rYmu, rMassMatrix},

	rYR = randomYDirac[p, q, YukawaRegime, zmin, zmax];
	rYmu = randomYMajorana[q, YukawaRegime, zmin, zmax];
	
	rYD = yukawaDiracCasasIbarra[hierarchy, mlight, p, q, scales, Omega, rYR, rYmu];
	
	rMassMatrix = massMatrixISSpq[p, q, scales, rYD, rYR, rYmu];
	Return[rMassMatrix]
]

(**********************************************************************)

massGenData::usage = "Generates array of mass spectrum with `ngen` generations for random Yukawa constants Y_R and diagonal Y_mu (uniform, 0.1 < |Ymu| < 1) and Casas-Ibarra-parametrization for Y_D"
massGenData[hierarchy_, Omega_, mlight_:0, p_:issP, q_:issQ, scales_:issScales, YukawaRegime_:"log-uniform", zmin_:zMIN, zmax_:zMAX, ngen_:50]:= 
Module[
	{randomYtab},

	randomYtab = Table[ {
		randomYDirac[p, q, YukawaRegime, zmin, zmax], 
		DiagonalMatrix[ RandomComplex[{0.1 + 0.1*I, 1 + 1*I}, 3] ]
		}, {ngen} ];

	Table[
		takagiDecomposition[
			massMatrixISSpq[
							p,
							q, 
							scales,
							yukawaDiracCasasIbarra[
								hierarchy,
								mlight,
								p,
								q,
								scales,
								Omega,
								randomYtab[[i,1]], (* random YR,  i-th generation *)
								randomYtab[[i,2]]  (* random Ymu, i-th generation *)
							], 
							randomYtab[[i,1]],
							randomYtab[[i,2]]
			]
		][[1, 3 - Min[p,q,3] + 1 ;; p+q+3 ]], (* [3 - min(p,q,3)]- number of zero act nu masses *)
		{i, 1, ngen} 
	]
]
(**********************************************************************)

plotRandomMasses::usage = "Draws the plot with mass spectrum [GeV] of ISS(p,q) model using Casas-Ibarra parametrization"
plotRandomMasses[hierarchy_, Omega_, mlight_:0, p_:issP, q_:issQ, scales_:issScales, YukawaRegime_:"log-uniform", zmin_:zMIN, zmax_:zMAX, ngen_:100, plotrange_:Automatic ]:= 
Module[
	{massesFixedScales, tab},

	massesFixedScales = massGenData[hierarchy, Omega, mlight, p, q, scales, YukawaRegime, zmin, zmax, ngen];
	
	tab = Table[ {i + RandomReal[{-0.2, 0.2}], #} & /@ massesFixedScales[[i]], {i, Length[massesFixedScales]} ];	
	
	ListLogPlot[
		Transpose[tab],

		PlotStyle -> AbsolutePointSize[6],
		Frame -> True,
		FrameLabel -> {label["# of random Yukawa generation", largeFont],label["Mass spectrum [GeV]", largeFont]},
		FrameTicks -> {
			{Table[{10^k, If[ Mod[k,3] ==0, label[Superscript[10,ToString[k]]], None ]} , {k, -20, 20}], 
				{
					{issScales[[1]], colorlabel["EW-scale", Darker[Green], largeFont-2] },
					{issScales[[2]], colorlabel["\[CapitalLambda](\!\(\*SubscriptBox[\(M\), \(R\)]\))", Red, largeFont-2] },
					{issScales[[3]], colorlabel["\[CapitalLambda](\[Mu])", Blue, largeFont-2] },
					{Sqrt[massSquareDif[[1]]], colorlabel["\!\(\*SubscriptBox[\(m\), \(\[Nu]\)]\) solar", Darker[Red], largeFont-2] },
					{Sqrt[massSquareDif[[2]]], colorlabel["\!\(\*SubscriptBox[\(m\), \(\[Nu]\)]\) atm", Darker[Red], largeFont-2] }
				} 
			},
			{Table[{0+i*ngen/10, label[ToString[0+i*ngen/10]]}, {i, 0, 10}], None}},
		GridLines -> {
			All, 
			Join[
				Table[10^i, {i,-20,20}],
				{
						{issScales[[1]], Directive[Green, Thickness[0.008], Opacity[.5]]},
						{issScales[[2]], Directive[Red, Thickness[0.008], Opacity[.5]]},
						{issScales[[3]], Directive[Blue, Thickness[0.008], Opacity[.5]]},
						{Sqrt[massSquareDif[[1]]], Directive[DotDashed, Darker[Red], Thickness[0.004], Opacity[1.0]]},
						{Sqrt[massSquareDif[[2]]], Directive[Dashed, Darker[Red], Thickness[0.004], Opacity[1.0]]}
				}
			]
		},
		AspectRatio->1/GoldenRatio,
		ImageSize -> Large, (* largeFont for ImageSize = Large *) 
		PlotRange -> {1*^-14, 1*^6},
		PlotLabel -> label["[random \!\(\*SubscriptBox[\(Y\), \(R\)]\)], NH,   \[CapitalLambda](\!\(\*SubscriptBox[\(M\), \(R\)]\))=" <> ToString[ScientificForm[issScales[[2]], 2], TraditionalForm] <> 
				", \[CapitalLambda](\[Mu])=" <> ToString[ScientificForm[issScales[[3]], 2], TraditionalForm] <>
				" [GeV]", largeFont]
	]
]


(**********************************************************************)
(*******************       Custom plot style       ********************)
(**********************************************************************)

(* label and  legend style *)

myimagesize = {500,400};

largeFont = 16; 

label[text_, fontsize_:14]:= Text[Style[text, Black, FontFamily->"Helvetica", FontSize->fontsize]];

colorlabel[text_, color_,fontsize_:14]:= Text[Style[text, FontFamily->"Helvetica", FontSize->fontsize, FontColor->color]];

legend[legendNames_, coords_:{0.8,0.8}]:= 
	Placed[
		PointLegend[
			legendNames, 
			LegendFunction->(Framed[#1, Background -> Directive[White, Opacity[0.2]], RoundingRadius -> 3] & )
			],
	coords]

(**********************************************************************)
(**********************************************************************)
(**********************************************************************)

Print["iss_functions.m loaded successfully"]

EndPackage[]