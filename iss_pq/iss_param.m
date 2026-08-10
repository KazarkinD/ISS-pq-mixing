(* ::Package:: *)

(* 
============================================
	parameters.m 
	Global set of parameters 
	for Inverse seesaw ISS(P,Q) model 
============================================
*)
BeginPackage["ISS`"]


(**********************************************************************)
(*                              Constants                             *)
(**********************************************************************)
vEW = 246;     (* GeV *)


(**********************************************************************)
(*                   Default model parameters                         *)
(**********************************************************************)

(* Number of fields *)
issP = 2;
issQ = 3;

(* Energy scales *)
mDscale = vEW*0.7; (*GeV*)
mRscale = 1.0*^4; (*GeV*)
muScale = 10.0*1*^-6; (*GeV*)
issScales = {mDscale, mRscale, muScale};
massSquareDif = {7.5*1*^-5, 2.4*1*^-3}*1*^-18; (*GeV^2*)

(* Numerical parameters of random Yukawa range *)
zMIN = 0.000001 + 0.000001 * I; 
zMAX = 1 + 1 * I;

(* Random seed (if needed) *)
(*seed = 1234; *)


(**********************************************************************)
(*                 PMNS-matrix    (NuFIT best fit)                    *)
(**********************************************************************)

\[Theta]23 = 48.5 Degree;
\[Theta]13 = 8.5 Degree;
\[Theta]12 = 33.7 Degree;
\[Delta]CPwIH = 270*Pi/180; (* for inverted ordering *)
\[Delta]CPwNH = 177*Pi/180; (* for normal ordering   *)
(*\[Delta]CP=0;*)
\[Eta]1= 0;
\[Eta]2= 0;

(*PMNS-matrix: parametrization form*)
U23[\[Theta]23_] := {
	{ 1,          0,              0        },
	{ 0,  Cos[\[Theta]23], Sin[\[Theta]23] },
	{ 0, -Sin[\[Theta]23], Cos[\[Theta]23] }
	};

U13[\[Theta]13_,\[Delta]CP_] := {
	{       Cos[\[Theta]13],              0,  Sin[\[Theta]13]*Exp[-I*\[Delta]CP] },
	{             0,                      1,                    0                },
	{-Sin[\[Theta]13]*Exp[I*\[Delta]CP],  0,            Cos[\[Theta]13]          }
	};

U21[\[Theta]21_] := {
	{ Cos[\[Theta]21],  Sin[\[Theta]21],  0 },
	{-Sin[\[Theta]21],  Cos[\[Theta]21],  0 },
	{ 0,                0,                1 }
	};

MajoranaPhases[\[Eta]1_, \[Eta]2_]:={
	{ Exp[I*\[Eta]1],        0,           0 },
	{ 0,              Exp[I*\[Eta]2],     0 },
	{ 0,                     0,           1 }
	};

pmns[\[Theta]23_,\[Theta]13_,\[Theta]21_,\[Delta]CP_,\[Eta]1_,\[Eta]2_]:=
	U23[\[Theta]23] . U13[\[Theta]13,\[Delta]CP] . U21[\[Theta]21] . MajoranaPhases[\[Eta]1,\[Eta]2];


(*Upmns={{Ue1,Ue2,Ue3},{U\[Mu]1,U\[Mu]2,U\[Mu]3},{U\[Tau]1,U\[Tau]2,U\[Tau]3}}*)
(*TBMix={{2/3, 1/3, 0},{1/6, 1/3, 1/3},{1/6, 1/3, 1/2}}; *)
(*Upmns=TBMix;*)
UpmnsNH = pmns[\[Theta]23,\[Theta]13,\[Theta]12,\[Delta]CPwNH,\[Eta]1,\[Eta]2];
UpmnsIH = pmns[\[Theta]23,\[Theta]13,\[Theta]12,\[Delta]CPwIH,\[Eta]1,\[Eta]2];

mNH = {0, Sqrt[massSquareDif[[1]]], Sqrt[massSquareDif[[2]]+massSquareDif[[1]]]};
mIH = {Sqrt[massSquareDif[[2]]-massSquareDif[[1]]], Sqrt[massSquareDif[[2]]] ,0};


Print[
	"Default parameters are loaded: \n P=", issP, ", Q=", issQ,
	"\n Scales: \n",
	"EW-scale = ", mDscale, 
	" GeV, \n \[CapitalLambda](\!\(\*SubscriptBox[\(M\), \(R\)]\)) = ", mRscale, 
	" GeV, \n \[CapitalLambda](\[Mu]) = ", muScale, 
	" GeV, \n Default random range: \n",
	"zMIN = ", ScientificForm[Re[zMIN],2], "+ i*", ScientificForm[Im[zMIN],2], 
	"\n zMAX = ", ScientificForm[Re[zMAX],2], "+ i*", ScientificForm[Im[zMAX],2],
	"\n Active neutrino masses: \n",
	"NH:", mNH, "\n IH:", mIH
]

EndPackage[]