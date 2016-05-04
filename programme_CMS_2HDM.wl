(* ::Package:: *)

BeginPackage["THDMfull`"]

Begin["Private`"]

Print["Begin..."]



$HistoryLength=1

SetDirectory["/home/slecorre/Documents/CMS/2HDMC-1.7.0/"]
Print["Directory : /home/slecorre/Documents/CMS/2HDMC-1.7.0"] 

(* Loading the files *)
(* 2HDM TYPE I *)
(* load the data contained in the files to form Mathematica arrays. Each line of theses arrays contains the whole parameter set (kappa, masses, limits...) of ONE specific point of the 2HDM *)



(* select the folder corresponding to the data you are interested in *)
(*dirOutfile="all_range_sin_BA/";*) (* data from a large scan (m12=0, mh in 80 to 120GeV, sinBA between -1 and +1... *)
dirOutfile="all_range_sin_BA_m12_30GeV/"; (* data from a large scan (m12=30GeV, mh in 80 to 120GeV, sinBA between -1 and +1... *)
(*dirOutfile="sinBA_-04_-01/";*) (* data from a scan with sinBA between -0.4 and -0.1 *)
(*dirOutfile="m12_30_40GeV/";*) (* data from a scan with m12 between 30 and 40GeV ; all parameters with the whole range *)
(*dirOutfile="plan_tanB_m12_sinBA_-01/";*) (* data from scan with sinBA=-0.1, with only tanB and m12 varying *)
(*dirOutfile="plan_tanB_m12_sinBA_-03/";*) (* data from scan with sinBA=-0.3, with only tanB and m12 varying *)
(*dirOutfile="plan_tanB_sinBA/";*)

Print["The current files directory is : ", dirOutfile]

Array101=Get[StringJoin["data/",dirOutfile,"outfile101"]];Print["file 1"];
Array102=Get[StringJoin["data/",dirOutfile,"outfile102"]];Print["file 2"];
Array103=Get[StringJoin["data/",dirOutfile,"outfile103"]];Print["file 3"];
Array104=Get[StringJoin["data/",dirOutfile,"outfile104"]];Print["file 4"];
Array105=Get[StringJoin["data/",dirOutfile,"outfile105"]];Print["file 5"];
Array106=Get[StringJoin["data/",dirOutfile,"outfile106"]];Print["file 6"];
Array107=Get[StringJoin["data/",dirOutfile,"outfile107"]];Print["file 7"];
Array108=Get[StringJoin["data/",dirOutfile,"outfile108"]];Print["file 8"];
Array109=Get[StringJoin["data/",dirOutfile,"outfile109"]];Print["file 9"];
Array110=Get[StringJoin["data/",dirOutfile,"outfile110"]];Print["file 10"];
Array111=Get[StringJoin["data/",dirOutfile,"outfile111"]];Print["file 11"];
Array112=Get[StringJoin["data/",dirOutfile,"outfile112"]];Print["file 12"];
Array113=Get[StringJoin["data/",dirOutfile,"outfile113"]];Print["file 13"];
Array114=Get[StringJoin["data/",dirOutfile,"outfile114"]];Print["file 14"];
Array115=Get[StringJoin["data/",dirOutfile,"outfile115"]];Print["file 15"];
Array116=Get[StringJoin["data/",dirOutfile,"outfile116"]];Print["file 16"];
Array117=Get[StringJoin["data/",dirOutfile,"outfile117"]];Print["file 17"];
Array118=Get[StringJoin["data/",dirOutfile,"outfile118"]];Print["file 18"];
Array119=Get[StringJoin["data/",dirOutfile,"outfile119"]];Print["file 19"];
Array120=Get[StringJoin["data/",dirOutfile,"outfile120"]];Print["file 20"];

Print["end array THDM type 1"]



(* Implementation of the flavor and theoretical limites *)
(* The uncertainties written here are combination of experimental uncertainties and theoretical uncertainties (which plays the role of the uncertainty on the computed data in 2HDMC, SuperIso...) *)

(* Experimental limits on B\[Rule] Xs \[Gamma] *)
limBXsgam=3.43*10^-4;
limBXsgamUp=0.29*10^-4;(* upper uncertainty on the value at 1\[Sigma] *)
limBXsgamDown=0.29*10^-4;(* lower uncertainty on the value at 1\[Sigma]*)

(* Experimental limits on Bs\[Rule] \[Mu]\[Mu] at 1\[Sigma] *)
limBsmumu=2.9*10^-9;
limBsmumuUp=0.8*10^-9;
limBsmumuDown=0.8*10^-9;

(* Experimental limits on \[CapitalDelta]Md at 1\[Sigma] *)
limDeltaMd=0.510;
limDeltaMdUp=0.091;
limDeltaMdDown=0.091;

(* Experimental limit on \[CapitalDelta]0 (isospin asymmetry for B\[Rule]K*\[Gamma] *)
limDelta0=0.052;
limDelta0Up=0.030;
limDelta0Down=0.030;

(* Experimental limits on S parameter at 1\[Sigma] *)
limS=0.05;
limSUp=0.11;
limSDown=0.11;

(* Experimental limits on T parameter at 1\[Sigma] *)
limT=0.09;
limTUp=0.13;
limTDown=0.13;

(* Experimental limits on U parameter at 1\[Sigma] *)
limU=0.01;
limUUp=0.11;
limUDown=0.11;

(* Minimum of the khi2 distribution obtained with the programm khi2_min.nb with 6dof *)
minLHC=\!\(TraditionalForm\`1.52653\);
boundLHC2\[Sigma]=12.8489; (* bound at 2\[Sigma] for a khi2 with 6dof*)
boundLHC=boundLHC2\[Sigma]+minLHC; (* upper bound taking into account the min of the khi2 distribution (minLHC) *)




(* Initialization of variables and extraction of data *)



(* /!\ Depending on the type of 2HDM you want to study, uncomment one of the two blocks below. *)


(*block for a 2HDM TYPE I  *)
txt="graph/Type1/" (* usefull for the graph's exportation *)
type="Type1" (* usefull for the graph's exportation *)
Array2HDM=Join[Array101,Array102,Array103,Array104,Array105,Array106,Array107,Array108,Array109,Array110,Array111,Array112,Array113,Array114,Array115,Array116,Array117,Array118,Array119,Array120]; (*compact the multiple arrays into one *)

Print[Length[Transpose[Array2HDM]]];
Print[Length[Array2HDM]];

(*
(*block for a 2HDM TYPE II  *)
txt="graph/Type2/"
type="Type2"
Array2HDM=Join[Array201,Array202,Array203,Array204,Array205,Array206,Array207,Array208,Array209,Array210,Array211,Array212,Array213,Array214,Array215,Array216,Array217,Array218,Array219,Array220];
*)



(* extraction of some data *)
(* from the big array Array2HDM containing many informations, we extract some interesting parameters and put them into Mathematica vectors. *)


(* kappas for h1 (light higgs) /!\ it is kappa^2 !!! *)
kg2=Array2HDM[[All,87]];
k\[Gamma]2=Array2HDM[[All,85]];
(*kv2=Array2HDM[[All,89]];*)
kv2=Array2HDM[[All,124]]^2; (* here we take kv=sin(B-A) *)
kb2=Array2HDM[[All,93]];
kt2=Array2HDM[[All,105]];

(* BR and masses *)
BRhGaGa=Array2HDM[[All,65]];
BRhGaGaSM=Array2HDM[[All,45]];

(*WidthSMh1=Array2HDM[[All,15]]; (* total width for a SM-like higgs boson with mass mh1 *)
WidthTHDMh1=Array2HDM[[All,17]]; (* total width for a THDM-like lighs higgs boson with mass mh1 *)
*)
mh1=Array2HDM[[All,120]]; (* mass of the light higgs boson mh1 *)
mH=Array2HDM[[All,121]];(* mass of the heavy higgs boson mh2 (always equal to 125GeV) *)
mHpm=Array2HDM[[All,123]]; (* mass of the charged Higgs H\:207a\:207b *)
mA=Array2HDM[[All,122]]; (* mass of the pseudo-scalar Higgs mA *)
m122=Array2HDM[[All,127]]; (* mass m12\.b2 *)

(* others *)
tanB=Array2HDM[[All,128]]; (* values of tan(\[Beta]) *)
sinBA=Array2HDM[[All,124]]; (* values of sin(\[Beta]-\[Alpha]) *)

(* Limits *)
S=Array2HDM[[All,107]];
T=Array2HDM[[All,108]];
U=Array2HDM[[All,109]];
stab=Array2HDM[[All,119]]; (* stability *)
pert=Array2HDM[[All,117]]; (* perturbativity *)
unit=Array2HDM[[All,118]]; (* unitarity *)
\[CapitalDelta]a\[Mu]=Array2HDM[[All,112]]; (* muon anomalous magnetic moment *)
Delta0=Array2HDM[[All,113]];
Bsmumu=Array2HDM[[All,114]]; (* Bs\[Rule]\[Mu]\[Mu] *)
BXsgam=Array2HDM[[All,115]]; (* B\[Rule] Subscript[X, s] \[Gamma] *)
DeltaMd=Array2HDM[[All,116]];
LEPconst=Array2HDM[[All,110]]; (* constraints from the LEP *)
LHCconst=Array2HDM[[All,111]]; (* constraints from the LHC *)

(* Computation of the signal strength Subscript[\[Mu], Subscript[h, 1]\[Rule]\[Gamma]\[Gamma]] *)
muggGamGam=kg2*BRhGaGa/BRhGaGaSM;
muVBFGamGam=kv2*BRhGaGa/BRhGaGaSM;
muBBGamGam=kb2*BRhGaGa/BRhGaGaSM;
muWHGamGam=muVBFGamGam;
muZHGamGam=muVBFGamGam;
muTTHGamGam=kt2*BRhGaGa/BRhGaGaSM;

(* dimensions of Array2HDM *)
lArray=Length[Array2HDM];

Print["Initializations DONE"]




(*(* "Heavy" higgs at 125 GeV *) 

BRhGaGa2=Array2HDM[[All,66]];
BRhgg2=Array2HDM[[All,68]];
BRhWW2=Array2HDM[[All,70]];
BRhZZ2=Array2HDM[[All,72]];
BRhbb2=Array2HDM[[All,74]];
BRh\[Tau]\[Tau]2=Array2HDM[[All,76]];
BRhcc2=Array2HDM[[All,78]];
BRhss2=Array2HDM[[All,80]];
BRh\[Mu]\[Mu]2=Array2HDM[[All,82]];
BRhZ\[Gamma]2=Array2HDM[[All,84]];

BRhGaGa2SM=Array2HDM[[All,46]];
BRhgg2SM=Array2HDM[[All,48]];
BRhWW2SM=Array2HDM[[All,50]];
BRhZZ2SM=Array2HDM[[All,52]];
BRhbb2SM=Array2HDM[[All,54]];
BRh\[Tau]\[Tau]2SM=Array2HDM[[All,56]];
BRhcc2SM=Array2HDM[[All,58]];
BRhss2SM=Array2HDM[[All,60]];
BRh\[Mu]\[Mu]2SM=Array2HDM[[All,62]];
BRhZ\[Gamma]2SM=Array2HDM[[All,64]];

k\[Gamma]22=Array2HDM[[All,86]];
kg22=Array2HDM[[All,88]];
kV22=Array2HDM[[All,90]];
kb22=Array2HDM[[All,94]];
k\[Tau]22=Array2HDM[[All,96]];
kc22=Array2HDM[[All,98]];
ks22=Array2HDM[[All,100]];
k\[Mu]22=Array2HDM[[All,102]];
kZ\[Gamma]22=Array2HDM[[All,104]];
kt22=Array2HDM[[All,106]];
Print["Done !"]

*)


(*chiWW=Block[{a,b,c,mug0,muV0,R}, a=22.53542693935295;
b=2.382300635300683;
c=6.972321500447737;
mug0=Table[1.0199730149215416,lArray];
muV0 = Table[1.4062238304760846,lArray];
(*R=kV22/(k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*BRh\[Tau]\[Tau]2SM+kc22*BRhcc2SM+ks22*BRhss2SM+k\[Mu]22*BRh\[Mu]\[Mu]2SM+kZ\[Gamma]22*BRhZ\[Gamma]2SM);*)
(* without s and some kappas *)
R=kV22/(BRhss2SM+k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*(BRh\[Tau]\[Tau]2SM+BRh\[Mu]\[Mu]2SM)+BRhcc2SM+BRhZ\[Gamma]2SM);
(a*(R*kg22-mug0)^2+2*b*(R*kg22-mug0)*(R*kV22-muV0)+c*(R*kV22-muV0)^2)
];
*)


(*chiZZ=Block[{a,b,c,mug0,muV0,R}, a=11.065910085978555;
b=1.8729205285669688;
c=1.079937160636133;
mug0=Table[1.4388133781554047,lArray];
muV0 = Table[0.8076327247400148,lArray];
(*R=kV22/(k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*BRh\[Tau]\[Tau]2SM+kc22*BRhcc2SM+ks22*BRhss2SM+k\[Mu]22*BRh\[Mu]\[Mu]2SM+kZ\[Gamma]22*BRhZ\[Gamma]2SM);*)
R=kV22/(BRhss2SM+k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*(BRh\[Tau]\[Tau]2SM+BRh\[Mu]\[Mu]2SM)+BRhcc2SM+BRhZ\[Gamma]2SM);
(a*(R*kg22-mug0)^2+2*b*(R*kg22-mug0)*(R*kV22-muV0)+c*(R*kV22-muV0)^2)
];
*)


(*
chi\[Tau]\[Tau]=Block[{a,b,c,mug0,muV0,R}, a=3.496073648864013;
b=2.581303927312333;
c=9.782178330965264;
mug0=Table[1.127715071184127,lArray];
muV0 = Table[1.1391447471474951,lArray];
(*R=k\[Tau]22/(k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*BRh\[Tau]\[Tau]2SM+kc22*BRhcc2SM+ks22*BRhss2SM+k\[Mu]22*BRh\[Mu]\[Mu]2SM+kZ\[Gamma]22*BRhZ\[Gamma]2SM);*)
R=k\[Tau]22/(BRhss2SM+k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*(BRh\[Tau]\[Tau]2SM+BRh\[Mu]\[Mu]2SM)+BRhcc2SM+BRhZ\[Gamma]2SM);
(a*(R*kg22-mug0)^2+2*b*(R*kg22-mug0)*(R*kV22-muV0)+c*(R*kV22-muV0)^2)
];
*)


(*
chi\[Gamma]\[Gamma]=Block[{a,b,c,mug0,muV0,R}, a=16.287512957235908;
b=3.067864180159719;
c=6.175671609346068;
mug0=Table[1.208066041934841,lArray];
muV0 = Table[1.073331195078713,lArray];
(*R=k\[Gamma]22/(k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*BRh\[Tau]\[Tau]2SM+kc22*BRhcc2SM+ks22*BRhss2SM+k\[Mu]22*BRh\[Mu]\[Mu]2SM+kZ\[Gamma]22*BRhZ\[Gamma]2SM);*)
R=k\[Gamma]22/(BRhss2SM+k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*(BRh\[Tau]\[Tau]2SM+BRh\[Mu]\[Mu]2SM)+BRhcc2SM+BRhZ\[Gamma]2SM);
(a*(R*kg22-mug0)^2+2*b*(R*kg22-mug0)*(R*kV22-muV0)+c*(R*kV22-muV0)^2)
];
*)


(*
chibb=Block[{a,b,c,mug0,muV0,R}, a=1.0812183083586497;
b=-0.0018041338479032067;
c=11.364289003398953;
mug0=Table[1.1218581135173131,lArray];
muV0 = Table[0.6628887587435617,lArray];
(*R=kb22/(k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*BRh\[Tau]\[Tau]2SM+kc22*BRhcc2SM+ks22*BRhss2SM+k\[Mu]22*BRh\[Mu]\[Mu]2SM+kZ\[Gamma]22*BRhZ\[Gamma]2SM);*)
R=kb22/(BRhss2SM+k\[Gamma]22*BRhGaGa2SM+kg22*BRhgg2SM+kV22*BRhZZ2SM+kV22*BRhWW2SM+kb22*BRhbb2SM+k\[Tau]22*(BRh\[Tau]\[Tau]2SM+BRh\[Mu]\[Mu]2SM)+BRhcc2SM+BRhZ\[Gamma]2SM);
(a*(R*kg22-mug0)^2+2*b*(R*kg22-mug0)*(R*kV22-muV0)+c*(R*kV22-muV0)^2)
];
*)


(*
chiJB=chiWW+chiZZ+chi\[Tau]\[Tau]+chi\[Gamma]\[Gamma]+chibb
LHCconst

100*(chi-LHCconst)/LHCconst
*)


(* Main program *)


(* Computation of \[Sigma]*BRh\[Gamma]\[Gamma] *)
(* We have computed the signal strenght Subscript[\[Mu], Subscript[h, 1]\[Rule]\[Gamma]\[Gamma]] for different production channels, but for some plots we need \[Sigma]*Subscript[BR, h\[Gamma]\[Gamma]]. In order to do this, we use official data for production and decay width of SM Higgs bosons for different masses avalaible on Atlas or CMS website. We fit them and then compute \[Sigma]*Subscript[BR, h\[Gamma]\[Gamma]](BSM)=Subscript[\[Mu], Subscript[h, 1]\[Rule]\[Gamma]\[Gamma]]*\[Sigma]*Subscript[BR, h\[Gamma]\[Gamma]](SM) *)

(* Fit of the SM cross section at 8TeV *)

(* Load cross sections from LHC official data : Subscript[\[Sigma], 8TeV]=f(Subscript[m, h]) *)
TabGGCrossSection8=Get["cross_section_data/CrossSectionGG8Tev.txt"];
TabVBFCrossSection8=Get["cross_section_data/CrossSectionVBF8Tev.txt"];
TabBBCrossSection8=Get["cross_section_data/CrossSectionBBH8Tev.txt"];
TabWHCrossSection8=Get["cross_section_data/CrossSectionWH8Tev.txt"];
TabZHCrossSection8=Get["cross_section_data/CrossSectionZH8Tev.txt"];
TabTTHCrossSection8=Get["cross_section_data/CrossSectionTTH8Tev.txt"];

(* Cut the data in order to fit the data for masses going from 80 to ~120 GeV *)
TabGGCrossSection8Cut=TabGGCrossSection8[[1;;100,All]]; 
TabVBFCrossSection8Cut=TabVBFCrossSection8[[1;;100,All]];
TabBBCrossSection8Cut=TabBBCrossSection8[[1;;100,All]];
TabWHCrossSection8Cut=TabWHCrossSection8[[1;;100,All]];
TabZHCrossSection8Cut=TabZHCrossSection8[[1;;100,All]];
TabTTHCrossSection8Cut=TabTTHCrossSection8[[1;;100,All]];

(* Fit the data and implement the fit function for each production channel *)
ParameterFitGGCrossSection8=FindFit[TabGGCrossSection8Cut,a+b*x+c*x^2+d*Log[x],{a,b,c,d},x] ;FitGGCrossSection8[x_]:=a+b*x+c*x^2+d*Log[x]/.ParameterFitGGCrossSection8 (* definition of the fit function which can be used later *)

ParameterFitVBFCrossSection8=FindFit[TabVBFCrossSection8Cut,a+b*x+c*x^2+d*Log[x],{a,b,c,d},x];
FitVBFCrossSection8[x_]:=a+b*x+c*x^2+d*Log[x]/.ParameterFitVBFCrossSection8

ParameterFitBBCrossSection8=FindFit[TabBBCrossSection8Cut,a+b*x+c*x^2+d*Log[x],{a,b,c,d},x];
FitBBCrossSection8[x_]:=a+b*x+c*x^2+d*Log[x]/.ParameterFitBBCrossSection8

ParameterFitWHCrossSection8=FindFit[TabWHCrossSection8Cut,a+b*x+c*x^2+d*Log[x],{a,b,c,d},x];
FitWHCrossSection8[x_]:=a+b*x+c*x^2+d*Log[x]/.ParameterFitWHCrossSection8

ParameterFitZHCrossSection8=FindFit[TabZHCrossSection8Cut,a+b*x+c*x^2+d*Log[x],{a,b,c,d},x];
FitZHCrossSection8[x_]:=a+b*x+c*x^2+d*Log[x]/.ParameterFitZHCrossSection8

ParameterFitTTHCrossSection8=FindFit[TabTTHCrossSection8Cut,a+b*x+c*x^2+d*Log[x],{a,b,c,d},x];
FitTTHCrossSection8[x_]:=a+b*x+c*x^2+d*Log[x]/.ParameterFitTTHCrossSection8

(*Plots for check*)
(*
PlotFitGG=Plot[FitGGCrossSection8[x],{x,50,130}];
PlotDataGG=ListPlot[TabGGCrossSection8,PlotRange\[Rule]{{50,120},All}];
Show[PlotDataGG,PlotFitGG]

PlotFitVBF=Plot[FitVBFCrossSection8[x],{x,50,130}];
PlotDataVBF=ListPlot[TabVBFCrossSection8,PlotRange\[Rule]{{50,120},All},PlotStyle\[Rule]Red];
Show[PlotDataVBF,PlotFitVBF]

PlotFitBB=Plot[FitBBCrossSection8[x],{x,50,130},PlotStyle\[Rule]Red];
PlotDataBB=ListPlot[TabBBCrossSection8,PlotRange\[Rule]{{50,120},All}];
Show[PlotDataBB,PlotFitBB]

PlotFitWH=Plot[FitWHCrossSection8[x],{x,50,130},PlotStyle\[Rule]Red];
PlotDataWH=ListPlot[TabWHCrossSection8,PlotRange\[Rule]{{50,120},All}];
Show[PlotDataWH,PlotFitWH]

PlotFitZH=Plot[FitZHCrossSection8[x],{x,50,130},PlotStyle\[Rule]Red];
PlotDataZH=ListPlot[TabZHCrossSection8,PlotRange\[Rule]{{50,120},All}];
Show[PlotDataZH,PlotFitZH]

PlotFitTTH=Plot[FitTTHCrossSection8[x],{x,50,130},PlotStyle\[Rule]Red];
PlotDataTTH=ListPlot[TabTTHCrossSection8,PlotRange\[Rule]{{50,120},All}];
Show[PlotDataTTH,PlotFitTTH]
*)

Print["Production cross section fits DONE"]

(* Computation of the \[Sigma]*BR_BSM *)
(* Now that we have fits functions for \[Sigma] and Subscript[BR, h\[Rule]\[Gamma]\[Gamma]], we can compute the value of \[Sigma]*BR(BSM) *)

(* initialization of the arrays *)
sigmaBRgghGamGam=Array[0,lArray,1]; (* \[Sigma]ggh*BR\[Gamma]\[Gamma]_BSM *)
sigmaBRvbfvhGamGam=Array[0,lArray,1]; (* (\[Sigma]VBF+\[Sigma]VH)*BR\[Gamma]\[Gamma]_BSM *)
sigmaBRhGamGam=Array[0,lArray,1]; (* \[Sigma]_tot*BR\[Gamma]\[Gamma]_BSM *)

sigmaGGHsm=Array[0,lArray,1];(*\[Sigma]GGH in SM*)
sigmaVBFsm=Array[0,lArray,1];(*\[Sigma]VBF in SM*)
sigmaVHsm=Array[0,lArray,1];(*\[Sigma]VH in SM*)

sigmaGGHthdm=Array[0,lArray,1];(*\[Sigma]GGH in THDM*)
sigmaVBFthdm=Array[0,lArray,1];(*\[Sigma]VBF in THDM*)
sigmaVHthdm=Array[0,lArray,1];(*\[Sigma]VH in THDM*)

(* BRhGaGaTHDM=Array[0,lArray,1];(*\BR_h->gamma gamma in THDM*) *)

sigmaTotSM=Array[0,lArray,1]; (* \[Sigma]_tot_SM *)
muTotGamGam=Array[0,lArray,1]; (* \[Mu]\[Gamma]\[Gamma] for all production channels *)
(* We make a loop which computes the value of \[Sigma]*BR_BSM(mh1) = Subscript[\[Mu], gg\[Rule]h\[Rule]\[Gamma]\[Gamma]]*\[Sigma]*BR_SM(mh1) for different production channels. First, we save the mass mh1 (mh1[[iloop]]), then we compute \[Sigma]ggh*BR_BSM thanks to the functions defined above *)

Clear[mloop] (* mass mh1 of a given point *)
For[iloop=1,iloop<=lArray,iloop++,mloop=mh1[[iloop]];sigmaBRgghGamGam[[iloop]]=kg2[[iloop]]*FitGGCrossSection8[mloop]*BRhGaGa[[iloop]];sigmaBRvbfvhGamGam[[iloop]]=kv2[[iloop]]*(FitVBFCrossSection8[mloop]+FitWHCrossSection8[mloop]+FitZHCrossSection8[mloop])*BRhGaGa[[iloop]];sigmaBRhGamGam[[iloop]]=sigmaBRgghGamGam[[iloop]]+sigmaBRvbfvhGamGam[[iloop]]+(kb2[[iloop]]*FitBBCrossSection8[mloop]+kt2[[iloop]]*FitTTHCrossSection8[mloop])*BRhGaGa[[iloop]];sigmaTotSM[[iloop]]=FitGGCrossSection8[mloop]+FitVBFCrossSection8[mloop]+FitBBCrossSection8[mloop]+FitWHCrossSection8[mloop]+FitZHCrossSection8[mloop]+FitTTHCrossSection8[mloop];muTotGamGam[[iloop]]=sigmaBRhGamGam[[iloop]]/(sigmaTotSM[[iloop]]*BRhGaGaSM[[iloop]]);sigmaGGHsm[[iloop]]=FitGGCrossSection8[mloop];sigmaVBFsm[[iloop]]=FitVBFCrossSection8[mloop];sigmaVHsm[[iloop]]=FitWHCrossSection8[mloop]+FitZHCrossSection8[mloop];sigmaGGHthdm[[iloop]]=kg2[[iloop]]*sigmaGGHsm[[iloop]];sigmaVBFthdm[[iloop]]=kv2[[iloop]]*sigmaVBFsm[[iloop]];sigmaVHthdm[[iloop]]=kv2[[iloop]]*sigmaVHsm[[iloop]]]
Print["sigma*BR computation DONE"]






(* convert to GnuPlot files (all data) *)

GnuTabAll=Transpose[{Join[{"mh1"},mh1],Join[{"mH"},mH],Join[{"mA"},mA],Join[{"mHpm"},mHpm],Join[{"kv2"},kv2],Join[{"kt2"},kt2],Join[{"kb2"},kb2],Join[{"tanB"},tanB],Join[{"sinBA"},sinBA],Join[{"muggh"},muggGamGam],Join[{"muVBF"},muVBFGamGam],Join[{"muTot"},muTotGamGam],Join[{"sigmaggh"},sigmaBRgghGamGam],Join[{"sigmaVBF_VF"},sigmaBRvbfvhGamGam],Join[{"sigma_tot_"},sigmaBRhGamGam],Join[{"m12"},Sqrt[m122]],Join[{"sigma_ggh"},sigmaGGHthdm],Join[{"BRhGaGa"},BRhGaGa]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/All_points.txt"],GnuTabAll,"Table"]
Print["Conversion Done"]




(* now that we have computed all the variables we need, we can check the different limits (EW, theoretical, LEP, LHC constraints) that we can put on the set of points ! *)


(* points passing EW and theoretical constraints *)

(* Steps :
1.For each points generated by 2HDMC and loaded through Array2HDM, we check if they pass EW and theoretical constraints. We fill a preliminary array (TabTH) with 1 if the point is ok, 0 else. We also count the number of points passing these constraints in order to later generate arrays of the right size containing only points passing EW and TH constraints.
2. We initialize the different arrays containing only points passing the constraints.
3. We make a loop in order to fill these arrays with the appropriate data. *)

TabTH=Array[0,lArray,1]; (* array with 1 if the point pass TH constraints, 0 else *)
Clear[sizeTH];
sizeTH=0; (* count the number of points passing theoretical constraints *)

(* loop which fills the array defined above *)
For[iloop=1,iloop<=lArray,iloop++,If[(limS-2*limSDown)<=S[[iloop]]<=(limS+2*limSUp)&&(limT-2*limTDown)<=T[[iloop]]<=(limT+2*limTUp)&&(limU-2*limUDown)<= U[[iloop]]<=(limU+2*limUUp)&&stab[[iloop]]==1&&pert[[iloop]]==1&&unit[[iloop]]==1&&(limBXsgam-2*limBXsgamDown)<=BXsgam[[iloop]]<=(limBXsgam +2*limBXsgamUp)&&(limDeltaMd-2*limDeltaMdDown)<=DeltaMd[[iloop]]<=(limDeltaMd+2*limDeltaMdUp)&&(limDelta0-2*limDelta0Down)<=Delta0[[iloop]]<=(limDelta0+2*limDelta0Up)&&(limBsmumu-2*limBsmumuDown)<=Bsmumu[[iloop]]<=(limBsmumu+2*limBsmumuUp)(*&&17.5*10^-10\[LessEqual]Abs[\[CapitalDelta]a\[Mu][[iloop]]]\[LessEqual]33.5*10^-10*),TabTH[[iloop]]=1;sizeTH++,TabTH[[iloop]]=0]]

Print["Loop for sizeTH DONE"]
(* initialization of the mass arrays. They will be filled with the value of the mass for a given point if this point passes TH constraints *)
mh1TH=Array[0,sizeTH,1];
mHTH=Array[0,sizeTH,1];
mHpmTH=Array[0,sizeTH,1];
mATH=Array[0,sizeTH,1];
m122TH=Array[0,sizeTH,1];

(* \[Mu] arrays *)
muggGamGamTH=Array[0,sizeTH,1]; 
muVBFGamGamTH=Array[0,sizeTH,1];
muTotGamGamTH=Array[0,sizeTH,1];

(* \[Kappa] arrays *)
kvTH=Array[0,sizeTH,1]; 
ktTH=Array[0,sizeTH,1]; 
kbTH=Array[0,sizeTH,1]; 

(* angles *)
tanBTH=Array[0,sizeTH,1];
sinBATH=Array[0,sizeTH,1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamTH=Array[0,sizeTH,1];
sigmaBRvbfvhGamGamTH=Array[0,sizeTH,1];
sigmaBRhGamGamTH=Array[0,sizeTH,1];

(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop]
jloop=0; (* index of the "arraysTH" *)

For[iloop=1,iloop<=lArray,iloop++,If[TabTH[[iloop]]==1,jloop++;muggGamGamTH[[jloop]]=muggGamGam[[iloop]];muVBFGamGamTH[[jloop]]=muVBFGamGam[[iloop]];kvTH[[jloop]]=kv2[[iloop]];ktTH[[jloop]]=kt2[[iloop]];kbTH[[jloop]]=kb2[[iloop]];mh1TH[[jloop]]=mh1[[iloop]];mHTH[[jloop]]=mH[[iloop]];mHpmTH[[jloop]]=mHpm[[iloop]];mATH[[jloop]]=mA[[iloop]];
m122TH[[jloop]]=m122[[iloop]];tanBTH[[jloop]]=tanB[[iloop]];sinBATH[[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamTH[[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamTH[[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamTH[[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamTH[[jloop]]=muTotGamGam[[iloop]]]]

If[jloop!=sizeTH,Print["Error in the tab for TH constraints !"]]

Print["EW constraints DONE"]






Print[sizeTH]



(* convert to GnuPlot files *)

GnuTabTH=Transpose[{Join[{"mh1TH"},mh1TH],Join[{"mHTH"},mHTH],Join[{"mATH"},mATH],Join[{"mHpmTH"},mHpmTH],Join[{"kvTH"},kvTH],Join[{"ktTH"},ktTH],Join[{"kbTH"},kbTH],Join[{"tanBTH"},tanBTH],Join[{"sinBATH"},sinBATH],Join[{"mugghTH"},muggGamGamTH],Join[{"muVBFTH"},muVBFGamGamTH],Join[{"muTotTH"},muTotGamGamTH],Join[{"sigmagghTH"},sigmaBRgghGamGamTH],Join[{"sigmaVBF_VFTH"},sigmaBRvbfvhGamGamTH],Join[{"sigma_tot_TH"},sigmaBRhGamGamTH],Join[{"m12TH"},Sqrt[m122TH]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_constraints.txt"],GnuTabTH,"Table"]
Print["Conversion Done"]



(* points passing LEP constraints *)
(* We do exactly the same things as before, with the same steps. *)

(* array with 1 if the point pass TH AND LEP constraints, 0 else *)
TabLEP=Array[0,lArray,1]; 

Clear[sizeLEP];
sizeLEP=0; (* count the number of points passing LEP constraints *)

(* loop which fills the array defined above : TabLep =1  if and only if the point iloop passes TH constraints (TabTH==1) AND passes LEP constraints (LEPconst\[Equal]1) *)
For[iloop=1,iloop<=lArray,iloop++,If[TabTH[[iloop]]==1&&LEPconst[[iloop]]==1,TabLEP[[iloop]]=1;sizeLEP++,TabLEP[[iloop]]=0]]

Print["Loop for sizeLEP DONE"]

(* initialisation of the mass array *)
mh1LEP=Array[0,sizeLEP,1];
mHLEP=Array[0,sizeLEP,1];
mHpmLEP=Array[0,sizeLEP,1];
mALEP=Array[0,sizeLEP,1];
m122LEP=Array[0,sizeLEP,1];

(* \[Mu] arrays *)
muggGamGamLEP=Array[0,sizeLEP,1]; 
muVBFGamGamLEP=Array[0,sizeLEP,1];
muTotGamGamLEP=Array[0,sizeLEP,1];
(* \[Kappa] arrays *)
kvLEP=Array[0,sizeLEP,1]; 
ktLEP=Array[0,sizeLEP,1]; 
kbLEP=Array[0,sizeLEP,1]; 

(* angles *)
tanBLEP=Array[0,sizeLEP,1];
sinBALEP=Array[0,sizeLEP,1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamLEP=Array[0,sizeLEP,1];
sigmaBRvbfvhGamGamLEP=Array[0,sizeLEP,1];
sigmaBRhGamGamLEP=Array[0,sizeLEP,1];
(* Loop which fills the arrays defined just before with the points passing TH and LEP constraints*)
Clear[jloop]
jloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabLEP[[iloop]]==1,jloop++;muggGamGamLEP[[jloop]]=muggGamGam[[iloop]];muVBFGamGamLEP[[jloop]]=muVBFGamGam[[iloop]];kvLEP[[jloop]]=kv2[[iloop]];ktLEP[[jloop]]=kt2[[iloop]];kbLEP[[jloop]]=kb2[[iloop]];mh1LEP[[jloop]]=mh1[[iloop]];mHLEP[[jloop]]=mH[[iloop]];mHpmLEP[[jloop]]=mHpm[[iloop]];mALEP[[jloop]]=mA[[iloop]];m122LEP[[jloop]]=m122[[iloop]];tanBLEP[[jloop]]=tanB[[iloop]];sinBALEP[[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamLEP[[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamLEP[[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamLEP[[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamLEP[[jloop]]=muTotGamGam[[iloop]]]]

If[jloop!=sizeLEP,Print["Error in the tab for TH constraints !"]]

Print["LEP constraintes DONE"]





Print[sizeLEP]



(* convert to GnuPlot files *)

GnuTabLEP=Transpose[{Join[{"mh1LEP"},mh1LEP],Join[{"mHLEP"},mHLEP],Join[{"mALEP"},mALEP],Join[{"mHpmLEP"},mHpmLEP],Join[{"kvLEP"},kvLEP],Join[{"ktLEP"},ktLEP],Join[{"kbLEP"},kbLEP],Join[{"tanBLEP"},tanBLEP],Join[{"sinBALEP"},sinBALEP],Join[{"mugghLEP"},muggGamGamLEP],Join[{"muVBFLEP"},muVBFGamGamLEP],Join[{"muTotLEP"},muTotGamGamLEP],Join[{"sigmagghLEP"},sigmaBRgghGamGamLEP],Join[{"sigmaVBF_VFLEP"},sigmaBRvbfvhGamGamLEP],Join[{"sigma_tot_LEP"},sigmaBRhGamGamLEP],Join[{"m12LEP"},Sqrt[m122LEP]]}];
Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/LEP_constraints.txt"],GnuTabLEP,"Table"]
Print["Conversion Done"]



(* points passing LHC constraints *)

(* Again the same method... *)


(* Array with 1 if the point pass TH, LEP AND LHC constraints, 0 else *)
TabLHC=Array[0,lArray,1]; 

Clear[sizeLHC];
sizeLHC=0; (* count the number of points passing theoretical constraints *)
(* loop which fills the array defined above : TabLHC =1 if and only if the point iloop pass TH and LEP constraints (TabLEP\[NotEqual]0 ) and LHC constraints (LHCconst\[LessEqual]12.89)*)
For[iloop=1,iloop<=lArray,iloop++,If[TabLEP[[iloop]]==1&&LHCconst[[iloop]]<=boundLHC,TabLHC[[iloop]]=1;sizeLHC++,TabLHC[[iloop]]=0]]

Print["Loop for sizeLHC DONE"]

(* initialisation of the mass arrays *)
mh1LHC=Array[0,sizeLHC,1];
mHLHC=Array[0,sizeLHC,1];
mHpmLHC=Array[0,sizeLHC,1];
mALHC=Array[0,sizeLHC,1];
m122LHC=Array[0,sizeLHC,1];

(* \[Mu] arrays *)
muggGamGamLHC=Array[0,sizeLHC,1]; 
muVBFGamGamLHC=Array[0,sizeLHC,1];
muTotGamGamLHC=Array[0,sizeLHC,1];

(* \[Kappa] arrays *)
kvLHC=Array[0,sizeLHC,1]; 
ktLHC=Array[0,sizeLHC,1]; 
kbLHC=Array[0,sizeLHC,1]; 

(* angles *)
tanBLHC=Array[0,sizeLHC,1];
sinBALHC=Array[0,sizeLHC,1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamLHC=Array[0,sizeLHC,1];
sigmaBRvbfvhGamGamLHC=Array[0,sizeLHC,1];
sigmaBRhGamGamLHC=Array[0,sizeLHC,1];

(* Loop which fill the arrays defined just before with the points passing TH and LEP constraints*)
Clear[jloop]
jloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabLHC[[iloop]]==1,jloop++;muggGamGamLHC[[jloop]]=muggGamGam[[iloop]];muVBFGamGamLHC[[jloop]]=muVBFGamGam[[iloop]];kvLHC[[jloop]]=kv2[[iloop]];ktLHC[[jloop]]=kt2[[iloop]];kbLHC[[jloop]]=kb2[[iloop]];mh1LHC[[jloop]]=mh1[[iloop]];mHLHC[[jloop]]=mH[[iloop]];mHpmLHC[[jloop]]=mHpm[[iloop]];mALHC[[jloop]]=mA[[iloop]];m122LHC[[jloop]]=m122[[iloop]];tanBLHC[[jloop]]=tanB[[iloop]];sinBALHC[[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamLHC[[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamLHC[[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamLHC[[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamLHC[[jloop]]=muTotGamGam[[iloop]]]]

If[jloop!=sizeLHC,Print["Error in the tab for TH constraints !"]]
Print["LHC constraints DONE"]

Print["END !"]




Print[sizeLHC]



(* convert to GnuPlot files *)

GnuTabLHC=Transpose[{Join[{"mh1LHC"},mh1LHC],Join[{"mHLHC"},mHLHC],Join[{"mALHC"},mALHC],Join[{"mHpmLHC"},mHpmLHC],Join[{"kvLHC"},kvLHC],Join[{"ktLHC"},ktLHC],Join[{"kbLHC"},kbLHC],Join[{"tanBLHC"},tanBLHC],Join[{"sinBALHC"},sinBALHC],Join[{"mugghLHC"},muggGamGamLHC],Join[{"muVBFLHC"},muVBFGamGamLHC],Join[{"muTotLHC"},muTotGamGamLHC],Join[{"sigmagghLHC"},sigmaBRgghGamGamLHC],Join[{"sigmaVBF_VFLHC"},sigmaBRvbfvhGamGamLHC],Join[{"sigma_tot_LHC"},sigmaBRhGamGamLHC],Join[{"m12LHC"},Sqrt[m122LHC]]}];
Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/LHC_constraints.txt"],GnuTabLHC,"Table"]
Print["Conversion Done"]



(* Histogram *)
(* We want to see how is distributed the ratio \[Sigma]_ggh/\[Sigma]_VBF and see if it is different from the one of the SM *)

(*We want to see how is distributed the ratio \[Sigma]_ggh/\[Sigma]_VBF and see if it is different from the one of the SM*)(*We can obtain the cross section ratio through (\[Sigma]_ggh/\[Sigma]_VBF)BSM=\[Kappa]g\.b2/\[Kappa]V\.b2*(\[Sigma]_ggh/\[Sigma]_VBF)SM and through (\[Sigma]_ggh/\[Sigma]_VBF)BSM=\[Sigma]_ggh*BR(h\[Rule]\[Gamma]\[Gamma])/(\[Sigma]_VBF*BR(h\[Rule]\[Gamma]\[Gamma])).The two methods give the same results.*)
(*
ratioSM=sigmaGGHsm/(sigmaVBFsm+sigmaVHsm);

ratioBSM=sigmaBRgghGamGam/sigmaBRvbfvhGamGam;

ratioBSMLHC=sigmaBRgghGamGamLHC/sigmaBRvbfvhGamGamLHC;

arrayRab=Table[-5,lArray];
arrayRabLHC=Table[-5,sizeLHC];
*)



(*convert to GnuPlot files*)

(*GnuTabRatioTHDM=Transpose[{Join[{"ratio BSM : sigmaggh/sigmaVBF (all points)"},ratioBSM],Join[{"ratio BSM : sigmaVBF/sigmaggh"},1/ratioBSM],Join[{"ratio SM"},ratioSM],Join[{"mh1"},mh1],Join[{"gnuTrick"},arrayRab]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/ratio_cross_sections_THDM.txt"],GnuTabRatioTHDM,"Table"]
*)


(*
GnuTabRatioTHDMLHC=Transpose[{Join[{"ratio BSM : sigmaggh/sigmaVBF (LHC points)"},ratioBSMLHC],Join[{"ratio BSM : sigmaVBF/sigmaggh"},1/ratioBSMLHC],Join[{"mh1"},mh1LHC],Join[{"gnuTrick"},arrayRabLHC]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/ratio_cross_sections_THDM_LHC.txt"],GnuTabRatioTHDMLHC,"Table"]

Print["Conversion Done"]
*)



(* Details on the different theoretical limits *)




(* points passing S,T,U constraints *)

TabSTU=Array[0,lArray,1]; (* array with 1 if the point pass TH constraints, 0 else *)

Clear[sizeSTU];
sizeSTU=0; (* count the number of points passing STU constraints *)

(* loop which fills the array defined above with 1 if it passes the constraints, 0 else ; limit at 2\[Sigma] *)
For[iloop=1,iloop<=lArray,iloop++,If[(limS-2*limSDown)<=S[[iloop]]<=(limS+2*limSUp)&&(limT-2*limTDown)<=T[[iloop]]<=(limT+2*limTUp)&&(limU-2*limUDown)<= U[[iloop]]<=(limU+2*limUUp),TabSTU[[iloop]]=1;sizeSTU++,TabSTU[[iloop]]=0]]

Print["Loop for sizeSTU DONE"]

(* initialisation of the mass arrays *)
mh1STU=Array[0,sizeSTU,1];
mHSTU=Array[0,sizeSTU,1];
mHpmSTU=Array[0,sizeSTU,1];
mASTU=Array[0,sizeSTU,1];
m122STU=Array[0,sizeSTU,1];

(* \[Mu] arrays *)
muggGamGamSTU=Array[0,sizeSTU,1]; 
muVBFGamGamSTU=Array[0,sizeSTU,1];
muTotGamGamSTU=Array[0,sizeSTU,1];

(* \[Kappa] arrays *)
kvSTU=Array[0,sizeSTU,1]; 
ktSTU=Array[0,sizeSTU,1]; 
kbSTU=Array[0,sizeSTU,1]; 

(* angles *)
tanBSTU=Array[0,sizeSTU,1];
sinBASTU=Array[0,sizeSTU,1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamSTU=Array[0,sizeSTU,1];
sigmaBRvbfvhGamGamSTU=Array[0,sizeSTU,1];
sigmaBRhGamGamSTU=Array[0,sizeSTU,1];
(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop]
jloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabSTU[[iloop]]==1,jloop++;muggGamGamSTU[[jloop]]=muggGamGam[[iloop]];muVBFGamGamSTU[[jloop]]=muVBFGamGam[[iloop]];kvSTU[[jloop]]=kv2[[iloop]];ktSTU[[jloop]]=kt2[[iloop]];kbSTU[[jloop]]=kb2[[iloop]];mh1STU[[jloop]]=mh1[[iloop]];mHSTU[[jloop]]=mH[[iloop]];mHpmSTU[[jloop]]=mHpm[[iloop]];mASTU[[jloop]]=mA[[iloop]];m122STU[[jloop]]=m122[[iloop]];tanBSTU[[jloop]]=tanB[[iloop]];sinBASTU[[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamSTU[[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamSTU[[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamSTU[[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamSTU[[jloop]]=muTotGamGam[[iloop]]]]

Print["Loop for STU DONE"]




sizeSTU


(* convert to GnuPlot files *)
GnuTabSTU=Transpose[{Join[{"mh1STU"},mh1STU],Join[{"mHSTU"},mHSTU],Join[{"mASTU"},mASTU],Join[{"mHpmSTU"},mHpmSTU],Join[{"kvSTU"},kvSTU],Join[{"ktSTU"},ktSTU],Join[{"kbSTU"},kbSTU],Join[{"tanBSTU"},tanBSTU],Join[{"sinBASTU"},sinBASTU],Join[{"mugghSTU"},muggGamGamSTU],Join[{"muVBFSTU"},muVBFGamGamSTU],Join[{"muTotSTU"},muTotGamGamSTU],Join[{"sigmagghSTU"},sigmaBRgghGamGamSTU],Join[{"sigmaVBF_VFSTU"},sigmaBRvbfvhGamGamSTU],Join[{"sigma_tot_STU"},sigmaBRhGamGamSTU],Join[{"m12STU"},Sqrt[m122STU]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_STU_limits.txt"],GnuTabSTU,"Table"]

Print["Conversion Done"]




(* points passing pert, unit, stab constraints *)

TabStab=Array[0,lArray,1]; (* array with 1 if the point pass TH AND pert., unit., stab. constraints, 0 else *)

Clear[sizeStab];
sizeStab=0; (* count the number of points passing Stab constraints *)

(* loop which fills the array defined above *)
For[iloop=1,iloop<=lArray,iloop++,If[TabSTU[[iloop]]==1&&stab[[iloop]]==1&&pert[[iloop]]==1&&unit[[iloop]]==1,TabStab[[iloop]]=1;sizeStab++,TabStab[[iloop]]=0]]

Print["Loop for sizeStab DONE"]

(* initialisation of the mass arrays *)
mh1Stab=Array[0,sizeStab,1];
mHStab=Array[0,sizeStab,1];
mHpmStab=Array[0,sizeStab,1];
mAStab=Array[0,sizeStab,1];
m122Stab=Array[0,sizeStab,1];

(* \[Mu] arrays *)
muggGamGamStab=Array[0,sizeStab,1]; 
muVBFGamGamStab=Array[0,sizeStab,1];
muTotGamGamStab=Array[0,sizeStab,1];

(* \[Kappa] arrays *)
kvStab=Array[0,sizeStab,1]; 
ktStab=Array[0,sizeStab,1]; 
kbStab=Array[0,sizeStab,1]; 

(* angles *)
tanBStab=Array[0,sizeStab,1];
sinBAStab=Array[0,sizeStab,1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamStab=Array[0,sizeStab,1];
sigmaBRvbfvhGamGamStab=Array[0,sizeStab,1];
sigmaBRhGamGamStab=Array[0,sizeStab,1];
(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop]
jloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabStab[[iloop]]==1,jloop++;muggGamGamStab[[jloop]]=muggGamGam[[iloop]];muVBFGamGamStab[[jloop]]=muVBFGamGam[[iloop]];kvStab[[jloop]]=kv2[[iloop]];ktStab[[jloop]]=kt2[[iloop]];kbStab[[jloop]]=kb2[[iloop]];mh1Stab[[jloop]]=mh1[[iloop]];mHStab[[jloop]]=mH[[iloop]];mHpmStab[[jloop]]=mHpm[[iloop]];mAStab[[jloop]]=mA[[iloop]];m122Stab[[jloop]]=m122[[iloop]];tanBStab[[jloop]]=tanB[[iloop]];sinBAStab[[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamStab[[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamStab[[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamStab[[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamStab[[jloop]]=muTotGamGam[[iloop]]]]
Print["Stab, unit, pert DONE"]



sizeStab


(* convert to GnuPlot files *)

GnuTabStab=Transpose[{Join[{"mh1Stab"},mh1Stab],Join[{"mHStab"},mHStab],Join[{"mAStab"},mAStab],Join[{"mHpmStab"},mHpmStab],Join[{"kvStab"},kvStab],Join[{"ktStab"},ktStab],Join[{"kbStab"},kbStab],Join[{"tanBStab"},tanBStab],Join[{"sinBAStab"},sinBAStab],Join[{"mugghStab"},muggGamGamStab],Join[{"muVBFStab"},muVBFGamGamStab],Join[{"muTotStab"},muTotGamGamStab],Join[{"sigmagghStab"},sigmaBRgghGamGamStab],Join[{"sigmaVBF_VFStab"},sigmaBRvbfvhGamGamStab],Join[{"sigma_tot_Stab"},sigmaBRhGamGamStab],Join[{"m12Stab"},Sqrt[m122Stab]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_StabPertUnit_limits.txt"],GnuTabStab,"Table"]

Print["Conversion Done"]


(* points passing B\[Rule]Xs \[Gamma] constraints (+ Stab, stab, pert, unit) *)

(* This time, we will have two arrays : one for points passing the constraints at 1\[Sigma], the second for points passing the constraints at 2\[Sigma] *)
TabBxs2\[Sigma]=Array[0,lArray,1]; 
TabBxs1\[Sigma]=Array[0,lArray,1];(* array with 1 if the point pass TH constraints, 0 else *)
Clear[sizeBxs1\[Sigma],sizeBxs2\[Sigma]];
sizeBxs1\[Sigma]=0; (* count the number of points passing Stab constraints *)
sizeBxs2\[Sigma]=0;
(* loop which fills the array defined above *)
For[iloop=1,iloop<=lArray,iloop++,If[TabStab[[iloop]]==1,If[(limBXsgam-limBXsgamDown)<=BXsgam[[iloop]]<=(limBXsgam+limBXsgamUp),TabBxs1\[Sigma][[iloop]]=1;sizeBxs1\[Sigma]++;TabBxs2\[Sigma][[iloop]]=1;sizeBxs2\[Sigma]++,If[(limBXsgam-2*limBXsgamDown)<=BXsgam[[iloop]]<=(limBXsgam+2*limBXsgamUp),TabBxs2\[Sigma][[iloop]]=1;sizeBxs2\[Sigma]++,TabBxs2\[Sigma][[iloop]]=0]]]]
Print["Loop for sizeBxs DONE"]
(* initialisation of the mass arrays *)
mh1Bxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
mHBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
mHpmBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
mABxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
m122Bxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];

mh1Bxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
mHBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
mHpmBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
mABxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
m122Bxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];

(* \[Mu] arrays *)
muggGamGamBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1]; 
muVBFGamGamBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
muTotGamGamBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];

muggGamGamBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1]; 
muVBFGamGamBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
muTotGamGamBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];

(* \[Kappa] arrays *)
kvBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1]; 
ktBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1]; 
kbBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1]; 

kvBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1]; 
ktBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1]; 
kbBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];

(* angles *)
tanBBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
sinBABxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];

tanBBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
sinBABxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
sigmaBRvbfvhGamGamBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];
sigmaBRhGamGamBxs1\[Sigma]=Array[0,sizeBxs1\[Sigma],1];

sigmaBRgghGamGamBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
sigmaBRvbfvhGamGamBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];
sigmaBRhGamGamBxs2\[Sigma]=Array[0,sizeBxs2\[Sigma],1];

(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop,kloop]
jloop=0;
kloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabBxs1\[Sigma][[iloop]]==1,jloop++;muggGamGamBxs1\[Sigma][[jloop]]=muggGamGam[[iloop]];muVBFGamGamBxs1\[Sigma][[jloop]]=muVBFGamGam[[iloop]];kvBxs1\[Sigma][[jloop]]=kv2[[iloop]];ktBxs1\[Sigma][[jloop]]=kt2[[iloop]];kbBxs1\[Sigma][[jloop]]=kb2[[iloop]];mh1Bxs1\[Sigma][[jloop]]=mh1[[iloop]];mHBxs1\[Sigma][[jloop]]=mH[[iloop]];mHpmBxs1\[Sigma][[jloop]]=mHpm[[iloop]];mABxs1\[Sigma][[jloop]]=mA[[iloop]];m122Bxs1\[Sigma][[jloop]]=m122[[iloop]];tanBBxs1\[Sigma][[jloop]]=tanB[[iloop]];sinBABxs1\[Sigma][[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamBxs1\[Sigma][[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamBxs1\[Sigma][[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamBxs1\[Sigma][[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamBxs1\[Sigma][[jloop]]=muTotGamGam[[iloop]]];
If[TabBxs2\[Sigma][[iloop]]==1,kloop++;muggGamGamBxs2\[Sigma][[kloop]]=muggGamGam[[iloop]];muVBFGamGamBxs2\[Sigma][[kloop]]=muVBFGamGam[[iloop]];kvBxs2\[Sigma][[kloop]]=kv2[[iloop]];ktBxs2\[Sigma][[kloop]]=kt2[[iloop]];kbBxs2\[Sigma][[kloop]]=kb2[[iloop]];mh1Bxs2\[Sigma][[kloop]]=mh1[[iloop]];mHBxs2\[Sigma][[kloop]]=mH[[iloop]];mHpmBxs2\[Sigma][[kloop]]=mHpm[[iloop]];mABxs2\[Sigma][[kloop]]=mA[[iloop]];m122Bxs2\[Sigma][[kloop]]=m122[[iloop]];tanBBxs2\[Sigma][[kloop]]=tanB[[iloop]];sinBABxs2\[Sigma][[kloop]]=sinBA[[iloop]];sigmaBRgghGamGamBxs2\[Sigma][[kloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamBxs2\[Sigma][[kloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamBxs2\[Sigma][[kloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamBxs2\[Sigma][[kloop]]=muTotGamGam[[iloop]]]]

Print["Loop for Bxs DONE"]




sizeBxs2\[Sigma]


(* convert to GnuPlot files *)

GnuTabBxs1\[Sigma]=Transpose[{Join[{"mh1Bxs1\[Sigma]"},mh1Bxs1\[Sigma]],Join[{"mHBxs1\[Sigma]"},mHBxs1\[Sigma]],Join[{"mABxs1\[Sigma]"},mABxs1\[Sigma]],Join[{"mHpmBxs1\[Sigma]"},mHpmBxs1\[Sigma]],Join[{"kvBxs1\[Sigma]"},kvBxs1\[Sigma]],Join[{"ktBxs1\[Sigma]"},ktBxs1\[Sigma]],Join[{"kbBxs1\[Sigma]"},kbBxs1\[Sigma]],Join[{"tanBBxs1\[Sigma]"},tanBBxs1\[Sigma]],Join[{"sinBABxs1\[Sigma]"},sinBABxs1\[Sigma]],Join[{"mugghBxs1\[Sigma]"},muggGamGamBxs1\[Sigma]],Join[{"muVBFBxs1\[Sigma]"},muVBFGamGamBxs1\[Sigma]],Join[{"muTotBxs1\[Sigma]"},muTotGamGamBxs1\[Sigma]],Join[{"sigmagghBxs1\[Sigma]"},sigmaBRgghGamGamBxs1\[Sigma]],Join[{"sigmaVBF_VFBxs1\[Sigma]"},sigmaBRvbfvhGamGamBxs1\[Sigma]],Join[{"sigma_tot_Bxs1\[Sigma]"},sigmaBRhGamGamBxs1\[Sigma]],Join[{"m12Bxs1\[Sigma]"},Sqrt[m122Bxs1\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_Bxs_1sigma_limits.txt"],GnuTabBxs1\[Sigma],"Table"]

GnuTabBxs2\[Sigma]=Transpose[{Join[{"mh1Bxs2\[Sigma]"},mh1Bxs2\[Sigma]],Join[{"mHBxs2\[Sigma]"},mHBxs2\[Sigma]],Join[{"mABxs2\[Sigma]"},mABxs2\[Sigma]],Join[{"mHpmBxs2\[Sigma]"},mHpmBxs2\[Sigma]],Join[{"kvBxs2\[Sigma]"},kvBxs2\[Sigma]],Join[{"ktBxs2\[Sigma]"},ktBxs2\[Sigma]],Join[{"kbBxs2\[Sigma]"},kbBxs2\[Sigma]],Join[{"tanBBxs2\[Sigma]"},tanBBxs2\[Sigma]],Join[{"sinBABxs2\[Sigma]"},sinBABxs2\[Sigma]],Join[{"mugghBxs2\[Sigma]"},muggGamGamBxs2\[Sigma]],Join[{"muVBFBxs2\[Sigma]"},muVBFGamGamBxs2\[Sigma]],Join[{"muTotBxs2\[Sigma]"},muTotGamGamBxs2\[Sigma]],Join[{"sigmagghBxs2\[Sigma]"},sigmaBRgghGamGamBxs2\[Sigma]],Join[{"sigmaVBF_VFBxs2\[Sigma]"},sigmaBRvbfvhGamGamBxs2\[Sigma]],Join[{"sigma_tot_Bxs2\[Sigma]"},sigmaBRhGamGamBxs2\[Sigma]],Join[{"m12Bxs2\[Sigma]"},Sqrt[m122Bxs2\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_Bxs_2sigma_limits.txt"],GnuTabBxs2\[Sigma],"Table"]
Print["Conversion Done"]




(* points passing Bs\[Rule]\[Mu]\[Mu] constraints (+ Stab, stab. unit.pert., B\[Rule]Xs \[Gamma] ) *)

(* same thing : we consider limits at 1\[Sigma] and 2\[Sigma] *)
TabBs2\[Sigma]=Array[0,lArray,1]; 
TabBs1\[Sigma]=Array[0,lArray,1];(* array with 1 if the point pass TH constraints, 0 else *)

Clear[sizeBs1\[Sigma],sizeBs2\[Sigma]];
sizeBs1\[Sigma]=0; (* count the number of points passing Stab constraints *)
sizeBs2\[Sigma]=0;

(* loop which fills the array defined above *)
For[iloop=1,iloop<=lArray,iloop++,If[TabBxs1\[Sigma][[iloop]]==1&&(limBsmumu-limBsmumuDown)<=Bsmumu[[iloop]]<=(limBsmumu+limBsmumuUp),TabBs1\[Sigma][[iloop]]=1;sizeBs1\[Sigma]++;TabBs2\[Sigma][[iloop]]=1;sizeBs2\[Sigma]++];If[(TabBxs1\[Sigma][[iloop]]==1||TabBxs2\[Sigma][[iloop]]==1)&&(limBsmumu-2*limBsmumuDown)<=Bsmumu[[iloop]]<=(limBsmumu+2*limBsmumuUp),TabBs2\[Sigma][[iloop]]=1;sizeBs2\[Sigma]++]]

Print["Loop for sizeBs DONE"]
(* initialisation of the mass arrays *)
mh1Bs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
mHBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
mHpmBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
mABs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
m122Bs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];

mh1Bs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
mHBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
mHpmBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
mABs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
m122Bs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];

(* \[Mu] arrays *)
muggGamGamBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1]; 
muVBFGamGamBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
muTotGamGamBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];

muggGamGamBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1]; 
muVBFGamGamBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
muTotGamGamBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];

(* \[Kappa] arrays *)
kvBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1]; 
ktBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1]; 
kbBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1]; 

kvBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1]; 
ktBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1]; 
kbBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];

(* angles *)
tanBBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
sinBABs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];

tanBBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
sinBABs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
sigmaBRvbfvhGamGamBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];
sigmaBRhGamGamBs1\[Sigma]=Array[0,sizeBs1\[Sigma],1];

sigmaBRgghGamGamBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
sigmaBRvbfvhGamGamBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
sigmaBRhGamGamBs2\[Sigma]=Array[0,sizeBs2\[Sigma],1];
(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop,kloop]
jloop=0;
kloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabBs1\[Sigma][[iloop]]==1,jloop++;muggGamGamBs1\[Sigma][[jloop]]=muggGamGam[[iloop]];muVBFGamGamBs1\[Sigma][[jloop]]=muVBFGamGam[[iloop]];kvBs1\[Sigma][[jloop]]=kv2[[iloop]];ktBs1\[Sigma][[jloop]]=kt2[[iloop]];kbBs1\[Sigma][[jloop]]=kb2[[iloop]];mh1Bs1\[Sigma][[jloop]]=mh1[[iloop]];mHBs1\[Sigma][[jloop]]=mH[[iloop]];mHpmBs1\[Sigma][[jloop]]=mHpm[[iloop]];mABs1\[Sigma][[jloop]]=mA[[iloop]];m122Bs1\[Sigma][[jloop]]=m122[[iloop]];tanBBs1\[Sigma][[jloop]]=tanB[[iloop]];sinBABs1\[Sigma][[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamBs1\[Sigma][[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamBs1\[Sigma][[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamBs1\[Sigma][[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamBs1\[Sigma][[jloop]]=muTotGamGam[[iloop]]];
If[TabBs2\[Sigma][[iloop]]==1,kloop++;muggGamGamBs2\[Sigma][[kloop]]=muggGamGam[[iloop]];muVBFGamGamBs2\[Sigma][[kloop]]=muVBFGamGam[[iloop]];kvBs2\[Sigma][[kloop]]=kv2[[iloop]];ktBs2\[Sigma][[kloop]]=kt2[[iloop]];kbBs2\[Sigma][[kloop]]=kb2[[iloop]];mh1Bs2\[Sigma][[kloop]]=mh1[[iloop]];mHBs2\[Sigma][[kloop]]=mH[[iloop]];mHpmBs2\[Sigma][[kloop]]=mHpm[[iloop]];mABs2\[Sigma][[kloop]]=mA[[iloop]];m122Bs2\[Sigma][[kloop]]=m122[[iloop]];tanBBs2\[Sigma][[kloop]]=tanB[[iloop]];sinBABs2\[Sigma][[kloop]]=sinBA[[iloop]];sigmaBRgghGamGamBs2\[Sigma][[kloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamBs2\[Sigma][[kloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamBs2\[Sigma][[kloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamBs2\[Sigma][[kloop]]=muTotGamGam[[iloop]]]]

Print["Loop for Bs DONE"]



sizeBxs2\[Sigma]
sizeBs2\[Sigma]
sizeBs1\[Sigma]
sizeTH


(* convert to GnuPlot files *)

GnuTabBs1\[Sigma]=Transpose[{Join[{"mh1Bs1\[Sigma]"},mh1Bs1\[Sigma]],Join[{"mHBs1\[Sigma]"},mHBs1\[Sigma]],Join[{"mABs1\[Sigma]"},mABs1\[Sigma]],Join[{"mHpmBs1\[Sigma]"},mHpmBs1\[Sigma]],Join[{"kvBs1\[Sigma]"},kvBs1\[Sigma]],Join[{"ktBs1\[Sigma]"},ktBs1\[Sigma]],Join[{"kbBs1\[Sigma]"},kbBs1\[Sigma]],Join[{"tanBBs1\[Sigma]"},tanBBs1\[Sigma]],Join[{"sinBABs1\[Sigma]"},sinBABs1\[Sigma]],Join[{"mugghBs1\[Sigma]"},muggGamGamBs1\[Sigma]],Join[{"muVBFBs1\[Sigma]"},muVBFGamGamBs1\[Sigma]],Join[{"muTotBs1\[Sigma]"},muTotGamGamBs1\[Sigma]],Join[{"sigmagghBs1\[Sigma]"},sigmaBRgghGamGamBs1\[Sigma]],Join[{"sigmaVBF_VFBs1\[Sigma]"},sigmaBRvbfvhGamGamBs1\[Sigma]],Join[{"sigma_tot_Bs1\[Sigma]"},sigmaBRhGamGamBs1\[Sigma]],Join[{"m12Bs1\[Sigma]"},Sqrt[m122Bs1\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_Bs_1sigma_limits.txt"],GnuTabBs1\[Sigma],"Table"]

GnuTabBs2\[Sigma]=Transpose[{Join[{"mh1Bs2\[Sigma]"},mh1Bs2\[Sigma]],Join[{"mHBs2\[Sigma]"},mHBs2\[Sigma]],Join[{"mABs2\[Sigma]"},mABs2\[Sigma]],Join[{"mHpmBs2\[Sigma]"},mHpmBs2\[Sigma]],Join[{"kvBs2\[Sigma]"},kvBs2\[Sigma]],Join[{"ktBs2\[Sigma]"},ktBs2\[Sigma]],Join[{"kbBs2\[Sigma]"},kbBs2\[Sigma]],Join[{"tanBBs2\[Sigma]"},tanBBs2\[Sigma]],Join[{"sinBABs2\[Sigma]"},sinBABs2\[Sigma]],Join[{"mugghBs2\[Sigma]"},muggGamGamBs2\[Sigma]],Join[{"muVBFBs2\[Sigma]"},muVBFGamGamBs2\[Sigma]],Join[{"muTotBs2\[Sigma]"},muTotGamGamBs2\[Sigma]],Join[{"sigmagghBs2\[Sigma]"},sigmaBRgghGamGamBs2\[Sigma]],Join[{"sigmaVBF_VFBs2\[Sigma]"},sigmaBRvbfvhGamGamBs2\[Sigma]],Join[{"sigma_tot_Bs2\[Sigma]"},sigmaBRhGamGamBs2\[Sigma]],Join[{"m12Bs2\[Sigma]"},Sqrt[m122Bs2\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_Bs_2sigma_limits.txt"],GnuTabBs2\[Sigma],"Table"]
Print["Conversion Done"]




(* points passing \[CapitalDelta]Md constraints  (+ Stab, stab. unit.pert., B\[Rule]Xs \[Gamma] but NOT Bs\[Rule]\[Mu]\[Mu] !! ) *)

Tab\[CapitalDelta]M2\[Sigma]=Array[0,lArray,1]; 
Tab\[CapitalDelta]M1\[Sigma]=Array[0,lArray,1];(* array with 1 if the point pass TH constraints, 0 else *)

Clear[size\[CapitalDelta]M1\[Sigma],size\[CapitalDelta]M2\[Sigma]];
size\[CapitalDelta]M1\[Sigma]=0; (* count the number of points passing Stab constraints *)
size\[CapitalDelta]M2\[Sigma]=0;

(* loop which fills the array defined above *)
For[iloop=1,iloop<=lArray,iloop++,If[TabBxs1\[Sigma][[iloop]]==1&&(limDeltaMd-limDeltaMdDown)<=DeltaMd[[iloop]]<=(limDeltaMd+limDeltaMdUp),Tab\[CapitalDelta]M1\[Sigma][[iloop]]=1;size\[CapitalDelta]M1\[Sigma]++;Tab\[CapitalDelta]M2\[Sigma][[iloop]]=1;size\[CapitalDelta]M2\[Sigma]++];If[(TabBxs1\[Sigma][[iloop]]==1||TabBxs2\[Sigma][[iloop]]==1)&&(limDeltaMd-2*limDeltaMdDown)<=DeltaMd[[iloop]]<=(limDeltaMd+2*limDeltaMdUp),Tab\[CapitalDelta]M2\[Sigma][[iloop]]=1;size\[CapitalDelta]M2\[Sigma]++]]

Print["Loop for size \[CapitalDelta]Md DONE"]
(* initialisation of the mass arrays *)
mh1\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
mH\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
mHpm\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
mA\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
m122\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];

mh1\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
mH\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
mHpm\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
mA\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
m122\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];

(* \[Mu] arrays *)
muggGamGam\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1]; 
muVBFGamGam\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
muTotGamGam\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];

muggGamGam\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1]; 
muVBFGamGam\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
muTotGamGam\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];

(* \[Kappa] arrays *)
kv\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1]; 
kt\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1]; 
kb\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1]; 

kv\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1]; 
kt\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1]; 
kb\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];

(* angles *)
tanB\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
sinBA\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];

tanB\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
sinBA\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];

(* \[Sigma]*BR *)
sigmaBRgghGamGam\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
sigmaBRvbfvhGamGam\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];
sigmaBRhGamGam\[CapitalDelta]M1\[Sigma]=Array[0,size\[CapitalDelta]M1\[Sigma],1];

sigmaBRgghGamGam\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
sigmaBRvbfvhGamGam\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
sigmaBRhGamGam\[CapitalDelta]M2\[Sigma]=Array[0,size\[CapitalDelta]M2\[Sigma],1];
(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop,kloop]
jloop=0;
kloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[Tab\[CapitalDelta]M1\[Sigma][[iloop]]==1,jloop++;muggGamGam\[CapitalDelta]M1\[Sigma][[jloop]]=muggGamGam[[iloop]];muVBFGamGam\[CapitalDelta]M1\[Sigma][[jloop]]=muVBFGamGam[[iloop]];kv\[CapitalDelta]M1\[Sigma][[jloop]]=kv2[[iloop]];kt\[CapitalDelta]M1\[Sigma][[jloop]]=kt2[[iloop]];kb\[CapitalDelta]M1\[Sigma][[jloop]]=kb2[[iloop]];mh1\[CapitalDelta]M1\[Sigma][[jloop]]=mh1[[iloop]];mH\[CapitalDelta]M1\[Sigma][[jloop]]=mH[[iloop]];mHpm\[CapitalDelta]M1\[Sigma][[jloop]]=mHpm[[iloop]];mA\[CapitalDelta]M1\[Sigma][[jloop]]=mA[[iloop]];m122\[CapitalDelta]M1\[Sigma][[jloop]]=m122[[iloop]];tanB\[CapitalDelta]M1\[Sigma][[jloop]]=tanB[[iloop]];sinBA\[CapitalDelta]M1\[Sigma][[jloop]]=sinBA[[iloop]];sigmaBRgghGamGam\[CapitalDelta]M1\[Sigma][[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGam\[CapitalDelta]M1\[Sigma][[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGam\[CapitalDelta]M1\[Sigma][[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGam\[CapitalDelta]M1\[Sigma][[jloop]]=muTotGamGam[[iloop]]];
If[Tab\[CapitalDelta]M2\[Sigma][[iloop]]==1,kloop++;muggGamGam\[CapitalDelta]M2\[Sigma][[kloop]]=muggGamGam[[iloop]];muVBFGamGam\[CapitalDelta]M2\[Sigma][[kloop]]=muVBFGamGam[[iloop]];kv\[CapitalDelta]M2\[Sigma][[kloop]]=kv2[[iloop]];kt\[CapitalDelta]M2\[Sigma][[kloop]]=kt2[[iloop]];kb\[CapitalDelta]M2\[Sigma][[kloop]]=kb2[[iloop]];mh1\[CapitalDelta]M2\[Sigma][[kloop]]=mh1[[iloop]];mH\[CapitalDelta]M2\[Sigma][[kloop]]=mH[[iloop]];mHpm\[CapitalDelta]M2\[Sigma][[kloop]]=mHpm[[iloop]];mA\[CapitalDelta]M2\[Sigma][[kloop]]=mA[[iloop]];m122\[CapitalDelta]M2\[Sigma][[kloop]]=m122[[iloop]];tanB\[CapitalDelta]M2\[Sigma][[kloop]]=tanB[[iloop]];sinBA\[CapitalDelta]M2\[Sigma][[kloop]]=sinBA[[iloop]];sigmaBRgghGamGam\[CapitalDelta]M2\[Sigma][[kloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGam\[CapitalDelta]M2\[Sigma][[kloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGam\[CapitalDelta]M2\[Sigma][[kloop]]=sigmaBRhGamGam[[iloop]];muTotGamGam\[CapitalDelta]M2\[Sigma][[kloop]]=muTotGamGam[[iloop]]]]

Print["Loop for \[CapitalDelta]M DONE"]





size\[CapitalDelta]M2\[Sigma]


(* convert to GnuPlot files *)
GnuTab\[CapitalDelta]M1\[Sigma]=Transpose[{Join[{"mh1\[CapitalDelta]M1\[Sigma]"},mh1\[CapitalDelta]M1\[Sigma]],Join[{"mH\[CapitalDelta]M1\[Sigma]"},mH\[CapitalDelta]M1\[Sigma]],Join[{"mA\[CapitalDelta]M1\[Sigma]"},mA\[CapitalDelta]M1\[Sigma]],Join[{"mHpm\[CapitalDelta]M1\[Sigma]"},mHpm\[CapitalDelta]M1\[Sigma]],Join[{"kv\[CapitalDelta]M1\[Sigma]"},kv\[CapitalDelta]M1\[Sigma]],Join[{"kt\[CapitalDelta]M1\[Sigma]"},kt\[CapitalDelta]M1\[Sigma]],Join[{"kb\[CapitalDelta]M1\[Sigma]"},kb\[CapitalDelta]M1\[Sigma]],Join[{"tanB\[CapitalDelta]M1\[Sigma]"},tanB\[CapitalDelta]M1\[Sigma]],Join[{"sinBA\[CapitalDelta]M1\[Sigma]"},sinBA\[CapitalDelta]M1\[Sigma]],Join[{"muggh\[CapitalDelta]M1\[Sigma]"},muggGamGam\[CapitalDelta]M1\[Sigma]],Join[{"muVBF\[CapitalDelta]M1\[Sigma]"},muVBFGamGam\[CapitalDelta]M1\[Sigma]],Join[{"muTot\[CapitalDelta]M1\[Sigma]"},muTotGamGam\[CapitalDelta]M1\[Sigma]],Join[{"sigmaggh\[CapitalDelta]M1\[Sigma]"},sigmaBRgghGamGam\[CapitalDelta]M1\[Sigma]],Join[{"sigmaVBF_VF\[CapitalDelta]M1\[Sigma]"},sigmaBRvbfvhGamGam\[CapitalDelta]M1\[Sigma]],Join[{"sigma_tot_\[CapitalDelta]M1\[Sigma]"},sigmaBRhGamGam\[CapitalDelta]M1\[Sigma]],Join[{"m12\[CapitalDelta]M1\[Sigma]"},Sqrt[m122\[CapitalDelta]M1\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_DeltaM_1sigma_limits.txt"],GnuTab\[CapitalDelta]M1\[Sigma],"Table"]

GnuTab\[CapitalDelta]M2\[Sigma]=Transpose[{Join[{"mh1\[CapitalDelta]M2\[Sigma]"},mh1\[CapitalDelta]M2\[Sigma]],Join[{"mH\[CapitalDelta]M2\[Sigma]"},mH\[CapitalDelta]M2\[Sigma]],Join[{"mA\[CapitalDelta]M2\[Sigma]"},mA\[CapitalDelta]M2\[Sigma]],Join[{"mHpm\[CapitalDelta]M2\[Sigma]"},mHpm\[CapitalDelta]M2\[Sigma]],Join[{"kv\[CapitalDelta]M2\[Sigma]"},kv\[CapitalDelta]M2\[Sigma]],Join[{"kt\[CapitalDelta]M2\[Sigma]"},kt\[CapitalDelta]M2\[Sigma]],Join[{"kb\[CapitalDelta]M2\[Sigma]"},kb\[CapitalDelta]M2\[Sigma]],Join[{"tanB\[CapitalDelta]M2\[Sigma]"},tanB\[CapitalDelta]M2\[Sigma]],Join[{"sinBA\[CapitalDelta]M2\[Sigma]"},sinBA\[CapitalDelta]M2\[Sigma]],Join[{"muggh\[CapitalDelta]M2\[Sigma]"},muggGamGam\[CapitalDelta]M2\[Sigma]],Join[{"muVBF\[CapitalDelta]M2\[Sigma]"},muVBFGamGam\[CapitalDelta]M2\[Sigma]],Join[{"muTot\[CapitalDelta]M2\[Sigma]"},muTotGamGam\[CapitalDelta]M2\[Sigma]],Join[{"sigmaggh\[CapitalDelta]M2\[Sigma]"},sigmaBRgghGamGam\[CapitalDelta]M2\[Sigma]],Join[{"sigmaVBF_VF\[CapitalDelta]M2\[Sigma]"},sigmaBRvbfvhGamGam\[CapitalDelta]M2\[Sigma]],Join[{"sigma_tot_\[CapitalDelta]M2\[Sigma]"},sigmaBRhGamGam\[CapitalDelta]M2\[Sigma]],Join[{"m12\[CapitalDelta]M2\[Sigma]"},Sqrt[m122\[CapitalDelta]M2\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_DeltaM_2sigma_limits.txt"],GnuTab\[CapitalDelta]M2\[Sigma],"Table"]
Print["Conversion Done"]






(* points passing \[CapitalDelta]0 constraints  (+ Stab, stab. unit.pert., B\[Rule]Xs \[Gamma] but NOT Bs\[Rule]\[Mu]\[Mu] nor \[CapitalDelta]Md !! ) *)

Tab\[CapitalDelta]02\[Sigma]=Array[0,lArray,1]; 
Tab\[CapitalDelta]01\[Sigma]=Array[0,lArray,1];(* array with 1 if the point pass TH constraints, 0 else *)

Clear[size\[CapitalDelta]01\[Sigma],size\[CapitalDelta]02\[Sigma]];
size\[CapitalDelta]01\[Sigma]=0; (* count the number of points passing Stab constraints *)
size\[CapitalDelta]02\[Sigma]=0;
(* loop which fills the array defined above *)

For[iloop=1,iloop<=lArray,iloop++,If[TabBxs1\[Sigma][[iloop]]==1&&(limDelta0-limDelta0Down)<=Delta0[[iloop]]<=(limDelta0+limDelta0Up),Tab\[CapitalDelta]01\[Sigma][[iloop]]=1;size\[CapitalDelta]01\[Sigma]++;Tab\[CapitalDelta]02\[Sigma][[iloop]]=1;size\[CapitalDelta]02\[Sigma]++];If[(TabBxs1\[Sigma][[iloop]]==1||TabBxs2\[Sigma][[iloop]]==1)&&(limDelta0-2*limDelta0Down)<=Delta0[[iloop]]<=(limDelta0+2*limDelta0Up),Tab\[CapitalDelta]02\[Sigma][[iloop]]=1;size\[CapitalDelta]02\[Sigma]++]]
Print["Loop for size \[CapitalDelta]0 DONE"]

(* initialisation of the mass arrays *)
mh1\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
mH\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
mHpm\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
mA\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
m122\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];

mh1\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
mH\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
mHpm\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
mA\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
m122\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];

(* \[Mu] arrays *)
muggGamGam\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1]; 
muVBFGamGam\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
muTotGamGam\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];

muggGamGam\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1]; 
muVBFGamGam\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
muTotGamGam\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];

(* \[Kappa] arrays *)
kv\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1]; 
kt\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1]; 
kb\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1]; 

kv\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1]; 
kt\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1]; 
kb\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];

(* angles *)
tanB\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
sinBA\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];

tanB\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
sinBA\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];

(* \[Sigma]*BR *)
sigmaBRgghGamGam\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
sigmaBRvbfvhGamGam\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];
sigmaBRhGamGam\[CapitalDelta]01\[Sigma]=Array[0,size\[CapitalDelta]01\[Sigma],1];

sigmaBRgghGamGam\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
sigmaBRvbfvhGamGam\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
sigmaBRhGamGam\[CapitalDelta]02\[Sigma]=Array[0,size\[CapitalDelta]02\[Sigma],1];
(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop,kloop]
jloop=0;
kloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[Tab\[CapitalDelta]01\[Sigma][[iloop]]==1,jloop++;muggGamGam\[CapitalDelta]01\[Sigma][[jloop]]=muggGamGam[[iloop]];muVBFGamGam\[CapitalDelta]01\[Sigma][[jloop]]=muVBFGamGam[[iloop]];kv\[CapitalDelta]01\[Sigma][[jloop]]=kv2[[iloop]];kt\[CapitalDelta]01\[Sigma][[jloop]]=kt2[[iloop]];kb\[CapitalDelta]01\[Sigma][[jloop]]=kb2[[iloop]];mh1\[CapitalDelta]01\[Sigma][[jloop]]=mh1[[iloop]];mH\[CapitalDelta]01\[Sigma][[jloop]]=mH[[iloop]];mHpm\[CapitalDelta]01\[Sigma][[jloop]]=mHpm[[iloop]];mA\[CapitalDelta]01\[Sigma][[jloop]]=mA[[iloop]];m122\[CapitalDelta]01\[Sigma][[jloop]]=m122[[iloop]];tanB\[CapitalDelta]01\[Sigma][[jloop]]=tanB[[iloop]];sinBA\[CapitalDelta]01\[Sigma][[jloop]]=sinBA[[iloop]];sigmaBRgghGamGam\[CapitalDelta]01\[Sigma][[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGam\[CapitalDelta]01\[Sigma][[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGam\[CapitalDelta]01\[Sigma][[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGam\[CapitalDelta]01\[Sigma][[jloop]]=muTotGamGam[[iloop]]];
If[Tab\[CapitalDelta]02\[Sigma][[iloop]]==1,kloop++;muggGamGam\[CapitalDelta]02\[Sigma][[kloop]]=muggGamGam[[iloop]];muVBFGamGam\[CapitalDelta]02\[Sigma][[kloop]]=muVBFGamGam[[iloop]];kv\[CapitalDelta]02\[Sigma][[kloop]]=kv2[[iloop]];kt\[CapitalDelta]02\[Sigma][[kloop]]=kt2[[iloop]];kb\[CapitalDelta]02\[Sigma][[kloop]]=kb2[[iloop]];mh1\[CapitalDelta]02\[Sigma][[kloop]]=mh1[[iloop]];mH\[CapitalDelta]02\[Sigma][[kloop]]=mH[[iloop]];mHpm\[CapitalDelta]02\[Sigma][[kloop]]=mHpm[[iloop]];mA\[CapitalDelta]02\[Sigma][[kloop]]=mA[[iloop]];m122\[CapitalDelta]02\[Sigma][[kloop]]=m122[[iloop]];tanB\[CapitalDelta]02\[Sigma][[kloop]]=tanB[[iloop]];sinBA\[CapitalDelta]02\[Sigma][[kloop]]=sinBA[[iloop]];sigmaBRgghGamGam\[CapitalDelta]02\[Sigma][[kloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGam\[CapitalDelta]02\[Sigma][[kloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGam\[CapitalDelta]02\[Sigma][[kloop]]=sigmaBRhGamGam[[iloop]];muTotGamGam\[CapitalDelta]02\[Sigma][[kloop]]=muTotGamGam[[iloop]]]]
Print["Loop for \[CapitalDelta]0 DONE"]




size\[CapitalDelta]02\[Sigma]


(* convert to GnuPlot files *)

GnuTab\[CapitalDelta]01\[Sigma]=Transpose[{Join[{"mh1\[CapitalDelta]01\[Sigma]"},mh1\[CapitalDelta]01\[Sigma]],Join[{"mH\[CapitalDelta]01\[Sigma]"},mH\[CapitalDelta]01\[Sigma]],Join[{"mA\[CapitalDelta]01\[Sigma]"},mA\[CapitalDelta]01\[Sigma]],Join[{"mHpm\[CapitalDelta]01\[Sigma]"},mHpm\[CapitalDelta]01\[Sigma]],Join[{"kv\[CapitalDelta]01\[Sigma]"},kv\[CapitalDelta]01\[Sigma]],Join[{"kt\[CapitalDelta]01\[Sigma]"},kt\[CapitalDelta]01\[Sigma]],Join[{"kb\[CapitalDelta]01\[Sigma]"},kb\[CapitalDelta]01\[Sigma]],Join[{"tanB\[CapitalDelta]01\[Sigma]"},tanB\[CapitalDelta]01\[Sigma]],Join[{"sinBA\[CapitalDelta]01\[Sigma]"},sinBA\[CapitalDelta]01\[Sigma]],Join[{"muggh\[CapitalDelta]01\[Sigma]"},muggGamGam\[CapitalDelta]01\[Sigma]],Join[{"muVBF\[CapitalDelta]01\[Sigma]"},muVBFGamGam\[CapitalDelta]01\[Sigma]],Join[{"muTot\[CapitalDelta]01\[Sigma]"},muTotGamGam\[CapitalDelta]01\[Sigma]],Join[{"sigmaggh\[CapitalDelta]01\[Sigma]"},sigmaBRgghGamGam\[CapitalDelta]01\[Sigma]],Join[{"sigmaVBF_VF\[CapitalDelta]01\[Sigma]"},sigmaBRvbfvhGamGam\[CapitalDelta]01\[Sigma]],Join[{"sigma_tot_\[CapitalDelta]01\[Sigma]"},sigmaBRhGamGam\[CapitalDelta]01\[Sigma]],Join[{"m12\[CapitalDelta]01\[Sigma]"},Sqrt[m122\[CapitalDelta]01\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_Delta0_1sigma_limits.txt"],GnuTab\[CapitalDelta]01\[Sigma],"Table"]

GnuTab\[CapitalDelta]02\[Sigma]=Transpose[{Join[{"mh1\[CapitalDelta]02\[Sigma]"},mh1\[CapitalDelta]02\[Sigma]],Join[{"mH\[CapitalDelta]02\[Sigma]"},mH\[CapitalDelta]02\[Sigma]],Join[{"mA\[CapitalDelta]02\[Sigma]"},mA\[CapitalDelta]02\[Sigma]],Join[{"mHpm\[CapitalDelta]02\[Sigma]"},mHpm\[CapitalDelta]02\[Sigma]],Join[{"kv\[CapitalDelta]02\[Sigma]"},kv\[CapitalDelta]02\[Sigma]],Join[{"kt\[CapitalDelta]02\[Sigma]"},kt\[CapitalDelta]02\[Sigma]],Join[{"kb\[CapitalDelta]02\[Sigma]"},kb\[CapitalDelta]02\[Sigma]],Join[{"tanB\[CapitalDelta]02\[Sigma]"},tanB\[CapitalDelta]02\[Sigma]],Join[{"sinBA\[CapitalDelta]02\[Sigma]"},sinBA\[CapitalDelta]02\[Sigma]],Join[{"muggh\[CapitalDelta]02\[Sigma]"},muggGamGam\[CapitalDelta]02\[Sigma]],Join[{"muVBF\[CapitalDelta]02\[Sigma]"},muVBFGamGam\[CapitalDelta]02\[Sigma]],Join[{"muTot\[CapitalDelta]02\[Sigma]"},muTotGamGam\[CapitalDelta]02\[Sigma]],Join[{"sigmaggh\[CapitalDelta]02\[Sigma]"},sigmaBRgghGamGam\[CapitalDelta]02\[Sigma]],Join[{"sigmaVBF_VF\[CapitalDelta]02\[Sigma]"},sigmaBRvbfvhGamGam\[CapitalDelta]02\[Sigma]],Join[{"sigma_tot_\[CapitalDelta]02\[Sigma]"},sigmaBRhGamGam\[CapitalDelta]02\[Sigma]],Join[{"m12\[CapitalDelta]02\[Sigma]"},Sqrt[m122\[CapitalDelta]02\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_Delta0_2sigma_limits.txt"],GnuTab\[CapitalDelta]02\[Sigma],"Table"]
Print["Conversion Done"]


(* points passing S, T OR U constraints *)

TabS1\[Sigma]=Array[0,lArray,1]; 
TabS2\[Sigma]=Array[0,lArray,1]; 
TabT1\[Sigma]=Array[0,lArray,1];
TabT2\[Sigma]=Array[0,lArray,1];
TabU1\[Sigma]=Array[0,lArray,1];
TabU2\[Sigma]=Array[0,lArray,1];

Clear[sizeS1\[Sigma],sizeS2\[Sigma],sizeT1\[Sigma],sizeT2\[Sigma],sizeU1\[Sigma],sizeU2\[Sigma]];

sizeS1\[Sigma]=0; (* count the number of points passing Stab constraints *)
sizeT1\[Sigma]=0;
sizeU1\[Sigma]=0;

sizeS2\[Sigma]=0; (* count the number of points passing Stab constraints *)
sizeT2\[Sigma]=0;
sizeU2\[Sigma]=0;

(* loop which fills the array defined above *)

For[iloop=1,iloop<=lArray,iloop++,If[(limS-limSDown)<=S[[iloop]]<=(limS+limSUp),TabS1\[Sigma][[iloop]]=1;sizeS1\[Sigma]++;TabS2\[Sigma][[iloop]]=1;sizeS2\[Sigma]++,If[(limS-2*limSDown)<=S[[iloop]]<=(limS+2*limSUp),TabS2\[Sigma][[iloop]]=1;sizeS2\[Sigma]++]];If[(limT-limTDown)<=T[[iloop]]<=(limT+limTUp),TabT1\[Sigma][[iloop]]=1;sizeT1\[Sigma]++;TabT2\[Sigma][[iloop]]=1;sizeT2\[Sigma]++,If[(limT-2*limTDown)<=T[[iloop]]<=(limT+2*limTUp),TabT2\[Sigma][[iloop]]=1;sizeT2\[Sigma]++]];If[(limU-limUDown)<=U[[iloop]]<=(limU+limUUp),TabU1\[Sigma][[iloop]]=1;sizeU1\[Sigma]++;TabU2\[Sigma][[iloop]]=1;sizeU2\[Sigma]++,If[(limU-2*limUDown)<=U[[iloop]]<=(limU+2*limUUp),TabU2\[Sigma][[iloop]]=1;sizeU2\[Sigma]++]]]

Print["Loop for sizeStab DONE"]

(* initialisation of the mass arrays *)
mh1S1\[Sigma]=Array[0,sizeS1\[Sigma],1];
mHS1\[Sigma]=Array[0,sizeS1\[Sigma],1];
mHpmS1\[Sigma]=Array[0,sizeS1\[Sigma],1];
mAS1\[Sigma]=Array[0,sizeS1\[Sigma],1];
m122S1\[Sigma]=Array[0,sizeS1\[Sigma],1];

mh1S2\[Sigma]=Array[0,sizeS2\[Sigma],1];
mHS2\[Sigma]=Array[0,sizeS2\[Sigma],1];
mHpmS2\[Sigma]=Array[0,sizeS2\[Sigma],1];
mAS2\[Sigma]=Array[0,sizeS2\[Sigma],1];
m122S2\[Sigma]=Array[0,sizeS2\[Sigma],1];

mh1T1\[Sigma]=Array[0,sizeT1\[Sigma],1];
mHT1\[Sigma]=Array[0,sizeT1\[Sigma],1];
mHpmT1\[Sigma]=Array[0,sizeT1\[Sigma],1];
mAT1\[Sigma]=Array[0,sizeT1\[Sigma],1];
m122T1\[Sigma]=Array[0,sizeT1\[Sigma],1];

mh1T2\[Sigma]=Array[0,sizeT2\[Sigma],1];
mHT2\[Sigma]=Array[0,sizeT2\[Sigma],1];
mHpmT2\[Sigma]=Array[0,sizeT2\[Sigma],1];
mAT2\[Sigma]=Array[0,sizeT2\[Sigma],1];
m122T2\[Sigma]=Array[0,sizeT2\[Sigma],1];

mh1U1\[Sigma]=Array[0,sizeU1\[Sigma],1];
mHU1\[Sigma]=Array[0,sizeU1\[Sigma],1];
mHpmU1\[Sigma]=Array[0,sizeU1\[Sigma],1];
mAU1\[Sigma]=Array[0,sizeU1\[Sigma],1];
m122U1\[Sigma]=Array[0,sizeU1\[Sigma],1];

mh1U2\[Sigma]=Array[0,sizeU2\[Sigma],1];
mHU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
mHpmU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
mAU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
m122U2\[Sigma]=Array[0,sizeU2\[Sigma],1];

(* \[Mu] arrays *)
muggGamGamS1\[Sigma]=Array[0,sizeS1\[Sigma],1]; 
muVBFGamGamS1\[Sigma]=Array[0,sizeS1\[Sigma],1];
muTotGamGamS1\[Sigma]=Array[0,sizeS1\[Sigma],1];

muggGamGamS2\[Sigma]=Array[0,sizeS2\[Sigma],1]; 
muVBFGamGamS2\[Sigma]=Array[0,sizeS2\[Sigma],1];
muTotGamGamS2\[Sigma]=Array[0,sizeS2\[Sigma],1];

muggGamGamT1\[Sigma]=Array[0,sizeT1\[Sigma],1]; 
muVBFGamGamT1\[Sigma]=Array[0,sizeT1\[Sigma],1];
muTotGamGamT1\[Sigma]=Array[0,sizeT1\[Sigma],1];

muggGamGamT2\[Sigma]=Array[0,sizeT2\[Sigma],1]; 
muVBFGamGamT2\[Sigma]=Array[0,sizeT2\[Sigma],1];
muTotGamGamT2\[Sigma]=Array[0,sizeT2\[Sigma],1];

muggGamGamU1\[Sigma]=Array[0,sizeU1\[Sigma],1]; 
muVBFGamGamU1\[Sigma]=Array[0,sizeU1\[Sigma],1];
muTotGamGamU1\[Sigma]=Array[0,sizeU1\[Sigma],1];

muggGamGamU2\[Sigma]=Array[0,sizeU2\[Sigma],1]; 
muVBFGamGamU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
muTotGamGamU2\[Sigma]=Array[0,sizeU2\[Sigma],1];

(* \[Kappa] arrays *)
kvS1\[Sigma]=Array[0,sizeS1\[Sigma],1]; 
ktS1\[Sigma]=Array[0,sizeS1\[Sigma],1]; 
kbS1\[Sigma]=Array[0,sizeS1\[Sigma],1]; 

kvS2\[Sigma]=Array[0,sizeS2\[Sigma],1]; 
ktS2\[Sigma]=Array[0,sizeS2\[Sigma],1]; 
kbS2\[Sigma]=Array[0,sizeS2\[Sigma],1];

kvT1\[Sigma]=Array[0,sizeT1\[Sigma],1]; 
ktT1\[Sigma]=Array[0,sizeT1\[Sigma],1]; 
kbT1\[Sigma]=Array[0,sizeT1\[Sigma],1]; 

kvT2\[Sigma]=Array[0,sizeT2\[Sigma],1]; 
ktT2\[Sigma]=Array[0,sizeT2\[Sigma],1]; 
kbT2\[Sigma]=Array[0,sizeT2\[Sigma],1];

kvU1\[Sigma]=Array[0,sizeU1\[Sigma],1]; 
ktU1\[Sigma]=Array[0,sizeU1\[Sigma],1]; 
kbU1\[Sigma]=Array[0,sizeU1\[Sigma],1]; 

kvU2\[Sigma]=Array[0,sizeU2\[Sigma],1]; 
ktU2\[Sigma]=Array[0,sizeU2\[Sigma],1]; 
kbU2\[Sigma]=Array[0,sizeU2\[Sigma],1];

(* angles *)
tanBS1\[Sigma]=Array[0,sizeS1\[Sigma],1];
sinBAS1\[Sigma]=Array[0,sizeS1\[Sigma],1];

tanBS2\[Sigma]=Array[0,sizeS2\[Sigma],1];
sinBAS2\[Sigma]=Array[0,sizeS2\[Sigma],1];

tanBT1\[Sigma]=Array[0,sizeT1\[Sigma],1];
sinBAT1\[Sigma]=Array[0,sizeT1\[Sigma],1];

tanBT2\[Sigma]=Array[0,sizeT2\[Sigma],1];
sinBAT2\[Sigma]=Array[0,sizeT2\[Sigma],1];

tanBU1\[Sigma]=Array[0,sizeU1\[Sigma],1];
sinBAU1\[Sigma]=Array[0,sizeU1\[Sigma],1];

tanBU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
sinBAU2\[Sigma]=Array[0,sizeU2\[Sigma],1];

(* \[Sigma]*BR *)
sigmaBRgghGamGamS1\[Sigma]=Array[0,sizeS1\[Sigma],1];
sigmaBRvbfvhGamGamS1\[Sigma]=Array[0,sizeS1\[Sigma],1];
sigmaBRhGamGamS1\[Sigma]=Array[0,sizeS1\[Sigma],1];

sigmaBRgghGamGamS2\[Sigma]=Array[0,sizeS2\[Sigma],1];
sigmaBRvbfvhGamGamS2\[Sigma]=Array[0,sizeS2\[Sigma],1];
sigmaBRhGamGamS2\[Sigma]=Array[0,sizeS2\[Sigma],1];

sigmaBRgghGamGamT1\[Sigma]=Array[0,sizeT1\[Sigma],1];
sigmaBRvbfvhGamGamT1\[Sigma]=Array[0,sizeT1\[Sigma],1];
sigmaBRhGamGamT1\[Sigma]=Array[0,sizeT1\[Sigma],1];

sigmaBRgghGamGamT2\[Sigma]=Array[0,sizeT2\[Sigma],1];
sigmaBRvbfvhGamGamT2\[Sigma]=Array[0,sizeT2\[Sigma],1];
sigmaBRhGamGamT2\[Sigma]=Array[0,sizeT2\[Sigma],1];

sigmaBRgghGamGamU1\[Sigma]=Array[0,sizeU1\[Sigma],1];
sigmaBRvbfvhGamGamU1\[Sigma]=Array[0,sizeU1\[Sigma],1];
sigmaBRhGamGamU1\[Sigma]=Array[0,sizeU1\[Sigma],1];

sigmaBRgghGamGamU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
sigmaBRvbfvhGamGamU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
sigmaBRhGamGamU2\[Sigma]=Array[0,sizeU2\[Sigma],1];
(* Loop which fills the arrays defined just before with the points passing TH constraints*)
Clear[jloop,kloop]
jloop=0;
kloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabS1\[Sigma][[iloop]]==1,jloop++;muggGamGamS1\[Sigma][[jloop]]=muggGamGam[[iloop]];muVBFGamGamS1\[Sigma][[jloop]]=muVBFGamGam[[iloop]];kvS1\[Sigma][[jloop]]=kv2[[iloop]];ktS1\[Sigma][[jloop]]=kt2[[iloop]];kbS1\[Sigma][[jloop]]=kb2[[iloop]];mh1S1\[Sigma][[jloop]]=mh1[[iloop]];mHS1\[Sigma][[jloop]]=mH[[iloop]];mHpmS1\[Sigma][[jloop]]=mHpm[[iloop]];mAS1\[Sigma][[jloop]]=mA[[iloop]];m122S1\[Sigma][[jloop]]=m122[[iloop]];tanBS1\[Sigma][[jloop]]=tanB[[iloop]];sinBAS1\[Sigma][[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamS1\[Sigma][[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamS1\[Sigma][[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamS1\[Sigma][[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamS1\[Sigma][[jloop]]=muTotGamGam[[iloop]]];
If[TabS2\[Sigma][[iloop]]==1,kloop++;muggGamGamS2\[Sigma][[kloop]]=muggGamGam[[iloop]];muVBFGamGamS2\[Sigma][[kloop]]=muVBFGamGam[[iloop]];kvS2\[Sigma][[kloop]]=kv2[[iloop]];ktS2\[Sigma][[kloop]]=kt2[[iloop]];kbS2\[Sigma][[kloop]]=kb2[[iloop]];mh1S2\[Sigma][[kloop]]=mh1[[iloop]];mHS2\[Sigma][[kloop]]=mH[[iloop]];mHpmS2\[Sigma][[kloop]]=mHpm[[iloop]];mAS2\[Sigma][[kloop]]=mA[[iloop]];m122S2\[Sigma][[kloop]]=m122[[iloop]];tanBS2\[Sigma][[kloop]]=tanB[[iloop]];sinBAS2\[Sigma][[kloop]]=sinBA[[iloop]];sigmaBRgghGamGamS2\[Sigma][[kloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamS2\[Sigma][[kloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamS2\[Sigma][[kloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamS2\[Sigma][[kloop]]=muTotGamGam[[iloop]]]]
Clear[jloop,kloop]
jloop=0;
kloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabT1\[Sigma][[iloop]]==1,jloop++;muggGamGamT1\[Sigma][[jloop]]=muggGamGam[[iloop]];muVBFGamGamT1\[Sigma][[jloop]]=muVBFGamGam[[iloop]];kvT1\[Sigma][[jloop]]=kv2[[iloop]];ktT1\[Sigma][[jloop]]=kt2[[iloop]];kbT1\[Sigma][[jloop]]=kb2[[iloop]];mh1T1\[Sigma][[jloop]]=mh1[[iloop]];mHT1\[Sigma][[jloop]]=mH[[iloop]];mHpmT1\[Sigma][[jloop]]=mHpm[[iloop]];mAT1\[Sigma][[jloop]]=mA[[iloop]];m122T1\[Sigma][[jloop]]=m122[[iloop]];tanBT1\[Sigma][[jloop]]=tanB[[iloop]];sinBAT1\[Sigma][[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamT1\[Sigma][[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamT1\[Sigma][[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamT1\[Sigma][[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamT1\[Sigma][[jloop]]=muTotGamGam[[iloop]]];
If[TabT2\[Sigma][[iloop]]==1,kloop++;muggGamGamT2\[Sigma][[kloop]]=muggGamGam[[iloop]];muVBFGamGamT2\[Sigma][[kloop]]=muVBFGamGam[[iloop]];kvT2\[Sigma][[kloop]]=kv2[[iloop]];ktT2\[Sigma][[kloop]]=kt2[[iloop]];kbT2\[Sigma][[kloop]]=kb2[[iloop]];mh1T2\[Sigma][[kloop]]=mh1[[iloop]];mHT2\[Sigma][[kloop]]=mH[[iloop]];mHpmT2\[Sigma][[kloop]]=mHpm[[iloop]];mAT2\[Sigma][[kloop]]=mA[[iloop]];m122T2\[Sigma][[kloop]]=m122[[iloop]];tanBT2\[Sigma][[kloop]]=tanB[[iloop]];sinBAT2\[Sigma][[kloop]]=sinBA[[iloop]];sigmaBRgghGamGamT2\[Sigma][[kloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamT2\[Sigma][[kloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamT2\[Sigma][[kloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamT2\[Sigma][[kloop]]=muTotGamGam[[iloop]]]]
Clear[jloop,kloop]
jloop=0;
kloop=0;

For[iloop=1,iloop<=lArray,iloop++,If[TabU1\[Sigma][[iloop]]==1,jloop++;muggGamGamU1\[Sigma][[jloop]]=muggGamGam[[iloop]];muVBFGamGamU1\[Sigma][[jloop]]=muVBFGamGam[[iloop]];kvU1\[Sigma][[jloop]]=kv2[[iloop]];ktU1\[Sigma][[jloop]]=kt2[[iloop]];kbU1\[Sigma][[jloop]]=kb2[[iloop]];mh1U1\[Sigma][[jloop]]=mh1[[iloop]];mHU1\[Sigma][[jloop]]=mH[[iloop]];mHpmU1\[Sigma][[jloop]]=mHpm[[iloop]];mAU1\[Sigma][[jloop]]=mA[[iloop]];m122U1\[Sigma][[jloop]]=m122[[iloop]];tanBU1\[Sigma][[jloop]]=tanB[[iloop]];sinBAU1\[Sigma][[jloop]]=sinBA[[iloop]];sigmaBRgghGamGamU1\[Sigma][[jloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamU1\[Sigma][[jloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamU1\[Sigma][[jloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamU1\[Sigma][[jloop]]=muTotGamGam[[iloop]]];
If[TabU2\[Sigma][[iloop]]==1,kloop++;muggGamGamU2\[Sigma][[kloop]]=muggGamGam[[iloop]];muVBFGamGamU2\[Sigma][[kloop]]=muVBFGamGam[[iloop]];kvU2\[Sigma][[kloop]]=kv2[[iloop]];ktU2\[Sigma][[kloop]]=kt2[[iloop]];kbU2\[Sigma][[kloop]]=kb2[[iloop]];mh1U2\[Sigma][[kloop]]=mh1[[iloop]];mHU2\[Sigma][[kloop]]=mH[[iloop]];mHpmU2\[Sigma][[kloop]]=mHpm[[iloop]];mAU2\[Sigma][[kloop]]=mA[[iloop]];m122U2\[Sigma][[kloop]]=m122[[iloop]];tanBU2\[Sigma][[kloop]]=tanB[[iloop]];sinBAU2\[Sigma][[kloop]]=sinBA[[iloop]];sigmaBRgghGamGamU2\[Sigma][[kloop]]=sigmaBRgghGamGam[[iloop]];sigmaBRvbfvhGamGamU2\[Sigma][[kloop]]=sigmaBRvbfvhGamGam[[iloop]];sigmaBRhGamGamU2\[Sigma][[kloop]]=sigmaBRhGamGam[[iloop]];muTotGamGamU2\[Sigma][[kloop]]=muTotGamGam[[iloop]]]]
Print["Loop for Stab DONE"]




(* convert to GnuPlot files *)

GnuTabS1\[Sigma]=Transpose[{Join[{"mh1S1\[Sigma]"},mh1S1\[Sigma]],Join[{"mHS1\[Sigma]"},mHS1\[Sigma]],Join[{"mAS1\[Sigma]"},mAS1\[Sigma]],Join[{"mHpmS1\[Sigma]"},mHpmS1\[Sigma]],Join[{"kvS1\[Sigma]"},kvS1\[Sigma]],Join[{"ktS1\[Sigma]"},ktS1\[Sigma]],Join[{"kbS1\[Sigma]"},kbS1\[Sigma]],Join[{"tanBS1\[Sigma]"},tanBS1\[Sigma]],Join[{"sinBAS1\[Sigma]"},sinBAS1\[Sigma]],Join[{"mugghS1\[Sigma]"},muggGamGamS1\[Sigma]],Join[{"muVBFS1\[Sigma]"},muVBFGamGamS1\[Sigma]],Join[{"muTotS1\[Sigma]"},muTotGamGamS1\[Sigma]],Join[{"sigmagghS1\[Sigma]"},sigmaBRgghGamGamS1\[Sigma]],Join[{"sigmaVBF_VFS1\[Sigma]"},sigmaBRvbfvhGamGamS1\[Sigma]],Join[{"sigma_tot_S1\[Sigma]"},sigmaBRhGamGamS1\[Sigma]],Join[{"m12S1\[Sigma]"},Sqrt[m122S1\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_S_1sigma_limits.txt"],GnuTabS1\[Sigma],"Table"]

GnuTabS2\[Sigma]=Transpose[{Join[{"mh1S2\[Sigma]"},mh1S2\[Sigma]],Join[{"mHS2\[Sigma]"},mHS2\[Sigma]],Join[{"mAS2\[Sigma]"},mAS2\[Sigma]],Join[{"mHpmS2\[Sigma]"},mHpmS2\[Sigma]],Join[{"kvS2\[Sigma]"},kvS2\[Sigma]],Join[{"ktS2\[Sigma]"},ktS2\[Sigma]],Join[{"kbS2\[Sigma]"},kbS2\[Sigma]],Join[{"tanBS2\[Sigma]"},tanBS2\[Sigma]],Join[{"sinBAS2\[Sigma]"},sinBAS2\[Sigma]],Join[{"mugghS2\[Sigma]"},muggGamGamS2\[Sigma]],Join[{"muVBFS2\[Sigma]"},muVBFGamGamS2\[Sigma]],Join[{"muTotS2\[Sigma]"},muTotGamGamS2\[Sigma]],Join[{"sigmagghS2\[Sigma]"},sigmaBRgghGamGamS2\[Sigma]],Join[{"sigmaVBF_VFS2\[Sigma]"},sigmaBRvbfvhGamGamS2\[Sigma]],Join[{"sigma_tot_S2\[Sigma]"},sigmaBRhGamGamS2\[Sigma]],Join[{"m12S2\[Sigma]"},Sqrt[m122S2\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_S_2sigma_limits.txt"],GnuTabS2\[Sigma],"Table"]


GnuTabT1\[Sigma]=Transpose[{Join[{"mh1T1\[Sigma]"},mh1T1\[Sigma]],Join[{"mHT1\[Sigma]"},mHT1\[Sigma]],Join[{"mAT1\[Sigma]"},mAT1\[Sigma]],Join[{"mHpmT1\[Sigma]"},mHpmT1\[Sigma]],Join[{"kvT1\[Sigma]"},kvT1\[Sigma]],Join[{"ktT1\[Sigma]"},ktT1\[Sigma]],Join[{"kbT1\[Sigma]"},kbT1\[Sigma]],Join[{"tanBT1\[Sigma]"},tanBT1\[Sigma]],Join[{"sinBAT1\[Sigma]"},sinBAT1\[Sigma]],Join[{"mugghT1\[Sigma]"},muggGamGamT1\[Sigma]],Join[{"muVBFT1\[Sigma]"},muVBFGamGamT1\[Sigma]],Join[{"muTotT1\[Sigma]"},muTotGamGamT1\[Sigma]],Join[{"sigmagghT1\[Sigma]"},sigmaBRgghGamGamT1\[Sigma]],Join[{"sigmaVBF_VFT1\[Sigma]"},sigmaBRvbfvhGamGamT1\[Sigma]],Join[{"sigma_tot_T1\[Sigma]"},sigmaBRhGamGamT1\[Sigma]],Join[{"m12T1\[Sigma]"},Sqrt[m122T1\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_T_1sigma_limits.txt"],GnuTabT1\[Sigma],"Table"]

GnuTabT2\[Sigma]=Transpose[{Join[{"mh1T2\[Sigma]"},mh1T2\[Sigma]],Join[{"mHT2\[Sigma]"},mHT2\[Sigma]],Join[{"mAT2\[Sigma]"},mAT2\[Sigma]],Join[{"mHpmT2\[Sigma]"},mHpmT2\[Sigma]],Join[{"kvT2\[Sigma]"},kvT2\[Sigma]],Join[{"ktT2\[Sigma]"},ktT2\[Sigma]],Join[{"kbT2\[Sigma]"},kbT2\[Sigma]],Join[{"tanBT2\[Sigma]"},tanBT2\[Sigma]],Join[{"sinBAT2\[Sigma]"},sinBAT2\[Sigma]],Join[{"mugghT2\[Sigma]"},muggGamGamT2\[Sigma]],Join[{"muVBFT2\[Sigma]"},muVBFGamGamT2\[Sigma]],Join[{"muTotT2\[Sigma]"},muTotGamGamT2\[Sigma]],Join[{"sigmagghT2\[Sigma]"},sigmaBRgghGamGamT2\[Sigma]],Join[{"sigmaVBF_VFT2\[Sigma]"},sigmaBRvbfvhGamGamT2\[Sigma]],Join[{"sigma_tot_T2\[Sigma]"},sigmaBRhGamGamT2\[Sigma]],Join[{"m12T2\[Sigma]"},Sqrt[m122T2\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_T_2sigma_limits.txt"],GnuTabT2\[Sigma],"Table"]


GnuTabU1\[Sigma]=Transpose[{Join[{"mh1U1\[Sigma]"},mh1U1\[Sigma]],Join[{"mHU1\[Sigma]"},mHU1\[Sigma]],Join[{"mAU1\[Sigma]"},mAU1\[Sigma]],Join[{"mHpmU1\[Sigma]"},mHpmU1\[Sigma]],Join[{"kvU1\[Sigma]"},kvU1\[Sigma]],Join[{"ktU1\[Sigma]"},ktU1\[Sigma]],Join[{"kbU1\[Sigma]"},kbU1\[Sigma]],Join[{"tanBU1\[Sigma]"},tanBU1\[Sigma]],Join[{"sinBAU1\[Sigma]"},sinBAU1\[Sigma]],Join[{"mugghU1\[Sigma]"},muggGamGamU1\[Sigma]],Join[{"muVBFU1\[Sigma]"},muVBFGamGamU1\[Sigma]],Join[{"muTotU1\[Sigma]"},muTotGamGamU1\[Sigma]],Join[{"sigmagghU1\[Sigma]"},sigmaBRgghGamGamU1\[Sigma]],Join[{"sigmaVBF_VFU1\[Sigma]"},sigmaBRvbfvhGamGamU1\[Sigma]],Join[{"sigma_tot_U1\[Sigma]"},sigmaBRhGamGamU1\[Sigma]],Join[{"m12U1\[Sigma]"},Sqrt[m122U1\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_U_1sigma_limits.txt"],GnuTabU1\[Sigma],"Table"]

GnuTabU2\[Sigma]=Transpose[{Join[{"mh1U2\[Sigma]"},mh1U2\[Sigma]],Join[{"mHU2\[Sigma]"},mHU2\[Sigma]],Join[{"mAU2\[Sigma]"},mAU2\[Sigma]],Join[{"mHpmU2\[Sigma]"},mHpmU2\[Sigma]],Join[{"kvU2\[Sigma]"},kvU2\[Sigma]],Join[{"ktU2\[Sigma]"},ktU2\[Sigma]],Join[{"kbU2\[Sigma]"},kbU2\[Sigma]],Join[{"tanBU2\[Sigma]"},tanBU2\[Sigma]],Join[{"sinBAU2\[Sigma]"},sinBAU2\[Sigma]],Join[{"mugghU2\[Sigma]"},muggGamGamU2\[Sigma]],Join[{"muVBFU2\[Sigma]"},muVBFGamGamU2\[Sigma]],Join[{"muTotU2\[Sigma]"},muTotGamGamU2\[Sigma]],Join[{"sigmagghU2\[Sigma]"},sigmaBRgghGamGamU2\[Sigma]],Join[{"sigmaVBF_VFU2\[Sigma]"},sigmaBRvbfvhGamGamU2\[Sigma]],Join[{"sigma_tot_U2\[Sigma]"},sigmaBRhGamGamU2\[Sigma]],Join[{"m12U2\[Sigma]"},Sqrt[m122U2\[Sigma]]]}];

Export[StringJoin[txt,dirOutfile,"1-GnuPlotFiles/TH_U_2sigma_limits.txt"],GnuTabU2\[Sigma],"Table"]
Print["Conversion Done"]




(* calculations to obtain the whole set of parameters for a given point*)



(*we want to find all the properties of some points having Subscript[\[Mu], gg\[Rule]h\[Rule]\[Gamma]\[Gamma]] above 0.8 and mh1 between 100GeV and 110 GeV *)
(*tabHistoryPointsLHCMu={{"\[Micro]_ggh>0.8","sigma_ggh*BR_gamgam","100<mh<110","mH","mA","mH+","tanB","sinBA"}}; (* Subscript[\[Mu], gg\[Rule]h\[Rule]\[Gamma]\[Gamma]], Subscript[\[Sigma], ggh]Subscript[BR, \[Gamma]\[Gamma]], mh,mH,mA,mH+,tan\[Beta],sin(\[Beta]-\[Alpha]) *)
For[iloop=1,iloop<=sizeLHC,iloop++,If[100<=mh1LHC[[iloop]]<=110&&0.8<=muggGamGamLHC[[iloop]],tabHistoryPointsLHCMu=Join[tabHistoryPointsLHCMu,{{muggGamGamLHC[[iloop]],sigmaBRgghGamGamLHC[[iloop]],mh1LHC[[iloop]],mHLHC[[iloop]],mALHC[[iloop]],mHpmLHC[[iloop]],tanBLHC[[iloop]],sinBALHC[[iloop]]}}]]]

Export[StringJoin[txt,dirOutfile,"history_point_muggh_LHC.txt"],tabHistoryPointsLHCMu,"Table"]

(*same thing but with restrictions on \[Sigma]_ggh and not \[Mu] *)

tabHistoryPointsLHC\[Sigma]={{"\[Micro]_ggh","sigma_ggh*BR_gamgam>0.5","100<mh<110","mH","mA","mH+","tanB","sinBA"}}; (*Subscript[\[Mu], gg\[Rule]h\[Rule]\[Gamma]\[Gamma]],Subscript[\[Sigma], ggh]Subscript[BR, \[Gamma]\[Gamma]],mh,mH,mA,mH+,tan\[Beta],sin(\[Beta]-\[Alpha]) *)
For[iloop=1,iloop<=sizeLHC,iloop++,If[100<=mh1LHC[[iloop]]<=110&&0.05<=sigmaBRgghGamGamLHC[[iloop]],tabHistoryPointsLHC\[Sigma]=Join[tabHistoryPointsLHC\[Sigma],{{muggGamGamLHC[[iloop]],sigmaBRgghGamGamLHC[[iloop]],mh1LHC[[iloop]],mHLHC[[iloop]],mALHC[[iloop]],mHpmLHC[[iloop]],tanBLHC[[iloop]],sinBALHC[[iloop]]}}]]]

Export[StringJoin[txt,dirOutfile,"history_point_sigma_ggh_LHC.txt"],tabHistoryPointsLHC\[Sigma],"Table"]*)


End[]
 
EndPackage[]







