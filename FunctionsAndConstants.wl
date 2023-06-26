(* ::Package:: *)

BeginPackage["Functions`"]


(* First define conversion factors between various units*)
CPrec        =500; (* numerical precission for later calculations*)
constants=SetPrecision[{hbar->1.05457*^-34,c0->2.997925*^8,ec->1.602177*^-19},5 CPrec];
(* Note that ec is the unitless(!!) value of the electric charge, which is used to convert eV to J. ev=ec*J *)

m2invMeV  =1/(hbar c0/(1*^6*ec))/.constants; (*to be read as convert meter to MeV^-1*)
kg2MeV       =1/(ec 10^6/(c0^2))/.constants; 
s2invMeV  =1/(hbar/(1*^6 ec))/.constants; 
N2invMeV2=kg2MeV m2invMeV/s2invMeV^2; (*Newton to inverse MeV^2*)
kgm32MeV4=kg2MeV/m2invMeV^3; (*to be read as kg/m^3 converted to MeV^4 / density to natural units*)

(*Next define various constants / parameters that are use in the LLR analysis*)

rhoV=SetPrecision[1*^-15,5 CPrec]; (* vacuum density in MeV^4 *)
rhoE=SetPrecision[0.0000236851,5 CPrec]; (*mean density of the earth in MeV^4 *)
rhoS=SetPrecision[1408*kgm32MeV4,5 CPrec]; (* mean density of the sun in MeV^4 *)
rhoM=SetPrecision[0.0000143877,5 CPrec]; (*mean density of the moon in MeV^4, Pitschmannvalue*)
rhoN=SetPrecision[1.435550104*10^10,5 CPrec]; (*neutron density in MeV^4, Pitschmannvalue*)

RE=SetPrecision[6378100*m2invMeV,5 CPrec]; (*equatorial radius of the earth in MeV*)
RS=SetPrecision[695700000*m2invMeV,5 CPrec]; (*equatorial radius of the sun in MeV*)
RM=SetPrecision[1738100*m2invMeV,5 CPrec]; (*equatorial radius of the moon in MeV*)
RN=SetPrecision[0.0025,5 CPrec]; (*Neutron radius MeV^-1 - Pitschmann value*)

ME=SetPrecision[5.9724*10^(24)*kg2MeV,5 CPrec]; (*total mass of the earth in MeV*)
MS=SetPrecision[1.9891*10^(30)*kg2MeV,5 CPrec]; (*total mass of the sun in MeV*)
MN=SetPrecision[939.565346,5 CPrec]; (* neutron mass in MeV, Pitschmann value *)

REM=SetPrecision[385000558.4*m2invMeV,5 CPrec]; (*Earth-Moon distance in MeV^-1*)   
REMax=SetPrecision[405000000*m2invMeV,5 CPrec]; (*Max. Earth-Moon distance in MeV^-1*) 

\[Rho]Cosmos=2.512792571065065102211739674794535033948495662732724260254001020163981835118*10^-35;

\[Rho]V = rhoV;
\[Rho]E = rhoE;
\[Rho]S = rhoS;
\[Rho]M = rhoM;
\[Rho]N = rhoN;

\[Rho]M2= SetPrecision[2514.0`*kgm32MeV4,5 CPrec]; (*mirror density cannex*)

ddr=SetPrecision[5.0677398985414185`*^7,1000]; (*10 \[Mu]m in MeV*)


RES=SetPrecision[149597870700*m2invMeV,5 CPrec]; (*Earth-Sun distance also known as AE or AU, Subscript[r, AU] in the paper*)

aG=SetPrecision[1.30199*10^(-32),5 CPrec]; (*acceleration of earth towards sun in MeV, Pitschmann value*)

GN=SetPrecision[6.70861*^-45,5 CPrec]; (*newtons constant in MeV^(-2), Pitschmann value*)
mpl=SetPrecision[2.4353635268*^21,5 CPrec](*planck mass in MeV*)

H0 := SetPrecision[197.3269788/(1.374*10^41),5CPrec];
\[CapitalOmega] := SetPrecision[0.73,5CPrec];
VC= 3 \[CapitalOmega] (mpl^2) (H0^2)         (* SetPrecision[6.73304*^-34,5 CPrec]; *)(*in MeV^4, this is the effective potential for the vacuum density and the vaccuum expactation value of \[Phi]*)


(*See DerivationOfNumericallySaveExpressions *)
(* ProductLogRhoMSafe[\[Lambda]_?NumericQ,A2_,VC_,rhoM_]:=If[(2 VC \[Lambda]^2)/(A2 rhoV)>130000,SetPrecision[1/(60 (-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))])^5) (60 (-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))])^6-60 (-2+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]) (1+(-2+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]) (-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))])) (Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]+(-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))])^2) Log[-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]]+30 (-10+(-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]) (6+(-4+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]) (-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]))) Log[-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]]^2+10 (35+(-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]) (-13+2 Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+2 Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))])) Log[-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]]^3+5 (-28+3 Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+3 Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]) Log[-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]]^4+12 Log[-1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]+Log[(2 VC \[Lambda]^2)/(rhoM (A2+Sqrt[A2 (A2+(2 VC \[Lambda]^2)/rhoV)]))]]^5),5 CPrec],SetPrecision[ProductLog[(2 E^((2 VC \[Lambda]^2)/(A2 rhoV (1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]))) VC \[Lambda]^2)/(A2 rhoM (1+Sqrt[1+(2 VC \[Lambda]^2)/(A2 rhoV)]))],5 CPrec]];*) 

V0[A_,\[Lambda]_]:=SetPrecision[(2 VC/(Sqrt[1 + (2 \[Lambda]^2)/(A \[Rho]Cosmos) VC] +1)) Exp[SetPrecision[(2 \[Lambda]^2)/(A \[Rho]Cosmos) VC/(Sqrt[1 + (2 \[Lambda]^2)/(A \[Rho]Cosmos) VC] +1),5CPrec]],5CPrec];

phirho[A_,\[Lambda]_,\[Beta]_,\[Rho]_]:=SetPrecision[mpl/\[Lambda] If[2*Log10[\[Lambda]]+\[Beta]-Log10[A]-Log10[\[Rho]]<1000000,ProductLog[(\[Lambda]^2 10^\[Beta])/(A \[Rho])],Log[\[Lambda]^2/(A \[Rho])]+1/(6 (\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])])^3) (6 \[Beta] Log[10] (\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])])^3-6 (-1+\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]) (1+\[Beta]^2 Log[10]^2+2 \[Beta] Log[10] Log[\[Lambda]^2/(A \[Rho])]+Log[\[Lambda]^2/(A \[Rho])]^2) Log[\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]]+3 (-3+\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]) Log[\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]]^2+2 Log[\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]]^3)],1000] 

murho[\[Lambda]_,\[Rho]_,A2_]:=SetPrecision[1/mpl Sqrt[\[Lambda]^2 V0[A2,\[Lambda]]Exp[SetPrecision[-\[Lambda] phirho[\[Lambda],\[Rho],A2]/mpl,5CPrec]]+A2 \[Rho]],5CPrec];

Veff[\[Lambda]_,\[Rho]_,\[Phi]_,A2_]:=SetPrecision[V0[A2,\[Lambda]] Exp[SetPrecision[-\[Lambda] \[Phi]/mpl,5CPrec]]+(A2 \[Rho])/(2mpl^2) \[Phi]^2,5CPrec];

j[x_]:=3/x^2 (1-Tanh[x]/x)

js[x_]:=1-(2 x^2)/5+(17 x^4)/105-(62 x^6)/945+(1382 x^8)/51975

SC[\[Lambda]_,\[Rho]_,A2_,R_] :=If[murho[\[Lambda],\[Rho],A2] R<=0.01,SetPrecision[js[murho[\[Lambda],\[Rho],A2] R]/(1+murho[\[Lambda],rhoV,A2]/murho[\[Lambda],\[Rho],A2] Tanh[murho[\[Lambda],\[Rho],A2] R]),5 CPrec],SetPrecision[j[murho[\[Lambda],\[Rho],A2] R]/(1+murho[\[Lambda],rhoV,A2]/murho[\[Lambda],\[Rho],A2] Tanh[murho[\[Lambda],\[Rho],A2] R]),5 CPrec]] (*screening charge*)

Q[\[Lambda]_,\[Rho]_,A2_,R_] :=SC[\[Lambda],\[Rho],A2,R]

\[Mu][A_,\[Lambda]_,\[Beta]_,\[Rho]_] :=SetPrecision[ Sqrt[SetPrecision[\[Lambda]^2  E^SetPrecision[-\[Lambda] \[Phi]M[A,\[Lambda],\[Beta],\[Rho]]/mpl+\[Beta]/Log10[E],1000] + A \[Rho],1000]]/mpl,1000]

\[Phi]M[A_,\[Lambda]_,\[Beta]_,\[Rho]_]:=phirho[A,\[Lambda],\[Beta],\[Rho]]

\[Mu]0[\[Lambda]_?NumberQ,A2_?NumberQ,\[Beta]_?NumberQ,\[Phi]0_?NumberQ,\[Rho]V2_?NumberQ]:=SetPrecision[1/mpl Sqrt[\[Lambda]^2 E^SetPrecision[-\[Lambda] \[Phi]0/mpl+\[Beta]/Log10[E],1000]+A2 \[Rho]V2],1000]

D0[\[Lambda]_?NumberQ,A2_?NumberQ,\[Beta]_?NumberQ,\[Phi]0_?NumberQ,\[Rho]V2_?NumberQ]:=SetPrecision[\[Lambda] /mpl E^SetPrecision[-\[Lambda] \[Phi]0/mpl+\[Beta]/Log10[E],1000]-(A2 \[Rho]V2)/mpl^2 \[Phi]0,1000]

\[Phi]d[\[Lambda]_,A2_,\[Beta]_,\[Phi]0_,d_,\[Rho]V2_]:=SetPrecision[\[Phi]0+D0[\[Lambda],A2,\[Beta],\[Phi]0,\[Rho]V2]/\[Mu]0[\[Lambda],A2,\[Beta],\[Phi]0,\[Rho]V2]^2 (1-Cosh[SetPrecision[\[Mu]0[\[Lambda],A2,\[Beta],\[Phi]0,\[Rho]V2]*d,5CPrec]]),5CPrec]

\[Phi]2[\[Lambda]_,A2_,\[Beta]_,\[Phi]0_,d_,z_,\[Rho]M2_,\[Rho]V2_]:=SetPrecision[Piecewise[{{\[Phi]0+D0[\[Lambda],A2,\[Beta],\[Phi]0,\[Rho]V2]/\[Mu]0[\[Lambda],A2,\[Beta],\[Phi]0,\[Rho]V2]^2 (1-Cosh[SetPrecision[\[Mu]0[\[Lambda],A2,\[Beta],\[Phi]0,\[Rho]V2]*(z),5CPrec]]),d-Abs[z]>0},
{\[Phi]M[A2,\[Lambda],\[Beta],\[Rho]M2]+(\[Phi]d[\[Lambda],A2,\[Beta],\[Phi]0,d,\[Rho]V2]-\[Phi]M[A2,\[Lambda],\[Beta],\[Rho]M2])E^SetPrecision[-\[Mu][A2,\[Lambda],\[Beta],\[Rho]M2](Abs[z]-d),5CPrec],Abs[z]-d>0}}],5CPrec]

\[Phi]\[Rho][A_,\[Lambda]_,\[Beta]_,\[Rho]_]:=SetPrecision[mpl/\[Lambda] If[2*Log10[\[Lambda]]+\[Beta]-Log10[A]-Log10[\[Rho]]<1000000,ProductLog[(\[Lambda]^2 10^\[Beta])/(A \[Rho])],Log[\[Lambda]^2/(A \[Rho])]+1/(6 (\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])])^3) (6 \[Beta] Log[10] (\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])])^3-6 (-1+\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]) (1+\[Beta]^2 Log[10]^2+2 \[Beta] Log[10] Log[\[Lambda]^2/(A \[Rho])]+Log[\[Lambda]^2/(A \[Rho])]^2) Log[\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]]+3 (-3+\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]) Log[\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]]^2+2 Log[\[Beta] Log[10]+Log[\[Lambda]^2/(A \[Rho])]]^3)],1000]
EQN\[Phi]0[\[Lambda]_,A2_,\[Beta]_,\[Phi]0_,\[Rho]_,\[Rho]V2_,d_]:=SetPrecision[-\[Phi]0+((mpl \[Lambda]-A2 E^SetPrecision[(\[Lambda] \[Phi]0)/mpl-\[Beta] Log[10],1000] \[Rho]V2 \[Phi]0) (-1+Cosh[SetPrecision[(d Sqrt[E^SetPrecision[-((\[Lambda] \[Phi]0)/mpl)+\[Beta] Log[10],1000] \[Lambda]^2+A2 \[Rho]V2])/mpl,1000]]))/(\[Lambda]^2+A2 E^SetPrecision[(\[Lambda] \[Phi]0)/mpl-\[Beta] Log[10],1000] \[Rho]V2)+\[Phi]\[Rho][A2,\[Lambda],\[Beta] ,\[Rho]]+((E^SetPrecision[-((\[Lambda] \[Phi]0)/mpl)+\[Beta] Log[10],1000] mpl \[Lambda]-A2 \[Rho]V2 \[Phi]0) Sinh[(d Sqrt[E^SetPrecision[-((\[Lambda] \[Phi]0)/mpl)+\[Beta] Log[10],1000] \[Lambda]^2+A2 \[Rho]V2])/mpl])/(Sqrt[E^SetPrecision[-((\[Lambda] \[Phi]0)/mpl)+\[Beta] Log[10],1000] \[Lambda]^2+A2 \[Rho]V2] Sqrt[A2 \[Rho] (1+(\[Lambda] \[Phi]\[Rho][A2,\[Lambda],\[Beta] ,\[Rho]])/mpl)]),1000]
\[Phi]0[\[Lambda]_?NumericQ,A2_,\[Beta]_,d_,\[Rho]_,\[Rho]V2_]:=Block[{res,startp},
startp=Abs[SetPrecision[\[Phi]\[Rho][A2,\[Lambda],\[Beta],\[Rho]V2],5 CPrec]];
res=N[(phi0/.SetPrecision[FindRoot[Re[SetPrecision[EQN\[Phi]0[SetPrecision[\[Lambda],1000],SetPrecision[A2,1000],SetPrecision[\[Beta],1000],SetPrecision[phi0,1000],SetPrecision[\[Rho],1000],SetPrecision[\[Rho]V2,1000],SetPrecision[d,1000]],5CPrec]],{phi0,startp,0,SetPrecision[1.5,5CPrec] startp},AccuracyGoal->5CPrec,PrecisionGoal->5CPrec,WorkingPrecision->5 CPrec],5 CPrec]),5CPrec];
test=checkval[res];
If[test!=NaN,
If[Im[test]/Re[test]>1*^-5,Print["Warning: getphi0i resulted in non-negligible Complex value at \[Lambda]=",N[\[Lambda],6]," and A2=",N[A2,6]]]];
Abs[res]];


\[Phi][r_,R_,\[Lambda]_,A2_,\[Rho]s_]:=If[r<=R,SetPrecision[\[Phi]\[Rho][\[Lambda],\[Rho]s,A2]+(\[Phi]\[Rho][\[Lambda],\[Rho]V,A2]-\[Phi]\[Rho][\[Lambda],\[Rho]s,A2])/Cosh[\[Mu][\[Lambda],\[Rho]s,A2]R] (1+\[Mu][\[Lambda],\[Rho]V,A2]R)/(\[Mu][\[Lambda],\[Rho]s,A2]+\[Mu][\[Lambda],\[Rho]V,A2]Tanh[\[Mu][\[Lambda],\[Rho]s,A2]R]) Sinh[\[Mu][\[Lambda],\[Rho]s,A2]r]/r,5 CPrec],SetPrecision[\[Phi]\[Rho][\[Lambda],\[Rho]V,A2]-Q[\[Lambda],\[Rho]s,A2,R] (\[Mu][\[Lambda],\[Rho]s,A2]^2 R^3)/3 (\[Phi]\[Rho][\[Lambda],\[Rho]V,A2]-\[Phi]\[Rho][\[Lambda],\[Rho]s,A2]) E^(-\[Mu][\[Lambda],\[Rho]V,A2](r-R))/r,5CPrec]]

\[Phi]Even[r_,R_,\[Lambda]_,A2_,\[Rho]s_]:=If[r>=0,\[Phi][r,R,\[Lambda],A2,\[Rho]s],\[Phi][-r,R,\[Lambda],A2,\[Rho]s]]

EQNLLR[\[Lambda]_,A2_]:=D[D[phi[r],r],r]+2/r D[phi[r],r]==HeavisideTheta[RS-Abs[r]]murho[\[Lambda],rhoS,A2]^2 (phi[r]-phirho[\[Lambda],rhoS,A2])+HeavisideTheta[Abs[r]-RS](murho[\[Lambda],rhoV,A2]^2) (phi[r]-phirho[\[Lambda],rhoV,A2]) (*this is the linearized DE for LLR*)

EQNLLRMoon[\[Lambda]_,A2_]:=D[D[phi[r],r],r]+2/r D[phi[r],r]==HeavisideTheta[RE-Abs[r]]murho[\[Lambda],rhoE,A2]^2 (phi[r]-phirho[\[Lambda],rhoE,A2])+HeavisideTheta[Abs[r]-RE](murho[\[Lambda],rhoV,A2]^2) (phi[r]-phirho[\[Lambda],rhoV,A2])


EndPackage[];
