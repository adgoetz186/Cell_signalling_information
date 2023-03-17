(* ::Package:: *)

(* ::Title:: *)
(*Basic functions*)


(* ::Chapter::Closed:: *)
(*Transformation of variables and moments*)


makeMomsTime=Subscript[y, i__]->Subscript[y, i][t];
makeMomsTimeDiff=Subscript[y, i__]->Subscript[y, i]'[t];
makeMomsNoTime=Subscript[y, i__][t]->Subscript[y, i];


makeCentralMomsTime=Subscript[z, i__]->Subscript[z, i][t];
makeCentralMomsTimeDiff=Subscript[z, i__]->Subscript[z, i]'[t];


trafoMoms[nSS_]:=Moment[i_]:>Subscript@@Prepend[Flatten[Table[m,{m,nSS},{n,i[[m]]}]],y]


trafoCentralMoms[nSS_]:={CentralMoment[i_]:>Subscript@@Prepend[Flatten[Table[m,{m,nSS},{n,i[[m]]}]],z],Moment[i_]:>Subscript@@Prepend[Flatten[Table[m,{m,nSS},{n,i[[m]]}]],z]}


repCentrMom[nSS_]:=Subscript[y, i__]:>If[Length[List[i]]==1,Subscript[z, i],MomentConvert[Moment[BinCounts[{i},{1,nSS+1,1}]],"CentralMoment"]]


centralMomentToMoment=Subscript[z, i__]:>If[Length[List[i]]==1,Subscript[y, i],MomentConvert[CentralMoment[BinCounts[{i},{1,nS+1,1}]],"Moment"]];


cumulantToMoment=Subscript[c, i__]:>MomentConvert[Cumulant[BinCounts[{i},{1,nS+1,1}]],"Moment"]/.trafoMoms[nS];


(* ::Chapter::Closed:: *)
(*Variable and moment lists*)


vars[nSS_]:=Table[Subscript[x, i],{i,nSS}]


tupList[order_,nSS_]:=DeleteDuplicates[Map[Sort[#]&,Tuples[Range[nSS],{order}]]]


momsList[order_,nSS_]:=Flatten[Table[Subscript@@Prepend[tup,y],{ord,order},{tup,tupList[ord,nSS]}]]


momsCentralList[order_,nSS_]:=momsList[order,nSS]/.Subscript[y, i__]->Subscript[z, i]
momsCentralListTime[order_,nSS_]:=momsCentralList[order,nSS]/.makeCentralMomsTime
momsCentralListTimeDiff[order_,nSS_]:=momsCentralList[order,nSS]/.makeCentralMomsTimeDiff


momsCentralListDiff[order_,nSS_]:=momsCentralList[order,nSS]/.Subscript[z, i__]->Subscript[z, i]'


(* ::Chapter::Closed:: *)
(*Moment equations from CME*)


momentEquation[stoch_,prop_,ind_]:=Block[{nR=Length[prop]},Expand[Sum[Product[(Subscript[x, i]+stoch[[i,r]]),{i,ind}]prop[[r]],{r,nR}]-Sum[Product[Subscript[x, i],{i,ind}]prop[[r]],{r,nR}]]]


momentCentralEquation[stoch_,prop_,ind_]:=Block[{nR=Length[prop],indtraf},
indtraf=BinCounts[ind,{1,nS+1,1}];
If[Length[ind]==1,
momentEquation[stoch,prop,ind],
D[MomentConvert[CentralMoment[indtraf],"Moment"]/.trafoMoms[nS]/.makeMomsTime,t]/.makeMomsNoTime/.Subscript[y, i__]'[t]:>(momentEquation[stoch,prop,{i}])]]


momentCentralEquationList[order_,stoch_,prop_]:=Flatten[Table[momentCentralEquation[stoch,prop,tup],{ord,order},{tup,tupList[ord,Length[stoch]]}]]


(* ::Chapter::Closed:: *)
(*Averaging of expressions*)


integrRule[rul_,nSS_]:=Block[{trans},
trans=Flatten[Table[i,{i,nSS},{j,First[rul][[i]]}]];
If[Length[trans]==0,1,
Subscript@@Prepend[trans,y]]]


integrate[term_,vars_,nSS_]:=Block[{rule},
rule=CoefficientRules[term,vars];
Sum[Last[rule[[i]]] integrRule[rule[[i]],nSS] ,{i,Length[rule]}]]


(* ::Title:: *)
(*Closure*)


(* ::Chapter:: *)
(*Define different moment-closure methods*)


clNumbToString[numb_,ord_,nSS_]:=Which[
numb==2,poisson[ord,nSS],
numb==1,normal[ord,nSS],
numb==3,lognormal[ord,nSS],
numb==4,cmn[ord,nSS],
numb==-1,definedMA]


clStringToNumb[method_]:=Which[
method=="Poisson",2,
method=="normal",1,
method=="log-normal",3,
method=="CMN",4,
method=="self-defined",-1]


(* ::Subchapter:: *)
(*Normal and Poisson*)


normal[ord_,nSS_]:=Subscript[y, i__]:>If[Length[{i}]!=ord,Subscript[y, i],
Expand[Subscript[y, i]-MomentConvert[Cumulant[BinCounts[{i},{1,nSS+1,1}]],"Moment"]]]


cmn[ord_,nSS_]:=Subscript[y, i__]:>If[Length[{i}]!=ord,Subscript[y, i],
Expand[Subscript[y, i]-MomentConvert[CentralMoment[BinCounts[{i},{1,nSS+1,1}]],"Moment"]]]


poisson[ord_,nSS_]:=Subscript[y, i__]:>If[Length[{i}]!=ord,Subscript[y, i],If[Length[Gather[{i}]]==1,Expand[
\!\(\*SubscriptBox[\(y\), \({i}[\([\)\(1\)\(]\)]\)]\) +Subscript[y, i]-MomentConvert[Cumulant[BinCounts[{i},{1,nSS+1,1}]],"Moment"]],
Expand[Subscript[y, i]-MomentConvert[Cumulant[BinCounts[{i},{1,nSS+1,1}]],"Moment"]]]]


(* ::Subchapter:: *)
(*Lognormal*)


lnSigDiag[i_]:=Log[1+(Subscript[y, i,i]-Subscript[y, i]^2)/Subscript[y, i]^2]


lnMuEl[i_]:=Log[Subscript[y, i]]-1/2 lnSigDiag[i]


lnMu[dim_]:=Map[lnMuEl[#]&,Range[dim]]


lnSig[dim_]:=Simplify[Table[Log[(Subscript[y, i,j]-Subscript[y, i] Subscript[y, j])/Exp[lnMuEl[i]+lnMuEl[j]+1/2 (lnSigDiag[i]+lnSigDiag[j])]+1],{i,dim},{j,dim}]]


lognormal[ord_,nSS_]:=Subscript[y, i__]:>Block[{in=BinCounts[{i},{1,nSS+1,1}]},If[Length[{i}]!=ord,Subscript[y, i],Exp[in.lnMu[nSS]+1/2 in.lnSig[nSS].in]]]


(* ::Chapter::Closed:: *)
(*Perform moment-closure*)


performMA[method_,order_,term_,vars_,nSS_]:=
Expand[Simplify[Fold[ReleaseHold[#1/.clNumbToString[method,#2,nSS]]/.trafoMoms[nSS]&,integrate[term,vars,nSS],Range[Max[Map[Total[#[[1]]]&,CoefficientRules[term,vars]]],order+1,-1]]]]


opsMAcentral[method_,order_,stoch_,prop_]:=Map[performMA[method,order,#,vars[Length[stoch]],Length[stoch]]&,momentCentralEquationList[order,stoch,prop]]/.repCentrMom[Length[stoch]]/.trafoCentralMoms[Length[stoch]];


opsMAcentralNE[method_,order_,so_,stoch_,prop_]:=Block[{expand,rules,rulestruncate,propseries},
expand=Normal[Series@@Join[{prop},Transpose[{vars[nS],momsList[1,nS],ConstantArray[so,nS]}]]];
rules=CoefficientRules[expand,vars[nS]];
rulestruncate=Map[Select[#,Total[First[#]]<=so&]&,rules];
propseries=Map[Sum[i[[2]]Times@@(vars[nS]^i[[1]]),{i,#}]&,rulestruncate];
Expand@opsMAcentral[method,order,stoch,propseries]]


eqMAcentralNE[method_?NumberQ,order_?NumberQ,so_?NumberQ,stoch_,prop_]:=
Block[{ops,dmoms},
ops=opsMAcentralNE[method,order,so,stoch,prop];
dmoms=momsCentralListDiff[order,nS];
Table[dmoms[[i]]==ops[[i]],{i,Length[dmoms]}]]


steadyMAcentralNE[method_,order_,so_,stoch_,prop_]:=
Map[0==#&,opsMAcentralNE[method,order,so,stoch,prop]]


timeMAcentralNE[method_,order_,so_,stoch_,prop_]:=Block[{ops,dmoms},
ops=opsMAcentralNE[method,order,so,stoch,prop]/.makeCentralMomsTime;
dmoms=momsCentralListTimeDiff[order,nS];
Table[dmoms[[i]]==ops[[i]],{i,Length[dmoms]}]]


(* ::Title::Closed:: *)
(*Evaluation of MAs*)


(* ::Chapter:: *)
(*Variable lists*)


posIndList[ord_,nSS_]:=Map[Flatten[#]&,Map[# {1,1}&,DeleteDuplicates[Map[Sort[#]&,Tuples[Table[Range[nSS],{ii,ord}]]]],{2}],1]


momsPositive[order_,nSS_]:=Flatten[{Table[Subscript[z, i],{i,nSS}],Table[Subscript@@Prepend[sub,z],{ord2,2,order,2},{sub,posIndList[ord2/2,nSS]}]}]


(* ::Chapter:: *)
(*Fixed points and eigenvalues*)


jacobSysMANE[method_,order_,s0_,stoch_,prop_]:=Table[D[f,v],{f,opsMAcentralNE[method,order,s0,stoch,prop]},{v,momsCentralList[order,Length[stoch]]}]


fpsMAallNE[method_,order_,s0_,stoch_,prop_] := NSolve[steadyMAcentralNE[method,order,s0,stoch,prop], momsCentralList[order,Length[stoch]]]


fpsMAposNE[method_,order_,s0_,stoch_,prop_] :=  Select[fpsMAallNE[method,order,s0,stoch,prop], And @@ NonNegative[momsPositive[order,Length[stoch]] /. #] &]


(* ::Title:: *)
(*GUI*)


(* ::Chapter:: *)
(*Self defined moment closure*)


MakeSingleRule[def_]:=Block[{ind,pat},
ind=Drop[List@@First[def],1];
pat=Map[If[NumberQ[#],#,Pattern[#,_]]&,ind];
(*pat=Map[Pattern[#,_]&,ind];*)
Subscript@@Prepend[pat,y]->Subscript@@Prepend[ind,y]-MomentConvert[Cumulant[Hold[BinCounts][ind,{1,nS+1,1}]],"Moment"]+Last[def]/.cumulantToMoment/.centralMomentToMoment/.trafoMoms[nS]
]


definedMA:=Map[MakeSingleRule[#]&,defineMomentclosure]


(* ::Chapter:: *)
(*For SS analysis*)


fpsMAposStableManipulateNE[method_,order_,so_,prop_] = Hold[Select[fpsMAposNE[method,order,so,stochMatrix,prop], And @@ Negative[Re[Eigenvalues[jacobSysMANE[method,order,so,stochMatrix,prop] /. #]]] &]];


fpsMAposStableManipulatePrint[method_,order_,series_,prop_] = Hold[Block[{res},res=Select[fpsMAposNE[method,order,series,stochMatrix,prop], And @@ Negative[Re[Eigenvalues[jacobSysMANE[method,order,series,stochMatrix,prop] /. #]]] &];
Map[NumberForm[Chop@#,4]&,If[Length[res]==0,ConstantArray["-",Length[momsCentralList[order,nS]]],If[Length[res]==1,momsCentralList[order,nS]/.res[[1]],Transpose[momsCentralList[order,nS]/.res]]]]]];


(* ::Chapter:: *)
(*For time trajectories*)


ndsolveNE[method_,order_,so_,prop_,ini_, tf_] :=ndsolveNE[method,order,so,prop,ini, tf]=NDSolve[Flatten[{Expand[timeMAcentralNE[method,order,so,stochMatrix,prop]], (momsCentralListTime[order,nS] /. t -> 0) == Flatten[{ini, 0 Range[Length[momsCentralList[order,nS]]-nS]}]}], momsCentralList[order,nS], {t, 0, tf}]


ndsolveMAplotManipulateNE[method_,order_,series_,prop_,ini_, tf_,plotorder_] = Hold[Plot[Evaluate[momsCentralListTime[plotorder,nS]/.ndsolveNE[method,order,series,prop,ini, tf]/.t->tt],{tt,0,tf},
PlotLegends->momsCentralList[plotorder,nS],Frame->True,FrameLabel->{Style["time",17],Style["central moments",17]},
ImageSize->300,
BaseStyle->{FontSize->13},
PlotRange->All]];


ndsolveMAplotDataManipulateNE[method_,order_,series_,prop_,ini_, tf_,dt_] = Hold[Block[{res,tlist},
tlist=Range[0,tf,dt];
res=momsCentralListTime[order,nS]/.ndsolveNE[method,order,series,prop,ini, tf][[1]];
Map[NumberForm[Chop[#],4]&,Transpose@Join[{Prepend[momsCentralList[order,nS]/.Subscript[z, i__]:>"z"<>StringJoin@Map[ToString@#&,{i}],t]},Map[Prepend[res,#]/.t->#&,tlist]],{2}]
]];


(* ::Title::Closed:: *)
(*Interface*)


DeriveEquations:=Manipulate[
Evaluate[Column[eqMAcentralNE[clStringToNumb[method],order,expansion,stochMatrix,propensity
]]],
Style["Method specifications",Bold,Medium],
Evaluate[If[ValueQ[defineMomentclosure],
{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN","self-defined"}},{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN"}}]],
{{order,2,"closure order"}, Range[2, 5]},
{{expansion,2,"expansion order"}, Range[ 10]},
Delimiter,
Style["Parameters",Bold,Medium],
Evaluate[Sequence@@Map[{#,ToString[#]}&,parameters]],
ControllerLinking->False,
ContinuousAction->None];


SteadyState:=ReleaseHold[Manipulate@@Join[
{Column@fpsMAposStableManipulateNE[clStringToNumb[method],order,expansion,propensity]},
{
Style["Method specifications",Bold,Medium],
Evaluate[If[ValueQ[defineMomentclosure],
{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN","self-defined"}},{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN"}}]],
{{order,2,"closure order"}, Range[2, 5]},
{{expansion,2,"expansion order"}, Range[ 10]},
Delimiter,
Style["Parameters",Bold,Medium],
Evaluate[Sequence@@Table[{parameters[[i]],1},{i,Length[parameters]}]],
ControllerLinking->False,
ContinuousAction->None}]];


SteadyStateVaryParameter:=Manipulate[ReleaseHold@With[{min=parMin,max=parMax,sp=parSpacing,vary=parameter,column=Prepend[Range[parMin,parMax,parSpacing],parameter]},
Manipulate@@Join[
{Column[{Grid[Hold[Transpose]@Hold[Join][{column},Hold[Transpose]@Join[{Hold@(momsCentralList[order,nS])},Map[fpsMAposStableManipulatePrint[clStringToNumb[method],order,series,propensity/.vary->#]&,Range[min,max,sp]]]],Alignment->Left,Spacings->{2,1},Frame->All,Background->{{Gray,None},{LightGray,None}}]
,
Hold[Button]["Save",Hold[Export][ToString[file]<>".txt",Hold[Transpose]@Hold[Join][{column},Hold[Transpose]@Join[{Hold@(momsCentralList[order,nS]/.Subscript[z, i__]:>"z"<>StringJoin@Map[ToString@#&,{i}])},Map[fpsMAposStableManipulatePrint[clStringToNumb[method],order,series,propensity/.vary->#]&,Range[min,max,sp]]]],"CSV"]]}]},
{
Style["Method specifications",Bold,Medium],
Evaluate[If[ValueQ[defineMomentclosure],
{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN","self-defined"}},{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN"}}]],
{{order,2,"closure order"}, Range[2, 5]},
{{series,2,"expansion order"}, Range[ 10]},
Delimiter,
Style["Fixed parameters",Bold,Medium],
Evaluate[Sequence@@Table[{Complement[parameters,{parameter}][[i]],1},{i,Length[parameters]-1}]]
,
ControllerLinking->False,
ContinuousAction->None}]]
,
Style["Paremeter scan specification",Bold,Medium],
{{parameter,parameters[[1]],"vary parameter"},parameters},
{{parMin,1,"minimal value"},1},
{{parMax,3,"maximal value"},3},
{{parSpacing,1,"grid spacing"},1},
Delimiter,
Style["Export results",Bold,Medium],
{{file,filename,"exported file's name"},filename}];


TimeTrajectory:=Manipulate[ReleaseHold@With[{},
Manipulate@@Join[
{Column[{ndsolveMAplotManipulateNE[clStringToNumb[method],order,expansion,propensity,momsCentralList[1,nS],tfinal,plotorder],
Row[{Hold[Button]["Save figure",Hold[Export][ToString[file]<>"."<>ToString[format],ndsolveMAplotManipulateNE[clStringToNumb[method],order,expansion,propensity,momsCentralList[1,nS],tfinal,plotorder]]],Hold[Button]["Save data",Hold[Export][ToString[file]<>".txt",Transpose@ndsolveMAplotDataManipulateNE[clStringToNumb[method],order,expansion,propensity,momsCentralList[1,nS],tfinal,dt],"CSV"]]}]},Alignment->Center]},
{
Style["Method specifications",Bold,Medium],
Evaluate[If[ValueQ[defineMomentclosure],
{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN","self-defined"}},{{method,"normal","closure method"},{"normal","Poisson","log-normal","CMN"}}]],
{{order,2,"closure order"}, Range[2, 5]},
{{expansion,2,"expansion order"}, Range[ 10]},
Delimiter,
Style["Parameters",Bold,Medium],
Evaluate[Sequence@@Map[{#,1}&,parameters]],
Delimiter,
Style["Initial conditions",Bold,Medium],
Evaluate[Sequence@@Map[{#,1}&,momsCentralList[1,nS]]],
Delimiter,
Style["Plot specifications",Bold,Medium],
{{tfinal,1,"final time"},1},
{{plotorder,1,"plot order"},Range[1, 5]}
,
ControllerLinking->False,
ContinuousAction->None}]]
,
Style["Export figure",Bold,Medium],
{{file,filename,"exported file's name"},filename},
{{format,gif,"export format"},{gif,jpg,eps,pdf}},
{{dt,0.1,"time spacing"},0.1}];


(* ::Title::Closed:: *)
(*Self coding*)


DeriveEquationsSelf[method_?NumberQ,order_?NumberQ,so_?NumberQ,pars_]:=
Block[{ops,dmoms},
ops=opsMAcentralNE[method,order,so,stochMatrix,propensityFct@@pars];
dmoms=momsCentralListDiff[order,nS];
Table[dmoms[[i]]==ops[[i]],{i,Length[dmoms]}]]


SteadyStateSelf[method_,order_,so_,pars_]:= Select[fpsMAposNE[method,order,so,stochMatrix,propensityFct@@pars], And @@ Negative[Re[Eigenvalues[jacobSysMANE[method,order,so,stochMatrix,propensityFct@@pars] /. #]]] &]


SolveMAsInTimeSelf[method_,order_,so_,pars_,ini_, tf_] :=NDSolve[Flatten[{Expand[timeMAcentralNE[method,order,so,stochMatrix,propensityFct@@pars]], (momsCentralListTime[order,nS] /. t -> 0) == Flatten[{ini, 0 Range[Length[momsCentralList[order,nS]]-Length[ini]]}]}], momsCentralList[order,nS], {t, 0, tf}]


PlotTimeTrajectorySelf[method_,order_,series_,pars_,ini_, tf_,plotorder_] := 
Block[{res},
res=NDSolve[Flatten[{Expand[timeMAcentralNE[method,order,series,stochMatrix,propensityFct@@pars]], (momsCentralListTime[order,nS] /. t -> 0) == Flatten[{ini, 0 Range[Length[momsCentralList[order,nS]]-Length[ini]]}]}], momsCentralList[order,nS], {t, 0, tf}];
Plot[Evaluate[momsCentralListTime[plotorder,nS]/.res[[1]]/.t->tt],{tt,0,tf},
PlotLegends->momsCentralList[plotorder,nS],Frame->True,FrameLabel->{Style["time",17],Style["central moments",17]},
ImageSize->300,
BaseStyle->{FontSize->13},
PlotRange->All]
]


(* ::Title::Closed:: *)
(*Rate equations*)


(* ::Chapter:: *)
(*Vars*)


makeREtime=Subscript[z, i_]->Subscript[z, i][t];


varsRE[nS_]:=Table[Subscript[z, i],{i,nS}];
varsREtime[nS_]:=Table[Subscript[z, i][t],{i,nS}];
varsREtimeDiff[nS_]:=Table[Subscript[z, i]'[t],{i,nS}];
varsREDiff[nS_]:=Table[Subscript[z, i]',{i,nS}];


(* ::Chapter:: *)
(*Equations*)


opsRE[stoch_,propMac_]:=Flatten[stoch.propMac]


ssRE[stoch_,propMac_]:=Block[{dr=opsRE[stoch,propMac]},Table[0==dr[[i]],{i,Length[stoch]}]]


timeRE[stoch_,propMac_]:=Block[{dvars=varsREtimeDiff[Length[stoch]],dr=opsRE[stoch,propMac]/.makeREtime},Table[dvars[[i]]==dr[[i]],{i,Length[stoch]}]]


(* ::Chapter:: *)
(*FP and eigenvalues*)


jacobSysRE[stoch_,propMac_]:=Table[D[f,v],{f,opsRE[stoch,propMac]},{v,varsRE[Length[stoch]]}]


fpsREall[stoch_,propMac_]:=fpsREall[stoch,propMac]=NSolve[ssRE[stoch,propMac],varsRE[Length[stoch]]]


fpsREreal[stoch_,propMac_]:=Select[fpsREall[stoch,propMac],And@@(Im[varsRE[Length[stoch]]/.#]==0 varsRE[Length[stoch]])&]


fpsREpos[stoch_,propMac_]:=Select[fpsREall[stoch,propMac],And@@Positive[varsRE[Length[stoch]]/.#]&]


fpsREposStable[stoch_,propMac_] := Select[fpsREpos[stoch,propMac], And @@ Negative[Re[Eigenvalues[jacobSysRE[stoch,propMac] /. #]]] &]


eigsREall[stoch_,propMac_]:=Map[Eigenvalues[jacobSysRE[stoch,propMac]/.#]&,fpsREall[stoch,propMac]]


eigsREposStable[stoch_,propMac_]:=Select[Map[Eigenvalues[jacobSysRE[stoch,propMac]/.#]&,fpsREpos[stoch,propMac]],And@@Negative[Re[#]]&]


eigsRErealStable[stoch_,propMac_]:=Select[Map[Eigenvalues[jacobSysRE[stoch,propMac]/.#]&,fpsREreal[stoch,propMac]],And@@Negative[Re[#]]&]


(* ::Chapter:: *)
(*Solve in time*)


solveTimeRE[stoch_,propMac_,ini_,tf_]:=NDSolve[Flatten[{timeRE[stoch,propMac],(varsREtime[Length[stoch]]/.t->0)==ini}],varsRE[Length[stoch]],{t,0,tf}]


(* ::Chapter:: *)
(*For manipulate*)


eqRE[stoch_,propRE_]:=Block[{dvars=varsREDiff[Length[stoch]],dr=opsRE[stoch,propRE]},Table[dvars[[i]]==dr[[i]],{i,Length[stoch]}]]


fpsREposStableManipulate[prop_] = Hold[Select[fpsREpos[stochMatrix,prop], And @@ Negative[Re[Eigenvalues[jacobSysRE[stochMatrix,prop] /. #]]] &]];


fpsREposStableManipulatePrint[prop_] = Hold[Block[{res},res=Select[fpsREposStable[stochMatrix,prop], And @@ Negative[Re[Eigenvalues[jacobSysRE[stochMatrix,prop] /. #]]] &];
Map[NumberForm[Chop@#,4]&,If[Length[res]==0,ConstantArray[res,nS],If[Length[res]==1,
varsRE[nS]/.res[[1]],Transpose[varsRE[nS]/.res]]]]]];


ndsolveREplotManipulate[prop_,ini_, tf_] = Hold[Plot[Evaluate[varsREtime[nS]/.solveTimeRE[stochMatrix,prop,ini,tf]/.t->tt],{tt,0,tf},
PlotLegends->varsRE[nS],Frame->True,FrameLabel->{Style["time",17],Style["central moments",17]},
ImageSize->300,
BaseStyle->{FontSize->13},
PlotRange->All]];


ndsolveMAplotManipulateNE[method_,order_,series_,prop_,ini_, tf_,plotorder_] = Hold[Plot[Evaluate[momsCentralListTime[plotorder,nS]/.ndsolveNE[method,order,series,prop,ini, tf]/.t->tt],{tt,0,tf},
PlotLegends->momsCentralList[plotorder,nS],Frame->True,FrameLabel->{Style["time",17],Style["central moments",17]},
ImageSize->300,
BaseStyle->{FontSize->13},
PlotRange->All]];


ndsolveREplotDataManipulate[prop_,ini_, tf_,dt_] = Hold[Block[{res,tlist},
tlist=Range[0,tf,dt];
res=varsREtime[nS]/.solveTimeRE[stochMatrix,prop,ini,tf][[1]];
Map[NumberForm[Chop[#],4]&,Transpose@Join[{Prepend[varsRE[nS]/.Subscript[z, i__]:>"z"<>StringJoin@Map[ToString@#&,{i}],t]},Map[Prepend[res,#]/.t->#&,tlist]],{2}]
]];


ndsolveMAplotDataManipulateNE[method_,order_,series_,prop_,ini_, tf_,dt_] = Hold[Block[{res,tlist},
tlist=Range[0,tf,dt];
res=momsCentralListTime[order,nS]/.ndsolveNE[method,order,series,prop,ini, tf][[1]];
Map[NumberForm[Chop[#],4]&,Transpose@Join[{Prepend[momsCentralList[order,nS]/.Subscript[z, i__]:>"z"<>StringJoin@Map[ToString@#&,{i}],t]},Map[Prepend[res,#]/.t->#&,tlist]],{2}]
]];


(* ::Chapter:: *)
(*Manipulate*)


DeriveEquationsRE:=Manipulate[
Evaluate[Column[eqRE[stochMatrix,propensityRE]]],
Delimiter,
Style["Parameters",Bold,Medium],
Evaluate[Sequence@@Map[{#,ToString[#]}&,parameters]],
ControllerLinking->False,
ContinuousAction->None];


SteadyStateRE:=ReleaseHold[Manipulate@@Join[
{Column@fpsREposStableManipulate[propensityRE]},
{
Delimiter,
Style["Parameters",Bold,Medium],
Evaluate[Sequence@@Table[{parameters[[i]],1},{i,Length[parameters]}]],
ControllerLinking->False,
ContinuousAction->None}]];


SteadyStateVaryParameterRE:=Manipulate[ReleaseHold@With[{min=parMin,max=parMax,sp=parSpacing,vary=parameter,column=Prepend[Range[parMin,parMax,parSpacing],parameter]},
Manipulate@@Join[
{Column[{Grid[Hold[Transpose]@Hold[Join][{column},Hold[Transpose]@Join[{Hold@(varsRE[nS])},Map[fpsREposStableManipulatePrint[propensityRE/.vary->#]&,Range[min,max,sp]]]],Alignment->Left,Spacings->{2,1},Frame->All,Background->{{Gray,None},{LightGray,None}}]
,
Hold[Button]["Save",Hold[Export][ToString[file]<>".txt",Hold[Transpose]@Hold[Join][{column},Hold[Transpose]@Join[{Hold@(varsRE[nS]/.Subscript[z, i__]:>"z"<>StringJoin@Map[ToString@#&,{i}])},Map[fpsREposStableManipulatePrint[propensityRE/.vary->#]&,Range[min,max,sp]]]],"CSV"]]}]},
{
Style["Parameters",Bold,Medium],
Evaluate[Sequence@@Table[{Complement[parameters,{parameter}][[i]],1},{i,Length[parameters]-1}]]
,
ControllerLinking->False,
ContinuousAction->None}]]
,
{parameter,parameters},
{parMin,1},
{parMax,2},
{parSpacing,0.5},
Delimiter,
Style["Export results",Bold,Medium],
{file,filename}(*,
ContinuousAction\[Rule]None*)];


TimeTrajectoryRE:=Manipulate[ReleaseHold@With[{},
Manipulate@@Join[
{Column[{ndsolveREplotManipulate[propensityRE,Table[Subscript[y, i],{i,nS}],tfinal],
Row[{Hold[Button]["Save figure",Hold[Export][ToString[file]<>"."<>ToString[format],ndsolveREplotManipulate[propensityRE,Table[Subscript[y, i],{i,nS}],tfinal]]],Hold[Button]["Save data",Hold[Export][ToString[file]<>".txt",Transpose@ndsolveREplotDataManipulate[propensityRE,Table[Subscript[y, i],{i,nS}],tfinal,dt],"CSV"]]}]},Alignment->Center]},
{
Style["Initial conditions",Bold,Medium],
Evaluate[Sequence@@Table[{{Subscript[y, i],1,Subscript[z, i]}},{i,nS}]],
Style["Parameters",Bold,Medium],
Evaluate[Sequence@@Map[{#,1}&,parameters]],
Style["Final time",Bold,Medium],
{tfinal,1},
ControllerLinking->False,
ContinuousAction->None}]]
,
Style["Export figure",Bold,Medium],
{file,filename},
{format,{gif,jpg,eps,pdf}},
{dt,0.1}];


(* ::Chapter:: *)
(*Specifically for self coding*)


DeriveEquationsSelfRE[pars_]:=
eqRE[stochMatrix,propensityFctRE@@pars]


SteadyStateSelfRE[pars_]:= fpsREposStable[stochMatrix,propensityFctRE@@pars]


SolveMAsInTimeSelfRE[pars_,ini_, tf_] :=solveTimeRE[stochMatrix,propensityFctRE@@pars,ini,tf]


PlotTimeTrajectorySelfRE[pars_,ini_, tf_] := Plot[Evaluate[varsREtime[nS]/.solveTimeRE[stochMatrix,propensityFctRE@@pars,ini,tf]/.t->tt],{tt,0,tf},
PlotLegends->varsRE[nS],Frame->True,FrameLabel->{Style["time",17],Style["central moments",17]},
ImageSize->300,
BaseStyle->{FontSize->13},
PlotRange->All]
