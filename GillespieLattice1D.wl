(* ::Package:: *)

(* ::Title:: *)
(*Public part*)


BeginPackage["GillespieLattice1D`"]


fReactionUnimol::usage="fReactionUnimol[name,kr,r,p] makes a unimolecular reaction 
with name and reaction rate kr. r is a list of length 1, with names of molecules 
as strings. p is a list of length 0, 1, or 2 with names of molecules as strings.";


fReactionUnimol::usage="fReactionUnimol[name,kr,r,p] makes a unimolecular reaction 
with name and reaction rate kr. r is a list of length 1, with names of molecules as 
strings. p is a list of length 0, 1, or 2 with names of molecules as strings.";


fReactionBimol::usage="fReactionBimol[name,prob,r,p] makes a unimolecular reaction 
with name and probability prob. r is a list of length 2, with names of molecules as 
strings. p is a list of length 0, 1, or 2 with names of molecules as strings. Every 
time the two reagents meet, they react with probability prob.";


fGillespie::usage="fGillespie[tMax,dt,lattice0,uniRxnList,biRxnList,verbose] runs the 
Gillespie code for tMax with timestep dt from an initial lattice lattice0. uniRxnList 
and biRxnList are the lists fof reactions. The verbose flag is False by default - 
warning: True is very verbose!";


latticeSpacing::usage="latticeSpacing[dc,dt] returns sqrt(2*dc*dt). This is the lattice 
spacing in units as diffusion constant. Note that the mols take a step at every iteration 
dt.";


fMean::usage="fMean[lattice,s] returns the number of species s (strings) in the lattice. 
fMean[lattice] gives the number of non-zero entries in lattice.";


fNN::usage="fNN[lattice,s1,s2] returns the number of NN s1,s2 + s2,s1 (strings) in the 
lattice. fNN[lattice] gives the total number NNs in lattice.";


fNNN::usage="fNNN[lattice,s1,s2] returns the number of NNN s1,s2 + s2,s1 (strings) in 
the lattice. fNNN[lattice] gives the total number NNNs in lattice.";


fTriplet::usage="fTriplet[lattice,s1,s2,s3] returns the number of triplets s1,s2,s3 + 
permutations (strings) in the lattice. fTriplet[lattice] gives the total number triplets 
in lattice.";


fGenLattice::usage="fGenLattice[hDict,jDict,nChain,nIter,trackMoments:False] generates a 
lattice of length nChain by Gibbs sampling. hDict is a dict from the species name 
(string) to the h value. jDict is a dict from the list of length 2 of species names to 
the j value (note: only supply A,B or B,A but not both). nIter iterations are run. 
trackMoments will also keep track of the means and NNs in the lattice, returned as 
{lattice, meanDict, nnDict}.";


fTESTrxnDictUniFromList::usage="INTERNAL TEST to generate a dictionary from a list 
of unimolecular reactions.";


fTESTrxnDictBiFromList::usage="INTERNAL TEST to generate a dictionary from a list
of bimolecular reactions.";


fTESTpredictDiff::usage="INTERNAL TEST to take a single step, where the molecules
diffuse and do the bimolecular reactions.";


fTESTscheduleUnimolRxns::usage="INTERNAL TEST to schedule unimolecular reactions.";


fTESTdoUniRxn::usage="INTERNAL TEST to do the scheduled unimolecular reaction.";


fTESTmoments::usage="INTERNAL TEST to calculate moments from a lattice.";


fTEST1DGenLattice::usage="INTERNAL TEST to generate a lattice by Gibbs sampling. 
One species.";


fTEST2DGenLattice::usage="INTERNAL TEST to generate a lattice by Gibbs sampling. 
Two species.";


(* ::Chapter:: *)
(*Private part*)


Begin["Private`"]


(* ::Subsection:: *)
(*Function to make a reaction structure*)


(* ::Input::Initialization:: *)
fReactionUnimol[name_,kr_,r_,p_]:=Return[<|"name"->name,"kr"->kr,"r"->r,"p"->p|>];


(* ::Input::Initialization:: *)
fReactionBimol[name_,prob_,r_,p_]:=Return[<|"name"->name,"prob"->prob,"r"->r,"p"->p|>];


(* ::Subsection:: *)
(*Reaction dict of unimol reactions from list*)


(* ::Input::Initialization:: *)
rxnDictUniFromList[uniRxnList_]:=Module[
{uniRxnDict,reactantName}
,
(* Go over all reactions *)
uniRxnDict=<||>;
Do[
reactant=rxn["r"];
If[Length[reactant]==0,reactantName={};];
If[Length[reactant]==1,reactantName=rxn["r"][[1]];];
If[
KeyExistsQ[uniRxnDict,reactantName]
,
AppendTo[uniRxnDict[reactantName],rxn]
,
uniRxnDict[reactantName]={rxn};
];
,{rxn,uniRxnList}];

Return[uniRxnDict];
];


(* ::Subsubsection:: *)
(*Test function*)


(* ::Input::Initialization:: *)
fTESTrxnDictUniFromList[]:=Module[
{rxnList,rxnDict}
,
(* Reactions: A\[Rule]B *)
Print["> Making reaction A->B "];
rxnList={fReactionUnimol["r1",0.5,{"A"},{"B"}]};
Print[rxnList];

Print["> Making dictionary:"];
rxnDict=rxnDictUniFromList[rxnList];
Print[rxnDict];
];


(* ::Subsection:: *)
(*Reaction dict of bimol reactions from list*)


(* ::Input::Initialization:: *)
rxnDictBiFromList[biRxnList_]:=Module[
{biRxnDict,reactantName}
,
(* Go over all reactions *)
biRxnDict=<||>;
Do[
reactantName=Sort[{rxn["r"][[1]],rxn["r"][[2]]}];
If[
KeyExistsQ[biRxnDict,reactantName]
,
AppendTo[biRxnDict[reactantName],rxn]
,
biRxnDict[reactantName]={rxn};
];
,{rxn,biRxnList}];

Return[biRxnDict];
];


(* ::Subsubsection:: *)
(*Test function*)


(* ::Input::Initialization:: *)
fTESTrxnDictBiFromList[]:=Module[
{rxnList,rxnDict}
,
(* Reactions: A+B\[Rule]A *)
Print["> Making reaction A+B->A"];
rxnList={fReactionBimol["r1",0.5,{"A","B"},{"A"}]};
Print[rxnList];

(* Dict *)
Print["> Making dict:"];
rxnDict=rxnDictBiFromList[rxnList];
Print[rxnDict];
];


(* ::Subsection:: *)
(*Predict the diffusions, bimol reactions*)


(* ::Input::Initialization:: *)
predictDiff[lattice_,biRxnDict_,verbose_:False]:=Module[
{bimolFlag,latticeOld,latticeNew,molPosList,molList,speciesOld,molListUse,posOld,posNew,speciesCollide,rxnName,rxnsPossible,rxnsPossibleSuccess,rxn,pos1,pos2,acceptMolMove,pStr}
,
(* Flag: did any bimol reactions occur? *)
bimolFlag=False;

(* Copy over old lattice *)
latticeOld=lattice;

(* New empty lattice *)
latticeNew=ConstantArray[0,Length[latticeOld]];

(* Generate the list of molecules and positions from the lattice *)
molPosList=Select[Range[1,Length[latticeOld]],latticeOld[[#]]=!=0&];
molList=Table[{ipos,latticeOld[[ipos]]},{ipos,molPosList}];

(* Generate a permutation of molList *)
molListUse=RandomSample[molList];

(* Go through all mols *)
While[
Length[molListUse]!=0
,
(* The mol to use *)
{posOld,speciesOld}=molListUse[[1]];

(* The new position *)
posNew=Min[Length[latticeOld],Max[1,posOld+RandomChoice[{-1,1}]]];

If[verbose,
Print["Moving mol at pos = ",posOld," to possible new pos: ",posNew];
];

(* Check the lattice - is this position occupied? *)
If[latticeOld[[posNew]]=!=0,
speciesCollide=latticeOld[[posNew]],
If[latticeNew[[posNew]]=!=0,
speciesCollide=latticeNew[[posNew]],
speciesCollide=0
];];

(* Is this move ok? For now, say yes *)
acceptMolMove=True;

(* Did we collide with a mol? *)
If[
speciesCollide=!=0
,
(* Yes, collided *)
If[verbose,
Print["Collided with a mol"];
];

(* For now, this is not acceptable, unless a reaction makes it ok *)
acceptMolMove=False;

(* Yes - check if a reaction is possible *)
rxnName=Sort[{speciesOld,speciesCollide}];
rxnsPossible={};
If[KeyExistsQ[biRxnDict,rxnName],
rxnsPossible=biRxnDict[rxnName]];

If[verbose,
Print["Possible reactions: ",Length[rxnsPossible]];
];

(* Was there a reaction? *)
If[rxnsPossible!={},
(* Check all reactions *)
rxnsPossibleSuccess={};
Do[
If[RandomReal[{0,1}]<=rxn["prob"],
AppendTo[rxnsPossibleSuccess,rxn];
];
,{rxn,rxnsPossible}];

If[verbose,
Print["Possible successful reactions: ",Length[rxnsPossibleSuccess]];
];

(* Did any reaction succeed? *)
If[rxnsPossibleSuccess!={},
(* Choose a random reaction to do *)
rxn=RandomChoice[rxnsPossibleSuccess];

If[verbose,
pStr="";
If[Length[rxn["p"]]==0,pStr="0"];
If[Length[rxn["p"]]==1,pStr=rxn["p"][[1]]];
If[Length[rxn["p"]]==2,pStr=rxn["p"][[1]]<>" + "<>rxn["p"][[2]]];
Print["Doing reaction: ",rxn["r"][[1]]," + ",rxn["r"][[2]]," -> ",pStr];
];

(* Remove the mol from the old lattice *)
latticeOld[[posOld]]=0;
If[verbose,
Print["Removed old mol at pos: ",posOld];
];
(* Remove the mol we collided with *)
latticeOld[[posNew]]=0;
If[verbose,
Print["Removed mol we collided with from old lattice at pos: ",posNew];
];
latticeNew[[posNew]]=0;
If[verbose,
Print["Removed mol we collided with from new lattice at pos: ",posNew];
];
(* If the collision mol is in the mol list, remove that too *)
Do[
If[molListUse[[ic]][[1]]==posNew,molListUse=Delete[molListUse,ic];Break[];];
,{ic,1,Length[molListUse]}];

(* Figure out where to place the mols *)
If[Length[rxn["p"]]==1,
(* If there is one, just put at the new pos = empty *)
latticeNew[[posNew]]=rxn["p"][[1]];
];
If[Length[rxn["p"]]==2,
(* If there are two, just put at the old pos = also empty *)
{pos1,pos2}=RandomChoice[{{posNew,posOld},{posOld,posNew}}];
latticeNew[[pos1]]=rxn["p"][[1]];
latticeNew[[pos2]]=rxn["p"][[2]];
];
(* Forbid the third case... No more than two products accepted *)

(* These molecule moves are accepted *)
acceptMolMove=True;

(* A bimol reaction has occurred *)
bimolFlag=True;
];
];
,
(* No collision with another molecule - just commit the move *)
latticeNew[[posNew]]=speciesOld;
(* Remove the mol from the old lattice *)
latticeOld[[posOld]]=0;
];

(* If this move was rejected, just keep the mol at the same place *)
If[!acceptMolMove,
(* Keep mol at the same place *)
latticeNew[[posOld]]=speciesOld;
(* Remove the mol from the old lattice, so it is not counted in future collisions *)
latticeOld[[posOld]]=0;
];

If[verbose,
Print["Finished mol at pos = ",posOld];
];

(* Next mol by deleting this one, if there are any left *)
If[Length[molListUse]!=0,
molListUse=Delete[molListUse,1];
];
];

Return[{latticeNew,bimolFlag}]
];


(* ::Subsubsection:: *)
(*Test function*)


(* ::Input::Initialization:: *)
fTESTpredictDiff[]:=Module[
{lattice,rxnList,rxnDict,TMP}
,

(* Make an empty lattice *)
Print["> Making empty lattice and populating:"];
lattice=ConstantArray[0,10];

(* Populate the lattice *)
lattice[[2]]="A";
lattice[[3]]="A";
Print[lattice];

(* Reactions: A+A \[Rule] A, A+A \[Rule] 0 *)
Print["> Making reactions A+A->A and A+A->0:"];
rxnList={
fReactionBimol["r1",0.8,{"A","A"},{"A"}],
fReactionBimol["r2",0.8,{"A","A"},{}]
};
Print[rxnList];

(* Dict *)
Print["> Making dict:"];
rxnDict=rxnDictBiFromList[rxnList];
Print[rxnDict];

(* Diffusion *)
Print["> Taking a diffusion step and doing bimol reactions:"];
{lattice,TMP}=predictDiff[lattice,rxnDict];
Print[lattice];
Print["> Has a bimol reaction occurred?"];
Print[TMP];
];


(* ::Subsection:: *)
(*Schedule unimolecular reactions*)


(* ::Input::Initialization:: *)
scheduleUnimolRxns[lattice_,uniRxnDict_,verbose_:False]:=Module[
{props,rxns,propsCum,ctr,rxnChoose,dt}
,
(* Propensities *)
props={};
rxns={};
propsCum=0.0;

(* Go through all possible reagants *)
Do[

(* Count how many mols there currently are of this species *)
ctr=0;
Do[
If[mol===rName,
ctr+=1;
];
,{mol,lattice}];

If[verbose,
Print["scheduleUnimolRxns: # mols for: ",rName," is ",ctr];
];

(* Propensities *)
Do[
propsCum+=ctr*rxn["kr"];
AppendTo[props,ctr*rxn["kr"]];
AppendTo[rxns,rxn];
,{rxn,uniRxnDict[rName]}];

,{rName,Keys[uniRxnDict]}];

If[verbose,
Print["scheduleUnimolRxns: propensities: ",props];
Print["scheduleUnimolRxns: rxns: ",rxns];
];

If[
Total[props]>0.0
,
(* Choose a reaction *)
rxnChoose=RandomChoice[props->rxns];

(* Time of next reaction *)
dt=Log[1.0/RandomReal[{0,1}]]/propsCum;

Return[{dt,rxnChoose}];
,
Return[{Infinity,Infinity}];
];
];


(* ::Subsubsection:: *)
(*Test function*)


(* ::Input::Initialization:: *)
fTESTscheduleUnimolRxns[]:=Module[
{lattice,rxnList,rxnDict,res}
,
(* Make an empty lattice *)
Print["> Making and populating a lattice:"];
lattice=ConstantArray[0,10];
(* Populate the lattice *)
lattice[[2]]="A";
lattice[[3]]="A";

(* Reactions: A \[Rule] A, A \[Rule] 0 *)
Print["> Making reactions:"];
rxnList={
fReactionUnimol["r1",1.0,{"A"},{"A"}],
fReactionUnimol["r2",1.0,{"A"},{}]
};
Print[rxnList];

(* Dict *)
Print["> Making dict:"];
rxnDict=rxnDictUniFromList[rxnList];
Print[rxnDict];

(* Rxn *)
Print["> Scheduling unimolecular reactions:"];res=scheduleUnimolRxns[lattice,rxnDict];
Print[res];
];


(* ::Subsection:: *)
(*Do a uni reaction*)


(* ::Input::Initialization:: *)
doUniRxn[rxn_,lattice_,verbose_:False]:=Module[
{latticeNew,molPossible,reactionDone,pos,species,ipos1,ipos2,il,ir,i1,i2},
latticeNew=lattice;

(* Go over all possible reagants *)
molPossible={};
Do[
If[lattice[[ipos]]=!=0,
If[lattice[[ipos]]==rxn["r"][[1]],AppendTo[molPossible,{ipos,lattice[[ipos]]}];
];];
,{ipos,1,Length[lattice]}];

reactionDone=False;
While[Length[molPossible]!=0&&!reactionDone,

(* Pick one at random *)
{pos,species}=RandomChoice[molPossible];
If[verbose,
Print["Chose one from # possible mols: ",Length[molPossible]];
Print["No products = ",Length[rxn["p"]]];
];

(* Do the reaction *)

(* Zero products *)
If[
Length[rxn["p"]]==0,
(* Remove the old *)
latticeNew[[pos]]=0;

(* Finished *)
reactionDone=True;
If[verbose,
Print["Result: ",pos," -> no new products"];
];
];

(* One product *)
If[
Length[rxn["p"]]==1,
(* Remove the old *)
latticeNew[[pos]]=0;

(* Place the new *)
latticeNew[[pos]]=rxn["p"][[1]];

(* Finished *)
reactionDone=True;
If[verbose,
Print["Result: ",pos," -> one new product"];
];
];

(* Two products *)
If[
Length[rxn["p"]]==2,

(* First check there is room *)
ipos1=pos;
ipos2=-1;
il=Max[1,ipos1-1];
ir=Min[Length[lattice],ipos1+1];
{i1,i2}=RandomChoice[{{il,ir},{ir,il}}];
If[lattice[[i1]]===0,
ipos2=i1,
If[lattice[[i2]]===0,
ipos2=i2
];
];

(* There is room? *)
If[
ipos2!=-1
,
(* Yes room - do the reaction *)
(* Remove the old *)
latticeNew[[pos]]=0;

(* Do it *)
latticeNew[[ipos1]]=rxn["p"][[1]];
latticeNew[[ipos2]]=rxn["p"][[2]];

(* Finished *)
reactionDone=True;
If[verbose,
Print["Result: ",pos," -> two new at: ",ipos1," ",ipos2];
];
,
(* No room - make sure this mol cant be chosen again *)
molPossible=Delete[molPossible,Position[molPossible,{pos,species}]];
];
];
];

Return[latticeNew];
];


(* ::Subsubsection:: *)
(*Test function*)


(* ::Input::Initialization:: *)
fTESTdoUniRxn[]:=Module[
{lattice,rxnList,rxnDict,dtrxn,rxn}
,
(* Make an empty lattice *)
Print["> Making and populating a lattice:"];
lattice=ConstantArray[0,10];
(* Populate the lattice *)
lattice[[2]]="A";
lattice[[3]]="A";

(* Reactions: A \[Rule] 0, A \[Rule] A+A *)
Print["> Making reactions:"];
rxnList={
fReactionUnimol["r1",1.0,{"A"},{"A","A"}],
fReactionUnimol["r2",1.0,{"A"},{}]
};
Print[rxnList];

(* Dict *)
Print["> Making dict:"];
rxnDict=rxnDictUniFromList[rxnList];
Print[rxnDict];

(* Schedule *)
Print["> Scheduling unimolecular reactions:"];{dtrxn,rxn}=scheduleUnimolRxns[lattice,rxnDict];
Print[{dtrxn,rxn}];

(* Do the reaction *)
Print["> Doing the reactions:"];
lattice=doUniRxn[rxn,lattice];
Print[lattice];
];


(* ::Subsection:: *)
(*Lattice spacing from dt, dc*)


(* ::Text:: *)
(*D = cm^2 / s*)
(*dt = s*)
(*step size = a = lattice grid spacing*)
(*a = sqrt(2*D*dt)*)
(*Unimolecular reaction rates in 1/s*)


(* ::Input::Initialization:: *)
latticeSpacing[dc_,dt_]:=Sqrt[2*dc*dt];


(* ::Subsection:: *)
(*Functions to get moments of a lattice - note they are symmetrized*)


(* ::Input::Initialization:: *)
fMean[lattice_,s_]:=Sum[If[l===s,1,0],{l,lattice}];


(* ::Input::Initialization:: *)
fNN[lattice_,s1_,s2_]:=Sum[
If[(lattice[[il]]===s1&&lattice[[il+1]]===s2)||(lattice[[il]]===s2&&lattice[[il+1]]===s1),
1,
0]
,{il,1,Length[lattice]-1}];


(* ::Input::Initialization:: *)
fNNN[lattice_,s1_,s2_]:=Sum[
If[(lattice[[il]]===s1&&lattice[[il+2]]===s2)||(lattice[[il]]===s2&&lattice[[il+2]]===s1),
1,0]
,{il,1,Length[lattice]-2}];


(* ::Input::Initialization:: *)
fTriplet[lattice_,s1_,s2_,s3_]:=Sum[
If[(lattice[[il]]===s1&&lattice[[il+1]]===s2&&lattice[[il+2]]===s3)||(lattice[[il]]===s1&&lattice[[il+1]]===s3&&lattice[[il+2]]===s2)||(lattice[[il]]===s2&&lattice[[il+1]]===s1&&lattice[[il+2]]===s3)||(lattice[[il]]===s2&&lattice[[il+1]]===s3&lattice[[il+2]]===s1)||
(lattice[[il]]===s3&&lattice[[il+1]]===s1&&lattice[[il+2]]===s2)||(lattice[[il]]===s3&&lattice[[il+1]]===s2&&lattice[[il+2]]===s1),
1,0]
,{il,1,Length[lattice]-2}];


(* ::Input::Initialization:: *)
fMean[lattice_]:=Total[Table[If[l=!=0,1,0],{l,lattice}]];


(* ::Input::Initialization:: *)
fNN[lattice_]:=Total[Table[If[lattice[[il]]=!=0&&lattice[[il+1]]=!=0,1,0],{il,1,Length[lattice]-1}]];


(* ::Input::Initialization:: *)
fNNN[lattice_]:=Total[Table[If[lattice[[il]]=!=0&&lattice[[il+2]]=!=0,1,0],{il,1,Length[lattice]-2}]];


(* ::Input::Initialization:: *)
fTriplet[lattice_]:=Total[Table[If[lattice[[il]]=!=0&&lattice[[il+1]]=!=0&&lattice[[il+2]]=!=0,1,0],{il,1,Length[lattice]-2}]];


(* ::Subsubsection:: *)
(*Test function*)


(* ::Input::Initialization:: *)
fTESTmoments[]:=Module[
{lattice}
,
(* Make an empty lattice *)
Print["> Making and populating a lattice:"];
lattice=ConstantArray[0,10];
(* Populate the lattice *)
lattice[[2]]="A";
lattice[[3]]="A";
lattice[[4]]="B";
lattice[[5]]="B";
Print[lattice];

Print["> Mean no particles:"];
Print[fMean[lattice]];
Print["> Mean no A:"];
Print[fMean[lattice,"A"]];
Print["> Mean no B:"];
Print[fMean[lattice,"B"]];

Print["> NN particles:"];
Print[fNN[lattice]];
Print["> NN AA:"];
Print[fNN[lattice,"A","A"]];
Print["> NN BB:"];
Print[fNN[lattice,"B","B"]];
Print["> NN AB + BA:"];
Print[fNN[lattice,"A","B"]];
];


(* ::Subsection:: *)
(*Function to generate an initial condition by gibbs sampling*)


(* ::Text:: *)
(*eg. http://www.cs.utoronto.ca/~yueli/CSC321_UTM_2014_files/tut9.pdf*)


(* ::Input::Initialization:: *)
fGenLattice[hDict_,jDict_,nChain_,nIter_,trackMoments_:False]:=Module[
{spNames,spNamesZ,spNNs,lattice,posAvailable,posSt,pos0,latticeNew,posFlip,oldVal,newVal,hOld,jOld,hNew,jNew,eDiff,posLeft,posRight,meanSt,nnSt}
,
(* Species names *)
spNames=Keys[hDict];
spNamesZ=Join[{0},spNames];

(* NN combinations *)
If[trackMoments,
spNNs={};
Do[
Do[
AppendTo[spNNs,{spNames[[isp1]],spNames[[isp2]]}];
,{isp2,isp1,Length[spNames]}];
,{isp1,1,Length[spNames]}];
];

(* Generate an empty lattice *)
lattice=ConstantArray[0,nChain];

(* Positions available *)
posAvailable=Range[1,nChain];

(* Store positions of each species *)
posSt=Association[];

(* Generate a random initial state *)
pos0=RandomSample[posAvailable,RandomInteger[{0,nChain}]];
Do[
lattice[[p0]]=RandomChoice[spNames];
,{p0,pos0}];

(* Store the mean, nn of the lattice *)
If[trackMoments,
meanSt=Association[];
nnSt=Association[];
Do[
meanSt[sp]={fMean[lattice,sp]};
,{sp,spNames}];
Do[
nnSt[sps]={fNN[lattice,sps[[1]],sps[[2]]]};
,{sps,spNNs}];
];

(* Do the sampling *)
Do[
(* Calculate a new test state by changing one of the spins *)
latticeNew=lattice;
posFlip=RandomInteger[{1,nChain}];
oldVal=lattice[[posFlip]];
newVal=RandomChoice[spNamesZ];
While[newVal===oldVal,newVal=RandomChoice[spNamesZ];];
latticeNew[[posFlip]]=newVal;

(* Calculate the energy difference: new - old *)
If[oldVal===0,
hOld=0;
jOld=0;
,
hOld=-hDict[oldVal];
jOld=0;
(* Coupling on the left *)
posLeft=posFlip-1;
If[posLeft>=1,
If[lattice[[posLeft]]=!=0,
jOld-=jDict[Sort[{oldVal,lattice[[posLeft]]}]];
];];
(* Coupling on the right *)
posRight=posFlip+1;
If[posRight<=nChain,
If[lattice[[posRight]]=!=0,
jOld-=jDict[Sort[{oldVal,lattice[[posRight]]}]];
];];
];
If[newVal===0,
hNew=0;
jNew=0;
,
hNew=-hDict[newVal];
jNew=0;
(* Coupling on the left *)
posLeft=posFlip-1;
If[posLeft>=1,
If[lattice[[posLeft]]=!=0,
jNew-=jDict[Sort[{newVal,lattice[[posLeft]]}]];
];];
(* Coupling on the right *)
posRight=posFlip+1;
If[posRight<=nChain,
If[lattice[[posRight]]=!=0,
jNew-=jDict[Sort[{newVal,lattice[[posRight]]}]];
];];
];
eDiff=hNew+jNew-hOld-jOld;

(* Accept spin flip? *)
If[eDiff<0||Exp[-eDiff]>RandomReal[{0,1}],
lattice=latticeNew;
];

(* Store new mean and nn *)
If[trackMoments,
Do[
AppendTo[meanSt[sp],fMean[lattice,sp]];
,{sp,spNames}];
Do[
AppendTo[nnSt[sps],fNN[lattice,sps[[1]],sps[[2]]]];
,{sps,spNNs}];
];

,{iIter,1,nIter}];

If[trackMoments,
Return[{lattice,meanSt,nnSt}];
,
Return[lattice];
];
];


(* ::Input::Initialization:: *)
fMomentsFromhj[h_,j_,nChain_]:={
(E^-j (2+E^j (-1+E^(h+j) (1+Sqrt[E^(-2 (h+j)) (1+E^h (4+E^j (-2+E^(h+j))))]))) nChain)/(Sqrt[E^(-2 (h+j)) (1+E^h (4+E^j (-2+E^(h+j))))] (1+E^(h+j) (1+Sqrt[E^(-2 (h+j)) (1+E^h (4+E^j (-2+E^(h+j))))])))
,
((-1+E^(h+j) (1+Sqrt[E^(-2 (h+j)) (1+E^h (4+E^j (-2+E^(h+j))))])) nChain)/(Sqrt[E^(-2 (h+j)) (1+E^h (4+E^j (-2+E^(h+j))))] (1+E^(h+j) (1+Sqrt[E^(-2 (h+j)) (1+E^h (4+E^j (-2+E^(h+j))))])))
}//N;


(* ::Subsubsection:: *)
(*Test function*)


(* ::Input::Initialization:: *)
fTEST1DGenLattice[]:=Module[
{}
,
Print["> Test parameters:"];
hTest=0.1;
jTest=0.3;
nTest=100;
Print["h=",hTest," j=",jTest," n=",nTest];

Print["> Correct moments:"];
{m1,m2}=fMomentsFromhj[hTest,jTest,nTest];
Print["Mean: ",m1," NN: ",m2];

Print["> Making dicts:"];
hDict=<|"A"->hTest|>;
jDict=<|{"A","A"}->jTest|>;
Print[hDict];
Print[jDict];

Print["> Generating:"];
{lattice,meanSt,nnSt}=fGenLattice[hDict,jDict,nTest,1000,True];

Return[{
Show[
ListLinePlot[meanSt["A"],PlotLabel->"Mean"]
,
ListLinePlot[{{0,m1},{1000,m1}},PlotStyle->Gray]
],
Show[
ListLinePlot[nnSt[{"A","A"}],PlotLabel->"NN"]
,
ListLinePlot[{{0,m2},{1000,m2}},PlotStyle->Gray]
]
}];
];


(* ::Subsubsection:: *)
(*Test - 2D*)


(* ::Input::Initialization:: *)
fTEST2DGenLattice[]:=Module[
{}
,
Print["> Test parameters:"];
hATest=0.1;
hBTest=0.2;
jAATest=0.3;
jBBTest=0;
jABTest=-0.1;
nTest=100;
Print["hA=",hATest,"hB=",hBTest," jAA=",jAATest," jBB=",jBBTest," jAB=",jABTest," n=",nTest];

Print["> Making dicts:"];
hDict=<|"A"->0.1,"B"->0.2|>;
jDict=<|{"A","A"}->0.3,{"B","B"}->0,{"A","B"}->-0.1|>;
Print[hDict];
Print[jDict];

Print["> Generating:"];
{lattice,meanSt,nnSt}=fGenLattice[hDict,jDict,nTest,1000,True];

Return[{
ListLinePlot[meanSt["A"],PlotLabel->"Mean A"]
,
ListLinePlot[meanSt["B"],PlotLabel->"Mean B"]
,
ListLinePlot[nnSt[{"A","A"}],PlotLabel->"NN A,A"]
,
ListLinePlot[nnSt[{"B","B"}],PlotLabel->"NN B,B"]
,
ListLinePlot[nnSt[{"A","B"}],PlotLabel->"NN A,B"]
}];
]


(* ::Subsection:: *)
(*Main Gillespie function*)


(* ::Input::Initialization:: *)
fGillespie[tMax_,dt_,lattice0_,uniRxnList_,biRxnList_,verbose_:False]:=Module[
{biRxnDict,uniRxnDict,lattice,dtNextURxn,nextURxn,checkUniRxns,tNextURxn,t,bimolFlag,
latticeSt}
,

(* Store the lattice *)
latticeSt=Association[];

(* Generate a reaction dict form the reaction lists *)
biRxnDict=rxnDictBiFromList[biRxnList];
uniRxnDict=rxnDictUniFromList[uniRxnList];

(* Initial condition *)
lattice=lattice0;
latticeSt[0]=lattice;

(* Are there any unimolecular reactions? *)
If[Length[uniRxnList]!=0,
checkUniRxns=True,
checkUniRxns=False
];

(* Schedule the next unimolecular reaction *)
If[checkUniRxns,
{dtNextURxn,nextURxn}=scheduleUnimolRxns[lattice,uniRxnDict];
If[dtNextURxn!=Infinity,tNextURxn=dtNextURxn,tNextURxn=Infinity];
];

If[verbose,
Print["First scheduled uni rxn: ",tNextURxn," ",nextURxn];
];

(* Iterate over all time *)
t=0.0;
While[t<tMax,
If[verbose,
Print["----------"];
Print["Time = ",t];
];

(* Do all unimolecular reactions *)
If[checkUniRxns,
While[tNextURxn<=t+dt,
If[verbose,
Print["Doing unimol rxn at time = ",tNextURxn];
];

(* Do the reaction *)
lattice=doUniRxn[nextURxn,lattice];

(* Next unimolecular reaction *)
{dtNextURxn,nextURxn}=scheduleUnimolRxns[lattice,uniRxnDict];
If[dtNextURxn!=Infinity,tNextURxn+=dtNextURxn,tNextURxn=Infinity];
If[verbose,
Print["Next uni rxn: ",tNextURxn," ",nextURxn]
];
];

If[verbose,
Print["Finished uni rxns"];
];
];

(* Move the mols, doing bimol reactions *)
{lattice,bimolFlag}=predictDiff[lattice,biRxnDict,verbose];
If[verbose,
Print["Did diffusion step. Did any bimol reactions happen? ",bimolFlag];
];

If[checkUniRxns&&bimolFlag,
(* Reschedule the next unimol reaction once more *)
{dtNextURxn,nextURxn}=scheduleUnimolRxns[lattice,uniRxnDict];
If[dtNextURxn!=Infinity,tNextURxn=t+dt+dtNextURxn,tNextURxn=Infinity];

If[verbose,
Print["Rescheduled uni rxn: ",tNextURxn," ",nextURxn];
];
];

(* Increment time *)
t+=dt;

(* Store *)
latticeSt[t]=lattice;
];

If[verbose,
Print["----------"];
Print["Finished sim"];
];

Return[latticeSt];
];


(* ::Chapter:: *)
(*End*)


End[]


EndPackage[]
