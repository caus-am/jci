% Copyright (c) 2014-2020, Antti Hyttinen, Frederick Eberhardt, Matti Järvisalo, Joris M. Mooij
% All rights reserved.
%
% Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
%
% ---
%
%%% Modified version of new_wmaxsat_acyclic.pl in the code of HEJ 2014
%%%
%%% - added ancestor predicate
%%% - removed before predicate
%%% - added comments

% guess (generate all possible directed mixed graphs)
{ edge(X,Y) } :- node(X), node(Y), X != Y.
{ conf(X,Y) } :- node(X), node(Y), X < Y.

% source node in computation graph
th(X,Y,0,0,0) :- edge(X,Y).
hh(X,Y,0,0,0) :- conf(X,Y).

% ancestral relations after perfect (surgical) intervention on J
ancestor(X,Y,J) :- th(X,Y,0,J,0), X!=Y, node(X), node(Y), set(J).
ancestor(X,Y,J) :- ancestor(X,Z,J), ancestor(Z,Y,J), 
                   X!=Y, X!=Z, Y!=Z, 
                   node(X), node(Y), node(Z), set(J).
% the following lines allow one to use (non)ancestral relations as input
:- not th(X,Y,0,J,M), ancestor(X,Y,J), node(X), node(Y), set(J), X!=Y, set(M), M + 2**(X) + 2**(Y) == 2**(nrnodes)-1.
:- th(X,Y,0,J,M), not ancestor(X,Y,J), node(X), node(Y), set(J), X!=Y, set(M), M + 2**(X) + 2**(Y) == 2**(nrnodes)-1.

% acyclicity
:- ancestor(X,Y,0), ancestor(Y,X,0), node(X), node(Y), X!=Y.

%%%%%%%%%% INTERVENTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt(X,Y,C,J,M) :- tt(X,Y,C,Jsub,M), 
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 intervene(C,Jsub,Z,J,M).

th(X,Y,C,J,M) :- th(X,Y,C,Jsub,M), 
                 X != Y,
                 not ismember(J,Y),
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 intervene(C,Jsub,Z,J,M).

hh(X,Y,C,J,M) :- hh(X,Y,C,Jsub,M), 
                 X < Y,
                 not ismember(J,Y), not ismember(J,X),
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 intervene(C,Jsub,Z,J,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% CONDITIONING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% X---Y => X---Y
tt(X,Y,C,J,M) :- tt(X,Y,Csub,J,M), 
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 condition(Csub,Z,C,J,M).

%% X-->Z<--Y => X---Y
tt(X,Y,C,J,M) :- th(X,Z,Csub,J,M),th(Y,Z,Csub,J,M), 
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 condition(Csub,Z,C,J,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% X-->Y => X-->Y
th(X,Y,C,J,M) :- th(X,Y,Csub,J,M), 
                 X != Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 condition(Csub,Z,C,J,M).

%% X-->Z<->Y => X-->Y
th(X,Y,C,J,M) :- th(X,Z,Csub,J,M),
                 { hh(Z,Y,Csub,J,M); hh(Y,Z,Csub,J,M) } >= 1,
                 X != Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 condition(Csub,Z,C,J,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% X<->Y => X<->Y
hh(X,Y,C,J,M) :- hh(X,Y,Csub,J,M), 
                 X < Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 condition(Csub,Z,C,J,M).

%% X<->Z<->Y => X<->Y
hh(X,Y,C,J,M) :- { hh(Z,X,Csub,J,M); hh(X,Z,Csub,J,M) } >= 1,
                 { hh(Z,Y,Csub,J,M); hh(Y,Z,Csub,J,M) } >= 1,
                 X < Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 condition(Csub,Z,C,J,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% MARGINALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% X---Y => X---Y
tt(X,Y,C,J,M) :- tt(X,Y,C,J,Msub), 
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X-->Z---Y => X---Y
tt(X,Y,C,J,M) :- th(X,Z,C,J,Msub), 
                 { tt(Z,Y,C,J,Msub); tt(Y,Z,C,J,Msub) } >= 1,
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X---Z<--Y => X---Y
tt(X,Y,C,J,M) :- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1, 
                 th(Y,Z,C,J,Msub),
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X---Z---Y => X---Y
tt(X,Y,C,J,M) :- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1, 
                 { tt(Z,Y,C,J,Msub); tt(Y,Z,C,J,Msub) } >= 1,
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X-->Z---Z<--Y => X---Y
tt(X,Y,C,J,M) :- th(X,Z,C,J,Msub),th(Y,Z,C,J,Msub),tt(Z,Z,C,J,Msub), 
                 X <= Y, 
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% X-->Y => X-->Y
th(X,Y,C,J,M) :- th(X,Y,C,J,Msub), 
                 X != Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X-->Z-->Y => X-->Y
th(X,Y,C,J,M) :- th(X,Z,C,J,Msub),th(Z,Y,C,J,Msub), 
                 X != Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X---Z<->Y => X-->Y
th(X,Y,C,J,M) :- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1,
                 { hh(Z,Y,C,J,Msub); hh(Y,Z,C,J,Msub) } >= 1, 
                 X != Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X---Z-->Y => X-->Y
th(X,Y,C,J,M) :- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1,
                 th(Z,Y,C,J,Msub), 
                 X != Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X-->Z---Z<->Y => X-->Y
th(X,Y,C,J,M) :- th(X,Z,C,J,Msub), 
                 tt(Z,Z,C,J,Msub),
                 { hh(Z,Y,C,J,Msub); hh(Y,Z,C,J,Msub) } >= 1, 
                 X != Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% X<->Y => X<->Y
hh(X,Y,C,J,M) :- hh(X,Y,C,J,Msub),
                 X < Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X<->Z-->Y => X<->Y
hh(X,Y,C,J,M) :- { hh(X,Z,C,J,Msub); hh(Z,X,C,J,Msub) } >= 1,
                 th(Z,Y,C,J,Msub),
                 X < Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X<--Z-->Y => X<->Y
hh(X,Y,C,J,M) :- th(Z,X,C,J,Msub),
                 th(Z,Y,C,J,Msub),
                 X < Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),                                
                 marginalize(C,J,Msub,Z,M).

%% X<--Z<->Y => X<->Y
hh(X,Y,C,J,M) :- th(Z,X,C,J,Msub),
                 { hh(Y,Z,C,J,Msub); hh(Z,Y,C,J,Msub) } >= 1,
                 X < Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).

%% X<->Z---Z<->Y => X<->Y
hh(X,Y,C,J,M) :- { hh(X,Z,C,J,Msub); hh(Z,X,C,J,Msub) } >= 1,
                 { hh(Y,Z,C,J,Msub); hh(Z,Y,C,J,Msub) } >= 1,
                 tt(Z,Z,C,J,Msub),
                 X < Y,
                 not ismember(C,X), not ismember(C,Y),
                 not ismember(M,X), not ismember(M,Y),
                 node(X),node(Y),
                 marginalize(C,J,Msub,Z,M).
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% LOSS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fail(X,Y,C,J,M,W) :- th(X,Y,C,J,M), indep(X,Y,C,J,M,W), X<Y.
fail(X,Y,C,J,M,W) :- th(Y,X,C,J,M), indep(X,Y,C,J,M,W), X<Y.
fail(X,Y,C,J,M,W) :- hh(X,Y,C,J,M), indep(X,Y,C,J,M,W), X<Y.
fail(X,Y,C,J,M,W) :- tt(X,Y,C,J,M), indep(X,Y,C,J,M,W), X<Y.

fail(X,Y,C,J,M,W) :- not th(X,Y,C,J,M), 
                     not th(Y,X,C,J,M), 
                     not hh(X,Y,C,J,M), 
                     not tt(X,Y,C,J,M), 
                     dep(X,Y,C,J,M,W), X<Y.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% OPTIMIZATION PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:~fail(X,Y,C,J,M,W). [W,X,Y,C,J,M]

%with penalty on edges:
%#:~fail(X,Y,C,J,M,W). [W,0,X,Y,C,J,M]
%:~edge(X,Y). [edgepenalty,1,X,Y,0,0,0]
%:~conf(X,Y). [confpenalty,2,X,Y,0,0,0]

%#minimize{W,X,Y,C,J,M:fail(X,Y,C,J,M,W) }.
