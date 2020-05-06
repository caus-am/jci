% Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
%
% ---
%
% Automatic generation of full computation tree for d-connection graphs
%
% NOTE: this one uses a different convention than in the HEJ2014 code:
% node labels start counting from 0!
%
% first intervene, then marginalize, then condition;
%   always remove the node with the highest possible label first
%
node(0..nrnodes-1).
set(0..2**nrnodes-1).
ismember(M,Z) :- set(M), node(Z), M & (2**Z) != 0.
intervene(0,Jsub,Z,J,0)   :- node(Z), set(J), ismember(J,Z), set(Jsub), not ismember(Jsub,Z), (Jsub + 2**Z) == J, 2**Z > Jsub.
marginalize(0,J,Msub,Z,M) :- node(Z), set(M), set(J), ismember(M,Z), set(Msub), not ismember(Msub,Z), (Msub + 2**Z) == M, 2**Z > Msub.
condition(Csub,Z,C,J,M)   :- node(Z), set(M), set(C), set(J), (C & M) == 0, ismember(C,Z), set(Csub), not ismember(Csub,Z), (Csub + 2**Z) == C, 2**Z > Csub.
