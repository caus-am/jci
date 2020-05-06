% Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

% Automatic generation of full computation tree for d-connection graphs
% NOTE: this one uses a different convention than Antti's R code, i.e.,
% node labels start counting from 0!
node(0..nrnodes-1).
set(0..2**nrnodes-1).
ismember(M,Z) :- set(M), node(Z), M & (2**Z) != 0.
intervene(0,Jsub,Z,J,0)   :- node(Z), set(J), ismember(J,Z), set(Jsub), not ismember(Jsub,Z), (Jsub + 2**Z) == J.
marginalize(C,J,Msub,Z,M) :- node(Z), set(M), set(C), set(J), (C & M) == 0, (M & J) == 0, (J & C) == 0, not ismember(J,Z), ismember(M,Z), not ismember(C,Z), set(Msub), not ismember(Msub,Z), (Msub + 2**Z) == M.
condition(Csub,Z,C,J,M)   :- node(Z), set(M), set(C), set(J), (C & M) == 0, (M & J) == 0, (J & C) == 0, not ismember(J,Z), not ismember(M,Z), ismember(C,Z), set(Csub), not ismember(Csub,Z), (Csub + 2**Z) == C.
