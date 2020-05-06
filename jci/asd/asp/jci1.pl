% Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

% JCI background knowledge, assumption 1 (exogenous context)
%
% assumes nodes and inodes have been defined

% context nodes have no incoming edges from system variables
:-inode(A), node(B), not inode(B), edge(B,A).
