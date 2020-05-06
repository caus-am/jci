% Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

% JCI background knowledge, assumption 3 (purely confounded context nodes)
%
% assumes nodes and inodes have been defined

% all pairs of context nodes are confounded
conf(A,B) :- inode(A), inode(B), A<B.
% context nodes have no directed edges between them
:-inode(A), inode(B), edge(B,A), A!=B.
