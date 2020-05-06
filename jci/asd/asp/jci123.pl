% Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

% JCI background knowledge
%
% assumes nodes and inodes have been defined
%
% context nodes have no incoming edges
:-inode(A), node(B), edge(B,A), A!=B.
% no context node is confounded with a system node
:-conf(A,B), inode(A), node(B), not inode(B), A<B.
:-conf(B,A), inode(A), node(B), not inode(B), B<A.
% all pairs of context nodes are confounded
conf(A,B) :- inode(A), inode(B), A<B.
