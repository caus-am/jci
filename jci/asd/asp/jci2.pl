% Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

% JCI background knowledge, assumption 2 (randomization)
%
% assumes nodes and inodes have been defined

% no context node is confounded with a system node
:-conf(A,B), inode(A), node(B), not inode(B), A<B.
:-conf(B,A), inode(A), node(B), not inode(B), B<A.
