% test
clc;
global FEM
% generate 8-noded CQUAD8 for a unit square panel

% 1) the current code works for even number element, will update it for all elements
% number, sorry for inconvenience.


xmesh_size = 8; % use even value
ymesh_size = 8; % use even value

% generate mesh
FEM = mesh_QUAD8_v2(xmesh_size,ymesh_size);

% plot
patch_plot( FEM.elementNodes,FEM.nodesCord,101,'skin')