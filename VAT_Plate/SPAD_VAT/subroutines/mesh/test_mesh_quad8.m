% test
clc;
global FEM
% generate 8-noded CQUAD8 for a unit square panel

xmesh_size = 12;
ymesh_size = 12;

% generate mesh
FEM = mesh_QUAD8_v2(xmesh_size,ymesh_size);

% plot
patch_plot( FEM.elementNodes,FEM.nodesCord,101,'skin')