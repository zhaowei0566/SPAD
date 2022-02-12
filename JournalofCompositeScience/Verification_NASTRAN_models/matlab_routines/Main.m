%% Read buckling mode shapes from NASTRAN output file in *.pch
% The stress for each layer described in the global coordinate system is
% exported to file *.rpt. This is done in PATRAN
%
% Provided to readers for cheecking the data/figure reported in the paper
%% 1. Read NASTRAN elements and nodes
% The original BASE bdf file is given in plate_with_hole_stiffeners.bdf
%
%
filename = [pwd filesep 'plate_with_hole_stiffeners.bdf'];

 
previous_folder_path = fileparts(pwd);

 
[QUAD4_Elem,TRIA3_Elem,BEAM3_Elem,bdf_node,bdf_cord]=read_bdf_v2(filename);

bdf_elem = QUAD4_Elem;

BASE_element_connectivity = bdf_elem(:,3:6);
BASE_node_cords = bdf_node;
BASE_node_cords(:,2) = BASE_node_cords (:,2)/0.15*0.254/2 + 0.254/2; % The original side length is 0.3, made slightly changed to 0.254 later.
BASE_node_cords(:,3) = BASE_node_cords (:,3)/0.15*0.254/2 + 0.254/2;

patch_plot(BASE_element_connectivity, BASE_node_cords,201,'skin',BASE_node_cords);colorbar off;title('NASTRAN MESH');

NASTRAN.grid_cords = BASE_node_cords;
NASTRAN.elements  = BASE_element_connectivity;


%% 2. Read buckling mode shapes


pchfname  = [previous_folder_path filesep 'plate_hole_4_stiffeners_bend\sol105_stiffened_plate_cutout_bend.pch'];

minodeid  = BASE_node_cords(1,1);
maxnodeid = BASE_node_cords(end,1);

Modelimit = 10; % maximum mode number

NASTRAN_ModeShape = readmodeshape_v2(pchfname,minodeid,maxnodeid,Modelimit);



modeNo = 7; % chosen for plotting

undeformedX = NASTRAN.grid_cords(:,2);
undeformedY = NASTRAN.grid_cords(:,3);
undeformedZ = NASTRAN.grid_cords(:,4);




deformUZ = NASTRAN_ModeShape.modeshape(:,4,modeNo);

scalefactor =  1/max(abs(deformUZ(:))); % normalize the maximum magnitude to 1.

deformedW = NASTRAN.grid_cords;
deformedW(:,4) =   deformedW(:,4) +  deformUZ* scalefactor;

figure1 = figure(1000+modeNo);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
patch_plot(NASTRAN.elements,NASTRAN.grid_cords,1000+modeNo,'modeshape',deformedW);axis image;
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'FontSize',16);
colorbar('FontSize',16)
colormap(jet(20));
hold off


%% 3. Read stress for each layer described in the global coordinate system
% global_cord_stress.rpt include stress described in the global coordinate
% system, which is obtained using PATRAN

minElmid = 1; % minimum element ID
maxElmid = size(NASTRAN.elements,1); % maximum element ID

Layerlimit = 8; % number of layers

rpt_fname = [previous_folder_path '\plate_hole_4_stiffeners_bend\global_cord_stress.rpt'];
Stress=read_stress_exported_from_patran(rpt_fname,minElmid,maxElmid,Layerlimit);

% interpolate

for ii = 1:maxElmid
    
    Resultant_Nxx(ii) =  sum(Stress(ii,2,:))*1.272e-4;
    Resultant_Nyy(ii) =  sum(Stress(ii,3,:))*1.272e-4;
    Resultant_Nxy(ii) =  sum(Stress(ii,4,:))*1.272e-4;
    
    elem_center_nodes(ii,1) = ii;
    elem_center_nodes(ii,2) = sum(NASTRAN.grid_cords(NASTRAN.elements(ii,:),2))/4;
    elem_center_nodes(ii,3) = sum(NASTRAN.grid_cords(NASTRAN.elements(ii,:),3))/4;
    elem_center_nodes(ii,4) = sum(NASTRAN.grid_cords(NASTRAN.elements(ii,:),4))/4;
    
end

% 

H_1D_3D = rbf_interface(elem_center_nodes(:,2:4),NASTRAN.grid_cords(:,2:4) ,2, 0.2,0);

Model_nodes_Nxx = H_1D_3D*Resultant_Nxx';

Model_nodes_Nyy = H_1D_3D*Resultant_Nyy';

Model_nodes_Nxy = H_1D_3D*Resultant_Nxy';

%%  Plot Nxx
deformed_stress = NASTRAN.grid_cords;
deformed_stress(:,4) = -Model_nodes_Nxx; % for plot in patch_plot

figure1 = figure(2001);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
patch_plot(NASTRAN.elements,NASTRAN.grid_cords,2001,'modeshape',deformed_stress);axis image;
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'FontSize',16);
colorbar('FontSize',16)
colormap(jet(20)); caxis auto

% Similar for Nyy

deformed_stress = NASTRAN.grid_cords;
deformed_stress(:,4) = -Model_nodes_Nyy; % for plot

figure1 = figure(2002);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
patch_plot(NASTRAN.elements,NASTRAN.grid_cords,2002,'modeshape',deformed_stress);axis image;
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'FontSize',16);
colorbar('FontSize',16)
colormap(jet(20)); caxis auto


% Similar for Nxy

deformed_stress = NASTRAN.grid_cords;
deformed_stress(:,4) = -Model_nodes_Nxy; % for plot

figure1 = figure(2003);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
patch_plot(NASTRAN.elements,NASTRAN.grid_cords,2003,'modeshape',deformed_stress);axis image;
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'FontSize',16);
colorbar('FontSize',16)
colormap(jet(20)); caxis auto

%% End of Code

