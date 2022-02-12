% read stress in stiffeners

%% 

filename = [pwd filesep 'plate_with_hole_stiffeners.bdf'];

 
previous_folder_path = fileparts(pwd);

 
[QUAD4_Elem,TRIA3_Elem,BEAM3_Elem,BEAM2_Elem,bdf_node,bdf_cord]=read_bdf_v2(filename);

bdf_elem = QUAD4_Elem;

BASE_element_connectivity = bdf_elem(:,3:6);
BASE_node_cords = bdf_node;
BASE_node_cords(:,2) = BASE_node_cords (:,2)/0.15*0.254/2 + 0.254/2; % The original side length is 0.3, made slightly changed to 0.254 later.
BASE_node_cords(:,3) = BASE_node_cords (:,3)/0.15*0.254/2 + 0.254/2;

patch_plot(BASE_element_connectivity, BASE_node_cords,201,'skin',BASE_node_cords);colorbar off;title('NASTRAN MESH');

NASTRAN.grid_cords = BASE_node_cords;
NASTRAN.elements  = BASE_element_connectivity;



%% stiffener 1 axial stress
bdf_name ='C:\Users\wzhao\Documents\GitHub\SPAD\JournalofCompositeScience\Verification_NASTRAN_models\plate_hole_4_stiffeners_shear\stiffener_1.dat';

[~,~,~,BEAM2_Elem,~,~,~]=read_bdf_v2(bdf_name);


% read stress in f06

% pch_filename = 'C:\Users\wzhao\Documents\GitHub\SPAD\JournalofCompositeScience\Verification_NASTRAN_models\plate_hole_4_stiffeners_axial\stress_punch_file\sol105_plate_axial_stiffeners.pch';

% pch_filename = 'C:\Users\wzhao\Documents\GitHub\SPAD\JournalofCompositeScience\Verification_NASTRAN_models\plate_hole_4_stiffeners_bend\stress_punch_file\sol105_stiffened_plate_cutout_bend.pch';

pch_filename ='C:\Users\wzhao\Documents\GitHub\SPAD\JournalofCompositeScience\Verification_NASTRAN_models\plate_hole_4_stiffeners_shear\stress_punch_file\sol105_stiffened_plate_cutout_shear.pch';

minElemid = min(BEAM2_Elem(:,1));
maxElemid = max(BEAM2_Elem(:,1));

Beam_stress(:,:,1) =read_cbeam_stress(pch_filename,minElemid,maxElemid);


stress_location_x = BASE_node_cords(Beam_stress(:,2,1),2);

figure; plot(stress_location_x,Beam_stress(:,4,1))

% present stiffener stress


% present_stiff_1_loc_x = Stiffener.PointsCoord(1:end-1,2);

present_stiff_loc_x = Stiffener.nodes_coordinates_all(Stiffener.elements_all(:,2,1),2);


present_stress = Stiffener.stress(:,1,1);%interp1(present_stiff_loc_x,Stiffener.stress(:,1,1),stress_location_x,'nearest')

hold on;

plot(present_stiff_loc_x ,present_stress ,'o')





%% stiffener 2



bdf_name ='C:\Users\wzhao\Documents\GitHub\SPAD\JournalofCompositeScience\Verification_NASTRAN_models\plate_hole_4_stiffeners_shear\stiffener_2.dat';

[~,~,~,BEAM2_Elem,~,~,~]=read_bdf_v2(bdf_name)


stiffno = 2;
% read stress in f06

% pch_filename = 'C:\Users\wzhao\Documents\GitHub\SPAD\JournalofCompositeScience\Verification_NASTRAN_models\plate_hole_4_stiffeners_axial\stress_punch_file\sol105_plate_axial_stiffeners.pch';



minElemid = min(BEAM2_Elem(:,1));
maxElemid = max(BEAM2_Elem(:,1));

Beam_stress(:,:,stiffno)  =read_cbeam_stress(pch_filename,minElemid,maxElemid);


stress_location_x = BASE_node_cords(Beam_stress(:,2,stiffno),3);

figure(12);hold on; plot(stress_location_x,Beam_stress(:,4,stiffno))

% present stiffener stress


% present_stiff_loc_x = Stiffener.PointsCoord(1:end-1,2);

present_stiff_loc_x = Stiffener.nodes_coordinates_all(Stiffener.elements_all(:,2,stiffno),3,3);

% present_stress = interp1(present_stiff_loc_x,Stiffener.stress(:,1,2),stress_location_x,'spline')
hold on;
plot(present_stiff_loc_x,Stiffener.stress(:,1,3),'o');









