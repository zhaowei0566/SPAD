%% Buckling anlaysis of varying angle tow laminated plate under inplane shortning
% Step 1: Linear static analysis
% Step 2: Buckling analysis

clear all;warning off;format long;
close all;

addpath(genpath([pwd filesep 'subroutines']));
%

%% Plate dimension
Plate.length = .254;
Plate.width  = .254;


%% Finite element mesh generation

FEM = mesh_QUAD8_v2(32,32); % mesh in natural space

FEM.numberElements = size(FEM.elementNodes ,1);

% Transform mesh in natural space to physical space
FEM.nodesCord(:,2) = FEM.nodesCord(:,2)/max(FEM.nodesCord(:,2))*Plate.length ;
FEM.nodesCord(:,3) = FEM.nodesCord(:,3)/max(FEM.nodesCord(:,3))*Plate.width ;

Xcoord = FEM.nodesCord(:,2);
Ycoord = FEM.nodesCord(:,3);
Zcoord = FEM.nodesCord(:,4);

%%


FEM.PlateNodeDof = 5;%% number of degrees of freedom per node for the plate

FEM.nodeCoordinates = FEM.nodesCord(:,2:3);

FEM.GDof = FEM.PlateNodeDof*size(FEM.nodeCoordinates,1);

FEM.typeplate='CQUAD8'; % 8 node for each plate element

FEM.Dimension='2D';

% selected integration
FEM.GaussPointShear='2by2';
FEM.GaussPointBend='3by3';


FEM.NodeNumber=size(FEM.nodesCord,1);
FEM.elementNumber=size(FEM.elementNodes,1);

FEM.nodeCoordinates_label = zeros(size(FEM.nodeCoordinates,1),4);

FEM.nodeCoordinates_label(:,2:3) = FEM.nodeCoordinates;



HH = figure(200);hold on;
plot(FEM.nodeCoordinates(:,1),FEM.nodeCoordinates(:,2),'ko')
FEM.nodeCoordinates_label(:,1)  = 1:size(FEM.nodeCoordinates,1);

patch_plot(FEM.elementNodes,FEM.nodeCoordinates_label,200,'skin');axis image;

%% Composite Material Properties

% MatProperty;
psi2pa=6894.75729;
% psipa=1;
Mat.kappa=5/6;
Mat.E1=181e9;
Mat.E2=10.273e9;
Mat.G12=7.1705e9;
Mat.G13=4e9;
Mat.G23=4e9;
Mat.v12=0.28;
Mat.v21=Mat.E2/Mat.E1*Mat.v12;
Mat.density=1800;
Mat.alpha1=-0.04e-6;
Mat.alpha2=16.7e-6;


% to from a 3D orthotropic material
Mat.E3  = Mat.E1;
Mat.v13 = Mat.v12;
Mat.v31 = Mat.v13/Mat.E1*Mat.E3;
Mat.v23 = Mat.v12;
Mat.v32 = Mat.v23*Mat.E3/Mat.E2;


%% VAT laminates

VAT.VAT_type = 'LV';
flag = 'SYM';

half_number = 4;

T01 =[45 0]; % Theta_0 and Theta_1

T0T1 = VAT_fiber(T01,half_number,flag);



% T0T1 = VAT_fiber(T01,half_number,flag);

Laminate.layer_thickness = 1.272e-4;
Laminate.layer= half_number*2;


Laminate.thickness=Laminate.layer_thickness*ones(1,Laminate.layer);%% Uniform thickness
%
Plate.thickness= sum(Laminate.thickness);

center =[Plate.length Plate.width]/2;

physical_length = Plate.length; % this is the length of the plate for fiber path angle calculation


% Mass matrix
[MassPlate,MelemPlate] = LinearMassLaminatedPlate(FEM,Mat,Plate);

% Stiffness matrix
[Kplate,Kelemp] = LinearStiffnessLaminatedPlate_VAT_v2_X_constant_angle(Mat, Plate,FEM,Laminate,T0T1,center(1),physical_length);

%% Static analysis
% Boundary condition used in static analysis


FEM.BCtype='SFSF'; %% v is free in four edges - Wu's example

[ActiveDof,Constrained_Dof]=InPlaneBCPlate5Dof(FEM); % set degrees of freedom to zero
%

xx=FEM.nodeCoordinates(:,1);
yy=FEM.nodeCoordinates(:,2);

nodeNum=size(FEM.nodeCoordinates,1);

SideNodesNumList_RHS = find(xx==Plate.length)';

SideNodesNumList_LHS = find(xx==0)';

% specify boundary conditions x=a, u = -1e-4;
% note the present degrees of freedom order is: w,theta_x,theta_y, u, v 

Dof_fixed_disp_LHS =  SideNodesNumList_LHS + 3*nodeNum ;
Dof_fixed_disp_RHS =  SideNodesNumList_RHS + 3*nodeNum;

ActiveDof_not_fixed  = setdiff(ActiveDof,[Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);

%        
fixed_u =[ .5e-4*ones(length(Dof_fixed_disp_LHS),1) ; -.5e-4*ones(length(Dof_fixed_disp_RHS),1) ];

K21 = Kplate( ActiveDof_not_fixed ,  [Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);
K22 = Kplate( ActiveDof_not_fixed ,  ActiveDof_not_fixed);

final_unkowns = -K22\(K21*fixed_u);

% discrete analytical stiffness matrix
% [ K11 K12 ]{X1}  {F1}
% [         ]    =
% [ K21 K22 ]{X2}  {0}
%
% K21 X1 + K22 X2 = 0; BASED ON X1 -> X2
%
%
FEM.displacement=zeros(1,FEM.GDof);

FEM.displacement( [Dof_fixed_disp_LHS Dof_fixed_disp_RHS] ) = fixed_u';
FEM.displacement( ActiveDof_not_fixed) =  final_unkowns';

% plot deformation of the wing

postprocessISOTROPIC;

[FEM.stress,FEM.strain]=StressRecoveryPlate_VAT_center_average_v2_constant_angle(FEM,Laminate,Mat,Plate,T0T1,center,Plate.length);

calculate_stress_mat_cord

KGplate = GeometryStiffnessPlate_stress_recovery(FEM,Plate,Laminate);

%% Buckling eigenvalue anlaysis 
FEM.BCtype = 'SSSS-3'; %% Four simple supported sides

[ActiveDof,SideDof] = EssentialBCPlate5Dof(FEM);

Solver = 'buckling';
tic
sparsKplate=  sparse(Kplate(ActiveDof,ActiveDof));
sparsKGplate=  sparse(KGplate(ActiveDof,ActiveDof));
[V,D] = eigs(sparsKplate,-sparsKGplate,20,'sm');toc


[DD,modeNo]=sort((diag(D)));
loadfactor=DD(1:10);
VVsort=V(:,modeNo);
frequency=loadfactor;
%-------------minimal load factor---------------
mineiglabel=find(min(abs(DD))==-DD);
if isempty(mineiglabel)==1
    mineiglabel=find(min(abs(DD))==DD);
end
plotmodenNo=mineiglabel;

modeshapeplate;

critical_loadfactor = DD(mineiglabel);

%%
% find out the stress for elements near edge, x=b/2
elements_of_interest = [];
num_interest = 1;
for elem = 1:FEM.elementNumber

    nodes_4_element = FEM.elementNodes(elem,:);
    
    element_yes = 0;
    
    for node_num = 1:length(nodes_4_element)

        cord_temp =  label2cord(nodes_4_element(node_num),FEM.nodesCord);

        if abs((cord_temp(2) - Plate.length)/Plate.length)<1e-5  % x=a, RHS

            element_yes = 1;

        end

    end

    if   element_yes ==1
        elements_of_interest(num_interest) = elem;
        
        num_interest = num_interest+1;
    end
    
end

% force under given displacement
FORCE = Kplate*FEM.displacement';

% critical force applied in the right hand side
force_rhs = FORCE(Dof_fixed_disp_RHS )*critical_loadfactor;
%
% resultant force applied in the right hand side, N/m
Nxcr = sum(force_rhs)/Plate.width

% buckling load factor, Kcr
Kcr = Nxcr*Plate.width^2/(Mat.E1*Plate.thickness^3)

% ==============  END OF CODE  ============== %
