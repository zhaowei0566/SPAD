% MAIN CODE FOR

% Example III.A in the paper

clear all;warning off;format long;
close all

%% ==== add subroutines ======
addpath(genpath([pwd filesep 'subroutines']));
%
%
%%

Plate.length = .4;
Plate.width  = .3;

%% Generate plate mesh

FEM = mesh_QUAD8_v2(40,30);

patch_plot(FEM.elementNodes,FEM.nodesCord,101,'skin')

FEM.nodesCord(:,2)=FEM.nodesCord(:,2)/max(FEM.nodesCord(:,2))*Plate.length;
FEM.nodesCord(:,3)=FEM.nodesCord(:,3)/max(FEM.nodesCord(:,3))*Plate.width;

Xcoord=FEM.nodesCord(:,2);
Ycoord=FEM.nodesCord(:,3);
Zcoord=FEM.nodesCord(:,4);
%%
%-----------------------------------------
Stru.length=max(abs(Xcoord));
Stru.width=max(abs(Ycoord));


%-----------D.O.F for each Node--------------
FEM.PlateNodeDof=5;%% for plate
Stiffener.nodedof=5; %% for stiffener
FEM.nodeCoordinates=FEM.nodesCord(:,2:3);

FEM.GDof = FEM.PlateNodeDof*size(FEM.nodeCoordinates,1);
FEM.typeplate='CQUAD8';
FEM.Dimension='2D';
switch FEM.typeplate
    case 'CQUAD4'
        FEM.GaussPointShear='1by1';
        FEM.GaussPointBend='2by2';
    case 'CQUAD8'
        FEM.GaussPointShear='2by2';
        FEM.GaussPointBend='3by3';
end
%---------------------------------------
FEM.NodeNumber=size(FEM.nodesCord,1);
FEM.elementNumber=size(FEM.elementNodes,1);


FEM.nodeCoordinates_label = zeros(size(FEM.nodeCoordinates,1),4);

FEM.nodeCoordinates_label(:,2:3) = FEM.nodeCoordinates;
figure(200);hold on;
plot(FEM.nodeCoordinates(:,1),FEM.nodeCoordinates(:,2),'ko')
FEM.nodeCoordinates_label(:,1)  = 1:size(FEM.nodeCoordinates,1);

patch_plot(FEM.elementNodes,FEM.nodeCoordinates_label,200,'skin')

%% ===========-Composite Material Properties ==========================

% MatProperty;
psi2pa=6894.75729;
% psipa=1;
Mat.kappa=5/6;
Mat.E1=126.3E9;
Mat.E2=8.765e9;
Mat.G12=4.92e9;
Mat.G13=4.92E9;
% Mat.G23=3.35E9;%Mat.G13;
Mat.G23=Mat.G13;
Mat.v12=0.334;
Mat.v21=Mat.E2/Mat.E1*Mat.v12;
Mat.density=1580;
Mat.alpha1=-0.04e-6;
Mat.alpha2=16.7e-6;

%% Parameteric studies

number_of_plies= 51;

Theta = linspace(-90,90,number_of_plies);

Critical_loadfactor = zeros(number_of_plies,number_of_plies);

Buckling_parameter = zeros(number_of_plies,number_of_plies);

T0T1 = [30 10
    -30 -10
    30 10
    -30 -10
    90 90
    90 90
    -30 -10
    30 10
    -30 -10
    30 10];

Laminate.layer_thickness = 0.00175/size(T0T1,1);%% 1.27e-4;

Laminate.layer= size(T0T1,1);

%%
VAT.VAT_type = 'LV';
disp('----------Laminates for Plate----------');


Laminate.thickness=Laminate.layer_thickness*ones(1,Laminate.layer);%% Uniform thickness
%
Stru.thickness= sum(Laminate.thickness);


FEM.typeplate='CQUAD8';


center = [Plate.length Plate.width]/2;
physical_length  = Plate.length;

[Kplate,Kelemp]=LinearStiffnessLaminatedPlate_VAT_v2_X_constant_angle(Mat, Stru,FEM,Laminate,T0T1,center,physical_length);


[MassPlate,MelemPlate] = LinearMassLaminatedPlate(FEM,Mat,Stru);


Ktotal = Kplate;

Mtotal = MassPlate;


%% Vibration ANALYSIS
 
SOL_flag = '103';% 102 - pure plate; 105 - buckling; 103 - free vibration; 106 - presstressed vibration

Solver='vibration';

ActiveDof = 1:FEM.GDof; %% Free vibration

[V,D]=eigs(Ktotal(ActiveDof,ActiveDof),Mtotal(ActiveDof,ActiveDof),40,'sm');

[DD,modeNo]=sort(diag(D));
eigenvalue=DD(1:20);
VVsort=V(:,modeNo);
cycleFreq=sqrt(eigenvalue);
frequency=cycleFreq/2/pi;

disp('==== The real part is:')
real(frequency)


%-------------minimal load factor---------------
mineiglabel=find(min(abs(DD))==-DD);
if isempty(mineiglabel)==1
    mineiglabel=find(min(abs(DD))==DD);
end

plotmodenNo=mineiglabel;

modeshapeplate;

%% End of code

