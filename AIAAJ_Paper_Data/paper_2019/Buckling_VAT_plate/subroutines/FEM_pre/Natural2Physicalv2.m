% Mesh in natural space
clear all;warning off;format long;

global FEM;
global Stru;

%% Geometry Coordinates
Stru.bdfname='naturalspace12by12.bdf';
Stru.pathfile='Z:\Paper\Stiffened Plate\StiffenedPlateCode\InputBDF';
FEM=geometryinfoplate(Stru);

Xcoord=FEM.nodesCord(:,2);
Ycoord=FEM.nodesCord(:,3);
Zcoord=FEM.nodesCord(:,4);

% For natural space, it's -1 at the first x coordinate
[xlength,ylength]=meshgrid(-1:2/((length(find(Xcoord==-1))-1)/2):1,...
    -1:2/((length(find(Xcoord==-1))-1)/2):1);
zlength=xlength.*0+ylength.*0;
figure(100)
mesh(xlength,ylength,zlength);hold on;view(2)
plot(Xcoord,Ycoord,'k.');hold on;
title('Nodes of Stiffened Plate','FontSize',15);axis image;hold on;

% map parameter
DesignStepub=0.95;
DesignSteplr=0.4;
% 4 points coordinates
panel_point_1=[0,0];
panel_point_2=[0.8,0.0];
panel_point_3=[0.8,1.75];
panel_point_4=[0.25,1.95];
%%  curve 1
% point 1 and point 3
XI1=[-0.5 -0.5]; % coordinate xi for two points at the side of the natural space
ETA1=[-1 1]; % coordinate eta for two points at the side of the natural space
xi11=XI1(1);eta11=ETA1(1);
xi13=XI1(2);eta13=ETA1(2);

direction_alpha_point1=[1,0];
curvature_d_point1=DesignStepub;
% the point 2
xi12=0.5*(xi11+xi13)+direction_alpha_point1(1)*curvature_d_point1;
eta12=0.5*(eta11+eta13)+direction_alpha_point1(2)*curvature_d_point1;

point1_xi=[xi11 xi12 xi13];
point1_eta=[eta11 eta12 eta13];

% spline points

Exacteta1=linspace(eta11,eta13,31);
xi1=interp1(point1_eta,point1_xi,Exacteta1,'spline');

figure(100);plot(point1_xi,point1_eta,'bo');hold on;
plot(xi1,Exacteta1,'k-','LineWidth',2);axis image;axis([-1 1 -1 1]);hold on;
plot([xi11 xi13],[eta11 eta13],'b-');

%%  curve 2
% point 1 and point 3
XI2=[0.4  0.4]; % coordinate xi for two points at the side of the natural space
ETA2=[-1 1]; % coordinate eta for two points at the side of the natural space
xi21=XI2(1);eta21=ETA2(1);
xi23=XI2(2);eta23=ETA2(2);

direction_alpha_point2=[-1,0];
curvature_d_point2=DesignStepub;
% the point 2
xi22=0.5*(xi21+xi23)+direction_alpha_point2(1)*curvature_d_point2;
eta22=0.5*(eta21+eta23)+direction_alpha_point2(2)*curvature_d_point2;

point2_xi=[xi21 xi22 xi23];
point2_eta=[eta21 eta22 eta23];

% spline points

Exacteta2=linspace(eta21,eta23,31);
xi2=interp1(point2_eta,point2_xi,Exacteta2,'spline');

figure(100);hold on;plot(point2_xi,point2_eta,'ro');hold on;
plot(xi2,Exacteta2,'k-','LineWidth',2);axis image;axis([-1 1 -1 1]);hold on;
plot([xi21 xi23],[eta21 eta23],'r-');

%% Curve 3
% point 1 and point 3
XI3=[-1  1]; % coordinate xi for two points at the side of the natural space
ETA3=[-0.6 -0.6]; % coordinate eta for two points at the side of the natural space
xi31=XI3(1);eta31=ETA3(1);
xi33=XI3(2);eta33=ETA3(2);

direction_alpha_point3=[0,1];
curvature_d_point3=DesignSteplr;
% the point 2
xi32=0.5*(xi31+xi33)+direction_alpha_point3(1)*curvature_d_point3;
eta32=0.5*(eta31+eta33)+direction_alpha_point3(2)*curvature_d_point3;

point3_xi=[xi31 xi32 xi33];
point3_eta=[eta31 eta32 eta33];

% spline points

Exactxi3=linspace(xi31,xi33,31);
eta3=interp1(point3_xi,point3_eta,Exactxi3,'spline');

figure(100);hold on;plot(point3_xi,point3_eta,'go');hold on;
plot(Exactxi3,eta3,'k-','LineWidth',2);axis image;axis([-1 1 -1 1]);hold on;
plot([xi31 xi33],[eta31 eta33],'r-');
%% Curve 4
%% Curve 3
% point 1 and point 3
XI4=[-1  1]; % coordinate xi for two points at the side of the natural space
ETA4=[0.5 0.5]; % coordinate eta for two points at the side of the natural space
xi41=XI4(1);eta41=ETA4(1);
xi43=XI4(2);eta43=ETA4(2);

direction_alpha_point4=[0,-1];
curvature_d_point4=DesignSteplr;
% the point 2
xi42=0.5*(xi41+xi43)+direction_alpha_point4(1)*curvature_d_point4;
eta42=0.5*(eta41+eta43)+direction_alpha_point4(2)*curvature_d_point4;

point4_xi=[xi41 xi42 xi43];
point4_eta=[eta41 eta42 eta43];

% spline points

Exactxi4=linspace(xi41,xi43,31);
eta4=interp1(point4_xi,point4_eta,Exactxi4,'spline');

figure(100);hold on;plot(point4_xi,point4_eta,'go');hold on;
plot(Exactxi4,eta4,'k-','LineWidth',2);axis image;axis([-1 1 -1 1]);hold on;
plot([xi41 xi43],[eta41 eta43],'r-');
%% Transformation through mapping functions
clear XXstiffener YYstiffener

% stiffener node coordinate in physical space
for node_num=1:length(xi1)
    
    %% stiffener - 1
    XXstiffener(1,node_num)=...
        panel_point_1(1)*1/4*(1-Exacteta1(node_num))*(1-xi1(node_num))+...
        panel_point_2(1)*1/4*(1-Exacteta1(node_num))*(1+xi1(node_num))+...
        panel_point_3(1)*1/4*(1+Exacteta1(node_num))*(1+xi1(node_num))+...
        panel_point_4(1)*1/4*(1+Exacteta1(node_num))*(1-xi1(node_num));
    
    YYstiffener(1,node_num)=...
        panel_point_1(2)*1/4*(1-Exacteta1(node_num))*(1-xi1(node_num))+...
        panel_point_2(2)*1/4*(1-Exacteta1(node_num))*(1+xi1(node_num))+...
        panel_point_3(2)*1/4*(1+Exacteta1(node_num))*(1+xi1(node_num))+...
        panel_point_4(2)*1/4*(1+Exacteta1(node_num))*(1-xi1(node_num));
    
    %% stiffener - 2
    XXstiffener(2,node_num)=...
        panel_point_1(1)*1/4*(1-Exacteta2(node_num))*(1-xi2(node_num))+...
        panel_point_2(1)*1/4*(1-Exacteta2(node_num))*(1+xi2(node_num))+...
        panel_point_3(1)*1/4*(1+Exacteta2(node_num))*(1+xi2(node_num))+...
        panel_point_4(1)*1/4*(1+Exacteta2(node_num))*(1-xi2(node_num));
    
    YYstiffener(2,node_num)=...
        panel_point_1(2)*1/4*(1-Exacteta2(node_num))*(1-xi2(node_num))+...
        panel_point_2(2)*1/4*(1-Exacteta2(node_num))*(1+xi2(node_num))+...
        panel_point_3(2)*1/4*(1+Exacteta2(node_num))*(1+xi2(node_num))+...
        panel_point_4(2)*1/4*(1+Exacteta2(node_num))*(1-xi2(node_num));
    %% stiffener - 3
    XXstiffener(3,node_num)=...
        panel_point_1(1)*1/4*(1-eta3(node_num))*(1-Exactxi3(node_num))+...
        panel_point_2(1)*1/4*(1-eta3(node_num))*(1+Exactxi3(node_num))+...
        panel_point_3(1)*1/4*(1+eta3(node_num))*(1+Exactxi3(node_num))+...
        panel_point_4(1)*1/4*(1+eta3(node_num))*(1-Exactxi3(node_num));
    
    YYstiffener(3,node_num)=...
        panel_point_1(2)*1/4*(1-eta3(node_num))*(1-Exactxi3(node_num))+...
        panel_point_2(2)*1/4*(1-eta3(node_num))*(1+Exactxi3(node_num))+...
        panel_point_3(2)*1/4*(1+eta3(node_num))*(1+Exactxi3(node_num))+...
        panel_point_4(2)*1/4*(1+eta3(node_num))*(1-Exactxi3(node_num));
    
    %% stiffener - 4
    XXstiffener(4,node_num)=...
        panel_point_1(1)*1/4*(1-eta4(node_num))*(1-Exactxi4(node_num))+...
        panel_point_2(1)*1/4*(1-eta4(node_num))*(1+Exactxi4(node_num))+...
        panel_point_3(1)*1/4*(1+eta4(node_num))*(1+Exactxi4(node_num))+...
        panel_point_4(1)*1/4*(1+eta4(node_num))*(1-Exactxi4(node_num));
    
    YYstiffener(4,node_num)=...
        panel_point_1(2)*1/4*(1-eta4(node_num))*(1-Exactxi4(node_num))+...
        panel_point_2(2)*1/4*(1-eta4(node_num))*(1+Exactxi4(node_num))+...
        panel_point_3(2)*1/4*(1+eta4(node_num))*(1+Exactxi4(node_num))+...
        panel_point_4(2)*1/4*(1+eta4(node_num))*(1-Exactxi4(node_num));
end

% plate node coordinate in physical space
for platenodeNum=1:length(Xcoord)
    
    platenodeX(platenodeNum)=...
        panel_point_1(1)*1/4*(1-Xcoord(platenodeNum))*(1-Ycoord(platenodeNum))+...
        panel_point_2(1)*1/4*(1+Xcoord(platenodeNum))*(1-Ycoord(platenodeNum))+...
        panel_point_3(1)*1/4*(1+Xcoord(platenodeNum))*(1+Ycoord(platenodeNum))+...
        panel_point_4(1)*1/4*(1-Xcoord(platenodeNum))*(1+Ycoord(platenodeNum));
    
    platenodeY(platenodeNum)=...
        panel_point_1(2)*1/4*(1-Xcoord(platenodeNum))*(1-Ycoord(platenodeNum))+...
        panel_point_2(2)*1/4*(1+Xcoord(platenodeNum))*(1-Ycoord(platenodeNum))+...
        panel_point_3(2)*1/4*(1+Xcoord(platenodeNum))*(1+Ycoord(platenodeNum))+...
        panel_point_4(2)*1/4*(1-Xcoord(platenodeNum))*(1+Ycoord(platenodeNum));
    
end

% figure;
% [xplate,yplate]=meshgrid(platenodeX,platenodeY);
% zplate=0.*xplate;
% surf(xplate,yplate,zplate);view(2)
% title('Nodes of Stiffened Plate','FontSize',15);

figure(200);
[xlength,ylength]=meshgrid(0:panel_point_2(1)/12:panel_point_2(1),...
    0:panel_point_3(2)/12:panel_point_3(2));
zlength=xlength.*0+ylength.*0;
mesh(xlength,ylength,zlength);hold on;view(2)
plot(platenodeX,platenodeY,'MarkerSize',4,'Marker','square','LineStyle','none',...
    'Color',[0 0 0]);axis image; hold on;
plot(XXstiffener(1,:),YYstiffener(1,:),'k-','LineWidth',2);hold on;
plot(XXstiffener(2,:),YYstiffener(2,:),'b-','LineWidth',2);axis image;
plot(XXstiffener(3,:),YYstiffener(3,:),'b-','LineWidth',2);axis image;
plot(XXstiffener(4,:),YYstiffener(4,:),'b-','LineWidth',2);axis image;
figure(200);hold on;
% rectangle('Position',[0,0,0.600,0.800]);hold on;axis image;
Xcord=[panel_point_1(1) panel_point_2(1) panel_point_3(1) panel_point_4(1)];
Ycord=[panel_point_1(2) panel_point_2(2) panel_point_3(2) panel_point_4(2)];
fill(Xcord,Ycord,'r','LineWidth',3,'FaceColor','none','EdgeColor',[1,0,0]);