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

% map
ETA=[-1/2 1/2 ];
eta11=ETA(1);xi11=1;
eta12=ETA(2);xi12=0;
eta13=ETA(1);xi13=-1;

% right curve -- -1~1

eta21=1-abs(eta11); xi21=1;
eta22=1-abs(eta12); xi22=0;
eta23=1-abs(eta13); xi23=-1;


deta1=[eta11 eta12 eta13];
dxi1=[xi11 xi12 xi13];

xi1=linspace(-1,1,31);
eta1=interp1(dxi1,deta1,xi1,'spline');
figure(100);hold on;plot(eta1,xi1,'r-');

xi2=xi1;
% for jj=1:length(eta1)
%
%     if eta1(jj)<0
%         eta2=
%
% end
% eta2=1-abs(1-abs(eta1));
eta2=-eta1;
figure(100);hold on;plot(eta2,xi2,'b-');axis image;box on;axis([-1 1 -1 1])


panel_point_1=[0,0];
panel_point_2=[0.6,0];
panel_point_3=[0.6,0.8];
panel_point_4=[0.0,0.8];

clear XXstiffener YYstiffener

% stiffener node coordinate in physical space
for node_num=1:length(xi1)
    
    XXstiffener(1,node_num)=panel_point_1(1)*1/4*(1-eta1(node_num))*(1-xi1(node_num))+...
        panel_point_2(1)*1/4*(1+eta1(node_num))*(1-xi1(node_num))+...
        panel_point_3(1)*1/4*(1+eta1(node_num))*(1+xi1(node_num))+...
        panel_point_4(1)*1/4*(1-eta1(node_num))*(1+xi1(node_num));
    
    XXstiffener(2,node_num)=panel_point_1(1)*1/4*(1-eta2(node_num))*(1-xi2(node_num))+...
        panel_point_2(1)*1/4*(1+eta2(node_num))*(1-xi2(node_num))+...
        panel_point_3(1)*1/4*(1+eta2(node_num))*(1+xi2(node_num))+...
        panel_point_4(1)*1/4*(1-eta2(node_num))*(1+xi2(node_num));
    
    YYstiffener(1,node_num)=panel_point_1(2)*1/4*(1-eta1(node_num))*(1-xi1(node_num))+...
        panel_point_2(2)*1/4*(1+eta1(node_num))*(1-xi1(node_num))+...
        panel_point_3(2)*1/4*(1+eta1(node_num))*(1+xi1(node_num))+...
        panel_point_4(2)*1/4*(1-eta1(node_num))*(1+xi1(node_num));
    
    YYstiffener(2,node_num)=panel_point_1(2)*1/4*(1-eta2(node_num))*(1-xi2(node_num))+...
        panel_point_2(2)*1/4*(1+eta2(node_num))*(1-xi2(node_num))+...
        panel_point_3(2)*1/4*(1+eta2(node_num))*(1+xi2(node_num))+...
        panel_point_4(2)*1/4*(1-eta2(node_num))*(1+xi2(node_num));
    
end

% plate node coordinate in physical space
for platenodeNum=1:length(Xcoord)
    
    platenodeX(platenodeNum)=panel_point_1(1)*1/4*(1-Xcoord(platenodeNum))*(1-Ycoord(platenodeNum))+...
        panel_point_2(1)*1/4*(1+Xcoord(platenodeNum))*(1-Ycoord(platenodeNum))+...
        panel_point_3(1)*1/4*(1+Xcoord(platenodeNum))*(1+Ycoord(platenodeNum))+...
        panel_point_4(1)*1/4*(1-Xcoord(platenodeNum))*(1+Ycoord(platenodeNum));
    
    platenodeY(platenodeNum)=panel_point_1(2)*1/4*(1-Xcoord(platenodeNum))*(1-Ycoord(platenodeNum))+...
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
% rectangle('Position',[0,0,0.600,0.800]);hold on;axis image;
% Xcord=[panel_point_1(1) panel_point_2(1) panel_point_3(1) panel_point_4(1)];
% Ycord=[panel_point_1(2) panel_point_2(2) panel_point_3(2) panel_point_4(2)];
% fill(Xcord,Ycord,'r','LineWidth',2);hold on;
[xlength,ylength]=meshgrid(0:panel_point_2(1)/12:panel_point_2(1),...
    0:panel_point_3(2)/12:panel_point_3(2));
zlength=xlength.*0+ylength.*0;
mesh(xlength,ylength,zlength);hold on;view(2)
figure(200);hold on;
plot(platenodeX,platenodeY,'k.');axis image; hold on;
plot(XXstiffener(1,:),YYstiffener(1,:),'k-','LineWidth',2);hold on;
plot(XXstiffener(2,:),YYstiffener(2,:),'b-','LineWidth',2);hold off;axis image;
