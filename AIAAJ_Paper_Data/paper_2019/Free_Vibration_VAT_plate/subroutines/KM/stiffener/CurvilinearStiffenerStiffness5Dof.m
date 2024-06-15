function [Kbeam,Kelem,StiffenerLength]=CurvilinearStiffenerStiffness5Dof(FEM,Mat,ebar,Stiffener,gausspoint)
% ONLY 1 GAUSS POINT FOR INTEGRATION
% The B matrix for the stiffener was modified. No curvature existed in the
% outplane deformation
% FOR CURVED BEAM, beam stiffness can be directly assembled when using
% nonlinear beam element.
disp('-------------Material Stiffness-------------');
% Rewrite the Nshape Function
format long;
FEM.typeplate='CBAR3';
E=Mat.E;
G=Mat.G;
kappa=Mat.kappa;

As=Stiffener.width*Stiffener.height;
Istiffener11=1/12*Stiffener.width*Stiffener.height^3+ebar^2*As;
Jstiffener=1/3*Stiffener.height*Stiffener.width^3;
% Jstiffener=StiffenerTorsionConstant1D(Stiffener);
% Jstiffener=0;
Gausspointbend=gausspoint;
% Kbeam---initialize to zeros matrix
% 1 Gauss Point instead
% [Kbeam]=zeros(FEM.Stiffener.gdof,FEM.Stiffener.gdof);
% [Mbeam]=zeros(FEM.Stiffener.gdof,FEM.Stiffener.gdof);
[Kbeam]=zeros(Stiffener.nodedof*Stiffener.nodenum,Stiffener.nodedof*Stiffener.nodenum);
[Mbeam]=zeros(Stiffener.nodedof*Stiffener.nodenum,Stiffener.nodedof*Stiffener.nodenum);


numberElements=size(FEM.Stiffener.element,1);
% DispFlag=zeros(numberElements, FEM.GDof/5);
% numberNodes=size(FEM.Stiffener.nodecord,1);
StiffenerLength=0;


for e=1:numberElements%% number of stiffener element

    disp(['-------------stiffener No#' num2str(e) '-------------']);
    
    stiffenerNODE=FEM.Stiffener.element(e,2:end); %% Node NO. for one element
    XYZ= Stiffener.PointsCoord(stiffenerNODE-4,2:3);
    
    elementDof2=stiffenerNODE-4;
            
    elementDof=[elementDof2 elementDof2+Stiffener.nodenum  elementDof2+Stiffener.nodenum*2 ...
        elementDof2+Stiffener.nodenum*3 elementDof2+Stiffener.nodenum*4];
    
    Kelem=zeros(length(elementDof),length(elementDof));
    Melem=zeros(size(elementDof2,2)*5,size(elementDof2,2)*5);
    %
    % %------ Shear Strain Energy ------------------
    [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointbend);
    
    for ee=1:size(GaussWeights,1) %% Gaussian Points number
        
        GaussPoint=GaussLocations(ee);
        xi=GaussPoint;
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);
        
        %         [Jacob,invJacob,XYderivatives]=...
        %             JacobianFEMCurved(naturalderivatives,XYZ);
        detJ=sqrt((naturalderivatives*XYZ(:,1))^2+(naturalderivatives*XYZ(:,2))^2);
        Curvature=StiffenerCurvature(XYZ,detJ, naturalderivatives, d2Nds2);% Curvature
        
        Es=E;
        An=kappa*As;
        Ab=kappa*As;
        Gs=G;
        Jt=Jstiffener;
        In=Istiffener11;
        
        %% ------------------------------------------------------
        Lst=zeros(5,15);%% Ali's paper
        dNds=naturalderivatives/detJ;
        N=shape';

        Lst(1,10:12)= dNds;
        Lst(1,13:15)= N*Curvature;
        
        Lst(2,10:12)= - N*Curvature;
        Lst(2,13:15)= dNds;
        
        Lst(3,1:3) =  dNds;
        Lst(3,4:6) = N;
        
        Lst(4,4:6) = dNds;
        Lst(4,7:9) = N*Curvature;
        
        
        Lst(5,4:6)=  - N*Curvature;
        Lst(5,7:9)=  dNds ;
   
        
        Dst=...
            [Es*As,       0,       0,  Es*As*ebar,   0;
            0,          Gs*An,   0,      0,    Gs*An*ebar;
            0,            0,      Gs*Ab,   0,        0;
            Es*As*ebar,   0,       0,    Es*In,      0;
            0,       Gs*An*ebar,  0,      0,      Gs*Jt];
        
        
        
        Kelem=Kelem+Lst'*Dst'*Lst*detJ*GaussWeights(ee);
        
    end
    %% SHEAR STRAIN ENERGY
    Kbeam(elementDof,elementDof)=Kbeam(elementDof,elementDof)+Kelem;
    StiffenerLength=StiffenerLength+2*detJ;
    %     disp('Tangential Angle is: ')
    %     alpha/pi*180
end
%% END OF DIARY