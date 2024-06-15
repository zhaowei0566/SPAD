function [Kbeam,Kelem,StiffenerLength]=...
    LaminatedStiffenerStiffness5Dof(FEM,Dst,Stiffener,Gausspointbend)
% 2 GAUSS POINTs FOR BOTH INTEGRATIONS
% The B matrix for the stiffener was modified. No curvature existed in the
% outplane deformation
% FOR CURVED BEAM, beam stiffness can be directly assembled when using
% nonlinear beam element.

disp('------------ ASSEMBLING STIFFENER STIFFNESS -------------');
FEM.typestiff='CBAR3';

% Kbeam---initialize to zeros matrix
 
% [Kbeam]=zeros(FEM.Stiffener.gdof,FEM.Stiffener.gdof);
% [Mbeam]=zeros(FEM.Stiffener.gdof,FEM.Stiffener.gdof);
Kbeam=zeros(Stiffener.nodedof*Stiffener.nodenum,Stiffener.nodedof*Stiffener.nodenum);
numberElements=size(FEM.Stiffener.element,1);
% DispFlag=zeros(numberElements, FEM.GDof/5);
% numberNodes=size(FEM.Stiffener.nodecord,1);
StiffenerLength=0;

for elem=1:numberElements%% number of stiffener element

%     disp(['-------------stiffener No#' num2str(e) '-------------']);
    
    stiffenerNODE=FEM.Stiffener.element(elem,2:end); %% Node NO. for one element
    XY = Stiffener.PointsCoord(stiffenerNODE,2:3);
    
    elementDof2=stiffenerNODE;
            
    elementDof=[elementDof2 elementDof2+Stiffener.nodenum  elementDof2+Stiffener.nodenum*2 ...
        elementDof2+Stiffener.nodenum*3 elementDof2+Stiffener.nodenum*4];
    
    Kelem=zeros(length(elementDof),length(elementDof));

    % %------ Shear Strain Energy ------------------
    [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointbend);
    
    for ee=1:size(GaussWeights,1) %% Gaussian Points number
        
        GaussPoint=GaussLocations(ee);
        xi=GaussPoint;
        
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);
        
        %         [Jacob,invJacob,XYderivatives]=...
        %             JacobianFEMCurved(naturalderivatives,XYZ);
        detJ=sqrt((naturalderivatives*XY(:,1))^2+(naturalderivatives*XY(:,2))^2);
        
        Curvature = StiffenerCurvature(XY,detJ, naturalderivatives, d2Nds2);% Curvature
        
        
        
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
        

        Kelem=Kelem+Lst'*Dst'*Lst*detJ*GaussWeights(ee);
        
    end
    
    
    
    %% SHEAR STRAIN ENERGY
    Kbeam(elementDof,elementDof)=Kbeam(elementDof,elementDof)+Kelem;
    
    StiffenerLength=StiffenerLength+2*detJ;
    %     disp('Tangential Angle is: ')
    %     alpha/pi*180
end
%% END OF DIARY