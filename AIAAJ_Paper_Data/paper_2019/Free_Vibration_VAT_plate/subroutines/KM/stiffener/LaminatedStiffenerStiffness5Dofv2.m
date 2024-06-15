function [Kbeam,Kelem,Mbeam,Melem,StiffenerLength]=LaminatedStiffenerStiffness5Dofv2(Stiffene,FEM,Mat,Dst,ebar,Stiffener)
% 2 GAUSSIAN POINTS FOR BENDING, AND 1 FOR SHEAR
% The B matrix for the stiffener was modified. No curvature existed in the
% outplane deformation
% FOR CURVED BEAM, beam stiffness can be directly assembled when using
% nonlinear beam element.


disp('-------------Material Stiffness-------------');
FEM.typestiff='CBAR3';

Gausspointbend='2';
Gausspointshear='1';
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

kappa=Mat.kappa;

      Amatrix= Stiffene.Amatrix;
       Bmatrix= Stiffene.Bmatrix;
       Dmatrix= Stiffene.Dmatrix;
       Ashear= Stiffene.Ashear;
       Qreduced=Stiffene.Qreduced;
       
        hs=Stiffener.height;
        A11=Amatrix(1,1);A16=Amatrix(1,2);A66=Amatrix(2,2);
        B11=Bmatrix(1,1);B16=Bmatrix(1,2);B66=Bmatrix(2,2);
        D11=Dmatrix(1,1);D16=Dmatrix(1,2);D66=Dmatrix(2,2);
        A55=Ashear;
        
        
        
         DsBend=[hs*A11,             0,          0             hs*ebar*A11,              -hs*B16;
        
                      0,              0   ,      0,                  0,                      0
        
                      0,              0,         0      ,           0,                        0;
        
                    hs*ebar*A11,      0,        0  ,     (hs^3/12+hs*ebar^2)*A11,     -hs*ebar*B16;
        
                    -hs*B16,          0  ,      0   ,        -hs*ebar*B16,        1/6*(Qreduced(2,2)+Qreduced(3,3))*Stiffener.width^3*hs];
                
                
                  
        DsShear=[0,                 0,             hs*A16*kappa,               0,                0;
         
                0,               hs*A55*6/5,         0,                       0,                 hs*ebar*A55;
        
             hs*A16*kappa,         0,              hs*A66*kappa,           hs*ebar*A16*kappa,      -hs*B66*kappa;
        
               0,                  0,               hs*ebar*A16*kappa,          0,                    0;
        
               0,             hs*ebar*A55,          -hs*B66*kappa ,             0,                    0];
        
        

for e=1:numberElements%% number of stiffener element

    disp(['-------------stiffener No#' num2str(e) '-------------']);
    
    stiffenerNODE=FEM.Stiffener.element(e,2:end); %% Node NO. for one element
    XYZ= Stiffener.PointsCoord(stiffenerNODE-4,2:3);
    
    elementDof2=stiffenerNODE-4;
            
    elementDof=[elementDof2 elementDof2+Stiffener.nodenum  elementDof2+Stiffener.nodenum*2 ...
        elementDof2+Stiffener.nodenum*3 elementDof2+Stiffener.nodenum*4];
    
    Kelem=zeros(length(elementDof),length(elementDof));
    Melem=zeros(size(elementDof2,2)*5,size(elementDof2,2)*5);
    
    
    %% ------ Bending Strain Energy ------------------
    [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointbend);
    
    for ee=1:size(GaussWeights,1) %% Gaussian Points number
        
        GaussPoint=GaussLocations(ee);
        xi=GaussPoint;
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);
        
        %         [Jacob,invJacob,XYderivatives]=...
        %             JacobianFEMCurved(naturalderivatives,XYZ);
        detJ=sqrt((naturalderivatives*XYZ(:,1))^2+(naturalderivatives*XYZ(:,2))^2);
        Curvature=StiffenerCurvature(XYZ,detJ, naturalderivatives, d2Nds2);% Curvature
        
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
        

        
   
        Kelem=Kelem+Lst'*DsBend'*Lst*detJ*GaussWeights(ee);
        
    end
    
    %% ------ Bending Strain Energy ------------------
    [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointshear);
    
    for ee=1:size(GaussWeights,1) %% Gaussian Points number
        
        GaussPoint=GaussLocations(ee);
        xi=GaussPoint;
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);
        
        %         [Jacob,invJacob,XYderivatives]=...
        %             JacobianFEMCurved(naturalderivatives,XYZ);
        detJ=sqrt((naturalderivatives*XYZ(:,1))^2+(naturalderivatives*XYZ(:,2))^2);
        Curvature=StiffenerCurvature(XYZ,detJ, naturalderivatives, d2Nds2);% Curvature
        
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
   
        Kelem=Kelem+Lst'*DsShear'*Lst*detJ*GaussWeights(ee);
        
    end
    
    
    %% SHEAR STRAIN ENERGY
    Kbeam(elementDof,elementDof)=Kbeam(elementDof,elementDof)+Kelem;
    StiffenerLength=StiffenerLength+2*detJ;
    %     disp('Tangential Angle is: ')
    %     alpha/pi*180
end
%% END OF DIARY