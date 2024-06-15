function [Kbeam,Kelem,StiffenerLength]=LaminatedStiffenerParallelV1(FEM,Ds,Stiffener)
% Stiffness Matrix for Curvilinear-Stiffener
% Cross section - rectanglar
%

FEM.typeplate='CBAR3';
%% ----Torsional Constant of Rect Section From NASTRAN----
% bb=Stiffener.width/2; aa=Stiffener.height/2;
% Jstiffener=aa*bb^3*(16/3-3.36*bb/aa*(1-bb^4/12/aa^4));
% Jstiffener=0;
%% ----------
nodenum=size(FEM.nodeCoordinates,1);
Gausspointbend='2';
% Kbeam---initialize to zeros matrix
% 1 Gauss Point instead
% [Kbeam]=zeros(FEM.Stiffener.gdof,FEM.Stiffener.gdof);
% [Mbeam]=zeros(FEM.Stiffener.gdof,FEM.Stiffener.gdof);
[Kbeam]=zeros(FEM.GDof,FEM.GDof);
[Mbeam]=zeros(FEM.GDof,FEM.GDof);



numberElements=size(FEM.Stiffener.element,1);
% DispFlag=zeros(numberElements, FEM.GDof/5);
% numberNodes=size(FEM.Stiffener.nodecord,1);
StiffenerLength=0;
DispFlagAll=[];
for e=1:numberElements%% number of stiffener element
    %
    disp(['-------------stiffener No#' num2str(e) '-------------']);
    
    stiffenerNODE=FEM.Stiffener.element(e,2:end); %% Node NO. for one element
    NshapeElem=zeros(3,FEM.GDof/5);
    clear XYZ;
    for eee=1:3 %% Three nodes in one beam element
        nodestifflag=stiffenerNODE(eee);
        TargetPlate=Stiffener.TargetPlate(nodestifflag-4);
        TargetGaussian=Stiffener.TargetGaussian(nodestifflag-4,:);
        xi=TargetGaussian(1);
        eta=TargetGaussian(2);
        
        Nshape(eee,:)=[1/4*(1-xi)*(1-eta)*(-xi-eta-1);
            1/4*(1+xi)*(1-eta)*(xi-eta-1);
            1/4*(1+xi)*(1+eta)*(xi+eta-1);
            1/4*(1-xi)*(1+eta)*(-xi+eta-1);
            1/2*(1-eta)*(1+xi)*(1-xi);
            1/2*(1+xi)*(1+eta)*(1-eta);
            1/2*(1+eta)*(1+xi)*(1-xi);
            1/2*(1-xi)*(1+eta)*(1-eta)]';
        
        NodeIndices(1+(eee-1)*8:eee*8)= FEM.elementNodes(TargetPlate,:);
        [aa,bb]=sort(NodeIndices(1+(eee-1)*8:eee*8));
        NshapeElem(eee,aa)=Nshape(eee,bb);
       
        XYZ(eee,:)=Stiffener.PointsCoord(nodestifflag-4,2:end);
        
        if eee==1
            dispflag1=find(NshapeElem(1,:)~=0);
        elseif eee==2
            dispflag2=find(NshapeElem(2,:)~=0);
        elseif eee==3
            dispflag3=find(NshapeElem(3,:)~=0);
        end
        
        %         eval(['DispFields' num2str(eee) '=find(NshapeElem(eee,:)~=0)']);
        %         eval(['DispFields=DispFields'  num2str(eee)]);
        %         eval(['dispflag' num2str(eee) '=DispFields']);
        
    end
    MatrixLength1=length(DispFlagAll);
    MatrixLength2=length(dispflag1)+length(dispflag2)+length(dispflag3);
    
    DispFlagAll(MatrixLength1+1:MatrixLength1+MatrixLength2)=[dispflag1 dispflag2 dispflag3]; 
    
    DispFlag=[dispflag1 dispflag2 dispflag3];
    elementDof2=unique(DispFlag);
    NshapeMatrixElem=NshapeElem(1:3,elementDof2);
    
    NshapeMatrix=zeros(15,size(elementDof2,2)*5);
    NshapeMatrix(1:3,1:size(elementDof2,2))= NshapeMatrixElem;
    NshapeMatrix(4:6,size(elementDof2,2)*1+1:size(elementDof2,2)*2)= NshapeMatrixElem;
    NshapeMatrix(7:9,size(elementDof2,2)*2+1:size(elementDof2,2)*3)= NshapeMatrixElem;
    NshapeMatrix(10:12,size(elementDof2,2)*3+1:size(elementDof2,2)*4)= NshapeMatrixElem;
    NshapeMatrix(13:15,size(elementDof2,2)*4+1:size(elementDof2,2)*5)= NshapeMatrixElem;
    
    elementDof=[elementDof2 elementDof2+FEM.NodeNumber  elementDof2+FEM.NodeNumber*2 ...
        elementDof2+FEM.NodeNumber*3 elementDof2+FEM.NodeNumber*4];
    
    
    Kelem=zeros(size(elementDof2,2)*5,size(elementDof2,2)*5);
    Melem=zeros(size(elementDof2,2)*5,size(elementDof2,2)*5);
    Kshear=zeros(size(elementDof2,2)*5,size(elementDof2,2)*5);
    Kbend=zeros(size(elementDof2,2)*5,size(elementDof2,2)*5);
    %
    % %------ Shear Strain Energy ------------------
    [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointbend);
    
    for ee=1:size(GaussWeights,1) %% Gaussian Points number
        
        FEM.typeplate='CBAR3';
        GaussPoint=GaussLocations(ee);
        xi=GaussPoint;
        eta=GaussPoint;
        [shape,naturalderivatives,d2Nds2]=shapefunction(xi,eta,FEM);
        
        %         [Jacob,invJacob,XYderivatives]=...
        %             JacobianFEM(naturalderivatives,FEM,NodeIndices);
        [Jacob,invJacob,XYderivatives]=...
            JacobianFEMCurved(naturalderivatives,XYZ);
        
        detJ=sqrt((naturalderivatives(1,:)*XYZ(:,1))^2+(naturalderivatives(2,:)*XYZ(:,2))^2);
        %         alpha=acos(naturalderivatives(1,:)*XYZ(:,1)/detJ);
        alpha=asin(naturalderivatives(2,:)*XYZ(:,2)/detJ);
        nx=cos(alpha);ny=sin(alpha);
        Curvature=StiffenerCurvature(XYZ,detJ, naturalderivatives, d2Nds2);% Curvature  
        %% ------------------------------------------------------
        Lst=zeros(5,15);%% Ali's paper
        dNds=naturalderivatives(1,:)/detJ;
        N=shape';
        
        Lst(1,10:12)= nx*dNds - N*Curvature*ny;
        Lst(1,13:15)= N*Curvature*nx + ny*dNds;
        
        Lst(2,10:12)= - N*Curvature*nx - ny*dNds;
        Lst(2,13:15)= nx*dNds - N*Curvature*ny;
        
        Lst(3,1:3)=dNds;
        Lst(3,4:6)= N*nx;
        Lst(3,7:9)= N*ny;
        
        Lst(4,4:6)= nx*dNds - N*Curvature*ny;
        Lst(4,7:9)= N*Curvature*nx + ny*dNds;
        
        Lst(5,4:6)=  - N*Curvature*nx - ny*dNds;
        Lst(5,7:9)=  nx*dNds - N*Curvature*ny;
        
        Dst=Ds;
        
        %         Dst=[Es*As,0,0,Es*As*ebar,0;
        %             0,Gs*An,0,0,Gs*An*ebar;
        %             0,0,Gs*Ab,0,0;
        %             Es*As*ebar,0,0,Es*In,0;
        %             0,Gs*As*ebar,0,0,Gs*Jt];
        Kelem=Kelem+NshapeMatrix'*Lst'*Dst'*Lst*NshapeMatrix*detJ*GaussWeights(ee);
        %         Kbeam(elementDof,elementDof)=Kbeam(elementDof,elementDof)+Kelem;
        %             NshapeMatrix'*Lst'*Dst*Lst*NshapeMatrix*detJ*GaussWeights(ee);
%         Kbend=Kbend+NshapeMatrix'*Lst'*Dst'*Lst*NshapeMatrix*detJ*GaussWeights(ee);
    end
    Kbeam(elementDof,elementDof)=Kbeam(elementDof,elementDof)+Kelem;
    StiffenerLength=StiffenerLength+2*detJ;
    disp('Tangential Angle is: ')
    alpha/pi*180
end
%% ---Output total stiffness matrix for beam ----
matrixNo=unique(DispFlagAll);
MatrixNo=[matrixNo matrixNo+FEM.NodeNumber  matrixNo+FEM.NodeNumber*2 ...
    matrixNo+FEM.NodeNumber*3 matrixNo+FEM.NodeNumber*4];
KK=Kbeam(MatrixNo,MatrixNo);
%% END OF DIARY