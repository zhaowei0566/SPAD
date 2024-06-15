function [KG]=GeometryStiffnessCurvedBeam(FEM,ebar,Stiffener, Gausspointbend)
% This function if to calculate the geometry stiffness of the curved
% stiffener
% Modified June 9, 2014
% Modified Aug, 9, 2014
% Modify axial stress for stiffener, May 23, 2015
disp(['-------------Assembling Stiffener Geometric Stiffness-------------']);
gDOF=Stiffener.nodedof*Stiffener.nodenum;
KG=zeros(gDOF,gDOF);
FEM.typeplate='CBAR3';
As=Stiffener.width*Stiffener.height;
% Istiffener11=1/12*Stiffener.width*Stiffener.height^3+ebar^2*As;
Istiffener11=Stiffener.Istiffener;
numberElements=size(FEM.Stiffener.element,1);
%-----Geometry Stiffness due to Tranverse deformation------
% Cycle for each element, element stiffness,
for e=1:numberElements
    %
    stiffenerNODE=FEM.Stiffener.element(e,2:end); %% Node NO. for one element
    XY= Stiffener.PointsCoord(stiffenerNODE,2:3);
    
    
    elementDof2=stiffenerNODE;
    
    elementDof=[elementDof2 elementDof2+Stiffener.nodenum  elementDof2+Stiffener.nodenum*2];
    %
    [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointbend);
    
    % Loop for Gauss Points
    for ee=1:size(GaussWeights,1)
        GaussPoint=GaussLocations(ee,:);
        xi=GaussPoint(1);
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);
        
        detJ=sqrt((naturalderivatives*XY(:,1))^2+(naturalderivatives*XY(:,2))^2);
        
        %% 
        if naturalderivatives*XY(:,1)>=0 && naturalderivatives*XY(:,2)>=0
            
            alpha=atan( (naturalderivatives*XY(:,2)) / (naturalderivatives*XY(:,1)) );
            
        elseif naturalderivatives*XY(:,1)<=0 && naturalderivatives*XY(:,2)>=0
            
            alpha=pi+atan( (naturalderivatives*XY(:,2)) / (naturalderivatives*XY(:,1)) );
            
        elseif naturalderivatives*XY(:,1)<0 && naturalderivatives*XY(:,2)<0
            
            alpha=-pi+atan( (naturalderivatives*XY(:,2)) / (naturalderivatives*XY(:,1)) );
            
        elseif naturalderivatives*XY(:,1)>=0 && naturalderivatives*XY(:,2)<=0
            
            alpha=atan( (naturalderivatives*XY(:,2)) / (naturalderivatives*XY(:,1)) );
            
        end
        
%         alpha/pi*180
        %%
        
        stress=(FEM.stress0(1,1)+FEM.stress0(2,2))/2+...
            (FEM.stress0(1,1)-FEM.stress0(2,2))/2*cos(2*alpha)+...
            FEM.stress0(1,2)*sin(2*alpha);
        
        stress00=[As*stress,0,0;
            0,(Istiffener11)*stress,0;
            0,0,(Istiffener11)*stress];
        
        Curvature=StiffenerCurvature(XY,detJ, naturalderivatives, d2Nds2);
        dNds=naturalderivatives/detJ;
        Lstiffener=zeros(3,9);
        
        Lstiffener(1,1:3)=dNds;
        
        Lstiffener(2,4:6)=dNds;
        Lstiffener(2,7:9)=Curvature*shape';
        
        Lstiffener(3,4:6)=-Curvature*shape';
        Lstiffener(3,7:9)=dNds;
        %         Lsbending=zeros(2,8);
        %         Lsbending(1,1:8)=XYderivatives(1,:);Lsbending(2,1:8)=XYderivatives(2,:);
        KG(elementDof,elementDof)=KG(elementDof,elementDof)+...
            Lstiffener'*stress00*Lstiffener*GaussWeights(ee)*detJ;
    end
end
