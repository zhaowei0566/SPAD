function [stress,strain,displacement,SdisplacementL]=...
    StressRecoveryStiffeners(FEM,Mat,Stiffener,Gausspointbend,ebar)
% ====================================================================
% recove the axial stress in the plate; Oct, 2015
% ====================================================================
% modified, Jan-3, 2016;
% Transform stiffener displacement in global
% coordinate system to local coordinate system and then calculate the axial
% stress in the stiffener for stiffener geometric stiffness;
%
Stiffener.thickness= Stiffener.width/Stiffener.layer*ones(1,Stiffener.layer);%% lamina along width
% stiffener displacement in global coordinate system, can be obtained by
% interpolation from surrounding plate node
Stiffener.displacement=zeros(Stiffener.GDof,size(Stiffener.TargetGaussian,2));

SdisplacementG=zeros(Stiffener.GDof,size(Stiffener.element,2));
% stiffener displacement in local coordinate system, obtained from
% transformation from global displacement;
SdisplacementL=zeros(Stiffener.GDof,size(Stiffener.element,2));

for stiffno=1:size(Stiffener.TargetGaussian,3)
    
    %% Calculate the nodal displacement for each stiffener;
    
    TargetPlate=Stiffener.TargetPlateID(:,stiffno);
    
    ShapeMatrixSP=Stiffener.ShapeMatrixSPTest(:,:,stiffno);
    
    for nodenum=1:Stiffener.nodenum
        
        plateElemID=TargetPlate(nodenum);
        
        plate_elem_node=FEM.elementNodes(plateElemID,:);
        
        plateElemDof1=[plate_elem_node plate_elem_node+FEM.NodeNumber plate_elem_node+2*FEM.NodeNumber...
            plate_elem_node+3*FEM.NodeNumber  plate_elem_node+4*FEM.NodeNumber];
        
        stiffenerDof=[nodenum nodenum+Stiffener.nodenum  nodenum+Stiffener.nodenum*2 ...
            nodenum+Stiffener.nodenum*3 nodenum+Stiffener.nodenum*4];
        
        
        SdisplacementG(stiffenerDof,stiffno)=...
            ShapeMatrixSP(stiffenerDof,plateElemDof1)*FEM.displacement(plateElemDof1)';
    end
    
    Stiffener.displacement=SdisplacementG;
    
    % == Need to transfer displacement in global coordinate system to ==
    % local coordinate system for axial stress, and then for stiffener
    % geometric stiffness. This is very important!!!
    
    SdisplacementL(1:Stiffener.GDof,stiffno)=Stiffener.Tmatrix(:,:,stiffno)*SdisplacementG(1:Stiffener.GDof,stiffno);
    
    %% For each laminate layer calculate stress
    
    for layer=1:length(Stiffener.thickness)
        disp('===================================')
        [AmatrixK,BmatrixK,DmatrixK,AshearK,Qreduced,Qbar,Q_bend_new,Q_shear_new]=...
            StiffenerCompositePerpendicular(Mat,Stiffener,layer);
        
        
        % perpendicular to the plate
        
        zcoordinate=ebar;
        
        for elem=1:Stiffener.elemnum

            
            stiffenerNODE=FEM.Stiffener.element(elem,2:end); %% Node NO. for one element
            
            XY= Stiffener.PointsCoord(stiffenerNODE,2:3);
            
            elementDof2=stiffenerNODE;
            
            elementDof=[elementDof2 elementDof2+Stiffener.nodenum  elementDof2+Stiffener.nodenum*2 ...
                elementDof2+Stiffener.nodenum*3 elementDof2+Stiffener.nodenum*4];
            
            NormalStrain=zeros(2,1); %
            Curvature0=zeros(2,1);
            ShearStrain=zeros(1,1);
            
            % Displacement indicies
            % each node has 5 DOFs
            % w, theta_x, theta_y, u, v
            
            
            [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointbend);
            
            %%  Loop for Gauss Points for in-plane stress \sigma_xx, \sigma_yy; and \sigma_{xy}
            for Gaussian_num=1:size(GaussWeights,1)
                % no need to sum them; strain at these Gaussian Points
                
                GaussPoint=GaussLocations(Gaussian_num,:);
                xi=GaussPoint(1);
                
                % shape function and derivatives
                [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);
                
                % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
                detJ=sqrt((naturalderivatives*XY(:,1))^2+(naturalderivatives*XY(:,2))^2);
                
                Curvature=StiffenerCurvature(XY,detJ, naturalderivatives, d2Nds2);% Curvature
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
                
                strainG(:,Gaussian_num)=Lst*SdisplacementL(elementDof)';
            end
            
            NormalStrain=sum(strainG(1:2,:),2)/Gaussian_num;
            
            Curvature0=sum(strainG(4:5,:),2)/Gaussian_num;
            ShearStrain=sum(strainG(3,:))/Gaussian_num;
            
            
            strain0=[NormalStrain(1)+zcoordinate*Curvature0(1),...
                NormalStrain(2)+zcoordinate*Curvature0(2)]';
            
            stress(elem,1:2,layer,stiffno)=Q_bend_new*strain0;
            strain(elem,1:2,layer,stiffno)=strain0';
        end
        
    end
end
displacement=Stiffener.displacement;