function [Kplate,LastPly]=LinearStiffnessLaminatedPlate_VAT_Lagrange_NLV(Mat, Stru,FEM,Laminate,center,VAT)
% elementype=1: CTRIA-3
% elementype=2: CTRIA-6
% elementype=3: CQUAD-4
% elementype=4: CQUAD-8
% nodedof=FEM.PlateNodeDof;

gDOF=FEM.GDof;

Kplate=zeros(gDOF,gDOF);

% Cycle for each element, element stiffness
for elem=1:FEM.numberElements
    %
    NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
    numberNodes=size(FEM.nodeCoordinates,1); %% how many nodes in FEM
    
    % Displacement indicies
    % each node has 5 DOFs
    % w, theta_x, theta_y, u, v
    elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes ...
        NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
    num_elem_nodes=length(NodeIndices); %% how many nodes for one element;
    
%     Kelem=zeros(num_elem_nodes*5,num_elem_nodes*5);%%% Element Stiffness Matrix;
    
    % ------------Bending Matrix-------------------
    [GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);
    % Loop for Gauss Points
    
    
    %%
    for ee=1:size(GaussWeights,1)
        
        
        GaussPoint=GaussLocations(ee,:);
        xi=GaussPoint(1);eta=GaussPoint(2);
        
        % shape function and derivatives
        [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
        % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
        [Jacob,invJacob,XYderivatives]=   JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
        %%
        
        Element_center_X   = sum(FEM.nodeCoordinates_label( NodeIndices,2))/length( NodeIndices);
        Element_center_Y   = sum(FEM.nodeCoordinates_label( NodeIndices,3))/length( NodeIndices);
        
%        
        
        Amatrix = zeros(3,3);
        Bmatrix = zeros(3,3);
        Dmatrix = zeros(3,3);
        Ashear  = zeros(2,2);
        
        for layer = 1:Laminate.layer
            
            %             T0 = T0T1(layer,1);
            %             T1 = T0T1(layer,2);
            
            %             theta = VAT_fiber_ply_angle_1D(T0,T1,Element_center_X,center,physical_length);
            
            
            cords = [abs(Element_center_X-center(1)) abs(Element_center_Y-center(2))];
            
            
            %             theta = VAT_fiber_ply_angle_Lagrangian(cords,T0T1(:,:,layer));
            
            
            theta = VAT_fiber_ply_angle_Lagrangian_2D(cords,VAT,layer);
            %%
            
            [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=LaminatedComposite_VAT(Mat,Stru,Laminate,layer,theta);
            
            Amatrix = Amatrix+AmatrixK;
            Bmatrix = Bmatrix+BmatrixK;
            Dmatrix = Dmatrix+DmatrixK;
            Ashear  = Ashear +AshearK;
            
            
        end
        
        %         Bmatrix
        %         Amatrix
        %         Dmatrix
        %% --------Membrane strain energy-------------
        Bstretch=zeros(3,num_elem_nodes*5);
        Bstretch(1,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(1,:);
        Bstretch(2,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(2,:);
        Bstretch(3,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(2,:);
        Bstretch(3,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(1,:);
        % B MATRIX FOR STRETCH
        
        Kplate(elementDof,elementDof)=Kplate(elementDof,elementDof)+...
            Bstretch'*Amatrix*Bstretch*GaussWeights(ee)*det(Jacob);%% Membrane strain energy
        
        %% -----Coupling extension-bending strain energy-------------
        Bbending=zeros(3,5*num_elem_nodes);
        Bbending(1,num_elem_nodes+1:num_elem_nodes*2)=XYderivatives(1,:);
        Bbending(2,num_elem_nodes*2+1:num_elem_nodes*3)=XYderivatives(2,:);
        Bbending(3,num_elem_nodes+1:num_elem_nodes*2)=XYderivatives(2,:);
        Bbending(3,num_elem_nodes*2+1:num_elem_nodes*3)=XYderivatives(1,:);
        %
        Kplate(elementDof,elementDof)=Kplate(elementDof,elementDof)+...
            Bstretch'*Bmatrix*Bbending*GaussWeights(ee)*det(Jacob);
        Kplate(elementDof,elementDof)=Kplate(elementDof,elementDof)+...
            Bbending'*Bmatrix*Bstretch*GaussWeights(ee)*det(Jacob);
        %% ---------Bending strain energy-----------------------
        Kplate(elementDof,elementDof)=Kplate(elementDof,elementDof)+...
            Bbending'*Dmatrix*Bbending*GaussWeights(ee)*det(Jacob);
    end
    
    %% --------------- For shear matrix----------------------
    [GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointShear);
    
    for ee=1:size(GaussWeights,1)
        
        
        GaussPoint=GaussLocations(ee,:);
        xi=GaussPoint(1);eta=GaussPoint(2);
        
        % shape function and derivatives
        [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
        
        % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
        [Jacob,invJacob,XYderivatives]=  JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
        %            
                        
        Amatrix = zeros(3,3);
        Bmatrix = zeros(3,3);
        Dmatrix = zeros(3,3);
        Ashear =zeros(2,2);
        
        for layer = 1:Laminate.layer
            
            %             T0 = T0T1(layer,1);
            %             T1 = T0T1(layer,2);
            %
            %             theta = VAT_fiber_ply_angle_1D(T0,T1,Element_center_X,center,physical_length);
            
            
            cords = [abs(Element_center_X-center(1)) abs(Element_center_Y-center(2))];
            
            theta = VAT_fiber_ply_angle_Lagrangian_2D(cords,VAT,layer);
            %             theta = VAT_fiber_ply_angle_Lagrangian(cords,T0T1(:,:,layer));
            
            LastPly(elem,:,layer) = [elem Element_center_X Element_center_Y theta];
            
            [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=LaminatedComposite_VAT(Mat,Stru,Laminate,layer,theta);
            
            Amatrix = Amatrix+AmatrixK;
            Bmatrix = Bmatrix+BmatrixK;
            Dmatrix = Dmatrix+DmatrixK;
            Ashear  = Ashear +AshearK;
            
            
        end
        
        % [B] matrix for bending
        Bshear=zeros(2,5*num_elem_nodes);
        Bshear(1,1:num_elem_nodes)=XYderivatives(1,:);
        Bshear(1,1+num_elem_nodes:num_elem_nodes*2)=shape;
        Bshear(2,1:num_elem_nodes)=XYderivatives(2,:);
        Bshear(2,1+num_elem_nodes*2:num_elem_nodes*3)=shape;
        %
        Kplate(elementDof,elementDof)=Kplate(elementDof,elementDof)+ Bshear'*Ashear*Bshear*GaussWeights(ee)*det(Jacob);
        
%         Kelem=Kelem+Kplate(elementDof,elementDof);
    end
    
%     det(Kelem);
    
    
end


