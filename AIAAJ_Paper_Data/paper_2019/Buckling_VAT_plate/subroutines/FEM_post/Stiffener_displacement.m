function Stiffener = Stiffener_displacement(FEM,Stiffener)

% displacement in both global and local coordinate system for stiffener
% nodes
% v0; works for one stiffener now


stiffeners_displacement_Global = zeros(Stiffener.GDof,Stiffener.stiffeners_num);


stiffeners_displacement_Local  = zeros(Stiffener.GDof,Stiffener.stiffeners_num);



for stiffno = 1:Stiffener.stiffeners_num
    
    
    % belong to which set of plate IDs
    
    TargetPlate=Stiffener.TargetPlateID(:,stiffno);
    
    % and the corresponding shape functions?
    ShapeMatrixSP=Stiffener.ShapeMatrixSPTest(:,:,stiffno);
    
    
    
    
    for nodenum=1:Stiffener.nodenum
        
        plateElemID=TargetPlate(nodenum);
        
        plate_elem_node=FEM.elementNodes(plateElemID,:);
        
        plateElemDof1=[plate_elem_node plate_elem_node+FEM.NodeNumber plate_elem_node+2*FEM.NodeNumber...
            plate_elem_node+3*FEM.NodeNumber  plate_elem_node+4*FEM.NodeNumber];
        
        stiffenerDof=[nodenum nodenum+Stiffener.nodenum  nodenum+Stiffener.nodenum*2 ...
            nodenum+Stiffener.nodenum*3 nodenum+Stiffener.nodenum*4];
        
        
        stiffeners_displacement_Global(stiffenerDof,stiffno)=...
            ShapeMatrixSP(stiffenerDof,plateElemDof1)*FEM.displacement(plateElemDof1)';
    end
    
    
    
    % Transformation matrix
    % global displacement to local displacement
    
    
    Tmatrix = Stiffener.Tmatrix(:,:,stiffno);
    
    
    stiffeners_displacement_Local(:,stiffno) = Tmatrix*stiffeners_displacement_Global(:,stiffno);
    

end


Stiffener.stiffeners_displacement_Global = stiffeners_displacement_Global;

Stiffener.stiffeners_displacement_Local = stiffeners_displacement_Local;
