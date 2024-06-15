function plateID  = plate_ID_stiffener_node(FEM,stiffener_node_cord)

% find out the plate ID in which the stiffener node located.

% works for 8-noded shell element


% element area


for stiff_node = 1:size(stiffener_node_cord,1)
    
    
    stiffener_pnt_temp = stiffener_node_cord(stiff_node,:);
    
    
    flag = 1;
    elem = 1;
    
    while elem <= FEM.elementNumber && flag
        
        
        element_nodes = FEM.elementNodes(elem,:);
        
        nodes_cord_temp  = zeros(4,4);
        
        % first 4 nodes to generate area
        for ii = 1:4
            nodes_cord_temp(ii,:) = label2cord(element_nodes(ii),FEM.nodeCoordinates_label);
            
        end

        total_area = polyarea(nodes_cord_temp(:,2),nodes_cord_temp(:,3));
        
        %%
        
        nodes_cord_temp(ii+1,:) = nodes_cord_temp(1,:);
        part_area = 0;
        for jj = 1:4
            
           XX = [ stiffener_pnt_temp(1) nodes_cord_temp(jj,2)  nodes_cord_temp(jj+1,2)] ;
            
             YY = [ stiffener_pnt_temp(2) nodes_cord_temp(jj,3)  nodes_cord_temp(jj+1,3)]; 
             
             
            part_area =  part_area + polyarea(XX,YY); 
        end
        
        
        %% compare 
        
        if abs((total_area-part_area)/total_area)<=1e-5
            flag = 0;
            
            plateID(stiff_node) = elem;
        else
            elem = elem + 1;
        end
        
    end
end