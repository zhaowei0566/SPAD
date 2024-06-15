function node_id = find_largest_difference_among_nodes(pntcords)


% pntcords = stiffener_cords';

node_id=[];

midpnt = pntcords((size(pntcords,1)+1)/2,:);
midpnt2 = pntcords((size(pntcords,1)+1)/2+1,:);

distance = norm(midpnt2 - midpnt);

for ii = 1:size(pntcords,1)-1
    
    cord1 = pntcords(ii,:);
    cord2 = pntcords(ii+1,:);
    
    dis = norm(cord2 - cord1);
    
    if dis/distance>1e2
        
        node_id = ii;
    end
    
    
    
end