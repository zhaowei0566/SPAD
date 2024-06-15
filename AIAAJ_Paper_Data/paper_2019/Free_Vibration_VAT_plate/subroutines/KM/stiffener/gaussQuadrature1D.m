function [GaussWeights,GaussLocations]=gaussQuadrature1D(points)

switch points

    case '1'
        
        GaussLocations=[0];
        GaussWeights=[2];
        
    case '2'
        GaussLocations=[-0.577350269189626;0.577350269189626];
        GaussWeights=ones(2,1);   
    case '3'
        
        GaussLocations=[-0.774596669241483;0.774596669241483;
            0];
        GaussWeights=[0.555555555555556;0.555555555555556;...
            0.888888888888889];
end