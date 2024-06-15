function theta = VAT_fiber_ply_angle_Lagrangian_v2(cords,ctrl_pnts_theta)

% use coordinate to identify the fiber ply orientation
%
%



Y = cords;

theta = 0;

% ctrl_ponts_theta1 = [0 68;
%     1/2*plate_width/2 55;
%     plate_width/2 19];

% LY = ones(size(ctrl_ponts_theta ,1),1);

for m = 1:size(ctrl_pnts_theta,1)
    
    LY=1;
    
    for ii = 1:size(ctrl_pnts_theta,1)
        
        if ii~=m
            temp_L = (Y - ctrl_pnts_theta(ii,1))/(ctrl_pnts_theta(m,1) - ctrl_pnts_theta(ii,1));
            
            LY = LY*temp_L;
            
        end
        
    end
    
    theta = theta+ctrl_pnts_theta(m,2)*LY;
end


% theta


theta =lagrange(Y,ctrl_pnts_theta(:,1)',ctrl_pnts_theta(:,2)');
