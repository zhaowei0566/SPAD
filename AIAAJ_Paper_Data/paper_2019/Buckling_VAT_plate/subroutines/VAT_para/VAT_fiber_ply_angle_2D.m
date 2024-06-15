function theta = VAT_fiber_ply_angle_2D(T0,T1,Phi,centerXY,center,physical_length)


% not working at this time.

% use coordinate to identify the fiber ply orientation
%
%




DeltaX = centerXY(1); - center(1);
DeltaY = centerXY(2);  - center(2);


xprime = DeltaX * cosd(Phi) + DeltaY * sind(Phi);


physical_length = sqrt(DeltaX^2 + DeltaY^2);
xbar = (xprime-center(1))/(physical_length);

% 
if xbar>=0
    
    theta = T0 + xbar*(T1-T0) ;
    
else
    
    
    theta = T0 - xbar*(T1-T0) ;
end



% theta = theta+90;




