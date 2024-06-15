function theta = VAT_fiber_ply_angle_1D_rotate(T0,T1,x,y,center,width,phi)
% all angles are in degrees

% use coordinate to identify the fiber ply orientation
%
%

d= width/2;

xo = center(1); % origin of x coordinate


% % xbar = (x-xo)/d;

% % if xbar>=0
% %     
% %     theta = T0 + (x-xo)/d*(T1-T0);
% %     
% % else
% %     
% %     
% %     theta = T0 - (x-xo)/d*(T1-T0);
% % end
% % 
% % 
% % 


xprime = cosd(phi)*x + sind(phi)*y;
xbar_prime = (xprime - xo)/d;




if xbar_prime>=0
    
    theta = T0 + phi+(xprime - xo)/d*(T1-T0);
    
else
    
    
    theta = T0  + phi-(xprime - xo)/d*(T1-T0);
end








