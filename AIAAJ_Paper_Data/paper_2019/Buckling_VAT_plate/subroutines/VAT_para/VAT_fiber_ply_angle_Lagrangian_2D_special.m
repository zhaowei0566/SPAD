function theta = VAT_fiber_ply_angle_Lagrangian_2D_special(cords,VAT,layer)

% the case for Wu design 3



T0T1_THETA1 = [   71.0000   49.5000   71.5000
                  67.0000   50.0000   51.0000
                  17.0000   12.0000   45.0000];

Theta1(:,:,1) = T0T1_THETA1 ;
Theta1(:,:,2) = -T0T1_THETA1 ;
Theta1(:,:,7) = -T0T1_THETA1 ;
Theta1(:,:,8) = T0T1_THETA1 ;

T0T1_THETA2 = -[89 67.5 64 65 82
    81 69 65 60 60
    80.5 66.5 58 54.5 59
    26 25 18 24 25
    8 -4.5 1 -5 0];

Theta2(:,:,3) = T0T1_THETA2 ;
Theta2(:,:,4) = -T0T1_THETA2 ;
Theta2(:,:,5) = -T0T1_THETA2 ;
Theta2(:,:,6) = T0T1_THETA2 ;

xi = cords(1);
yi = cords(2);

if layer>=3 && layer<=6
    
    
    
    
    
    ctrlT0T1 = Theta2(:,:,layer);
else
    
    
    
    
    ctrlT0T1 = Theta1(:,:,layer);
end

zd = ctrlT0T1;


ctrlpnt_X = linspace(0,0.127,size(zd,1));%(VAT.ctrlpnts_X);

ctrlpnt_Y = linspace(0,0.127,size(zd,2));%( VAT.ctrlpnts_Y);

%     ctrlT0T1 = VAT.T0T1(:,:,layer);




mx = length(ctrlpnt_X)-1;
my = length(ctrlpnt_Y)-1;

xd_1d = ctrlpnt_X ;
yd_1d = ctrlpnt_Y ;

%     zd = ctrlT0T1;

ni = 1;



theta = lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi );



