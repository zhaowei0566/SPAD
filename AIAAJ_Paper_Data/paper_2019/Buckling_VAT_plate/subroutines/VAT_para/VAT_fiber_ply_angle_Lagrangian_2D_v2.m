function theta = VAT_fiber_ply_angle_Lagrangian_2D_v2(cords,VAT,layer)

ctrlpnt_X = (VAT.ctrlpnts_X);

ctrlpnt_Y = fliplr(VAT.ctrlpnts_Y);

ctrlT0T1 = VAT.T0T1(:,:,layer);


xi = cords(1);
yi = cords(2);


mx = length(ctrlpnt_X)-1;
my = length(ctrlpnt_Y)-1;

xd_1d = ctrlpnt_X ;
yd_1d = ctrlpnt_Y ;

zd = ctrlT0T1;

ni = 1;

theta = lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi );


