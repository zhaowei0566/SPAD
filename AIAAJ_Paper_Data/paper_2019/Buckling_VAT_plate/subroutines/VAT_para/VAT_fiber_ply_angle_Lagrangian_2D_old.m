function theta = VAT_fiber_ply_angle_Lagrangian_2D_old(cords,VAT,layer)

% ctrlpnt_X = VAT.ctrlpnts_X;
% ctrlpnt_Y = VAT.ctrlpnts_Y;


ctrlpnt_X = VAT.ctrlpnts_X;

ctrlpnt_Y = VAT.ctrlpnts_Y;


ctrlT0T1 = (VAT.T0T1(:,:,layer));


X = cords(1);
Y = cords(2);

theta = 0;

for m = 1:length(ctrlpnt_X) % X-axis
    
    for n = 1:length(ctrlpnt_Y) % Y-axis

        Tmn = ctrlT0T1(n,m);
        

         LX = 1;
        
%         for m = 1:size(ctrl_ponts_theta,1)
            
            for ii = 1:length(ctrlpnt_X )
                
                if ii~=m
                    
                    temp_L = (X - ctrlpnt_X(ii))/(ctrlpnt_X(m) - ctrlpnt_X(ii));
                    
                    LX = LX*temp_L;
                    
                end
                
            end

        % 
         LY = 1;
        
%         for m = 1:size(ctrl_ponts_theta,1)
            
            
            for jj = 1:length(ctrlpnt_Y)
                
                if jj~=n
                    
                    temp_L = (Y - ctrlpnt_Y(jj))/(ctrlpnt_Y(n) - ctrlpnt_Y(jj));
                    
                    LY = LY*temp_L;
                    
                end
                
            end
            
        theta = theta + LX*LY*Tmn;
  
        
    end
end
