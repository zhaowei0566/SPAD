%% A quick tool to check the fiber path

layer = 1

[centerX,centerY]  = meshgrid(linspace(0,Plate.length,1001),linspace(0,Plate.width,1001));


if exist('VAT')
    
    flag = 'NLV';
else
    flag = 'LV';
end

flag = VAT.VAT_type

% [centerX,centerY]  = meshgrid(linspace(-0.5,0.5,101),linspace(-0.5,0.5,101));
for itemp = 1:size(centerX,1)
    
    for jtemp = 1:size(centerY,2)
        
        cordtemp =  [abs(centerX(itemp,jtemp)-center(1)) abs(centerY(itemp,jtemp)-center(2)) ] ;
        
        
        switch flag
            
            case 'NLV'
                
                theta = VAT_fiber_ply_angle_Lagrangian_2D( cordtemp ,VAT,layer);
             
                
            case  'LV'
                T0 = T0T1(layer,1);
                T1 = T0T1(layer,2);
                
                theta = VAT_fiber_ply_angle_1D(T0,T1,centerX(itemp,jtemp),center(1),physical_length);
                
              case  'LV-Y'
                T0 =  T0T1(layer,1);
                T1 = T0T1(layer,2);
                
                theta = VAT_fiber_ply_angle_1D(T0,T1,centerY(itemp,jtemp),center(2),physical_length);
                
            case 'LV-XY'

                T0 = T0T1(layer,1);
                T1 =  T0T1(layer,2);
                phi = VAT.rigid_rotational_angle;
                theta = VAT_fiber_ply_angle_1D_rotate(T0,T1,centerX(itemp,jtemp),centerY(itemp,jtemp),center(1),physical_length,phi);
                
                
        end
        
        
        test_Theta(itemp,jtemp) = theta;
        
        u(itemp,jtemp) = cosd(theta);
        v(itemp,jtemp) = sind(theta);
        
    end
    
end



vec_p = [1 0 0 0 ;
    2 Plate.length 0 0;
    3 Plate.length Plate.width 0;
    4 0 Plate.width 0];



vec_t = [1 2 3 4];
patch_plot(vec_t,vec_p,34,'skin')
figure(34);hold on;streamslice(centerX,centerY,u,v,2,'noarrows','nearest');axis image;box on;
set(gcf,'color','w')
axis([0 Plate.length 0 Plate.width]);
set(gca,'FontSize',16);
set(gcf,'color','w');

figure(34);hold on; line([0 Plate.length],[Plate.width/2 Plate.width/2])
figure(34);hold on; line([Plate.length/2 Plate.length/2],[Plate.width/2*0 Plate.width])

% 
% num = 20;
% figure;
% % startx = linspace(0,Plate.length,num);
% % starty = ones(1,length(startx))*Plate.width;
% % streamline(centerX,centerY,u,v,startx,starty);axis image;box on;hold on;
% % streamline(centerX,centerY,-u,v,startx,starty);axis image;box on;hold on;
% 
% 
% startx = linspace(0,Plate.length,num);
% starty = ones(1,num)*Plate.width*0;
% 
% streamline(centerX,centerY,u,v,startx,starty);axis image;box on;
% 
% starty = linspace(0,Plate.length,num);
% startx = ones(1,num)*Plate.width*0;
% 
% streamline(centerX,centerY,u,v,startx,starty);axis image;box on;
% 
% 
% % hold on;
% % streamline(centerX,centerY,-u,v,startx,starty);axis image;box on;
% 
% figure;
% num = 40;
% startx = linspace(0.,Plate.length,num);
% starty = 0.*ones(num,1);
% % starty = linspace(0,Plate.length,num);
% 
% streamline(centerX,centerY,u,v,startx,starty);axis image;box on;
% 
% num = 20;
% startx = linspace(0.,Plate.length,num)*0;
% starty = linspace(0.,Plate.length,num);
% % starty = linspace(0,Plate.length,num);
% 
% streamline(centerX,centerY,u,v,startx,starty);axis image;box on;

