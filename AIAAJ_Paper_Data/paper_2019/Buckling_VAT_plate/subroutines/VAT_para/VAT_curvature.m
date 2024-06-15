function curvature = VAT_curvature(FEM,VAT_type,VAT)



curvature  =zeros(size(FEM.nodeCoordinates_label,1),1);

switch VAT_type
    
    
    
    
    case 'LV-X'
        
        
        for pnt_num = 1:size(FEM.nodeCoordinates_label,1)
            
            
            
            
            x = FEM.nodeCoordinates_label(pnt_num,2);
            
            % only check first layer fiber path
            
            T0 = VAT.T0T1(1,1);
            T1 = VAT.T0T1(1,2);
            center = VAT.center;
            physical_length = VAT.physical_length;
            
            theta = VAT_fiber_ply_angle_1D(T0,T1,x,center(1),physical_length);
            
            theta_4_grad = VAT_fiber_ply_angle_1D(T0,T1,x+1e-5,center(1),physical_length);
            
            
            
            grad_theta = (theta_4_grad/180*pi - theta/180*pi)/1e-5;
            
            
            curvature(pnt_num) = cosd(theta)*grad_theta;
            
            %             approximate_curvature = 1/cosd(theta)^2* grad_theta
            
            
        end
        
        
        
        
    case 'LV-Y'
        
        
        for pnt_num = 1:size(FEM.nodeCoordinates_label,1)
            
            
            
            
            y = FEM.nodeCoordinates_label(pnt_num,3);
            
            % only check first layer fiber path
            
            T0 = VAT.T0T1(1,1);
            T1 = VAT.T0T1(1,2);
            center = VAT.center;
            physical_length = VAT.physical_length;
            
            theta = VAT_fiber_ply_angle_1D(T0,T1,y,center(2),physical_length);
            
            theta_4_grad = VAT_fiber_ply_angle_1D(T0,T1,y+1e-5,center(2),physical_length);
            
            
            
            grad_theta = (theta_4_grad/180*pi - theta/180*pi)/1e-5;
            
            
            curvature(pnt_num) = cosd(theta)*grad_theta;
            
            %             approximate_curvature = 1/cosd(theta)^2* grad_theta
            
            
        end
        
        
        
        
    case 'LV-XY'
        
        
        for pnt_num = 1:size(FEM.nodeCoordinates_label,1)
            
            phi = VAT.rigid_rotation_angle;
            
            x = FEM.nodeCoordinates_label(pnt_num,2);
            y = FEM.nodeCoordinates_label(pnt_num,3);
            
            % only check first layer fiber path % rotate does not change
            % the curvature
            
            T0 = VAT.T0T1(1,1);
            T1 = VAT.T0T1(1,2);
            center = VAT.center;
            physical_length = VAT.physical_length;
            
            %             theta = VAT_fiber_ply_angle_1D(T0,T1,y,center(2),physical_length);
            
            theta = VAT_fiber_ply_angle_1D_rotate(T0,T1,x,y,center(1),physical_length,phi);
            
            theta_4_grad = VAT_fiber_ply_angle_1D_rotate(T0,T1,x+1e-5,y,center(1),physical_length,phi);
            
            
            
            grad_theta = (theta_4_grad/180*pi - theta/180*pi)/1e-5;
            
            
            curvature(pnt_num) = cosd(theta)*grad_theta;
            
            %             approximate_curvature = 1/cosd(theta)^2* grad_theta
            
        end
        
        
    case 'NLV'
        
        
        for pnt_num = 1:size(FEM.nodeCoordinates_label,1)
            
            
            x = FEM.nodeCoordinates_label(pnt_num,2);
            y = FEM.nodeCoordinates_label(pnt_num,3);
            
            cordtemp =  [abs(x-VAT.center(1)) abs(y-VAT.center(2)) ] ;
            
            theta = VAT_fiber_ply_angle_Lagrangian_2D(cordtemp ,VAT,1);
            
            
            
            x = FEM.nodeCoordinates_label(pnt_num,2)+1e-5;
            
            
            cordtemp =  [abs(x-VAT.center(1)) abs(y-VAT.center(2)) ] ;
            
            theta_4_grad = VAT_fiber_ply_angle_Lagrangian_2D(cordtemp ,VAT,1);
            
            
            grad_theta = (theta_4_grad/180*pi - theta/180*pi)/1e-5;
            
            
            curvature(pnt_num) = cosd(theta)*grad_theta;
            
            
            
        end
        
        
end