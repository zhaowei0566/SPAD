function Stiffener = Stiffener_perpendicular_stress_recover(FEM,Stiffener)


% compute the axial stress along stiffener

% use displacement in local coordinate system


disp ( ' ======== S T I F F E N E R    S T R E S S      R E C O V E R Y ==========');



displacements_local = Stiffener.stiffeners_displacement_Local;

for stiffno = 1:Stiffener.stiffeners_num
    
    
    displacements_local_each = displacements_local(:,stiffno);
    
    for layer = 1:length(Stiffener.theta)
        
        for elem=1:Stiffener.elemnum
            
            % strain at median plane
            
            
            stiffenerNODE=Stiffener.element(elem,2:end);
            
            
            XY= Stiffener.PointsCoord(stiffenerNODE,2:3);
            
            elementDof=[elementDof2 elementDof2+Stiffener.nodenum  elementDof2+Stiffener.nodenum*2 ...
                           elementDof2+Stiffener.nodenum*3 elementDof2+Stiffener.nodenum*4];
            
            
            
            
            NormalStrain=zeros(2,1); %
            Curvature0=zeros(2,1);
            ShearStrain=zeros(1,1);
            
            % Displacement indicies
            % each node has 5 DOFs
            % w, theta_x, theta_y, u, v
            
            
            [GaussWeights,GaussLocations]=gaussQuadrature1D( Stiffener.StiffenerGaussPoint);
            
            
            
            
            for Gaussian_num=1:size(GaussWeights,1)
                
                
                GaussPoint=GaussLocations(Gaussian_num,:);
                xi=GaussPoint(1);
                
                [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,Stiffener);
                
                
                detJ=sqrt((naturalderivatives*XY(:,1))^2+(naturalderivatives*XY(:,2))^2);
                
                Curvature=StiffenerCurvature(XY,detJ, naturalderivatives, d2Nds2);% Curvature
                
                
                %% ------------------------------------------------------
                Lst=zeros(5,15);%% Ali's paper
                dNds=naturalderivatives/detJ;
                N=shape';
                
                Lst(1,10:12)= dNds;
                Lst(1,13:15)= N*Curvature;
                
                Lst(2,10:12)= - N*Curvature;
                Lst(2,13:15)= dNds;
                
                Lst(3,1:3) =  dNds;
                Lst(3,4:6) = N;
                
                Lst(4,4:6) = dNds;
                Lst(4,7:9) = N*Curvature;
                
                
                Lst(5,4:6)=  - N*Curvature;
                Lst(5,7:9)=  dNds ;
                
                strain(:,Gaussian_num)=Lst* displacements_local_each(elementDof)';
                
                
                
                
            end
            
            
            
            
            
            
            
        end
    end
end