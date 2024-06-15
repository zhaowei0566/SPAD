function [TargetGaussian,TargetPlateElement,Residual]=StiffenerNodesGaussianPointsLocationsv2(FEM,Stiffener)
% Gaussian location of nodes on stiffeners %
% find the nodes in which plate element and find out Gaussian Points
% Location
% Ali-1stiffener
% Ill-condition in solving the eqaution
% clc;
format long;
PointsCoord=Stiffener.PointsCoord*1;
disp('---------Computing Natural Coordinates of the stiffener nodes-----------');
%% Loop run to find the nodes where they located
TargetGaussian=zeros(size(PointsCoord,1),2);

for stiffenerID=1:size(PointsCoord,1)
    %     disp(['-----Stiffener Node #' num2str(stiffenerID)  '----------']);
    stiffenerX=PointsCoord(stiffenerID,2);
    stiffenerY=PointsCoord(stiffenerID,3);
    % Initialize some data for searching plate element
    plateID=1;
    flag=1;
    while plateID<=size(FEM.elementNodes,1) && flag %% Loop
        PlateNodeID=FEM.elementNodes(plateID,:); %% Nodes ID in the plate elements
        plateX=FEM.nodesCord(PlateNodeID,2)*1;
        plateY=FEM.nodesCord(PlateNodeID,3)*1;
        
        % Judge wheather the points in this plate element
        flagX=stiffenerX<=max( plateX) && stiffenerX>=min( plateX);
        flagY=stiffenerY<=max( plateY) && stiffenerY>=min( plateY);
        
        
        
        
        
        if flagX && flagY
            
            flag=0;
            syms xi eta real;
            shape=[
                1/4*(1-xi)*(1-eta)*(-xi-eta-1);
                1/4*(1+xi)*(1-eta)*(xi-eta-1);
                1/4*(1+xi)*(1+eta)*(xi+eta-1);
                1/4*(1-xi)*(1+eta)*(-xi+eta-1);
                1/2*(1-eta)*(1+xi)*(1-xi);
                1/2*(1+xi)*(1+eta)*(1-eta);
                1/2*(1+eta)*(1+xi)*(1-xi);
                1/2*(1-xi)*(1+eta)*(1-eta)];
            %             disp('----fun1 and fun2-----')
            fun1=vpa(shape'*plateX-stiffenerX,4);
            fun2=vpa(shape'*plateY-stiffenerY,4);
            
            %                         fun1=simplify(shape'*plateX-stiffenerX,4)
            %             fun2=shape'*plateY-stiffenerY,4
            
            solution=vpasolve([fun1==0,fun2==0],[xi,eta]);
            solution.xi;
            solution.eta;
            right_xi_id=find(round(solution.xi)<=1&round(solution.xi)>=-1);
            right_eta_id=find(round(solution.eta)<=1&round(solution.eta)>=-1);
            right_id=intersect(right_xi_id,right_eta_id);
            eta=solution.eta(right_id);
            xi=solution.xi(right_id);
            
            %             eta=solve(fun2==0,eta>=-1,eta<=1);
            %             eta=vpa(eta,4);
            %             xi=solve(fun1==0,xi>=-1, xi<=1);
            %             xi=vpa(xi,4);
            
            shape=[1/4*(1-xi)*(1-eta)*(-xi-eta-1);
                1/4*(1+xi)*(1-eta)*(xi-eta-1);
                1/4*(1+xi)*(1+eta)*(xi+eta-1);
                1/4*(1-xi)*(1+eta)*(-xi+eta-1);
                1/2*(1-eta)*(1+xi)*(1-xi);
                1/2*(1+xi)*(1+eta)*(1-eta);
                1/2*(1+eta)*(1+xi)*(1-xi);
                1/2*(1-xi)*(1+eta)*(1-eta)];
            
            residualX=(shape'*plateX-stiffenerX);
            residualY=(shape'*plateY-stiffenerY);
            % %             eval([sol_xi,sol_eta])
            TargetPlateElement(stiffenerID)=plateID;
            %             TargetGaussian(stiffenerID,:)=([xi,eta]);
            
            
            TargetGaussian(stiffenerID,1)=xi;
            TargetGaussian(stiffenerID,2)=eta;
            Residual(stiffenerID,:)=[residualX,residualY];
            
        end
        plateID=plateID+1;
    end
end

