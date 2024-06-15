function [ActiveDof,Constrained_Dof]=InPlaneBCPlate5Dof(FEM)

type=FEM.BCtype;

%% FIND the SidesDof through Code for Square plate or Rectangle Plate
xx=FEM.nodeCoordinates(:,1);
yy=FEM.nodeCoordinates(:,2);

nodeNum=size(FEM.nodeCoordinates,1);

%% Offer ActiveDOF manually
% w, betax, betay, u, v
%%

SideNodesNumList_LHS = find(xx==min(FEM.nodeCoordinates(:,1)))';
figure(200);hold on; plot(xx(SideNodesNumList_LHS ),yy(SideNodesNumList_LHS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
SideNodesNumList_RHS = find(xx>=0.99*max(FEM.nodeCoordinates(:,1)))';
figure(200);hold on; plot(xx(SideNodesNumList_RHS ),yy(SideNodesNumList_RHS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
SideNodesNumList_THS = find(yy>=0.99*max(FEM.nodeCoordinates(:,2)))';
figure(200);hold on; plot(xx(SideNodesNumList_THS ),yy(SideNodesNumList_THS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
SideNodesNumList_BHS = find(yy==min(FEM.nodeCoordinates(:,2)))';
figure(200);hold on; plot(xx(SideNodesNumList_BHS ),yy(SideNodesNumList_BHS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);

switch type
    
    % Order of degrees of freedom: 
    % w,theta_x,theta_y, u, v
    
    case 'SFSF'
        
        % left side
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        
        SidesDof_LHS = [];
        
        SidesDof_LHS=[SidesDof_LHS SideNodesNumList];%% w is constrainted
        
        % right side
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        
        % top side
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % bottom side
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        %          Constrained_Dof= unique([SidesDof_LHS  SidesDof_THS SidesDof_BHS]);
        
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        

end

