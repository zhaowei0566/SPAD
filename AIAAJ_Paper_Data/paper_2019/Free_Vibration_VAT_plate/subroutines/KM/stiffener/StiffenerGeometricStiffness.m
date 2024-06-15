for stiffno=[1:2] % constistent with the number in the code: stiffener_Gauss_points;
    
    disp(['----Stiffener #' num2str(stiffno) '-----'])
    
    xstiffener=XXstiffener(stiffno,:);
    ystiffener=YYstiffener(stiffno,:);
    
    Stiffener.PointsCoord=zeros(stiffenernodesnumber,4);
    Stiffener.PointsCoord(:,1)=[1:stiffenernodesnumber]';
    Stiffener.PointsCoord(:,2)=[xstiffener]';
    Stiffener.PointsCoord(:,3)=[ystiffener]';
    
    Stiffener.element=zeros(stiffenerelementnumber,4);
    Stiffener.element(:,1)=[101:100+stiffenerelementnumber]';
    Stiffener.element(:,2)=[1:2:2*stiffenerelementnumber-1]';
    Stiffener.element(:,3)=[3:2:2*stiffenerelementnumber+1]';
    Stiffener.element(:,4)=[2:2:2*stiffenerelementnumber]';
    
    % specify gaussian points for each stiffener
    
    Stiffener.TargetPlate=TargetPlateID(stiffno,:);
%     eval(['Stiffener.TargetGaussian=TargetGaussianNum' num2str(stiffno) ';'])
    
    Stiffener.TargetGaussian=Stiffener.elementnode_naturalcoord(:,:,stiffno);
    
    
    %% ==========Laminates for stiffeners==========================
    %%%% ============Laminates for Stiffener
    %%%%========================================================
    display('----------Laminates for Stiffeners-----------');
    
    Stiffener.theta=PlyAngleStiffener;
    Stiffener.layer=length(Stiffener.theta);
    
    Stiffener.width=Stru.thickness;

    Stiffener.height=depthratio(depthNO)*Stiffener.width;
       
    Stiffener.CrossArea=Stiffener.width*Stiffener.height;
    Stiffener.Istiffener=1/12*Stiffener.width*Stiffener.height^3+ebar^2*Stiffener.CrossArea;
    %
    FEM.Stiffener=Stiffener;
    
    StiffenerElemNodeNum=size(Stiffener.element(1,:),2)-1;
    Stiffener.nodeID=unique(Stiffener.element(:,2:StiffenerElemNodeNum+1));
    Stiffener.nodecord= FEM.nodesCord(Stiffener.nodeID,:);
    
    %% Stiffener stiffness calculations
    
    %%%% Perpendicular or parallel laminate %%%%%%%%%%%%%%%
%     Stiffener.orientation='perpendicular'; % 1-parallel, 2-perpendicular
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     GJflag=1; %% 1- effective GJ; 0 - GJ=0; 2- originial GJ
    % ------------Calculate Constitutive Matrix [Ds] for stiffener----
    [Ds,Qbar,Stiffene]=LaminatedBeamDs(Stiffener,FEM,ebar,MatS,GJflag);
    
    %% ------------Stiffeness for each stiffener --------------
%     StiffenerGaussPoint='2';%%
    
%     [Kbeam,Kelem,StiffenerLength1]=...
%         LaminatedStiffenerStiffness5Dof(FEM,Ds,Stiffener,StiffenerGaussPoint);
    
    %--------- Transform local displacement fields into global ones ----
    
    [Tmatrix,alpha]=TransformationMatrix(FEM,Stiffener);
%     Kbeam=Tmatrix'*Kbeam*Tmatrix;
    
    %% -------------- stiffener {d} to plate {d} -------------------
    % [ShapeMatrixSP,ShapeMatrixElem]=PlateShapeStiffenerv2(FEM,Stiffener);
%     eval(['ShapeMatrixSP=ShapeMatrixSP' num2str(stiffno) ';']);
    
    ShapeMatrixSP=Stiffener.ShapeMatrixSPTest(:,:,stiffno);
    % ----------------- Sum stiffness -------------------------------
%     Kstiffener=Kstiffener+ShapeMatrixSP'*Kbeam*ShapeMatrixSP;
    
    %% ----------Geometric Stiffness of the stiffener -----------------
    % displacement compatability condition for stiffener axial stress
%     KGbeam_Local=NonGeometryStiffnessCurvedBeam(FEM,ebar,Stiffener,StiffenerGaussPoint,stiffno);
    % stress relatin for stiffener axial stress;
    
    KGbeam_Local=...
        NonGeometryStiffnessCurvedBeam_Interp(FEM,ebar,Stiffener,StiffenerGaussPoint,stiffno);
    
    KGbeam=Tmatrix'*KGbeam_Local*Tmatrix;
    
    KGStiffener=KGStiffener+ShapeMatrixSP'*KGbeam*ShapeMatrixSP;
        
end