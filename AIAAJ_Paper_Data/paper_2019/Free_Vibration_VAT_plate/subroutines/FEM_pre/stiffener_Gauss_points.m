disp(['----Stiffener #' num2str(stiffnum) '-----'])

xxstiffener=XXstiffener(stiffnum,:);
yystiffener=YYstiffener(stiffnum,:);

Stiffener.PointsCoord=zeros(stiffenernodesnumber,4);
Stiffener.PointsCoord(:,1)=[1:stiffenernodesnumber]';
Stiffener.PointsCoord(:,2)=[xxstiffener]';
Stiffener.PointsCoord(:,3)=[yystiffener]';

Stiffener.element=zeros(stiffenerelementnumber,4);
Stiffener.element(:,1)=[101:100+stiffenerelementnumber]';
Stiffener.element(:,2)=[1:2:2*stiffenerelementnumber-1]';
Stiffener.element(:,3)=[3:2:2*stiffenerelementnumber+1]';
Stiffener.element(:,4)=[2:2:2*stiffenerelementnumber]';

[TargetGaussian,TargetPlateElement,Residual]=...
    StiffenerNodesGaussianPointsLocationsv2(FEM,Stiffener); %% To calculate the gaussian locations of the nodes at the stiffener in the plate located.
%     Stiffener.TargetPlate(stiffno,:)=TargetPlateElement;
%     Stiffener.TargetGaussian(stiffno,:)=TargetGaussian;

TargetPlateID(stiffnum,:)=TargetPlateElement;
Stiffener.TargetPlateID(:,stiffnum)=TargetPlateElement';


Stiffener.elementnode_naturalcoord(:,:,stiffnum)=TargetGaussian;

% eval(['TargetGaussianNum' num2str(stiffnum) '=TargetGaussian;'])

Stiffener.TargetGaussianTest(:,:,stiffnum)=TargetGaussian;

[ShapeMatrixSP,ShapeMatrixElem]=...
    PlateShapeStiffenerv3(FEM,Stiffener,TargetGaussian,TargetPlateElement);
%
%
% eval(['ShapeMatrixSP' num2str(stiffnum) '=ShapeMatrixSP;'])
Stiffener.ShapeMatrixSPTest(:,:,stiffnum)=ShapeMatrixSP;

