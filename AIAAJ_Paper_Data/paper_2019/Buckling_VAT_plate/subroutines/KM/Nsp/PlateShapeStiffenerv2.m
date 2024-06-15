function [ShapeMatrixSP,ShapeMatrixElem]=PlateShapeStiffenerv2(FEM,Stiffener)

disp(' -------- Calculate Shape Function for Stiffener ----------')


ShapeMatrixSP=zeros(Stiffener.GDof,FEM.GDof);

ShapeMatrixElem=zeros(Stiffener.nodenum,FEM.NodeNumber);


numberElements=Stiffener.nodenum;% number of beam elements for one stiffener

for ee=1:numberElements%% number of stiffener nodes
    %
%      disp(['-------------stiffener node No#' num2str(ee) '-------------']);

    
    TargetGaussian=Stiffener.TargetGaussian(ee,:); %% xi and eta
    plate=Stiffener.TargetPlate(ee);
    platenode=FEM.elementNodes(plate,:);
    xi=TargetGaussian(1);
    eta=TargetGaussian(2);
    NshapeElement(ee,platenode)=...
        [1/4*(1-xi)*(1-eta)*(-xi-eta-1);
        1/4*(1+xi)*(1-eta)*(xi-eta-1);
        1/4*(1+xi)*(1+eta)*(xi+eta-1);
        1/4*(1-xi)*(1+eta)*(-xi+eta-1);
        1/2*(1-eta)*(1+xi)*(1-xi);
        1/2*(1+xi)*(1+eta)*(1-eta);
        1/2*(1+eta)*(1+xi)*(1-xi);
        1/2*(1-xi)*(1+eta)*(1-eta)]';
    
    % w, thetaX, thetaY, u, v
    
    ShapeMatrixElem(ee,platenode)=ShapeMatrixElem(ee,platenode)+ ...
        NshapeElement(ee,platenode);
    
end

ShapeMatrixSP(1:Stiffener.nodenum,1:FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum:Stiffener.nodenum*2,1+FEM.NodeNumber:2*FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum*2:Stiffener.nodenum*3,1+FEM.NodeNumber*2:3*FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum*3:Stiffener.nodenum*4,1+FEM.NodeNumber*3:4*FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum*4:Stiffener.nodenum*5,1+FEM.NodeNumber*4:5*FEM.NodeNumber)=ShapeMatrixElem;