function [ShapeMatrixSP,ShapeMatrixElem]=PlateShapeStiffener(FEM,Stiffener)



ShapeMatrixSP=zeros(Stiffener.GDof,FEM.GDof);

ShapeMatrixElem=zeros(Stiffener.nodenum,FEM.NodeNumber);


numberElements=size(FEM.Stiffener.element,1);% number of beam elements for one stiffener

for e=1:numberElements%% number of stiffener element
    %
    disp(['-------------stiffener No#' num2str(e) '-------------']);
    
    stiffenerNODE=FEM.Stiffener.element(e,2:end); %% Node NO. for one element
    
    NshapeElement=zeros(3,FEM.GDof/5);
    for eee=1:3 %% Three nodes in one beam element 5 7 6
        nodestifflag=stiffenerNODE(eee);
        TargetPlate=Stiffener.TargetPlate(nodestifflag-4); %% node in which plate element
        platenode=FEM.elementNodes(TargetPlate,:);
        
        TargetGaussian=Stiffener.TargetGaussian(nodestifflag-4,:); %% xi and eta
        xi=TargetGaussian(1);
        eta=TargetGaussian(2);
        NshapeElement(eee,platenode)=...
            [1/4*(1-xi)*(1-eta)*(-xi-eta-1);
            1/4*(1+xi)*(1-eta)*(xi-eta-1);
            1/4*(1+xi)*(1+eta)*(xi+eta-1);
            1/4*(1-xi)*(1+eta)*(-xi+eta-1);
            1/2*(1-eta)*(1+xi)*(1-xi);
            1/2*(1+xi)*(1+eta)*(1-eta);
            1/2*(1+eta)*(1+xi)*(1-xi);
            1/2*(1-xi)*(1+eta)*(1-eta)]';
        
        % w, thetaX, thetaY, u, v
        stiffenernodeID=nodestifflag-4;
        
        ShapeMatrixElem(stiffenernodeID,platenode)=ShapeMatrixElem(stiffenernodeID,platenode)+ ...
            NshapeElement(eee,platenode);
        
    end
end
ShapeMatrixSP(1:Stiffener.nodenum,1:FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum:Stiffener.nodenum*2,1+FEM.NodeNumber:2*FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum*2:Stiffener.nodenum*3,1+FEM.NodeNumber*2:3*FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum*3:Stiffener.nodenum*4,1+FEM.NodeNumber*3:4*FEM.NodeNumber)=ShapeMatrixElem;

ShapeMatrixSP(1+Stiffener.nodenum*4:Stiffener.nodenum*5,1+FEM.NodeNumber*4:5*FEM.NodeNumber)=ShapeMatrixElem;