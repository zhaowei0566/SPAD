function [Tmatrix,alpha]=TransformationMatrix(FEM,Stiffener)
% Transform local coordinate system to global coordinate system

% Tmatrix=zeros(Stiffener.nodenum,Stiffener.nodenum);

alpha=zeros(Stiffener.nodenum,1);

XI=[-1,1,0];

for ii=1:Stiffener.elemnum
    
    NodeIndices=Stiffener.element(ii,2:end);
    XYZ=Stiffener.PointsCoord(NodeIndices,2:3);
    
    for ee=[1 2 3]
        
        xi=XI(ee);
        
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);
        
        if naturalderivatives*XYZ(:,1)>=0 && naturalderivatives*XYZ(:,2)>=0
            
            alpha(NodeIndices(ee))=atan( (naturalderivatives*XYZ(:,2)) / (naturalderivatives*XYZ(:,1)) );
            
        elseif naturalderivatives*XYZ(:,1)<=0 && naturalderivatives*XYZ(:,2)>=0
            
            alpha(NodeIndices(ee))=pi+atan( (naturalderivatives*XYZ(:,2)) / (naturalderivatives*XYZ(:,1)) );
            
        elseif naturalderivatives*XYZ(:,1)<0 && naturalderivatives*XYZ(:,2)<0
            
            alpha(NodeIndices(ee))=-pi+atan( (naturalderivatives*XYZ(:,2)) / (naturalderivatives*XYZ(:,1)) );
            
        elseif naturalderivatives*XYZ(:,1)>=0 && naturalderivatives*XYZ(:,2)<=0
            
            alpha(NodeIndices(ee))=atan( (naturalderivatives*XYZ(:,2)) / (naturalderivatives*XYZ(:,1)) );
            
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=Stiffener.nodenum;
Tmatrix1=eye(n);
Tmatrix2=diag(cos(alpha));
Tmatrix3=diag(sin(alpha));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w, theta_x, theta_y, u, v for the plate displacement fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmatrix(1:n,1:n)=Tmatrix1;

Tmatrix(1+n:2*n,1+n:2*n)=Tmatrix2;
Tmatrix(1+n:2*n,1+n*2:3*n)=Tmatrix3;

Tmatrix(1+n*2:3*n,1+n:2*n)=-Tmatrix3;
Tmatrix(1+n*2:3*n,1+n*2:3*n)=Tmatrix2;

Tmatrix(1+n*3:4*n,1+n*3:4*n)=Tmatrix2;
Tmatrix(1+n*3:4*n,1+n*4:5*n)=Tmatrix3;

Tmatrix(1+n*4:5*n,1+n*3:4*n)=-Tmatrix3;
Tmatrix(1+n*4:5*n,1+n*4:5*n)=Tmatrix2;