function[ KG_Stiffener,KGbeam_global ]= GeometryStiffnessStiffener_Orthotropic_main_v2(FEM,Stiffener,ebar,Ds)

% compute geometric stiffness for stiffness using stress obtained from
% stress recover


StiffenerGaussPoint = Stiffener.StiffenerGaussPoint;

KG_Stiffener=zeros(FEM.GDof,FEM.GDof); %

stiffeners_number = size(Stiffener.Tmatrix,3);

% % %  Stiffener Eccentricity % % %
% ebar=(Stru.thickness+Stiffener.height)/2;

%% Geometric stiffness for stiffeners


for stiffno = 1:stiffeners_number
    
    %
    Tmatrix = Stiffener.Tmatrix(:,:,stiffno);
    
    %
    ShapeMatrixSP=Stiffener.ShapeMatrixSPTest(:,:,stiffno);
    
    
    
   % Geometric stiffness for each stiffener described in local coordinate
   % system
    KGbeam_local = GeometryStiffnessStiffener_stress_recovery_v2(FEM,Stiffener,StiffenerGaussPoint,ebar,stiffno);
    
    
    % transform it to global coordinae system
    KGbeam_global = Tmatrix'*KGbeam_local*Tmatrix;
    
    
    
    % transform it to the panel's geometric stiffness
    KG_Stiffener = KG_Stiffener+ShapeMatrixSP'*KGbeam_global*ShapeMatrixSP;
    
    
end