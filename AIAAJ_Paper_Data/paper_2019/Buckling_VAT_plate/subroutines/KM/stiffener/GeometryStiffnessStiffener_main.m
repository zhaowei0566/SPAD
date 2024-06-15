function KG_Stiffener = GeometryStiffnessStiffener_main(FEM,Stiffener,Mat,StiffenerGaussPoint,ebar)

% compute geometric stiffness for stiffness using stress obtained from
% stress recover

KG_Stiffener=zeros(FEM.GDof,FEM.GDof); %


stiffeners_number = size(Stiffener.Tmatrix,3);



% % %  Stiffener Eccentricity % % %
% ebar=(Stru.thickness+Stiffener.height)/2;

%% Stress recovery for stiffeners

[Stiffener.stress,Stiffener.strain,Stiffener.displacement,Stiffener.Sdisplacement_Local]=...
    StressRecoveryStiffeners(FEM,Mat,Stiffener,StiffenerGaussPoint,ebar);




%% Geometric stiffness for stiffeners


for stiffno = 1:stiffeners_number
    
    Tmatrix = Stiffener.Tmatrix(:,:,stiffno);
    ShapeMatrixSP=Stiffener.ShapeMatrixSPTest(:,:,stiffno);
    
   
    KGbeam = GeometryStiffnessStiffener_stress_recovery_v1(FEM,ebar,Stiffener,StiffenerGaussPoint);
    KGbeam = Tmatrix'*KGbeam*Tmatrix;
    KG_Stiffener = KG_Stiffener+ShapeMatrixSP'*KGbeam*ShapeMatrixSP;
    
    
end