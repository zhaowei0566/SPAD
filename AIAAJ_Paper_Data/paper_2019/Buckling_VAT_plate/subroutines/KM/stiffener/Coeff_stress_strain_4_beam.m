function [C11bar,C55bar,C66bar,C] = Coeff_stress_strain_4_beam(Mat)

% reduce 6by6 to 3by3 matrix





Compliance_matrix = [1/Mat.E1 -Mat.v21/Mat.E2 -Mat.v31/Mat.E3 0 0 0;
                      -Mat.v12/Mat.E1 1/Mat.E2 -Mat.v32/Mat.E3  0 0 0;
                      -Mat.v13/Mat.E1 -Mat.v23/Mat.E2 1/Mat.E3 0 0 0;
                               0        0       0     1/Mat.G23    0           0;
                                0      0         0      0     1/Mat.G13       0;
                               0      0          0      0          0      1/Mat.G12];
                 
                 
C = inv(Compliance_matrix);

C11bar_n = C(2,3)^2 - C(2,2)*C(3,3);
C11bar_d = C(1,3)^2*C(2,2)-2*C(1,2)*C(1,3)*C(2,3) + C(1,1)*C(2,3)^2+C(1,2)^2*C(3,3) - C(1,1)*C(2,2)*C(3,3);

C11bar = C11bar_d/C11bar_n;



C55bar = C(5,5);

C66bar = C(6,6);


% CBAR = diag([C11bar,C66bar,C55bar]);




