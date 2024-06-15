%% Material Constant
%  Lamina parallel to the reference plane
function [Amatrix,Bmatrix,Dmatrix,Ashear,Qreduced,Qbar]=StiffenerCompositeParallel(Mat,Stiffener,layer,ebar)
theta=Stiffener.theta(layer);
thickness=Stiffener.thickness(layer);
display(['Laminated Composites Stiffener - Orthotropic Materials, Layer=' num2str(layer)])

Q11=Mat.E1/(1-Mat.v12*Mat.v21);
Q12=Mat.v12*Mat.E2/(1-Mat.v12*Mat.v21);
Q22=Mat.E2/(1-Mat.v12*Mat.v21);
Q66=Mat.G12;
Q44=Mat.G23;
Q55=Mat.G13;

%------Orthotropic Material--------------
c=cos(theta/180*pi);
s=sin(theta/180*pi);

Q11b=Q11*c^4+2*(Q12+2*Q66)*s^2*c^2+Q22*s^4;
Q12b=(Q11+Q22-4*Q66)*s^2*c^2+Q12*(s^4+c^4);
Q22b=Q11*s^4+2*(Q12+2*Q66)*s^2*c^2+Q22*c^4;
Q16b=(Q11-Q12-2*Q66)*s*c^3+(Q12-Q22+2*Q66)*s^3*c;
Q26b=(Q11-Q12-2*Q66)*c*s^3+(Q12-Q22+2*Q66)*c^3*s;
Q66b=(Q11+Q22-2*Q12-2*Q66)*s^2*c^2+Q66*(s^4+c^4);
Q44b=Q44*c^2+Q55*s^2;
Q45b=(Q55-Q44)*c*s;
Q55b=Q55*c^2+Q44*s^2;



Qbar=[Q11b,  Q12b,  0,   0,   Q16b;
      Q12b,  Q22b,  0,   0,   Q26b;
       0,      0,  Q44b, Q45b,  0;
       0,      0,  Q45b, Q55b,  0;
     Q16b,    Q26b, 0,    0,   Q66b];

bendingNO=[1,2,5];shearNO=[3,4];
D_bend=Qbar(bendingNO,bendingNO);
D_shear=Qbar(shearNO,shearNO);

bendNewNo=[1,3];shearNewNo=[2];
Qreduced_bend=D_bend(bendNewNo,bendNewNo)-D_bend(2,bendNewNo)'/D_bend(2,2)*D_bend(2,bendNewNo);
Qreduced_shear=D_shear(2,2)-D_shear(1,2)^2/D_shear(1,1);

zk1=-Stiffener.height/2+ebar+layer*thickness;
zk=-Stiffener.height/2+ebar+(layer-1)*thickness;
% zk=Stiffener.thickness/2-layer*thickness;
Amatrix=(zk1-zk)*Qreduced_bend;
Bmatrix=1/2*(zk1^2-zk^2)*Qreduced_bend;
Dmatrix=1/3*(zk1^3-zk^3)*Qreduced_bend;
Ashear=Mat.kappa*(zk1-zk)*Qreduced_shear;
Qreduced=[Q11, Q16b,0;
    Q16b, Q66, 0;
    0, 0, Q55];
%% End of diary