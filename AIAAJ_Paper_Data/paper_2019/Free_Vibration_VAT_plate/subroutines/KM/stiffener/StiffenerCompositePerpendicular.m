%% Material Constant
%  Lamina parallel to the reference plane
function [AmatrixReduced,BmatrixReduced,DmatrixReduced,AshearReduced,Qreduced,Qbar,Q_bend_new,Q_shear_new]=...
    StiffenerCompositePerpendicular(Mat,Stiffener,layer)
theta=Stiffener.theta(layer);
thickness=Stiffener.thickness(layer);
display(['Laminated Composites Stiffener - Orthotropic Materials, Layer=' num2str(layer)])
c=cos(theta/180*pi);
s=sin(theta/180*pi);
%--------Transform Matrix for Coordinate Transformation--------
Q11=Mat.E1/(1-Mat.v12*Mat.v21);
Q12=Mat.v12*Mat.E2/(1-Mat.v12*Mat.v21);
Q22=Mat.E2/(1-Mat.v12*Mat.v21);

Q66=Mat.G12;Q55=Mat.G13;
Q44=Mat.G23;


%------Orthotropic Material--------------
Q11b=Q11*c^4+2*(Q12+2*Q66)*s^2*c^2+Q22*s^4;
Q12b=(Q11+Q22-4*Q66)*s^2*c^2+Q12*(s^4+c^4);
Q22b=Q11*s^4+2*(Q12+2*Q66)*s^2*c^2+Q22*c^4;
Q16b=(Q11-Q12-2*Q66)*s*c^3+(Q12-Q22+2*Q66)*s^3*c;
Q26b=(Q11-Q12-2*Q66)*c*s^3+(Q12-Q22+2*Q66)*c^3*s;
Q66b=(Q11+Q22-2*Q12-2*Q66)*s^2*c^2+Q66*(s^4+c^4);
Q44b=Q44*c^2+Q55*s^2;
Q45b=(Q55-Q44)*c*s;
Q55b=Q55*c^2+Q44*s^2;


Qbar=[Q11b,  Q12b,  0,   0,  Q16b;
      Q12b,  Q22b,  0,   0,  Q26b;
       0,      0,  Q44b, Q45b, 0;
       0,      0,  Q45b, Q55b, 0;
      Q16b,   Q26b, 0,    0,   Q66b];

bendingNO=[1,2,5];
Q_bend=Qbar(bendingNO,bendingNO);
shearNO=[3,4];
Q_shear=Qbar(shearNO,shearNO);

bendNewNo=[1,3];

Q_bend_new=Q_bend(bendNewNo,bendNewNo)-Q_bend(2,bendNewNo)'/Q_bend(2,2)*Q_bend(2,bendNewNo);
Q_shear_new=Q_shear(2,2)-Q_shear(1,2)^2/Q_shear(1,1);


zk1=-Stiffener.width/2+layer*thickness;
zk=-Stiffener.width/2+(layer-1)*thickness;

AmatrixReduced=(zk1-zk)*Q_bend_new;
BmatrixReduced=1/2*(zk1^2-zk^2)*Q_bend_new;
DmatrixReduced=1/3*(zk1^3-zk^3)*Q_bend_new;
AshearReduced=Mat.kappa*(zk1-zk)*Q_shear_new;
Qreduced=[Q11b, Q16b,0;
    Q16b, Q66b, 0;
    0, 0, Q55b];

%% End of diary