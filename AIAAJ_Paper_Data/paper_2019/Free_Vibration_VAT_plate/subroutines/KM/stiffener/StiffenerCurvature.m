function  Curvature=StiffenerCurvature(XYZ,detJ, naturalderivatives, d2Nds2)

% =========================================================================
%    Be careful on the stiffener orientation
% =========================================================================

% X=XYZ(:,1);Y=XYZ(:,2);
% dxds=naturalderivatives*X;
% dyds=naturalderivatives*Y;
% 
% d2yds2=d2Nds2*Y;
% d2xds2=d2Nds2*X;
% Curvature =(dyds*d2xds2-dxds*d2yds2)/(detJ)^3*0;
% if abs(Curvature)<=eps
%     Curvature=0;
% end
% % % det(Jacob)=length/2;
% % dy2dx2=(naturalderivatives(1,:)*X*d2Nds2*Y-naturalderivatives(2,:)*Y*d2Nds2*X)...
% %     /(naturalderivatives(1,:)*X)^3;
% % Curvature=dy2dx2/(1+dy2dx2^2)^1.5;


X=XYZ(:,1);Y=XYZ(:,2);


dxds=naturalderivatives*X;
dyds=naturalderivatives*Y;
% dydx=dyds/dxds;

d2yds2=d2Nds2*Y;
d2xds2=d2Nds2*X;
% d2ydx2=-(dxds*d2yds2-dyds*d2xds2)/(dxds)^3;

Curvature= (dyds*d2xds2-dxds*d2yds2)/(detJ)^3;

