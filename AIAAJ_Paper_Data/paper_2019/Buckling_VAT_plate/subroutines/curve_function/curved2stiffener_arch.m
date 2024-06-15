clear xstiffener YYstiffener XXstiffener ystiffener;

%% Case I
% Xcenter=0.300;
% Ycenter=-0.100;
% radius=0.375;
% 
% theta1=atan((0.125+abs(Ycenter))/Xcenter);
% theta2=pi-theta1;
% 
% theta=linspace(theta1,theta2,stiffenernodesnumber);
% 
% 
% 
% 
% 
% % x-stiffener
% XXstiffener(1,:)=Xcenter-radius.*cos(theta);XXstiffener(1,1)=0;XXstiffener(1,end)=0.6;
% YYstiffener(1,:)=sqrt(radius^2-(XXstiffener(1,:)-0.300).^2)-0.100;
% 
% 
% % y-stiffener
% XXstiffener(2,:)=XXstiffener(1,:);
% YYstiffener(2,:)=-sqrt(0.375^2-(XXstiffener(1,:)-0.300).^2)+0.9;
% 
% % figure;plot(XXstiffener(1,:),YYstiffener(1,:),'ro',XXstiffener(2,:),YYstiffener(2,:),'ro')
% % axis([0 0.6 0 0.8]);axis image;

%% Case II

Xcenter=0.300;
Ycenter=-0.65;
radius=0.925;

theta1=atan((0.225+abs(Ycenter))/Xcenter);
theta2=pi-theta1;

theta=linspace(theta1,theta2,stiffenernodesnumber);

% x-stiffener
XXstiffener(1,:)=Xcenter-radius.*cos(theta);XXstiffener(1,1)=0;XXstiffener(1,end)=0.6;
YYstiffener(1,:)=sqrt(radius^2-(XXstiffener(1,:)-0.300).^2)-0.650;


% y-stiffener
XXstiffener(2,:)=XXstiffener(1,:);
YYstiffener(2,:)=-sqrt(radius^2-(XXstiffener(1,:)-0.300).^2)+1.45;
% figure;plot(XXstiffener(1,:),YYstiffener(1,:),'ro',XXstiffener(2,:),YYstiffener(2,:),'ro')
% hold on;axis([0 0.6 0 0.8]);axis image; hold off;