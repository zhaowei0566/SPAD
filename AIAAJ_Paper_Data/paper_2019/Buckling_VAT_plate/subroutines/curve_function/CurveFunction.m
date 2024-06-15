%--- Curved Stiffener Function -------------
clear xstiffener YYstiffener XXstiffener ystiffener;
% xstiffener=linspace(0,Stru.length,stiffenernodesnumber);
%
% xstiffenerAngle=linspace(0,pi/2,stiffenernodesnumber);
% % xstiffener=0.11-0.11.*cos(xstiffenerAngle)+0.01;
% %     ystiffener=linspace(0,Stru.length,stiffenernodesnumber);



%% x-stiffener and y-stiffener
% % x-stiffener
XXstiffener(1,:)=linspace(0,Stru.length,stiffenernodesnumber);
YYstiffener(1,:)=Stru.width/2*ones(1,stiffenernodesnumber);% y=b/2
% 
% y-stiffener
YYstiffener(2,:)=linspace(0,Stru.width,stiffenernodesnumber);
XXstiffener(2,:)=Stru.length/2*ones(1,stiffenernodesnumber);% x=a/2


%% 2 Parallel stiffeners
% 
%   xstiffener=linspace(0,Stru.length,stiffenernodesnumber);
%    
%    XXstiffener(1,:)=xstiffener;
%    XXstiffener(2,:)=xstiffener;
%    YYstiffener(1,:)=Stru.width*1/4*ones(1,size(xstiffener,2));% y=b/2
%    YYstiffener(2,:)=3*Stru.width/4*ones(1,size(xstiffener,2));% y=2b/3

