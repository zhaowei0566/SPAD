clear xstiffener YYstiffener XXstiffener ystiffener;
% x=0:0.001:0.12;
XXstiffener(1,:)=linspace(0,Stru.length,stiffenernodesnumber);
stiffener_xcord=XXstiffener(1,:);
YYstiffener(1,:)=0.00.*sin(2*pi/0.12.*stiffener_xcord)+0.06;
% figure;
% plot(x,y);
% axis([0 0.12 0 0.12]);
