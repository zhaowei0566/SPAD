function p_spl=Bspline_web(n,order,xinput,stiffenernodenumber)

% function spline(n,order)
%
% Plots the B-slpine-curve of n control-points.
% The control points can be chosen by clicking
% with the mouse on the figure.
%
% COMMAND:  spline(n,order)
% INPUT:    n     Number of Control-Points
%           order Order ob B-Splines
%                 Argnument is arbitrary
%                 default: order = 4
%
% Date:     2007-11-28
% Author:   Stefan Hüeber

% if (nargin ~= 2)
% 	order = 4;
% end
nplot = 100;

if (n < order)
    display([' !!! Error: Choose n >= order=',num2str(order),' !!!']);
    return;
end


% t = linspace(0,1,nplot);

for i = 1:n
%     title(['Choose ',num2str(i),' th. control point']);
    pointCord(i,:) = xinput(i,:);
%     hold off;
%     figure(1);hold on; box on;
%     set(gca,'Fontsize',16);
%     plot(pointCord(:,1),pointCord(:,2),'k:','LineWidth',2);hold on;
%     axis([-1 1 -1 1]);
%     hold on; box on;
    if (i  >= order)
        T = linspace(0,1,i-order+2);
        y = linspace(0,1,stiffenernodenumber);
        p_spl = DEBOOR(T,pointCord,y,order);
%         figure(1);hold on;
%         plot(p_spl(:,1),p_spl(:,2),'b-','LineWidth',4);
    end
%     plot(pointCord(:,1),pointCord(:,2),'ro','MarkerSize',10,'MarkerFaceColor','r');hold on;
end

% title(['B-Spline-curve with ',num2str(n),' control points of order ',num2str(order)]);
% hold on;plot(xinput(1,:),xinput(3,:),'k:','LineWidth',2)
