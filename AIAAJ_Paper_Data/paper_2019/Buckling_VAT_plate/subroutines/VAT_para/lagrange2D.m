% function lagrange2D
%% STARTVARIABLES
% 
clear all; close all;
[X,Y]=meshgrid(-1:1:1,-1:1:1);
Z = sin(Y)*sin(X);



[X,Y] = meshgrid(linspace(0,0.127,3),linspace(0,0.127,3));



VAT.Theta1 = [17 12 45
             67 50 51
             71 49.5 71.5 ];


VAT.Theta2 = [ 14 11.5  6
       -65 -54 -50.5
-72.5 -59 -59.5 ];
           
           

Z = VAT.Theta1;

step = 0.03;
%% INTERPOLATION
% erst kurven in x richtung ausrechnen
xCurves={};
for i=1:size(X,1)
    % x vector must be a column vector for the lagrange function
    x = X(i,:)';
    % z vector must be a column vector for the lagrange function
    z = Z(i,:)';
    p=[];
    for j = x(1):step:x(end)
        % interpolate for every parameter value j and add it to p
        p = [p,lagrange1d(x,z,j)];
    end
    % save curves in x direction
    xCurves{i} = p;
end

y = Y(:,1);
% matrix for the graphical outpu
A=[];
% interpolate in y-direction
for i=1:length(xCurves{1})
    p=[];
    z=[];
    for l=1:length(y)
        z = [z;xCurves{l}(i)];
    end
    for j = y(1):step:y(end)
         % interpolate for every parameter value j and add it to p
        p = [p;lagrange1d(y,z,j)];
    end
    A = [A,p];
end

%% GRAPHICAL OUTPUT
% plot surface
surf(x(1):(x(end)-x(1))/(size(A,1)-1):x(end),y(1):(y(end)-y(1))/(size(A,1)-1):y(end),A);
hold on;
% plot points
for i=1:size(X,1)
    for j=1:size(Y,1)
        p = plot3(X(i,j),Y(i,j),Z(i,j));
        set(p,'Marker','.');
        set(p,'MarkerSize',30);
    end
end

% end
