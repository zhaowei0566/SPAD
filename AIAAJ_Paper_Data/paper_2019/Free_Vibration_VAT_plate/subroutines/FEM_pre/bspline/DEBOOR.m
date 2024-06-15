function val = DEBOOR(T,p,y,order)

% function val = DEBOOR(T,p,y,order)
%
% INPUT:  T     Stützstellen
%         p     Kontrollpunkte (nx2-Matrix)
%         y     Auswertungspunkte (Spaltenvektor)
%         order Spline-Ordnung
%
% OUTPUT: val   Werte des B-Splines an y (mx2-Matrix)
%
% Date:   2007-11-27
% Author: Jonas Ballani
% p;
m = size(p,1);
n = length(y);
X = zeros(order,order);
Y = zeros(order,order);
a = T(1);
b = T(end);
T = [ones(1,order-1)*a,T,ones(1,order-1)*b];


for l = 1:n
    t0 = y(l);
    id = find(t0 >= T);
    k = id(end);
    if (k > m)
        return;
    end
    X(:,1) = p(k-order+1:k,1);
    Y(:,1) = p(k-order+1:k,2);
    
    for i = 2:order
        for j = i:order
            num = t0-T(k-order+j);
            if num == 0
                weight = 0;
            else
                s = T(k+j-i+1)-T(k-order+j);
                weight = num/s;
            end
            X(j,i) = (1-weight)*X(j-1,i-1) + weight*X(j,i-1);
            Y(j,i) = (1-weight)*Y(j-1,i-1) + weight*Y(j,i-1);
        end
    end
    val(l,1) = X(order,order);
    val(l,2) = Y(order,order);
end


