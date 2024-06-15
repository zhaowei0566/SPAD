%% nested functions
function res = lagrange1d(x,z,t)
% lagrange interpolation in 1d
% input: x...vektor, welcher die x (bzw. y) Werte beinhaltet

res = 0;
for i=1:length(x)
    % delete the component not needed
    temp = [x(1:i-1);x(i+1:end)];
    
    denominator = x(i)-temp;
    
    numerator = t-temp;
    res = res + z(i)*(prod(numerator))/(prod(denominator));
end
end