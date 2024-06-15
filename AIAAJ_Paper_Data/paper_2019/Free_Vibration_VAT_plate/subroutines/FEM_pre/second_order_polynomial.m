function XYstiffener = second_order_polynomial(two_pnt_cords,pnts_number,Length,Width)


x2=two_pnt_cords(1,1);  y2=two_pnt_cords(1,2);
x1=Length/2;            y1=Width/2;
x3=two_pnt_cords(2,1);  y3=two_pnt_cords(2,2);


a = (x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))/((x1 - x2) *(x1 - x3)* (x2 - x3));
b = (x3^2 *(y1 - y2) + x1^2 *(y2 - y3) + x2^2 *(-y1 + y3))/((x1 - x2) *(x1 - x3)* (x2 - x3));
c =( x3* (x2* (x2 - x3) *y1 + x1 *(-x1 + x3) *y2) + x1* (x1 - x2) *x2* y3)/((x1 - x2) *(x1 - x3)* (x2 - x3));

xmin = min([x1 x2 x3]);

xmax = max([x1 x2 x3]);
x = linspace(xmin,xmax,pnts_number);
y = a.*x.^2 + b.*x + c;
figure;
plot(x,y,'b-o');
hold on;
plot([x1 x2 x3],[y1 y2 y3],'ro');axis([0 Length 0 Width])


XYstiffener = [x;y];



polyfit([x1 x2 x3],[y1 y2 y3],2)