function [d_array,pdd_pdf,pdd_cdf] = f_formula_pdd_pdf_triangle(x,y,d_step)%(X,Y,Z,a,b,c)
% This function is to get the Rand2Rand distance distribution within a
% single triangle, including both pdf and cdf.
% Input:
%   (x,y): the 3 vertexes of the triangle.
%   d_step: A large number which determines the accuracy of the results.
% Output:
%   numerical results of pdf: [d_array,pdd_pdf] and cdf: [d_array,pdd_cdf]
%
% Author: Fei Tong
% Date: May. 10, 2016

if nargin < 2
    error 'This function needs at least 2 arguments.'
end
if nargin == 2
    d_step = 1000;
end
%% move and rotate
% M = [0 0] - [x(1) y(1)];
% for i = 1:length(x)
%     x(i) = x(i)+M(1);
%     y(i) = y(i)+M(2);
% end
% r = norm([x(2) y(2)]-[x(1) y(1)]);
% theta = acos(x(2)/r); % rotating angle
% rotation_x = [cos(theta) sin(theta)];%rotation matrix
% rotation_y = [-sin(theta) cos(theta)];
% v_x = rotation_x * [x;y];
% v_y = rotation_y * [x;y];
% x = v_x;
% y = v_y;
% 
% % text(x(1),y(1),'1');text(x(2),y(2),'2');text(x(3),y(3),'3');text(x(4),y(4),'4');text(x(5),y(5),'5');
% text(x(1),y(1),'1');text(x(2),y(2),'2');text(x(3),y(3),'3');
% line([x x(1)],[y y(1)],'color','r');
% figure;
%%
aa = norm([x(1),y(1)]-[x(2),y(2)]);
bb = norm([x(1),y(1)]-[x(3),y(3)]);
cc = norm([x(2),y(2)]-[x(3),y(3)]);

a = max(aa,max(bb,cc));
c = min(aa,min(bb,cc));
b = aa+bb+cc - a - c;

X = acos((b^2+c^2-a^2)/(2*b*c));
Y = acos((c^2+a^2-b^2)/(2*c*a));
Z = acos((a^2+b^2-c^2)/(2*a*b));

% d_step = 1000;
delta_d = 1/d_step; 
% delta=1/10000; 
d_array=0:delta_d:a; 
pdd_pdf=zeros(1,length(d_array));
pdd_cdf=zeros(1,length(d_array));
ff1 = zeros(1,length(d_array));
ff2 = zeros(1,length(d_array));
ff3 = zeros(1,length(d_array));

% lower = 0; upper = Z;%%% case2: %theta \in [0, Z]
% gg1 = b*cos(Z-upper)-a*cos(upper)-b*cos(Z-lower)+a*cos(lower);
% lower = Z; upper = pi-Y;%%% case4: theta \in [Z, pi-Y]
% gg2 = -c*cos(Y+upper)-b*cos(Z-upper)+c*cos(Y+lower)+b*cos(Z-lower);
% lower = pi-Y; upper = pi;%%% case6: theta \in [pi-Y, pi]
% gg3 = c*cos(Y+upper)-a*cos(upper)-c*cos(Y+lower)+a*cos(lower);

for i = 1:length(d_array) %theta \in [0, pi]
    d = d_array(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% case2: %theta \in [0, Z]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    lower = 0; upper = Z; flag = 1; 
    ff1(i) = f_specific_d_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% case4: theta \in [Z, pi-Y]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    lower = Z; upper = pi-Y; flag = 2;
    ff2(i) = f_specific_d_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% case6: theta \in [pi-Y, pi]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    lower = pi-Y; upper = pi; flag = 3;
    ff3(i) = f_specific_d_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag);
    
    pdd_pdf(i) = ff1(i)+ff2(i)+ff3(i);
end
% pdd_pdf = real(pdd_pdf);

for i = 1:length(pdd_cdf)
    pdd_cdf(i)=sum(pdd_pdf(1:i));
end
pdd_cdf = pdd_cdf/pdd_cdf(length(pdd_cdf));

pdd_pdf = pdd_pdf*(d_step/sum(pdd_pdf));
end