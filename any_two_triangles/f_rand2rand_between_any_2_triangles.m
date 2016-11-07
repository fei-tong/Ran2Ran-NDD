function [ d_array, pdf_array, cdf_array ] = f_rand2rand_between_any_2_triangles( t1,t2,d_step )
% This function is to get the Rand2Rand distance distribution between any
% two triangles.
%
% Input:
%   t1: triangle 1, t1 is like [x1 y1;x2 y2;x3 y3], where [xi yi] is a
%       vertex of t1.
%   t2: triangle 2, t2 is similar to t1;
%   d_step: A large number which determines the accuracy of the results.
%
% Output:
%   [d_array, cdf_array]: Numerical result of the Rand2Rand distance
%                         distribution between 2 triangles
% Author: Fei Tong
% Date: May. 11, 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin < 4
%     error 'This function needs at least 4 arguments.'
% end
% if nargin == 4
%     d_step = 1000;
% end
%% Move and rotate
% this section will be commented
% % clear;clc;
% % a = 1;
% % A_1 = 30*pi/180;
% % A_2 = 90*pi/180;
% % A_3 = 110*pi/180;
% % 
% % % A_1 = 36*pi/180;
% % % A_2 = 72*pi/180;
% % % A_3 = 108*pi/180;
% % % A_4 = (180-36)*pi/180;
% % b = a/sqrt(2*(1-cos(A_3)));
% % 
% % A = [0 b*sin(A_3)]; D = [-b*cos(A_3) 0]; B = [-b*cos(A_3)+a*cos(A_2) a*sin(A_2)];
% % C = [b-b*cos(A_3) 0]; Bp = [b-2*b*cos(A_3) a*sin(A_1)]; % Bp is B'
% % B = [0.6 0.4];
% % x = [D(1) C(1) Bp(1) B(1) A(1)];
% % y = [D(2) C(2) Bp(2) B(2) A(2)];
% % 
% % line([x x(1)],[y y(1)]);
% % text(D(1),D(2),'D');text(Bp(1),Bp(2),'Bp');text(B(1),B(2),'B');text(C(1),C(2),'C');text(A(1),A(2),'A');
% % 
% % x1 = [Bp(1) B(1) D(1)];
% % y1 = [Bp(2) B(2) D(2)];
% % x2 = [A(1) B(1) D(1)];
% % y2 = [A(2) B(2) D(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = t1(:,1)';
y1 = t1(:,2)';
x2 = t2(:,1)';
y2 = t2(:,2)';
vx = min([x1 x2]);
vy = min([y1 y2]);

M = [0 0] - [vx vy]; % moving vector
for i = 1:length(x1) % move all vertexes correspondingly
    x1(i) = x1(i) + M(1);
    x2(i) = x2(i) + M(1);
    y1(i) = y1(i) + M(2);
    y2(i) = y2(i) + M(2);
end
% figure;
% line([x1 x1(1)],[y1 y1(1)]);
% line([x2,x2(1)],[y2,y2(1)]);
%% prepare d_array and pdf_array
x = [x1 x2];
y = [y1 y2];
max_d = 0;
for i = 1:length(x)-1
    for j = i+1:length(x)
        d = norm([x(i) y(i)]-[x(j) y(j)]);
        if max_d < d
            max_d = d;
        end
    end
end
% d_step = 1000;
delta_d = 1/d_step; 
d_array = 0:delta_d:max_d;
pdf_array = zeros(1,length(d_array));
%% area
s1 = shoelace(x1,y1);
s2 = shoelace(x2,y2);
%%
delta_theta = pi/360;
theta_array = 0:delta_theta:pi;
for theta = theta_array
%     theta
    [ L ] = f_theta_func_between_any_2_triangles( theta,t1,t2 );
    if length(L)==1 && L == 0
        continue;
    end
    [n,~] = size(L);
    for i = 1:n
        L1 = L(i,1); L2 = L(i,2); L3 = L(i,3);
        
        L_low = L1+L2; L_up = L2+L3;
        L_s = L1;
        if L_low > L_up
            L_low = L2+L3; L_up = L1+L2;
            L_s = L3;
        end
        
        for j = 1:length(d_array)
            g = 0;
            if d_array(j)>=L2 && d_array(j)<=L_low
                g = 2*d_array(j)*(d_array(j)-L2);
            elseif d_array(j)>=L_low && d_array(j) <= L_up
                g = 2*d_array(j)*L_s;
            elseif d_array(j) <= L1+L2+L3 && d_array(j) >= L_up
                g = 2*d_array(j)*(L1+L2+L3 - d_array(j));
            end
            pdf_array(j) = pdf_array(j) + g/(s1*s2);%*u/(s1+s2)^2;
        end
    end
end

cdf_array = zeros(1,length(pdf_array));
for i = 1:length(cdf_array)
    cdf_array(i)=sum(pdf_array(1:i));
end
cdf_array = cdf_array/cdf_array(length(cdf_array));
pdf_array = pdf_array*(d_step/sum(pdf_array));
end