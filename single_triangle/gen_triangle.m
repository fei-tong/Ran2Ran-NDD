% the function for generating a triangle
% Input parameters: 3 angles ([0,180]): alpha,beta,gamma
%             the lenght of the opposite edge of angle bets: a
% Output: three decreasing angles: X,Y,Z
%         three edge lengths: a,b,c
%         three altitudes: h_a, h_b, h_c
%         three vertices: A,B,C
%
% Author: Fei T.
% Date: Nov. 05, 2013
% Last update: Nov. 06, 2013

function [X,Y,Z,a,b,c,h_a,h_b,h_c,A,B,C]=gen_triangle(alpha,beta,gamma,a)
%sort the angle: X >= Y >= Z
X = alpha; Y = beta; Z = gamma;
if beta > X
    X = beta; Y =alpha;
elseif gamma > X;
    X = gamma; Z = alpha;
end
if Y < Z
    temp = Y;
    Y=Z;
    Z=temp;
end

X=X*pi/180; Y=Y*pi/180; Z=Z*pi/180;
h_b=a*sin(Z); h_c=a*sin(Y);
if X <= pi/2
    b=h_c/sin(X); c=h_b/sin(X);
else
    b=h_c/sin(pi-X); c=h_b/sin(pi-X);
end
h_a = b*sin(Z); 
A=[b*cos(Z) h_a]; B=[a 0]; C = [0 0];
end