function [ A,B,C ] = get_line( X,Y,flag )
% Input:
%   if flag = 'a': X is the angle of the line (to be got) with regard to the
%                x-axis; 
%                Y is the point the line passes.
%   if flag = 'p' (default): X and Y are two points the line passes.
% Output:
%   [A,B,C]: A*x+B*y+C=0
if nargin == 2
    flag = 'p'; % default: X and Y are two points
end
epsilon = 1e-15;
if flag == 'a' % X is an angle, Y is a point
    theta = X;
    X = Y;
    
    if abs(theta-0)<epsilon || abs(theta - pi)<epsilon % theta == 0 || theta == pi
        A = 0;
        B = 1;
        C = -X(2);
    elseif abs(theta - pi/2) < epsilon % theta == pi/2
        A = 1;
        B = 0;
        C = -X(1);
    else
        k = tan(theta); % scope
        A = k;
        B = -1;
        C = X(2) - k*X(1);
    end
elseif flag == 'p' % X and Y are two points
    if abs(X(1)-Y(1))<epsilon && abs(X(2)-Y(2))<epsilon
        error 'Only one point cannot determine a line.'
    end
    if abs(X(1)-Y(1))<epsilon
        A = 1;
        B = 0;
        C = -X(1);
    elseif abs(X(2)-Y(2))<epsilon
        A = 0;
        B = 1;
        C = -X(2);
    else
        a = Y(1) - X(1);
        b = Y(2) - X(2);
        A = b;
        B = -a;
        C = a*X(2) - b*X(1);
    end
end

end

