function [d_array,pdd_cdf] = f_sim_rand2rand_tiered_polygon_23( x11,y11,x22,y22)
% Simulation: to get the Rand2Rand distance distribution between two
%             arbitrary polygons (x11,y11) and (x22,y22)
% Input: 
%   (x11,y11): the first polygon. 
%   (x22,y22): the second polygon.
%
% Output:
%   [ d_array, pdd_cdf ]: Rand2Rand distance distribution between the two
%       polygons.
% Author: Fei Tong
% Date: May. 18, 2016

if nargin < 1
    error 'This function needs at least 1 argument.'
end

max_iter = 10000;
count = 0;
distance_in_sim=zeros(1,max_iter); % distance of two random points with a triangle


x11_length = max(x11)-min(x11);
y11_length = max(y11)-min(y11);

x22_length = max(x22)-min(x22);
y22_length = max(y22)-min(y22);

while count < max_iter
    % generage the 1st random point
    while 1
        x1 = min(x11)+x11_length*rand(1);
        y1 = min(y11)+y11_length*rand(1); % h_a is the altitude from point A to edge BC
        IN = inpolygon(x1,y1,x11,y11);
        IN2 = inpolygon(x1,y1,x22,y22);
        if IN == 1 && IN2~=1
            break;
        end
    end
    
    % generage the 2nd random point
    while 1
        x2 = min(x22)+x22_length*rand(1);
        y2 = min(y22)+y22_length*rand(1);
        IN = inpolygon(x2,y2,x22,y22);
        if IN == 1
            break;
        end
    end
    
    dist = sqrt((x1-x2)^2+(y1-y2)^2);
    count = count + 1;
    distance_in_sim(count)=dist;
end

[cc,x] = ecdf(distance_in_sim);% from simulation
point_num = 40;
pdd_cdf = zeros(1,point_num);
d_array = zeros(1,point_num);
tick = int32(length(x)/point_num);
for i =1:point_num
    next = i*tick;
    if next >length(x)
        next = length(x);
        d_array(i) = x(next);
        pdd_cdf(i) = cc(next);
        break;
    end
    d_array(i) = x(next);
    pdd_cdf(i) = cc(next);
end

end