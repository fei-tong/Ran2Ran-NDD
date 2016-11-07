function [d_array,pdd_cdf,sim_density] = f_sim_rand2rand_arbitrary_polygon( triangle_cell,varargin)
% Simulation: to get the Rand2Rand distance distribution within an
%             arbitrary polygon triangulated into several triangules,
%             stored in triangle_cell
% Input: 
%   triangle_cell: Triangulated triangules of an arbitrary polygon. 
%
%   varargin contains the following one optional argument:
%       -> density array: (n)*1, containing "node densities" of each
%          triangle (n: the number of triangles).
% Output:
%   [ d_array, pdd_cdf, sim_density ]: Rand2Rand distance distribution within the
%                      polygon, with the simulated density of each triangle.
% Author: Fei Tong
% Date: May. 12, 2016

if nargin < 1
    error 'This function needs at least 1 argument.'
end

% s = shoelace(x,y); % total area
t_n = length(triangle_cell); % number of triangles: totally, there are n-2 triangles
if isempty(varargin{1})
    density_array = ones(t_n,1);
else
    density_array = varargin{1};
end
[q1,q2] = size(density_array);
if q1 == 1 && q2 > 1
    density_array = density_array';
end
q = length(density_array);
if q ~= t_n
    density_array = ones(t_n,1);
    disp('The size of density_array is not right. We will assume all areas have the same node density.');
end
if all(density_array==0)
    error 'All densities are zero.'
end

s = zeros(t_n,1); % area array
density_ratio = zeros(t_n,1);
counter_array = zeros(t_n,1);
%%
% Rand2Rand distance distribution within a single triangle
for i = 1:t_n % obtain triangles
    % vertexes (xv,yv) of a triangle
    xv = triangle_cell{i}(:,1)';
    yv = triangle_cell{i}(:,2)';
    s(i) = shoelace(xv,yv); % triangle area
end

for i = 1 : t_n
    density_ratio(i) = s(i)*density_array(i)/(s'*density_array);
    if i ~= 1
        density_ratio(i) = density_ratio(i) + density_ratio(i-1);
    end
end
%%
max_iter = 10000;
count = 0;
distance_in_sim=zeros(1,max_iter); % distance of two random points with a triangle

% x = triangle_cell{1}(:,1)';
% y = triangle_cell{1}(:,2)';
% for i = 2:t_n
%     x = [x triangle_cell{i}(:,1)'];
%     y = [y triangle_cell{i}(:,2)'];
% end
% x_length = max(x)-min(x);
% y_length = max(y)-min(y);


while count < max_iter
    % generage the 1st random point within the pentagon
    rand_density_ratio = rand;
    for i = 1:t_n
        if rand_density_ratio < density_ratio(i)
            xv = triangle_cell{i}(:,1)';
            yv = triangle_cell{i}(:,2)';
            x_length = max(xv)-min(xv);
            y_length = max(yv)-min(yv);
            while 1
                r_x1 = min(xv)+x_length*rand(1);
                r_y1 = min(yv)+y_length*rand(1);
                IN = inpolygon(r_x1,r_y1,xv,yv);
                if IN == 1
                    break;
                end
            end
            counter_array(i) = counter_array(i)+1;
            break;
        end
    end
    
    % generage the 2nd random point within the pentagon
    rand_density_ratio = rand;
    for i = 1:t_n
        if rand_density_ratio < density_ratio(i)
            xv = triangle_cell{i}(:,1)';
            yv = triangle_cell{i}(:,2)';
            x_length = max(xv)-min(xv);
            y_length = max(yv)-min(yv);
            while 1
                r_x2 = min(xv)+x_length*rand(1);
                r_y2 = min(yv)+y_length*rand(1);
                IN = inpolygon(r_x2,r_y2,xv,yv);
                if IN == 1
                    break;
                end
            end
            counter_array(i) = counter_array(i)+1;
            break;
        end
    end
%     hold on;
%     plot([r_x1 r_x2],[r_y1 r_y2],'k');
    dist = norm([r_x1 r_y1]-[r_x2 r_y2]);%sqrt((r_x1-r_x2)^2+(r_y1-r_y2)^2);
    count = count + 1;
    distance_in_sim(count)=dist;
end
sim_density = counter_array./s;

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