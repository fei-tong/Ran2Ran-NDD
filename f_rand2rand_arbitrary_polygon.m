function [ d_array, cdf_array ] = f_rand2rand_arbitrary_polygon( triangle_cell,varargin )
% This function is to get the numerical Rand2Rand distance distribution
%   within an arbitrary polygon.
%
% Input:
%   triangle_cell: Triangulated triangules of an arbitrary polygon. Each
%       cell element contains a triangle, which is like
%       [x1 y1;x2 y2;x3 y3], where [xi yi] is a vertex of the triangle.
%
%   varargin contains the following one optional argument:
%       -> density array: (n-2)*1, containing "node densities" of each
%          triangle.
% Output:
%   [ d_array, cdf_array ]: Rand2Rand distance distribution within the polygon
%
% Author: Fei Tong
% Date: May. 12, 2016

if nargin < 1
    error 'This function needs at least 1 arguments.'
end

t_n = length(triangle_cell); % number of triangles: totally, there are n-2 triangles
if isempty(varargin{1})
    density_array = ones(t_n,1);
else
    density_array = varargin{1};
end
[q1,q2] = size(density_array);
if q1 == 1 && q2 > 1
    density_array = density_array'; % change to n*1 from 1*n
end
q = length(density_array);
if q ~= t_n
    density_array = ones(t_n,1);
    disp('The size of density_array is not right. We will assume all areas have the same node density.');
end
if all(density_array==0)
    error 'All densities are zero.'
end

cdf_cell = cell(t_n,t_n); % the CDF results for each triangle
d_step = 1000;
%% prepare d_array and pdf_array
x = triangle_cell{1}(:,1)';
y = triangle_cell{1}(:,2)';
for i = 2:t_n
    x = [x triangle_cell{i}(:,1)'];
    y = [y triangle_cell{i}(:,2)'];
end
max_d = 0;
n=length(x);
for i = 1:n
    for j = i+1:n
        d = norm([x(i) y(i)]-[x(j) y(j)]);
        if max_d < d
            max_d = d;
        end
    end
end
% d_step = 1000;
delta_d = 1/d_step; 
d_array = 0:delta_d:max_d;
cdf_array = zeros(1,length(d_array));
s = zeros(t_n,1); % area array
%%
% Rand2Rand distance distribution within a single triangle
for i = 1:t_n % obtain triangles
    xv = triangle_cell{i}(:,1)';
    yv = triangle_cell{i}(:,2)';
    s(i) = shoelace(xv,yv); % triangle area
    [~,~,pdd_cdf] = f_formula_pdd_pdf_triangle(xv,yv,d_step);
    cdf_cell{i,i}=pdd_cdf;
end

% Rand2Rand distance distribution between two triangles
for i = 1:t_n-1
    for j = i+1:t_n
        % vertexes of a 5-gon, used to calculate the distance distribution
        % between two triangles:
%         xv = [triangle_cell{i}(1,:) triangle_cell{j}(1,2:end)];
%         yv = [triangle_cell{i}(2,:) triangle_cell{j}(2,2:end)];
        
        [ ~, ~, pdd_cdf ] = f_rand2rand_between_any_2_triangles( triangle_cell{i},triangle_cell{j},d_step );
        cdf_cell{i,j} = pdd_cdf;
        cdf_cell{j,i} = pdd_cdf;
    end
end

% probabilistic sum
F = zeros(t_n,t_n);
for k = 1:length(d_array)
    for i = 1:t_n
        for j = i:t_n
            if k > length(cdf_cell{i,j})
                F(i,j) = 1;
                F(j,i) = 1;
            else
                F(i,j) = cdf_cell{i,j}(k);
                F(j,i) = cdf_cell{i,j}(k);
            end
        end
    end
    
    for i = 1:t_n
        for j = 1:t_n
%             cdf_array(k) = cdf_array(k) + triangle_cell{i,2}/s*triangle_cell{j,2}/s*F(i,j);
            cdf_array(k) = cdf_array(k) + s(i)*density_array(i)/(s'*density_array) * s(j)*density_array(j)/(s'*density_array)*F(i,j);
        end
    end
end


end

