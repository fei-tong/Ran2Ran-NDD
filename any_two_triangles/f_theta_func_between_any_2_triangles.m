function [ L ] = f_theta_func_between_any_2_triangles( theta,t1,t2 )
% This function is to get the lengths of 3 segments which are generated due
% to the intersection between a line with angle theta (with regard to
% x-axis) and the two triangles. 
% The function will be called by the following function:
% function [ d_array, pdf_array, cdf_array ] = f_rand2rand_between_any_2_triangles( t1,t2,d_step )
%
% Input:
%   theta: the angle of an arbitrary line passed the origin.
%   t1: triangle 1, t1 is like [x1 y1;x2 y2;x3 y3], where [xi yi] is a
%       vertex of t1.
%   t2: triangle 2, t2 is similar to t1.
% Output:
%   L: n*3, Lengths of three segments. n is determined by the step distance
%       between two adjacent parallel lines with angle of theta.
% Author: Fei Tong
% Date: May. 12, 2016

%% Move and rotate
% this section will be commented
% clear;clc;
% a = 1;
% A_1 = 30*pi/180;
% A_2 = 90*pi/180;
% A_3 = 110*pi/180;
% 
% % A_1 = 36*pi/180;
% % A_2 = 72*pi/180;
% % A_3 = 108*pi/180;
% % A_4 = (180-36)*pi/180;
% b = a/sqrt(2*(1-cos(A_3)));
% 
% A = [0 b*sin(A_3)]; D = [-b*cos(A_3) 0]; B = [-b*cos(A_3)+a*cos(A_2) a*sin(A_2)];
% C = [b-b*cos(A_3) 0]; Bp = [b-2*b*cos(A_3) a*sin(A_1)]; % Bp is B'
% B = [0.6 0.4];
% x = [D(1) C(1) Bp(1) B(1) A(1)];
% y = [D(2) C(2) Bp(2) B(2) A(2)];
% fig = 1;
% figure(fig);fig=fig+1;
% line([x x(1)],[y y(1)]);
% text(D(1),D(2),'D');text(Bp(1),Bp(2),'Bp');text(B(1),B(2),'B');text(C(1),C(2),'C');text(A(1),A(2),'A');
% 
% x1 = [Bp(1) B(1) D(1)];
% y1 = [Bp(2) B(2) D(2)];
% x2 = [A(1) B(1) D(1)];
% y2 = [A(2) B(2) D(2)];
% 
% % x1 = [D(1) Bp(1) C(1) ];
% % y1 = [D(2) Bp(2) C(2) ];
% % x2 = [D(1) Bp(1) B(1) ];
% % y2 = [D(2) Bp(2) B(2) ];
% x1 = [A(1) B(1) Bp(1)];
% y1 = [A(2) B(2) Bp(2)];
% x2 = [D(1) 0.3 C(1)];
% y2 = [D(2) 0.2 C(2)];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % vx = min([x1 x2]);
% % vy = min([y1 y2]);
% % 
% % M = [0 0] - [vx vy]; % moving vector
% % for i = 1:length(x1) % move all vertexes correspondingly
% %     x1(i) = x1(i) + M(1);
% %     x2(i) = x2(i) + M(1);
% %     y1(i) = y1(i) + M(2);
% %     y2(i) = y2(i) + M(2);
% % end
% % % figure(fig);fig=fig+1;
% % figure;
% % line([x1 x1(1)],[y1 y1(1)]);
% % line([x2,x2(1)],[y2,y2(1)]);

%%
% theta = 1.8850
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = t1(:,1)';
y1 = t1(:,2)';
x2 = t2(:,1)';
y2 = t2(:,2)';
epsilon = 1e-10;
[ A,B,~ ] = get_line( theta,[0 0], 'a'); % the line passing the origin
% Find the vertex, so the line with theta and passing this vertex is at the
% same side of all vertexes:
x=[x1 x2];
y=[y1 y2];
for i = 1:length(x)
    C = -A*x(i)-B*y(i);
    pn = 0; % positive or negative
    sameside = 1; % assume yes
    for j = 1:length(x)
        if i==j
            continue;
        end
        temp = A*x(j) + B*y(j) + C;
        if abs(temp) < epsilon % temp==0: the vertex [x(j) y(j)] is on the line
            continue;
        end
        if pn == 0
            pn = sign(temp);
        elseif pn*sign(temp) < 0
            sameside = 0; % not at the same side
            break;
        end
    end
    if sameside == 1 % have found the vertex
        vx = x(i);
        vy = y(i);
        break;
    end
end
% text(vx,vy,'v');
[ LA,LB,LC ] = get_line( theta,[vx vy], 'a'); % the line passing the [vx,vy]
% T1 and T2 contain the following 6 columns: triangle (1 or 2); vertex index;
% LA LB LC (a line LA*x+LB*y+LC=0); distance (between the line passing a vertex and the line passing [vx,vy])
col = 6;
row = length(x1);
T1 = zeros(row,col);
T2 = T1;
for i = 1:row
    T1(i,1) = 1; %triangle 1
    T1(i,2) = i; % vertex index
    % line: LA*x+LB*y+LC=0:
    T1(i,3) = LA; 
    T1(i,4) = LB;
    T1(i,5) = -LA*x1(i)-LB*y1(i);
    % distance (from [vx yx] to the above line):
    T1(i,6) = abs( LA*vx+LB*vy+T1(i,5) )/sqrt(LA^2+LB^2);
end

for i = 1:row
    T2(i,1) = 2; % triangle 2
    T2(i,2) = i; % vertex index
    % line: LA*x+LB*y+LC=0:
    T2(i,3) = LA;
    T2(i,4) = LB;
    T2(i,5) = -LA*x2(i)-LB*y2(i);
    % distance (from the origin to the above line):
    T2(i,6) = abs(LA*vx+LB*vy+T2(i,5))/sqrt(LA^2+LB^2);
end
% T = [T1;T2];
T = sortrows([T1;T2],col);
% Each triangle has a flag array, which indicates the number of vertexes
% the line has hitted:
row = size(T,1);
T_flag = zeros(row,2);
temp = zeros(1,2);
j = 1;
T_index = zeros(1,1);
T_edges = zeros(4,2);%cell(2,1); % store the 4 edges of triangle 1 and 2 intersected with line G
T_edges_cell = cell(1,1);
k = 1;
for i = 1:row
    t = T(i,1); % the triangle with changes
    if t == 1 % the first triangle has a change
        T_flag(i,1) = temp(1) + 1;
        T_flag(i,2) = temp(2);
        temp(1) = temp(1) + 1;
    else % t == 2, the second triangle has a change
        T_flag(i,2) = temp(2) + 1;
        T_flag(i,1) = temp(1);
        temp(2) = temp(2) + 1;
    end
    if temp(1) * temp(2) == 0 && temp(t) ~= 3
%         if isempty(T_edges{t})
        if all(all(T_edges==0)) 
            if T(i,2) == 1
                T_edges(t*2-1,:) = [1 2];
                T_edges(t*2,:) = [1 3];
            elseif T(i,2) == 2
                T_edges(t*2-1,:) = [2 1];
                T_edges(t*2,:) = [2 3];
            elseif T(i,2) == 3
%                 T_edges{t} = [3 1;3 2];
                T_edges(t*2-1,:) = [3 1];
                T_edges(t*2,:) = [3 2];
            end
        else
%             if T_edges{t}(1,2) == T(i,2)
            if T_edges(t*2-1,2) == T(i,2)
                T_edges(t*2-1,:) = [T_edges(t*2-1,2) T_edges(t*2,2)];
            elseif T_edges(t*2,2) == T(i,2)
                T_edges(t*2,:) = [T_edges(t*2-1,2) T_edges(t*2,2)];
            end
        end
    end
    if temp(1) * temp(2) ~= 0
        if temp(t) == 1
%             if T(i,2) == 1
%                 T_edges{t} = [1 2;1 3];
%             elseif T(i,2) == 2
%                 T_edges{t} = [2 1;2 3];
%             elseif T(i,2) == 3
%                 T_edges{t} = [3 1;3 2];
%             end
            
            if T(i,2) == 1
                T_edges(t*2-1,:) = [1 2];
                T_edges(t*2,:) = [1 3];
            elseif T(i,2) == 2
                T_edges(t*2-1,:) = [2 1];
                T_edges(t*2,:) = [2 3];
            elseif T(i,2) == 3
%                 T_edges{t} = [3 1;3 2];
                T_edges(t*2-1,:) = [3 1];
                T_edges(t*2,:) = [3 2];
            end
            T_edges_cell{k}=T_edges;
            k = k+1;
        elseif temp(t) == 2
            if T_edges(t*2-1,2) == T(i,2)
                T_edges(t*2-1,:) = [T_edges(t*2-1,2) T_edges(t*2,2)];
            elseif T_edges(t*2,2) == T(i,2)
                T_edges(t*2,:) = [T_edges(t*2-1,2) T_edges(t*2,2)];
            end
            T_edges_cell{k}=T_edges;
            k = k+1;
        end
        T_index(j) = i;
        j = j+1;
    end
    if temp(1) == 3 || temp(2) == 3
        break;
    end
end
% [T T_flag]
if T_index == 0 % the two triangles have no intersection
    L = 0;
    return;
end
if T_index == 1
    error 'Error: T_index cannot be 1.'
end
% T_index
% T_edges_cell{1}
% T_edges_cell{2}
% T_edges_cell{3}
% T_edges_cell{4}
% T_edges{1}
% T_edges{2}

p_min = T(T_index(1),6);
p_max = T(T_index(end),6);
delta_p=1/1000;
p_array = p_min:delta_p:p_max;
L = zeros(length(p_array),3);
inter_points = zeros(4,2); % intersection points with triangle edges
for i = 1:length(p_array)
    var_p = p_array(i);
    % first, get the line with a distance var_p to the line which passes v1
    GA = LA;
    GB = LB;
    %         GC = LC(1)+var_p*sign(LC(3)-LC(1));
    GC = LC + var_p*sign(T(end,5)-LC)*sqrt(GA^2+GB^2);
    for q = 2:length(T_index)
        T_edges = T_edges_cell{q-1};
        if var_p >= T(T_index(q-1),6) && var_p <= T(T_index(q),6)
            for j = 1:4 % calculate the 4 intersection points
                if j == 1 || j == 2 %triangle 1
                    k = 0;
                else % triangle 2
                    k = 3;
                end
                [ A,B,C ] = get_line( [x(T_edges(j,1)+k) y(T_edges(j,1)+k)],[x(T_edges(j,2)+k) y(T_edges(j,2)+k)], 'p');
                temp = [GA,GB;A,B]\[-GC;-C]; % calculate the intersection points
                inter_points(j,1) = temp(1);
                inter_points(j,2) = temp(2);
            end
            break;
        end
    end
%     inter_points
    if abs(theta-pi/2)<epsilon % theta == pi/2, sort inter_points according to y-axis
        inter_points = sortrows(inter_points,2);
    else % sort inter_points according to x-axis
        inter_points = sortrows(inter_points,1);
    end
    
    for j = 2:4
        L(i,j-1) = norm([inter_points(j,1) inter_points(j,2)]-[inter_points(j-1,1) inter_points(j-1,2)]);
    end
%     L
end

end