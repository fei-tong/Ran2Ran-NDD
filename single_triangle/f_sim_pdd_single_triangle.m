% Simulation on the cdf of two random points within an arbitrary triangle
%
% Author: F. Tong
% Date: Nov. 13, 2013
% Last update: Nov. 13, 2013

function [d_array,pdd_cdf] = f_sim_pdd_single_triangle(xv,yv)
max_iter = 5000;
count = 0;
distance_in_sim=zeros(1,max_iter); % distance of two random points with a triangle

% the coordinate system is established following the way in the arxiv
% report, specifically,  Let ¡÷ ABC denote an arbitrary triangle with the 
% side lengths a, b, and c, the internal angles ¦Á, ¦Â, and ¦Ã, and the 
% altitudes h_a , h_b , and h_c , respectively (see an example shown in 
% Fig. 1). Without loss of generality, let a>=b>=c. Without loss of 
% generality, let vertex C be located at the origin, and side CB on the 
% positive x-axis.

x_length = max(xv)-min(xv);
y_length = max(yv)-min(yv);

while count < max_iter
    % generage the 1st random point
    while 1
        x1 = min(xv)+x_length*rand(1);
        y1 = min(yv)+y_length*rand(1); % h_a is the altitude from point A to edge BC
        IN = inpolygon(x1,y1,xv,yv);
        if IN == 1
            break;
        end
    end
    
    % generage the 2nd random point
    while 1
        x2 = min(xv)+x_length*rand(1);
        y2 = min(yv)+y_length*rand(1);
        IN = inpolygon(x2,y2,xv,yv);
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

% h = cdfplot(distance_in_sim);
% set( h, 'Color', 'r','LineStyle',':');

% [cc,x] = ecdf(distance_in_sim);% from simulation
% for i =1:50
% tick = int32(length(x)/50);
% plot( x(i*tick), cc(i*tick), '--rs','MarkerFaceColor','y','MarkerSize',3);
% end
% 
% box on;
% xlabel('Distance','fontsize',14);
% ylabel('CDF','fontsize',14);
% 
% AX = legend('Systematic Approach','Simulation');
% LEG = findobj(AX,'type','text');
% set(LEG,'FontSize',13);
% grid on;

end