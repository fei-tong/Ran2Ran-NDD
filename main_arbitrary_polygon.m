% Author: Fei Tong
% Date: May. 10, 2016

%% Prepare polygon vertexes
%   NOTE THAT:
%       1. The points in (x,y) are sorted in a counterclockwise order.
%       2. The polygon must can be triangulated starting from a single
%          vertex (After triangulation, all triangles are within the polygon).
clear;
clc;
%%
a = 1;
A_1 = 30*pi/180;
A_2 = 90*pi/180;
A_3 = 110*pi/180;

% A_1 = 36*pi/180;
% A_2 = 72*pi/180;
% A_3 = 108*pi/180;
% A_4 = (180-36)*pi/180;
b = a/sqrt(2*(1-cos(A_3)));

A = [0 b*sin(A_3)]; D = [-b*cos(A_3) 0]; B = [-b*cos(A_3)+a*cos(A_2) a*sin(A_2)];
C = [b-b*cos(A_3) 0]; Bp = [b-2*b*cos(A_3) a*sin(A_1)]; % Bp is B'
B = [0.6 0.4];
% D = [0 0]; C=[2 0]; Bp=[1 1];B=[-1 1];A=[-2 0];
x = [D(1) C(1) Bp(1) B(1) A(1)];
y = [D(2) C(2) Bp(2) B(2) A(2)];

F=Bp;
B = [0.4 0.4];
E = [0.6 0.3];
% D = [0 0]; C=[2 0]; Bp=[1 1];B=[-1 1];A=[-2 0];
x = [D(1) E(1) C(1) Bp(1) B(1) A(1)];
y = [D(2) E(2) C(2) Bp(2) B(2) A(2)];

% x = [D(1) C(1) Bp(1) B(1)];
% y = [D(2) C(2) Bp(2) B(2)];
% x = [C(1) Bp(1) B(1) A(1)];
% y = [C(2) Bp(2) B(2) A(2)];
% x = [ D(1) C(1) Bp(1) A(1)];
% y = [ D(2) C(2) Bp(2) A(2) ];
% 
% x = [A(1) D(1) C(1) Bp(1) ];
% y = [A(2) D(2) C(2) Bp(2) ];
% 
% x = [Bp(1) A(1) D(1) C(1) ];
% y = [Bp(2) A(2) D(2) C(2) ];
% x = [C(1) Bp(1) B(1)];
% y = [C(2) Bp(2) B(2)];
line([x x(1)],[y y(1)]);
% regular pentagon
text(D(1),D(2),'D');text(Bp(1),Bp(2),'F');text(B(1),B(2),'B');text(C(1),C(2),'C');text(A(1),A(2),'A');text(E(1),E(2),'E');
line([B(1) D(1)],[B(2) D(2)],'LineStyle','-.','color','r');
line([B(1) E(1)],[B(2) E(2)],'LineStyle','-.','color','r');
line([F(1) E(1)],[F(2) E(2)],'LineStyle','-.','color','r');
t1=[A;B;D];
t2=[B;E;D];
t3=[B;E;F];
t4=[E;F;C];
t_cell = cell(1,1);
t_cell{1}=t1;
t_cell{2}=t2;
t_cell{3}=t3;
t_cell{4}=t4;

% A=A-D;
% B=B-D;
% C=C-D;
% Bp=Bp-D;
% D=D-D;
% E=C;
% C=Bp;
% C = [(C(1)+D(1))/3 (C(2)+D(2))/3];
% E =C;
% C = B;
% B = [0.2 0.2];
% C= B;

% x = [D(1) E(1) C(1) B(1) A(1)];
% y = [D(2) E(2) C(2) B(2) A(2)];
% x = [D(1) E(1) C(1) C(1) B(1)];
% y = [D(2) E(2) C(2) C(2) B(2)];
% x = [D(1) C(1) B(1) B(1) A(1)];
% y = [D(2) C(2) B(2) B(2) A(2)];
% line([x x(1)],[y y(1)]);
% line([D(1) B(1)],[D(2) B(2)]);
% line([D(1) C(1)],[D(2) C(2)]);
% text(D(1),D(2),'D');text(E(1),E(2),'E');text(B(1),B(2),'B');text(C(1),C(2),'C');text(A(1),A(2),'A');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% approach and simulation

% the same density
d1 = 1;
d2 = 1;
d3 = 1;
d4 = 1;
density_array = [d1 d2 d3 d4];
[ d_array, cdf_array ] = f_rand2rand_arbitrary_polygon( t_cell,density_array );
figure;
analysis = plot(d_array,cdf_array);

hold on;
[sim_d_array,sim_pdd_cdf,sim_density] = f_sim_rand2rand_arbitrary_polygon( t_cell,density_array );
sim_density
% [d_array,pdd_cdf]=f_systematic_pdd_2_triangles_vertex_regular_pentagon(a,theta_array);
sim1 = plot(sim_d_array,sim_pdd_cdf,'r*');

%% approach and simulation

% node density ratio: 1:5:10
d1 = 1;
d2 = 5;
d3 = 10;
d4 = 15;
density_array = [d1 d2 d3 d4];
[ d_array, cdf_array ] = f_rand2rand_arbitrary_polygon( t_cell,density_array );
analysis = plot(d_array,cdf_array);

hold on;
[sim_d_array,sim_pdd_cdf,sim_density] = f_sim_rand2rand_arbitrary_polygon( t_cell,density_array );
sim_density
% [d_array,pdd_cdf]=f_systematic_pdd_2_triangles_vertex_regular_pentagon(a,theta_array);
sim2 = plot(sim_d_array,sim_pdd_cdf,'k*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box on;
% axis([0 1.0 0 1.0]);
xlabel('Distance','fontsize',16);
ylabel('CDF','fontsize',16);

% AX = legend([lkm_based existing sim_based],'LKM-based','Existing Result','Simulation',4);
AX = legend([sim1 sim2],'density ratio: 1:1:1:1','density ratio: 1:5:10:15',4);
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',16);
grid on;