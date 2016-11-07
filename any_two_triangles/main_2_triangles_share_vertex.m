%%
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
% x1 = [D(1) Bp(1) C(1) ];
% y1 = [D(2) Bp(2) C(2) ];
% x2 = [D(1) Bp(1) B(1) ];
% y2 = [D(2) Bp(2) B(2) ];
% x1 = [A(1) B(1) Bp(1)];
% y1 = [A(2) B(2) Bp(2)];
% x2 = [D(1) 0.3 C(1)];
% y2 = [D(2) 0.2 C(2)];
% figure(fig);fig=fig+1;
% line([x1 x1(1)],[y1 y1(1)]);
% line([x2 x2(1)],[y2 y2(1)]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
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
%%
d_step = 1000;
t1 = t3;
t2 = t4;
% [ d_array, pdf_array, cdf_array ] = f_rand2rand_between_2_triangles( xv,yv,d_step );
[ d_array, pdf_array, cdf_array ] = f_rand2rand_between_any_2_triangles( t1,t2,d_step );
figure;
analysis = plot(d_array,cdf_array);
hold on;

[sim_d_array,sim_pdd_cdf] = f_sim_pdd_2_triangles(t1,t2);
% [d_array,pdd_cdf]=f_systematic_pdd_2_triangles_vertex_regular_pentagon(a,theta_array);
sim = plot(sim_d_array,sim_pdd_cdf,'r*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box on;
% axis([0 1.0 0 1.0]);
xlabel('Distance','fontsize',16);
ylabel('CDF','fontsize',16);

% AX = legend([lkm_based existing sim_based],'LKM-based','Existing Result','Simulation',4);
AX = legend([analysis sim],'LKM-based','Simulation','Location','SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',16);
grid on;