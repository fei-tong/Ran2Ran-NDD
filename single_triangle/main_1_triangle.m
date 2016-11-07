%%
clear;clc;
format long
interval = 20;
lkm = 'b-';
cld = 'k+';
sim = 'r*';

a = 1;
A_1 = 30*pi/180;
A_2 = 90*pi/180;
A_3 = 110*pi/180;
% A_4 = (180-36)*pi/180;
b = a/sqrt(2*(1-cos(A_3)));

A = [0 b*sin(A_3)]; D = [-b*cos(A_3) 0]; B = [-b*cos(A_3)+a*cos(A_2) a*sin(A_2)];
C = [b-b*cos(A_3) 0]; Bp = [b-2*b*cos(A_3) a*sin(A_1)]; % Bp is B'
% A=A-D;
% B=B-D;
% C=C-D;
% Bp=Bp-D;
% D=D-D;
E=C;
C=Bp;
% C = [(C(1)+D(1))/3 (C(2)+D(2))/3];
% E =C;
% C = B;
% B = [0.2 0.2];
% C= B;

% x = [D(1) E(1) C(1) B(1) A(1)];
% y = [D(2) E(2) C(2) B(2) A(2)];
% x = [D(1) E(1) C(1) C(1) B(1)];
% y = [D(2) E(2) C(2) C(2) B(2)];
x = [D(1) C(1) B(1) B(1) A(1)];
y = [D(2) C(2) B(2) B(2) A(2)];
line([x x(1)],[y y(1)]);
line([D(1) B(1)],[D(2) B(2)]);
line([D(1) C(1)],[D(2) C(2)]);
text(D(1),D(2),'D');text(E(1),E(2),'E');text(B(1),B(2),'B');text(C(1),C(2),'C');text(A(1),A(2),'A');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = [D(1) C(1) B(1)];
y = [D(2) C(2) B(2)];

% x = [D(1) B(1) A(1)];
% y = [D(2) B(2) A(2)];
d_step = 1000;
[d_a,~,pdd_cdf] = f_formula_pdd_pdf_triangle(x,y,d_step);
lkm_based = plot(d_a,pdd_cdf,lkm);
hold on;
% % random pdd simulation 
[d_array,pdd_cdf] = f_sim_pdd_single_triangle(x,y);
sim_based=plot( d_array, pdd_cdf, sim);