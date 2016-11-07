% Called by: f_cld_cdf_for_specific_l
% call process: main -> formula -> specific l -> integral
% fg_flag: intigate whether the l is smaller than the base, if fg_flag ==
%           'f', l<= base, if fg_flag='p', l>base
function pdd_pdf = f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag,fg_flag)

u = a+b+c;
s = sqrt(u/2*(u/2-a)*(u/2-b)*(u/2-c));

if flag == 1 % case2: %theta \in [0, Z]
    if fg_flag == '11'
%         pdd_pdf_u = -(1/2)*d*(-2*sin(X)^2*b^2*(sin(Y)*cos(Z)+cos(Y)*sin(Z))*log(tan(Z-upper)*cos(Y)*cos(Z)-tan(Z-upper)*sin(Y)*sin(Z)-sin(Y)*cos(Z)-cos(Y)*sin(Z))+sin(X)^2*b^2*(sin(Y)*cos(Z)+cos(Y)*sin(Z))*log(tan(Z-upper)^2+1)-(1/2)*d^2*sin(2*upper+Y-Z)+2*sin(X)^2*b^2*(-cos(Y)*cos(Z)+sin(Y)*sin(Z))*atan(tan(Z-upper))+d*(upper*cos(Y+Z)*d+4*sin(X)*cos(Z-upper)*b))/sin(X);
%         pdd_pdf_l = -(1/2)*d*(-2*sin(X)^2*b^2*(sin(Y)*cos(Z)+cos(Y)*sin(Z))*log(tan(Z-lower)*cos(Y)*cos(Z)-tan(Z-lower)*sin(Y)*sin(Z)-sin(Y)*cos(Z)-cos(Y)*sin(Z))+sin(X)^2*b^2*(sin(Y)*cos(Z)+cos(Y)*sin(Z))*log(tan(Z-lower)^2+1)-(1/2)*d^2*sin(2*lower+Y-Z)+2*sin(X)^2*b^2*(-cos(Y)*cos(Z)+sin(Y)*sin(Z))*atan(tan(Z-lower))+d*(lower*cos(Y+Z)*d+4*sin(X)*cos(Z-lower)*b))/sin(X);
%         pdd_pdf = (pdd_pdf_u-pdd_pdf_l)/s^2;
        pdd_pdf = f_H_I(X,Y,Z,a,b,c,d,s,upper,1)-f_H_I(X,Y,Z,a,b,c,d,s,lower,1);
    elseif fg_flag == '12'
%         pdd_pdf_u = d*a*(b^2*sin(Y)*(-1+cos(X))*(cos(X)+1)*log(sin(upper)*cos(Y)+cos(upper)*sin(Y))+(b^2-b^2*cos(X)^2)*cos(Y)*atan(sin(upper)/cos(upper))+2*d*cos(upper)*b*sin(X)-(1/4)*d^2*sin(2*upper+Y)+(1/2)*d^2*upper*cos(Y))/(sin(X)*b);
%         pdd_pdf_l = d*a*(b^2*sin(Y)*(-1+cos(X))*(cos(X)+1)*log(sin(lower)*cos(Y)+cos(lower)*sin(Y))+(b^2-b^2*cos(X)^2)*cos(Y)*atan(sin(lower)/cos(lower))+2*d*cos(lower)*b*sin(X)-(1/4)*d^2*sin(2*lower+Y)+(1/2)*d^2*lower*cos(Y))/(sin(X)*b);
%         pdd_pdf = (pdd_pdf_u-pdd_pdf_l)/s^2;
        pdd_pdf = f_H_I(X,Y,Z,a,b,c,d,s,upper,2)-f_H_I(X,Y,Z,a,b,c,d,s,lower,2);
    else
        disp('ER');
    end
elseif flag == 2 % case4: theta \in [Z, pi-Y]
    if fg_flag == '21'
        pdd_pdf = f_H_II(X,Y,Z,a,b,c,d,s,upper,1)-f_H_II(X,Y,Z,a,b,c,d,s,lower,1);
%         pdd_pdf_u = d*c*sin(Y)*b*cos(Z)*upper-d*c*sin(Y)*b*sin(Z)*log(sin(upper))+(1/4)*d^3*b*sin(-2*upper+Z)/(c*sin(Y))+(1/2)*d^3*b*cos(Z)*upper/(c*sin(Y))+2*d^2*b*cos(-upper+Z);
%         pdd_pdf_l = d*c*sin(Y)*b*cos(Z)*lower-d*c*sin(Y)*b*sin(Z)*log(sin(lower))+(1/4)*d^3*b*sin(-2*lower+Z)/(c*sin(Y))+(1/2)*d^3*b*cos(Z)*lower/(c*sin(Y))+2*d^2*b*cos(-lower+Z);
%         pdd_pdf = (pdd_pdf_u-pdd_pdf_l)/s^2;
    elseif fg_flag == '22'
        pdd_pdf = f_H_II(X,Y,Z,a,b,c,d,s,upper,2)-f_H_II(X,Y,Z,a,b,c,d,s,lower,2);
%         pdd_pdf_u = d*sin(Y)*c^2*cos(Y)*upper+d*sin(Y)^2*c^2*log(sin(upper))+2*d^2*c*cos(upper+Y)+(1/2)*d^3*cos(Y)*upper/sin(Y)-(1/4)*d^3*sin(2*upper+Y)/sin(Y);
%         pdd_pdf_l = d*sin(Y)*c^2*cos(Y)*lower+d*sin(Y)^2*c^2*log(sin(lower))+2*d^2*c*cos(lower+Y)+(1/2)*d^3*cos(Y)*lower/sin(Y)-(1/4)*d^3*sin(2*lower+Y)/sin(Y);
%         pdd_pdf = (pdd_pdf_u-pdd_pdf_l)/s^2;
    else
        disp('ER');
    end
elseif flag == 3 % case6: theta \in [pi-Y, pi]
    if fg_flag == '31'
        pdd_pdf = f_H_III(X,Y,Z,a,b,c,d,s,upper,1)-f_H_III(X,Y,Z,a,b,c,d,s,lower,1);
%         pdd_pdf_u = -(c^2*sin(Z)*(cos(X)-1)*(cos(X)+1)*log(sin(upper)*cos(Z)-sin(Z)*cos(upper))+c^2*cos(Z)*(cos(X)-1)*(cos(X)+1)*atan(sin(upper)/cos(upper))-(1/4)*d^2*sin(-2*upper+Z)-(1/2)*d^2*upper*cos(Z)-2*d*cos(upper)*c*sin(X))*d*a/(c*sin(X));
%         pdd_pdf_l = -(c^2*sin(Z)*(cos(X)-1)*(cos(X)+1)*log(sin(lower)*cos(Z)-sin(Z)*cos(lower))+c^2*cos(Z)*(cos(X)-1)*(cos(X)+1)*atan(sin(lower)/cos(lower))-(1/4)*d^2*sin(-2*lower+Z)-(1/2)*d^2*lower*cos(Z)-2*d*cos(lower)*c*sin(X))*d*a/(c*sin(X));
%         pdd_pdf = (pdd_pdf_u-pdd_pdf_l)/s^2;
    elseif fg_flag == '32'
        pdd_pdf = f_H_III(X,Y,Z,a,b,c,d,s,upper,2)-f_H_III(X,Y,Z,a,b,c,d,s,lower,2);
%         pdd_pdf_u = -2*d*(-(1/8)*d^2*sin(2*upper-Z+Y)-(1/2)*c^2*(cos(X)-1)*(cos(X)+1)*(cos(Y)*sin(Z)+cos(Z)*sin(Y))*log(sin(-upper+Z))+d*c*cos(upper+Y)*sin(X)+(1/2)*upper*((1/2)*cos(Z+Y)*d^2+c^2*(cos(X)-1)*(cos(X)+1)*(-cos(Z)*cos(Y)+sin(Z)*sin(Y))))/sin(X);
%         pdd_pdf_l = -2*d*(-(1/8)*d^2*sin(2*lower-Z+Y)-(1/2)*c^2*(cos(X)-1)*(cos(X)+1)*(cos(Y)*sin(Z)+cos(Z)*sin(Y))*log(sin(-lower+Z))+d*c*cos(lower+Y)*sin(X)+(1/2)*lower*((1/2)*cos(Z+Y)*d^2+c^2*(cos(X)-1)*(cos(X)+1)*(-cos(Z)*cos(Y)+sin(Z)*sin(Y))))/sin(X);
%         pdd_pdf = (pdd_pdf_u-pdd_pdf_l)/s^2;
    else
        disp('ER');
    end
end    
end