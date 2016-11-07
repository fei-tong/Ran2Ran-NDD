% Call: f_integral_cld_cdf
% Called by: f_formula_cld_cdf
% call process: main -> formula -> specific d -> integral
% flag: indicate which case of three cases: 1. [0,Z]; 2. [Z, pi-Y]; 3.
%       [pi-Y, pi]
function pdd_pdf = f_specific_d_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag)

pdd_pdf = 0;
if flag == 1 % F1, case2: theta \in [0, Z]
%     lower = 0; upper = Z; 
%     f11
    if upper <= pi/2 - Y
        %f11 & f12
        if (b*sin(X))/d <= 1
            theta = asin((b*sin(X))/d) - Y;
        else
            theta = pi/2 -Y;
        end
        if theta < lower % d > b*sin(X)/sin(Y) && d > -b*sin(X)/sin(Y) %
            pdd_pdf =  0;
%             disp('error11');
        elseif  theta <= upper % d >= b*sin(X)/sin(Y+Z) && d >= - b*sin(X)/sin(Y+Z)%
            pdd_pdf = f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'11');
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'12');
        elseif theta > upper %l < b*sin(X)/sin(Y+Z) || l < -b*sin(X)/sin(Y+Z) %
            pdd_pdf = f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag,'11');
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag,'12');
        else
            disp('error12');
        end
    else % upper > pi/2 - Y
        %f12
        %divided into two parts£º [lower, pi/2-Y] ºÍ [pi/2-Y,upper]
        % part1: F121
%         theta = asin((b*sin(X))/d) - Y;
        if (b*sin(X))/d <= 1
            theta = asin((b*sin(X))/d) - Y;
        else
            theta = pi/2 -Y;
        end
        upper1 = pi/2-Y;
        if theta < lower % d > b*sin(X)/sin(Y) && d > -b*sin(X)/sin(Y) % % lower = 0; upper = Z; 
            pdd_pdf = 0;
%             disp('error21');
        elseif theta <= upper1 % d >= b*sin(X) || d >= -b*sin(X) 
            pdd_pdf = f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'11');
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'12');
        else
            disp('error22');
        end 
        % part2: F122
%         theta = pi - asin((b*sin(X))/d) - Y;
        if (b*sin(X))/d <= 1
            theta = pi - asin((b*sin(X))/d) - Y;
        else
            theta = pi/2 -Y;
        end
        if theta <= upper % l <= b % 
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'11');
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'12');
        elseif theta > upper
            pdd_pdf = pdd_pdf + 0;
        else
            disp('error23');
        end        
    end

elseif flag == 2 % F2, case4: theta \in [Z, pi-Y]
%     lower = Z; upper = pi-Y; 
    %divided into two parts£º [lower, pi/2] ºÍ [pi/2,upper]
    %part 1: F21
    if (c*sin(Y))/d<=1
        theta = asin((c*sin(Y))/d);
    else
        theta = pi/2;
    end
    upper1 = pi/2;

    if theta < lower
        pdd_pdf = 0;
    elseif theta <= upper1
        pdd_pdf =  f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'21');
        pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'22');
    else
        disp('error31');
    end
    %part 2: F22
    if (c*sin(Y))/d <= 1
        theta = pi - asin((c*sin(Y))/d);
    else
        theta = pi/2;
    end
    if theta <= upper
        pdd_pdf = pdd_pdf +  f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'21');
        pdd_pdf = pdd_pdf +  f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'22');
    elseif theta > upper
        pdd_pdf = pdd_pdf + 0;
    else
        disp('error32');
    end
elseif flag == 3 % F3, case6: theta \in [pi-Y, pi]
%     lower = pi-Y; upper = pi; 
    %     lower = pi-Y; upper = pi;
    if lower >= (pi/2+Z)
        %F31
        if (c*sin(X))/d <= 1
            theta = pi - asin((c*sin(X))/d) + Z;
        else
            theta = pi/2+Z;
        end
        if theta < lower
            pdd_pdf = f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag,'31');
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,upper,flag,'32');
        elseif theta <= upper
            pdd_pdf = f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'31');
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'32');
        elseif theta > upper
            pdd_pdf = 0;
        else
            disp('error34');
        end
    else % lower < pi/2+Z
        %F32
        %divided into two parts: [lower, pi/2+Z] ºÍ [pi/2+Z,upper]
        %part 1: F321
        if (c*sin(X))/d <= 1
            theta = asin((c*sin(X))/d) + Z;
        else
            theta = pi/2 + Z;
        end
        upper1 = pi/2+Z;
        if theta < lower
            pdd_pdf = 0;
        elseif theta <= upper1
            pdd_pdf = f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'31');
            pdd_pdf = pdd_pdf + f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,lower,theta,flag,'32');
        else
            disp('error35');
        end        
        %part 2: F322
        if (c*sin(X))/d<=1
            theta = pi - asin((c*sin(X))/d) + Z;
        else
            theta = pi/2 + Z;
        end
        if theta <= upper
            pdd_pdf = pdd_pdf +  f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'31');
            pdd_pdf = pdd_pdf +  f_integral_pdd_pdf_triangle(X,Y,Z,a,b,c,d,theta,upper,flag,'32');
        else
            disp('error37');
        end        
    end    
else
    disp('error0');
end

end