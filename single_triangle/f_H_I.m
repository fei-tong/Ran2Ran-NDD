% calculate H^I_1 and H^I_2
function I_value = f_H_I(X,Y,Z,a,b,c,d,s,theta,flag)
    if flag == 1
        I_value = (d*((d^2*sin(Y - Z + 2*theta))/2 - d*(4*b*sin(X)*cos(Z - theta) + d*theta*cos(Y + Z)) + (b^2*log(-sin(Y + theta)/cos(Z - theta))*(2*sin(Y + Z) - sin(2*X + Y + Z) + sin(2*X - Y - Z)))/2 + 2*b^2*atan(tan(Z - theta))*cos(Y + Z)*sin(X)^2 - b^2*log(tan(Z - theta)^2 + 1)*sin(Y + Z)*sin(X)^2))/(2*sin(X));
    elseif flag == 2
        I_value = (a*d*((d^2*theta*cos(Y))/2 - (d^2*sin(Y + 2*theta))/4 - b^2*atan(tan(theta))*cos(Y)*(cos(X)^2 - 1) + 2*b*d*sin(X)*cos(theta) + b^2*log(sin(Y + theta))*sin(Y)*(cos(X) - 1)*(cos(X) + 1)))/(b*sin(X));
    end
    I_value = I_value/s^2;
end