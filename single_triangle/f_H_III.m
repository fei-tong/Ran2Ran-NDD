% calculate H^I_1 and H^I_2
function I_value = f_H_III(X,Y,Z,a,b,c,d,s,theta,flag)
    if flag == 1
        I_value = (a*d*((d^2*sin(Z - 2*theta))/4 + (d^2*theta*cos(Z))/2 + 2*c*d*sin(X)*cos(theta) + c^2*theta*cos(Z)*sin(X)^2 + c^2*sin(Z)*log(-sin(Z - theta))*sin(X)^2))/(c*sin(X));
    elseif flag == 2
        I_value = (2*d*((d^2*sin(Y - Z + 2*theta))/8 - (theta*cos(Y + Z)*(2*c^2*sin(X)^2 + d^2))/4 - c*d*cos(Y + theta)*sin(X) - (c^2*log(sin(Z - theta))*sin(Y + Z)*sin(X)^2)/2))/sin(X);
    end
    I_value = I_value/s^2;
end