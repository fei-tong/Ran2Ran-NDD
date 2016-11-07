% calculate H^I_1 and H^I_2
function I_value = f_H_II(X,Y,Z,a,b,c,d,s,theta,flag)
    if flag == 1
        I_value = (b*d*(d^2*sin(Z - 2*theta) + 2*d^2*theta*cos(Z) - 4*c^2*log(sin(theta))*sin(Y)^2*sin(Z) + 4*c^2*theta*cos(Z)*sin(Y)^2 + 8*c*d*sin(Y)*cos(Z - theta)))/(4*c*sin(Y));
    elseif flag == 2
        I_value = (d*(2*d^2*theta*cos(Y) - d^2*sin(Y + 2*theta) + 4*c^2*log(sin(theta))*sin(Y)^3 + 4*c^2*theta*cos(Y)*sin(Y)^2 + 8*c*d*cos(Y + theta)*sin(Y)))/(4*sin(Y));
    end
    I_value = I_value/s^2;
end