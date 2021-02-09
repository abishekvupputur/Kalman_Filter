%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zpred = kf_calcHx(x) Calculates the output dynamics equation h(x,u,t) 
%   
%   Author: C.C. de Visser, Delft University of Technology, 2013
%   email: c.c.devisser@tudelft.nl
%   Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zpred = kf_calcHx(t, x, U)
   u = x(1);
    v = x(2);
    w = x(3);
    px = x(4);
    py = x(5);
    pz = x(6);
    phi = x(7);
    theta = x(8);
    psi = x(9);
    lambda_x = x(10);
    lambda_y = x(11);
    lambda_z = x(12);
    lambda_p = x(13);
    lambda_q = x(14);
    lambda_r = x(15);
    u_wind = x(16);
    v_wind = x(17);
    w_wind = x(18);
    
    Ax = U(1);
    Ay = U(2);
    Az = U(3);
    p = U(4);
    q = U(5);
    r = U(6);
    
    g=9.81;


    zpred(1) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi) - w*sin(phi))*sin(psi) + u_wind;
    zpred(2) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi) - w*sin(phi))*cos(psi) + v_wind;
    zpred(3) = -u*sin(theta) + (v*sin(phi) + w*cos(phi))*cos(theta) + w_wind;
    zpred(4:9) = x(4:9);
    zpred(10) = sqrt(u^2 + v^2 + w^2);
    zpred(11) = atan2(w,u);
    zpred(12) = atan2(v,sqrt(u^2 + w^2));
    zpred=zpred';
    end
    