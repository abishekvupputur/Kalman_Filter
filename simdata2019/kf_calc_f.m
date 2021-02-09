%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xdot = kf_calcFx(x) Calculates the system dynamics equation f(x,u,t) 
%   
%   Author: C.C. de Visser, Delft University of Technology, 2013
%   email: c.c.devisser@tudelft.nl
%   Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = kf_calcFx(t, x, U)
    %%Function F(x)
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

    xdot(1) = (Ax - lambda_x) -g*sin(theta) + (r - lambda_r)*v - (q - lambda_q)*w;
    xdot(2) = (Ay - lambda_y) +g*cos(theta)*sin(phi) + (p - lambda_p)*w - (r - lambda_r)*u;
    xdot(3) = (Az - lambda_z) +g*cos(theta)*cos(phi) + (q - lambda_q)*u - (p - lambda_p)*v;
    xdot(4) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi) - w*sin(phi))*sin(psi) + u_wind;
    xdot(5) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi) - w*sin(phi))*cos(psi) + v_wind;
    xdot(6) = -u*sin(theta) + (v*sin(phi) + w*cos(phi))*cos(theta) + w_wind;
    xdot(7) = (p - lambda_p) + (q - lambda_q)*sin(phi)*tan(theta) + (r - lambda_r)*cos(phi)*tan(theta);
    xdot(8) = (q - lambda_q)*cos(phi) - (r - lambda_r)*sin(phi);
    xdot(9) = (q - lambda_q)*sin(phi)/cos(theta) + (r - lambda_r)*cos(phi)/cos(theta);
    xdot(10) = 0;
    xdot(11) = 0;
    xdot(12) = 0;
    xdot(13) = 0;
    xdot(14) = 0;
    xdot(15) = 0;
    xdot(16) = 0;
    xdot(17) = 0;
    xdot(18) = 0;
    
    xdot =xdot';
    end
