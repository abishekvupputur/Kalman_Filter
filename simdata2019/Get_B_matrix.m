function B_matrix = Get_B_matrix(x,U)
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
    
    B_matrix = zeros(18,6);
    B_matrix(1,:) = [1 0 0 0 -w v];
    B_matrix(2,:) = [0 1 0 w 0 -u];
    B_matrix(3,:) = [0 0 1 -v u 0];
    B_matrix(7,:) = [0 0 0 1 sin(phi)*tan(theta) cos(phi)*tan(theta)];
    B_matrix(8,:) = [0 0 0 0 cos(phi) -sin(phi)];
    B_matrix(9,:) = [0 0 0 0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
end
