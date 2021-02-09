function G_matrix = Get_G_matrix(x)
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
    
    G_matrix=zeros(18,6);
    G_matrix(1,:) = -1*[1 0 0 0 -w v];
    G_matrix(2,:) = -1*[0 1 0 w 0 -u];
    G_matrix(3,:) = -1*[0 0 1 -v u 0];
    G_matrix(7,:) = -1*[0 0 0 1 sin(phi)*tan(theta) cos(phi)*tan(theta)];
    G_matrix(8,:) = -1*[0 0 0 0 cos(phi) -sin(phi)];
    G_matrix(9,:) = -1*[0 0 0 0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
end