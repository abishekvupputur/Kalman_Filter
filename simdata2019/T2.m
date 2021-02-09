%% Task 2
clc;
%% 2.1. AIRDATA_ESTIMATES
Files = ["da3211.mat" ; "de3211.mat" ; "dadoublet.mat" ; "drdoublet.mat"; "dr3211.mat"];
[a, ~] = size(Files);
for i=1:dataset
    load(Files(i))
    display("loading file "+Files(i));
    display(string);
    u=u_n + 1*10; %X direction Airvelocity corrected for zero wind and in Flat earth.
    v=v_n + 1*6; %Y Direction Airvelocity corrected for zero wind and in Flat earth.
    w=w_n + 1*1; %Z Direction Airvelocity corrected for zero wind and in Flat earth.
    t=t; %Log of time
    x=cumsum(u)*0.01;%cumtrapz(u,t);
    y=cumsum(v)*0.01;%cumtrapz(v,t);
    z=cumsum(w)*0.01;%cumtrapz(w,t);
%     %% TASK 2.1
%      figure;
%      subplot(3,1,1)
%      plot(x)
%      grid on
%      title(strcat("Flight path Generated with Airdata Estimates, Manoeuvre: ",Files(i)));
%      hold off
%      ylabel('x-axis [m]')
%      xlabel('Dataset')
%      xlim([0 12000])
%      
%      subplot(3,1,2)
%      plot(y)
%      grid on
%      title(strcat("Flight path Generated with Airdata Estimates, Manoeuvre: ",Files(i)));
%      hold off
%      ylabel('y-axis [m]')
%      xlabel('Dataset')
%      xlim([0 12000])
%      
%      subplot(3,1,3)
%      plot(z)
%      grid on
%      title(strcat("Flight path Generated with Airdata Estimates, Manoeuvre: ",Files(i)));
%      hold off
%      ylabel('z-axis [m]')
%      xlabel('Dataset')
%      xlim([0 12000])
%      
%      set(gcf,'units','points','position',[10,10,1600,800])
%     File_names = ["da3211" ; "de3211" ; "dadoublet" ; "drdoublet"; "dr3211"];
%     saveas(gcf,strcat('FlightPath_',File_names(dataset),'.png'));
% end

%% 2.2. Generate GPS, IMU, Airdata

u_w = 10; %x-direction wind speed m/s
v_w = 6; %y-direction wind speed m/s
w_w = 1; %z-direction wind speed m/s

W_xe = 0;
W_ye = 0;
W_ze = 0;

bias_x = 0.001; %m/s2
bias_y = 0.001; %m/s2
bias_z = 0.001; %m/s2

bias_p = deg2rad(0.001); %rad/s
bias_q = deg2rad(0.001); %rad/s
bias_r = deg2rad(0.001); %rad/s

%%%%% IMU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ax_imu = Ax + bias_x + 0.001*randn(size(Ax,1),1);
Ay_imu = Ay + bias_y + 0.001*randn(size(Ay,1),1);
Az_imu = Az + bias_z + 0.001*randn(size(Az,1),1);

p_imu = p + bias_p + deg2rad(0.001)*randn(size(p,1),1);
q_imu = q + bias_q + deg2rad(0.001)*randn(size(q,1),1);
r_imu = r + bias_r + deg2rad(0.001)*randn(size(r,1),1);

%%%% GPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_gps = x + 10*randn(size(x,1),1);
y_gps = y + 10*randn(size(y,1),1);
z_gps = z + 10*randn(size(z,1),1);

u_gps = u + 0.1*randn(size(u,1),1);
v_gps = v + 0.1*randn(size(v,1),1);
w_gps = w + 0.1*randn(size(w,1),1);

phi_gps = phi + deg2rad(0.1)*randn(size(phi,1),1);
theta_gps = theta + deg2rad(0.1)*randn(size(theta,1),1);
psi_gps = psi + deg2rad(0.1)*randn(size(psi,1),1);

%%%% Airdata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_b = u_n.*cos(theta).*cos(psi) + v_n.*sin(psi).*cos(theta) - w_n.*sin(theta);
v_b = u_n.*(sin(phi).*sin(theta).*cos(psi)-cos(phi).*sin(psi))+ ...
        v_n.*(sin(phi).*sin(theta).*sin(psi)+cos(phi).*cos(psi))+ ...
        w_n.*sin(phi).*cos(theta);
w_b = u_n.*(cos(phi).*sin(theta).*cos(psi)+sin(phi).*sin(psi))+ ...
        v_n.*(cos(phi).*sin(theta).*sin(psi)-sin(phi).*cos(psi))+ ...
        w_n.*cos(phi).*cos(theta);

u_ardata = 1*u_w + u_n;
v_ardata = 1*v_w + v_n;
w_ardata = 1*w_w + w_n;

V_ardta = vtas +  0.1*randn(size(u,1),1); 
alpha_ardta = alpha + deg2rad(0.1)*randn(size(alpha,1),1);
beta_ardta = beta + deg2rad(0.1)*randn(size(beta,1),1);

%% Task 2.2
    %GPS DATA PLOTS
% figure;
%     subplot(331)
%     plot(t, x_gps)
%     grid on
%     hold on
%     plot(t, x)
%     hold off
%     title(strcat("GPS: x-axis, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Integrated Velocity, Nav frame');
%     ylabel('x-axis [m]')
%     xlabel('time [s]')
%     
%     subplot(332)
%     plot(t, y_gps)
%     grid on
%     hold on
%     plot(t, y)
%     hold off
%     title(strcat("GPS: y-axis, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Integrated Velocity, Nav frame');
%     ylabel('y-axis [m]')
%     xlabel('time [s]')
%     
%     subplot(333)
%     plot(t, z_gps)
%     grid on
%     hold on
%     plot(t, z)
%     hold off
%     title(strcat("GPS: z-axis, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Integrated Velocity, Nav frame');
%     ylabel('z-axis [m]')
%     xlabel('time [s]')
% 
%     subplot(334)
%     plot(t, u_gps)
%     grid on
%     hold on
%     plot(t, u_n)
%     hold off
%     title(strcat("GPS: U, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Velocity, wind corrected');
%     ylabel('Velocity-x [m/s]')
%     xlabel('time [s]')
%     
%     subplot(335)
%     plot(t, v_gps)
%     grid on
%     hold on
%     plot(t, v_n)
%     hold off
%     title(strcat("GPS: V, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Velocity, wind corrected');
%     ylabel('Velocity-y [m/s]')
%     xlabel('time [s]')
% 
%     subplot(336)
%     plot(t, w_gps)
%     grid on
%     hold on
%     plot(t, w_n)
%     hold off
%     title(strcat("GPS: W, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Velocity, wind corrected');
%     ylabel('Velocity-z [m/s]')
%     xlabel('time [s]')
%     
%     subplot(337)
%     plot(t, phi_gps)
%     grid on
%     hold on
%     plot(t, phi)
%     hold off
%     title(strcat("GPS: Roll Angle, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Given Angles');
%     ylabel('Roll Angle [rad]')
%     xlabel('time [s]')
%     
%     subplot(338)
%     plot(t, theta_gps)
%     grid on
%     hold on
%     plot(t, theta)
%     hold off
%     title(strcat("GPS: Pitch Angle, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Given Angles');
%     ylabel('Pitch Angle [rad]')
%     xlabel('time [s]')
%     
%     subplot(339)
%     plot(t, psi_gps)
%     grid on
%     hold on
%     plot(t, psi)
%     hold off
%     title(strcat("GPS: Yaw Angle, Manoeuvre: ",Files(i)));
%     legend('GPS Generated Data', 'Given Angles');
%     ylabel('Yaw Angle [rad]')
%     xlabel('time [s]')
%     
    %% IMU and AIRDATA
%     figure;
%     subplot(331)
%     plot(t, Ax_imu)
%     grid on
%     hold on
%     plot(t, Ax)
%     hold off
%     title(strcat("IMU: x-axis, Manoeuvre: ",Files(i)));
%     legend('IMU Generated Data', 'Acceleration, body frame');
%     ylabel('x-axis [m/s2]')
%     xlabel('time [s]')
%     
%     subplot(332)
%     plot(t, Ay_imu)
%     grid on
%     hold on
%     plot(t, Ay)
%     hold off
%     title(strcat("IMU: y-axis, Manoeuvre: ",Files(i)));
%     legend('IMU Generated Data', 'Acceleration, body frame');
%     ylabel('y-axis [m/s2]')
%     xlabel('time [s]')
%     
%     subplot(333)
%     plot(t, Az_imu)
%     grid on
%     hold on
%     plot(t, Az)
%     hold off
%     title(strcat("IMU: z-axis, Manoeuvre: ",Files(i)));
%     legend('IMU Generated Data', 'Acceleration, body frame');
%     ylabel('z-axis [m/s2]')
%     xlabel('time [s]')
% 
%     subplot(334)
%     plot(t, p_imu)
%     grid on
%     hold on
%     plot(t, p)
%     hold off
%     title(strcat("Roll rate, Manoeuvre: ",Files(i)));
%     legend('IMU Generated Data', 'Given Angular rate');
%     ylabel('Roll Rate [rad/s]')
%     xlabel('time [s]')
%     
%     subplot(335)
%     plot(t, q_imu)
%     grid on
%     hold on
%     plot(t, q)
%     hold off
%     title(strcat("Pitch rate, Manoeuvre: ",Files(i)));
%     legend('IMU Generated Data', 'Given Angular rate');
%     ylabel('Pitch Rate [rad/s]')
%     xlabel('time [s]')
% 
%     subplot(336)
%     plot(t, r_imu)
%     grid on
%     hold on
%     plot(t, r)
%     hold off
%     title(strcat("Yaw rate, Manoeuvre: ",Files(i)));
%     legend('IMU Generated Data', 'Given Angular rate');
%     ylabel('Yaw Rate [rad/s]')
%     xlabel('time [s]')
%     
%     subplot(337)
%     plot(t, V_ardta)
%     grid on
%     title(strcat("Air Data: Air Velocity, Manoeuvre: ",Files(i)));
%     legend('Pitot Data');
%     ylabel('Velocity [m/s]')
%     xlabel('time [s]')
%     
%     subplot(338)
%     plot(t, alpha_ardta)
%     grid on
%     title(strcat("Air Data: Angle of Attack, Manoeuvre: ",Files(i)));
%     legend('Angle of Attack Sensor');
%     ylabel('Angle of attack [rad]')
%     xlabel('time [s]')
%     
%     subplot(339)
%     plot(t, beta_ardta)
%     grid on
%     title(strcat("Air Data: Side slip Angle, Manoeuvre: ",Files(i)));
%     ylabel('Side slip Angle [rad]')
%     legend('Side slip Sensor')
%     xlabel('time [s]')
%     
   end
%% Generate Measurement and Input
U_k = [Ax_imu, Ay_imu, Az_imu, p_imu, q_imu, r_imu];
Z_k = [u_gps, v_gps, w_gps, x_gps, y_gps, z_gps, phi_gps, theta_gps, psi_gps, V_ardta, alpha_ardta, beta_ardta];
raw = [u, v, w, x, y, z, phi, theta, psi, vtas, alpha, beta];
