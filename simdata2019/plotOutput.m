%%Plot output
t = 0:dt:120;
figure;
    subplot(4,3,1)
    plot(t,Z_k(1,:))
    hold on
    plot(t(2:end),z_pred(1,:))
    grid on
    title('GPS Velocity X')
    ylabel('U_{GPS} [m/s]')
    xlabel('time [s]')
    
    subplot(4,3,2)
    plot(t,Z_k(2,:))
    hold on
    plot(t(2:end),z_pred(2,:))
    grid on
    title('GPS Velocity Y')
    ylabel('V_{GPS} [m/s]')
    xlabel('time [s]')
    
    subplot(4,3,3)
    plot(t,Z_k(3,:))
    hold on
    plot(t(2:end),z_pred(3,:))
    grid on
    title('GPS Velocity Z')
    ylabel('W_{GPS} [m/s]')
    xlabel('time [s]')
    
    subplot(4,3,4)
    plot(t,Z_k(4,:))
    hold on
    plot(t(2:end),z_pred(4,:))
    grid on
    title('GPS Position X')
    ylabel('X_{GPS} [m]')
    xlabel('time [s]')
    
    subplot(4,3,5)
    plot(t,Z_k(5,:))
    hold on
    plot(t(2:end),z_pred(5,:))
    grid on
    title('GPS Position Y')
    ylabel('Y_{GPS} [m]')
    xlabel('time [s]')
    
    subplot(4,3,6)
    plot(t,Z_k(6,:))
    hold on
    plot(t(2:end),z_pred(6,:))
    grid on
    title('GPS Position Z')
    ylabel('Z_{GPS} [m]')
    xlabel('time [s]')
    
    subplot(4,3,7)
    plot(t,Z_k(7,:))
    hold on
    plot(t(2:end),z_pred(7,:))
    grid on
    title('GPS Roll Angle')
    ylabel('Phi_{GPS} [rad]')
    xlabel('time [s]')
    
    subplot(4,3,8)
    plot(t,Z_k(8,:))
    hold on
    plot(t(2:end),z_pred(8,:))
    grid on
    title('GPS Pitch Angle')
    ylabel('Theta_{GPS} [rad]')
    xlabel('time [s]')
    
    subplot(4,3,9)
    plot(t,Z_k(9,:))
    hold on
    plot(t(2:end),z_pred(9,:))
    grid on
    title('GPS Yaw Angle')
    ylabel('Psi_{GPS} [rad]')
    xlabel('time [s]')
    
    subplot(4,3,10)
    plot(t,Z_k(10,:))
    hold on
    plot(t(2:end),z_pred(10,:))
    grid on
    title('True Air Speed')
    ylabel('Vtas [m/s]')
    xlabel('time [s]')
    
    subplot(4,3,11)
    plot(t,Z_k(11,:))
    hold on
    plot(t(2:end),z_pred(11,:))
    grid on
    title('Angle of attack')
    ylabel('AoA [rad]')
    xlabel('time [s]')
    
    subplot(4,3,12)
    plot(t,Z_k(12,:))
    hold on
    plot(t(2:end),z_pred(12,:))
    grid on
    title('Side Slip Angle')
    ylabel('Beta [rad]')
    xlabel('time [s]')