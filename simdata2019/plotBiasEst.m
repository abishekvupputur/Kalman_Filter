%% TASK 3.2
useful_data = 1;
t = 0:dt:120-dt;
ref = 0.001*ones(size(t));
figure;
    subplot(2,3,1)
    plot(t(useful_data:end),XX_k1k1(10,useful_data:end))
    hold on
    plot(t(useful_data:end), ref(useful_data:end))
    legend('Estimated Value','Groundtruth')
    grid on
    title('Acceleration Bias X')
    xlabel('time [s]')
    ylabel('\lambda_X [m/s2]')
    ylim([-0.03 0.03])
    
    subplot(2,3,2)
    plot(t(useful_data:end),XX_k1k1(11,useful_data:end))
    hold on
    plot(t(useful_data:end), ref(useful_data:end))
    legend('Estimated Value','Groundtruth')
    grid on
    title('Acceleration Bias Y')
    xlabel('time [s]')
    ylabel('\lambda_Y [m/s2]')
    ylim([-0.03 0.03])
    
    subplot(2,3,3)
    plot(t(useful_data:end),XX_k1k1(12,useful_data:end))
    hold on
    plot(t(useful_data:end), ref(useful_data:end))
    legend('Estimated Value','Groundtruth')
    grid on
    title('Acceleration Bias Z')
    xlabel('time [s]')
    ylabel('lambda_Z [m/s2]')
    ylim([-0.03 0.03])
    
    subplot(2,3,4)
    plot(t(useful_data:end),rad2deg(XX_k1k1(13,useful_data:end)))
    hold on
    plot(t(useful_data:end), ref(useful_data:end))
    legend('Estimated Value','Groundtruth')
    grid on
    title('Gyro Drift Roll rate')
    xlabel('time [s]')
    ylabel('\lambda_p [deg/s]')
    ylim([-0.03 0.03])
    
    subplot(2,3,5)
    plot(t(useful_data:end),rad2deg(XX_k1k1(14,useful_data:end)))
    hold on
    plot(t(useful_data:end), ref(useful_data:end))
    legend('Estimated Value','Groundtruth')
    grid on
    title('Gyro Drift Pitch rate')
    xlabel('time [s]')
    ylabel('\lambda_q [deg/s]')
    ylim([-0.03 0.03])
    
    subplot(2,3,6)
    plot(t(useful_data:end),rad2deg(XX_k1k1(15,useful_data:end)))
    hold on
    plot(t(useful_data:end), ref(useful_data:end))
    legend('Estimated Value','Groundtruth')
    grid on
    title('Gyro Drift Yaw rate')
    xlabel('time [s]')
    ylabel('\lambda_r [deg/s]')
    ylim([-0.03 0.03])
set(gcf,'units','points','position',[10,10,1000,800])
File_names = ["da3211" ; "de3211" ; "dadoublet" ; "drdoublet"; "dr3211"];
saveas(gcf,strcat('biasEst_',File_names(dataset),'.png'));