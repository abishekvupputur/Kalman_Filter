%figure;
%heatmap(P_k_1k_1)
for i=2:N+1
   Z_k1k1(1:12,i)=kf_calc_h(0, XX_k1k1(:,i-1), zeros(6,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
%title(Files(dataset));
for i=1:12
    i
    subplot(3,4, i)
    plot(t,Z_k(i,2:end));
    hold on
    plot(t,z_pred(i,1:end));
    grid on
    switch i
        case 1
            title('GPS U');
            ylabel('Velocity [m/s]');
        case 2
            title('GPS V');
            ylabel('Velocity [m/s]');
        case 3
            title('GPS W');
            ylabel('Velocity [m/s]');
        case 4
            title('GPS X');
            ylabel('Position [m]');
        case 5
            title('GPS Y');
            ylabel('Position [m]');
        case 6
            title('GPS Z');
            ylabel('Position [m]');
        case 7
            title('GPS \Phi');
            ylabel('Angle [rad]');
        case 8
            title('GPS \Theta');
            ylabel('Angle [rad]');
        case 9
            title('GPS \Psi');
            ylabel('Angle [rad]');
        case 10
            title('Airdata V_{tas}');
            ylabel('Velocity [m/s]');
        case 11
            title('Airdata \alpha');
            ylabel('Angle [rad]');
        case 12
            title('Airdata \beta');
            ylabel('Angle [rad]');
    end
    xlabel('Time [s]')
    legend('Measured value','Estimated Value');
end    
%suptitle(char(Files(dataset)));
set(gcf,'units','points','position',[10,10,1600,800])
File_names = ["da3211" ; "de3211" ; "dadoublet" ; "drdoublet"; "dr3211"];
saveas(gcf,strcat('IEKF_',File_names(dataset),'.png'));

% figure;
% for i=1:12
%     subplot(3,4, i)
%     plot(STDx_cor(i,:));
%     grid on
% end
%%
%LSE
%%Plot Bias Estimates
%plotBiasEst
%% Plot Wind States
x_wind = 10*ones(size(t));
y_wind = 6*ones(size(t));
z_wind = 1*ones(size(t));

figure;
subplot(3,1,1)
plot(t,XX_k1k1(16,:))
hold on
plot(t,x_wind)
grid on
title('Windspeed Estimate in X-axis');
ylabel('Velocity [m/s]');
xlabel('Time [s]');
legend('Estimated value','True Value');

subplot(3,1,2)
plot(t,XX_k1k1(17,:))
hold on
plot(t,y_wind)
grid on
title('Windspeed Estimate in Y-axis');
ylabel('Velocity [m/s]');
xlabel('Time [s]');
legend('Estimated value','True Value');

subplot(3,1,3)
plot(t,XX_k1k1(18,:))
hold on
plot(t,z_wind)
grid on
title('Windspeed Estimate in Z-axis');
ylabel('Velocity [m/s]');
xlabel('Time [s]');
legend('Estimated value','True Value');

set(gcf,'units','points','position',[10,10,1000,800])
File_names = ["da3211" ; "de3211" ; "dadoublet" ; "drdoublet"; "dr3211"];
saveas(gcf,strcat('Wind_',File_names(dataset),'.png'));