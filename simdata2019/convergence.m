for i=2:N+1
   Z_k1k1(1:12,i)=kf_calc_h(0, XX_k1k1(:,i-1), zeros(6,1));
end
converged_value = 3000;
Residuals = Z_k(:,converged_value:end) - Z_k1k1(:,converged_value:end);
Noise = Z_k(:,converged_value:end) - raw(converged_value:end,:)';

%%
figure;
for i=1:12
    subplot(3,4, i)
    pd(i) = fitdist(Residuals(i,:)','Normal');
    [f,xi] = ksdensity(Residuals(i,:));
    plot(xi,f,'LineWidth',2);
    hold on
    [fx,xxi] = ksdensity(Noise(i,:));
    plot(xxi,fx,'LineWidth',1);
    switch i
        case 1
            title('GPS U');
            xlabel('Residual_{Velocity} [m/s]');
        case 2
            title('GPS V');
            xlabel('Residual_{Velocity} [m/s]');
        case 3
            title('GPS W');
            xlabel('Residual_{Velocity} [m/s]');
        case 4
            title('GPS X');
            xlabel('Residual_{Position} [m]');
        case 5
            title('GPS Y');
            xlabel('Residual_{Position} [m]');
        case 6
            title('GPS Z');
            xlabel('Residual_{Position} [m]');
        case 7
            title('GPS \Phi');
            xlabel('Residual_{Angle} [rad]');
        case 8
            title('GPS \Theta');
            xlabel('Residual_{Angle} [rad]');
        case 9
            title('GPS \Psi');
            xlabel('Residual_{Angle} [rad]');
        case 10
            title('Airdata V_{tas}');
            xlabel('Residual_{Velocity} [m/s]');
        case 11
            title('Airdata \alpha');
            xlabel('Residual_{Angle} [rad]');
        case 12
            title('Airdata \beta');
            xlabel('Residual_{Angle} [rad]');
    end
    grid on
    hold on
    xline(pd(i).mu,'-',{'Mean:',num2str(pd(i).mu)},'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    hold on
    xline(pd(i).mu + pd(i).sigma,'-',{'STD:',num2str(pd(i).sigma + pd(i).mu)},'LabelHorizontalAlignment','right','LabelVerticalAlignment','top');
    hold on
    xline(pd(i).mu-1*pd(i).sigma,'-',{'STD:',num2str(pd(i).mu-1*pd(i).sigma)},'LabelHorizontalAlignment','left','LabelVerticalAlignment','top');
    legend('Noise Estimated', 'Noise Added','location','southeast','Orientation','horizontal')
    [h(i), pp(i)] = ztest(Residuals(i,:), mean(Noise(i,:)), std(Noise(i,:)));
    htt = normalitytest(Residuals(i,:));
    hti(i,:) = htt(:,3);
    pti(i,:) = htt(:,2);
end
set(gcf,'units','points','position',[10,10,1600,800])
File_names = ["da3211" ; "de3211" ; "dadoublet" ; "drdoublet"; "dr3211"];
saveas(gcf,strcat('Convergence_',File_names(dataset),'.png'));

St = ["GPS U";"GPS V";"GPS W";"GPS X";"GPS Y";"GPS Z";"GPS \Phi";"GPS \Theta";"GPS \Psi";"Airdata Vtas";"Airdata \alpha";"Airdata \beta"];
T = table(St, (mean(hti(:,1:3),2) > 0.5),pti(:,1), pti(:,2), pti(:,3), ~h', pp')
writetable(T,strcat('Convergence_',File_names(dataset),'.txt'),'Delimiter',';');
