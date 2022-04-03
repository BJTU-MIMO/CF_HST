
LSFD_Scalable_th10 = squeeze(min(mean(SE_n_LSFD_Scalable_th10(:,:,:,:),2)));
LSFD_Scalable_Power_th10 = squeeze(min(mean(SE_n_LSFD_Scalable_Power_th10(:,:,:,:),2)));
MF_Scalable_MaxMin_LSFD_th10 = squeeze(min(mean(SE_n_MF_Scalable_MaxMin_LSFD_th10(:,:,:,:),2)));
LSFD_Scalable_th30 = squeeze(min(mean(SE_n_LSFD_Scalable_th30(:,:,:,:),2)));
LSFD_Scalable_Power_th30 = squeeze(min(mean(SE_n_LSFD_Scalable_Power_th30(:,:,:,:),2)));
MF_Scalable_MaxMin_LSFD_th30 = squeeze(min(mean(SE_n_MF_Scalable_MaxMin_LSFD_th30(:,:,:,:),2)));

s = 0;
Number_new = Number-2*s;
ss = (s+1):(Number-s);

%%

hold on; box on;

plot(speed,squeeze(mean(MF_Scalable_MaxMin_LSFD_th30(ss,:),1)),'b:s','LineWidth',2);
plot(speed,squeeze(mean(LSFD_Scalable_Power_th30(ss,:),1)),'b:>','LineWidth',2);
plot(speed,squeeze(mean(LSFD_Scalable_th30(ss,:),1)),'b:o','LineWidth',2);
plot(speed,squeeze(mean(MF_Scalable_MaxMin_LSFD_th10(ss,:),1)),'r-s','LineWidth',2);
plot(speed,squeeze(mean(LSFD_Scalable_Power_th10(ss,:),1)),'r->','LineWidth',2);
plot(speed,squeeze(mean(LSFD_Scalable_th10(ss,:),1)),'r-o','LineWidth',2);

xlabel('Speed (km/h)','Interpreter','latex');
ylabel('Average SE of the worst TA (bit/s/Hz)','Interpreter','latex');
legend('$\Theta = 30$ dB (max-min)','$\Theta = 30$ dB (fractional)','$\Theta = 30$ dB (full)','$\Theta = 10$ dB (max-min)','$\Theta = 10$ dB (fractional)','$\Theta = 10$ dB (full)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
axis([100 600 1.6 6.5]);





