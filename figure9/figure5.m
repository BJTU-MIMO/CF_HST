

LSFD_cor01 = squeeze(mean(mean(SE_n_LSFD_cor01(:,:,:,:),2),1));
LSFD_cor1000 = squeeze(mean(mean(SE_n_LSFD_cor1000(:,:,:,:),2),1));
LSFD_un01 = squeeze(mean(mean(SE_n_LSFD_un01(:,:,:,:),2),1));
LSFD_un1000 = squeeze(mean(mean(SE_n_LSFD_un1000(:,:,:,:),2),1));

%%

hold on; box on;

s = 0;
%Number_new = Number-2*s;
ss = (s+1):(Number-s);
plot(speed,squeeze(mean(LSFD_un1000(ss,:),1)),'b--o','LineWidth',2);
plot(speed,squeeze(mean(LSFD_cor1000(ss,:),1)),'m:+','LineWidth',2);
plot(speed,squeeze(mean(LSFD_un01(ss,:),1)),'r-s','LineWidth',2);
plot(speed,squeeze(mean(LSFD_cor01(ss,:),1)),'k-.>','LineWidth',2);

xlabel('Speed (km/h)','Interpreter','latex');
ylabel('Average SE (bit/s/Hz)','Interpreter','latex');
legend('$\bar{K} = 30$ dB (uncorrelated)','$\bar{K} = 30$ dB (ASD = $10^\mathbf{o}$)','$\bar{K} = -10$ dB (uncorrelated)','$\bar{K} = -10$ dB (ASD = $10^\mathbf{o}$)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
%axis([100 600 0 10]);

%%
magnify;
