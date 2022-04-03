

LSFD_or = squeeze(mean(mean(SE_n_LSFD_or(:,:,:,:),2),1));
LSFD_K6 = squeeze(mean(mean(SE_n_LSFD_K6(:,:,:,:),2),1));
LSFD_N6 = squeeze(mean(mean(SE_n_LSFD_N6(:,:,:,:),2),1));

%%

hold on; box on;

s = 0;
%Number_new = Number-2*s;
ss = (s+1):(Number-s);
plot(speed,squeeze(mean(LSFD_K6(ss,:),1)),'m:>','LineWidth',2);
plot(speed,squeeze(mean(LSFD_N6(ss,:),1)),'k-.o','LineWidth',2);

plot(speed,squeeze(mean(LSFD_or(ss,:),1)),'r-s','LineWidth',2);


xlabel('Speed (km/h)','Interpreter','latex');
ylabel('Average SE (bit/s/Hz)','Interpreter','latex');
legend('$K = 6$, $N = 4$','$K = 8$, $N = 6$','$K = 8$, $N = 4$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
axis([100 600 6 8]);




