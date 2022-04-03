
LSFD_L20 = squeeze(mean(mean(SE_n_LSFD_L20(:,:,:,:),2),1));
LSFD_L40 = squeeze(mean(mean(SE_n_LSFD_L40(:,:,:,:),2),1));

%%

s = 0;
%Number_new = Number-2*s;
ss = (s+1):(Number-s);

d_L20 = squeeze(mean(LSFD_L20(ss,:),1));
d_L40 = squeeze(mean(LSFD_L40(ss,:),1));

%%
hold on; box on;

plot(d_ve,d_L40,'b--s','LineWidth',2);
plot(d_ve,d_L20,'r-o','LineWidth',2);

xlabel('$d_\mathrm{ve}$ (m)','Interpreter','latex');
ylabel('Average SE (bit/s/Hz)','Interpreter','latex');
legend('$L=40$','$L=20$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
axis([10 250 5.5 9]);



