subplot(2,1,1);
hold on; box on;

s = 0;
Number_new = Number-2*s;
ss = (s+1):(Number-s);

t = 1;

plot(sort(reshape(SE_n_LSFD_Scalable(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'r-','LineWidth',2);
plot(sort(reshape(SE_n_LSFD_Scalable_Fractional(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'b--','LineWidth',2);
plot(sort(reshape(SE_n_MF_Scalable_MaxMin_LSFD(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'k-.','LineWidth',2);
plot(sort(reshape(SE_n_MF_Scalable_SumSE_LSFD(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'m:','LineWidth',2);

xlabel('SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend('full','fractional','max-min','sum SE','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
axis([2 6 0 1]);

%------------------------------------------
subplot(2,1,2);
hold on; box on;

LSFD_Scalable = squeeze(sum(mean(mean(SE_n_LSFD_Scalable(:,:,:),2),3),1));
LSFD_Scalable_Fractional = squeeze(sum(mean(mean(SE_n_LSFD_Scalable_Fractional(:,:,:),2),3),1));
MF_Scalable_MaxMin_LSFD = squeeze(sum(mean(mean(SE_n_MF_Scalable_MaxMin_LSFD(:,:,:),2),3),1));
MF_Scalable_SumSE_LSFD = squeeze(sum(mean(mean(SE_n_MF_Scalable_SumSE_LSFD(:,:,:),2),3),1));


data = [LSFD_Scalable; LSFD_Scalable_Fractional; MF_Scalable_MaxMin_LSFD; MF_Scalable_SumSE_LSFD];
%data = [8.9938 8.6810 8.3631 7.8827 6.6549 3.9520];
x = 1:4;

y = bar(data,0.25);
plot(x,data,'r-o','LineWidth',2);

set(gca,'xtick',(1:4));
set(gca,'XTickLabel',{'full';'fractional';'max-min';'sum-SE'});

ylabel('Sum SE (bit/s/Hz)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);

legend('power control','values','Interpreter','latex');

grid on;
axis([0.5 4.5 28 33]) 





