

LSFD_L10 = squeeze(mean(mean(SE_n_LSFD_L10(:,:,:),2),1));
LSFD_L20 = squeeze(mean(mean(SE_n_LSFD_L20(:,:,:),2),1));
LSFD_L40 = squeeze(mean(mean(SE_n_LSFD_L40(:,:,:),2),1));
SmallCell_L10 = squeeze(mean(mean(SE_n_SmallCell_L10(:,:,:),2),1));
SmallCell_L20 = squeeze(mean(mean(SE_n_SmallCell_L20(:,:,:),2),1));
SmallCell_L40 = squeeze(mean(mean(SE_n_SmallCell_L40(:,:,:),2),1));
Cellular_MMSE_L10 = squeeze(mean(mean(SE_n_Cellular_MMSE_L10(:,:,:),2),1));
Cellular_MMSE_L20 = squeeze(mean(mean(SE_n_Cellular_MMSE_L20(:,:,:),2),1));
Cellular_MMSE_L40 = squeeze(mean(mean(SE_n_Cellular_MMSE_L40(:,:,:),2),1));


hold on; box on;

s = 0;
%Number_new = Number-2*s;
ss = (s+1):(Number-s);

LSFD_L10_mean = mean(LSFD_L10(ss));
LSFD_L20_mean = mean(LSFD_L20(ss));
LSFD_L40_mean = mean(LSFD_L40(ss));

SmallCell_L10_mean = mean(SmallCell_L10(ss));
SmallCell_L20_mean = mean(SmallCell_L20(ss));
SmallCell_L40_mean = mean(SmallCell_L40(ss));

Cellular_MMSE_L10_mean = mean(Cellular_MMSE_L10(ss));
Cellular_MMSE_L20_mean = mean(Cellular_MMSE_L20(ss));
Cellular_MMSE_L40_mean = mean(Cellular_MMSE_L40(ss));

%data = [LSFD_L10_mean,LSFD_L20_mean,LSFD_L30_mean,LSFD_L40_mean; SmallCell_L10_mean,SmallCell_L20_mean,SmallCell_L30_mean,SmallCell_L40_mean; Cellular_L10_mean,Cellular_L20_mean,Cellular_L30_mean,Cellular_L40_mean];
data = [LSFD_L10_mean,LSFD_L20_mean,LSFD_L40_mean; SmallCell_L10_mean,SmallCell_L20_mean,SmallCell_L40_mean; Cellular_MMSE_L10_mean,Cellular_MMSE_L20_mean,Cellular_MMSE_L40_mean];

y = bar(data);

%ch = get(y,'children');
%xlabel('Dataset','FontSize',12);
set(gca,'xtick',(1:3));
set(gca,'XTickLabel',{'LSFD';'Small Cell';'Cellular'});

ylabel('Average SE (bit/s/Hz)','Interpreter','latex');

legend('$L = 10$','$L = 20$','$L = 40$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

grid on;
axis([0.5 3.5 0 9]) 
