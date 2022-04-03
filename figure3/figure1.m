
Centralized_MMSE = squeeze(mean(mean(SE_n_Centralized_MMSE(:,:,:),2),1));
SmallCell_MMSE = squeeze(mean(mean(SE_n_SmallCell_MMSE(:,:,:),2),1));
LSFD_MMSE = squeeze(mean(mean(SE_n_LSFD_MMSE(:,:,:),2),1));

%%
subplot(2,1,1);
hold on; box on;

plot(mobile(1:Number),SmallCell_MMSE(1:Number),'r-','LineWidth',2);
plot(mobile(1:Number),LSFD_MMSE(1:Number),'b-','LineWidth',2);
plot(mobile(1:Number),Centralized_MMSE(1:Number),'m-','LineWidth',2);

xlabel('Position (m)','Interpreter','latex');
ylabel('SE (bit/s/Hz)','Interpreter','latex');
%legend('MF','LSFD','Centralized','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

grid on;
axis([-300 300 0 10]);

%%
%magnify;

%% CDF

s = 0;
Number_new = Number-2*s;
ss = (s+1):(Number-s);

x = [7.714 7.968 8.177 8.351 8.621];
y = [0.1 0.3 0.5 0.7 0.9];

subplot(2,1,2);
hold on; box on;

plot(sort(SmallCell_MMSE(ss)),linspace(0,1,Number_new),'r-','LineWidth',2);
plot(sort(LSFD_MMSE(ss)),linspace(0,1,Number_new),'b--','LineWidth',2);
plot(sort(Centralized_MMSE(ss)),linspace(0,1,Number_new),'m-.','LineWidth',2);
plot(x,y,'ks','LineWidth',1);

xlabel('SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend('Small Cell (MMSE)','LSFD (MMSE)','Centralized (MMSE)','Sequential','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

grid on;
axis([0 10 0 1]);



