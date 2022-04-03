

LSFD = squeeze(mean(mean(SE_n_LSFD(:,:,:),2),1));
SmallCell = squeeze(mean(mean(SE_n_SmallCell(:,:,:),2),1));
Cellular_MMSE = squeeze(mean(mean(SE_n_Cellular_MMSE(:,:,:),2),1));

%%
subplot(2,1,1);
hold on; box on;

plot(mobile(1:Number),SmallCell(1:Number),'r-','LineWidth',2);
plot(mobile(1:Number),LSFD(1:Number),'b-','LineWidth',2);
plot(mobile(1:Number),Cellular_MMSE(1:Number),'m-','LineWidth',2);

xlabel('Position (m)','Interpreter','latex');
ylabel('SE (bit/s/Hz)','Interpreter','latex');
%legend('MF','LSFD','Small Cell','Cellular','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

grid on;
%axis([-300 300 0 4.5]);

%%
%magnify;

%%
s = 0;
Number_new = Number-2*s;
ss = (s+1):(Number-s);

subplot(2,1,2);
hold on; box on;

x = [3.897 4.054 4.149 4.237 4.348];
y = [0.1 0.3 0.5 0.7 0.9];

plot(sort(SmallCell(ss)),linspace(0,1,Number_new),'r-.','LineWidth',2);
plot(sort(LSFD(ss)),linspace(0,1,Number_new),'b--','LineWidth',2);
plot(sort(Cellular_MMSE(ss)),linspace(0,1,Number_new),'m:','LineWidth',2);
plot(x,y,'ko','LineWidth',1);

xlabel('SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend('Small Cell','LSFD','Cellular','Simulation','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

grid on;
%axis([0.9 1.6 0 1]);