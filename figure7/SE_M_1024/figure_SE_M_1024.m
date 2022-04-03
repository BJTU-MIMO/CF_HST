
hold on; box on;

plot(M,LSFD100,'r-o','LineWidth',2);
plot(M,LSFD300,'b--s','LineWidth',2);
plot(M,LSFD600,'k-.>','LineWidth',2);

xlabel('Number of subcarriers ($M$)','Interpreter','latex');
ylabel('Average SE at initial position (bit/s/Hz)','Interpreter','latex');
legend('$v = 100$ km/h','$v = 300$ km/h','$v = 600$ km/h','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
axis([0 64 6.5 7.2]);

