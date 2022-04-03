hold on; box on;

data = [8.9938; 8.6810; 8.3631; 7.8827; 6.6549; 3.9520];
%data = [8.9938 8.6810 8.3631 7.8827 6.6549 3.9520];
x = 1:6;

y = bar(data,0.3);
plot(x,data,'r-o','LineWidth',2);

set(gca,'xtick',(1:6));
set(gca,'XTickLabel',{'160\times1';'80\times2';'40\times4';'20\times8';'10\times16';'5\times32'});

ylabel('Average SE (bit/s/Hz)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

legend('Parameters: $L \times N$','Values','Interpreter','latex');

grid on;
axis([0.5 6.5 0 10]) 