clc;
clear
inversiondc=load('DC.txt');
depth=inversiondc(:,1);
obsdc=inversiondc(:,2);
caldc=inversiondc(:,3);
loglog(depth,obsdc,'*r',depth,caldc,'ob');%˫��������ʾ

%legend(strcat('��ʵģ�ͦ�=',num2str(realres),'��.m'),strcat('��ʼģ�ͦ�=',num2str(initialres),'��.m'),strcat('����ģ�ͦ�=',num2str(calres),'��.m'));
% legend(['����ģ�ͦ�=',num2str(realmodel)],['����ģ�ͦ�=',num2str(gbest)]);
legend('DCģ������ֵ','DCģ�ͷ���ֵ');
xlabel('AB/2(m)','FontSize',13);
ylabel('��(��*m)','FontSize',13);