clc;
clear
inversiondc=load('DC.txt');
depth=inversiondc(:,1);
obsdc=inversiondc(:,2);
caldc=inversiondc(:,3);
loglog(depth,obsdc,'*r',depth,caldc,'ob');%双对数域显示

%legend(strcat('真实模型ρ=',num2str(realres),'Ω.m'),strcat('初始模型ρ=',num2str(initialres),'Ω.m'),strcat('反演模型ρ=',num2str(calres),'Ω.m'));
% legend(['理论模型ρ=',num2str(realmodel)],['反演模型ρ=',num2str(gbest)]);
legend('DC模型理论值','DC模型反演值');
xlabel('AB/2(m)','FontSize',13);
ylabel('ρ(Ω*m)','FontSize',13);