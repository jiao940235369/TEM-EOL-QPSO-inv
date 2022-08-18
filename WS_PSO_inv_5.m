 clc
 clear
 close all
tic

nmodel=200;  %%粒子数
nter=50;     %%最大迭代次数
%%粒子群优化反演

[obstime,obshz,calhz,njter,ymin,yaver,iter,gbest] = INVpsoN(nmodel,nter);

obscal=[obstime',obshz',calhz'];  %%%输出实测值
save TEM.txt obscal -ascii ;  %%改文件编号


% save('F:\\基金与论文\\课题组电磁法硕士论文\\毕业交接-周金\\3软件源程序\\一维反演\\PSO\\da.txt','shuchu');
wucha=[njter',ymin'];
save 拟合误差.txt wucha -ascii ;
gbest=gbest';
save gbest.txt gbest -ascii;
%%画拟合图和误差下降曲线图
toc;%计时结束
figure(1);   
loglog(obstime,obshz,'-k',obstime,calhz,'*k');%双对数域显示
% legend(strcat('真实模型ρ=',num2str(realres),'Ω.m'),strcat('初始模型ρ=',num2str(initialres),'Ω.m'),strcat('反演模型ρ=',num2str(calres),'Ω.m'));
% legend(['理论模型ρ=',num2str(realmodel)],['反演模型ρ=',num2str(gbest)]);
legend('TEM模型理论值','TEM模型QPSO反演值');
xlabel('时间(s)','FontSize',13);
ylabel('hz (A/m)','FontSize',13);
grid on
hold on
% title('瞬变响应','FontSize',18);
loglog(obstime,obshz,'-k',obstime,calhz1,'ok',obstime,calhz,'*k');
legend('TEM模型理论值','TEM模型PSO反演值','TEM模型QPSO反演值');


% figure(2);
% semilogy(njter,ymin,'b:o');%绘制最优值随迭代次数变化图像
% xlabel('联合反演迭代次数','Fontsize',13);
% ylabel('联合反演拟合误差','Fontsize',13);
% 

gbest
iter
ymin(end)
% dzl'