 clc
 clear
 close all
tic

nmodel=200;  %%������
nter=50;     %%����������
%%����Ⱥ�Ż�����

[obstime,obshz,calhz,njter,ymin,yaver,iter,gbest] = INVpsoN(nmodel,nter);

obscal=[obstime',obshz',calhz'];  %%%���ʵ��ֵ
save TEM.txt obscal -ascii ;  %%���ļ����


% save('F:\\����������\\�������ŷ�˶ʿ����\\��ҵ����-�ܽ�\\3���Դ����\\һά����\\PSO\\da.txt','shuchu');
wucha=[njter',ymin'];
save ������.txt wucha -ascii ;
gbest=gbest';
save gbest.txt gbest -ascii;
%%�����ͼ������½�����ͼ
toc;%��ʱ����
figure(1);   
loglog(obstime,obshz,'-k',obstime,calhz,'*k');%˫��������ʾ
% legend(strcat('��ʵģ�ͦ�=',num2str(realres),'��.m'),strcat('��ʼģ�ͦ�=',num2str(initialres),'��.m'),strcat('����ģ�ͦ�=',num2str(calres),'��.m'));
% legend(['����ģ�ͦ�=',num2str(realmodel)],['����ģ�ͦ�=',num2str(gbest)]);
legend('TEMģ������ֵ','TEMģ��QPSO����ֵ');
xlabel('ʱ��(s)','FontSize',13);
ylabel('hz (A/m)','FontSize',13);
grid on
hold on
% title('˲����Ӧ','FontSize',18);
loglog(obstime,obshz,'-k',obstime,calhz1,'ok',obstime,calhz,'*k');
legend('TEMģ������ֵ','TEMģ��PSO����ֵ','TEMģ��QPSO����ֵ');


% figure(2);
% semilogy(njter,ymin,'b:o');%��������ֵ����������仯ͼ��
% xlabel('���Ϸ��ݵ�������','Fontsize',13);
% ylabel('���Ϸ���������','Fontsize',13);
% 

gbest
iter
ymin(end)
% dzl'