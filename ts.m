clc;clear
 % �ó���Ϊһά��˲��������ģ�����

 p=[50 200 100];   % �ײ������
 h=[100 50];        % �ز���
 resh=[p,h];
 ninv=length(resh);
 N=(ninv+1)/2;
 thk=resh(N+1:2*N-1);
 [tt,bbz]=TEMloop(ninv,resh)