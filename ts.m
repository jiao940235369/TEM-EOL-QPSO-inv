clc;clear
 % 该程序为一维的瞬变电磁正演模拟程序

 p=[50 200 100];   % 底层电阻率
 h=[100 50];        % 地层厚度
 resh=[p,h];
 ninv=length(resh);
 N=(ninv+1)/2;
 thk=resh(N+1:2*N-1);
 [tt,bbz]=TEMloop(ninv,resh)