  function  [tt,bbz]=gdduocengzhengyan(toff,pou)
 % 该程序为一维的瞬变电磁正演模拟程序
 p=[50 200 100];   % 底层电阻率
 h=[100 50];        % 地层厚度
 resh=[p,h];
 ninv=length(resh);
%  a=100/sqrt(pi);             % 发射半径

 I0=1;
%  toff=0.0003;  %%改关断时间
 bbz=zeros(41,1);
%  kk=zeros(41,1);
 tt=logspace(-6,-1,41);   % 接收时间
%  uu=tt;
 tic       
 nter=length(tt);
 hwait=waitbar(0,'开始计算');
 step=nter/100;
 for i=1:nter
     Str=fix(i/step);sstr=['正在运行中',num2str(Str),'%'];waitbar(i/nter,hwait,sstr); pause(0.05);
      sum=0.0;
         for j=1-pou:1:0
             kk=tt(i)-(1*j)*toff/pou;
%              tt(i)=kk(i); 
               bz=TEMforward(ninv,resh,kk)*I0/pou;
             sum=sum+bz;
%              tt(i)=uu(i);
         end
     bbz(i)=sum;
 end
 close(hwait);
 toc;
 
%  A=tt';
%  B=bz';
%  for i=1:41;
%  dzl(i)=(I*u0^(5/2)*Sr*pi*a*a/(20*sqrt(pi)*B(i)*A(i)^(5/2)))^(2/3);
%  end
% C=[A B dzl'];
% bbz=bz';
  kkk=[tt',bbz];
  figure(1)
  loglog(tt',bbz)
  save g300.txt  kkk  -ascii;    % 保存计算结果 第一列是时间 第二列是磁场的导数
 