  function  [tt,bbz]=gdduocengzhengyan(toff,pou)
 % �ó���Ϊһά��˲��������ģ�����
 p=[50 200 100];   % �ײ������
 h=[100 50];        % �ز���
 resh=[p,h];
 ninv=length(resh);
%  a=100/sqrt(pi);             % ����뾶

 I0=1;
%  toff=0.0003;  %%�Ĺض�ʱ��
 bbz=zeros(41,1);
%  kk=zeros(41,1);
 tt=logspace(-6,-1,41);   % ����ʱ��
%  uu=tt;
 tic       
 nter=length(tt);
 hwait=waitbar(0,'��ʼ����');
 step=nter/100;
 for i=1:nter
     Str=fix(i/step);sstr=['����������',num2str(Str),'%'];waitbar(i/nter,hwait,sstr); pause(0.05);
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
  save g300.txt  kkk  -ascii;    % ��������� ��һ����ʱ�� �ڶ����Ǵų��ĵ���
 