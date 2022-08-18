% 中心回线中心点感应电动势解析解
% H=0是地面瞬变电磁法 均匀大地解析解
clc;clear;
% a=100/(pi^0.5);
a=100/sqrt(pi);
I0=1;
u0=4*pi*10^-7;
cc=load('TEMtoff300.txt');
t=cc(:,1);
shice=cc(:,2);
m=length(t);
toff=0.0003;
sr=1;
% bbz=zeros(m,1);
quanqi=zeros(m,1);
W=zeros(m,1); %数据的位数
 sum1=0.0;
 tic
 nter=length(t);
 hwait=waitbar(0,'开始计算');
 step=nter/100;

for i=1:m
    Str=fix(i/step);sstr=['正在运行中',num2str(Str),'%'];waitbar(i/nter,hwait,sstr); pause(0.05);
    p1=0.1;
    p2=1000;
    cishu=1;
     for j=-6:12  %-4
        sc=shice(i)*10^j;
        if (1<sc && sc<10)
            W(i,1)=j;
        end
     end
    wucha=10^(-W(i,1)-4); %变误差精度
%   wucha=0.001;
    kk=2;
     while(kk>wucha)
         p=(p1+p2)/2; 
         for j=-49:1:0
             k(i)=t(i)-(1*j)*toff/50;
             u=a*(u0/(p*k(i)))^0.5/2;     
             Hz=I0*(3/sqrt(pi)/u*exp(-u^2)+(1-3/2/u^2)*erf(u))/2/a/50;%%        %均匀大地磁场
%              Vt=(I0*p1*sr/a^3)*(3*erf(u)-2*u*(3+2*u^2)*exp(-u^2)/sqrt(pi)); %均匀大地磁场对时间的导数
             sum1=sum1+Hz;
%              sum2=sum2+Vt;
         end
         bbz=sum1;
         
         kk=abs(bbz-shice(i));
             if(bbz<shice(i))
                 p2=p;
             elseif(bbz>shice(i))
                 p1=p;
             else 
                 break;
             end
             
             sum1=0.0;
      end
         quanqi(i)=p;
                              
end
close(hwait);
toc;
figure(3);
loglog(t,quanqi,'.-b','LineWidth',0.1)    % 解析解
xlabel('时间t/ms');
ylabel('电阻率Ω*m');
title('考虑time-off的全期视电阻率');
axis([10.^-5 1 0.1 1000]);
grid on
% set(gca,'xgrid','on');
% set(gca,'ygrid','on');
% % set(gca,'xminorgrid','off');  %去辅刻度        
% dzl1=u0/pi.*t.^(-1)*(70*70*1/30).^(2/3).*(bbz).^(-2/3);
% dzl2=u0/(4*pi).*t.^(-1)*(2*u0*70*70*sr/5).^(2/3).*t.^(-2/3).*vvz.^(-2/3);
%%此处晚期公式采用电动势计算出的电阻率为3.7，而磁场计算出的电阻率为50；前者为视电阻率
aa=[t,quanqi];
save qq1.txt aa -ascii;
 
