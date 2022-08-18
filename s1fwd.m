function [depth,rhos] = s1fwd(ninv,resh)

%
% c系数，来自程志平1997，p
%

c=[0.003042,-0.001198,0.01284,0.0235,0.08688,0.2374,0.6194,1.1817,0.4248,-3.4507...
    2.7044,-1.1324,0.393,-0.1436,0.05812,-0.02521,0.01125,-0.004978,0.002072,-0.000318];
ab2=10:10:300;
% ab2 = [
% 10
% 20
% 30
% 40
% 50
% 60
% 70
% 80
% 90
% 100
% 110
% 120
% 130
% 140
% 150
% 160
% 170
% 180
% 190
% 200
%    ];
nab2 = length(ab2);
depth=ab2;
N=(ninv+1)/2;
res=resh(1:N);
thick=resh(N+1:2*N-1);
nlayer=length(res);

sum = 0;

for iab2 = 1:nab2
    r = ab2(iab2);

for k=1:20
    %
    % 求电阻率转换函数t(k)
    %
    T = res(nlayer);
    m= exp(k*log(10)/6-2.1719)/r;
    for i = nlayer-1:-1:1
        tt1 = 1 - exp(-2*m*thick(i));
        tt2 = 1 + exp(-2*m*thick(i));
        T   = res(i)*(res(i)*tt1+T*tt2)/(res(i)*tt2+T*tt1);
    end
    t(k) = T;
    sum = sum + t(k)*c(k);
end
    rhos(iab2) =sum;
    sum=0;
end
return