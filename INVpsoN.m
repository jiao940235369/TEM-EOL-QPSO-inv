
function [obstime,obshz,calhz,njter,ymin,yaver,iter,gbest] = INVpsoN(nmodel,nter)

clear obstime obshz calhz  njter ymin yaver iter gbest;

%%%%%%%%%%%%%%%%%%%以下为TEM初始模型%%%%%%%%%%%%%%%%
% minres=[1 1 1 1];%每层最低电阻率
% mindeep=[5 5 5];%每层最小厚度（除去最上一层、最下一层）
% maxres=[300 300 200 200];%每层最高电阻率
% maxdeep=  [100 100 100];%每层最大厚度（除去最下一层）
% 
% realres=[200 50 300 150];%各层电阻率，从最上一层到最下一层
% realdeep= [50 20 50];%除最上一层和最下一层的厚度
% realmodel=[realres,realdeep];%实际模型
% ninv=length(realmodel);% 反演参数个数


%%%%%%%%%%%%%%%%% 输入模型空间取值范围%%%%%%%%%%%%%%
% minres=repmat(10,1,15);
% mindeep=repmat(20,1,14);
% maxres=repmat(500,1,15);
% maxdeep=repmat(20,1,14);
minres=[10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];%每层最低电阻率
mindeep=  [20 20 20 20 20 20 20 20 20 20 20 20 20 20];%每层最小厚度（除去最上一层、最下一层）
maxres=[1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];%每层最高电阻率
maxdeep=  [20 20 20 20 20 20 20 20 20 20 20 20 20 20];%每层最大厚度（除去最上一层、最下一层）
mmin=[minres,mindeep];%最小模型参数
mmax=[maxres,maxdeep];%最大模型参数
ninv=length(mmin);
realres=[50 200 50];%各层电阻率，从最内一层到最外一层
realdeep= [100 100];%除最外一层的厚度
realmodel=[realres,realdeep];%实际模型
ninv=length(realmodel);% 反演参数个数
% mmin=realmodel*0.5;%最小模型参数
% mmax=realmodel*1.5;%最大模型参数


%%%%%%%%%%%%%%%%%%%%%%%%%%以下为DC初始模型%%%%%%%%%%%%%%%%%%%%%%%
% minres=[50 50 1 1 1 50 50 50 50 100 100 100];%每层最低电阻率
% mindeep=  [5 5 5 5 5 5 5 5 5 5 5 5];%每层最小厚度（除去最上一层、最下一层）
% maxres=[300 300 100 100 100 200 200 200 200 400 400 400];%每层最高电阻率
% maxdeep=  [20 20 10 10 10 20 20 20 20 10 10 10];%每层最大厚度（除去最下一层）

mmin=[minres,mindeep];%最小模型参数
mmax=[maxres,maxdeep];%最大模型参数
ninv=length(mmin);% 反演参数个数
tic;

% 设置速度限制
ResDeepSpan=(mmax - mmin);%[电阻率&层厚]尺度
vmax=ResDeepSpan*0.2;%最大追踪速度  所有参数
vmin = -1 * vmax;%最小追踪速度
lamd_z=0.0000001;  %正则化因子

k=2;
[Lz]=gnerate_Lz_w(ninv);
% 初始化速度
v=vmax*0.5;
c1 = 2; 
c2 = 2;
eps = 1.0e-6;%相对误差因子
beta=1.5;
% 读取实测值
shiceTEM=load('3cengmodelKnoisy.txt');
obstime=shiceTEM(:,1);
obshz=shiceTEM(:,2)*1;%利用正演函数过程计算实际模型的 响应时间及磁场强度,代表实测值
obshz=obshz';
obstime=obstime';
ntime=length(obstime);%时间门数量

% obshz=(floor(rand(1,5)*2)*2-1)*obshz*0.1+obshz;
% for i=1:ntime     %%加噪声
%     obshz(i)=(floor(rand(1,1)*2)*2-1)*0.2*obshz(i)+obshz(i);
% end

% 随机产生nmodel个初始模型   将每个随机模型的误差都存入xxerror中，个体最优误差存入pberror中
for imodel = 1 : nmodel  
%     v(imodel,:) = v(1,:);%初始速度
    
    % 产生随机模型并初始化粒子群，初始化个体极值
  %     [res,deep] = creatmodel(nlayer,minres,maxres,mindeep,maxdeep);%产生初始模型[电阻率，层厚]=creatmodel[层数，最小电阻率，最大，最小层厚，最大]
         model(imodel,:)=creatmodel(ninv,mmin,mmax);%产生初始模型,赋值给一个粒子
         alpha=rand;
         fxmodel(imodel,:)=alpha*(mmin+mmax)-model(imodel,:);%%%在此加反向学习 dy
        for i=1:ninv 
         if fxmodel(imodel,i)<mmin(i)
%              fxmodel(imodel,i)=mmin(i)+(mmax(i)-mmin(i))*alpha;
            fxmodel(imodel,i)=mmin(i)+mmax(i)-model(imodel,i);
         elseif fxmodel(imodel,i)>mmax(i)
             fxmodel(imodel,i)=mmin(i)+model(imodel,i);
         end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%二阶微分
%    meme(imodel,:)=model(imodel,:);
    pbest(imodel,:) = model(imodel,:);  %个体最优
    fxpbest(imodel,:)=fxmodel(imodel,:);
%  对生成的模型进行正演
    [~,calhz]=TEMloop(ninv,model(imodel,:));%正演计算  初始模型理论
    [~,fxcalhz]=TEMloop(ninv,fxmodel(imodel,:));%正演计算  初始模型理论
   
%     计算拟和误差
    xerror=(sum(((obshz-calhz)./obshz).^2))/(ntime)+lamd_z*sqrt((Lz*model(imodel,:)')'*(Lz*model(imodel,:)'));
    fxxerror=(sum(((obshz-fxcalhz)./obshz).^2))/(ntime)+lamd_z*sqrt((Lz*fxmodel(imodel,:)')'*(Lz*fxmodel(imodel,:)'));
    if xerror<=fxxerror
          xxerror(imodel) = xerror;
          pbest(imodel,:)=model(imodel,:);
    else
        xxerror(imodel) = fxxerror;
%         model(imodel,:)=fxmodel(imodel,:);
        pbest(imodel,:)=fxpbest(imodel,:);
    end   
    pberror(imodel) = xxerror(imodel);
    v(imodel,:) = pbest(imodel,:);%初始速度
end  % end of imodel

% 求取全局极值
jter       = 1;
njter(1)   = 1;
yaver(1) = mean(xxerror);%第一次迭代误差均值
[gg,sub]   = min(xxerror);%第一次迭代误差最小值，sub为最小值索引 第n个模型最小
yerr(1)   = gg;%最小值 误差值
% bestmodel  = model;%各个粒子的个体最优，因只计算了一次，故最优值也即当前已有值  20个粒子的
ymin(jter) = pberror(sub);%第jter次迭代误差最小值，sub为最小值粒子序号索引
gbest(1,:) = pbest(sub,:);%第一次迭代最优模型（最优粒子）  全局最优化

hwait=waitbar(0,'开始计算');
 step=nter/100;
% 粒子群优化算法计算
for iter = 2 : nter
        Str=fix(iter/step);sstr=['正在运行中',num2str(Str),'%'];waitbar(iter/nter,hwait,sstr); pause(0.05);
        r1 = rand(1,ninv);   %%随机生成ninv个0-1的数
        r2 = rand(1,ninv);
         fai1=c1*r1;
         fai2=c2*r2;
        
    for imodel = 1 : nmodel  
        ssum=0.0;
        for i=1:imodel
            ssum=ssum+pbest(i,:);
        end
        Mbest=ssum/imodel;
      
        

% v(imodel,:) = (0.99^iter*rand(1)/2+0.1) .* v(imodel,:) + fai1.*(pbest(imodel,:) - model(imodel,:)) + fai2.*(gbest(1,:) - pbest(imodel,:));  
v(imodel,:) =r1.*( pbest(imodel,:) ) + (1-r1).*(gbest(1,:) ) ;%./ ( fai1+fai2 )
%   v(imodel,:) = (1+fai1+fai2-sqrt((fai1+fai2)*(2+fai1+fai2))) .* v(imodel,:) + fai1.*(pbest(imodel,:) - model(imodel,:)) + fai2.*(gbest(1,:) - model(imodel,:));   %%%粒子群论文
%  v(imodel,:) = (1+fai1+fai2-sqrt((fai1+fai2)*(2+fai1+fai2))) .* v(imodel,:) + abs( ( fai1*( pbest(imodel,:) ) - fai2*( gbest(1,:)) )) / ( fai1+fai2 );
v(imodel,:) = ( (v(imodel,:) <= vmin).*vmin  ) + ( (v(imodel,:) > vmin).*v(imodel,:) );%速度超出限制范围的量置为范围值，未超限的仍为原值
        v(imodel,:) = ( (v(imodel,:) >= vmax).*vmax  ) + ( (v(imodel,:) < vmax).*v(imodel,:) ); 
        
       %%%%%%%%%%%%%%%以下为量子粒子群位置更新代码%%%%%%%%%%%%%
        u = rand;
        if u>0.5
            model(imodel,:)=v(imodel,:)+( 0.5+0.5*(nter-iter)/nter ).*(v(imodel,:)-model(imodel,:)) .*log(u)+( 0.5+0.5*(nter-iter)/nter ).*(Mbest-model(imodel,:)) .*randn;%pbest(imodel,:)-
        else
            model(imodel,:)=v(imodel,:)-( 0.5+0.5*(nter-iter)/nter ).*(v(imodel,:)-model(imodel,:)) .*log(u)+( 0.5+0.5*(nter-iter)/nter ).*(Mbest-model(imodel,:)) .*randn;%pbest(imodel,:)-
        end

        
%          model(imodel,:) = model(imodel,:) + v(imodel,:);%更新粒子位置  对所有参数均更新位置，而不是某个参数     在此计算出粒子适应度，决定是否改变粒子速度（位置）董
%        model(imodel,:) = pbest(imodel,:) + v(imodel,:);%更新粒子位置  对所有参数均更新位置，而不是某个参数     在此计算出粒子适应度，决定是否改变粒子速度（位置）董 
        
        min_throwaway = model(imodel,:) <= mmin;  % 小于最小限制的位置 为 1
        min_keep      = model(imodel,:) >  mmin;  % 大于最小限制的位置 为 1
        max_throwaway = model(imodel,:) >= mmax;  % 大于最大限制的位置 为 1
        max_keep      = model(imodel,:) <  mmax;  % 小于最大限制的位置 为 1
        
        model(imodel,:) = ( min_throwaway.*mmin ) + ( min_keep.*model(imodel,:) ); %粒子位置超出限制范围的量置为范围值，未超限的仍为原值
        model(imodel,:) = ( max_throwaway.*mmax ) + ( max_keep.*model(imodel,:) );

%        model(imodel,:)=0.1*meme(imodel,:)+0.9*model(imodel,:);
   
        v(imodel,:) = (v(imodel,:).*min_keep) + (-v(imodel,:).*min_throwaway);%对到达边界的粒子速度变号，使其返回继续搜寻
        v(imodel,:) = (v(imodel,:).*max_keep) + (-v(imodel,:).*max_throwaway);
     
        % 正演计算
%         res=model(imodel,1:nlayer);
%         deep=model(imodel,nlayer+1:(2*nlayer-1));
        
        [~,calhz]=TEMloop(ninv,model(imodel,:));%正演计算

        % 计算拟和误差，更新粒子的个体极值         在这里判断是否有必要在更新变动下一个更新变动上一个位置的极值
        xerror=(sum(((obshz-calhz)./obshz).^2))/(ntime)+lamd_z*sqrt((Lz*model(imodel,:)')'*(Lz*model(imodel,:)'));
        xxerror(imodel) = xerror;
        if  xxerror(imodel) < pberror(imodel)
            pberror(imodel) = xxerror(imodel);
            pbest(imodel,:) = model(imodel,:);
        end
    end %end of imodel  
    yaver(iter) = mean(xxerror);
       
    % 求取全局极值
    [ggg,sub] = min(xxerror);
    yerr(iter) = ggg;
    if ggg < gg
        gg = ggg;
        bestmodel = model;
        jter = jter + 1;
        njter(jter) = iter;
        ymin(jter) = xxerror(sub);%群体最优适应度值
        gbest(1,:) = model(sub,:);%群体最优粒子       
        lamd_z=lamd_z*k;
    else
        lamd_z=lamd_z/k;
    end    
    
    % 运算终止条件
    if ymin(jter) < eps
        break;
    end

end %end of iter

 
% 对最优模型进行正演计算
[~,calhz]=TEMloop(ninv,gbest(1,:));
close(hwait);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%显示各个部分占用时间长短
% profile report
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%