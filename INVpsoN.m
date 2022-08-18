
function [obstime,obshz,calhz,njter,ymin,yaver,iter,gbest] = INVpsoN(nmodel,nter)

clear obstime obshz calhz  njter ymin yaver iter gbest;

%%%%%%%%%%%%%%%%%%%����ΪTEM��ʼģ��%%%%%%%%%%%%%%%%
% minres=[1 1 1 1];%ÿ����͵�����
% mindeep=[5 5 5];%ÿ����С��ȣ���ȥ����һ�㡢����һ�㣩
% maxres=[300 300 200 200];%ÿ����ߵ�����
% maxdeep=  [100 100 100];%ÿ������ȣ���ȥ����һ�㣩
% 
% realres=[200 50 300 150];%��������ʣ�������һ�㵽����һ��
% realdeep= [50 20 50];%������һ�������һ��ĺ��
% realmodel=[realres,realdeep];%ʵ��ģ��
% ninv=length(realmodel);% ���ݲ�������


%%%%%%%%%%%%%%%%% ����ģ�Ϳռ�ȡֵ��Χ%%%%%%%%%%%%%%
% minres=repmat(10,1,15);
% mindeep=repmat(20,1,14);
% maxres=repmat(500,1,15);
% maxdeep=repmat(20,1,14);
minres=[10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];%ÿ����͵�����
mindeep=  [20 20 20 20 20 20 20 20 20 20 20 20 20 20];%ÿ����С��ȣ���ȥ����һ�㡢����һ�㣩
maxres=[1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];%ÿ����ߵ�����
maxdeep=  [20 20 20 20 20 20 20 20 20 20 20 20 20 20];%ÿ������ȣ���ȥ����һ�㡢����һ�㣩
mmin=[minres,mindeep];%��Сģ�Ͳ���
mmax=[maxres,maxdeep];%���ģ�Ͳ���
ninv=length(mmin);
realres=[50 200 50];%��������ʣ�������һ�㵽����һ��
realdeep= [100 100];%������һ��ĺ��
realmodel=[realres,realdeep];%ʵ��ģ��
ninv=length(realmodel);% ���ݲ�������
% mmin=realmodel*0.5;%��Сģ�Ͳ���
% mmax=realmodel*1.5;%���ģ�Ͳ���


%%%%%%%%%%%%%%%%%%%%%%%%%%����ΪDC��ʼģ��%%%%%%%%%%%%%%%%%%%%%%%
% minres=[50 50 1 1 1 50 50 50 50 100 100 100];%ÿ����͵�����
% mindeep=  [5 5 5 5 5 5 5 5 5 5 5 5];%ÿ����С��ȣ���ȥ����һ�㡢����һ�㣩
% maxres=[300 300 100 100 100 200 200 200 200 400 400 400];%ÿ����ߵ�����
% maxdeep=  [20 20 10 10 10 20 20 20 20 10 10 10];%ÿ������ȣ���ȥ����һ�㣩

mmin=[minres,mindeep];%��Сģ�Ͳ���
mmax=[maxres,maxdeep];%���ģ�Ͳ���
ninv=length(mmin);% ���ݲ�������
tic;

% �����ٶ�����
ResDeepSpan=(mmax - mmin);%[������&���]�߶�
vmax=ResDeepSpan*0.2;%���׷���ٶ�  ���в���
vmin = -1 * vmax;%��С׷���ٶ�
lamd_z=0.0000001;  %��������

k=2;
[Lz]=gnerate_Lz_w(ninv);
% ��ʼ���ٶ�
v=vmax*0.5;
c1 = 2; 
c2 = 2;
eps = 1.0e-6;%����������
beta=1.5;
% ��ȡʵ��ֵ
shiceTEM=load('3cengmodelKnoisy.txt');
obstime=shiceTEM(:,1);
obshz=shiceTEM(:,2)*1;%�������ݺ������̼���ʵ��ģ�͵� ��Ӧʱ�估�ų�ǿ��,����ʵ��ֵ
obshz=obshz';
obstime=obstime';
ntime=length(obstime);%ʱ��������

% obshz=(floor(rand(1,5)*2)*2-1)*obshz*0.1+obshz;
% for i=1:ntime     %%������
%     obshz(i)=(floor(rand(1,1)*2)*2-1)*0.2*obshz(i)+obshz(i);
% end

% �������nmodel����ʼģ��   ��ÿ�����ģ�͵�������xxerror�У���������������pberror��
for imodel = 1 : nmodel  
%     v(imodel,:) = v(1,:);%��ʼ�ٶ�
    
    % �������ģ�Ͳ���ʼ������Ⱥ����ʼ�����弫ֵ
  %     [res,deep] = creatmodel(nlayer,minres,maxres,mindeep,maxdeep);%������ʼģ��[�����ʣ����]=creatmodel[��������С�����ʣ������С������]
         model(imodel,:)=creatmodel(ninv,mmin,mmax);%������ʼģ��,��ֵ��һ������
         alpha=rand;
         fxmodel(imodel,:)=alpha*(mmin+mmax)-model(imodel,:);%%%�ڴ˼ӷ���ѧϰ dy
        for i=1:ninv 
         if fxmodel(imodel,i)<mmin(i)
%              fxmodel(imodel,i)=mmin(i)+(mmax(i)-mmin(i))*alpha;
            fxmodel(imodel,i)=mmin(i)+mmax(i)-model(imodel,i);
         elseif fxmodel(imodel,i)>mmax(i)
             fxmodel(imodel,i)=mmin(i)+model(imodel,i);
         end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����΢��
%    meme(imodel,:)=model(imodel,:);
    pbest(imodel,:) = model(imodel,:);  %��������
    fxpbest(imodel,:)=fxmodel(imodel,:);
%  �����ɵ�ģ�ͽ�������
    [~,calhz]=TEMloop(ninv,model(imodel,:));%���ݼ���  ��ʼģ������
    [~,fxcalhz]=TEMloop(ninv,fxmodel(imodel,:));%���ݼ���  ��ʼģ������
   
%     ����������
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
    v(imodel,:) = pbest(imodel,:);%��ʼ�ٶ�
end  % end of imodel

% ��ȡȫ�ּ�ֵ
jter       = 1;
njter(1)   = 1;
yaver(1) = mean(xxerror);%��һ�ε�������ֵ
[gg,sub]   = min(xxerror);%��һ�ε��������Сֵ��subΪ��Сֵ���� ��n��ģ����С
yerr(1)   = gg;%��Сֵ ���ֵ
% bestmodel  = model;%�������ӵĸ������ţ���ֻ������һ�Σ�������ֵҲ����ǰ����ֵ  20�����ӵ�
ymin(jter) = pberror(sub);%��jter�ε��������Сֵ��subΪ��Сֵ�����������
gbest(1,:) = pbest(sub,:);%��һ�ε�������ģ�ͣ��������ӣ�  ȫ�����Ż�

hwait=waitbar(0,'��ʼ����');
 step=nter/100;
% ����Ⱥ�Ż��㷨����
for iter = 2 : nter
        Str=fix(iter/step);sstr=['����������',num2str(Str),'%'];waitbar(iter/nter,hwait,sstr); pause(0.05);
        r1 = rand(1,ninv);   %%�������ninv��0-1����
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
%   v(imodel,:) = (1+fai1+fai2-sqrt((fai1+fai2)*(2+fai1+fai2))) .* v(imodel,:) + fai1.*(pbest(imodel,:) - model(imodel,:)) + fai2.*(gbest(1,:) - model(imodel,:));   %%%����Ⱥ����
%  v(imodel,:) = (1+fai1+fai2-sqrt((fai1+fai2)*(2+fai1+fai2))) .* v(imodel,:) + abs( ( fai1*( pbest(imodel,:) ) - fai2*( gbest(1,:)) )) / ( fai1+fai2 );
v(imodel,:) = ( (v(imodel,:) <= vmin).*vmin  ) + ( (v(imodel,:) > vmin).*v(imodel,:) );%�ٶȳ������Ʒ�Χ������Ϊ��Χֵ��δ���޵���Ϊԭֵ
        v(imodel,:) = ( (v(imodel,:) >= vmax).*vmax  ) + ( (v(imodel,:) < vmax).*v(imodel,:) ); 
        
       %%%%%%%%%%%%%%%����Ϊ��������Ⱥλ�ø��´���%%%%%%%%%%%%%
        u = rand;
        if u>0.5
            model(imodel,:)=v(imodel,:)+( 0.5+0.5*(nter-iter)/nter ).*(v(imodel,:)-model(imodel,:)) .*log(u)+( 0.5+0.5*(nter-iter)/nter ).*(Mbest-model(imodel,:)) .*randn;%pbest(imodel,:)-
        else
            model(imodel,:)=v(imodel,:)-( 0.5+0.5*(nter-iter)/nter ).*(v(imodel,:)-model(imodel,:)) .*log(u)+( 0.5+0.5*(nter-iter)/nter ).*(Mbest-model(imodel,:)) .*randn;%pbest(imodel,:)-
        end

        
%          model(imodel,:) = model(imodel,:) + v(imodel,:);%��������λ��  �����в���������λ�ã�������ĳ������     �ڴ˼����������Ӧ�ȣ������Ƿ�ı������ٶȣ�λ�ã���
%        model(imodel,:) = pbest(imodel,:) + v(imodel,:);%��������λ��  �����в���������λ�ã�������ĳ������     �ڴ˼����������Ӧ�ȣ������Ƿ�ı������ٶȣ�λ�ã��� 
        
        min_throwaway = model(imodel,:) <= mmin;  % С����С���Ƶ�λ�� Ϊ 1
        min_keep      = model(imodel,:) >  mmin;  % ������С���Ƶ�λ�� Ϊ 1
        max_throwaway = model(imodel,:) >= mmax;  % ����������Ƶ�λ�� Ϊ 1
        max_keep      = model(imodel,:) <  mmax;  % С��������Ƶ�λ�� Ϊ 1
        
        model(imodel,:) = ( min_throwaway.*mmin ) + ( min_keep.*model(imodel,:) ); %����λ�ó������Ʒ�Χ������Ϊ��Χֵ��δ���޵���Ϊԭֵ
        model(imodel,:) = ( max_throwaway.*mmax ) + ( max_keep.*model(imodel,:) );

%        model(imodel,:)=0.1*meme(imodel,:)+0.9*model(imodel,:);
   
        v(imodel,:) = (v(imodel,:).*min_keep) + (-v(imodel,:).*min_throwaway);%�Ե���߽�������ٶȱ�ţ�ʹ�䷵�ؼ�����Ѱ
        v(imodel,:) = (v(imodel,:).*max_keep) + (-v(imodel,:).*max_throwaway);
     
        % ���ݼ���
%         res=model(imodel,1:nlayer);
%         deep=model(imodel,nlayer+1:(2*nlayer-1));
        
        [~,calhz]=TEMloop(ninv,model(imodel,:));%���ݼ���

        % ����������������ӵĸ��弫ֵ         �������ж��Ƿ��б�Ҫ�ڸ��±䶯��һ�����±䶯��һ��λ�õļ�ֵ
        xerror=(sum(((obshz-calhz)./obshz).^2))/(ntime)+lamd_z*sqrt((Lz*model(imodel,:)')'*(Lz*model(imodel,:)'));
        xxerror(imodel) = xerror;
        if  xxerror(imodel) < pberror(imodel)
            pberror(imodel) = xxerror(imodel);
            pbest(imodel,:) = model(imodel,:);
        end
    end %end of imodel  
    yaver(iter) = mean(xxerror);
       
    % ��ȡȫ�ּ�ֵ
    [ggg,sub] = min(xxerror);
    yerr(iter) = ggg;
    if ggg < gg
        gg = ggg;
        bestmodel = model;
        jter = jter + 1;
        njter(jter) = iter;
        ymin(jter) = xxerror(sub);%Ⱥ��������Ӧ��ֵ
        gbest(1,:) = model(sub,:);%Ⱥ����������       
        lamd_z=lamd_z*k;
    else
        lamd_z=lamd_z/k;
    end    
    
    % ������ֹ����
    if ymin(jter) < eps
        break;
    end

end %end of iter

 
% ������ģ�ͽ������ݼ���
[~,calhz]=TEMloop(ninv,gbest(1,:));
close(hwait);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%��ʾ��������ռ��ʱ�䳤��
% profile report
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%