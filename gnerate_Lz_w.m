function [Lz]=gnerate_Lz_w(numodel)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%��Ϊxģ������xzģ��ÿһ����ɣ�����xz��z������������Լ���൱��x����������Լ��
%%%%%%%%%%%%%%%%%%%%%%%%%
%P=nmodel_average
Lz=eye(numodel)*0;
for i=2:1:numodel
    Lz(i,i)=1;
    Lz(i,i-1)=-1;
end
for i=numodel+1:numodel:numodel%ģ��ÿ�е����һ������Լ��
    Lz(i,i)=0;
    Lz(i,i-1)=0;
end
