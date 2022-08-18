function [Lz]=gnerate_Lz_w(numodel)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%因为x模型是由xz模型每一列组成，所以xz中z方向相邻网格约束相当于x中相邻网格约束
%%%%%%%%%%%%%%%%%%%%%%%%%
%P=nmodel_average
Lz=eye(numodel)*0;
for i=2:1:numodel
    Lz(i,i)=1;
    Lz(i,i-1)=-1;
end
for i=numodel+1:numodel:numodel%模型每行的最后一个网格不约束
    Lz(i,i)=0;
    Lz(i,i-1)=0;
end
