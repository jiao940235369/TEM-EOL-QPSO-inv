function inimodel = creatmodel(nlayer,mmin,mmax)

resdeep=zeros(nlayer,1);
for ilayer = 1:nlayer
    res0 = mmin(ilayer) + rand(1) * (mmax(ilayer) - mmin(ilayer));
%res0 = mmin(ilayer) + betarnd(0.8,0.8,1,1) * (mmax(ilayer) - mmin(ilayer));
    resdeep(ilayer) = floor(res0);
end
inimodel=resdeep;
