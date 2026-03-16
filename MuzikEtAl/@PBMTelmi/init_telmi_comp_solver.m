function init_telmi_comp_solver(obj,file_pka)
%BUILD_UNIV_COMP_INTERPOLANT Summary of this function goes here

pka=readmatrix(file_pka,'sheet','pka');
T=pka(:,1)+273.15;
pka=pka(:,2:end);

Kafun_T=@(Tq) 10.^(-interp1(T,pka,Tq));

s1 = @(a) (1-a(:,1)).*a(:,2).*a(:,3); %1: HT-N-N (0)
s2 = @(a) prod(a,2); %2 : -T-N-N  (-1)
s3= @(a) (1-a(:,1)).*a(:,2).*(1-a(:,3));  %3 : HT-N-NH (+1)
s4= @(a) (1-a(:,1)).*(1-a(:,2)).*a(:,3); %4 : HT-NH-N (+1)
s5= @(a) a(:,1).*a(:,2).*(1-a(:,3)); %5  : -T-N-NH (0)
s6= @(a) a(:,1).*(1-a(:,2)).*a(:,3); %6  : -T-NH-N (0)
s7= @(a) prod(1-a,2); %7  : HT-NH-NH (+2)
s8= @(a) a(:,1).*(1-a(:,2)).*(1-a(:,3)); %8  : -T-NH-NH (+1)

ind={};
ind.ineutral=[1,5,6];
ind.iHdonnated=[2];
ind.nHdonnated=[1];
ind.iHaccepted=[3,4,7,8];
ind.nHaccepted=[1,1,2,1];
ind.charge=[0,-1,1,1,0,0,2,1];


afun=@(Ka,cH) Ka./(cH+Ka);
compositionfun_= @(a) [s1(a),s2(a),s3(a),s4(a),s5(a),s6(a),s7(a),s8(a)];
compfun_KacH= @(Ka,cH) compositionfun_(afun(Ka,cH));

obj.Kafun_T=Kafun_T;
obj.compfun_KacH=compfun_KacH;
obj.ind=ind;
end

