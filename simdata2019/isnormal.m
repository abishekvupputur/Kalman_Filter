function Results=isnormal(x)
alpha=0.05;
n=length(x);
i=1:n;
y=sort(x);
fx=normcdf(zscore(y));
dplus=max(abs(fx-i/n));
dminus=max(abs(fx-(i-1)/n));
Dn=max(dplus,dminus);
KSz=sqrt(n)*Dn;
s=-20:1:20;
a=(-1).^s.*exp(-2*(s.*KSz).^2); 
pvalue=1-sum(a);
Results(1,1)=KSz;
Results(1,2)=pvalue;
 
% KOLMOGOROV-SMIRNOV TEST - STEPHENS MODIFICATION
 
dKSz=Dn*(sqrt(n)-0.01+0.85/sqrt(n));
 
if dKSz<0.775
    pvalue=0.15;
elseif dKSz<0.819
    pvalue=((0.10-0.15)/(0.819-0.775))*(dKSz-0.775)+0.15;
elseif dKSz<0.895
    pvalue=((0.05-0.10)/(0.895-0.819))*(dKSz-0.819)+0.10;
elseif dKSz<0.995
    pvalue=((0.025-0.05)/(0.995-0.895))*(dKSz-0.895)+0.05;
elseif dKSz<1.035
    pvalue=((0.01-0.025)/(1.035-0.995))*(dKSz-0.995)+0.025;
else
    pvalue=0.01;
end
Results(2,1)=dKSz;
Results(2,2)=pvalue;
 
% KOLMOGOROV-SMIRNOV TEST - MARSAGLIA METHOD
 
k=ceil(n*Dn);
m=2*k-1;
h=k-n*Dn;
 
Hmatrix=zeros(m,m);
 
for i=1:m-1
   for j=2:m
      if i-j+1>=0
      Hmatrix(i,j)=1/factorial(i-j+1);
    else
      Hmatrix(i,j)=0;
    end
    end
end
 
for i=1:m-1
    Hmatrix(i,1)=(1-h^i)/factorial(i);
end
 
Hmatrix(m,:)=fliplr(Hmatrix(:,1)');
 
if h<=0.5
Hmatrix(m,1)=(1 - 2*h^m)/factorial(m);
else
Hmatrix(m,1)=(1 - 2*h^m + max(0,2*h-1)^m)/factorial(m);
end
    lmax = max(eig(Hmatrix));
    Hmatrix = (Hmatrix./lmax)^n;
    pvalue = (1 - exp(gammaln(n+1) + n*log(lmax) - n*log(n)) * Hmatrix(k,k));
Results(3,1)=KSz;
Results(3,2)=pvalue;

% Compare p-value to alpha
for i=1:3
    if Results(i,2)>alpha
        Results(i,3)=1;
    else
        Results(i,3)=0;
    end
end

disp(' ') 
disp('Test Name                  Test Statistic   p-value   Normality (1:Normal,0:Not Normal)')
disp('-----------------------    --------------  ---------  --------------------------------')
fprintf('KS Limiting Form               %6.4f \t     %6.4f                 %1.0f \r',KSz,Results(1,2),Results(1,3))
fprintf('KS Stephens Modification       %6.4f \t     %6.4f                 %1.0f \r',dKSz,Results(2,2),Results(2,3))
fprintf('KS Marsaglia Method            %6.4f \t     %6.4f                 %1.0f \r',KSz,Results(3,2),Results(3,3))

 