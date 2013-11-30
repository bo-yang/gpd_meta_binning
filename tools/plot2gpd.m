function ret=plot2gpd(l1,l2,m1,m2)

x=[0:200];
y=gpd_distr(x,l1,l2);
figure;
plot(x,y,'b');
hold on

y=gpd_distr(x,m1,m2);
plot(x,y,'r');

anno1=sprintf('\\lambda_1=%0.5g,\\lambda_2=%0.5g',l1,l2);
anno2=sprintf('\\lambda_1=%0.5g,\\lambda_2=%0.5g',m1,m2);
legend(anno1,anno2)
title('GPD')

xlim([0 30])
