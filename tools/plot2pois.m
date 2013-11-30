function ret=plot2pois(l1,l2)

% plot two poisson distributions, lambda=l1 and lambda=l2

x=[0:100];
y=poisspdf(x,l1);
figure;
plot(x,y,'b');
hold on

y=poisspdf(x,l2);
plot(x,y,'r');

anno1=sprintf('\\lambda_1=%0.5g',l1);
anno2=sprintf('\\lambda_1=%0.5g',l2);
legend(anno1,anno2)
title('Poisson')

xlim([0 30])
