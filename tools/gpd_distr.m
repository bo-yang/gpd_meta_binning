function y=gpd_distr(x,l1,l2)

y=zeros(size(x));

for i=1:size(x')
	y(i)=gpdpdf(x(i),l1,l2);
end
