function c=gpdcdf(x,l1,l2)

c=0;
for k=0:x
	c=c+gpdpdf(k,l1,l2);
end
