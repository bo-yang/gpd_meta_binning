function [kmer,x,l] = load_data( file )
%load_data: read kmer occurrance file.
%   input data should be in the format of "k-mer occurrance".

f=fopen(file,'r');
a=textscan(f,'%s%f');
kmer=a{1};
[l,x]=hist(a{2},unique(a{2}));

%l=l(2:numel(l));
%x=x(2:numel(x));
figure
bar(l)

end