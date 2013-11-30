function [k,l] = load_kmer( file )
%load_data: read kmer occurrance file.
%   input data should be in the format of "k-mer occurrance".
%	output:
%		k - unqiue k-mers
%		l - occurrences of k-mer

f=fopen(file,'r');
a=textscan(f,'%s%f');
k=a{1};
l=a{2};

end