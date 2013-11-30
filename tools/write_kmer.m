function [ n ] = write_kmer( file, kmer, cluster )
%write_kmer: write kmers and corresponding clusters into file.
%   Detailed explanation goes here

f=fopen(file,'w');
for i=1:numel(cluster)
	fprintf(f,'%s\t%d\n',kmer{i},cluster(i));
end
fclose(f);

end

