
% normalized heterogeneity index from Hu, Wang 2008

function H = mea_degree_heterogeneity_hu_wang(d)
n=length(d);
H=sum(sum(abs(repmat(d',n,1)-repmat(d,1,n))))/(2*(n^2)*mean(d));