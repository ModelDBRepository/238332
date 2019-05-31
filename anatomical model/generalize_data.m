function output_data=generalize_data(input_data)
%input_data=pri_axon_length;
[p,x] = ecdf(input_data);
w=rand;
r=find(w>p(1:end-1) & w<=p(2:end));
if p(r+1)-p(r) ~= 0
    output_data=x(r)+(x(r+1)-x(r))*(w-p(r))/(p(r+1)-p(r));
else
    output_data=x(r);
end;
clear x;
clear p;