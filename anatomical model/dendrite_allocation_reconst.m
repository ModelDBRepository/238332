function [d, v] = dendrite_allocation_reconst(bn,an,sd_d,sd_v,ro,lim_fact)

limiting_factor=lim_fact-1+2*rand;
lolimit=min(bn-an)-limiting_factor; % limits shortest dendrites
limiting_factor=lim_fact-1+2*rand;
dmax=max(bn)+limiting_factor; % limits maximum dorsal position
dmax=min(112,dmax);
limiting_factor=lim_fact-1+2*rand;
vmax=max(an)+limiting_factor; % limits maximum ventral position
vmax=min(112,dmax);
len_lim=6.+2*rand;
n=length(bn);

i=floor(rand*n)+1;
Z = mvnrnd([an(i) bn(i)], [sd_v^2 ro*sd_v*sd_d; ro*sd_v*sd_d sd_d^2], 1);
v=Z(1,1); d=Z(1,2);

while (d<v || d>dmax || d<0 || v>vmax  || v<0 || d-v<lolimit || d-v<len_lim)  % limits "extreme" dendrites
    Z = mvnrnd([an(i) bn(i)], [sd_v^2 ro*sd_v*sd_d; ro*sd_v*sd_d sd_d^2], 1);
    v=Z(1,1); d=Z(1,2);
end
if d < v
    error('Problem allocating dendritic field!');
end