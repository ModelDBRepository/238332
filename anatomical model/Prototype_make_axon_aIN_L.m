function [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_aIN_L(rc, dorsal_dendrite, ventral_dendrite, cell_types, this_cell_index,cell_pos) % true/is ascending
synapse_indices_asc = [];
synapse_depths_asc = [];
synapse_xs_asc = [];
synapse_indices_desc = [];
synapse_depths_desc = [];
synapse_xs_desc = [];
%%
global cell_colours;
global dendwidth;
global total_number_of_cells;
global prob_syn_low;
%
%% Environment/gradient parameters
%
del=1; % step size
betaR=0; % set to give longitudinal polarity rather than gradient+
betaV=log(10)/30;
betaD=log(10)/30;
%
%% aIN Primary projection Parameters
alphap=0.09; % weak stochasticity PRIMARY
%
gRp=0.02; % Initial values of the sensitivity values
gVp=0.02;
gDp=0.03;

gfRp=0.0539; % final values of gradient sensitivities
gfVp=0.1326;
gfDp=0.0378;
%
sbetaRp=log(10)/30; % gradient sensitivity slopes (/50->20�m; /1000->435�m
sbetaVp=log(10)/100;
sbetaDp=log(10)/100;
%
%% aIN Secondary projection Parameters
%
alphas=0.090; % stochasticity SECONDARY
%
gRs=0.0539;%Initial values of the sensitivities
gVs=0.1326;
gDs=0.0378;
%
gfRs=0.0539; % final values of gradient sensitivities
gfVs=0.1326;
gfDs=0.0378;
%
sbetaRs=log(10)/30; % gradient sensitivity slope
sbetaVs=log(10)/100; % gradient sensitivity slope
sbetaDs=log(10)/100; % gradient sensitivity slope
%
%% aIN data
%
soma_RC=rc; % fixed position for testing
soma_DV_position=[24 26 31 36 36 38 40 40 40 44 44 45 45 46 47 50 51 51 51 54 55 55 56 56 56 56 56 58 59 61 61 61 63 74]; % aIN values
pri_axon_length=[397 397 461 474 515 576 589 678 691 691 781 781 794 819 832 832 845 896 947 960 986 1024 1024 1062 1139 1152 1152 1216 1216 1229 1280 1293 1344 1344 1459 1536 1536 1843 ]; % updated 17Sep2013 from aINs-DC2D model standard.xls columnAJ%
sec_axon_length=[0  576  730  256  1088  1024  922  1306  166  90  1306  947  192  973  282  0  1344  154  691  486  346  909  0  422  384  115  614  230  128  384  1178  154  704  538  589  602  141  269  320  154  397  0  166  154  ];

%%
%
axproj=-1; % ipsilateral
axdir=-1;  % primary ascending
initial_angle=normrnd(-93,31);
branch=abs(normrnd(70,23)); % mean and SD but negative values made positive
branch_angle=normrnd(39,42);
%
ventraledgesoma_random=generalize_data(soma_DV_position);
axonlen_random_pri=floor(generalize_data(pri_axon_length))+1;
%
axonlen_random_sec=floor(generalize_data(sec_axon_length))+1;
%% Additional axon initial angle calculation
%
thetap(1)=initial_angle*pi/180;
thetas=branch_angle*pi/180;
%
%% Adjusting angles according to primary direction and side
if axproj == -1  % ipsilateral
    if axdir == -1 % primary ascending
        thetap=(sign(thetap)*pi-thetap) - ((1-abs(sign(thetap)))*pi);
    else           % primary descending
        thetas=(sign(thetas)*pi-thetas) - ((1-abs(sign(thetas)))*pi);
    end
else             % contralateral
    
    if axdir == -1 % primary ascending
        thetap=(sign(thetap)*pi-thetap) - ((1-abs(sign(thetap)))*pi);
        thetas=-thetas;
    else           % primary descending
        thetap=(sign(thetap)*pi-thetap) - ((1-abs(sign(thetap)))*pi);
        thetas=(sign(thetas)*pi-thetas) - ((1-abs(sign(thetas)))*pi);
        thetas=-thetas;
    end
end
%%
xp(1)=rc;
yp(1)= ventraledgesoma_random+25;  % offset by width of floor plate (25�m each side)
emerge=1;
np=axonlen_random_pri; %
%%
for i=2:np
    yp(i)=yp(i-1)+del*sin(thetap(i-1));    
    if abs(yp(i))>=145 || (abs(yp(i))>=137 && xp(i-1)>=495) || (abs(yp(i))<=127 && abs(yp(i)) >=125 && xp(i-1)>=700) || abs(yp(i))<=25 % the y value has crossed any barrier
        xt=xp(i-1)+axdir*del*cos(thetap(i-1));  % what xp(i) would be
        xp(i)=xp(i-1)+(del*sign(xt-xp(i-1)))+(del*(1-abs(sign(xt-xp(i-1))))); % longitudinal direction determined by direction of approach to barrier (last step)
        yp(i)=yp(i-1);                           % yp reset to previous value (before barrier crossing)
        thetap(i-1)=pi/2-(axdir*pi/2*(sign(xt-xp(i-1))))-(axdir*pi/2*(1-abs(sign(xt-xp(i-1))))); % longitudinal direction determined by direction of approach to barrier (last five steps)
    else
        xp(i)=xp(i-1)+axdir*del*cos(thetap(i-1));
    end
    thetap(i)=thetap(i-1)-((gRp(1)-gfRp(1))*exp(-sbetaRp*((i-1)-emerge))+gfRp(1))*exp(-betaR*(xp(i-1)-500))*sin(thetap(i-1))+ (sign(yp(i)))*((gVp(1)-gfVp(1))*exp(-sbetaVp*((i-1)-emerge))+gfVp(1))*exp(-betaV*(abs(yp(i-1))-5))*cos(thetap(i-1))-(sign(yp(i)))*((gDp(1)-gfDp(1))*exp(-sbetaDp*((i-1)-emerge))+gfDp(1))*exp(betaD*(abs(yp(i-1))-145))*cos(thetap(i-1))+(-alphap+2*alphap*rand);
    w=rand;
    real_x=xp(i-1);
    dis10=abs(rc-real_x);
    y0=yp(i-1);
    [distance_nearest_dend,nearest_dend_index]=min(abs(real_x-cell_pos(1:total_number_of_cells)));
    nearest_dend_pos=cell_pos(nearest_dend_index);
    
    cont_prob=prob_syn_low;%0.46
    if dis10>=10
        if nearest_dend_index >= 1 && nearest_dend_index <= total_number_of_cells
            if real_x > (nearest_dend_pos-dendwidth/2) && real_x < (nearest_dend_pos+dendwidth/2)
                if y0 <= dorsal_dendrite(nearest_dend_index) && y0 >= ventral_dendrite(nearest_dend_index)
                    if (w < cont_prob)
                        ww1=find(synapse_indices_asc == nearest_dend_index);
                        if ( isempty(ww1)& this_cell_index ~= nearest_dend_index )%isempty(ww1)
                            synapse_indices_asc = [synapse_indices_asc nearest_dend_index];
                            synapse_depths_asc = [synapse_depths_asc y0];
                            synapse_xs_asc = [synapse_xs_asc nearest_dend_pos];
                        end
                    end
                end
            end
        end
    end
end
%% plotting
lq=length(find(xp>=500));
coord_asc(1:2*lq)=reshape([xp(1:lq)' yp(1:lq)']',1,2*lq);
% plot(xp(1:lq),yp(1:lq),'Color',cell_colours(3,1:3));
% plot(xp(1),yp(1),'*');
% for i=1:length(synapse_indices_asc)
%     if i ~= this_cell_index
%         rectangle('Position',[synapse_xs_asc(i)-0.5,synapse_depths_asc(i)-0.5,1,1],'FaceColor',cell_colours(3,1:3),'Curvature',[1,1]);
%     end
% end
%% Secondary projection Parameters
%
ns=axonlen_random_sec; % secondary axon length (from data)
%
if axproj==1    % if secondary projection is commissural
    len_rand=floor(emerge+branch); % len_rand is distance to branch point (after emergence from floor plate)
else
    len_rand=floor(branch); % generates a branch point at distance from soma (from cell_param)
end
if 0<len_rand & len_rand<(np-1)    % Provided the proposed branch point is within the length of the primary axon
    xs(1)=xp(len_rand);      %initial rostro-caudal position of branch / randomly picked along the longitudinal axis
    ys(1)=yp(len_rand);      %initial dorso-ventral position of branch
elseif np==emerge   % if there is no primary axon, secondary starts at emergence of initial growth from floor plate
    xs(1)=xp(emerge);
    ys(1)=yp(emerge);
else    % If the proposed branch point is beyond the length of the primary axon
    ns=0; % no branches on very short primary axons
    xs(1)=xp(1);
    ys(1)=yp(1);
end
%
%% Secondary axon growth

for i=2:ns
    
    ys(i)=ys(i-1)+del*sin(thetas(i-1));    %
    if abs(ys(i))>=145 || (abs(ys(i))>=137 && xs(i-1)>=495) || (abs(ys(i))<=127 && abs(ys(i)) >=125 && xs(i-1)>=700) || abs(ys(i))<=25 % the y value has crossed any barrier
        xtt=xs(i-1)-axdir*del*cos(thetas(i-1));  % what xp(i) would be
        xs(i)=xs(i-1)+(del*sign(xtt-xs(i-1)))+(del*(1-abs(sign(xtt-xs(i-1))))); % longitudinal direction determined by direction of approach to barrier (last five steps)
        ys(i)=ys(i-1);                           % ys reset to previous value (before barrier crossing)
        thetas(i-1)=pi/2+(axdir*pi/2*(sign(xtt-xs(i-1))))+(axdir*pi/2*(1-abs(sign(xtt-xs(i-1))))); % longitudinal direction determined by direction of approach to barrier (last five steps)
    else
        xs(i)=xs(i-1)-axdir*del*cos(thetas(i-1));
    end
    thetas(i)=thetas(i-1)-((gRs(1)-gfRs(1))*exp(-sbetaRs*(i-1))+gfRs(1))*exp(-betaR*(xs(i-1)-500))*sin(thetas(i-1))+ (sign(ys(i)))*((gVs(1)-gfVs(1))*exp(-sbetaVs*(i-1))+gfVs(1))*exp(-betaV*(abs(ys(i-1))-5))*cos(thetas(i-1))-(sign(ys(i)))*((gDs(1)-gfDs(1))*exp(-sbetaDs*(i-1))+gfDs(1))*exp(betaD*(abs(ys(i-1))-145))*cos(thetas(i-1))+(-alphas+2*alphas*rand);
    w=rand;
    real_x=xs(i-1);
    dis10=abs(rc-real_x);
    y0=ys(i-1);
    
    [distance_nearest_dend,nearest_dend_index]=min(abs(real_x-cell_pos(1:total_number_of_cells)));
    nearest_dend_pos=cell_pos(nearest_dend_index);
    
    cont_prob=prob_syn_low;%0.46;
    if dis10>=10
        if nearest_dend_index >= 1 && nearest_dend_index <= total_number_of_cells
            if real_x > (nearest_dend_pos-dendwidth/2) && real_x < (nearest_dend_pos+dendwidth/2)
                if y0 <= dorsal_dendrite(nearest_dend_index) && y0 >= ventral_dendrite(nearest_dend_index)
                    if (w < cont_prob)
                        ww1=find(synapse_indices_desc == nearest_dend_index);
                        if ( isempty(ww1)& this_cell_index ~= nearest_dend_index )
                            synapse_indices_desc = [synapse_indices_desc nearest_dend_index];
                            synapse_depths_desc = [synapse_depths_desc y0];
                            synapse_xs_desc = [synapse_xs_desc nearest_dend_pos];
                        end
                    end
                end
            end
        end
    end
end
lqq=length(find(xs<=2000));
coord_desc(1:2*lqq)=reshape([xs(1:lqq)' ys(1:lqq)']',1,2*lqq);
% plot(xs(1:lqq),ys(1:lqq),'Color',cell_colours(2,1:3));
% plot(xs(1),ys(1),'*r');
% for i=1:length(synapse_indices_desc)
%     if i ~= this_cell_index
%         rectangle('Position',[synapse_xs_desc(i)-0.5,synapse_depths_desc(i)-0.5,1,1],'FaceColor',cell_colours(2,1:3),'Curvature',[1,1]);
%     end
% end
