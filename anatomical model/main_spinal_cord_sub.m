function main_spinal_cord_sub(simulation_number)
% 13 December 2012
% This is the main routine to generate the Connectome
% This routine additionally contains dendritic field data
% In this version of 24 Oct 2012 a new routine Generalize_data was used
% also probability of synaptic contacts can be change in the main routine:
% low one is universal and high is for synapses from RB to dla and dlc
%
%%
%clear all
close all
%
%simulation_number should be specified from a calling routine;
global cell_colours;
global gap_between_cells;
global dendwidth;
global total_number_of_cells;
global side_shift;
global prob_syn_low;
global prob_syn_high;
prob_syn_low=0.46;% probability of synapse for all cells exceptfrom RB to dl, ususally 0.46
prob_syn_high=0.63;% probability of synapse from RB to dl, usually 0.63
%
%%
% from Bristol
% The RGB values in the current figures are: % divide 255
rb_colour = [255/255,210/255,50/255];
dlc_colour = [255/255,0/255,0/255];
ain_colour = [70/255,70/255,180/255];
cin_colour = [0/255,170/255,220/255];
din_colour = [150/255,80/255,30/255];
mn_colour = [0/255,150/255,60/255];
dla_colour = [255/255,170/255,140/255]; %dla colour code
%
cell_colours(1,1:3) = rb_colour;
cell_colours(2,1:3) = dlc_colour;
cell_colours(3,1:3) = ain_colour;
cell_colours(4,1:3) = cin_colour;
cell_colours(5,1:3) = din_colour;
cell_colours(6,1:3) = mn_colour;
cell_colours(7,1:3) = dla_colour; %dla colour call
%
%%
section_length = 2000;			     % Length of cord to simulate (in microns) - from 0.6 to 4000
gap_between_cells =1.5;	% It is a gap between cells
pos_per_bin=floor(100/gap_between_cells);
dendwidth = 1;
%
dorsal_dendrite = [];
ventral_dendrite = [];
%
%%
n_bins=15;
for j=1:7 %we do it for all cell types
    for i=1:n_bins  % we consider the body from 500 microns till 2000 microns, therefore 34 bins of 100 microns lengths
        D=0.6+0.1*(i-1);
        switch j
            %Distribution of the number of RB cells
            case 1
                if D<0.9
                    number_of_cells(i,1)=floor(15.78*D-7.92);
                elseif D>=0.9 & D<3.5
                    number_of_cells(i,1)=floor(7.11*exp(-0.53*D)+1.67);
                elseif D>=3.5 & D<3.9
                    number_of_cells(i,1)=floor(-5.87*D+22.87);
                else
                    number_of_cells(i,1)=0;
                end;
                
                %Distribution of the number of dlc cells
            case 2
                if D<1.9
                    number_of_cells(i,2)=floor(-142.06*exp(-6.0*D)+4.54);
                elseif D>=1.9 & D<3.2
                    number_of_cells(i,2)=floor(-0.0244*exp(1.665*D)+4.99);
                else
                    number_of_cells(i,2)=0;
                end;
                
                %Distribution of the number of aIN cells
            case 3
                if D<0.8
                    number_of_cells(i,3)=0;
                elseif D>=0.8 & D<3.0
                    number_of_cells(i,3)=floor(28.58*exp(-1.78*D)+2.94);
                elseif D>=3.2 & D<3.5
                    number_of_cells(i,3)=floor(-5.59*D+19.6);
                else
                    number_of_cells(i,3)=0;
                end;
                
                %Distribution of the number of cIN cells
            case 4
                if D<1.5
                    number_of_cells(i,4)=floor(-33.45*exp(-2.961*D)+15.816);
                else
                    number_of_cells(i,4)=floor(45.23*exp(-0.3047*D)-13.746);
                end;
                
                %Distribution of the number of dIN cells
            case 5
                
                %        number_of_cells(i,5)=floor(16*exp(-0.8*D))+1;
                number_of_cells(i,5)=floor(11.5-3.2*D)+1;
                
                %Distribution of the number of mn cells
            case 6
                if D<1.2
                    number_of_cells(i,6)=floor(-40.11*exp(-3.598*D)+13.46);
                elseif D<3.2
                    number_of_cells(i,6)=floor(-0.065*exp(1.665*D)+13.3);
                else
                    number_of_cells(i,6)=0;
                end;
                %Distribution of the number of dla cells
            case 7
                if D<1.2
                    number_of_cells(i,7)=0;
                elseif D>=1.2 & D<1.6
                    number_of_cells(i,7)=floor(-2.74*10^6*exp(-11.8*D)+4.95);
                elseif D>=1.6 & D<2.0
                    number_of_cells(i,7)=floor(-5.87*10^(-8)*exp(9.2*D)+5.04);
                else
                    number_of_cells(i,7)=0;
                end;
        end;
    end;
end;

total_number_of_cells = sum(sum(number_of_cells));% This is a total number of cells on one side
final_cell_number=total_number_of_cells+total_number_of_cells; % this id total number of cells on both left and rigth sides
side_shift=total_number_of_cells;  %This is decrement which we will add to cell numbers on the rigt side,
%%
% All this stuff for left side
cell_types1(1:n_bins,1:pos_per_bin)=0;
ct(1:7)=randperm(7);
for i=1:n_bins
    for j=1:7
        i1=ct(j);
        if  number_of_cells(i,i1) ~= 0
            del=(pos_per_bin/(number_of_cells(i,i1)));
            del_h=(del/2);
            for j1=1:number_of_cells(i,i1)
                pos=round(del_h+(j1-1)*del);
                in=pos;
                in1=in;
                emp=cell_types1(i,in);
                if (emp~=0)
                    k=0;
                    while emp~=0
                        k=k+1;
                        in1=in+k;
                        if in1<pos_per_bin+1
                            emp=cell_types1(i,in1);
                        end;
                        if emp~=0
                            in1=in-k;
                            if in1>0
                                emp=cell_types1(i,in1);
                            end;
                        end;
                    end
                end;%------------------------------------
                cell_types1(i,in1)=i1;
            end;
        end;
    end
end;

total_positions=n_bins*pos_per_bin;
%%
%
cell_types1=reshape(cell_types1',total_positions,1);

rc_l=find(cell_types1);
cell_types=cell_types1(rc_l);
rc_l=500+rc_l*gap_between_cells;
%
%%

% Now we repeat the same for the right side
clear cell_types1;
cell_types1(1:n_bins,1:pos_per_bin)=0;
ct(1:7)=randperm(7);%[1 4 2 6 3 7];%[4 1 2 3 6 7 ];
for i=1:n_bins
    for j=1:7
        i1=ct(j);
        if  number_of_cells(i,i1) ~= 0
            del=(pos_per_bin/(number_of_cells(i,i1)));
            del_h=(del/2);
            for j1=1:number_of_cells(i,i1)
                pos=round(del_h+(j1-1)*del);
                in=pos;%----------------------------empty? if not - adjust and put to in1
                in1=in;
                emp=cell_types1(i,in);
                if (emp~=0)
                    k=0;
                    while emp~=0
                        k=k+1;
                        in1=in+k;
                        if in1<pos_per_bin+1
                            emp=cell_types1(i,in1);
                        end;
                        if emp~=0
                            in1=in-k;
                            if in1>0
                                emp=cell_types1(i,in1);
                            end;
                        end;
                    end
                end;%------------------------------------
                cell_types1(i,in1)=i1;
            end;
        end;
    end
end;
total_positions=n_bins*pos_per_bin;
cell_types1=reshape(cell_types1',total_positions,1);
%
rc_r=find(cell_types1);
%
cell_types(1+side_shift:final_cell_number)=cell_types1(rc_r);
rc_r=500+rc_r*gap_between_cells;
%
rc=[rc_l; rc_r]';
cell_types=cell_types';
%
%%
dorsal_den_exp5=[70.3232 61.2352 77.3632 72.9856 76.672 99.456 97.024 94.9504 66.56 64.896];
ventral_den_exp5=[32.9728 50.0736 51.584 67.2256 61.5168 85.6832 71.5264 74.24 44.8512 40.1664];
dorsal_dendrites_dlc=dorsal_den_exp5';%./dlc_DSC_DVsize'*100;    %normalized
ventral_dendrites_dlc=ventral_den_exp5';%./dlc_DSC_DVsize'*100;  %normalized
%
%
%%
dorsal_den_exp4=[47.2576 54.1696 66.56 37.0176 53.76 50.176 76.8 40.6528 53.76 60.16];
ventral_den_exp4=[1.5104 4.992 4.1216 2.944 8.192 0.256 30.72 2 1.28 12.8];
dorsal_dendrites_aIN=dorsal_den_exp4';%./aIN_DSC_DVsize'*100;    %normalized
ventral_dendrites_aIN=ventral_den_exp4';%./aIN_DSC_DVsize'*100;  %normalized
%
%
%%
dorsal_den_exp2=[93 79 79 77 76 76 73 73 70 70 69 67 67 65 65 63 58 55 52 52 52 50 50 45 38 38 37 37 36 36 35 35 33 19];
ventral_den_exp2=[50 44 35 41 24 24 22 49 31 36 19 26 33 26 26 45 15 13 12 19 9 12 38 24 24 20 35 35 15 15 26 22 19 14];
dorsal_dendrites_cIN=dorsal_den_exp2';%./cIN_DSC_DVsize'*100;    %normalized
ventral_dendrites_cIN=ventral_den_exp2';%./cIN_DSC_DVsize'*100;  %normalized
%
%
%%
%
dIN_proby_table=[];
dorsal_den_gr1=[77.7700 105.2300 65.4100  92.6700 94.5700 70.4000 82.5700  48.6100   67.9700   32.7100   91.3600   89.2300   80.6000   75.5200  12.8 61.7600   61.2100   78.1000   54.3900];
ventral_den_gr1=[25.1600   53.2300   0   40.4100   0   19.2000    9.4400    4.1700   30.5300    2.3000  0   38.7000   29.4000    7.6800 5.1  38.3700  35.4900  23.0000  1.2200];
%
dorsal_den_gr2=[34.6000 54.4600   72.9600   72.7600   62.4700   76.8000   52.0000   40.3400   61.4000   54.8400   64.5300   59.9800  49.1700   74.2000   59.8800];
ventral_den_gr2=[2.6000   28.5400    2.5600   56.4300   35.3000   30.7000   20.0000    1.3600    2.6000   16.7000   13.3300   15.2400   24.1700   58.9000 28.7100];
%
dorsal_den_gr3=[87.0000   56.3000   75.5000   64.6000   61.4500   54.4100   35.8000   58.4000   37.1000   74.2000   94.7000   51.2000   32.4200   49.4900   66.6000   89.6000  32.6200   64.0000   69.1000];
ventral_den_gr3=[46.1000 10.2000 52.5000 18.2400 41.9100 37.6900 29.4000 0 24.3000   33.3000 43.5000 48.6000 15.2300 10.9800 35.8000 69.1000 14.5400 33.3000 25.6000];
%
dorsal_dendrites_dIN_gr1=dorsal_den_gr1;
dorsal_dendrites_dIN_gr2=dorsal_den_gr2;
dorsal_dendrites_dIN_gr3=dorsal_den_gr3;
%
ventral_dendrites_dIN_gr1=ventral_den_gr1;
ventral_dendrites_dIN_gr2=ventral_den_gr2;
ventral_dendrites_dIN_gr3=ventral_den_gr3;
%
%
%%
dorsal_den_exp6=[36.00291 39.40073 51.47927 46.89455 32.60509 55.22618 64.13964 48.2304 23.5008 51.072 45.184];
ventral_den_exp6=[1.536 10.30982 10.40291 12.89309 5.888 10.03055 3.095273 18.8416 1.152 1.9968 1.6128];
dorsal_dendrites_mn=dorsal_den_exp6'; %normalized
ventral_dendrites_mn=ventral_den_exp6'; %normalized
%
%%
dorsal_dendrites_dla=[59 71 86 72 68 72 64 64 68];
ventral_dendrites_dla=[46 67 75 54 41 48 61 41 54];
%
%
%% Determine dendritic fields on the left side
%
limiting_factor=5.;
%
for i=1:final_cell_number
    ct=cell_types(i);
    sd_v=15.; % Standard deviation for ventral end
    sd_d=15.;  % Standard deviation for dorsal end
    ro=0.8;
    %
    switch ct
        case(1)
            dorsal_dendrite(i) = 0;
            ventral_dendrite(i) = 0;
        case{2}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dlc,ventral_dendrites_dlc,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  137;
            ventral_dendrite(i) = (137-(dd+25))+ vd+25;
            
        case{3}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_aIN,ventral_dendrites_aIN,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  dd+25;
            ventral_dendrite(i) =  vd+25;
            
        case{4}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_cIN,ventral_dendrites_cIN,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  dd+25;
            ventral_dendrite(i) =  vd+25;
            
            
        case(5)
            if rc(i)<=850
                [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dIN_gr1,ventral_dendrites_dIN_gr1,sd_d,sd_v,ro,limiting_factor);
                dorsal_dendrite(i) =  dd+25;
                ventral_dendrite(i) =  vd+25;
            end
            
            if rc(i)>850 & rc(i)<=1400
                [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dIN_gr2,ventral_dendrites_dIN_gr2,sd_d,sd_v,ro,limiting_factor);
                dorsal_dendrite(i) =  dd+25;
                ventral_dendrite(i) =  vd+25;
            end
            
            if rc(i)>1400
                [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dIN_gr3,ventral_dendrites_dIN_gr3,sd_d,sd_v,ro,limiting_factor);
                dorsal_dendrite(i) =  dd+25;
                ventral_dendrite(i) =  vd+25;
            end
            
        case{6}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_mn,ventral_dendrites_mn,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  dd+25;
            ventral_dendrite(i) =  vd+25;
            
        case{7}
            %dorsal_dendrites_dla=137;
            %sd_d=0.;
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dla,ventral_dendrites_dla,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  137;
            ventral_dendrite(i) =  (137-(dd+25))+vd+25;
            
    end
end
%

%% Display dendrites on the left side and right side
%%  left side
% figure(1);
% hold on;
% xlabel('Distance from midbrain (\mum)','FontSize',14);
% ylabel('D-V Axis (\mum)','FontSize',14);
% title('Generated Tadpole Spinal Connectome','FontSize',14);
% caud=2000;  % gives caudal extent of graph
% axis([000 caud -145 145]);
% %set(fig2,'OuterPosition',pos2)
% rectangle('Position',[0, -145, caud, 290],'FaceColor',[255/255,255/255,204/255]);
% rectangle('Position',[0, -25, caud, 50],'FaceColor',[153/204 204/204 153/204]);
% rectangle('Position',[700, 125, caud, 2],'FaceColor',[255/255,0/255,0/255]);
% rectangle('Position',[700, -127, caud, 2],'FaceColor',[255/255,0/255,0/255]);
% rectangle('Position',[500, 137, caud, 8],'FaceColor',[255/255,210/255,50/255]);
% rectangle('Position',[500, -145, caud, 8],'FaceColor',[255/255,210/255,50/255]);
% grid off;
% %
% for i=1:total_number_of_cells
%     if cell_types(i)> 1
%         rectangle('Position',[rc(i)-dendwidth/2, ventral_dendrite(i), dendwidth, dorsal_dendrite(i) - ventral_dendrite(i)],'FaceColor',cell_colours(cell_types(i),1:3));
%     end
% end
% %
%% Determine dendritic fields on the right side
for i=side_shift+1:final_cell_number
    ct=cell_types(i);
    sd_v=15.; % Standard deviation for ventral end
    sd_d=15.;  % Standard deviation for dorsal end
    ro=0.8;
    switch ct
        case(1)
            dorsal_dendrite(i) = 0;
            ventral_dendrite(i) = 0;
        case{2}
            % sd_d=0.;
            % dorsal_dendrites_dlc=137;
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dlc,ventral_dendrites_dlc,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  137;
            ventral_dendrite(i) = (137-(dd+25))+(vd+25);
        case{3}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_aIN,ventral_dendrites_aIN,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  (dd+25);
            ventral_dendrite(i) =  (vd+25);
        case{4}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_cIN,ventral_dendrites_cIN,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  (dd+25);
            ventral_dendrite(i) =  (vd+25);
        case(5)
            if rc(i)<=850
                [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dIN_gr1,ventral_dendrites_dIN_gr1,sd_d,sd_v,ro,limiting_factor);
                dorsal_dendrite(i) =  (dd+25);
                ventral_dendrite(i) = (vd+25);
            end
            if rc(i)>850 & rc(i)<=1400
                [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dIN_gr2,ventral_dendrites_dIN_gr2,sd_d,sd_v,ro,limiting_factor);
                dorsal_dendrite(i) =  (dd+25);
                ventral_dendrite(i) =  (vd+25);
            end
            if rc(i)>1400
                [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dIN_gr3,ventral_dendrites_dIN_gr3,sd_d,sd_v,ro,limiting_factor);
                dorsal_dendrite(i) =  (dd+25);
                ventral_dendrite(i) =  (vd+25);
            end
        case{6}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_mn,ventral_dendrites_mn,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) =  (dd+25);
            ventral_dendrite(i) =  (vd+25);
        case{7}
            [dd, vd] = dendrite_allocation_reconst(dorsal_dendrites_dla,ventral_dendrites_dla,sd_d,sd_v,ro,limiting_factor);
            dorsal_dendrite(i) = 137;
            ventral_dendrite(i) = (137-(dd+25))+ (vd+25);
    end
end
%%  right side
% Generate axons and make synapses on the right side
% for i=side_shift+1:final_cell_number   %total_number_of_cells
%     if cell_types(i)> 1
%         rectangle('Position',[rc(i)-dendwidth/2, -dorsal_dendrite(i), dendwidth, (dorsal_dendrite(i) - ventral_dendrite(i))],'FaceColor',cell_colours(cell_types(i),1:3));
%     end
% end

%% ***********************************Generate axons & record synapses for both Left and Right **************************
%
[syn_tab_asc_L syn_tab_desc_L axon_tab_asc_L axon_tab_desc_L] = generate_axons_and_make_synapses_L(cell_types, rc, dorsal_dendrite, ventral_dendrite);
[syn_tab_asc_R syn_tab_desc_R axon_tab_asc_R axon_tab_desc_R] = generate_axons_and_make_synapses_R(cell_types, rc, dorsal_dendrite, ventral_dendrite);
%
all_syn=[syn_tab_asc_L;syn_tab_desc_L;syn_tab_asc_R;syn_tab_desc_R];
%
i=1:1:final_cell_number; %for i=1:final_cell_number
s=[i; cell_types(i); rc(i); dorsal_dendrite(i); ventral_dendrite(i)]';
save(['connectome files/DendriteGrad' num2str(simulation_number) '.txt'], 's', '-ascii', '-tabs');

[q1,q2]=size(all_syn);
in_con(1:final_cell_number,1:800)=0;
ks(1:final_cell_number)=-1;
for i=1:q1
    i_cur=all_syn(i,3);
    i_pre=all_syn(i,1);
    if (i_cur ~= i_pre)
        ks(i_cur)= ks(i_cur)+2;
        in_con(i_cur,ks(i_cur))=all_syn(i,1);
        in_con(i_cur,ks(i_cur)+1)=all_syn(i,2);
    end;
end;
for i=1:final_cell_number
    clear q3;
    q3=find(in_con(i,:));
    b(i,1)=length(q3)+1;
    b(i,2:ks(i)+2)=in_con(i,q3);
end;
%
[qq1,maxrow_len]=size(b);
%
c=zeros(1,maxrow_len);
c(1,1)=maxrow_len;
c=[c;b];

fid=fopen(['connectome files/inc_connectGrad' num2str(simulation_number) '.txt'],'wt');

[q1,q2]=size(c);
for i=1:q1
    fprintf(fid,'% 5.0f',c(i,1:q2));
    fprintf(fid,'\n');
end;
fclose(fid);


