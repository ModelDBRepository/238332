function [syn_tab_asc syn_tab_desc axon_tab_asc axon_tab_desc] = generate_axons_and_make_synapses_L(cell_types, rc, dorsal_dendrite, ventral_dendrite)
%
global side_shift;
global total_number_of_cells;
%
syn_tab_asc=[];
syn_tab_desc=[];
axon_tab_asc=[];
axon_tab_desc=[];
%
asc=0;
desc=1;
als=2; % axon on left side
ars=3; % axon on right side
%
tt(1:8000)=0;
%
for i=1:total_number_of_cells
    cell_type = cell_types(i);
    switch cell_type
        case 1
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_rb_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if ~isempty(synapse_indices_asc)
                syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if ~isempty(synapse_indices_desc)
                syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
        case 2
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_dlc_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if ~isempty(synapse_indices_asc)
                syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' (synapse_indices_asc+side_shift)' cell_types(synapse_indices_asc+side_shift)' synapse_xs_asc' synapse_depths_asc']];
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if ~isempty(synapse_indices_desc)
                syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' (synapse_indices_desc+side_shift)' cell_types(synapse_indices_desc+side_shift)' synapse_xs_desc' synapse_depths_desc']];
            end;
        case 3
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_aIN_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if ~isempty(synapse_indices_asc)
                syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if ~isempty(synapse_indices_desc)
                syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
        case 4
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_cIN_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if ~isempty(synapse_indices_asc)
                syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' (synapse_indices_asc+side_shift)' cell_types(synapse_indices_asc+side_shift)' synapse_xs_asc' synapse_depths_asc']];
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if ~isempty(synapse_indices_desc)
                syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' (synapse_indices_desc+side_shift)' cell_types(synapse_indices_desc+side_shift)' synapse_xs_desc' synapse_depths_desc']];
            end;
        case 5
            flag_asc=0;
            if rc(i) <=850
                work_space=rand;
                if work_space <0.85
                    flag_asc=1;
                end;
            else
                if (rc(i) >850 && rc(i)<=1400)
                    work_space=rand;
                    if work_space <0.6
                        flag_asc=1;
                    end;
                end;
            end;
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_dIN_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc,flag_asc); % true/is ascending
            tt(1:8000)=0;
            tt1=[i,desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if ~isempty(synapse_indices_desc)
                syn_tab_desc= [syn_tab_desc;[i*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
            if flag_asc==1
                tt(1:8000)=0;
                tt1=[i, asc, cell_types(i), als, length(coord_asc), coord_asc];
                tt(1:length(tt1))=tt1;
                axon_tab_asc=[axon_tab_asc; tt];
                if ~isempty(synapse_indices_asc)
                    syn_tab_asc= [syn_tab_asc;[i*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];
                end;
            end;
        case 6
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_mn_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if ~isempty(synapse_indices_asc)
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if ~isempty(synapse_indices_desc)
                syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
        case 7
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_dla_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc);
            tt(1:8000)=0;
            tt1=[i, asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if ~isempty(synapse_indices_asc)
                syn_tab_asc= [syn_tab_asc;[i*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if ~isempty(synapse_indices_desc)
                syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
    end
end
%% Truncate axon tables (delete zeros)
%
axon_asc_L_nmax=max(axon_tab_asc(:,5))+5;
axon_tab_asc=axon_tab_asc(:,1:axon_asc_L_nmax);
%
axon_desc_L_nmax=max(axon_tab_desc(:,5))+5;
axon_tab_desc=axon_tab_desc(:,1:axon_desc_L_nmax);
%
end
