% A code to measure spontanous activity. Code needs the input of cells, and
% the input of backgroung, estimitad from green channel 

%this is code 2/2

close all
clear all
clc

exp_mice=[17,50,200,52,54];
saline_cont=[1000 1001 1002 1003 1004];
cno_cont=[10102 1011];


mice=exp_mice;
% mice=50;

% load(['BG_green_STRUCT'])  % if runing not for a group


    if length(mice)==length(exp_mice);
    if sum(mice-exp_mice)==0;
        condition='exp_mice';
%         BG_green_STRUCT=load(['BG_green_STRUCT_' condition])
       load(['BG_green_STRUCT_' condition])

    end
    end

    if length(mice)==length(saline_cont);
    if sum(mice-saline_cont)==0;
        condition='saline_cont';
%         BG_green_STRUCT=
       load(['BG_green_STRUCT_' condition])

    end
    end

    if length(mice)==length(cno_cont);
    if sum(mice-cno_cont)==0;
        condition='cno_cont';
%         BG_green_STRUCT=
        load(['BG_green_STRUCT_' condition])

    end
    end


    




for i=1:length(mice);

mice_table{1,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results before']); %create a struct- rows for conditions, coloumns for mice
mice_table{2,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results after']); 

end




%vec_for_locations_in_cells_per_field
for i=1:size(mice_table{1,1}.fields_for_analysis,2)
loc_vec(i)=str2num(mice_table{1,1}.fields_for_analysis(i).name(6:end))
end

field_quantities=cell(1,size(mice_table,2));

%% F_over_bg
for mouse=mice;
        ll=find(mice==mouse);
        field_ind=find(loc_vec==mouse);
        cells_per_field=mice_table{1, 1}.fields_for_analysis(field_ind).cells_per_field;
       %         cells_per_field=mice_table{1,ll}.fields_for_analysis(field_ind).cells_per_field; check
       %         this if the is a [rpblem

        field_quantities{ll}=cells_per_field;
   
end

for i=1:length(mice)
BG_green_STRUCT(1,i).cells_per_field=field_quantities{i}
end

raw_f_cells=cell(2,size(mice_table,2));
c=0;
for mouse=mice;
    c=c+1;
        ll=find(mice==mouse);
        field_ind=find(loc_vec==mouse);
      temp_data_bef=mice_table{1, 1}.fields_for_analysis(field_ind).data_before;
            temp_data_aft=mice_table{1, 1}.fields_for_analysis(field_ind).data_after;

      raw_f_cells{1,c}=(temp_data_bef)
      raw_f_cells{2,c}=(temp_data_aft)

end



four_d_mat_raw_f_cells=cell(2,size(mice_table,2));
HZ= 7;
numOfOdors = 12;  %including blank
numOfTrails = 5;
lengthOfOdor=112;
lengthOfTrails=lengthOfOdor*numOfOdors;
    

for mouse=mice;
    ll=find(mice==mouse);
    cells_raw_data_before=raw_f_cells{1,ll};
    cells_raw_data_after=raw_f_cells{2,ll};

    


    
    numOfNeurons = size( cells_raw_data_before,2); %same like data after   
    
    DataMat_bef=zeros(numOfNeurons,numOfTrails,numOfOdors,lengthOfOdor); 
    DataMat_aft=DataMat_bef;
    ii = 0; %index for the loop to allow overcoming the difference due to the fact that first 5 lines doesn't contain Neurons
   
   
    
    %% now this is one loop that organizes all the data in a 4-dim
    %matrix (with Ami's help)+ transforming to delta_f_over_f(normalizing)+ ing
    
    startOfNeurons=1;   %because it is cleaned at the begining
    for neuron = startOfNeurons:numOfNeurons + startOfNeurons - 1;
        ii = ii + 1;
        
        
        for kk = 1:numOfOdors
            for jj = 1:numOfTrails;
                startInd = (jj-1) * lengthOfTrails;%lengthOfOdor ;
                tempind = (kk-1) * lengthOfOdor + 1 : kk * lengthOfOdor;
                ind = startInd + tempind;
                DataMat_bef(ii,jj,kk,:) = cells_raw_data_before(ind,neuron);
                DataMat_aft(ii,jj,kk,:) = cells_raw_data_after(ind,neuron);

            end
           
        end
        
    end
four_d_mat_raw_f_cells{1,ll}=DataMat_bef;
four_d_mat_raw_f_cells{2,ll}=DataMat_aft;

 end




all_mice_aug_bg_bef_af=[];
all_mice_aug_bg_bef_af=cell(2,size(mice_table,2));

    for mouse=mice;
    total_bg_per_mouse_bef=[];
     total_bg_per_mouse_aft=[];
      ll=find(mice==mouse);
%       bg_data_bef=BG_green_STRUCT(ll).raw_data_mat_before;
%       bg_data_aft=BG_green_STRUCT(ll).raw_data_mat_after;
         

bg_data_bef=BG_green_STRUCT(ll).raw_data_mat_before
bg_data_aft=BG_green_STRUCT(ll).raw_data_mat_after


      cell_count=field_quantities{ll};
    for field =1:length(field_quantities{ll});
   
        aug_bg_data_bef=bg_data_bef(field,:,:,:);
        aug_bg_data_bef=repelem(aug_bg_data_bef,cell_count(field),1);
        total_bg_per_mouse_bef=[total_bg_per_mouse_bef;aug_bg_data_bef];
        
        aug_bg_data_aft=bg_data_aft(field,:,:,:);
        aug_bg_data_aft=repelem(aug_bg_data_aft,cell_count(field),1);
        total_bg_per_mouse_aft=[total_bg_per_mouse_aft;aug_bg_data_aft];
       
    end
all_mice_aug_bg_bef_af{1,ll}= total_bg_per_mouse_bef;
all_mice_aug_bg_bef_af{2,ll}= total_bg_per_mouse_aft;
    end
    
    %cells_per_field=mice_table{1, ll}.fields_for_analysis(field_ind).cells_per_field;
        
Cells_over_bg_all_mice=cell(2,size(mice_table,2));

a=(1:length(mice));
b=(1:length(mice));
for i=1:length(mice);
    a(i)=size(mice_table{1, i}.SmoothDff,1);     %Check that a and b are correct together
    b(i)=size(all_mice_aug_bg_bef_af{1, i},1);
    
end

%  FramesForF0start =prestim_and_odor_to_final*HZ-HZ*5; % use this!!!
%  FramesForF0end =prestim_and_odor_to_final*HZ-2*HZ; %change to check
 
 FramesForF0start =14; % use this!!!
 FramesForF0end =35; %change to check

Baseline=cell(2,size(mice_table,2));
std_Baseline=cell(2,size(mice_table,2));
% std_no_bg=cell(2,size(mice_table,2));
 all_bsln_bef=[];
 all_bsln_aft=[];
 
 all_bsln_subtraction_bef=[];
 all_bsln_subtraction_aft=[];
 
 all_std_Baseline_bef=[];
 all_std_Baseline_aft=[];
 
%  std_no_bg_bef=[];
%  std_no_bg_aft=[];

for i=1:size(BG_green_STRUCT,2);
   
    Cells_over_bg_all_mice{1,i}= four_d_mat_raw_f_cells{1, i}./all_mice_aug_bg_bef_af{1,i};
    Cells_over_bg_all_mice{2,i}= four_d_mat_raw_f_cells{2, i}./all_mice_aug_bg_bef_af{2,i};
    
    
    
    Baseline{1,i}=squeeze(mean(mean(mean(Cells_over_bg_all_mice{1,i}(:,:,:, FramesForF0start:FramesForF0end),4),3),2));
    Baseline{2,i}=squeeze(mean(mean(mean(Cells_over_bg_all_mice{2,i}(:,:,:, FramesForF0start:FramesForF0end),4),3),2));
    
        Cells_over_bg_all_mice_spont{1,i} =Cells_over_bg_all_mice{1,i}(:,:,:, FramesForF0start:FramesForF0end)
        Cells_over_bg_all_mice_spont{2,i} =Cells_over_bg_all_mice{2,i}(:,:,:, FramesForF0start:FramesForF0end)

    std_Baseline{1,i}=squeeze(mean(mean(std(Cells_over_bg_all_mice_spont{1,i},0,[4]),3),2));
    std_Baseline{2,i}=squeeze(mean(mean(std(Cells_over_bg_all_mice_spont{2,i},0,[4]),3),2));

%     std_no_bg{1,i}=squeeze(mean(mean(std(four_d_mat_raw_f_cells{1,i},0,[4]),3),2));
%     std_no_bg{2,i}=squeeze(mean(mean(std(four_d_mat_raw_f_cells{2,i},0,[4]),3),2));
%     
%     std_no_bg_bef=[std_no_bg_bef;std_no_bg{1,i}]
%      std_no_bg_aft=[std_no_bg_aft;std_no_bg{2,i}]


    
    all_bsln_bef=[all_bsln_bef;Baseline{1,i}];
    all_bsln_aft=[all_bsln_aft;Baseline{2,i}];
    
    all_std_Baseline_bef=[ all_std_Baseline_bef;std_Baseline{1,i}]
    all_std_Baseline_aft=[ all_std_Baseline_aft;std_Baseline{2,i}]

    
    Cells_minus_bg_all_mice{1,i}= four_d_mat_raw_f_cells{1, i}-all_mice_aug_bg_bef_af{1,i};
    Cells_minus_bg_all_mice{2,i}= four_d_mat_raw_f_cells{2, i}-all_mice_aug_bg_bef_af{2,i};
       
    
    Baseline_subtracion{1,i}=squeeze(mean(mean(mean(Cells_minus_bg_all_mice{1,i}(:,:,:, FramesForF0start:FramesForF0end),4),3),2));
    Baseline_subtracion{2,i}=squeeze(mean(mean(mean(Cells_minus_bg_all_mice{2,i}(:,:,:, FramesForF0start:FramesForF0end),4),3),2));
    
    all_bsln_subtraction_bef=[all_bsln_subtraction_bef;Baseline_subtracion{1,i}];
    all_bsln_subtraction_aft=[all_bsln_subtraction_aft;Baseline_subtracion{2,i}];
    
end

mean(all_bsln_bef)
mean(all_bsln_aft)

[h,p]=ttest(all_bsln_bef,all_bsln_aft)
[h,p]=signrank(all_bsln_bef,all_bsln_aft)


mean(all_bsln_subtraction_bef)
mean(all_bsln_subtraction_aft)

mean(all_std_Baseline_bef)
mean(all_std_Baseline_aft)

  save(['Spont_analysis_' condition],'all_bsln_bef','all_bsln_aft')

