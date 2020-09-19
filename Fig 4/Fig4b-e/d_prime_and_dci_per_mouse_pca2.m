clc
close all
clear all

%change here according to the group you would like to see

type = 'Awake';
%  type= 'mice_awake_control';
%      type='first_timelapse'
%      type='sec_timelapse'
%     type='third_timelapse'
%    type='Control'
%   type='Experiment'
%type='Wt_cno_cont'

filename = ['Mice_table_',type,'.mat'];
load(filename)

condition=type; % for saving results
%%
AllFieldsStruct = mice_table{1,1}.fields_for_analysis;
nMice = length(mice_table);

stimframes = 49 + (1:28);

CellsPerAnimal  = zeros(1,nMice)
for ani=1:nMice;
    CellsPerAnimal(ani) = sum(AllFieldsStruct(ani).cells_per_field);   
end


%nOdors=11;  % skip Blank
%odor_inds = [1:4,6:12];

nOdors=12; odor_inds=1:12; % include Blank
nAllCells = sum(CellsPerAnimal);

All_Before = zeros(nAllCells,5,nOdors);
All_After = zeros(nAllCells,5,nOdors);

TrialAvg_Before = zeros(nOdors,nAllCells);
TrialAvg_After = zeros(nOdors,nAllCells);


Dist_within_per_mouse_Before = zeros(nOdors,nMice);
Dist_within_per_mouse_After = zeros(nOdors,nMice);
Dist_btw_per_mouse_Before = zeros(nOdors,nOdors, nMice);
Dist_btw_per_mouse_After = zeros(nOdors,nOdors, nMice);


dist_measure = @(x,y) sqrt(sum((x-y).^2));


Names = cell(1,nMice);


last_cell = 0;
for ani=1:nMice


    FieldStruct = AllFieldsStruct(ani);
    Names{ani} = FieldStruct.name;
  
    nMouseCells = CellsPerAnimal(ani);
    
     CurrentCells = last_cell + (1:nMouseCells);
     
     
    valves_names_12 ={'Valeraldehyde','Methyl propanoate','Ethyl acetate','Butyraldehyde','Blank','TMT','Female Pee', 'Male Pee', 'Peanut butter','Ethil tiglate','Propanal', 'Pups Bedding'}
   
    valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;    %%new order for valves
    ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
    
    for gg=1:length( ploting)
    k_for_plot(gg)=find( valves_12==ploting(gg));
    end 
     
    Before_raw = mice_table{1,ani}.SmoothDff;   %addition to organize the odors consistently with other places in the paper
    Before=(Before_raw(:,:, k_for_plot,:));
    [nCells,nTrials,~,nT] = size(Before);

    BeforeTrialAvg = mean(Before,2);
    BeforeStim = Before(:,:,:,stimframes);
    BeforeStimAmp = mean(BeforeStim,4); % average over stimulus period
    BeforeStimTrialAvg = squeeze(mean(BeforeStimAmp,2)); % average over trials
 

    After_raw = mice_table{2,ani}.SmoothDff;
    After=(After_raw(:,:, k_for_plot,:));
    AfterTrialAvg = mean(After,2);
    AfterStim = After(:,:,:,stimframes);
    AfterStimAmp = mean(AfterStim,4);
    AfterStimTrialAvg = squeeze(mean(AfterStimAmp,2));
    
   
    %%

    % All Responses 
    All_Before(CurrentCells,:,:) = (squeeze(BeforeStimAmp(:,:,odor_inds)));
    All_After(CurrentCells,:,:) = (squeeze(AfterStimAmp(:,:,odor_inds)));

    % Trial Avg Response per Odor, size: num Cells in Field * nOdors
    trialAvg_Odor_Resp_Before = BeforeStimTrialAvg(:,odor_inds);
    trialAvg_Odor_Resp_After = AfterStimTrialAvg(:,odor_inds);

    TrialAvg_Before(:,CurrentCells) = trialAvg_Odor_Resp_Before';
    TrialAvg_After(:,CurrentCells) = trialAvg_Odor_Resp_After';


    last_cell = last_cell + nMouseCells;

%     field_counter=field_counter+1;

    for oi=1:nOdors
        dpertrialBefore = zeros(1,5);
        dpertrialAfter = zeros(1,5);
        for ti=1:5
            % Distance within odor oi, per trial ti
            dpertrialBefore(ti) = dist_measure(BeforeStimAmp(:,ti,oi),BeforeStimTrialAvg(:,oi));
            dpertrialAfter(ti) = dist_measure(AfterStimAmp(:,ti,oi),AfterStimTrialAvg(:,oi));
        end
        
        % Avg distance within odor oi, for animal ani
        Dist_within_per_mouse_Before(oi,ani) = mean(dpertrialBefore);
        Dist_within_per_mouse_After(oi,ani) = mean(dpertrialAfter);
        
        % Distances btw pairs of odors oi, wi. For animal ani
        for wi=1:nOdors
            Dist_btw_per_mouse_Before(oi,wi,ani) = dist_measure(BeforeStimTrialAvg(:,oi),BeforeStimTrialAvg(:,wi));
            Dist_btw_per_mouse_After(oi,wi,ani) = dist_measure(AfterStimTrialAvg(:,oi),AfterStimTrialAvg(:,wi));
        end
    end
end



%% Naive cell-avg distances
Cell_avg_Dist_btw_odors_Before = zeros(nOdors,nOdors);
Cell_avg_Dist_btw_odors_After = zeros(nOdors,nOdors);

Cel_avg_Dist_within_odor_Before = zeros(nOdors,1);
Cell_avg_Dist_within_odor_After = zeros(nOdors,1);

for oi=1:nOdors
    dpertrialBefore = zeros(1,5);
    dpertrialAfter = zeros(1,5);
    for ti=1:5
        dpertrialBefore(ti) = dist_measure(All_Before(:,ti,oi),TrialAvg_Before(oi,:)');
        dpertrialAfter(ti) = dist_measure(All_After(:,ti,oi),TrialAvg_After(oi,:)');
    end
    Cel_avg_Dist_within_odor_Before(oi) = mean(dpertrialBefore);
    Cell_avg_Dist_within_odor_After(oi) = mean(dpertrialAfter);
    for wi=1:nOdors
        Cell_avg_Dist_btw_odors_Before(oi,wi) = dist_measure(TrialAvg_Before(oi,:),TrialAvg_Before(wi,:));
        Cell_avg_Dist_btw_odors_After(oi,wi) = dist_measure(TrialAvg_After(oi,:),TrialAvg_After(wi,:));
    end
end



Cell_avg_NormDist_Before = Cell_avg_Dist_btw_odors_Before./sqrt(Cel_avg_Dist_within_odor_Before*Cel_avg_Dist_within_odor_Before');
Cell_avg_NormDist_After = Cell_avg_Dist_btw_odors_After./sqrt(Cell_avg_Dist_within_odor_After*Cell_avg_Dist_within_odor_After');

figure;imagesc(Cell_avg_NormDist_Before); title('Normalized Dist Before, avg all cells')
caxis([min([Cell_avg_NormDist_Before(Cell_avg_NormDist_Before>0);Cell_avg_NormDist_After(Cell_avg_NormDist_After>0)]),max([Cell_avg_NormDist_Before(:);Cell_avg_NormDist_After(:)])])
colormap(gray)
colorbar()
pbaspect([1 1 1])
figure;imagesc(Cell_avg_NormDist_After); title('Normalized Dist After, avg all cells')
caxis([min([Cell_avg_NormDist_Before(Cell_avg_NormDist_Before>0);Cell_avg_NormDist_After(Cell_avg_NormDist_After>0)]),max([Cell_avg_NormDist_Before(:);Cell_avg_NormDist_After(:)])])
colormap(gray)
colorbar()
pbaspect([1 1 1])
%% Normalized per Animal Distances

% compute distance within for pairs of odors
matrix_of_dists_within_per_mouse_before = zeros(nOdors,nOdors, nMice);
matrix_of_dists_within_per_mouse_after = zeros(nOdors,nOdors, nMice);
for ffi=1:nMice
    % compute sqrt of product of dists within for each pair of odors
    matrix_of_dists_within_per_mouse_before(:,:,ffi) = sqrt(Dist_within_per_mouse_Before(:,ffi)*Dist_within_per_mouse_Before(:,ffi)');
    matrix_of_dists_within_per_mouse_after(:,:,ffi) = sqrt(Dist_within_per_mouse_After(:,ffi)*Dist_within_per_mouse_After(:,ffi)');
    
end

NormDist_per_mouse_Before = Dist_btw_per_mouse_Before./matrix_of_dists_within_per_mouse_before;
NormDist_per_mouse_After = Dist_btw_per_mouse_After./matrix_of_dists_within_per_mouse_after;

%% independent part- make CI for all mice seperately

CIs_all_mice= (NormDist_per_mouse_Before-NormDist_per_mouse_After)./(NormDist_per_mouse_Before+NormDist_per_mouse_After);
CIs_all_mice_reg_mean=nanmean(CIs_all_mice,3);
%create weighted average
temp=ones(size(CIs_all_mice));
for i=1:nMice
    temp(:,:,i)=temp(:,:,i)*CellsPerAnimal(i);
end
temp2=CIs_all_mice.*temp; temp2=nansum(temp2,3)/sum(CellsPerAnimal); %to get cis weighted by num of cells for everything
temp3=squeeze(nanmean(CIs_all_mice,1)); 
CIs_all_mice_weighted_mean=temp2;

%% these temps here are just for the per odor plot
temp4=[];
for i=1:nMice
    temp4(:,i)=temp3(:,i)*CellsPerAnimal(i);
end

weighted_CIS_stds_PER_ODOR=[];
for j=1:nOdors;
weighted_CIS_stds_PER_ODOR(j)=sqrt(var(temp3(j,:),CellsPerAnimal));
end



weighted_CIS_MEANS_PER_ODOR=sum(temp4,2)/sum(CellsPerAnimal);
weighted_CIS_sems_PER_ODOR=weighted_CIS_stds_PER_ODOR/sqrt(nMice);


%%
%mean_CIs_weighted=mean(CIs_all_mice_weighted_mean(CIs_all_mice_weighted_mean~=0));
vec1=triu(CIs_all_mice_weighted_mean); vec1=vec1(vec1~=0);
mean_weighted=mean(vec1);
sem_weighted=std(vec1)/sqrt(length(vec1));

vec2=triu(CIs_all_mice_reg_mean); vec2=vec2(vec2~=0); vec2=vec2(~isnan(vec2));
mean_non_weighted=mean(vec2);
sem_non_weighted=std(vec2)/sqrt(length(vec2));




non_weighted_CI=vec2; %order is oposite to be consistent with other plots, all is ok
weighted_CI=vec1;

x1 = linspace(4,5,size(non_weighted_CI,1));
x2 = linspace(6,7,size(weighted_CI,1));
mean_vec_non_weighted=mean(non_weighted_CI)*ones(size(x1))
mean_vec_weighted=mean(weighted_CI)*ones(size(x2))
x3=3:8;
y=zeros(size(x3))

figure
xlim([3 8])
hold on
plot(x1,non_weighted_CI,'k*')
plot(x2,weighted_CI,'ko')
plot(x1,mean_vec_non_weighted,'r--','linewidth',2)
plot(x2,mean_vec_weighted,'r--','linewidth',2)
plot(x3,y,'k--')
XTick=[]
ylim([-0.3 0.3])

legend('non weighted CI', 'weighted CI')
title([num2str(type) ' CI Calculated per mouse'])
sem_non_weighted_CI=std(non_weighted_CI)/sqrt(length(non_weighted_CI))
sem_weighted_CI=std(weighted_CI)/sqrt(length(weighted_CI))

x_er=[4.5 6.5];
y_er=[mean(non_weighted_CI) mean(weighted_CI)];
sems=[sem_non_weighted_CI sem_weighted_CI];
errorbar(x_er,y_er,sems,'r.','linewidth',2)


box off

[h,p]=ttest2(weighted_CI,non_weighted_CI)
[h,p]=ranksum(weighted_CI,non_weighted_CI)


%% Plot per animal Normalized Dist

 mean_d_prime_per_odor_before=zeros(nMice,nOdors);
  mean_d_prime_per_odor_after=zeros(nMice,nOdors);
  
%no place for weighted average data when checking significance according to mouse!

for ffi=1:nMice
   
    values_before = NormDist_per_mouse_Before(:,:,ffi);
    values_after = NormDist_per_mouse_After(:,:,ffi);
    all_vals = [values_before(values_before>0);values_after(values_after>0)];
    
    cmin = min(all_vals);
    cmax = max(all_vals);   
    
    
    figure;imagesc(NormDist_per_mouse_Before(:,:,ffi)); title(['Mouse ', Names{ffi}, ' Normalized Dist Before'])
    caxis([cmin, cmax])
    colorbar()
    colormap(gray)
    pbaspect([1 1 1])
    
figure;imagesc(NormDist_per_mouse_After(:,:,ffi)); title(['Mouse ', Names{ffi}, ' Normalized Dist After'])
    caxis([cmin, cmax])
    colorbar()
      colormap(gray)
      pbaspect([1 1 1])
      
      %plot DCI per mouse
      Dif=NormDist_per_mouse_After(:,:,ffi)-NormDist_per_mouse_Before(:,:,ffi);
      Baseline=NormDist_per_mouse_After(:,:,ffi)+NormDist_per_mouse_Before(:,:,ffi);
      
      DCI=(Dif./Baseline);
      idx    = isnan(DCI);
      DCI(idx) = 0;
    
      low=min(DCI(:));
      high=max(DCI(:));
      
      limit=max(abs(low),abs(high))
      
      figure;imagesc(DCI); title(['Mouse ', Names{ffi}, ' DCI'])
    caxis([-limit, limit])
    colorbar()
      colormap(gray)
      pbaspect([1 1 1])
      
   % Create stats for odors:
   for i=1:length(values_before)
  temp1=values_before(i,:); temp1=mean(temp1(temp1>0));
  temp2=values_after(i,:); temp2=mean(temp2(temp2>0));

  mean_d_prime_per_odor_before(ffi,i)=  temp1;
  mean_d_prime_per_odor_after(ffi,i)=  temp2;

  
   end
    

end

 Deltas_per_odor_all_mice=mean_d_prime_per_odor_before-mean_d_prime_per_odor_after; 
 CI_per_odor_all_mice=(mean_d_prime_per_odor_before-mean_d_prime_per_odor_after)./(mean_d_prime_per_odor_before+mean_d_prime_per_odor_after);
 
 SEMs_Deltas_per_odor=std(Deltas_per_odor_all_mice)/sqrt(size(Deltas_per_odor_all_mice,1));
 means_Deltas_per_odor=mean(Deltas_per_odor_all_mice,1);
 
 SEMs_CIs_per_odor=std(CI_per_odor_all_mice)/sqrt(size(CI_per_odor_all_mice,1));
 means_CIs_per_odor=mean(CI_per_odor_all_mice,1);
 
 
%here!
%% plot results for norm dist small sems are per mouse big sems are per odors (12)
  % plot per odor:

x1 = linspace(4,5,size(means_Deltas_per_odor,2));
x2 = linspace(6,7,size(means_CIs_per_odor,2));
mean_vec_Deltas=mean(means_Deltas_per_odor)*ones(size(x1))
mean_vec_CIs=mean(means_CIs_per_odor)*ones(size(x2))
x3=3:8;
y=zeros(size(x3))

figure
xlim([3 8])
hold on
plot(x1,means_Deltas_per_odor,'k*')
plot(x2,means_CIs_per_odor,'ko')
plot(x1,mean_vec_Deltas,'r--','linewidth',2)
plot(x2,mean_vec_CIs,'r--','linewidth',2)
plot(x3,y,'k--')
XTick=[]
sem_odor_deltas=SEMs_Deltas_per_odor;
sem_odor_CIs=SEMs_CIs_per_odor;
errorbar(x1,means_Deltas_per_odor,sem_odor_deltas,'r.','linewidth',1)
errorbar(x2,means_CIs_per_odor,sem_odor_CIs,'r.','linewidth',1)
%ylim([-0.3 0.3])

legend('DELTAS PER ODOR', 'CIs PER ODOR')
title([num2str(type)])
sem_means_Deltas_per_odor=std(means_Deltas_per_odor)/sqrt(length(means_Deltas_per_odor));
sem_means_CIs_per_odor=std(means_CIs_per_odor)/sqrt(length(means_CIs_per_odor));

x_er=[4.5 6.5];
y_er=[mean(means_Deltas_per_odor) mean(means_CIs_per_odor)];
sems=[sem_means_Deltas_per_odor sem_means_CIs_per_odor];
errorbar(x_er,y_er,sems,'r.','linewidth',2)


box off

% [h,p]=ttest2(weighted_CI,non_weighted_CI)
% [h,p]=ranksum(weighted_CI,non_weighted_CI)

%% weighted cis per odor , noise by animals, not completed
weighted_means_CIs_per_odor=mean(CI_per_odor_all_mice,1);
 weighted_SEMs_CIs_per_odor=std(CI_per_odor_all_mice)/sqrt(size(CI_per_odor_all_mice,1));

 n_cells_per_mouse_bef=zeros(1,ani);
 n_cells_per_mouse_aft=zeros(1,ani);


  
for m=1:ani
n_cells_per_mouse_bef(m)=size(mice_table{1,m}.SmoothDff,1);
n_cells_per_mouse_aft(m)=size(mice_table{2,m}.SmoothDff,1);
end



%% this part is independent. here I try to seperate artificial odors, mixtures TMT and blank  sems are for mice

Dist_blank_per_mouse_before=NormDist_per_mouse_Before(1,:,:); 
Dist_artificial_per_mouse_before=NormDist_per_mouse_Before(2:8,2:8,:); 
%Dist_TMT_per_mouse_before=NormDist_per_mouse_Before(8,:,:); 
Dist_mixtures_per_mouse_before=NormDist_per_mouse_Before(9:12,9:12,:);

Dist_blank_per_mouse_after=NormDist_per_mouse_After(1,:,:); 
Dist_artificial_per_mouse_after=NormDist_per_mouse_After(2:8,2:8,:); 
%Dist_TMT_per_mouse_after=NormDist_per_mouse_After(8,:,:); 
Dist_mixtures_per_mouse_after=NormDist_per_mouse_After(9:12,9:12,:);



for mouse=1:nMice
data_blank_bef=Dist_blank_per_mouse_before(:,:,mouse);
data_art_bef=Dist_artificial_per_mouse_before(:,:,mouse);
%data_tmt_bef=Dist_TMT_per_mouse_before(:,:,mouse);
data_mix_bef=Dist_mixtures_per_mouse_before(:,:,mouse);

data_blank_aft=Dist_blank_per_mouse_after(:,:,mouse);
data_art_aft=Dist_artificial_per_mouse_after(:,:,mouse);
%data_tmt_aft=Dist_TMT_per_mouse_after(:,:,mouse);
data_mix_aft=Dist_mixtures_per_mouse_after(:,:,mouse);

Avg_blank_per_mouse_bef(mouse)=mean(data_blank_bef(data_blank_bef>0));   %to ignore the diagonal
Avg_art_per_mouse_bef(mouse)=mean(data_art_bef(data_art_bef>0));
%Avg_tmt_per_mouse_bef(mouse)=mean(data_tmt_bef(data_tmt_bef>0));
Avg_mix_per_mouse_bef(mouse)=mean(data_mix_bef(data_mix_bef>0));

Avg_blank_per_mouse_aft(mouse)=mean(data_blank_aft(data_blank_aft>0));
Avg_art_per_mouse_aft(mouse)=mean(data_art_aft(data_art_aft>0));
%Avg_tmt_per_mouse_aft(mouse)=mean(data_tmt_aft(data_tmt_aft>0));
Avg_mix_per_mouse_aft(mouse)=mean(data_mix_aft(data_mix_aft>0));

end

Deltas_blank=Avg_blank_per_mouse_bef-Avg_blank_per_mouse_aft;
Deltas_art=Avg_art_per_mouse_bef-Avg_art_per_mouse_aft;
%Deltas_tmt=Avg_tmt_per_mouse_bef-Avg_tmt_per_mouse_aft;
Deltas_mix=Avg_mix_per_mouse_bef-Avg_mix_per_mouse_aft;

mean_deltas_blank=mean(Deltas_blank);
mean_deltas_art=mean(Deltas_art);
%mean_deltas_tmt=mean(Deltas_tmt);
mean_deltas_mix=mean(Deltas_mix);

sem_deltas_blank=std(Deltas_blank)/sqrt(nMice);
sem_deltas_art=std(Deltas_art)/sqrt(nMice);
%sem_deltas_tmt=std(Deltas_tmt)/sqrt(nMice);
sem_deltas_mix=std(Deltas_mix)/sqrt(nMice);

%x=1:4
%means=[mean_deltas_blank mean_deltas_art mean_deltas_tmt mean_deltas_mix];
%sems=[sem_deltas_blank sem_deltas_art sem_deltas_tmt sem_deltas_mix];

x=1:3;
means=[mean_deltas_blank mean_deltas_art  mean_deltas_mix];
sems=[sem_deltas_blank sem_deltas_art  sem_deltas_mix];


figure
hold on
plot(x,means,'ko')
errorbar(x,means,sems)
%title('Blank Artificial Odors TMT Mixtures')
title('Blank Artificial Odors  Mixtures')

%% independent as well- grand avg and cis (2 types) per mouse

CI_CALCULATED_PER_MOUSE=struct
grand_Avg_per_mouse_before=ones(1,nMice);
grand_Avg_per_mouse_after=ones(1,nMice);
CI_per_mouse_calculated_individually=ones(1,nMice);

for i=1:nMice
 temp1=NormDist_per_mouse_Before(:,:,i)
 temp2=NormDist_per_mouse_After(:,:,i)
 
 grand_Avg_per_mouse_before(i)=mean(temp1(temp1~=0));
 grand_Avg_per_mouse_after(i)=mean(temp2(temp2~=0));
 
 CI_CALCULATED_PER_MOUSE(i).CIs=(temp1-temp2)./(temp1+temp2)                            %here every aquare in the matrix is caculated by ci
 CI_CALCULATED_PER_MOUSE(i).Avg_CIs=mean(nanmean(CI_CALCULATED_PER_MOUSE(i).CIs))
 CI_per_mouse_calculated_individually(i)=CI_CALCULATED_PER_MOUSE(i).Avg_CIs;
end

grand_Avg_deltas=grand_Avg_per_mouse_before-grand_Avg_per_mouse_after;
CI_per_mouse_calculated_generally=(grand_Avg_per_mouse_before-grand_Avg_per_mouse_after)./(grand_Avg_per_mouse_before+grand_Avg_per_mouse_after);

grand_mean_deltas_per_mice=mean(grand_Avg_deltas);
sem_deltas_per_mice=std(grand_Avg_deltas)./sqrt(nMice);

grand_mean_CI_per_mouse_calculated_generally=mean(CI_per_mouse_calculated_generally);
sem_CI_per_mouse_calculated_generally=std(CI_per_mouse_calculated_generally)./sqrt(nMice);

grand_mean_CI_per_mouse_calculated_individually=mean(CI_per_mouse_calculated_individually);
sem_CI_per_mouse_calculated_individually=std(CI_per_mouse_calculated_individually)./sqrt(nMice);

x=1:3;
y=zeros(size(x))
means=[grand_mean_deltas_per_mice grand_mean_CI_per_mouse_calculated_generally grand_mean_CI_per_mouse_calculated_individually];
sems=[sem_deltas_per_mice sem_CI_per_mouse_calculated_generally sem_CI_per_mouse_calculated_individually];

figure
hold on
title([num2str(type) 'Grand Avgs per mouse: Deltas  CI general calculation CI individual calculation'])

xlim([0.5 3.5])
ylim([-0.3 0.3])

errorbar(x, means, sems)
plot(x,y,'k--')

snr=means./sems
%%
% Dist_blank_per_mouse_before=mean(Dist_blank_per_mouse_before) 
% Dist_artificial_per_mouse_before=NormDist_per_mouse_Before(2:7,2:7,:); 
% Dist_TMT_per_mouse_before=NormDist_per_mouse_Before(8,:,:); 
% Dist_mixtures_per_mouse_before=NormDist_per_mouse_Before(9:12,9:12,:);
%     
%     
% Deltas_per_odor_all_mice=mean_d_prime_per_odor_before-mean_d_prime_per_odor_after; 
%  CI_per_odor_all_mice=(mean_d_prime_per_odor_before-mean_d_prime_per_odor_after)./(mean_d_prime_per_odor_before+mean_d_prime_per_odor_after);
%  
%  SEMs_Deltas_per_odor=std(Deltas_per_odor_all_mice)/sqrt(size(Deltas_per_odor_all_mice,1));
%  means_Deltas_per_odor=mean(Deltas_per_odor_all_mice,1);
%  
%  SEMs_CIs_per_odor=std(CI_per_odor_all_mice)/sqrt(size(CI_per_odor_all_mice,1));
%  means_CIs_per_odor=mean(CI_per_odor_all_mice,1);
%  
%   
%  Dist_blank_per_mouse_before=mean(Dist_blank_per_mouse_before(Dist_blank_per_mouse_before>0))
  
  %%

%% Avg over Animals Normalized Distance

Avg_NormDist_Before = mean(NormDist_per_mouse_Before,3);
Avg_NormDist_After = mean(NormDist_per_mouse_After,3);
figure;imagesc(Avg_NormDist_Before); title('Normalized Dist Before, avg over mice')
 caxis([min([Avg_NormDist_Before(Avg_NormDist_Before>0);Avg_NormDist_After(Avg_NormDist_After>0)]),...
max([Avg_NormDist_Before(:);Avg_NormDist_After(:)])])
 colorbar
  colormap(gray)
  pbaspect([1 1 1])
  
figure;imagesc(Avg_NormDist_After); title('Normalized Dist After, avg over mice')
 caxis([min([Avg_NormDist_Before(Avg_NormDist_Before>0);Avg_NormDist_After(Avg_NormDist_After>0)]),...
     max([Avg_NormDist_Before(:);Avg_NormDist_After(:)])])
 colormap(gray)
  colorbar()
  pbaspect([1 1 1])
    

%% weighted Avg- give distances matrices power with accordance to the number of cells

NormDist_per_mouse_Before_weighted=[];
NormDist_per_mouse_After_weighted=[];

for ffi=1:nMice
NormDist_per_mouse_Before_weighted(:,:,ffi)=NormDist_per_mouse_Before(:,:,ffi).*CellsPerAnimal(ffi);
NormDist_per_mouse_After_weighted(:,:,ffi)=NormDist_per_mouse_After(:,:,ffi).*CellsPerAnimal(ffi);
end


Avg_NormDist_Before_weighted = (sum(NormDist_per_mouse_Before_weighted,3))/sum(CellsPerAnimal);
Avg_NormDist_After_weighted =(sum(NormDist_per_mouse_After_weighted,3))/sum(CellsPerAnimal);


figure;imagesc(Avg_NormDist_Before_weighted); title('Normalized Dist Before, weighted avg over mice')
 caxis([min([Avg_NormDist_Before_weighted(Avg_NormDist_Before_weighted>0);Avg_NormDist_After_weighted(Avg_NormDist_After_weighted>0)]),...
     max([Avg_NormDist_Before_weighted(:);Avg_NormDist_After_weighted(:)])])
 colorbar
 colormap(gray)
 pbaspect([1 1 1])
figure;imagesc(Avg_NormDist_After_weighted); title('Normalized Dist After, weighted avg over mice')
 caxis([min([Avg_NormDist_Before_weighted(Avg_NormDist_Before_weighted>0);Avg_NormDist_After_weighted(Avg_NormDist_After_weighted>0)]),...
     max([Avg_NormDist_Before_weighted(:);Avg_NormDist_After_weighted(:)])])
  colorbar()
 colormap(gray)
  pbaspect([1 1 1])
  %% Create stats for odors:
  
  mean_d_prime_per_odor_before=zeros(1,length(Avg_NormDist_Before));
  mean_d_prime_per_odor_after=zeros(1,length(Avg_NormDist_Before));
  
  mean_d_prime_per_odor_weighted_before=zeros(1,length(Avg_NormDist_Before));
  mean_d_prime_per_odor_weighted_after=zeros(1,length(Avg_NormDist_Before));
  
  for i=1:length(Avg_NormDist_Before);
  temp1=Avg_NormDist_Before(i,:); temp1=mean(temp1(temp1>0));
  temp2=Avg_NormDist_After(i,:); temp2=mean(temp2(temp2>0));
  temp3=Avg_NormDist_Before_weighted(i,:); temp3=mean(temp3(temp3>0));
  temp4=Avg_NormDist_After_weighted(i,:); temp4=mean(temp4(temp4>0));
  
  mean_d_prime_per_odor_before(i)=  temp1;
  mean_d_prime_per_odor_after(i)=  temp2;
  mean_d_prime_per_odor_weighted_before(i)=  temp3;
  mean_d_prime_per_odor_weighted_after(i)=  temp4; 
  
  end
   
  Deltas_per_odor=mean_d_prime_per_odor_before-mean_d_prime_per_odor_after; 
  CI_per_odor=(mean_d_prime_per_odor_before-mean_d_prime_per_odor_after)./(mean_d_prime_per_odor_before+mean_d_prime_per_odor_after);
 
  Deltas_per_odor_weighted=mean_d_prime_per_odor_weighted_before-mean_d_prime_per_odor_weighted_after;
  CI_per_odor_weighted=(mean_d_prime_per_odor_weighted_before-mean_d_prime_per_odor_weighted_after)./(mean_d_prime_per_odor_weighted_before+mean_d_prime_per_odor_weighted_after);
  
  %plot:

 %% %% test for all distances, first for weighted distances 

A=triu(Avg_NormDist_Before_weighted)
A(A==0)=[];
nans=isnan(A);
A(nans==1)=[];

B=triu(Avg_NormDist_After_weighted)
B(B==0)=[];
nans=isnan(B);
B(nans==1)=[];

%tests for all
[p,h]=signrank(A,B);
[h,p]=ttest(A,B);


%signrank for blank only, paired ttest for all
[h,p]=ttest(Avg_NormDist_Before_weighted,Avg_NormDist_After_weighted);
[p,h]=signrank(Avg_NormDist_Before_weighted(1,:),Avg_NormDist_After_weighted(1,:))

 %% %% test for all distances, now for regular distances 

C=triu(Avg_NormDist_Before)
C(C==0)=[];
nans=isnan(C);
C(nans==1)=[];



D=triu(Avg_NormDist_After)
D(D==0)=[];
nans=isnan(D);
D(nans==1)=[];

%tests for all
[p,h]=signrank(C,D);
[h,p]=ttest(C,D);


%signrank for blank only, paired ttest for all
[h,p]=ttest(Avg_NormDist_Before,Avg_NormDist_After);
[p,h]=signrank(Avg_NormDist_Before(1,:),Avg_NormDist_After(1,:))



%% create all sort of comparisons FOR WEIGHTED AVERAGE FIRST

Deltas_Dist_mat_weighted=Avg_NormDist_Before_weighted-Avg_NormDist_After_weighted;   ,
Deltas_Dist_vec_weighted=A-B;

Deltas_Dist_mat_norm_weighted=Deltas_Dist_mat_weighted./((Avg_NormDist_Before_weighted+Avg_NormDist_After_weighted));

nans=isnan(Deltas_Dist_mat_norm_weighted);
Deltas_Dist_mat_norm_weighted(nans==1)=0;


Deltas_Dist_vec_norm_weighted=(A-B)./((A+B));
Deltas_Dist_vec_norm=(C-D)./(C+D);

Deltas_weighted_no_norm=A-B;
Deltas_not_weighted_no_norm=C-D;

Dist_mat_ratios_weighted=Avg_NormDist_Before_weighted./Avg_NormDist_After_weighted;
Dist_vec_ratios_weighted=A./B;

Dist_vec_ratios_log=log2(Dist_vec_ratios_weighted);

% %% plot deltas
% figure
% hold on
% title('just Deltas')
% clims=[-1 1];
% imagesc(Deltas_Dist_mat,clims)
% 
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
% 
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
% axis tight
% %% plot normalized deltas
% 
% width = 2;     % Width in inches
%  height = 2;    % Height in inches
%  alw = 1;    % AxesLineWidth
%  fsz = 10;      % Fontsize
%  lw = 2;      % LineWidth
%  msz = 8;       % MarkerSize
% 
% figure
% hold on
% title('normalized Deltas')
% clims=[-0.3 0.3];
% imagesc(flip(Deltas_Dist_mat_norm_weighted),clims)
% 
% colorbar                             %cancel this for figure square plot
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
% 
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
% axis tight
% % CAL_MIR_PCA.m
% % Displaying CAL_MIR_PCA.m.
% 
% % %% plot ratios
% figure
% hold on
% title('ratios')
% clims=[0 2];
% imagesc((Dist_mat_ratios),clims)
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
% 
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
% axis tight
%%
%% per odor, noise evaluated by matrices   %CIs_all_mice_weighted_mean may be found above
odor_dist_vec_mean=zeros(size(1,nOdors));
odor_dist_vec_sems=zeros(size(1,nOdors));

CIs_general=zeros(size(Avg_NormDist_Before_weighted));
CIs_general=(Avg_NormDist_After_weighted-Avg_NormDist_Before_weighted)./(Avg_NormDist_After_weighted+Avg_NormDist_Before_weighted);


for i=1:nOdors
% odor_dist_vec_mean(i)=mean(CIs_all_mice_weighted_mean(i,CIs_all_mice_weighted_mean(i,:)~=0));
% odor_dist_vec_sems(i)=std(CIs_all_mice_weighted_mean(i,CIs_all_mice_weighted_mean(i,:)~=0))/sqrt(nOdors-1);

 odor_dist_vec_mean(i)=nanmean(CIs_general(i,:));
 odor_dist_vec_sems(i)=nanstd(CIs_general(i,:))/sqrt(nOdors-1);
 
end


%% plot log ratios
% figure
% hold on
% title('log ratios')
% clims=[-1 1];
% imagesc((log2(Dist_mat_ratios)),clims)
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
% 
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
% axis tight
%%


dis_vec_weighted=triu(Deltas_Dist_mat_norm_weighted);%get rid of zeros and nans
dis_vec_weighted=reshape(dis_vec_weighted,1,144);
dis_vec_weighted(dis_vec_weighted==0)=[];
nans=isnan(dis_vec_weighted);
dis_vec_weighted(nans==1)=[];

x=1:length(dis_vec_weighted);
y=zeros(size(x));
% figure
% hold on
% plot(x,dis_vec,'k','linewidth',2)
% plot(x,y,'r--','linewidth',2)
% xlim([1 length(dis_vec)])

%%
dis_vec_ratios_weighted=triu(Dist_mat_ratios_weighted);%get rid of zeros and nans
dis_vec_ratios_weighted=reshape(dis_vec_ratios_weighted,1,144);
dis_vec_ratios_weighted(dis_vec_ratios_weighted==0)=[];
nans=isnan(dis_vec_ratios_weighted);
dis_vec_ratios_weighted(nans==1)=[];

x=1:length(dis_vec_ratios_weighted);
y=ones(size(x));
% figure
% hold on
% plot(x,dis_vec_ratios,'k','linewidth',2)
% plot(x,y,'r--','linewidth',2)
% xlim([1 length(dis_vec_ratios)])
% 
% sum(dis_vec_ratios>=1)
% length(dis_vec_ratios)

% Deltas_Dist_vec_norm_cont=Deltas_Dist_vec_norm;   %mark as coment 

%% plot results for norm dist

weighted_CI=Deltas_Dist_vec_norm_weighted
non_weighted_CI=Deltas_Dist_vec_norm;%Deltas_Dist_vec_norm_cont       %need to fix this issue not to do this manually

x1 = linspace(4,5,size(non_weighted_CI,2));
x2 = linspace(6,7,size(weighted_CI,2));
mean_vec_cont=mean(non_weighted_CI)*ones(size(x1))
mean_vec_exp=mean(weighted_CI)*ones(size(x2))
x3=3:8;
y=zeros(size(x3))

figure
xlim([3 8])
hold on
plot(x1,non_weighted_CI,'k*')
plot(x2,weighted_CI,'ko')
plot(x1,mean_vec_cont,'r--','linewidth',2)
plot(x2,mean_vec_exp,'r--','linewidth',2)
plot(x3,y,'k--')
XTick=[]
ylim([-0.3 0.3])

legend('non weighted CI', 'weighted CI')
title([num2str(type) 'CI calculated for all together'])
sem_non_weighted_CI=std(non_weighted_CI)/sqrt(length(non_weighted_CI))
sem_weighted_CI=std(weighted_CI)/sqrt(length(weighted_CI))

x_er=[4.5 6.5];
y_er=[mean(non_weighted_CI) mean(weighted_CI)];
sems=[sem_non_weighted_CI sem_weighted_CI];
errorbar(x_er,y_er,sems,'r.','linewidth',2)


box off

[h,p]=ttest2(weighted_CI,non_weighted_CI)
[h,p]=ranksum(weighted_CI,non_weighted_CI)


%% exactly the same, but for deltas and not for CIs

weighted_deltas=Deltas_weighted_no_norm;
non_weighted_deltas=Deltas_not_weighted_no_norm;%Deltas_Dist_vec_norm_cont       %need to fix this issue not to do this manually

x1 = linspace(4,5,size(non_weighted_deltas,2));
x2 = linspace(6,7,size(weighted_deltas,2));
mean_vec_cont=mean(non_weighted_deltas)*ones(size(x1));
mean_vec_exp=mean(weighted_deltas)*ones(size(x2));
x3=3:8;
y=zeros(size(x3));

figure
xlim([3 8])
hold on
plot(x1,non_weighted_deltas,'k*')
plot(x2,weighted_deltas,'ko')
plot(x1,mean_vec_cont,'r--','linewidth',2)
plot(x2,mean_vec_exp,'r--','linewidth',2)
plot(x3,y,'k--')
XTick=[]
ylim([-1 1])
legend('non weighted Deltas', 'weighted Deltas')
title([num2str(type)])
sem_non_weighted=std(non_weighted_deltas)/sqrt(length(non_weighted_deltas))
sem_weighted=std(weighted_deltas)/sqrt(length(weighted_deltas))

x_er=[4.5 6.5];
y_er=[mean(non_weighted_deltas) mean(weighted_deltas)];
sems=[sem_non_weighted sem_weighted];
errorbar(x_er,y_er,sems,'r.','linewidth',2)

box off

[h,p]=ttest2(weighted_deltas,non_weighted_deltas)
[h,p]=ranksum(weighted_deltas,non_weighted_deltas)



if condition(1)=='C'
save(['population_results_control'],'weighted_CI','CIs_all_mice_weighted_mean','odor_dist_vec_mean','odor_dist_vec_sems','weighted_CIS_MEANS_PER_ODOR','weighted_CIS_sems_PER_ODOR','-v7.3')%,'string')
end

if condition(1)=='E'
save(['population_results_experiment'],'weighted_CI','CIs_all_mice_weighted_mean','odor_dist_vec_mean','odor_dist_vec_sems','weighted_CIS_MEANS_PER_ODOR','weighted_CIS_sems_PER_ODOR','-v7.3')%,'string')
end

if condition(1)=='A'
save(['population_results_awake'],'weighted_CI','CIs_all_mice_weighted_mean','odor_dist_vec_mean','odor_dist_vec_sems','weighted_CIS_MEANS_PER_ODOR','weighted_CIS_sems_PER_ODOR','-v7.3')%,'string')
end

if condition(1)=='m'
save(['population_results_awake_cont'],'weighted_CI','CIs_all_mice_weighted_mean','odor_dist_vec_mean','odor_dist_vec_sems','weighted_CIS_MEANS_PER_ODOR','weighted_CIS_sems_PER_ODOR','-v7.3')%,'string')
end

if condition(1)=='f'
save(['population_results_first_timelapse'],'weighted_CI','CIs_all_mice_weighted_mean','odor_dist_vec_mean','odor_dist_vec_sems','weighted_CIS_MEANS_PER_ODOR','weighted_CIS_sems_PER_ODOR','-v7.3')%,'string')
end

if condition(1)=='s'
save(['population_results_sec_timelapse'],'weighted_CI','CIs_all_mice_weighted_mean','odor_dist_vec_mean','odor_dist_vec_sems','weighted_CIS_MEANS_PER_ODOR','weighted_CIS_sems_PER_ODOR','-v7.3')%,'string')
end

if condition(1)=='t'
save(['population_results_third_timelapse'],'weighted_CI','CIs_all_mice_weighted_mean','odor_dist_vec_mean','odor_dist_vec_sems','weighted_CIS_MEANS_PER_ODOR','weighted_CIS_sems_PER_ODOR','-v7.3')%,'string')
end


