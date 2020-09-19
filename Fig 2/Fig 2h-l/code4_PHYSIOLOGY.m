
%Comaprisons_control_vs_experiment
%ONE SHOULD COMMENT/UNCOMMENT THE CORRECT GROUP AT THE BEGINING IN ORDER TO
%PROCEED. THE FILES YOU NEED TO LOAD ARE IN THE PHYSIOLOGY FOLDER

clc
clear all
close all
%%

stats=struct;

%for anes                                              %Fig 2
data_set=['anes']
Cont_data=load('Mice_table_Control.mat');
Exp_data=load('Mice_table_Experiment.mat');
general_exp=load('General_results_experiment.mat');
general_cont=load('General_results_control.mat');

%                                                     %Fig 3
% % % %for awake
% data_set=['awake']
% Cont_data=load('Mice_table_mice_awake_control.mat');
% Exp_data=load('Mice_table_awake.mat');
% general_exp=load('General_results_awake.mat');
% general_cont=load('General_results_awake_cont.mat');
                                                          %Fig S3
%for cno control
% Cont_data=load('Mice_table_Wt_cno_cont.mat');
% Exp_data=load('Mice_table_awake.mat');
% general_exp=load('General_results_awake.mat');
% general_cont=load('General_results_Wt_cno_cont.mat');


                                                          %Fig 5
%% DONT USE THIS, BUT THE OTHER CODE THAT IS PROVIDED FOR TIME LAPSE.
% %for time lapse
% data_set=['time_lapse']
% Cont_data=load('Mice_table_sec_timelapse.mat');
% Exp_data=load('Mice_table_first_timelapse.mat');
% general_exp=load('General_results_first_timelapse.mat');
% general_cont=load('General_results_sec_timelapse.mat');
% 
% % %for time lapse
% Cont_data=load('Mice_table_third_timelapse.mat');
% Exp_data=load('Mice_table_sec_timelapse.mat');
% general_exp=load('General_results_sec_timelapse.mat');
% general_cont=load('General_results_third_timelapse.mat');

%%

temp_for_cont=size(Cont_data.mice_table);
Num_of_cont_mice=temp_for_cont(2);
Num_of_cont_conditions=temp_for_cont(1);
Num_of_odors=size(Exp_data.mice_table{1,1}.meta_desicion,2);

temp_for_exp=size(Exp_data.mice_table);
Num_of_exp_mice=temp_for_exp(2);
Num_of_exp_conditions=temp_for_exp(1);

 valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;    %%new orfer for valves
    ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
    
    for gg=1:length( ploting)
    k_for_plot(gg)=find( valves_12==ploting(gg));
    end


%% Histogramas for number of responses per odor
exp_resp_per_odor_before=zeros(1,Num_of_odors);
exp_resp_per_odor_after=zeros(1,Num_of_odors);
signif_exp_matrix_before=[];
siginf_exp_matrix_after=[];
% 
% prob_vec_bef=[];
% prob_vec_aft=[];

for j=1:Num_of_exp_mice
    
a=sum(Exp_data.mice_table{1,j}.meta_desicion~=0);
b=sum(Exp_data.mice_table{2,j}.meta_desicion~=0);
signif_exp_matrix_before=[signif_exp_matrix_before;Exp_data.mice_table{1,j}.meta_desicion~=0];
siginf_exp_matrix_after=[siginf_exp_matrix_after;Exp_data.mice_table{2,j}.meta_desicion~=0];
exp_resp_per_odor_before=exp_resp_per_odor_before+a;
exp_resp_per_odor_after=exp_resp_per_odor_after+b;
% prob_to_resp_per_mouse_bef=sum(a)/size(Exp_data.mice_table{1,j}.meta_desicion,1);
% prob_to_resp_per_mouse_aft=sum(b)/size(Exp_data.mice_table{2,j}.meta_desicion,1);
% prob_vec_bef=[prob_vec_bef;prob_to_resp_per_mouse_bef];
% prob_vec_aft=[prob_vec_aft;prob_to_resp_per_mouse_aft];
end

% figure
% hold on
% x=(1:2)
% y=[mean(prob_vec_bef),mean(prob_vec_aft)];
% sem_bef=std(prob_vec_bef)/sqrt(size(prob_vec_bef,1));
% sem_aft=std(prob_vec_aft)/sqrt(size(prob_vec_aft,1));
% sems=[sem_bef sem_aft]
% bar(x(1),y(1),'b')
% bar(x(2),y(2),'r')
% errorbar(x,y,sems,'.')
% %errorbar(x(2),sems(2),'.')
% ranksum(prob_vec_bef,prob_vec_aft)


exp_matrix_before_organized=signif_exp_matrix_before(:,k_for_plot);
exp_matrix_after_organized=siginf_exp_matrix_after(:,k_for_plot);
exp_resp_per_odor_before_organized=exp_resp_per_odor_before(:,k_for_plot);
exp_resp_per_odor_after_organized=exp_resp_per_odor_after(:,k_for_plot);
exp_resp_per_odor_before_organized(:,1)=[];
exp_resp_per_odor_after_organized(:,1)=[];
exp_matrix_before_organized(:,1)=[];
exp_matrix_after_organized(:,1)=[];
figure
bar(exp_resp_per_odor_before_organized,0.6,'histc','b')
hold on
bar (exp_resp_per_odor_after_organized,0.6,'r')
box off
xlim([0.5 12])
ylim([0 350])
title('Number of responses per odor-experiment')

figure   %normalized- prob
bar(exp_resp_per_odor_before_organized./size(exp_matrix_before_organized,1),0.6,'histc','b')
hold on
bar (exp_resp_per_odor_after_organized./size(exp_matrix_after_organized,1),0.6,'r')
box off
xlim([0.5 12])
ylim([0 0.7])
title('Proportion of responses per odor-experiment')

% title('Number of responses per odor-Experiment')
% legend('Before', 'After')
% xlabel('Number of odor')
% ylabel('Number of responses')

%Same for Control:
cont_resp_per_odor_before=zeros(1,Num_of_odors);
cont_resp_per_odor_after=zeros(1,Num_of_odors);
signif_cont_matrix_before=[];
signif_cont_matrix_after=[];
for j=1:Num_of_cont_mice
    
a=sum(Cont_data.mice_table{1,j}.meta_desicion~=0);
b=sum(Cont_data.mice_table{2,j}.meta_desicion~=0);
signif_cont_matrix_before=[signif_cont_matrix_before;Cont_data.mice_table{1,j}.meta_desicion~=0];
signif_cont_matrix_after=[signif_cont_matrix_after;Cont_data.mice_table{2,j}.meta_desicion~=0];
cont_resp_per_odor_before=cont_resp_per_odor_before+a;
cont_resp_per_odor_after=cont_resp_per_odor_after+b;
end

cont_matrix_before_organized=signif_cont_matrix_before(:,k_for_plot);
cont_matrix_after_organized=signif_cont_matrix_after(:,k_for_plot);
cont_resp_per_odor_before_organized=cont_resp_per_odor_before(:,k_for_plot);
cont_resp_per_odor_after_organized=cont_resp_per_odor_after(:,k_for_plot);
cont_resp_per_odor_before_organized(:,1)=[];
cont_resp_per_odor_after_organized(:,1)=[];
figure
bar(cont_resp_per_odor_before_organized,0.6,'histc','b')
hold on
bar (cont_resp_per_odor_after_organized,0.6,'r')
xlim([0.5 12])
ylim([0 350])

title('Number of responses per odor-control')
% title('Number of responses per odor-Control')
% legend('Before', 'After')
% xlabel('Number of odor')
% ylabel('Number of responses')
box off


figure   %normalized- prob
bar(cont_resp_per_odor_before_organized./size(cont_matrix_before_organized,1),0.6,'histc','b')
hold on
bar (cont_resp_per_odor_after_organized./size(cont_matrix_after_organized,1),0.6,'r')
box off
xlim([0.5 12])
ylim([0 0.7])
title('Proportion of responses per odor-control')

%% all trials magnitudes together histograms for magnitude- sorted and not sorted by oddor strength
Total_exp_mat_for_trials_magnitudes_before=[];
Total_exp_mat_for_trials_magnitudes_after=[];

Total_cont_mat_for_trials_magnitudes_before=[];
Total_cont_mat_for_trials_magnitudes_after=[];

for j=1:Num_of_exp_mice
 Total_exp_mat_for_trials_magnitudes_before=[Total_exp_mat_for_trials_magnitudes_before;Exp_data.mice_table{1,j}.mat_for_diffrence_test] ; 
 Total_exp_mat_for_trials_magnitudes_after=[Total_exp_mat_for_trials_magnitudes_after;Exp_data.mice_table{2,j}.mat_for_diffrence_test];  
end

for j=1:Num_of_cont_mice
 Total_cont_mat_for_trials_magnitudes_before=[Total_cont_mat_for_trials_magnitudes_before;Cont_data.mice_table{1,j}.mat_for_diffrence_test];  
 Total_cont_mat_for_trials_magnitudes_after=[Total_cont_mat_for_trials_magnitudes_after;Cont_data.mice_table{2,j}.mat_for_diffrence_test];  
end

%added here down : for each not to reduce dimension

 Total_exp_mat_for_trials_magnitudes_before_organized= Total_exp_mat_for_trials_magnitudes_before(:,:,k_for_plot);
 Total_exp_mat_for_trials_magnitudes_after_organized= Total_exp_mat_for_trials_magnitudes_after(:,:,k_for_plot);
 Total_cont_mat_for_trials_magnitudes_before_organized= Total_cont_mat_for_trials_magnitudes_before(:,:,k_for_plot);
 Total_cont_mat_for_trials_magnitudes_after_organized= Total_cont_mat_for_trials_magnitudes_after(:,:,k_for_plot);
 
 Total_exp_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
 Total_exp_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
 Total_cont_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
 Total_cont_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
 
% Total_exp_mat_for_trials_magnitudes_before_organized=abs(Total_exp_mat_for_trials_magnitudes_before_organized);  %delete only this part to get them not as abs
% Total_exp_mat_for_trials_magnitudes_after_organized=abs(Total_exp_mat_for_trials_magnitudes_after_organized);
% Total_cont_mat_for_trials_magnitudes_before_organized=abs(Total_cont_mat_for_trials_magnitudes_before_organized);
% Total_cont_mat_for_trials_magnitudes_after_organized=abs(Total_cont_mat_for_trials_magnitudes_after_organized);
 %%for exp not sorted by magnitude for each cell:
 
  x=1:11;
  y1=(mean(Total_exp_mat_for_trials_magnitudes_before_organized,2));
  std_y1=squeeze(std(y1,1));
  sem_y1=std_y1/sqrt(length(y1));
  y1=mean(y1,1);
  y1=squeeze(y1);

  y2=(mean(Total_exp_mat_for_trials_magnitudes_after_organized,2))
  std_y2=squeeze(std(y2,1));
  sem_y2=std_y1/sqrt(length(y2));
  y2=mean(y2,1);
  y2=squeeze(y2);
  
 mono_dif=y1(1:6)-y2(1:6);
 mono_dif_prop=mono_dif./y1(1:6);
 nat_dif=y1(7:11)-y2(7:11);
 nat_dif_prop=nat_dif./y1(7:11);
 
 [h,p]=ranksum(mono_dif_prop,nat_dif_prop);
 [h,p]=ttest2(mono_dif_prop,nat_dif_prop);
 
 figure
 errorbar(x,y1,sem_y1,'b','linewidth',2)
 hold on
 errorbar(x,y2,sem_y2,'r','linewidth',2)
 title('magnitude of response per odor-experiment')
 xlabel('Odor identity')
 ylabel('response magnitude')
  box off
  
  ylim([-0.1 0.1])
% ylim([-0.05 0.4])
 %ylim([0 0.2])  %for integral
 %%for cont not sorted by magnitude for each cell:
 
 x=1:11;
  y3=(mean(Total_cont_mat_for_trials_magnitudes_before_organized,2));
  std_y3=squeeze(std(y3,1));
  sem_y3=std_y3/sqrt(length(y3));
  y3=mean(y3,1);
  y3=squeeze(y3);

  y4=(mean(Total_cont_mat_for_trials_magnitudes_after_organized,2));
  std_y4=squeeze(std(y4,1));
  sem_y4=std_y4/sqrt(length(y4));
  y4=mean(y4,1);
  y4=squeeze(y4);
  

 
 
  figure
 errorbar(x,y3,sem_y3,'b','linewidth',2)
 hold on
 errorbar(x,y4,sem_y4,'r','linewidth',2)
 box off
 %ylim([-0.05 0.4])
  ylim([-0.1 0.1])
%  ylim([0 0.2])  %for integral
 title('magnitude of response per odor-control')
  xlabel('Odor identity')
 ylabel('response magnitude')
 % now sorted by magnitude:
 
 %  %sort by magnitude:
 sort_exp_bef= Total_exp_mat_for_trials_magnitudes_before_organized;
 sort_exp_aft=Total_exp_mat_for_trials_magnitudes_after_organized;
 sort_cont_bef=Total_cont_mat_for_trials_magnitudes_before_organized;
 sort_cont_aft=Total_cont_mat_for_trials_magnitudes_after_organized;
 
 [sort_exp_bef,Indices]=sort(squeeze(mean(sort_exp_bef,2)),2);
  sort_exp_aft=sort(squeeze(mean( sort_exp_aft,2)),2);
  sort_exp_aft_paired=sort_exp_aft;
  sort_exp_aft_paired(Indices)= sort_exp_aft_paired;
  
  [sort_cont_bef,Indices]=sort(squeeze(mean(sort_cont_bef,2)),2);
  sort_cont_aft=sort(squeeze(mean( sort_cont_aft,2)),2);
 sort_cont_aft_paired=sort_cont_aft;
  sort_cont_aft_paired(Indices)= sort_cont_aft_paired;
  
  %%added at 8.6.20 to relate to change is spontanous activity
%  save(['responses_per_cell_exp'], 'sort_exp_bef','sort_exp_aft','sort_exp_aft_paired')
%  %save(['responses_per_cell_cont_saline'], 'sort_cont_bef','sort_cont_aft','sort_cont_aft_paired')
%  save(['responses_per_cell_cont_cno'], 'sort_cont_bef','sort_cont_aft','sort_cont_aft_paired')

  
  
 mean_sort_exp_bef=mean(sort_exp_bef,1);
 mean_sort_exp_aft=mean(sort_exp_aft,1);
 mean_sort_cont_bef=mean(sort_cont_bef,1);
 mean_sort_cont_aft=mean(sort_cont_aft,1); 
%% 28.8.19 curves ratios analysis- stats used for paper

std( mean_sort_exp_bef)/ std( mean_sort_exp_aft); %itamars measure for curvature- the whole population
std_ratio_per_cell_exp=std(sort_exp_bef')./std(sort_exp_aft');  %a vector of sds for each cell
mean_curves_ratio_exp=mean(std(sort_exp_bef')./std(sort_exp_aft')); %mean curve;
sem_curves_ratio_exp=(std(std_ratio_per_cell_exp))/sqrt(length(std_ratio_per_cell_exp)); %sem for curves


std( mean_sort_cont_bef)/ std( mean_sort_cont_aft); %itamars measure for curvature- the whole population
std_ratio_per_cell_cont=std(sort_cont_bef')./std(sort_cont_aft');  %a vector of sds for each cell
mean_curves_ratio_cont=mean(std(sort_cont_bef')./std(sort_cont_aft')); %mean curve;
sem_curves_ratio_cont=(std(std_ratio_per_cell_cont))/sqrt(length(std_ratio_per_cell_cont)); %sem for curves

means=[mean_curves_ratio_exp mean_curves_ratio_cont]
sems=[sem_curves_ratio_exp sem_curves_ratio_cont]
figure
errorbar(means,sems)

[h,p,stat]=ranksum(std_ratio_per_cell_exp,std_ratio_per_cell_cont) %%%the diffrence between curves- DO NOT report this! 
[p,h,stat]=signrank(std_ratio_per_cell_exp,1)
[p,h,stat]=signrank(std_ratio_per_cell_cont,1)

%stats for paper
cont_bef=std(sort_cont_bef');        %REPORT THESE
cont_aft=std(sort_cont_aft');
[h,p,stat]=signrank(cont_bef,cont_aft)

exp_bef=std(sort_exp_bef');          %REPORT THESE
exp_aft=std(sort_exp_aft');
[h,p,stat]=signrank(exp_bef,exp_aft)
%[h,p,ci,stat]=ttest(exp_bef,exp_aft)

%only for cno control- c
% 
% [p,h,ci,stat]=ttest(log(std_ratio_per_cell_cont))
% [p,h,ci,stat]=ttest(log(std_ratio_per_cell_exp))

disp('STOP! here is the curve analysis')
%%
 
 std_sort_exp_bef=std(sort_exp_bef,1)/sqrt(length(sort_exp_bef));
 std_sort_exp_aft=std(sort_exp_aft,1)/sqrt(length(sort_exp_aft));
 std_sort_cont_bef=std(sort_cont_bef,1)/sqrt(length(sort_cont_bef));
 std_sort_cont_aft=std(sort_cont_aft,1)/sqrt(length(sort_cont_aft)); 
  
 
  %exp:
 figure
errorbar (x,mean_sort_exp_bef,std_sort_exp_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_exp_aft,std_sort_exp_aft,'r','linewidth',2)
title('responses sorted by magnitude-experiment')
 box off
 ylim([-0.15 0.5])
%  ylim([-0.1 0.45]) for anes
 small_dif=mean_sort_exp_bef(4:7)-mean_sort_exp_aft(4:7);
 small_dif_prop=small_dif./mean_sort_exp_bef(4:7);
 big_dif=mean_sort_exp_bef(8:11)-mean_sort_exp_aft(8:11);
 big_dif_prop=big_dif./mean_sort_exp_bef(8:11);
 
 [h,p]=ranksum(small_dif_prop,big_dif_prop);
 [h,p]=ttest2(small_dif_prop,big_dif_prop);
 
 %cont:
   figure
errorbar (x,mean_sort_cont_bef,std_sort_cont_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_cont_aft,std_sort_cont_aft,'r','linewidth',2)
 title('responses sorted by magnitude-control')
box off
 ylim([-0.15 0.5])
 
 

 
 
 
%  ylim([-0.1 0.45]) for anes
%   ylim([-0.15 0.25]) for awake

%% Deltas_before_after_for_strongest+weakest responses

exp_deltas_first=(sort_exp_aft(:,1)-sort_exp_bef(:,1));
exp_deltas_last=(sort_exp_aft(:,11)-sort_exp_bef(:,11));

cont_deltas_first=(sort_cont_aft(:,1)-sort_cont_bef(:,1));
cont_deltas_last=(sort_cont_aft(:,11)-sort_cont_bef(:,11));

Exp_CI_first=((sort_exp_aft(:,1)-sort_exp_bef(:,1)))./(abs(sort_exp_aft(:,1))+abs(sort_exp_bef(:,1)));
Exp_CI_last=((sort_exp_aft(:,11)-sort_exp_bef(:,11)))./(abs(sort_exp_aft(:,11))+abs(sort_exp_bef(:,11)));


cont_CI_first=((sort_cont_aft(:,1)-sort_cont_bef(:,1)))./(abs(sort_cont_aft(:,1))+abs(sort_cont_bef(:,1)));
cont_CI_last=((sort_cont_aft(:,11)-sort_cont_bef(:,11)))./(abs(sort_cont_aft(:,11))+abs(sort_cont_bef(:,11)));


sem_reg_exp_first=std(exp_deltas_first)/sqrt(length(exp_deltas_first));
sem_reg_cont_first=std(cont_deltas_first)/sqrt(length(cont_deltas_first));


sem_reg_exp_last=std(exp_deltas_last)/sqrt(length(exp_deltas_last));
sem_reg_cont_last=std(cont_deltas_last)/sqrt(length(cont_deltas_last));

sems_reg_first=[sem_reg_exp_first sem_reg_cont_first]
sems_reg_last=[sem_reg_exp_last sem_reg_cont_last]

means_reg_first=[mean(exp_deltas_first) mean(cont_deltas_first)];
means_reg_last=[mean(exp_deltas_last) mean(cont_deltas_last)];



figure
hold on

x=(1:2);

bar(x(1),means_reg_first(1),'k')
bar(x(2),means_reg_first(2),'k')
errorbar(x,means_reg_first,sems_reg_first,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

ylim([min([means_reg_first means_reg_last])-max([sems_reg_first sems_reg_last]) max([means_reg_first means_reg_last])+max([sems_reg_first sems_reg_last])])

title('Response magnitudes Deltas per cell, regular-sorted first. cont vs exp')

%
% x1 = linspace(4,5,size(exp_deltas_first,1));
% x2 = linspace(6,7,size(cont_deltas_first,1));
% 
% mean_vec_exp=mean(means_reg_first(1))*ones(size(x1))
% mean_vec_cont=mean(means_reg_first(2))*ones(size(x2))
% 
% x3=3:8;
% y=zeros(size(x3))
% 
% figure
% title('Norm deltas regular distances')
% xlim([3 8])
% hold on
% plot(x1,exp_deltas_first,'k*')
% plot(x2,cont_deltas_first,'ko')
% plot(x1,mean_vec_exp,'r--','linewidth',2)
% plot(x2,mean_vec_cont,'r--','linewidth',2)
% plot(x3,y,'k--')
% XTick=[]
% 
% sem_exp=sem_reg_exp_first
% sem_cont=sem_reg_cont_first
% 
% x_er=[4.5 6.5];
% y_er=[means_reg_first(1) means_reg_first(2)];
% sems=[sem_cont sem_exp];
% errorbar(x_er,y_er,sems,'r.','linewidth',2)
% 
% box off
%





figure
hold on


bar(x(1),means_reg_last(1),'k')
bar(x(2),means_reg_last(2),'k')
errorbar(x,means_reg_last,sems_reg_last,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]
ylabel('dff')
xlabel('Exp              Cont')
title('Response magnitudes Deltas per cell, regular-sorted last. cont vs exp')
ylim([min([means_reg_first means_reg_last])-max([sems_reg_first sems_reg_last]) max([means_reg_first means_reg_last])+max([sems_reg_first sems_reg_last])])

%now for CIs
Exp_CI_first=((sort_exp_aft(:,1)-sort_exp_bef(:,1)))./(abs(sort_exp_aft(:,1))+abs(sort_exp_bef(:,1)));
Exp_CI_last=((sort_exp_aft(:,11)-sort_exp_bef(:,11)))./(abs(sort_exp_aft(:,11))+abs(sort_exp_bef(:,11)));


cont_CI_first=((sort_cont_aft(:,1)-sort_cont_bef(:,1)))./(abs(sort_cont_aft(:,1))+abs(sort_cont_bef(:,1)));
cont_CI_last=((sort_cont_aft(:,11)-sort_cont_bef(:,11)))./(abs(sort_cont_aft(:,11))+abs(sort_cont_bef(:,11)));


sem_CI_exp_first=std(Exp_CI_first)/sqrt(length(Exp_CI_first));
sem_CI_cont_first=std(cont_CI_first)/sqrt(length(cont_CI_first));


sem_CI_exp_last=std(Exp_CI_last)/sqrt(length(Exp_CI_last));
sem_CI_cont_last=std(cont_CI_last)/sqrt(length(cont_CI_last));

sems_CI_first=[sem_CI_exp_first sem_CI_cont_first];
sems_CI_last=[sem_CI_exp_last sem_CI_cont_last];

means_CI_first=[mean(Exp_CI_first) mean(cont_CI_first)];
means_CI_last=[mean(Exp_CI_last) mean(cont_CI_last)];

figure
hold on

x=(1:2);

bar(x(1),means_CI_first(1),'k')
bar(x(2),means_CI_first(2),'k')
errorbar(x,means_CI_first,sems_CI_first,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

ylim([min([means_CI_first means_CI_last])-max([sems_CI_first sems_CI_last]) max([means_CI_first means_CI_last])+max([sems_CI_first sems_CI_last])])

title('Response magnitudes Deltas per cell, CIs-sorted first. cont vs exp')


figure
hold on


bar(x(1),means_CI_last(1),'k')
bar(x(2),means_CI_last(2),'k')
errorbar(x,means_CI_last,sems_CI_last,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

ylim([min([means_CI_first means_CI_last])-max([sems_CI_first sems_CI_last]) max([means_CI_first means_CI_last])+max([sems_CI_first sems_CI_last])])

title('Response magnitudes Deltas per cell, CIs-sorted last. cont vs exp')

%% same like last but by absolute values

%% all trials magnitudes together histograms for magnitude- sorted and not sorted by oddor strength
Total_exp_mat_for_trials_magnitudes_before=[];
Total_exp_mat_for_trials_magnitudes_after=[];

Total_cont_mat_for_trials_magnitudes_before=[];
Total_cont_mat_for_trials_magnitudes_after=[];

for j=1:Num_of_exp_mice
 Total_exp_mat_for_trials_magnitudes_before=[Total_exp_mat_for_trials_magnitudes_before;Exp_data.mice_table{1,j}.mat_for_diffrence_test] ; 
 Total_exp_mat_for_trials_magnitudes_after=[Total_exp_mat_for_trials_magnitudes_after;Exp_data.mice_table{2,j}.mat_for_diffrence_test];  
end

for j=1:Num_of_cont_mice
 Total_cont_mat_for_trials_magnitudes_before=[Total_cont_mat_for_trials_magnitudes_before;Cont_data.mice_table{1,j}.mat_for_diffrence_test];  
 Total_cont_mat_for_trials_magnitudes_after=[Total_cont_mat_for_trials_magnitudes_after;Cont_data.mice_table{2,j}.mat_for_diffrence_test];  
end

%added here down : for each not to reduce dimension

%  Total_exp_mat_for_trials_magnitudes_before_organized= Total_exp_mat_for_trials_magnitudes_before(:,:,k_for_plot);
%  Total_exp_mat_for_trials_magnitudes_after_organized= Total_exp_mat_for_trials_magnitudes_after(:,:,k_for_plot);
%  Total_cont_mat_for_trials_magnitudes_before_organized= Total_cont_mat_for_trials_magnitudes_before(:,:,k_for_plot);
%  Total_cont_mat_for_trials_magnitudes_after_organized= Total_cont_mat_for_trials_magnitudes_after(:,:,k_for_plot);
%  
%  Total_exp_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
%  Total_exp_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
%  Total_cont_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
%  Total_cont_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
 
Total_exp_mat_for_trials_magnitudes_before_organized_abs=[];  %delete only this part to get them not as abs
Total_exp_mat_for_trials_magnitudes_after_organized_abs=[];
Total_cont_mat_for_trials_magnitudes_before_organized_abs=[];
Total_cont_mat_for_trials_magnitudes_after_organized_abs=[];
 
Total_exp_mat_for_trials_magnitudes_before_organized_abs=abs(Total_exp_mat_for_trials_magnitudes_before_organized);  %delete only this part to get them not as abs note that this is important for the KATO analysis
Total_exp_mat_for_trials_magnitudes_after_organized_abs=abs(Total_exp_mat_for_trials_magnitudes_after_organized);
Total_cont_mat_for_trials_magnitudes_before_organized_abs=abs(Total_cont_mat_for_trials_magnitudes_before_organized);
Total_cont_mat_for_trials_magnitudes_after_organized_abs=abs(Total_cont_mat_for_trials_magnitudes_after_organized);
 %%for exp not sorted by magnitude for each cell:
 
  x=1:11;
  y1=(mean(Total_exp_mat_for_trials_magnitudes_before_organized_abs,2));
      exp_abs_bef=squeeze(y1);

  std_y1=squeeze(std(y1,1));
  sem_y1=std_y1/sqrt(length(y1));
  y1=mean(y1,1);
  y1=squeeze(y1);

  y2=(mean(Total_exp_mat_for_trials_magnitudes_after_organized_abs,2));
    exp_abs_aft=squeeze(y2);
  std_y2=squeeze(std(y2,1));
  sem_y2=std_y1/sqrt(length(y2));
  y2=mean(y2,1);
  y2=squeeze(y2);
  
 mono_dif=y1(1:6)-y2(1:6);
 mono_dif_prop=mono_dif./y1(1:6);
 nat_dif=y1(7:11)-y2(7:11);
 nat_dif_prop=nat_dif./y1(7:11);
 
 [h,p]=ranksum(mono_dif_prop,nat_dif_prop);
 [h,p]=ttest2(mono_dif_prop,nat_dif_prop);
 
 figure
 errorbar(x,y1,sem_y1,'b','linewidth',2)
 hold on
 errorbar(x,y2,sem_y2,'r','linewidth',2)
 title('magnitude of response per odor-experiment')
 xlabel('Odor identity')
 ylabel('response magnitude')
  box off
 ylim([-0.05 0.4])
 ylim([0 0.2])  %for integral
 %%for cont not sorted by magnitude for each cell:
 
 x=1:11;
  y3=(mean(Total_cont_mat_for_trials_magnitudes_before_organized_abs,2));
 cont_abs_bef=squeeze(y3);

  std_y3=squeeze(std(y3,1));
  sem_y3=std_y3/sqrt(length(y3));
  y3=mean(y3,1);
  y3=squeeze(y3);

  y4=(mean(Total_cont_mat_for_trials_magnitudes_after_organized_abs,2));
  cont_abs_aft=squeeze(y4);
  std_y4=squeeze(std(y4,1));
  sem_y4=std_y4/sqrt(length(y4));
  y4=mean(y4,1);
  y4=squeeze(y4);
  
% stats for paper
%  [h,p]=ttest(cont_abs_bef, cont_abs_aft);
%  [h,p]=ttest(exp_abs_bef, exp_abs_aft);
 
 cont_for_stat_temp=cont_abs_bef-cont_abs_aft;
  cont_for_stat=mean(cont_for_stat_temp,1);
    
%   cont_for_stat=mean(cont_for_stat_temp,2);   uncomment these to check by
%   cells and not odors
 % [h,p,ci,stat]=ttest(cont_for_stat);
   [h,p,stat]=signrank(cont_for_stat)  %this is reported eventually
  
  exp_for_stat_temp=exp_abs_bef-exp_abs_aft;
  exp_for_stat=mean(exp_for_stat_temp,1);
 % [h,p,ci,stat]=ttest(exp_for_stat)
   [h,p,stat]=signrank(exp_for_stat)  %this is reported eventually
 
  figure
 errorbar(x,y3,sem_y3,'b','linewidth',2)
 hold on
 errorbar(x,y4,sem_y4,'r','linewidth',2)
 box off
 ylim([-0.05 0.4])
  ylim([0 0.2])  %for integral
 title('magnitude of response per odor-control')
  xlabel('Odor identity')
 ylabel('response magnitude')
 % now sorted by magnitude:
 
 %  %sort by magnitude:
 sort_exp_bef_abs= Total_exp_mat_for_trials_magnitudes_before_organized_abs;
 sort_exp_aft_abs=Total_exp_mat_for_trials_magnitudes_after_organized_abs;
 sort_cont_bef_abs=Total_cont_mat_for_trials_magnitudes_before_organized_abs;
 sort_cont_aft_abs=Total_cont_mat_for_trials_magnitudes_after_organized_abs;
 
 [sort_exp_bef_abs,Indices]=sort(squeeze(mean(sort_exp_bef_abs,2)),2);
  sort_exp_aft_abs=sort(squeeze(mean( sort_exp_aft_abs,2)),2);
 % sort_exp_aft_paired=sort_exp_aft_abs;  %change here to have abs for kato
  sort_exp_aft_paired=sort_exp_aft;
  sort_exp_aft_paired(Indices)= sort_exp_aft_paired;
  
  [sort_cont_bef_abs,Indices]=sort(squeeze(mean( sort_cont_bef_abs,2)),2);
  sort_cont_aft_abs=sort(squeeze(mean( sort_cont_aft_abs,2)),2);
% sort_cont_aft_paired=sort_cont_aft_abs;  %change here to have abs for kato
sort_cont_aft_paired=sort_cont_aft
  sort_cont_aft_paired(Indices)= sort_cont_aft_paired;
  
  
 mean_sort_exp_bef=mean(sort_exp_bef_abs,1);
 mean_sort_exp_aft=mean(sort_exp_aft_abs,1);
 mean_sort_cont_bef=mean(sort_cont_bef_abs,1);
 mean_sort_cont_aft=mean(sort_cont_aft_abs,1); 
 
 std_sort_exp_bef=std(sort_exp_bef_abs,1)/sqrt(length(sort_exp_bef_abs));
 std_sort_exp_aft=std(sort_exp_aft_abs,1)/sqrt(length(sort_exp_aft_abs));
 std_sort_cont_bef=std(sort_cont_bef_abs,1)/sqrt(length(sort_cont_bef_abs));
 std_sort_cont_aft=std(sort_cont_aft_abs,1)/sqrt(length(sort_cont_aft_abs)); 
 


%%
 
 
 
  %exp:
 figure
errorbar (x,mean_sort_exp_bef,std_sort_exp_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_exp_aft,std_sort_exp_aft,'r','linewidth',2)
title('responses sorted by abs magnitude-experiment')
 box off
 ylim([-0.2 0.2])
  ylim([0 0.3])
 small_dif=mean_sort_exp_bef(4:7)-mean_sort_exp_aft(4:7);
 small_dif_prop=small_dif./mean_sort_exp_bef(4:7);
 big_dif=mean_sort_exp_bef(8:11)-mean_sort_exp_aft(8:11);
 big_dif_prop=big_dif./mean_sort_exp_bef(8:11);
 
 [h,p]=ranksum(small_dif_prop,big_dif_prop);
 [h,p]=ttest2(small_dif_prop,big_dif_prop);
 
 %cont:
   figure
errorbar (x,mean_sort_cont_bef,std_sort_cont_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_cont_aft,std_sort_cont_aft,'r','linewidth',2)
 title('responses sorted by abs magnitude-control')
box off
 ylim([-0.2 0.2])
  ylim([0 0.3])
 

%paper stats

 %% measure the change un curves for the abs:
 
 %% 28.8.19 curves ratios analysis for absolute value!
std( mean_sort_exp_bef)/ std( mean_sort_exp_aft); %itamars measure for curvature- the whole population
std_ratio_per_cell_exp=std(sort_exp_bef_abs')./std(sort_exp_aft_abs');  %a vector of sds for each cell
mean_curves_ratio_exp=mean(std(sort_exp_bef_abs')./std(sort_exp_aft_abs')); %mean curve;
sem_curves_ratio_exp=(std(std_ratio_per_cell_exp))/sqrt(length(std_ratio_per_cell_exp)); %sem for curves


std( mean_sort_cont_bef)/ std( mean_sort_cont_aft); %itamars measure for curvature- the whole population
std_ratio_per_cell_cont=std(sort_cont_bef_abs')./std(sort_cont_aft_abs');  %a vector of sds for each cell
mean_curves_ratio_cont=mean(std(sort_cont_bef_abs')./std(sort_cont_aft_abs')); %mean curve;
sem_curves_ratio_cont=(std(std_ratio_per_cell_cont))/sqrt(length(std_ratio_per_cell_cont)); %sem for curves

% means=[mean_curves_ratio_exp mean_curves_ratio_cont]
% sems=[sem_curves_ratio_exp sem_curves_ratio_cont]
% figure
% errorbar(means,sems)
% [h,p]=ranksum(std_ratio_per_cell_exp,std_ratio_per_cell_cont) %%%the diffrence between curves- dont report this!
%  
% [h,p]=ttest(log(std_ratio_per_cell_cont))
% [h,p]=ttest(log(std_ratio_per_cell_exp))

% 
% [p,h,stats_cont]=signrank(log(std_ratio_per_cell_cont))
% [p,h,stats_exp]=signrank(log(std_ratio_per_cell_exp))
% 
% 
% 
% [h,p]=ttest(log(std_ratio_per_cell_cont))
% [h,p]=ttest(log(std_ratio_per_cell_exp))
% 
% [p,h,stats_cont]=signrank(log(std_ratio_per_cell_cont))
% [p,h,stats_exp]=signrank(log(std_ratio_per_cell_exp))

%stats for paper- curves in absolute values
cont_bef=std(sort_cont_bef_abs');
cont_aft=std(sort_cont_aft_abs');
[h,p,stat]=signrank(cont_bef,cont_aft)

exp_bef=std(sort_exp_bef_abs');
exp_aft=std(sort_exp_aft_abs');
[h,p,stat]=signrank(exp_bef,exp_aft)

 %% kato:
 resp_thresh=3; %include only cells that had this number of exc responses at least (both before and after)
 
 sort_exp_bef_norm=sort_exp_bef;
 sort_exp_aft_norm=sort_exp_aft_paired;
 sort_cont_bef_norm=sort_cont_bef;
 sort_cont_aft_norm=sort_cont_aft_paired;
 
 for i=1:size(sort_exp_bef,1);
   sort_exp_bef_norm(i,:)= sort_exp_bef_norm(i,:)./sort_exp_bef(i,end);
   sort_exp_aft_norm(i,:)= sort_exp_aft_norm(i,:)./sort_exp_bef(i,end);
 end
 
 for i=1:size(sort_cont_bef,1);
   sort_cont_bef_norm(i,:)= sort_cont_bef_norm(i,:)./sort_cont_bef(i,end);
   sort_cont_aft_norm(i,:)= sort_cont_aft_norm(i,:)./sort_cont_bef(i,end);
 end
 
meta_bef_af_exp=general_exp.all_total_meta_des_bef_af_inh_ex;
meta_bef_af_cont=general_cont.all_total_meta_des_bef_af_inh_ex;

sum_meta_exp=sum(meta_bef_af_exp~=0,2)-sum(meta_bef_af_exp<=-2,2);    %dont include inhibition
sum_meta_cont=sum(meta_bef_af_cont~=0,2)-sum(meta_bef_af_cont<=-2,2);
 
sort_exp_bef_norm(sum_meta_exp<resp_thresh,:)=[];
sort_exp_aft_norm(sum_meta_exp<resp_thresh,:)=[];

sort_cont_bef_norm(sum_meta_cont<resp_thresh,:)=[];
sort_cont_aft_norm(sum_meta_cont<resp_thresh,:)=[];

sort_exp_bef_norm=fliplr(sort_exp_bef_norm);
sort_exp_aft_norm=fliplr(sort_exp_aft_norm);

sort_cont_bef_norm=fliplr(sort_cont_bef_norm);
sort_cont_aft_norm=fliplr(sort_cont_aft_norm);


x=[6 4 2 1 3 5 7];
x2=[4 3 5 2 6 1 7];
x3=1:7;

y_bef=mean(sort_exp_bef_norm);
y_bef(8:end)=[];
y_bef=y_bef(x);

sem_bef=std(sort_exp_bef_norm)./sqrt(size(sort_exp_bef_norm,1));
sem_bef=sem_bef(x);

y_aft=mean(sort_exp_aft_norm);
y_aft(8:end)=[];
y_aft=y_aft(x);

sem_aft=std(sort_exp_aft_norm)./sqrt(size(sort_exp_aft_norm,1));
sem_aft=sem_aft(x);

%% Now linear regression for each cell sperately:
aaa=(-1:5);
bbb=(-1:5);
xx=-1:length(aaa);
x_ax=zeros(length(aaa)+2,1);

xxx=[0 0 0];
yyy=[-10 0 10];

b1_for_all_exp=[];
b_for_all_exp=[];

figure
for i=1:size(sort_exp_bef_norm,1);
    %figure      % cancel this to plot some of the cells

    x= sort_exp_bef_norm(i,1:7)';
    y= sort_exp_aft_norm(i,1:7)';
    
    b1 = x\y;
yCalc1 = b1*x;
%scatter(x,y,'k','linewidth',3)
hold on
%plot(x,yCalc1)
xlabel('Response magnitude before')
ylabel('Response magnitude after')

grid on
 xlim([-0.5 1.5]) 
 ylim([-0.5 1.5]) 
X = [ones(length(x),1) x];
b = X\y;

yCalc2 = X*b;
% plot(x,yCalc2,'k--','linewidth',3)             uncomment these to see single cells
% 
% plot(aaa,bbb,'k --')
% plot(xxx,yyy,'k -')
% plot(xx,x_ax,'k -')
%legend('Data','Slope','Slope & Intercept','X=Y','Location','best'); also
%for only slope
%legend('Data','Slope & Intercept','X=Y','Location','best'); %only slope+intercept

b1_for_all_exp=[b1_for_all_exp;b1];
b_for_all_exp(1,i)=b(1);
b_for_all_exp(2,i)=b(2);

title([' Intercept=' num2str(b(1)) ' Slope= ' num2str(b(2))])
end

%analyse intercept and slopes
figure
hold on
yy1=zeros(size( sort_exp_bef_norm,1),1);
yy2=ones(size( sort_exp_bef_norm,1),1);
% x1=ones(size( sort_exp_bef_norm,1),1);
% x2=ones(size( sort_exp_bef_norm,1),1)*2;
x1 = linspace(1,3,size( sort_exp_bef_norm,1));

x2 = linspace(6,8,size( sort_exp_bef_norm,1));

plot(x1,b_for_all_exp(1,:),'*')
plot(x2,b_for_all_exp(2,:),'*')
plot(x1,yy1,'k--');
plot(x2,yy2,'k--');

xlim([0 9])



h_exp=mean(b_for_all_exp,2);

%kato
figure
hold on
title(['Experiment-Kato analysis' num2str(size( sort_exp_bef_norm,1)) 'cells'])
ylim([0 1.2])
errorbar(x3,y_bef,sem_bef,'b','linewidth',2)
errorbar(x3,y_aft,sem_aft,'r','linewidth',2)
y_linear_prediction=h_exp(1)+y_bef.*h_exp(2);
plot(x3,y_linear_prediction,'k--','linewidth',2)


%%

%same for control
xx=-1:length(aaa);
x_ax=zeros(length(aaa)+2,1);

xxx=[0 0 0];
yyy=[-10 0 10];

x=[6 4 2 1 3 5 7];

y_bef=mean(sort_cont_bef_norm);
y_bef(8:end)=[];
y_bef=y_bef(x);

sem_bef=std(sort_cont_bef_norm)./sqrt(size(sort_cont_bef_norm,1));
sem_bef=sem_bef(x);

y_aft=mean(sort_cont_aft_norm);
y_aft(8:end)=[];
y_aft=y_aft(x);

sem_aft=std(sort_cont_aft_norm)./sqrt(size(sort_cont_aft_norm,1));
sem_aft=sem_aft(x);

aaa=(-1:5);
bbb=(-1:5);
b1_for_all_cont=[];
b_for_all_cont=[];

figure
for i=1:size( sort_cont_bef_norm,1);
%    figure      % cancel this to plot some of the cells

    x= sort_cont_bef_norm(i,1:7)';
    y= sort_cont_aft_norm(i,1:7)';
    
    b1 = x\y;
yCalc1 = b1*x;
%scatter(x,y,'k','linewidth',3)
hold on
%plot(x,yCalc1)
xlabel('Response magnitude before')
ylabel('Response magnitude after')

grid on
 xlim([-0.5 3.5]) 
 ylim([-0.5 3.5]) 
X = [ones(length(x),1) x];
b = X\y;

yCalc2 = X*b;
% plot(x,yCalc2,'k--','linewidth',3)   uncomment these to see single cells
% plot(xxx,yyy,'k -')
% plot(xx,x_ax,'k -')
% 
% plot(aaa,bbb,'k --')
%legend('Data','Slope','Slope & Intercept','X=Y','Location','best');
    
b1_for_all_cont=[b1_for_all_cont;b1];
b_for_all_cont(1,i)=b(1);
b_for_all_cont(2,i)=b(2);

title([' Intercept=' num2str(b(1)) ' Slope= ' num2str(b(2))])
end

%analyse intercept and slopes
figure
xlim([0 9])
ylim([-2 4])


hold on
x1 = linspace(1,3,size( sort_cont_bef_norm,1));
x2 = linspace(6,8,size( sort_cont_bef_norm,1));
plot(x1,b_for_all_cont(1,:),'o')
plot(x2,b_for_all_cont(2,:),'o')


yy1=zeros(size( sort_cont_bef_norm,1),1);
yy2=ones(size( sort_cont_bef_norm,1),1);
plot(x1,yy1,'k--');
plot(x2,yy2,'k--');
ylim([-2 4])





%kato
h_cont=mean(b_for_all_cont,2)

figure
title(['Control-Kato analysis ' num2str(size( sort_cont_bef_norm,1)) 'cells'])
hold on
ylim([0 1.2])
errorbar(x3,y_bef,sem_bef,'b','linewidth',2)
errorbar(x3,y_aft,sem_aft,'r','linewidth',2)
y_linear_prediction=h_cont(1)+y_bef.*h_cont(2);
plot(x3,y_linear_prediction,'k--','linewidth',2)
hold off

disp('no')

%% intercepts and slopes sep exp cont combined:
figure
ylim([-1 1])
hold on
y1=zeros(size( sort_exp_bef_norm,1),1);
y2=ones(size( sort_exp_bef_norm,1),1);
xx1 = linspace(1,2,size( sort_exp_bef_norm,1));
xx2 = linspace(7,8,size( sort_exp_bef_norm,1));

x1 = linspace(4,5,size( sort_cont_bef_norm,1));
x2 = linspace(10,11,size( sort_cont_bef_norm,1));
yy1=zeros(size( sort_cont_bef_norm,1),1);
yy2=ones(size( sort_cont_bef_norm,1),1);


plot(xx1,b_for_all_exp(1,:),'k*')
plot(x1,b_for_all_cont(1,:),'ko')

plot(xx1,y1,'r--');
plot(x1,yy1,'r--');
xlim([0 6])

figure
hold on
ylim([0 4])
xlim([6 12])

plot(xx2,b_for_all_exp(2,:),'k*')
plot(x2,b_for_all_cont(2,:),'ko')
plot(xx2,y2,'r--');
plot(x2,yy2,'r--');




mean(b_for_all_cont(1,:))
mean(b_for_all_exp(1,:))
mean(b_for_all_cont(2,:))
mean(b_for_all_exp(2,:))

[p,h]=ranksum(b_for_all_cont(1,:),b_for_all_exp(1,:))
[p,h]=ranksum(b_for_all_cont(2,:),b_for_all_exp(2,:))



%%
 sort_exp_bef_norm=sort_exp_bef;
 sort_exp_aft_norm=sort_exp_aft_paired;
 sort_cont_bef_norm=sort_cont_bef;
 sort_cont_aft_norm=sort_cont_aft_paired;
 

 %% Variance between trials
 STD_mat_exp_before=zeros(size(exp_matrix_before_organized,1),Num_of_odors);
 STD_mat_exp_after=zeros(size(exp_matrix_before_organized,1),Num_of_odors);
 STD_mat_cont_before=zeros(size(cont_matrix_before_organized,1),Num_of_odors);
 STD_mat_cont_after=zeros(size(cont_matrix_before_organized,1),Num_of_odors);
 
 for i=1:size(exp_matrix_before_organized,1);
     for j=1:Num_of_odors;
          STD_mat_exp_before(i,j)=std(Total_exp_mat_for_trials_magnitudes_before(i,:,j));
          STD_mat_exp_after(i,j)=std(Total_exp_mat_for_trials_magnitudes_after(i,:,j));
          
     end
 end
a= mean(STD_mat_exp_before,2);
b=mean(STD_mat_exp_after,2);
 
 STD_mat_cont_before=zeros(size(cont_matrix_before_organized,1),Num_of_odors);
 STD_mat_cont_after=zeros(size(cont_matrix_before_organized,1),Num_of_odors);
 
 for i=1:size(cont_matrix_before_organized,1);
     for j=1:Num_of_odors
          STD_mat_cont_before(i,j)=std(Total_cont_mat_for_trials_magnitudes_before(i,:,j));
          STD_mat_cont_after(i,j)=std(Total_cont_mat_for_trials_magnitudes_after(i,:,j));
          
     end
 end
c= mean(STD_mat_cont_before,2);
d=mean(STD_mat_cont_after,2);
%  
%% 
Total_signif_desicion_exp=[];

for j=1:Num_of_exp_mice;
% temp1=Exp_data.mice_table{1,j}.meta_desicion>1;
% temp2=Exp_data.mice_table{2,j}.meta_desicion>1;
 temp1=Exp_data.mice_table{1,j}.meta_desicion~=0;
 temp2=Exp_data.mice_table{2,j}.meta_desicion~=0;

temp3=temp1+temp2;
Total_signif_desicion_exp=[Total_signif_desicion_exp;temp3];
end


Total_signif_desicion_exp_organized=Total_signif_desicion_exp(:,k_for_plot);
Total_signif_desicion_exp_organized(:,1)=[];

Total_signif_desicion_cont=[];
for j=1:Num_of_cont_mice;
temp1=Cont_data.mice_table{1,j}.meta_desicion~=0;
temp2=Cont_data.mice_table{2,j}.meta_desicion~=0;
% temp1=Cont_data.mice_table{1,j}.meta_desicion>1;
% temp2=Cont_data.mice_table{2,j}.meta_desicion>1;

temp3=temp1+temp2;
Total_signif_desicion_cont=[Total_signif_desicion_cont;temp3];
end


Total_signif_desicion_cont_organized=Total_signif_desicion_cont(:,k_for_plot);
Total_signif_desicion_cont_organized(:,1)=[];
%%
a=STD_mat_exp_before(Total_signif_desicion_exp_organized>0);
b=STD_mat_exp_after(Total_signif_desicion_exp_organized>0);

c=STD_mat_cont_before(Total_signif_desicion_cont_organized>0);
d=STD_mat_cont_after(Total_signif_desicion_cont_organized>0);

[h,p]=ttest(a,b);
[h,p]=ttest(c,d);
j=a-b;
k=c-d;
[h,p]=ttest2(j,k);
x=[1 2];
y=[mean(a) mean(b)];
z=[mean(c) mean(d)];

e1=[std(a)/sqrt(length(a)) std(b)/sqrt(length(b))]
e2=[std(c)/sqrt(length(c)) std(d)/sqrt(length(d))]

figure
% plot(x,y);
% hold on
% plot(x,z);
hold on 
errorbar(x,y, e1)
 errorbar(x,z, e2)

title(['Variability between trials, control vs experiment. P=' num2str(p)])
xlabel('Before                                               After','FontSize',14)
ylabel('Variability between trials (SD)','FontSize',14)
legend('Exp','Cont', 'Location', 'northwest')

%%
 %% chi_test_for_dots_plot- proportion of cells
cont_observ=general_cont.stats_for_chi;
exp_observ=general_exp.stats_for_chi;

  %first all then all non signif then those that are strongerbefore then stronger after
       cont_non_sig = cont_observ(2); cont_bef=cont_observ(3);cont_af=cont_observ(4); all_resp_cont = cont_observ(1);   %even vs odd
       exp_non_sig = exp_observ(2); exp_bef=exp_observ(3);exp_af=exp_observ(4); all_resp_exp = exp_observ(1); 

        % Pooled estimate of proportion
       p_non_sig= (cont_non_sig+exp_non_sig)/(all_resp_cont+all_resp_exp);
       p_bef= (cont_bef+exp_bef)/(all_resp_cont+all_resp_exp);
       p_af=(exp_af+cont_af)/(all_resp_cont+all_resp_exp);
       
       % Expected counts under H0 (null hypothesis)
       p_non_sig_cont_expected= p_non_sig*all_resp_cont;
       p_bef_cont_expected=p_bef*all_resp_cont;
       p_af_cont_expected=p_af*all_resp_cont;
       
       p_non_exp_exp_expected= p_non_sig*all_resp_exp;
       p_bef_exp_expected=p_bef*all_resp_exp;
       p_af_exp_expected=p_af*all_resp_exp;
        
       a_c=p_non_sig_cont_expected;  %just for comfort
       b_c= p_bef_cont_expected;
       c_c=p_af_cont_expected;
       
       a_e=p_non_exp_exp_expected;
       b_e= p_bef_exp_expected;
       c_e= p_af_exp_expected;
       
       expected=[a_c b_c c_c a_e b_e c_e]
       observed=[cont_observ(2:end) exp_observ(2:end)]
    
       % Chi-square test, by hand
      
       chi2stat = sum((observed-expected).^2 ./ expected)
       p = 1 - chi2cdf(chi2stat,2) %degrees of freedom is (rows-1)*(colmuns-1)
       
       %normalization
       cont_observ_norm=cont_observ/cont_observ(1)
       exp_observ_norm=exp_observ/exp_observ(1)
       
       %plot
       x=1:6;
       figure
       
     
       bar(x(1:2:5),cont_observ(2:end),0.45)
       hold on
       bar(x(2:2:6),exp_observ(2:end),0.45,'k')
      
       ax=gca;
       ax.FontSize=14;
       ax.XTick=[];
      
      
      
      title(['Changes in response magnitude. P=' num2str(p)])
      xlabel('Non significant      Stronger Before        Stronger after','FontSize',14)
      ylabel('Number of responses','FontSize',14)
      legend('Control','Experiment', 'Location', 'northeast')


       %plot normalized bar
       x=1:6;
       figure
       
     
       bar(x(1:2:5),cont_observ_norm(2:end),0.45)
       hold on
       bar(x(2:2:6),exp_observ_norm(2:end),0.45,'k')
      
       ax=gca;
       ax.FontSize=14;
       ax.XTick=[];
      
      
      
      title(['Changes in response magnitude. P=' num2str(p)])
      xlabel('Non significant      Stronger Before        Stronger after','FontSize',14)
      ylabel('Number of responses','FontSize',14)
      legend('Control','Experiment', 'Location', 'northeast')
%% mice deltas

deltas_cont=general_cont.mice_deltas;
deltas_exp=general_exp.mice_deltas;



figure
for i=1:length(deltas_cont)
    plot(2, deltas_cont(i),'ko','linewidth', 2)
    hold on
end


for i=1:length(deltas_exp)
    plot(1, deltas_exp(i),'ro','linewidth', 2)
    hold on
end
y=zeros(1,4);
x=(0:3);
plot(x,y,'k--')

xlim([0 3])
ylim([-0.2 0.2])


mean_cont=mean(deltas_cont);
mean_exp=mean(deltas_exp);

means=[mean_exp,mean_cont];

sem_cont=std(deltas_cont)/sqrt(length(deltas_cont));
sem_exp=std(deltas_exp)/sqrt(length(deltas_exp));

sems=[sem_exp,sem_cont];

errorbar(means, sems,'.','linewidth', 2);

xx=1:2;
if length(deltas_cont)==length(deltas_exp)
for i=1:length(deltas_cont)
    yy=[deltas_exp(i) deltas_cont(i)];
    plot(xx,yy,'k')
    hold on
end
end

[h,p]=ttest(deltas_exp,deltas_cont,'Tail','right')  %for paired- one tail assumption. use for awake repeated measures
[h,p,stats]=signrank(deltas_exp,deltas_cont,'Tail','right')

[h,p]=ttest2(deltas_exp,deltas_cont,'Tail','right')  %for unpaired. use for different subjects as exp-cont
[h,p,stats]=ranksum(deltas_exp,deltas_cont,'Tail','right') 

title('abs deltas per mice- dff')
%%
%% prop (ratios)

deltas_cont_prop=general_cont.prop_for_comparison;
deltas_exp_prop=general_exp.prop_for_comparison;


figure
for i=1:length(deltas_cont_prop);
    plot(1, deltas_exp_prop(i),'ko','linewidth', 2)
    hold on
end


for i=1:length(deltas_exp_prop);
    plot(2, deltas_cont_prop(i),'ro','linewidth', 2);
    hold on
end

y=ones(1,4);
x=(0:3);

hold on
plot(x,y,'k--','linewidth', 2)

% xlim([0 3])
 ylim([0.1 10])

ax = gca;
ax.YScale = 'log'
ax.FontSize = 14;
box off;
set(gca,'xtick',[]);

mean_cont=mean(deltas_cont_prop);
mean_exp=mean(deltas_exp_prop);

means=[mean_exp,mean_cont];

sem_cont=std(deltas_cont_prop)/sqrt(length(deltas_cont_prop));
sem_exp=std(deltas_exp_prop)/sqrt(length(deltas_exp_prop));

sems=[sem_exp,sem_cont];

errorbar(means, sems,'.','linewidth', 2);


xx=1:2;
if length(deltas_cont_prop)==length(deltas_exp_prop)
for i=1:length(deltas_cont_prop)
    yy=[deltas_exp_prop(i) deltas_cont_prop(i)];
    plot(xx,yy,'k')
    ax = gca;
ax.YScale = 'log'
    hold on
end
end

[h,p]=ttest(deltas_exp_prop,deltas_cont_prop)  %for paired
[h,p]=signrank(deltas_exp_prop,deltas_cont_prop)

[h,p]=ttest2(deltas_exp_prop,deltas_cont_prop)  %for unpaired
[h,p]=ranksum(deltas_exp_prop,deltas_cont_prop) 
%ax.xtick=[];

% XTick off
% YTick off

%%
deltas_cont_prop=general_cont.prop_for_comparison;
deltas_exp_prop=general_exp.prop_for_comparison;


figure
for i=1:length(deltas_cont_prop)
    plot(1, deltas_cont_prop(i),'ko')
    hold on
end


for i=1:length(deltas_exp_prop);
    plot(2, deltas_exp_prop(i),'ro')
    hold on
end

y=ones(1,4);
x=(0:3);

hold on
plot(x,y,'k--','linewidth', 1.5)

% xlim([0 3])
 ylim([0.1 10])

ax = gca;
ax.YScale = 'log'
ax.FontSize = 14;
box off
set(gca,'xtick',[])

mean_cont=mean(deltas_cont_prop);
mean_exp=mean(deltas_exp_prop);

means=[mean_cont,mean_exp];

sem_cont=std(deltas_cont_prop)/length(deltas_cont_prop);
sem_exp=std(deltas_exp_prop)/length(deltas_exp_prop);

sems=[sem_cont,sem_exp];

errorbar(means, sems,'.');


title('abs deltas per mice- dff')

title('ratios per mice- dff')
%ax.xtick=[];

% XTick off
% YTick off
%%

% exp_deltas_bp=general_exp.deltas_for_boxplots;
% cont_deltas_bp=general_cont.deltas_for_boxplots;

exp_deltas_bp=general_exp.abs_deltas_for_boxplots;
cont_deltas_bp=general_cont.abs_deltas_for_boxplots;


%% boxplots
kx=[0:4];
ky=[0 0 0 0 0];

figure
%x=[all_significant_responses_integrals_before;all_significant_responses_integrals_after];
x=[cont_deltas_bp;exp_deltas_bp];

g = [ones(size(cont_deltas_bp)); 2*ones(size(exp_deltas_bp))];
% boxplot(x,g)
boxplot(x,g,'symbol','');

%title(['Before and After ' num2str(length(all_significant_responses_magnitudes_before)) ' significant responses ' ' Before Mean=' num2str(mean_before) ' SD=' num2str(SD_before) ' After Mean=' num2str(mean_after) ' SD=' num2str(SD_after) '    p=' num2str(p)],'FontSize',24);% '    p2=' num2str(p_for_all)]);
%ylabel('Response magnitude')
%xlabel('Control                                                       Experiment')
%ylim([0 0.6])
ylim([-0.3 0.3])
a=gca;
a.FontSize=20;
hold on
plot(kx,ky,'k--','linewidth',2)

[h,p]=ttest2(cont_deltas_bp,exp_deltas_bp)
[p,h]=ranksum(cont_deltas_bp,exp_deltas_bp)
% figure
% boxplot(all_significant_responses_magnitudes_before)
% title('Before')
% ylim([0 1.8])
% 
% figure
% boxplot(all_significant_responses_magnitudes_after)
% title('after')
% ylim([0 1.8])
box off

%% bars for comparison amount of responses
resp_per_cell_exp_before=general_exp.total_prop_of_cells_before;
resp_per_cell_exp_after=general_exp.total_prop_of_cells_after;

resp_per_cell_cont_before=general_cont.total_prop_of_cells_before;
resp_per_cell_cont_after=general_cont.total_prop_of_cells_after;

exp=resp_per_cell_exp_after-resp_per_cell_exp_before;
cont=resp_per_cell_cont_after-resp_per_cell_cont_before;
%stats used for the paper
[h,p,ci,stat]=ttest2(exp,cont)    %for unpiared data
if size(exp)==size(cont)
[h,p,ci,stat]=ttest(exp,cont)    %for paired data
end

 %used only for CNO control which is not compared to anything. report in
 %the paper!
 [h,p,ci,stat]=ttest(cont) 
 [h,p,stat]=signrank(cont) 


exp_norm_general=mean(resp_per_cell_exp_before);
cont_norm_general=mean(resp_per_cell_cont_before);

y=[mean(exp),mean(cont)];
sem_exp=std(exp)/sqrt(size(exp,1));
sem_cont=std(cont)/sqrt(size(cont,1));
sems=[sem_exp sem_cont]

figure
hold on
x=(1:2)

bar(x(1),y(1),'k')
bar(x(2),y(2),'k')
errorbar(x,y,sems,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]
ylim([-2.5 0.5])
title('Change in responsiveness # of responses per cell')

y=[mean(exp)/exp_norm_general,mean(cont)/cont_norm_general];
sem_exp=std(exp)/sqrt(size(exp,1));
sem_cont=std(cont)/sqrt(size(cont,1));
sems_norm=[sem_exp/exp_norm_general sem_cont/cont_norm_general];

figure
hold on
x=(1:2)

bar(x(1),y(1),'k')
bar(x(2),y(2),'k')
errorbar(x,y,sems_norm,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

title('Proportional Change in responsiveness-general Norm')

exp_prop=(resp_per_cell_exp_after-resp_per_cell_exp_before)./(resp_per_cell_exp_before+resp_per_cell_exp_after);
cont_prop=(resp_per_cell_cont_after-resp_per_cell_cont_before)./(resp_per_cell_cont_before+resp_per_cell_cont_after);
%cont_prop(cont_prop==inf)=1;
cont_prop(isnan(cont_prop))=0;
%exp_prop(exp_prop==inf)=1;
exp_prop(isnan(exp_prop))=0;

y=[mean(exp_prop) mean(cont_prop)];
sem_exp_prop=std(exp_prop)/sqrt(size(exp_prop,1));
sem_cont_prop=std(cont_prop)/sqrt(size(cont_prop,1));
sems_prop=[sem_exp_prop sem_cont_prop];

figure
hold on
x=(1:2)

bar(x(1),y(1),'k')
bar(x(2),y(2),'k')
errorbar(x,y,sems_prop,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

title('Proportional Change in responsiveness-individual Norm')

%% %% bars for comparison magnitude of responses per responses
resp_over_time_exp_ex=general_exp.deltas_over_time_ex;
resp_over_time_exp_in=general_exp.deltas_over_time_in;

resp_over_time_cont_ex=general_cont.deltas_over_time_ex;
resp_over_time_cont_in=general_cont.deltas_over_time_in;


exp_ex=mean(resp_over_time_exp_ex(:,49:77),2);
cont_ex=mean(resp_over_time_cont_ex(:,49:77),2);

length(exp_ex)
length(cont_ex)



exp_in=mean(resp_over_time_exp_in(:,49:77),2);
cont_in=mean(resp_over_time_cont_in(:,49:77),2);

length(exp_in)
length(cont_in)
%Stats for paper- assesing the difference for excitatory responses
[p,h,ci,stat]=ttest2(exp_ex,cont_ex)   %unpaired
[p,h,ci,stat]=ttest2(exp_in,cont_in) 

%only for cno control- compare to zero:
[p,h,ci,stat]=ttest(cont_ex)
[p,h,ci,stat]=ttest(cont_in)

if size(exp_ex)==size(cont_ex)
[p,h,ci,stat]=ttest(exp_ex,cont_ex) 
end


y_ex=[mean(exp_ex),mean(cont_ex)].*(-1);     %to make after-beofre (flip) 
sem_exp_ex=std(exp_ex)/sqrt(size(exp_ex,1));
sem_cont_ex=std(cont_ex)/sqrt(size(cont_ex,1));
sems_ex=[sem_exp_ex sem_cont_ex];

y_in=[mean(exp_in),mean(cont_in)].*(-1);              %to make after-beofre (flip) 
sem_exp_in=std(exp_in)/sqrt(size(exp_in,1));
sem_cont_in=std(cont_in)/sqrt(size(cont_in,1));
sems_in=[sem_exp_in sem_cont_in];


figure
hold on
x=(1:2);

bar(x(1),y_ex(1),'k')
bar(x(2),y_ex(2),'k')
errorbar(x,y_ex,sems_ex,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

title('Response magnitudes Deltas Excitatory')
ylim([min([y_in y_ex])-max([sems_ex sems_in]) max([y_in y_ex])+max([sems_ex sems_in])])
ylim([-0.05 0.08])


figure
hold on
x=(1:2);

bar(x(1),y_in(1),'k')
bar(x(2),y_in(2),'k')
errorbar(x,y_in,sems_in,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

ylim([min([y_in y_ex])-max([sems_ex sems_in]) max([y_in y_ex])+max([sems_ex sems_in])])
ylim([-0.05 0.08])
title('Response magnitudes Deltas Inhibitory')

%% statistics with tables. first for number of responses. this is for the anesthesized data, where bef-aft are not paired
if data_set(2)=='n'
A=zeros(1,size(exp,1)+size(cont,1));
A(1:size(exp,1))=1;
Conditions = [categorical(A,[1 0],{'exp' 'cont'})]'

all_n_resp_bef=[resp_per_cell_exp_before;resp_per_cell_cont_before];
all_n_resp_aft=[resp_per_cell_exp_after;resp_per_cell_cont_after];
deltas_n_resp=all_n_resp_bef-all_n_resp_aft

t_table_n_responses = table(Conditions,all_n_resp_bef,all_n_resp_aft,...
'VariableNames',{'Conditions','before','after'});
Times = table([1 2]','VariableNames',{'Measurements'});
Times.Measurements=categorical(Times.Measurements);
rm = fitrm(t_table_n_responses,'before-after~ Conditions','WithinDesign',Times);
 rm.Coefficients
rm.DFE
ranovatbl = ranova(rm,'WithinModel','Measurements');   %full model with the before after
%ranovatbl = ranova(rm,'WithinModel','Measurements')
stats.n_responses.general_results_full_model=ranovatbl; 
tbl_main = multcompare(rm,'Conditions'); %this shows no significant diff between cont and exp
stats.n_resp.multicompare_tukey_kramer_main=tbl_main;   % this is just the general effect accross
tbl_spec = multcompare(rm,'Measurements','By','Conditions');  %this is tukey kremer, specific effects for manipulation by group
stats.n_resp.multicompare_tukey_kramer_spec=tbl_spec;
[p,tbl]=anova1(deltas_n_resp,Conditions);   %one way anova on the deltas- same like the interaction term
stats.n_responses_deltas=tbl;    %this is the 1 way comparison for the deltas
disp('These were results for number of responses for anes') %meaning:  (Intercept):Measurements is for the effect of bef-af    Conditions:Measurements  is for the interaction

%% %%now for magnitudes, all for cells: (first deltas)

sss=zeros(22,1); sss(12:end)=1;
Times = [categorical(sss,[0 1],{'before' 'after'})];

aaaa=sort_cont_bef-sort_cont_aft;
bbbb=sort_exp_bef-sort_exp_aft;
% dummy=NaN(abs(size(aaaa,1)-size(bbbb,1)),11)
% bbbb=[bbbb;dummy]
cccc=[aaaa;bbbb]
[d,p,stats] = manova1(cccc,Conditions)
%% now the full data set


abs(size(aaaa,1)-size(bbbb,1))

sort_bef_exp_cont=[sort_exp_bef;sort_cont_bef];
sort_aft_exp_cont=[sort_exp_aft;sort_cont_aft];
deltas_sort=sort_bef_exp_cont-sort_aft_exp_cont;


t_table_sorted_mag= table();
t_table_sorted_mag{:,1} = Conditions;
t_table_sorted_mag{:,2:12} = sort_bef_exp_cont;
t_table_sorted_mag{:,13:23} = sort_aft_exp_cont;
t_table_sorted_mag.Properties.VariableNames={'Conditions','odor1b','odor2b','odor3b','odor4b','odor5b','odor6b','odor7b','odor8b','odor9b','odor10b','odor11b','odor1a','odor2a','odor3a','odor4a','odor5a','odor6a','odor7a','odor8a','odor9a','odor10a','odor11a'};

within= table([1:11 1:11]',Times,'VariableNames',{'Ranks','Times'});
%Times.Measurements=categorical(Times.Measurements)     %rank has meaning
%here, so don't change
rm_mag = fitrm(t_table_sorted_mag,'odor1b-odor11a~ Conditions','WithinDesign',within,'WithinModel','Ranks+Times');
% manovatbl = manova(rm_mag,'WithinModel','Times')
manova(rm_mag,'By','Conditions','WithinModel','Times')   %here mannova is used because there are so many means. this is the right comp
tbl_spec = multcompare(rm_mag,'Times','By','Conditions')


% ranovatbl = ranova(rm_mag,'WithinModel','Times')


% tbl_main = multcompare(rm_mag,'Conditions'); %this shows  significant diff between cont and exp
% tbl_main = multcompare(rm_mag,'Times'); %this shows  significant diff between cont and exp
% tbl_main = multcompare(rm_mag,'Ranks')

% stats.sorted_mag_resp.multicompare_tukey_kramer_main=tbl_main;   % this is just the general effect accross

% tbl_spec = multcompare(rm_mag,'Times','By','Conditions');  %this is tukey kremer, specific effects for manipulation by group
% tbl_spec = multcompare(rm_mag,'Conditions','By','Times')

% stats.sorted_mag_resp.multicompare_tukey_kramer_spec=tbl_spec;




% [p,tbl]=anovan(deltas_n_resp,Conditions);   %one way anova on the deltas- same like the interaction term
% stats.n_responses_deltas=tbl;    %this is the 1 way comparison for the deltas

%NOW MODEL FOR THE DELTAS::

% 
% t_table_sorted_mag_deltas= table();
% t_table_sorted_mag_deltas{:,1} = Conditions;
% t_table_sorted_mag_deltas{:,2:12} = deltas_sort;
% t_table_sorted_mag_deltas.Properties.VariableNames={'Conditions','odor1','odor2','odor3','odor4','odor5','odor6','odor7','odor8','odor9','odor10','odor11'};
% 
% within= table([1:11]','VariableNames',{'Ranks'});
% %Times.Measurements=categorical(Times.Measurements)     %rank has meaning
% %here, so don't change
% rm__mag_deltas = fitrm(t_table_sorted_mag_deltas,'odor1-odor11~ Conditions','WithinDesign',within);
% 
% 
% %Times.Measurements=categorical(Times.Measurements)     %rank has meaning
% %here, so don't change
% manovatbl = manova(rm__mag_deltas,'WithinModel','Ranks')
% ranovatbl = manova(rm__mag_deltas)
% 
% 
% stats.mag_sorted_deltas.general_results=ranovatbl;   %measurenents are the ranks, cond foe exp-cont

% figure
% plot(rm__mag_deltas)
% tbl = multcompare(rm__mag_deltas,'Conditions')
% 
% 
% % tbl = multcompare(rm_deltas,'Measurements')
% 
% stats.mag_sorted_deltas.main=tbl;   %measurenents are the ranks, cond foe exp-cont

% tbl = multcompare(rm_deltas,'Conditions','By','Measurements')
%% at the response level: 

end


%% read the end of this cell!!!

% Feilds_per_mouse_cont=zeros(1,size(cont,1))
% Feilds_per_mouse_exp=zeros(1,size(cont,1))
% last_cont=1;
% last_exp=1;
% n_mouse=0;
% mouse_name_exp={};
% mouse_name_cont={};
% 
% for i=1:Num_of_cont_mice;
% n_mouse=n_mouse+1;
% cells_per_mouse_exp=sum(Exp_data.mice_table{1,1}.fields_for_analysis(i).cells_per_field);
% cells_per_mouse_cont=sum(Cont_data.mice_table{1,1}.fields_for_analysis(i).cells_per_field);
% mouse_name_exp{i}=Exp_data.mice_table{1,1}.fields_for_analysis(i).name
% mouse_name_cont{i}=Cont_data.mice_table{1,1}.fields_for_analysis(i).name
% Feilds_per_mouse_exp(last_exp:last_exp+cells_per_mouse_exp-1)=n_mouse;
% Feilds_per_mouse_cont(last_cont:last_cont+cells_per_mouse_cont-1)=n_mouse;
% last_exp=last_exp+cells_per_mouse_exp;
% last_cont=last_cont+cells_per_mouse_cont;
% end
% 
% Feilds_per_mouse_exp(Feilds_per_mouse_exp==0)=[]
% Feilds_per_mouse_cont(Feilds_per_mouse_cont==0)=[]
% 
% Mice_exp = [categorical(Feilds_per_mouse_exp,[1:10],mouse_name_exp)]'
% Mice_cont = [categorical(Feilds_per_mouse_cont,[1:10],mouse_name_cont)]'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%commented area is good and ready to add mice identity!!!!!!!!!1
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if data_set(2)=='n'
A=zeros(1,size(exp,1)+size(cont,1));
A(1:size(exp,1))=1;
Conditions = [categorical(A,[1 0],{'exp' 'cont'})]'

all_n_resp_bef=[resp_per_cell_exp_before;resp_per_cell_cont_before];
all_n_resp_aft=[resp_per_cell_exp_after;resp_per_cell_cont_after];
deltas_n_resp=all_n_resp_bef-all_n_resp_aft

t_table_n_responses = table(Conditions,all_n_resp_bef,all_n_resp_aft,...
'VariableNames',{'Conditions','before','after'});
Times = table([1 2]','VariableNames',{'Measurements'});
Times.Measurements=categorical(Times.Measurements);
rm = fitrm(t_table_n_responses,'before-after~ Conditions','WithinDesign',Times);
% rm.Coefficients
rm.DFE
ranovatbl = ranova(rm);   %full model with the before after
%ranovatbl = ranova(rm,'WithinModel','Measurements')
stats.n_responses.general_results_full_model=ranovatbl; 
tbl_main = multcompare(rm,'Conditions'); %this shows no significant diff between cont and exp
stats.n_resp.multicompare_tukey_kramer_main=tbl_main;   % this is just the general effect accross
tbl_spec = multcompare(rm,'Measurements','By','Conditions');  %this is tukey kremer, specific effects for manipulation by group
stats.n_resp.multicompare_tukey_kramer_spec=tbl_spec;
[p,tbl]=anova1(deltas_n_resp,Conditions);   %one way anova on the deltas- same like the interaction term
stats.n_responses_deltas=tbl;    %this is the 1 way comparison for the deltas
disp('These were results for number of responses for anes') %meaning:  (Intercept):Measurements is for the effect of bef-af    Conditions:Measurements  is for the interaction

%% %%now for magnitudes, all for cells: (first deltas)

sss=zeros(22,1); sss(12:end)=1;
Times = [categorical(sss,[0 1],{'before' 'after'})];

aaaa=sort_cont_bef-sort_cont_aft;
bbbb=sort_exp_bef-sort_exp_aft;
% dummy=NaN(abs(size(aaaa,1)-size(bbbb,1)),11)
% bbbb=[bbbb;dummy]
cccc=[aaaa;bbbb]
[d,p,stats] = manova1(cccc,Conditions)
%% now the full data set


abs(size(aaaa,1)-size(bbbb,1))

sort_bef_exp_cont=[sort_exp_bef;sort_cont_bef];
sort_aft_exp_cont=[sort_exp_aft;sort_cont_aft];
deltas_sort=sort_bef_exp_cont-sort_aft_exp_cont;


t_table_sorted_mag= table();
t_table_sorted_mag{:,1} = Conditions;
t_table_sorted_mag{:,2:12} = sort_bef_exp_cont;
t_table_sorted_mag{:,13:23} = sort_aft_exp_cont;
t_table_sorted_mag.Properties.VariableNames={'Conditions','odor1b','odor2b','odor3b','odor4b','odor5b','odor6b','odor7b','odor8b','odor9b','odor10b','odor11b','odor1a','odor2a','odor3a','odor4a','odor5a','odor6a','odor7a','odor8a','odor9a','odor10a','odor11a'};

within= table([1:11 1:11]',Times,'VariableNames',{'Ranks','Times'});
%Times.Measurements=categorical(Times.Measurements)     %rank has meaning
%here, so don't change
rm_mag = fitrm(t_table_sorted_mag,'odor1b-odor11a~ Conditions','WithinDesign',within,'WithinModel','Ranks+Times');
% manovatbl = manova(rm_mag,'WithinModel','Times')
manova(rm_mag,'By','Conditions','WithinModel','Times')   %here mannova is used because there are so many means. this is the right comp
tbl_spec = multcompare(rm_mag,'Times','By','Conditions')


%

end