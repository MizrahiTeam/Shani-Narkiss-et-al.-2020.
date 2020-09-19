%%IMPORTANT:::    STARTED TO WRITE THIS CODE, BUT IT IS NOT WORKING. NOT
%%NECCESAY, NOT IN USE

%Comaprisons_thirdrol_vs_firsteriment
clc
clear all
close all
%%
% %for time lapse Fig 5

first_data=load('Mice_table_first_timelapse.mat');
sec_data=load('Mice_table_sec_timelapse.mat');
third_data=load('Mice_table_third_timelapse.mat');

general_first=load('General_results_first_timelapse.mat');
general_sec=load('General_results_sec_timelapse.mat');
general_third=load('General_results_third_timelapse.mat');

Num_of_odors=size(first_data.mice_table{1,1}.meta_desicion,2);

temp_for_first=size(first_data.mice_table);
Num_of_first_mice=temp_for_first(2);
Num_of_first_conditions=temp_for_first(1);

temp_for_sec=size(sec_data.mice_table);
Num_of_sec_mice=temp_for_sec(2);
Num_of_sec_conditions=temp_for_sec(1);


temp_for_third=size(third_data.mice_table);
Num_of_third_mice=temp_for_third(2);
Num_of_third_conditions=temp_for_third(1);

 valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;    %%new orfer for valves
    ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
    
    for gg=1:length( ploting)
    k_for_plot(gg)=find( valves_12==ploting(gg));
    end




%% all trials magnitudes together histograms for magnitude- sorted and not sorted by oddor strength
Total_first_mat_for_trials_magnitudes_before=[];
Total_first_mat_for_trials_magnitudes_after=[];

Total_sec_mat_for_trials_magnitudes_before=[];
Total_sec_mat_for_trials_magnitudes_after=[];

Total_third_mat_for_trials_magnitudes_before=[];
Total_third_mat_for_trials_magnitudes_after=[];

for j=1:Num_of_first_mice
 Total_first_mat_for_trials_magnitudes_before=[Total_first_mat_for_trials_magnitudes_before;first_data.mice_table{1,j}.mat_for_diffrence_test] ; 
 Total_first_mat_for_trials_magnitudes_after=[Total_first_mat_for_trials_magnitudes_after;first_data.mice_table{2,j}.mat_for_diffrence_test];  
end

for j=1:Num_of_sec_mice
 Total_sec_mat_for_trials_magnitudes_before=[Total_sec_mat_for_trials_magnitudes_before;sec_data.mice_table{1,j}.mat_for_diffrence_test] ; 
 Total_sec_mat_for_trials_magnitudes_after=[Total_sec_mat_for_trials_magnitudes_after;sec_data.mice_table{2,j}.mat_for_diffrence_test];  
end

for j=1:Num_of_third_mice
 Total_third_mat_for_trials_magnitudes_before=[Total_third_mat_for_trials_magnitudes_before;third_data.mice_table{1,j}.mat_for_diffrence_test];  
 Total_third_mat_for_trials_magnitudes_after=[Total_third_mat_for_trials_magnitudes_after;third_data.mice_table{2,j}.mat_for_diffrence_test];  
end

%added here down : for each not to reduce dimension

 Total_first_mat_for_trials_magnitudes_before_organized= Total_first_mat_for_trials_magnitudes_before(:,:,k_for_plot);
 Total_first_mat_for_trials_magnitudes_after_organized= Total_first_mat_for_trials_magnitudes_after(:,:,k_for_plot);
 Total_sec_mat_for_trials_magnitudes_before_organized= Total_sec_mat_for_trials_magnitudes_before(:,:,k_for_plot);
 Total_sec_mat_for_trials_magnitudes_after_organized= Total_sec_mat_for_trials_magnitudes_after(:,:,k_for_plot);
 Total_third_mat_for_trials_magnitudes_before_organized= Total_third_mat_for_trials_magnitudes_before(:,:,k_for_plot);
 Total_third_mat_for_trials_magnitudes_after_organized= Total_third_mat_for_trials_magnitudes_after(:,:,k_for_plot);
 
 Total_first_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
 Total_first_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
 Total_sec_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
 Total_sec_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
 Total_third_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
 Total_third_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
 

 
  x=1:11;
  y1=(mean(Total_first_mat_for_trials_magnitudes_before_organized,2));
  std_y1=squeeze(std(y1,1));
  sem_y1=std_y1/sqrt(length(y1));
  y1=mean(y1,1);
  y1=squeeze(y1);

  y2=(mean(Total_first_mat_for_trials_magnitudes_after_organized,2))
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
 title('magnitude of response per odor-firsteriment')
 xlabel('Odor identity')
 ylabel('response magnitude')
  box off
  
  ylim([-0.1 0.1])
% ylim([-0.05 0.4])
 %ylim([0 0.2])  %for integral
 %%for third not sorted by magnitude for each cell:
 
 x=1:11;
  y3=(mean(Total_third_mat_for_trials_magnitudes_before_organized,2));
  std_y3=squeeze(std(y3,1));
  sem_y3=std_y3/sqrt(length(y3));
  y3=mean(y3,1);
  y3=squeeze(y3);

  y4=(mean(Total_third_mat_for_trials_magnitudes_after_organized,2));
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
 title('magnitude of response per odor-thirdrol')
  xlabel('Odor identity')
 ylabel('response magnitude')
 % now sorted by magnitude:
 
 %  %sort by magnitude:
 sort_first_bef= Total_first_mat_for_trials_magnitudes_before_organized;
 sort_first_aft=Total_first_mat_for_trials_magnitudes_after_organized;
 sort_sec_bef= Total_sec_mat_for_trials_magnitudes_before_organized;
 sort_sec_aft=Total_sec_mat_for_trials_magnitudes_after_organized;
 sort_third_bef=Total_third_mat_for_trials_magnitudes_before_organized;
 sort_third_aft=Total_third_mat_for_trials_magnitudes_after_organized;
 
 [sort_first_bef,Indices]=sort(squeeze(mean(sort_first_bef,2)),2);
  sort_first_aft=sort(squeeze(mean( sort_first_aft,2)),2);
  sort_first_aft_paired=sort_first_aft;
  sort_first_aft_paired(Indices)= sort_first_aft_paired;
  
  [sort_sec_bef,Indices]=sort(squeeze(mean(sort_sec_bef,2)),2);
  sort_sec_aft=sort(squeeze(mean( sort_sec_aft,2)),2);
  sort_sec_aft_paired=sort_sec_aft;
  sort_sec_aft_paired(Indices)= sort_sec_aft_paired;
  
  [sort_third_bef,Indices]=sort(squeeze(mean( sort_third_bef,2)),2);
  sort_third_aft=sort(squeeze(mean( sort_third_aft,2)),2);
 sort_third_aft_paired=sort_third_aft;
  sort_third_aft_paired(Indices)= sort_third_aft_paired;
  
  
 mean_sort_first_bef=mean(sort_first_bef,1);
 mean_sort_first_aft=mean(sort_first_aft,1);
 mean_sort_sec_bef=mean(sort_sec_bef,1);
 mean_sort_sec_aft=mean(sort_sec_aft,1);
 mean_sort_third_bef=mean(sort_third_bef,1);
 mean_sort_third_aft=mean(sort_third_aft,1); 
 
 std_sort_first_bef=std(sort_first_bef,1)/sqrt(length(sort_first_bef));
 std_sort_first_aft=std(sort_first_aft,1)/sqrt(length(sort_first_aft));
  std_sort_sec_bef=std(sort_sec_bef,1)/sqrt(length(sort_sec_bef));
 std_sort_sec_aft=std(sort_sec_aft,1)/sqrt(length(sort_sec_aft));
 std_sort_third_bef=std(sort_third_bef,1)/sqrt(length(sort_third_bef));
 std_sort_third_aft=std(sort_third_aft,1)/sqrt(length(sort_third_aft)); 
  
 
  %first:
 figure
errorbar (x,mean_sort_first_bef,std_sort_first_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_first_aft,std_sort_first_aft,'r','linewidth',2)
title('responses sorted by magnitude-firsteriment')
 box off
 ylim([-0.15 0.5])
%  ylim([-0.1 0.45]) for anes
%  small_dif=mean_sort_first_bef(4:7)-mean_sort_first_aft(4:7);
%  small_dif_prop=small_dif./mean_sort_first_bef(4:7);
%  big_dif=mean_sort_first_bef(8:11)-mean_sort_first_aft(8:11);
%  big_dif_prop=big_dif./mean_sort_first_bef(8:11);
%  
%  [h,p]=ranksum(small_dif_prop,big_dif_prop);
%  [h,p]=ttest2(small_dif_prop,big_dif_prop);
 
 %sec
  figure
errorbar (x,mean_sort_sec_bef,std_sort_sec_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_sec_aft,std_sort_sec_aft,'r','linewidth',2)
title('responses sorted by magnitude-sec')
 box off
 ylim([-0.15 0.5])
 %third:
   figure
errorbar (x,mean_sort_third_bef,std_sort_third_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_third_aft,std_sort_third_aft,'r','linewidth',2)
 title('responses sorted by magnitude-thirdrol')
box off
 ylim([-0.15 0.5])
%  ylim([-0.1 0.45]) for anes
 

%% Deltas_before_after_for_strongest+weakest responses

first_deltas_first=(sort_first_aft(:,1)-sort_first_bef(:,1));
first_deltas_last=(sort_first_aft(:,11)-sort_first_bef(:,11));

sec_deltas_first=(sort_first_aft(:,1)-sort_first_bef(:,1));
sec_deltas_last=(sort_first_aft(:,11)-sort_first_bef(:,11));

third_deltas_first=(sort_third_aft(:,1)-sort_third_bef(:,1));
third_deltas_last=(sort_third_aft(:,11)-sort_third_bef(:,11));

% first_CI_first=((sort_first_aft(:,1)-sort_first_bef(:,1)))./(abs(sort_first_aft(:,1))+abs(sort_first_bef(:,1)));
% first_CI_last=((sort_first_aft(:,11)-sort_first_bef(:,11)))./(abs(sort_first_aft(:,11))+abs(sort_first_bef(:,11)));
% 
% 
% third_CI_first=((sort_third_aft(:,1)-sort_third_bef(:,1)))./(abs(sort_third_aft(:,1))+abs(sort_third_bef(:,1)));
% third_CI_last=((sort_third_aft(:,11)-sort_third_bef(:,11)))./(abs(sort_third_aft(:,11))+abs(sort_third_bef(:,11)));


sem_reg_first_first=std(first_deltas_first)/sqrt(length(first_deltas_first));
sem_reg_sec_first=std(sec_deltas_first)/sqrt(length(sec_deltas_first));
sem_reg_third_first=std(third_deltas_first)/sqrt(length(third_deltas_first));


sem_reg_first_last=std(first_deltas_last)/sqrt(length(first_deltas_last));
sem_reg_sec_last=std(sec_deltas_last)/sqrt(length(sec_deltas_last));
sem_reg_third_last=std(third_deltas_last)/sqrt(length(third_deltas_last));

sems_reg_first=[sem_reg_first_first sem_reg_sec_first sem_reg_third_first]
sems_reg_last=[sem_reg_first_last sem_reg_sec_last sem_reg_third_last]

means_reg_first=[mean(first_deltas_first) mean(sec_deltas_first) mean(third_deltas_first)];
means_reg_last=[mean(first_deltas_last) mean(sec_deltas_last) mean(third_deltas_last)];
y=6;
x=1:3;

x=1:3;
figure
hold on



bar(x(1),means_reg_first(1),'k')
bar(x(2),means_reg_first(2),'k')
bar(x(3),means_reg_first(3),'k')
errorbar(x,means_reg_first,sems_reg_first,'r.')

ax=gca;
ax.FontSize=14;
ax.XTick= []
ax.XTickLabels=[]

ylim([min([means_reg_first means_reg_last])-max([sems_reg_first sems_reg_last]) max([means_reg_first means_reg_last])+max([sems_reg_first sems_reg_last])])

title('Response magnitudes Deltas per cell, regular-sorted first. third vs first')
%%
%
% x1 = linspace(4,5,size(first_deltas_first,1));
% x2 = linspace(6,7,size(third_deltas_first,1));
% 
% mean_vec_first=mean(means_reg_first(1))*ones(size(x1))
% mean_vec_third=mean(means_reg_first(2))*ones(size(x2))
% 
% x3=3:8;
% y=zeros(size(x3))
% 
% figure
% title('Norm deltas regular distances')
% xlim([3 8])
% hold on
% plot(x1,first_deltas_first,'k*')
% plot(x2,third_deltas_first,'ko')
% plot(x1,mean_vec_first,'r--','linewidth',2)
% plot(x2,mean_vec_third,'r--','linewidth',2)
% plot(x3,y,'k--')
% XTick=[]
% 
% sem_first=sem_reg_first_first
% sem_third=sem_reg_third_first
% 
% x_er=[4.5 6.5];
% y_er=[means_reg_first(1) means_reg_first(2)];
% sems=[sem_third sem_first];
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
xlabel('first              third')
title('Response magnitudes Deltas per cell, regular-sorted last. third vs first')
ylim([min([means_reg_first means_reg_last])-max([sems_reg_first sems_reg_last]) max([means_reg_first means_reg_last])+max([sems_reg_first sems_reg_last])])

%now for CIs
first_CI_first=((sort_first_aft(:,1)-sort_first_bef(:,1)))./(abs(sort_first_aft(:,1))+abs(sort_first_bef(:,1)));
first_CI_last=((sort_first_aft(:,11)-sort_first_bef(:,11)))./(abs(sort_first_aft(:,11))+abs(sort_first_bef(:,11)));


third_CI_first=((sort_third_aft(:,1)-sort_third_bef(:,1)))./(abs(sort_third_aft(:,1))+abs(sort_third_bef(:,1)));
third_CI_last=((sort_third_aft(:,11)-sort_third_bef(:,11)))./(abs(sort_third_aft(:,11))+abs(sort_third_bef(:,11)));


sem_CI_first_first=std(first_CI_first)/sqrt(length(first_CI_first));
sem_CI_third_first=std(third_CI_first)/sqrt(length(third_CI_first));


sem_CI_first_last=std(first_CI_last)/sqrt(length(first_CI_last));
sem_CI_third_last=std(third_CI_last)/sqrt(length(third_CI_last));

sems_CI_first=[sem_CI_first_first sem_CI_third_first];
sems_CI_last=[sem_CI_first_last sem_CI_third_last];

means_CI_first=[mean(first_CI_first) mean(third_CI_first)];
means_CI_last=[mean(first_CI_last) mean(third_CI_last)];

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

title('Response magnitudes Deltas per cell, CIs-sorted first. third vs first')


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

title('Response magnitudes Deltas per cell, CIs-sorted last. third vs first')

%% same like last but by absolute values

%% all trials magnitudes together histograms for magnitude- sorted and not sorted by oddor strength
Total_first_mat_for_trials_magnitudes_before=[];
Total_first_mat_for_trials_magnitudes_after=[];

Total_third_mat_for_trials_magnitudes_before=[];
Total_third_mat_for_trials_magnitudes_after=[];

for j=1:Num_of_first_mice
 Total_first_mat_for_trials_magnitudes_before=[Total_first_mat_for_trials_magnitudes_before;first_data.mice_table{1,j}.mat_for_diffrence_test] ; 
 Total_first_mat_for_trials_magnitudes_after=[Total_first_mat_for_trials_magnitudes_after;first_data.mice_table{2,j}.mat_for_diffrence_test];  
end

for j=1:Num_of_third_mice
 Total_third_mat_for_trials_magnitudes_before=[Total_third_mat_for_trials_magnitudes_before;sec_data.mice_table{1,j}.mat_for_diffrence_test];  
 Total_third_mat_for_trials_magnitudes_after=[Total_third_mat_for_trials_magnitudes_after;sec_data.mice_table{2,j}.mat_for_diffrence_test];  
end

%added here down : for each not to reduce dimension

%  Total_first_mat_for_trials_magnitudes_before_organized= Total_first_mat_for_trials_magnitudes_before(:,:,k_for_plot);
%  Total_first_mat_for_trials_magnitudes_after_organized= Total_first_mat_for_trials_magnitudes_after(:,:,k_for_plot);
%  Total_third_mat_for_trials_magnitudes_before_organized= Total_third_mat_for_trials_magnitudes_before(:,:,k_for_plot);
%  Total_third_mat_for_trials_magnitudes_after_organized= Total_third_mat_for_trials_magnitudes_after(:,:,k_for_plot);
%  
%  Total_first_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
%  Total_first_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
%  Total_third_mat_for_trials_magnitudes_before_organized(:,:,1)=[];
%  Total_third_mat_for_trials_magnitudes_after_organized(:,:,1)=[];
 
Total_first_mat_for_trials_magnitudes_before_organized_abs=[];  %delete only this part to get them not as abs
Total_first_mat_for_trials_magnitudes_after_organized_abs=[];
Total_third_mat_for_trials_magnitudes_before_organized_abs=[];
Total_third_mat_for_trials_magnitudes_after_organized_abs=[];
 
Total_first_mat_for_trials_magnitudes_before_organized_abs=abs(Total_first_mat_for_trials_magnitudes_before_organized);  %delete only this part to get them not as abs
Total_first_mat_for_trials_magnitudes_after_organized_abs=abs(Total_first_mat_for_trials_magnitudes_after_organized);
Total_third_mat_for_trials_magnitudes_before_organized_abs=abs(Total_third_mat_for_trials_magnitudes_before_organized);
Total_third_mat_for_trials_magnitudes_after_organized_abs=abs(Total_third_mat_for_trials_magnitudes_after_organized);
 %%for first not sorted by magnitude for each cell:
 
  x=1:11;
  y1=(mean(Total_first_mat_for_trials_magnitudes_before_organized_abs,2));
  std_y1=squeeze(std(y1,1));
  sem_y1=std_y1/sqrt(length(y1));
  y1=mean(y1,1);
  y1=squeeze(y1);

  y2=(mean(Total_first_mat_for_trials_magnitudes_after_organized_abs,2))
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
 title('magnitude of response per odor-firsteriment')
 xlabel('Odor identity')
 ylabel('response magnitude')
  box off
 ylim([-0.05 0.4])
 ylim([0 0.2])  %for integral
 %%for third not sorted by magnitude for each cell:
 
 x=1:11;
  y3=(mean(Total_third_mat_for_trials_magnitudes_before_organized_abs,2));
  std_y3=squeeze(std(y3,1));
  sem_y3=std_y3/sqrt(length(y3));
  y3=mean(y3,1);
  y3=squeeze(y3);

  y4=(mean(Total_third_mat_for_trials_magnitudes_after_organized_abs,2));
  std_y4=squeeze(std(y4,1));
  sem_y4=std_y4/sqrt(length(y4));
  y4=mean(y4,1);
  y4=squeeze(y4);
  

 
 
  figure
 errorbar(x,y3,sem_y3,'b','linewidth',2)
 hold on
 errorbar(x,y4,sem_y4,'r','linewidth',2)
 box off
 ylim([-0.05 0.4])
  ylim([0 0.2])  %for integral
 title('magnitude of response per odor-thirdrol')
  xlabel('Odor identity')
 ylabel('response magnitude')
 % now sorted by magnitude:
 
 %  %sort by magnitude:
 sort_first_bef= Total_first_mat_for_trials_magnitudes_before_organized_abs;
 sort_first_aft=Total_first_mat_for_trials_magnitudes_after_organized_abs;
 sort_third_bef=Total_third_mat_for_trials_magnitudes_before_organized_abs;
 sort_third_aft=Total_third_mat_for_trials_magnitudes_after_organized_abs;
 
 [sort_first_bef,Indices]=sort(squeeze(mean(sort_first_bef,2)),2);
  sort_first_aft=sort(squeeze(mean( sort_first_aft,2)),2);
  sort_first_aft_paired=sort_first_aft;
  sort_first_aft_paired(Indices)= sort_first_aft_paired;
  
  [sort_third_bef,Indices]=sort(squeeze(mean( sort_third_bef,2)),2);
  sort_third_aft=sort(squeeze(mean( sort_third_aft,2)),2);
 sort_third_aft_paired=sort_third_aft;
  sort_third_aft_paired(Indices)= sort_third_aft_paired;
  
  
 mean_sort_first_bef=mean(sort_first_bef,1);
 mean_sort_first_aft=mean(sort_first_aft,1);
 mean_sort_third_bef=mean(sort_third_bef,1);
 mean_sort_third_aft=mean(sort_third_aft,1); 
 
 std_sort_first_bef=std(sort_first_bef,1)/sqrt(length(sort_first_bef));
 std_sort_first_aft=std(sort_first_aft,1)/sqrt(length(sort_first_aft));
 std_sort_third_bef=std(sort_third_bef,1)/sqrt(length(sort_third_bef));
 std_sort_third_aft=std(sort_third_aft,1)/sqrt(length(sort_third_aft)); 
  
 
  %first:
 figure
errorbar (x,mean_sort_first_bef,std_sort_first_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_first_aft,std_sort_first_aft,'r','linewidth',2)
title('responses sorted by magnitude-firsteriment')
 box off
 ylim([-0.2 0.2])
  ylim([0 0.3])
 small_dif=mean_sort_first_bef(4:7)-mean_sort_first_aft(4:7);
 small_dif_prop=small_dif./mean_sort_first_bef(4:7);
 big_dif=mean_sort_first_bef(8:11)-mean_sort_first_aft(8:11);
 big_dif_prop=big_dif./mean_sort_first_bef(8:11);
 
 [h,p]=ranksum(small_dif_prop,big_dif_prop);
 [h,p]=ttest2(small_dif_prop,big_dif_prop);
 
 %third:
   figure
errorbar (x,mean_sort_third_bef,std_sort_third_bef,'b','linewidth',2)
 hold on
errorbar(x,mean_sort_third_aft,std_sort_third_aft,'r','linewidth',2)
 title('responses sorted by magnitude-thirdrol')
box off
 ylim([-0.2 0.2])
  ylim([0 0.3])
 



 %% kato:
 resp_thresh=3; %include only cells that had this number of exc responses at least (both before and after)
 
 sort_first_bef_norm=sort_first_bef;
 sort_first_aft_norm=sort_first_aft_paired;
 sort_third_bef_norm=sort_third_bef;
 sort_third_aft_norm=sort_third_aft_paired;
 
 for i=1:size(sort_first_bef,1);
   sort_first_bef_norm(i,:)= sort_first_bef_norm(i,:)./sort_first_bef(i,end);
   sort_first_aft_norm(i,:)= sort_first_aft_norm(i,:)./sort_first_bef(i,end);
 end
 
 for i=1:size(sort_third_bef,1);
   sort_third_bef_norm(i,:)= sort_third_bef_norm(i,:)./sort_third_bef(i,end);
   sort_third_aft_norm(i,:)= sort_third_aft_norm(i,:)./sort_third_bef(i,end);
 end
 
meta_bef_af_first=general_first.all_total_meta_des_bef_af_inh_ex;
meta_bef_af_third=general_sec.all_total_meta_des_bef_af_inh_ex;

sum_meta_first=sum(meta_bef_af_first~=0,2)-sum(meta_bef_af_first<=-2,2);    %dont include inhibition
sum_meta_third=sum(meta_bef_af_third~=0,2)-sum(meta_bef_af_third<=-2,2);
 
sort_first_bef_norm(sum_meta_first<resp_thresh,:)=[];
sort_first_aft_norm(sum_meta_first<resp_thresh,:)=[];

sort_third_bef_norm(sum_meta_third<resp_thresh,:)=[];
sort_third_aft_norm(sum_meta_third<resp_thresh,:)=[];

sort_first_bef_norm=fliplr(sort_first_bef_norm);
sort_first_aft_norm=fliplr(sort_first_aft_norm);

sort_third_bef_norm=fliplr(sort_third_bef_norm);
sort_third_aft_norm=fliplr(sort_third_aft_norm);


x=[6 4 2 1 3 5 7];
x2=[4 3 5 2 6 1 7];
x3=1:7;

y_bef=mean(sort_first_bef_norm);
y_bef(8:end)=[];
y_bef=y_bef(x);

sem_bef=std(sort_first_bef_norm)./sqrt(size(sort_first_bef_norm,1));
sem_bef=sem_bef(x);

y_aft=mean(sort_first_aft_norm);
y_aft(8:end)=[];
y_aft=y_aft(x);

sem_aft=std(sort_first_aft_norm)./sqrt(size(sort_first_aft_norm,1));
sem_aft=sem_aft(x);

%% Now linear regression for each cell sperately:
aaa=(-1:5);
bbb=(-1:5);
xx=-1:length(aaa);
x_ax=zeros(length(aaa)+2,1);

xxx=[0 0 0];
yyy=[-10 0 10];

b1_for_all_first=[];
b_for_all_first=[];

figure
for i=1:size(sort_first_bef_norm,1);
    %figure      % cancel this to plot some of the cells

    x= sort_first_bef_norm(i,1:7)';
    y= sort_first_aft_norm(i,1:7)';
    
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

b1_for_all_first=[b1_for_all_first;b1];
b_for_all_first(1,i)=b(1);
b_for_all_first(2,i)=b(2);

title([' Intercept=' num2str(b(1)) ' Slope= ' num2str(b(2))])
end

%analyse intercept and slopes
figure
hold on
yy1=zeros(size( sort_first_bef_norm,1),1);
yy2=ones(size( sort_first_bef_norm,1),1);
% x1=ones(size( sort_first_bef_norm,1),1);
% x2=ones(size( sort_first_bef_norm,1),1)*2;
x1 = linspace(1,3,size( sort_first_bef_norm,1));

x2 = linspace(6,8,size( sort_first_bef_norm,1));

plot(x1,b_for_all_first(1,:),'*')
plot(x2,b_for_all_first(2,:),'*')
plot(x1,yy1,'k--');
plot(x2,yy2,'k--');

xlim([0 9])



h_first=mean(b_for_all_first,2);

%kato
figure
hold on
title(['firsteriment-Kato analysis' num2str(size( sort_first_bef_norm,1)) 'cells'])
ylim([0 1.2])
errorbar(x3,y_bef,sem_bef,'b','linewidth',2)
errorbar(x3,y_aft,sem_aft,'r','linewidth',2)
y_linear_prediction=h_first(1)+y_bef.*h_first(2);
plot(x3,y_linear_prediction,'k--','linewidth',2)


%%

%same for thirdrol
xx=-1:length(aaa);
x_ax=zeros(length(aaa)+2,1);

xxx=[0 0 0];
yyy=[-10 0 10];

x=[6 4 2 1 3 5 7];

y_bef=mean(sort_third_bef_norm);
y_bef(8:end)=[];
y_bef=y_bef(x);

sem_bef=std(sort_third_bef_norm)./sqrt(size(sort_third_bef_norm,1));
sem_bef=sem_bef(x);

y_aft=mean(sort_third_aft_norm);
y_aft(8:end)=[];
y_aft=y_aft(x);

sem_aft=std(sort_third_aft_norm)./sqrt(size(sort_third_aft_norm,1));
sem_aft=sem_aft(x);

aaa=(-1:5);
bbb=(-1:5);
b1_for_all_third=[];
b_for_all_third=[];

figure
for i=1:size( sort_third_bef_norm,1);
%    figure      % cancel this to plot some of the cells

    x= sort_third_bef_norm(i,1:7)';
    y= sort_third_aft_norm(i,1:7)';
    
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
    
b1_for_all_third=[b1_for_all_third;b1];
b_for_all_third(1,i)=b(1);
b_for_all_third(2,i)=b(2);

title([' Intercept=' num2str(b(1)) ' Slope= ' num2str(b(2))])
end

%analyse intercept and slopes
figure
xlim([0 9])
ylim([-2 4])


hold on
x1 = linspace(1,3,size( sort_third_bef_norm,1));
x2 = linspace(6,8,size( sort_third_bef_norm,1));
plot(x1,b_for_all_third(1,:),'o')
plot(x2,b_for_all_third(2,:),'o')


yy1=zeros(size( sort_third_bef_norm,1),1);
yy2=ones(size( sort_third_bef_norm,1),1);
plot(x1,yy1,'k--');
plot(x2,yy2,'k--');
ylim([-2 4])





%kato
h_third=mean(b_for_all_third,2)

figure
title(['thirdrol-Kato analysis ' num2str(size( sort_third_bef_norm,1)) 'cells'])
hold on
ylim([0 1.2])
errorbar(x3,y_bef,sem_bef,'b','linewidth',2)
errorbar(x3,y_aft,sem_aft,'r','linewidth',2)
y_linear_prediction=h_third(1)+y_bef.*h_third(2);
plot(x3,y_linear_prediction,'k--','linewidth',2)
hold off

disp('no')

%% intercepts and slopes sep first third combined:
figure
ylim([-1 1])
hold on
y1=zeros(size( sort_first_bef_norm,1),1);
y2=ones(size( sort_first_bef_norm,1),1);
xx1 = linspace(1,2,size( sort_first_bef_norm,1));
xx2 = linspace(7,8,size( sort_first_bef_norm,1));

x1 = linspace(4,5,size( sort_third_bef_norm,1));
x2 = linspace(10,11,size( sort_third_bef_norm,1));
yy1=zeros(size( sort_third_bef_norm,1),1);
yy2=ones(size( sort_third_bef_norm,1),1);


plot(xx1,b_for_all_first(1,:),'k*')
plot(x1,b_for_all_third(1,:),'ko')

plot(xx1,y1,'r--');
plot(x1,yy1,'r--');
xlim([0 6])

figure
hold on
ylim([0 4])
xlim([6 12])

plot(xx2,b_for_all_first(2,:),'k*')
plot(x2,b_for_all_third(2,:),'ko')
plot(xx2,y2,'r--');
plot(x2,yy2,'r--');




mean(b_for_all_third(1,:))
mean(b_for_all_first(1,:))
mean(b_for_all_third(2,:))
mean(b_for_all_first(2,:))

[p,h]=ranksum(b_for_all_third(1,:),b_for_all_first(1,:))
[p,h]=ranksum(b_for_all_third(2,:),b_for_all_first(2,:))



%%
 sort_first_bef_norm
 sort_first_aft_norm=sort_first_aft_paired;
 sort_third_bef_norm=sort_third_bef;
 sort_third_aft_norm=sort_third_aft_paired;
 

 %% Variance between trials
%  STD_mat_first_before=zeros(size(first_matrix_before_organized,1),Num_of_odors);
%  STD_mat_first_after=zeros(size(first_matrix_before_organized,1),Num_of_odors);
%  STD_mat_third_before=zeros(size(third_matrix_before_organized,1),Num_of_odors);
%  STD_mat_third_after=zeros(size(third_matrix_before_organized,1),Num_of_odors);
%  
%  for i=1:size(first_matrix_before_organized,1);
%      for j=1:Num_of_odors;
%           STD_mat_first_before(i,j)=std(Total_first_mat_for_trials_magnitudes_before(i,:,j));
%           STD_mat_first_after(i,j)=std(Total_first_mat_for_trials_magnitudes_after(i,:,j));
%           
%      end
%  end
% a= mean(STD_mat_first_before,2);
% b=mean(STD_mat_first_after,2);
%  
%  STD_mat_third_before=zeros(size(third_matrix_before_organized,1),Num_of_odors);
%  STD_mat_third_after=zeros(size(third_matrix_before_organized,1),Num_of_odors);
%  
%  for i=1:size(third_matrix_before_organized,1);
%      for j=1:Num_of_odors
%           STD_mat_third_before(i,j)=std(Total_third_mat_for_trials_magnitudes_before(i,:,j));
%           STD_mat_third_after(i,j)=std(Total_third_mat_for_trials_magnitudes_after(i,:,j));
%           
%      end
%  end
% c= mean(STD_mat_third_before,2);
% d=mean(STD_mat_third_after,2);
%  
%% 
Total_signif_desicion_first=[];

for j=1:Num_of_first_mice;
% temp1=first_data.mice_table{1,j}.meta_desicion>1;
% temp2=first_data.mice_table{2,j}.meta_desicion>1;
 temp1=first_data.mice_table{1,j}.meta_desicion~=0;
 temp2=first_data.mice_table{2,j}.meta_desicion~=0;

temp3=temp1+temp2;
Total_signif_desicion_first=[Total_signif_desicion_first;temp3];
end


Total_signif_desicion_first_organized=Total_signif_desicion_first(:,k_for_plot);
Total_signif_desicion_first_organized(:,1)=[];

Total_signif_desicion_third=[];
for j=1:Num_of_third_mice;
temp1=sec_data.mice_table{1,j}.meta_desicion~=0;
temp2=sec_data.mice_table{2,j}.meta_desicion~=0;
% temp1=third_data.mice_table{1,j}.meta_desicion>1;
% temp2=third_data.mice_table{2,j}.meta_desicion>1;

temp3=temp1+temp2;
Total_signif_desicion_third=[Total_signif_desicion_third;temp3];
end


Total_signif_desicion_third_organized=Total_signif_desicion_third(:,k_for_plot);
Total_signif_desicion_third_organized(:,1)=[];
% %%
% a=STD_mat_first_before(Total_signif_desicion_first_organized>0);
% b=STD_mat_first_after(Total_signif_desicion_first_organized>0);
% 
% c=STD_mat_third_before(Total_signif_desicion_third_organized>0);
% d=STD_mat_third_after(Total_signif_desicion_third_organized>0);
% 
% [h,p]=ttest(a,b);
% [h,p]=ttest(c,d);
% j=a-b;
% k=c-d;
% [h,p]=ttest2(j,k);
% x=[1 2];
% y=[mean(a) mean(b)];
% z=[mean(c) mean(d)];
% 
% e1=[std(a)/sqrt(length(a)) std(b)/sqrt(length(b))]
% e2=[std(c)/sqrt(length(c)) std(d)/sqrt(length(d))]
% 
% figure
% % plot(x,y);
% % hold on
% % plot(x,z);
% hold on 
% errorbar(x,y, e1)
%  errorbar(x,z, e2)
% 
% title(['Variability between trials, thirdrol vs firsteriment. P=' num2str(p)])
% xlabel('Before                                               After','FontSize',14)
% ylabel('Variability between trials (SD)','FontSize',14)
% legend('first','third', 'Location', 'northwest')

%%
 %% chi_test_for_dots_plot- proportion of cells
third_observ=general_sec.stats_for_chi;
first_observ=general_first.stats_for_chi;

  %first all then all non signif then those that are strongerbefore then stronger after
       third_non_sig = third_observ(2); third_bef=third_observ(3);third_af=third_observ(4); all_resp_third = third_observ(1);   %even vs odd
       first_non_sig = first_observ(2); first_bef=first_observ(3);first_af=first_observ(4); all_resp_first = first_observ(1); 

        % Pooled estimate of proportion
       p_non_sig= (third_non_sig+first_non_sig)/(all_resp_third+all_resp_first);
       p_bef= (third_bef+first_bef)/(all_resp_third+all_resp_first);
       p_af=(first_af+third_af)/(all_resp_third+all_resp_first);
       
       % firstected counts under H0 (null hypothesis)
       p_non_sig_third_firstected= p_non_sig*all_resp_third;
       p_bef_third_firstected=p_bef*all_resp_third;
       p_af_third_firstected=p_af*all_resp_third;
       
       p_non_first_first_firstected= p_non_sig*all_resp_first;
       p_bef_first_firstected=p_bef*all_resp_first;
       p_af_first_firstected=p_af*all_resp_first;
        
       a_c=p_non_sig_third_firstected;  %just for comfort
       b_c= p_bef_third_firstected;
       c_c=p_af_third_firstected;
       
       a_e=p_non_first_first_firstected;
       b_e= p_bef_first_firstected;
       c_e= p_af_first_firstected;
       
       firstected=[a_c b_c c_c a_e b_e c_e]
       observed=[third_observ(2:end) first_observ(2:end)]
    
       % Chi-square test, by hand
      
       chi2stat = sum((observed-firstected).^2 ./ firstected)
       p = 1 - chi2cdf(chi2stat,2) %degrees of freedom is (rows-1)*(colmuns-1)
       
       %normalization
       third_observ_norm=third_observ/third_observ(1)
       first_observ_norm=first_observ/first_observ(1)
       
       %plot
       x=1:6;
       figure
       
     
       bar(x(1:2:5),third_observ(2:end),0.45)
       hold on
       bar(x(2:2:6),first_observ(2:end),0.45,'k')
      
       ax=gca;
       ax.FontSize=14;
       ax.XTick=[];
      
      
      
      title(['Changes in response magnitude. P=' num2str(p)])
      xlabel('Non significant      Stronger Before        Stronger after','FontSize',14)
      ylabel('Number of responses','FontSize',14)
      legend('thirdrol','firsteriment', 'Location', 'northeast')


       %plot normalized bar
       x=1:6;
       figure
       
     
       bar(x(1:2:5),third_observ_norm(2:end),0.45)
       hold on
       bar(x(2:2:6),first_observ_norm(2:end),0.45,'k')
      
       ax=gca;
       ax.FontSize=14;
       ax.XTick=[];
      
      
      
      title(['Changes in response magnitude. P=' num2str(p)])
      xlabel('Non significant      Stronger Before        Stronger after','FontSize',14)
      ylabel('Number of responses','FontSize',14)
      legend('thirdrol','firsteriment', 'Location', 'northeast')
%% mice deltas

deltas_third=general_sec.mice_deltas;
deltas_first=general_first.mice_deltas;



figure
for i=1:length(deltas_third)
    plot(1, deltas_third(i),'ko','linewidth', 2)
    hold on
end


for i=1:length(deltas_first)
    plot(2, deltas_first(i),'ro','linewidth', 2)
    hold on
end
y=zeros(1,4);
x=(0:3);
plot(x,y,'k--')

xlim([0 3])
ylim([-0.2 0.2])


mean_third=mean(deltas_third);
mean_first=mean(deltas_first);

means=[mean_third,mean_first];

sem_third=std(deltas_third)/sqrt(length(deltas_third));
sem_first=std(deltas_first)/sqrt(length(deltas_first));

sems=[sem_third,sem_first];

errorbar(means, sems,'.','linewidth', 2);


%%
%% prop (ratios)

deltas_third_prop=general_sec.prop_for_comparison;
deltas_first_prop=general_first.prop_for_comparison;


figure
for i=1:length(deltas_third_prop);
    plot(1, deltas_third_prop(i),'ko','linewidth', 2)
    hold on
end


for i=1:length(deltas_first_prop);
    plot(2, deltas_first_prop(i),'ro','linewidth', 2);
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

mean_third=mean(deltas_third_prop);
mean_first=mean(deltas_first_prop);

means=[mean_third,mean_first];

sem_third=std(deltas_third_prop)/sqrt(length(deltas_third_prop));
sem_first=std(deltas_first_prop)/sqrt(length(deltas_first_prop));

sems=[sem_third,sem_first];

errorbar(means, sems,'.','linewidth', 2);


%ax.xtick=[];

% XTick off
% YTick off

%%
deltas_third_prop=general_sec.prop_for_comparison;
deltas_first_prop=general_first.prop_for_comparison;


figure
for i=1:length(deltas_third_prop)
    plot(1, deltas_third_prop(i),'ko')
    hold on
end


for i=1:length(deltas_first_prop);
    plot(2, deltas_first_prop(i),'ro')
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

mean_third=mean(deltas_third_prop);
mean_first=mean(deltas_first_prop);

means=[mean_third,mean_first];

sem_third=std(deltas_third_prop)/length(deltas_third_prop);
sem_first=std(deltas_first_prop)/length(deltas_first_prop);

sems=[sem_third,sem_first];

errorbar(means, sems,'.')


%ax.xtick=[];

% XTick off
% YTick off
%%

% first_deltas_bp=general_first.deltas_for_boxplots;
% third_deltas_bp=general_third.deltas_for_boxplots;

first_deltas_bp=general_first.abs_deltas_for_boxplots;
third_deltas_bp=general_sec.abs_deltas_for_boxplots;


%% boxplots
kx=[0:4];
ky=[0 0 0 0 0];

figure
%x=[all_significant_responses_integrals_before;all_significant_responses_integrals_after];
x=[third_deltas_bp;first_deltas_bp];

g = [ones(size(third_deltas_bp)); 2*ones(size(first_deltas_bp))];
% boxplot(x,g)
boxplot(x,g,'symbol','');

%title(['Before and After ' num2str(length(all_significant_responses_magnitudes_before)) ' significant responses ' ' Before Mean=' num2str(mean_before) ' SD=' num2str(SD_before) ' After Mean=' num2str(mean_after) ' SD=' num2str(SD_after) '    p=' num2str(p)],'FontSize',24);% '    p2=' num2str(p_for_all)]);
%ylabel('Response magnitude')
%xlabel('thirdrol                                                       firsteriment')
%ylim([0 0.6])
ylim([-0.3 0.3])
a=gca;
a.FontSize=20;
hold on
plot(kx,ky,'k--','linewidth',2)

[h,p]=ttest2(third_deltas_bp,first_deltas_bp)
[p,h]=ranksum(third_deltas_bp,first_deltas_bp)
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
resp_per_cell_first_before=general_first.total_prop_of_cells_before;
resp_per_cell_first_after=general_first.total_prop_of_cells_after;

resp_per_cell_third_before=general_sec.total_prop_of_cells_before;
resp_per_cell_third_after=general_sec.total_prop_of_cells_after;

first=resp_per_cell_first_after-resp_per_cell_first_before;
third=resp_per_cell_third_after-resp_per_cell_third_before;

first_norm_general=mean(resp_per_cell_first_before);
third_norm_general=mean(resp_per_cell_third_before);

y=[mean(first),mean(third)];
sem_first=std(first)/sqrt(size(first,1));
sem_third=std(third)/sqrt(size(third,1));
sems=[sem_first sem_third]

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

title('Change in responsiveness # of responses per cell')

y=[mean(first)/first_norm_general,mean(third)/third_norm_general];
sem_first=std(first)/sqrt(size(first,1));
sem_third=std(third)/sqrt(size(third,1));
sems_norm=[sem_first/first_norm_general sem_third/third_norm_general];

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

first_prop=(resp_per_cell_first_after-resp_per_cell_first_before)./(resp_per_cell_first_before+resp_per_cell_first_after);
third_prop=(resp_per_cell_third_after-resp_per_cell_third_before)./(resp_per_cell_third_before+resp_per_cell_third_after);
%third_prop(third_prop==inf)=1;
third_prop(isnan(third_prop))=0;
%first_prop(first_prop==inf)=1;
first_prop(isnan(first_prop))=0;

y=[mean(first_prop) mean(third_prop)];
sem_first_prop=std(first_prop)/sqrt(size(first_prop,1));
sem_third_prop=std(third_prop)/sqrt(size(third_prop,1));
sems_prop=[sem_first_prop sem_third_prop];

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
resp_over_time_first_ex=general_first.deltas_over_time_ex;
resp_over_time_first_in=general_first.deltas_over_time_in;

resp_over_time_third_ex=general_sec.deltas_over_time_ex;
resp_over_time_third_in=general_sec.deltas_over_time_in;


first_ex=mean(resp_over_time_first_ex(:,49:77),2);
third_ex=mean(resp_over_time_third_ex(:,49:77),2);

first_in=mean(resp_over_time_first_in(:,49:77),2);
third_in=mean(resp_over_time_third_in(:,49:77),2);





y_ex=[mean(first_ex),mean(third_ex)].*(-1);     %to make after-beofre (flip) 
sem_first_ex=std(first_ex)/sqrt(size(first_ex,1));
sem_third_ex=std(third_ex)/sqrt(size(third_ex,1));
sems_ex=[sem_first_ex sem_third_ex];

y_in=[mean(first_in),mean(third_in)].*(-1);              %to make after-beofre (flip) 
sem_first_in=std(first_in)/sqrt(size(first_in,1));
sem_third_in=std(third_in)/sqrt(size(third_in,1));
sems_in=[sem_first_in sem_third_in];


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

title('Response magnitudes Deltas')
ylim([min([y_in y_ex])-max([sems_ex sems_in]) max([y_in y_ex])+max([sems_ex sems_in])])



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

%% chi test for the rasters:

%first all then all non signif then those that are strongerbefore then
%stronger after:
general_first.stats_for_chi(1)
N1=general_first.stats_for_chi(1); n1=general_first.stats_for_chi(3)+general_first.stats_for_chi(4)
N2= general_sec.stats_for_chi(1); n2=general_sec.stats_for_chi(3)+general_sec.stats_for_chi(4)
N3=general_third.stats_for_chi(1); n3=general_third.stats_for_chi(3)+general_third.stats_for_chi(4)

    p0 = (n1+n2+n3) / (N1+N2+N3)
       % firstected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       n30=N3*p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2 n3 N3-n3];
       firstected = [n10 N1-n10 n20 N2-n20 n30 N3-n30];
       chi2stat = sum((observed-firstected).^2 ./ firstected);
       p = 1 - chi2cdf(chi2stat,1)
     
%report this