clc
clear all
close all

stats=struct;
    
    
exp_mice=[1,2,5,9,10,39,12,13,32,36];   %Anesthesized
cont_mice=[30,34,35,37,38,11,7,15,16,40]; %Anesthesized
awake_mice=[17,50,200,52,54] ;
mice_awake_cont=[1000 1001 1002 1003 1004];
first_timelapse=[100 200 300 500 600];
sec_timelapse=[101 201 301 601];
third_timelapse=[102 202 302 502];
cno_cont_mice=[10102 1011];  

    

%%%%%%%%change only here:

% mice=cont_mice;        %Fig 2 CONT
% mice=exp_mice;         %Fig 2 EXP
% mice=awake_mice;       %Fig 3 EXP
% mice=mice_awake_cont;  %Fig 3 CONT
 mice=first_timelapse;  %Fig 5 8WPI
% mice=sec_timelapse;    %Fig 5 12WPI
% mice=third_timelapse   %Fig 5 16WPI
% mice=cno_cont_mice;    %Fig S3 

  valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;    %%new orfer for valves
    ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
         for gg=1:length( ploting);
     k_for_plot(gg)=find( valves_12==ploting(gg));
     end


if length(mice)==length(cont_mice);
    if sum(mice-cont_mice)==0;
        condition='Control'
    end
end

if length(mice)==length(exp_mice);
    if sum(mice-exp_mice)==0;
        condition='Experiment'
    end
end

if length(mice)==length(awake_mice);
    if sum(mice-awake_mice)==0;
        condition='awake'
    end
end

if length(mice)==length(mice_awake_cont);
    if sum(mice-mice_awake_cont)==0;
        condition='mice_awake_control'
    end
end

if length(mice)==length(first_timelapse);
    if sum(mice-first_timelapse)==0;
        condition='first_timelapse'
    end
end

if length(mice)==length(sec_timelapse);
    if sum(mice-sec_timelapse)==0;
        condition='sec_timelapse'
    end
end

if length(mice)==length(third_timelapse);
    if sum(mice-third_timelapse)==0;
        condition='third_timelapse'
    end
end


if length(mice)==length(cno_cont_mice);
    if sum(mice-cno_cont_mice)==0;
        condition='Wt_cno_cont'
    end
end

num_of_mice=length(mice); 


%for experiment only 
total_neurons_before=0;
total_neurons_after=0;
total_cells_resp_before=0;
total_cells_resp_after=0;
all_responses_magnituds_before=[];
all_responses_magnituds_after=[];
all_significant_responses_magnitudes_before=[];
all_significant_responses_magnitudes_after=[];
all_significant_responses_integrals_before=[];
all_significant_responses_integrals_after=[];
total_prop_of_cells_before=[];
total_prop_of_cells_after=[];
total_prop_of_cells_before_inhib=[];
total_prop_of_cells_after_inhib=[];
total_prop_of_cells_before_excit=[];
total_prop_of_cells_after_excit=[];
%%parameters to change
numOfOdors =12; %including blank
blank=5;
%% load data

for i=1:num_of_mice;

mice_table{1,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results before']); %create a struct- rows for conditions, coloumns for mice
mice_table{2,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results after']); 

end

num_of_mice=length(mice)
for i=1:num_of_mice;

% perm_table{1,i} = load(['perutation results- mouse' num2str(mice(i)) 'results before']); %create a struct- rows for conditions, coloumns for mice
% perm_table{2,i} = load(['perutation results- mouse' num2str(mice(i)) 'results after']); 
perm_table{3,i} = load(['perutation results- mouse' num2str(mice(i)) ' bef_vs_after']); 

end

 %% significant test for between after difference. so now a ttest for 2 samples is conducted on the means throughout the whole response window
 %exp_mice
 
 all_mice_p_values=[];%zeros(2000,numOfOdors);
for i=1:num_of_mice;
    
     p_values_for_dif_mat=perm_table{3,i}.perm_p_vals_bef_vs_aft
    
    
%      p_values_for_dif_mat=[];
%     group_before=mice_table{1,i}.mat_for_diffrence_test;
%      group_after=mice_table{2,i}.mat_for_diffrence_test;
%      for ii=1:size(group_before,1);
%          for kk=1:size(group_before,3);
%              [h,p]=ttest2(group_before(ii,:,kk),group_after(ii,:,kk));
%              p_values_for_dif_mat(ii,kk)=p;
%              
%          end
%      end
     mice_table{1,i}.p_values_for_dif_mat= p_values_for_dif_mat;
     mice_table{2,i}.p_values_for_dif_mat= p_values_for_dif_mat;
     all_mice_p_values= [all_mice_p_values;p_values_for_dif_mat];
end
         thresh=0.05;
         signif_dif_bef_af_all_mice=all_mice_p_values<thresh;
         

%%


num_of_responses_per_mouse=zeros(1,num_of_mice);
all_total_meta_des=[];
all_total_meta_des_changed=[];
all_signif_dif_bef_af_binary=[];
all_mice_sig_bef_both_bef_af_integral=[];
all_mice_sig_af_both_bef_af_integral=[];
all_total_meta_des_bef_af_inh_ex=[];
all_response_magnitude_by_integral_before=[];
all_response_magnitude_by_integral_after=[];
all_responses_magnituds_before=[];
all_responses_magnituds_after=[];

for i=1:num_of_mice;


total_prop_of_cells_before=[total_prop_of_cells_before;(mice_table{1,i}.prop_of_cells)];
total_prop_of_cells_after=[total_prop_of_cells_after;(mice_table{2,i}.prop_of_cells)];

total_prop_of_cells_before_inhib=[total_prop_of_cells_before_inhib;(mice_table{1,i}.prop_of_cells_inhib)];
total_prop_of_cells_after_inhib=[total_prop_of_cells_after_inhib;(mice_table{2,i}.prop_of_cells_inhib)];


total_prop_of_cells_before_excit=[total_prop_of_cells_before_excit;(mice_table{1,i}.prop_of_cells_excit)];
total_prop_of_cells_after_excit=[total_prop_of_cells_after_excit;(mice_table{2,i}.prop_of_cells_excit)];
%total_meta_des=(mice_table{1,i}.meta_desicion>1)+(mice_table{2,i}.meta_desicion>1);%13.9.18

total_meta_des=abs(mice_table{1,i}.meta_desicion~=0)+abs(mice_table{2,i}.meta_desicion~=0);
total_meta_des(:,blank)=0;
meta_bef=mice_table{1,i}.meta_desicion;
meta_af=mice_table{2,i}.meta_desicion;
meta_bef(meta_bef==-1)=-3;
meta_bef(meta_bef==1)=3;
total_meta_des_both_bef_aft_changed=meta_bef.*meta_af;  %this one includes only responses that were sig both bef and after
total_meta_des_bef_or_af_inh_ex=meta_bef+meta_af;
total_meta_des_bef_or_af_inh_ex(:,blank)=0;
total_meta_des_both_bef_aft_changed(:,blank)=0;

all_total_meta_des_bef_af_inh_ex=[all_total_meta_des_bef_af_inh_ex;total_meta_des_bef_or_af_inh_ex];
all_total_meta_des_changed=[all_total_meta_des_changed;total_meta_des_both_bef_aft_changed];
all_total_meta_des=[all_total_meta_des;total_meta_des];
%all_total_meta_des(:,blank)=0; 

significant_responses_before=(mice_table{1,i}.response_magnitude);   %problematic definition 29.1.19
significant_responses_after=(mice_table{2,i}.response_magnitude);

% significant_responses_before=(mice_table{1,i}.response_magnitude_by_integral);
% significant_responses_after=(mice_table{2,i}.response_magnitude_by_integral);


sig_bef_integral=mice_table{1,i}.response_magnitude_by_integral;
sig_af_integral=mice_table{2,i}.response_magnitude_by_integral;

p_values_for_dif_mat= mice_table{1,i}.p_values_for_dif_mat;
signif_dif_bef_af_binary=p_values_for_dif_mat<thresh;
%signif_dif_bef_af_binary(total_meta_des<1)=[];
signif_dif_bef_af_binary(total_meta_des==0)=[]; %13.9.18
% 
% significant_responses_before(total_meta_des<1)=[]; %13.9.18
% significant_responses_after(total_meta_des<1)=[];

significant_responses_before(total_meta_des==0)=[]; %13.9.18
significant_responses_after(total_meta_des==0)=[];

sig_bef_both_bef_af_integral=mice_table{1,i}.response_magnitude_by_integral;
sig_af_both_bef_af_integral=mice_table{2,i}.response_magnitude_by_integral;
sig_bef_both_bef_af_integral(total_meta_des_both_bef_aft_changed==0)=[];
sig_af_both_bef_af_integral(total_meta_des_both_bef_aft_changed==0)=[];
all_mice_sig_bef_both_bef_af_integral=[all_mice_sig_bef_both_bef_af_integral;sig_bef_both_bef_af_integral'];
all_mice_sig_af_both_bef_af_integral=[all_mice_sig_af_both_bef_af_integral;sig_af_both_bef_af_integral'];
sig_bef_integral(total_meta_des==0)=[];
sig_af_integral(total_meta_des==0)=[];



mice_table{1,i}.PAIRED_response_magnitude=significant_responses_before;
mice_table{2,i}.PAIRED_response_magnitude=significant_responses_after;

mice_table{1,i}.PAIRED_response_integral=sig_bef_integral;
mice_table{2,i}.PAIRED_response_integral=sig_af_integral;

num_of_responses_per_mouse(i)=size(significant_responses_before,2);

all_significant_responses_magnitudes_before=[all_significant_responses_magnitudes_before;significant_responses_before'];
all_significant_responses_magnitudes_after=[all_significant_responses_magnitudes_after;significant_responses_after'];

all_significant_responses_integrals_before=[all_significant_responses_integrals_before;sig_bef_integral'];
all_significant_responses_integrals_after=[all_significant_responses_integrals_after;sig_af_integral'];

all_signif_dif_bef_af_binary=[all_signif_dif_bef_af_binary;signif_dif_bef_af_binary'];

total_neurons_before=total_neurons_before+(mice_table{1,i}.numOfNeurons); %count the total number of neurons
total_neurons_after=total_neurons_after+(mice_table{2,i}.numOfNeurons);
total_cells_resp_before=total_cells_resp_before+(mice_table{1,i}.CELLS_RESP);
total_cells_resp_after=total_cells_resp_after+(mice_table{2,i}.CELLS_RESP);
all_responses_magnituds_before=[all_responses_magnituds_before;(mice_table{1,i}.response_magnitude)];
all_responses_magnituds_after=[all_responses_magnituds_after;(mice_table{2,i}.response_magnitude)];

all_response_magnitude_by_integral_before=[all_response_magnitude_by_integral_before;(mice_table{1,i}.response_magnitude_by_integral)];
all_response_magnitude_by_integral_after=[all_response_magnitude_by_integral_after;(mice_table{2,i}.response_magnitude_by_integral)];
% all_significant_responses_magnitudes_before=[all_significant_responses_magnitudes_before;(mice_table{1,i}.significant_responses_magnitudes)];
% all_significant_responses_magnitudes_after=[all_significant_responses_magnitudes_after;(mice_table{2,i}.significant_responses_magnitudes)];
end
all_total_meta_des(:,blank)=0;
all_total_meta_des_bef_af_inh_ex(:,blank)=0; %now use this for a pai chart

a=sum(sum(all_total_meta_des_bef_af_inh_ex==4)) %ex to ex
b=sum(sum(all_total_meta_des_bef_af_inh_ex==3)) %ex to null
c=sum(sum(all_total_meta_des_bef_af_inh_ex==2)) %ex to in

d=sum(sum(all_total_meta_des_bef_af_inh_ex==-4)) %in to in
e=sum(sum(all_total_meta_des_bef_af_inh_ex==-3)) %in to null
f=sum(sum(all_total_meta_des_bef_af_inh_ex==-2)) %in to ex

g=sum(sum(all_total_meta_des_bef_af_inh_ex==-1)) %null to in
h=sum(sum(all_total_meta_des_bef_af_inh_ex==1)) %null to ex


vec_for_pai_chart=[a b c d e f g h]

figure
pie(vec_for_pai_chart)
legend('ex to ex','ex to null','ex to in', 'in to in', 'in to null','in to ex','null to in', 'null to ex')



%exclude the blank

% %cont_mice
% num_of_responses_per_mouse=zeros(1,num_of_cont_mice);
% all_total_meta_des=[];
% all_signif_dif_bef_af_binary=[];
% for i=1:num_of_cont_mice;
% 
% 
% total_prop_of_cells_before=[total_prop_of_cells_before;(mice_table_cont{1,i}.prop_of_cells)];
% total_prop_of_cells_after=[total_prop_of_cells_after;(mice_table_cont{2,i}.prop_of_cells)];
% 
% total_meta_des=(mice_table_cont{1,i}.meta_desicion>1)+(mice_table_cont{2,i}.meta_desicion>1);
% total_meta_des(:,blank)=0;
% all_total_meta_des=[all_total_meta_des;total_meta_des];
% %all_total_meta_des(:,blank)=0;
% 
% significant_responses_before=(mice_table_cont{1,i}.response_magnitude);
% significant_responses_after=(mice_table_cont{2,i}.response_magnitude);
% 
% 
% p_values_for_dif_mat= mice_table_cont{1,i}.p_values_for_dif_mat;
% signif_dif_bef_af_binary=p_values_for_dif_mat<thresh;
% signif_dif_bef_af_binary(total_meta_des<1)=[];
% 
% significant_responses_before(total_meta_des<1)=[];
% significant_responses_after(total_meta_des<1)=[];
% 
% mice_table_cont{1,i}.PAIRED_response_magnitude=significant_responses_before;
% mice_table_cont{2,i}.PAIRED_response_magnitude=significant_responses_after;
% num_of_responses_per_mouse(i)=size(significant_responses_before,2);
% 
% all_significant_responses_magnitudes_before=[all_significant_responses_magnitudes_before;significant_responses_before'];
% all_significant_responses_magnitudes_after=[all_significant_responses_magnitudes_after;significant_responses_after'];
% all_signif_dif_bef_af_binary=[all_signif_dif_bef_af_binary;signif_dif_bef_af_binary'];
% 
% total_neurons_before=total_neurons_before+(mice_table_cont{1,i}.numOfNeurons); %count the total number of neurons
% total_neurons_after=total_neurons_after+(mice_table_cont{2,i}.numOfNeurons);
% total_cells_resp_before=total_cells_resp_before+(mice_table_cont{1,i}.CELLS_RESP);
% total_cells_resp_after=total_cells_resp_after+(mice_table_cont{2,i}.CELLS_RESP);
% all_responses_magnituds_before=[all_responses_magnituds_before;(mice_table_cont{1,i}.response_magnitude)];
% all_responses_magnituds_after=[all_responses_magnituds_after;(mice_table_cont{2,i}.response_magnitude)];
% % all_significant_responses_magnitudes_before=[all_significant_responses_magnitudes_before;(mice_table{1,i}.significant_responses_magnitudes)];
% % all_significant_responses_magnitudes_after=[all_significant_responses_magnitudes_after;(mice_table{2,i}.significant_responses_magnitudes)];
% end



figure
plot(all_mice_sig_bef_both_bef_af_integral,all_mice_sig_af_both_bef_af_integral,'*')
xlim([-0.5 1.5]);
ylim([-0.5 1.5]);
x=-5:5;
y=x;

hold on
title('only sig both bef and after by integral')
plot(x,y,'--')
      %% means

 mean_before_after=zeros(2,num_of_mice);
for i=1:size(mice_table,1);
    for j=1:num_of_mice;
      mean_before_after(i,j)=mean(abs(mice_table{i,j}.PAIRED_response_magnitude));
      %  mean_before_after(i,j)=mean(mice_table{i,j}.significant_responses_magnitudes);
    end
end

mice_deltas=mean_before_after(1,:)-mean_before_after(2,:);
[p,h]=signrank(mean_before_after(1,:),mean_before_after(2,:))
stats.mice_deltas_signrank_p_n=([p,length(mean_before_after(1,:))]);

%   
%     for j=1:num_of_exp_mice;
%     %   deltas_before_after_exp(j)=(mice_table_exp{1,j}.PAIRED_response_magnitude)-(mice_table_exp{2,j}.PAIRED_response_magnitude)
%      a=(mice_table_exp{1,j}.PAIRED_response_magnitude)-(mice_table_exp{2,j}.PAIRED_response_magnitude);
%        %  mean_before_after(i,j)=mean(mice_table{i,j}.significant_responses_magnitudes);
%      exp_vs_cont{1,j}.paired_deltas=a;
%        
%     end

figure 
plot( mean_before_after)
legend('M1', 'M2', 'M3', 'M4', 'M5','M6', 'M7', 'M8', 'M9','M10')
ylabel('mean response magnitude-before and after')
xlabel(['before' '                                                           '  'after'])
title('Mean response absolute magnitude')

prop_mean_before_after=mean_before_after;
for m=1:num_of_mice;
prop_mean_before_after(1,m)=  mean_before_after(1,m)/mean_before_after(1,m);
prop_mean_before_after(2,m)=  mean_before_after(2,m)/mean_before_after(1,m);
end

prop_for_comparison=prop_mean_before_after(2,:);

figure
plot( prop_mean_before_after,'linewidth',2)
legend('M1', 'M2', 'M3', 'M4', 'M5','M6','M7','M8','M9','M10');% Orientation,'horizontal')
%legend('Orientation','horizontal');
if numOfOdors==12;
legend('Location','southwest');
legend('boxoff');
title('Mean response absolute magnitude')

end

if numOfOdors==8;
legend('Location','northwest');
legend('boxoff');
end
ylim([0 2]);
ylabel('Normalized mean response magnitude')
xlabel(['before' '                                                           '  'after'])
 ax=gca;
 ax.FontSize=14;
ax.XTick=[];
   
  %% cdfs
  
prop_before=total_cells_resp_before/total_neurons_before;
prop_after=total_cells_resp_after/total_neurons_after;
% 

    for ii=1:numOfOdors; %so num of odors include blank so it feets
       CDF_BEFORE(ii)=sum(prop_before(1:ii));
       CDF_AFTER(ii)=sum(prop_after(1:ii));
       
       CDF_DELTAS(ii)=CDF_BEFORE(ii)-CDF_AFTER(ii)
    end
    
  
    
     deltas_for_ks=(total_prop_of_cells_before-total_prop_of_cells_after);
     deltas_for_ks_norm=deltas_for_ks./std(deltas_for_ks);
  [h,p,ksstat,cv]=kstest(deltas_for_ks_norm)
    %stats for paper, use the first line here:
    [h,p,ks2stat]=kstest2(total_prop_of_cells_before,total_prop_of_cells_after) %for the first time only
    stats.cdf_ks_test_2_tailed_p_ksstat_n=[p,ks2stat,length(total_prop_of_cells_before)]
    [h,p,ks2stat]=kstest2(total_prop_of_cells_before,total_prop_of_cells_after,'Tail','smaller') % for 2 and 3rd tl
        stats.cdf_ks_test_1_tail_p_ksstat_n=[p,ks2stat,length(total_prop_of_cells_before)]

      
  [h,p,ks2stat]=kstest(total_prop_of_cells_before-total_prop_of_cells_after)
    
    x=(0:numOfOdors-1);
    figure; %figure5
     set(gcf, 'Position', [350, 350, 350, 350])
    plot(fliplr(x),CDF_BEFORE,'b','LineWidth',3);
    set(gca,'XDir','reverse');
    hold on
    plot(fliplr(x),CDF_AFTER, 'r','LineWidth',3)
    legend(['Before-' num2str(total_neurons_before) 'cells'], ['After-' num2str(total_neurons_after) 'cells'],'FontSize',12)
    legend('Location','northwest');
    legend('boxoff');
    title(['CDF for odor responses- time lapse ,10  mice.  p=' num2str(p)]);
    %title('CDF for odor responses- the whole population of cells,control mouse')   %for control
    ylabel('proportion of responding cells');
    xlabel('odors');
    xlim([0 numOfOdors-1]);
    ylim([0 1]);
    ax = gca;
    ax.FontSize = 14;
    box off

    
    figure; %figure5
     set(gcf, 'Position', [350, 350, 350, 350])
    plot(fliplr(x),CDF_DELTAS,'K','LineWidth',3);
    set(gca,'XDir','reverse');
    hold on
   
    legend(['Before-' num2str(total_neurons_before) 'cells'], ['After-' num2str(total_neurons_after) 'cells'],'FontSize',12)
    legend('Location','northwest');
    legend('boxoff');
    title(['CDF DELTAS- time lapse ,10  mice.  p=' num2str(p)]);
    %title('CDF for odor responses- the whole population of cells,control mouse')   %for control
    ylabel('proportion of responding cells');
    xlabel('odors');
    xlim([0 numOfOdors-1]);
    ylim([-0.1 0.2]);
    ax = gca;
    ax.FontSize = 14;
    box off
    

    %% plot the means of the cdfs
figure
   
hold on
x=(1:2)
y=[mean(total_prop_of_cells_before),mean(total_prop_of_cells_after)];
sem_bef=std(total_prop_of_cells_before)/sqrt(size(total_prop_of_cells_before,1));
sem_aft=std(total_prop_of_cells_after)/sqrt(size(total_prop_of_cells_after,1));
sems=[sem_bef sem_aft]
bar(x(1),y(1),'b')
bar(x(2),y(2),'r')
errorbar(x,y,sems,'.')
%errorbar(x(2),sems(2),'.')
% [h,p]=signrank(total_prop_of_cells_before,total_prop_of_cells_after)
[h,p]=ttest(total_prop_of_cells_before,total_prop_of_cells_after);
stats.number_of_resp_after_vs_before_p_n=[p,length(total_prop_of_cells_after)];

ylim([0 5])      
title('Average number of responses per cell')

%deltas
deltas=(total_prop_of_cells_before-total_prop_of_cells_after)
sem_deltas=std(deltas)/sqrt(size(deltas,1));
sems=[sem_deltas sem_deltas]   %plot it all twicw to replace manually with cont/exp
y=[mean(deltas),mean(deltas)]
figure
hold on
bar(x(1),y(1),'k')
bar(x(2),y(2),'k')
errorbar(x,y,sems,'.')

   

figure
hold on
x=(1:2)
y=[mean(total_prop_of_cells_before_inhib),mean(total_prop_of_cells_after_inhib)];
sem_bef=std(total_prop_of_cells_before_inhib)/sqrt(size(total_prop_of_cells_before_inhib,1));
sem_aft=std(total_prop_of_cells_after_inhib)/sqrt(size(total_prop_of_cells_after_inhib,1));
sems=[sem_bef sem_aft]
bar(x(1),y(1),'b')
bar(x(2),y(2),'r')
errorbar(x,y,sems,'.')
%errorbar(x(2),sems(2),'.')
ranksum(total_prop_of_cells_before_inhib,total_prop_of_cells_after_inhib)
[h,p]=ttest(total_prop_of_cells_before_inhib,total_prop_of_cells_after_inhib)
ylim([0 5])      
title('Average number of inhibitory responses per cell')


figure
   
hold on
x=(1:2)
y=[mean(total_prop_of_cells_before_excit),mean(total_prop_of_cells_after_excit)];
sem_bef=std(total_prop_of_cells_before_excit)/sqrt(size(total_prop_of_cells_before_excit,1));
sem_aft=std(total_prop_of_cells_after_excit)/sqrt(size(total_prop_of_cells_after_excit,1));
sems=[sem_bef sem_aft]
bar(x(1),y(1),'b')
bar(x(2),y(2),'r')
errorbar(x,y,sems,'.')
%errorbar(x(2),sems(2),'.')
ranksum(total_prop_of_cells_before_excit,total_prop_of_cells_after_excit)
[h,p]=ttest(total_prop_of_cells_before_excit,total_prop_of_cells_after_excit)
ylim([0 5])      
title('Average number of excitatory responses per cell')
%  figure
%  hist(all_responses_magnituds_before)
%  %title(['all responses magnitudes before-' num2str(total_neurons_before) 'cells' ',three experiment mice' ])
%  title(['all responses magnitudes before-' num2str(total_neurons_before) 'cells' ',control mouse' ]) %for control
%  xlabel('Response Magnitude')
%  ylabel('Number of cells')
%  legend('Blank', 'Isoamyl acetate', 'Ethyl butyrate','Amyl Acetate','2-phenylethanol','Ethyl tiglate', 'Alpha Pinen','Geraniol') 
% 
%  figure
%  hist(all_responses_magnituds_after)
%  %title(['all responses magnitudes after-' num2str(total_neurons_after) 'cells' ',three experiment mice' ])
%  title(['all responses magnitudes after-' num2str(total_neurons_after) 'cells' ',control mouse' ]) %for control
%  xlabel('Response Magnitude')
%  ylabel('Number of cells')
%  legend('Blank', 'Isoamyl acetate', 'Ethyl butyrate','Amyl Acetate','2-phenylethanol','Ethyl tiglate', 'Alpha Pinen','Geraniol') 

%% boxplots   %at30.1.19 changed all to abs
% mean_before=mean(all_significant_responses_magnitudes_before);    %for peacks
% SD_before=std(all_significant_responses_magnitudes_before);
% mean_after=mean(all_significant_responses_magnitudes_after);
% SD_after=std(all_significant_responses_magnitudes_after);

mean_before=mean(abs(all_significant_responses_integrals_before));    %for integrals
SD_before=std(abs(all_significant_responses_integrals_before));
mean_after=mean(abs(all_significant_responses_integrals_after));
SD_after=std(abs(all_significant_responses_integrals_after));


resp_vec_before=[];
resp_vec_after=[];
resp_vec_before_integral=[];
resp_vec_after_integral=[];

for j=1:size(all_responses_magnituds_before,1)*size(all_responses_magnituds_before,2);
    resp_vec_before(j)=abs(all_responses_magnituds_before(j));
    resp_vec_after(j)=abs(all_responses_magnituds_after(j));
    
    resp_vec_before_integral(j)=abs(all_response_magnitude_by_integral_before(j));
    resp_vec_after_integral(j)=abs(all_response_magnitude_by_integral_after(j));
    
end
[h,p] = signrank( resp_vec_before, resp_vec_after)
p_for_all=p;
[h,p] = signrank( resp_vec_before_integral,resp_vec_after_integral)

[h,p] = ttest(all_significant_responses_magnitudes_after,all_significant_responses_magnitudes_before)
[psr,h]=signrank(all_significant_responses_magnitudes_after,all_significant_responses_magnitudes_before)
[h,p] = ttest( resp_vec_before, resp_vec_after)

[h_int,p_int] = signrank( resp_vec_before_integral, resp_vec_after_integral) %all responses- not only significant
p_for_all=p_int;
[h_int,p_int] = ttest(all_significant_responses_integrals_before,all_significant_responses_integrals_after)
% [psr_int,h_int]=signrank(all_significant_responses_integrals_before,all_significant_responses_integrals_after)
% [h_int,p_int] = ttest( resp_vec_before_integral, resp_vec_after_integral)

%% boxplots

figure
%x=[all_significant_responses_integrals_before;all_significant_responses_integrals_after];
x=[abs(all_significant_responses_magnitudes_before);abs(all_significant_responses_magnitudes_after)];
x2=[all_significant_responses_magnitudes_before;all_significant_responses_magnitudes_after];
%deltas_for_boxplots=[all_significant_responses_magnitudes_before-all_significant_responses_magnitudes_after];
%by peack
deltas_for_boxplots=[all_significant_responses_integrals_before-all_significant_responses_integrals_after]; %by integral

abs_deltas_for_boxplots=[abs(all_significant_responses_integrals_before)-abs(all_significant_responses_integrals_after)]; %added at 11.1.19


g = [ones(size(all_significant_responses_magnitudes_before)); 2*ones(size(all_significant_responses_magnitudes_after))];
% boxplot(x,g)
boxplot(x,g,'symbol','');
%for peack:
%title(['Before and After ' num2str(length(all_significant_responses_magnitudes_before)) ' significant responses ' ' Before Mean=' num2str(mean_before) ' SD=' num2str(SD_before) ' After Mean=' num2str(mean_after) ' SD=' num2str(SD_after) '    p=' num2str(p)],'FontSize',24);% '    p2=' num2str(p_for_all)]);
%for integral:
title(['Before and After abs ' num2str(length(all_significant_responses_magnitudes_before)) ' significant responses ' ' Before Mean=' num2str(mean_before) ' SD=' num2str(SD_before) ' After Mean=' num2str(mean_after) ' SD=' num2str(SD_after) '    p=' num2str(p)],'FontSize',24);% '    p2=' num2str(p_for_all)]);

ylabel('Response magnitude');
xlabel('Before                                                       after');
%ylim([0 0.6])
ylim([0 1]);
a=gca;
a.FontSize=20;

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

figure
boxplot(x2,g,'symbol','');
title(['Before and After' num2str(length(all_significant_responses_magnitudes_before)) ' significant responses ' ' Before Mean=' num2str(mean_before) ' SD=' num2str(SD_before) ' After Mean=' num2str(mean_after) ' SD=' num2str(SD_after) '    p=' num2str(p)],'FontSize',24);% '    p2=' num2str(p_for_all)]);

ylabel('Response magnitude');
xlabel('Before                                                       after');
%ylim([0 0.6])
ylim([-0.5 1]);
a=gca;
a.FontSize=20;
box off

%% histograms
 figure
 
 vect=(-0.8:0.02:1.5);
%  all_significant_responses_magnitudes_before_plot=all_significant_responses_magnitudes_before;
%  all_significant_responses_magnitudes_after_plot=all_significant_responses_magnitudes_after;
 
hist(all_significant_responses_magnitudes_before,vect,'linewidth',3) ;
title(['all significant' num2str(length(all_significant_responses_magnitudes_before)) ' responses peaks']);
xlabel('Response Magnitude - dF/F','FontSize',14);
ylabel('Number of Responses','FontSize',14);
%ylim([0 160])
xlim([-0.8 1.5]);
% set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
set(get(gca,'child'),'FaceColor','b','EdgeColor','b','linewidth',3);

%ylim([0 100]) %for control
%xlim([0 1.2]) %for control

hold on
hist(all_significant_responses_magnitudes_after,vect,'EdgeColor','r','FaceColor','r','linewidth',3);
set(get(gca,'child'),'linewidth',3);
xlim([-0.8 1.5]);
%legend('Before', 'After')
%title(['After- all significant' num2str(length(all_significant_responses_magnitudes_after)) 'responses' 'Mean=' num2str(mean_after) 'SD=' num2str(SD_after) ])
%xlabel('Response Magnitude')
%ylabel('Number of cells')

%ylim([0 100]) %for control
%xlim([0 1.2]) %for control
ax = gca;
ax.FontSize = 14;



%% histograms for abs values

 figure
 
 vect=(-0.5:0.02:2);
%  all_significant_responses_magnitudes_before_plot=all_significant_responses_magnitudes_before;
%  all_significant_responses_magnitudes_after_plot=all_significant_responses_magnitudes_after;
 
hist(abs(all_significant_responses_magnitudes_before),vect,'linewidth',3) ;
title(['all significant' num2str(length(all_significant_responses_magnitudes_before)) ' responses. ' 'peaks absolute values']);
xlabel('Response Magnitude - dF/F','FontSize',14);
ylabel('Number of Responses','FontSize',14);
%ylim([0 160])
xlim([-0.5 2]);
% set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
set(get(gca,'child'),'FaceColor','b','EdgeColor','b','linewidth',3);

%ylim([0 100]) %for control
%xlim([0 1.2]) %for control

hold on
hist(abs(all_significant_responses_magnitudes_after),vect,'EdgeColor','r','FaceColor','r','linewidth',3);
set(get(gca,'child'),'linewidth',3);
xlim([0 2]);

%legend('Before', 'After')
%title(['After- all significant' num2str(length(all_significant_responses_magnitudes_after)) 'responses' 'Mean=' num2str(mean_after) 'SD=' num2str(SD_after) ])
%xlabel('Response Magnitude')
%ylabel('Number of cells')

%ylim([0 100]) %for control
%xlim([0 1.2]) %for control
ax = gca;
ax.FontSize = 14;



%% histograms ny integral
 figure
 
 vect=(-0.7:0.02:0.7);
%  all_significant_responses_magnitudes_before_plot=all_significant_responses_magnitudes_before;
%  all_significant_responses_magnitudes_after_plot=all_significant_responses_magnitudes_after;
 
hist(all_significant_responses_integrals_before,vect,'linewidth',3); 
title(['all significant' num2str(length(all_significant_responses_integrals_before)) ' responses. integrals']);
xlabel('Response Magnitude - dF/F','FontSize',14);
ylabel('Number of Responses','FontSize',14);
%ylim([0 160])
xlim([-0.7 0.7])
% set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
set(get(gca,'child'),'FaceColor','none','EdgeColor','b','linewidth',3);



%ylim([0 100]) %for control
%xlim([0 1.2]) %for control

hold on
hist(all_significant_responses_integrals_after,vect,'EdgeColor','r','linewidth',3);
set(get(gca,'child'),'FaceColor','none','linewidth',3);
xlim([-0.7 0.7]);
legend('Before', 'After');
%title(['After- all significant' num2str(length(all_significant_responses_magnitudes_after)) 'responses' 'Mean=' num2str(mean_after) 'SD=' num2str(SD_after) ])
%xlabel('Response Magnitude')
%ylabel('Number of cells')

%ylim([0 100]) %for control
%xlim([0 1.2]) %for control
ax = gca;
ax.FontSize = 14;

%%
%% histograms by integral absolute values
 figure
 
 vect=(-0.5:0.02:2);
%  all_significant_responses_magnitudes_before_plot=all_significant_responses_magnitudes_before;
%  all_significant_responses_magnitudes_after_plot=all_significant_responses_magnitudes_after;
plot_bef=all_significant_responses_integrals_before;
plot_af=all_significant_responses_integrals_after;
plot_bef(abs(plot_bef)>=1)=1;
plot_af(abs(plot_af)>=1)=1;
 
hist(abs(plot_bef),vect,'linewidth',3); 
title(['all significant' num2str(length(all_significant_responses_integrals_before)) ' responses. integrals']);
xlabel('Response Magnitude - dF/F','FontSize',14);
ylabel('Number of Responses','FontSize',14);
%ylim([0 160])
xlim([-0.5 2])
% set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
set(get(gca,'child'),'FaceColor','none','EdgeColor','b','linewidth',3);



%ylim([0 100]) %for control
%xlim([0 1.2]) %for control

hold on
hist(abs(plot_af),vect,'EdgeColor','r','linewidth',3);
set(get(gca,'child'),'FaceColor','none','linewidth',3);
xlim([0 1]);
legend('Before', 'After');
%title(['After- all significant' num2str(length(all_significant_responses_magnitudes_after)) 'responses' 'Mean=' num2str(mean_after) 'SD=' num2str(SD_after) ])
%xlabel('Response Magnitude')
%ylabel('Number of cells')

%ylim([0 100]) %for control
%xlim([0 1.2]) %for control
ax = gca;
ax.FontSize = 14;


%% dots plot


x=(-1:4);
y=(-1:4);
%all_significant_responses_magnitudes_before_exp=all_significant_responses_magnitudes_before;
%all_significant_responses_magnitudes_after_exp=all_significant_responses_magnitudes_after;
figure
% plot( resp_vec_before, resp_vec_after, '*')
all_red_black=all_significant_responses_magnitudes_before-all_significant_responses_magnitudes_after;

a=all_significant_responses_magnitudes_before(all_signif_dif_bef_af_binary==0);
b=all_significant_responses_magnitudes_after(all_signif_dif_bef_af_binary==0);

c=all_significant_responses_magnitudes_before(all_signif_dif_bef_af_binary==1);
d=all_significant_responses_magnitudes_after(all_signif_dif_bef_af_binary==1);
e=c-d;
f=e(e>=0);
g=e(e<=0);

blacks=a-b;

up=length(blacks(blacks<=0));
down=length(blacks(blacks>0));
up_percent=up/(up+down)*100;
down_percent=down/(up+down)*100;
signif_up=length(g) ;
signif_up_percent= length(g)/(length(g)+length(f))*100;
signif_down=length(f);
signif_down_percent= length(f)/(length(g)+length(f))*100;

stats_for_chi=[length(a)+length(c) length(a) length(f) length(g)]; %first all then all non signif then those that are strongerbefore then stronger after
% plot(all_significant_responses_magnitudes_before, all_significant_responses_magnitudes_after, 'r *')
plot(a,b,'k o');
hold on
plot(c,d,'r o');
hold on

%plot(all_significant_responses_magnitudes_before_control, all_significant_responses_magnitudes_after_control,'b *')
legend([num2str(length(all_significant_responses_magnitudes_before)), ' responses'], [num2str(length(c)) ' significant differences'])
legend('Location','northwest');
set(gcf,'color','w');
  xlim([-0.8 3.5]) ;
  ylim([-0.8 3.5]); 
hold on
plot(x,y,'k --');

xlabel('Response magnitude before (dF/F)');
ylabel('Response magnitude after (dF/F)');
title('Response magnitude by peak')
axis square;
box off

%for combined plot:
%all_significant_responses_magnitudes_before_exp=all_significant_responses_magnitudes_before;
%all_significant_responses_magnitudes_after_exp=all_significant_responses_magnitudes_after;

ax = gca;
ax.FontSize = 14;


%% only sig both before and after

all_total_meta_des(:,blank)=0;

figure
plot(all_mice_sig_bef_both_bef_af_integral,all_mice_sig_af_both_bef_af_integral,'*');
xlim([-0.5 1.5]);
ylim([-0.5 1.5]);
x=-5:5;
y=x;

hold on
plot(x,y,'--')

%% by integral
%all paired responses, sorted in unpaired manner


%all_significant_responses_magnitudes_before_exp=all_significant_responses_magnitudes_before;
%all_significant_responses_magnitudes_after_exp=all_significant_responses_magnitudes_after;

% plot( resp_vec_before, resp_vec_after, '*')
all_red_black=all_significant_responses_integrals_before-all_significant_responses_integrals_after; 
bef_sig_sort=sort(all_significant_responses_integrals_before);
af_sig_sort=sort(all_significant_responses_integrals_after);
x=1:length(all_significant_responses_integrals_before);
figure
set(gcf, 'Position', [150, 150, 550, 150])
plot(x,bef_sig_sort,'b','linewidth',1)
hold on
title('all paired responses- both before and after')
plot(x,af_sig_sort,'r','linewidth',1)
ylim([-0.5 3.5]);
ylim([-0.5 3.5]);

ylim([-0.7 1]);
ylim([-0.7 1]);

xlim([1 length(all_significant_responses_integrals_before)])
xlabel('Response number')
ylabel('response magnitude')

%%
% dots plot by integral

figure

x=(-1:4);
y=(-1:4);
a=all_significant_responses_integrals_before(all_signif_dif_bef_af_binary==0);
b=all_significant_responses_integrals_after(all_signif_dif_bef_af_binary==0);


c=all_significant_responses_integrals_before(all_signif_dif_bef_af_binary==1);
d=all_significant_responses_integrals_after(all_signif_dif_bef_af_binary==1);

% count all inhibitory responses. only for both before and after inhibitory responses
c_inhib_sig=c(c<0);
c_inhib_af_sig=d(c<0);
all_inhib_sig=[c_inhib_sig,c_inhib_af_sig];


for i=1:length(c_inhib_af_sig);
    if all_inhib_sig(i,2)>0;
        all_inhib_sig(i,:)=7777;
    end
end

temp=all_inhib_sig(:,1);
temp2=all_inhib_sig(:,2);
temp(temp==7777)=[];
temp2(temp2==7777)=[];
all_inhib_sig=[temp, temp2];
dec_inhib_sig=sum(abs(temp)>=abs(temp2));
inc_inhib_sig=sum(abs(temp)<=abs(temp2));

[h,p]=ttest(all_inhib_sig(:,1),all_inhib_sig(:,2));

a_inhib_non_sig=a(a<0);
a_inhib_af_non_sig=b(a<0);
all_inhib_non_sig=[a_inhib_non_sig,a_inhib_af_non_sig];


for i=1:length(a_inhib_af_non_sig);
    if all_inhib_non_sig(i,2)>0;
        all_inhib_non_sig(i,:)=7777;
    end
end

temp=all_inhib_non_sig(:,1);
temp2=all_inhib_non_sig(:,2);
temp(temp==7777)=[];
temp2(temp2==7777)=[];
all_inhib_non_sig=[temp, temp2];
dec_inhib_non_sig=sum(abs(temp)>=abs(temp2));
inc_inhib_non_sig=sum(abs(temp)<=abs(temp2));

dec_inhib_percent=dec_inhib_non_sig/(dec_inhib_non_sig+inc_inhib_non_sig)*100;
inc_inhib_percent=inc_inhib_non_sig/(dec_inhib_non_sig+inc_inhib_non_sig)*100;

[h,p]=ttest(all_inhib_non_sig(:,1),all_inhib_non_sig(:,2))

all_inhib_bef=[a_inhib_non_sig;all_inhib_sig(:,1)];

all_inhib_aft=[a_inhib_af_non_sig;all_inhib_sig(:,2)];

[h,p]=ttest(all_inhib_bef,all_inhib_aft);


%signif
e=c-d;
f=e(e>=0);
g=e(e<=0);

%nonsignif
ee=a-b;
ff=ee(ee>=0);
gg=ee(ee<=0);

blacks=a-b;
blacks_pos=(a(a>0));

up=length(blacks(blacks<=0));
down=length(blacks(blacks>0));
up_percent=up/(up+down)*100;
down_percent=down/(up+down)*100;
signif_up=length(g);
signif_up_percent= length(g)/(length(g)+length(f))*100;
signif_down=length(f);
signif_down_percent= length(f)/(length(g)+length(f))*100;


non_signif_up=length(gg) ;
non_signif_up_percent= length(gg)/(length(gg)+length(ff))*100;
non_signif_down=length(ff);
non_signif_down_percent= length(ff)/(length(gg)+length(ff))*100;

stats_for_chi=[length(a)+length(c) length(a) length(f) length(g)]; %first all then all non signif then those that are strongerbefore then stronger after
% plot(all_significant_responses_magnitudes_before, all_significant_responses_magnitudes_after, 'r *')
plot(a,b,'k o')
hold on
plot(c,d,'r o')
hold on

%plot(all_significant_responses_magnitudes_before_control, all_significant_responses_magnitudes_after_control,'b *')
%legend([num2str(length(all_significant_responses_magnitudes_before)), ' responses'], [num2str(length(c)) ' significant differences'])
% legend('Location','northwest');
% legend('boxoff')
set(gcf,'color','w');
  xlim([-0.5 1.5]) 
 %  xylim([-0.5 1])%for awake
%   xlim([-0.5 0]) 
 ylim([-0.5 1.5]) 
 % ylim([-0.5 1]) for awake
 %  ylim([-0.5 0])
hold on
plot(x,y,'k --')
xx=[0 0]
yy=[-5 5]
plot(xx,yy,'k --')
xxx=[-5 5]
yyy=[0 0]
plot(xxx,yyy,'k --')

%title('by integral')
xlabel('Response magnitude before (dF/F)')
ylabel('Response magnitude after (dF/F)')
axis square
box off

%for combined plot:
%all_significant_responses_magnitudes_before_exp=all_significant_responses_magnitudes_before;
%all_significant_responses_magnitudes_after_exp=all_significant_responses_magnitudes_after;

ax = gca;
ax.FontSize = 14;

%deltas, naive calculation
%dist_from_diag=a-b;  %for all- black
dist_from_diag=c-d;   %for significant differences only- red 

sem_deltas=(std(dist_from_diag))/sqrt(length(dist_from_diag));

x=(1:2)
y=[mean(dist_from_diag),mean(dist_from_diag)];
figure
hold on
sems=[sem_deltas sem_deltas]
bar(x(1),y(1),'b')
bar(x(2),y(2),'r')
errorbar(x,y,sems,'.')

title('dist from diag plotted twice just to compare later')

red_percent=length(c)/(length(c)+length(a))*100;

%%
%close all

% delta_exp=zeros(size(all_responses_magnituds_after));
% delta_exp=abs(all_responses_magnituds_after-all_responses_magnituds_before);
% delta_exp=sum(delta_exp,2);
% delta_exp=delta_exp.*0.125;
% figure
% hist(delta_exp);
% title('experiment- average response change per cell. all responses included')
% mean(delta_exp);
% sem_exp=std((delta_exp)/sqrt(length(delta_exp)));
% figure
% orbar(mean(delta_exp),sem_exp);
% 
% delta_exp_sig=abs(all_significant_responses_magnitudes_after-all_significant_responses_magnitudes_before);
% figure
% hist(delta_exp_sig);
% title('experiment- distribution of deltas for cell-odor pairs with significant responses only')
% mean_exp_sig= mean(delta_exp_sig)
% sem_exp_sig=std((delta_exp_sig)/sqrt(length(delta_exp_sig)));
% figure
% errorbar(mean(delta_exp_sig),sem_exp_sig);

%% DELTAS by abs
%deltas=abs(all_significant_responses_magnitudes_before)-abs(all_significant_responses_magnitudes_after); %for peak
deltas=abs(all_significant_responses_integrals_before)-abs(all_significant_responses_integrals_after);   %for integral
mean_deltas=mean(deltas);
std_deltas=std(deltas);
deltas_for_plot=deltas;
limit_for_plot=[-0.5 0.5];
[h,p]=ttest(deltas);
for i=1:length(deltas_for_plot);
if deltas_for_plot(i)>=limit_for_plot(2)
    deltas_for_plot(i)=limit_for_plot(2);
end
if deltas_for_plot(i)<=limit_for_plot(1)
    deltas_for_plot(i)=limit_for_plot(1);
end
end   
figure
set(gcf, 'Position', [400, 400, 700, 400])
plot(deltas_for_plot,'k *')
ylim([limit_for_plot(1) limit_for_plot(2)])
color_vec={'blue', 'green', 'red', 'cyan', 'magenta', 'yellow'};
color_vec={'blue','red', 'blue', 'red','blue', 'red', 'blue', 'red'};
zer=zeros(length(deltas));

hold on
for i=1:num_of_mice;
    f=mod(i,2);
if f==0;
    d=[0.5 0.5 0.5];
else
d=[1 1 1];
end

%patch([sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i)] , [-2 -2 2 2]  ,char(color_vec(i)) , 'facealpha' , 0.2 , 'edgecolor' , 'none')
patch([sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i)] , [-2 -2 2 2]  ,d , 'facealpha' , 0.2 , 'edgecolor' , 'none')
end
xlim([0 length(deltas)]);
plot(zer,'r --', 'linewidth',1);
xlabel('Response Number','FontSize',16);
ylabel('Before-After (dF/F)','FontSize',16);
title('Deltas: responses before-Responses after','FontSize',24);
title(['Deltas for all significant' num2str(length(all_significant_responses_magnitudes_before)) ' responses. ' ' Deltas Mean=' num2str(mean_deltas) ' , STD='  num2str(std_deltas) ', p=' num2str(p)],'FontSize',20)
[p,h]=signrank(deltas);
ax = gca;
ax.FontSize = 14;
%ax.Color='r';
s=deltas<=0;
sum(s);

%% Ratios with abs
%for integrals only:
resp_bef=all_significant_responses_integrals_before;
resp_aft=all_significant_responses_integrals_after;

%for peak:
% resp_bef=all_significant_responses_magnitudes_before;
% resp_aft=all_significant_responses_magnitudes_after;

ratios=abs(resp_bef)./abs(resp_aft);
%ratios=abs(all_significant_responses_magnitudes_after)./abs(all_significant_responses_magnitudes_before);
ratios_forplot=ratios;
figure
%set(gcf, 'Position', [300, 300, 500, 300])
set(gcf, 'Position', [400, 400, 700, 400])
x=1:length(ratios);
y=x;
y(:)=1;

a=0.1;
a_op=10;
ex_ratios=[];
in_ratios=[];
hh=0;
for h=1:length(ratios)
if resp_bef(h)>0;
    if resp_aft(h)>0;
        hh=hh+1;
        if ratios(h)>=a_op;
            ratios(h)=a_op;
        end
        if ratios(h)<=a;
            ratios(h)=a;
        end
        plot(x(h),ratios(h),'k o','markersize',sqrt(mean([abs(resp_bef(h)),abs(resp_aft(h))]))*5,'linewidth',2);%); 
        ex_ratios(hh)=ratios(h);
        %plot(x(h),ratios(h),'k o','markersize',4,'linewidth',2);%); 
        %dots are being magnified sublinearly by sqrt
        %plot(x(h),ratios(h),'k *','markersize',4)%all_significant_responses_magnitudes_before(h)*5)
hold on
    end
end
end
hh=0;
for h=1:length(ratios);
if resp_bef(h)<0;
    if resp_aft(h)<0;
         hh=hh+1;
         if ratios_forplot(h)<=a;
            ratios_forplot(h)=a;
         end
        if ratios_forplot(h)>=a_op;
            ratios_forplot(h)=a_op;
        end
plot(x(h),ratios_forplot(h),'b o','markersize',sqrt(mean([abs(resp_bef(h)),abs(resp_aft(h))]))*10,'linewidth',2);%); 
in_ratios(hh)=ratios(h);
%plot(x(h),ratios_forplot(h),'b o','markersize',4,'linewidth',2);%); 
hold on
    end
end
end

ex_up=sum(ex_ratios>=1);
ex_down=sum(ex_ratios<=1);
in_up=sum(in_ratios>=1);
in_down=sum(in_ratios<=1);

ex_up_percent=ex_up/(ex_up+ex_down)*100
ex_down_percent=ex_down/(ex_up+ex_down)*100

in_up_percent=in_up/(in_up+in_down)*100
in_down_percent=in_down/(in_up+in_down)*100

% ylim([0 5.1])
% color_vec={'blue', 'green', 'red', 'cyan', 'magenta', 'yellow'};
% color_vec={'blue','red', 'blue', 'red','blue', 'red', 'blue', 'red'};

%this is for the patch, for backgrounds in different colors for different
%mice
d=[0.5 0.5 0.5];
e=[1 1 1];

   
zer=zeros(length(ratios_forplot));

hold on
for i=1:num_of_mice;
f=mod(i,2);
if f==0;
    d=[0.5 0.5 0.5];
else
d=[1 1 1];
end
%patch([sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i)] , [a a a_op a_op]  ,char(color_vec(i)) , 'facealpha' , 0.2 , 'edgecolor' , 'none')
patch([sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i)) sum(num_of_responses_per_mouse(1:i))-num_of_responses_per_mouse(i)] , [a a a_op a_op]  ,d , 'facealpha' , 0.2 , 'edgecolor' , 'none')
end
xlim([0 length(ratios_forplot)]);
plot(zer);
xlabel('Number of Response','FontSize',16);
ylabel('Before-After (dF/F)','FontSize',16);
title('RATIOS-Responses before/responses after','FontSize',24);
%title('Responses after/responses before','FontSize',24)
%legend('Excitatory','Inhibitory')

hold on
plot(x,y,'k --','linewidth',1)

ax = gca;
ax.YScale = 'log';
ax.FontSize = 14;
%ax.Color='r';


%% hist for deltas
figure   % commented below is an old version to do everything not normalized equally 
% set(get(gca,'child'),'FaceColor','none','EdgeColor','r','linewidth',3);
% a=hist(deltas);
% b=a./sum(a);
% %hist (,length(deltas)/10,'b')
%xlim([-1.5 1.5])
set(gcf, 'Position', [300, 300, 300, 300])

xx=[0 0 0];
yy=[-10 0 10];

deltas_plot=deltas;   
limits=[-0.5 0.5]%was already changed to show abs values
% deltas_plot(deltas_plot>=1.5)=1.5; %was good for peaks
% deltas_plot(deltas_plot<-1.5)=-1.5;
figure
xlim(limits)

deltas_plot(deltas_plot>=limits(2))=limits(2);
deltas_plot(deltas_plot<limits(1))=limits(1);  %for integrals


[h,p]=ztest(deltas,0,std(deltas)); 
[p,h]=signrank(deltas); 

histogram(deltas_plot,'EdgeColor','b','BinWidth',0.05); 
hold on;
%ylim([0 300]);


ylabel('Number of responses','FontSize',16);
xlabel('Before-After (dF/F)','FontSize',16);
title('Deltas: before-after','FontSize',24);
ax = gca;
ax.FontSize = 14;

figure
set(gcf, 'Position', [300, 300, 300, 300])
histogram(deltas_plot,'Normalization','probability','EdgeColor','b','BinWidth',0.05);
hold on
ylim([0 0.25]);
xlim(limits);
ylabel('Probability','FontSize',16);
xlabel('Before-After (dF/F)','FontSize',16);
title('Deltas: before-after','FontSize',24);
ax = gca;
ax.FontSize = 14;
box off
plot(xx,yy,'k--','linewidth',2)
%

% 
% hold on 
% histogram(b,'Normalization','probability','EdgeColor','r','BinWidth',0.001) 

%% check spontanous activity   good for same number bef and after. otherwise need to change
total_raw_f_mat_bef=[];
total_raw_f_mat_aft=[];
for i=1:num_of_mice;
   
     temp_bef= mice_table{1,i}.temp_raw_f_mat;
     temp_aft= mice_table{2,i}.temp_raw_f_mat;
     total_raw_f_mat_bef=[total_raw_f_mat_bef; temp_bef];
     total_raw_f_mat_aft=[total_raw_f_mat_aft; temp_aft];
end


total_size= size(total_raw_f_mat_aft,1)*size(total_raw_f_mat_aft,2)*size(total_raw_f_mat_aft,3);
f_bef=reshape(total_raw_f_mat_bef,[1,total_size]);
f_aft=reshape(total_raw_f_mat_aft,[1,total_size]);

means=[mean(f_bef) mean(f_aft)];
sems=[std(f_bef)/sqrt(length(total_raw_f_mat_bef)) std(f_bef)/sqrt(length(total_raw_f_mat_bef))];

[h,p]=ttest2(f_bef,f_aft);

x=(1:2);
figure
errorbar (x,means,sems);


%% check spontanous activity   good for same number bef and after. otherwise need to change. with f0 as a measure
total_raw_f_mat_bef=[];
total_raw_f_mat_aft=[];
for i=1:num_of_mice
   
     temp_bef= mice_table{1,i}.temp_raw_f_mat;
     temp_aft= mice_table{2,i}.temp_raw_f_mat;
     total_raw_f_mat_bef=[total_raw_f_mat_bef; temp_bef];
     total_raw_f_mat_aft=[total_raw_f_mat_aft; temp_aft];
end


total_size= size(total_raw_f_mat_aft,1)*size(total_raw_f_mat_aft,2)*size(total_raw_f_mat_aft,3);
f_bef=reshape(total_raw_f_mat_bef,[1,total_size]);
f_aft=reshape(total_raw_f_mat_aft,[1,total_size]);

means=[mean(f_bef) mean(f_aft)];
sems=[std(f_bef)/sqrt(length(total_raw_f_mat_bef)) std(f_bef)/sqrt(length(total_raw_f_mat_bef))];

[h,p]=ttest2(f_bef,f_aft);

x=(1:2);
figure
errorbar (x,means,sems);


%% check spontanous activity   good for same number bef and after. otherwise need to change. std of off reesponse window.
total_f0_std_bef=[];
total_f0_std_aft=[];
for i=1:num_of_mice;
   
     temp_bef= mice_table{1,i}.spont_est;
     temp_aft= mice_table{2,i}.spont_est;
     total_f0_std_bef=[total_f0_std_bef; temp_bef];
     total_f0_std_aft=[total_f0_std_aft; temp_aft];
end


total_size= size(total_f0_std_aft,1)*size(total_f0_std_aft,2)*size(total_f0_std_aft,3);
f0_std_bef=reshape(total_f0_std_bef,[1,total_size]);
f0_std_aft=reshape(total_f0_std_aft,[1,total_size]);

means=[mean(f0_std_bef) mean(f0_std_aft)];
sems=[std(f0_std_bef)/sqrt(length(total_f0_std_bef)) std(f0_std_bef)/sqrt(length(total_f0_std_bef))];

[h,p]=ttest2(f0_std_bef,f0_std_aft);

x=(1:2);
figure
errorbar (x,means,sems);
%% Ariels feedback. everything here is unbiased. responses for before and
%%after are undependent

all_inhib_traces_before=[];
all_excit_traces_before=[];
all_inhib_traces_after=[];
all_excit_traces_after=[];

all_prop_ex_bef=[];
all_prop_inh_bef=[];
all_prop_ex_af=[];
all_prop_inh_af=[];

total_exc_resp_per_odor_bef=zeros(1,numOfOdors);
total_exc_resp_per_odor_after=zeros(1,numOfOdors);
total_inh_resp_per_odor_bef=zeros(1,numOfOdors);
total_inh_resp_per_odor_after=zeros(1,numOfOdors);

for i=1:num_of_mice;
   
     trace_bef_sig_exc= mice_table{1,i}.sig_excit_trials_unpaired;
     trace_bef_sig_inh= mice_table{1,i}.sig_inhib_trials_unpaired;
     
     trace_af_sig_exc= mice_table{2,i}.sig_excit_trials_unpaired;
     trace_af_sig_inh= mice_table{2,i}.sig_inhib_trials_unpaired;
     
     all_excit_traces_before=[all_excit_traces_before;trace_bef_sig_exc];
     all_inhib_traces_before=[all_inhib_traces_before;trace_bef_sig_inh];
    
     all_excit_traces_after=[all_excit_traces_after;trace_af_sig_exc];
     all_inhib_traces_after=[all_inhib_traces_after;trace_af_sig_inh];
     
    exc_resp_per_odor_bef= mice_table{1,i}.exc_prop_for_odor;
    inh_resp_per_odor_bef= mice_table{1,i}.inhib_prop_for_odor;
    
    exc_resp_per_odor_after= mice_table{2,i}.exc_prop_for_odor;
    inh_resp_per_odor_after= mice_table{2,i}.inhib_prop_for_odor;
    
    exc_resp_per_odor_bef_prop=exc_resp_per_odor_bef./size(trace_bef_sig_exc,1);
    inh_resp_per_odor_bef_prop=inh_resp_per_odor_bef./size(trace_bef_sig_inh,1);
    
    exc_resp_per_odor_after_prop=exc_resp_per_odor_after./size(trace_af_sig_exc,1);
    inh_resp_per_odor_after_prop=inh_resp_per_odor_after./size(trace_af_sig_inh,1);
    
    total_exc_resp_per_odor_bef=total_exc_resp_per_odor_bef+exc_resp_per_odor_bef;
    total_exc_resp_per_odor_after=total_exc_resp_per_odor_after+exc_resp_per_odor_after;
    total_inh_resp_per_odor_bef=total_inh_resp_per_odor_bef+inh_resp_per_odor_bef;
    total_inh_resp_per_odor_after=total_inh_resp_per_odor_after+inh_resp_per_odor_after;
    
   
    all_prop_ex_bef=[all_prop_ex_bef;exc_resp_per_odor_bef_prop];
    all_prop_inh_bef=[all_prop_inh_bef;inh_resp_per_odor_bef_prop];
    all_prop_ex_af=[all_prop_ex_af;exc_resp_per_odor_after_prop];
    all_prop_inh_af=[all_prop_inh_af;inh_resp_per_odor_after_prop];
    
  
end

 all_prop_ex_bef(:,blank)=[]
   all_prop_inh_bef(:,blank)=[]
   all_prop_ex_af(:,blank)=[]
   all_prop_inh_af(:,blank)=[]

%% number of responses for excitatory responses only:
x=1:11;


y1=total_exc_resp_per_odor_bef(k_for_plot);
y1(1)=[];
y2=total_exc_resp_per_odor_after(k_for_plot);
y2(1)=[];

figure
title('ex resp per odor')
hold on
plot(x, y1,'b')
plot(x, y2,'r')


y3=total_inh_resp_per_odor_bef(k_for_plot);
y3(1)=[];
y4=total_inh_resp_per_odor_after(k_for_plot);
y4(1)=[];

figure
title('in resp per odor')
hold on
plot(x, y3,'b')
plot(x, y4,'r')

%% 
%here there is a problem with error bars
y1_prob=y1./total_neurons_before;
y2_prob=y2./total_neurons_before;
figure
hold on
title('excitatory unpaired, number of responses per odor');

avg_across_odors_ex_bef=mean(y1_prob,2);
avg_across_odors_ex_af=mean(y2_prob,2);
sem_all_odors_bef_ex=std(y1_prob)./sqrt(length(y1_prob));
sem_all_odors_af_ex=std(y2_prob)./sqrt(length(y2));
x=(1:2);

means=[mean(mean(y1_prob)) mean(mean(y2_prob))];
general_sems=[sem_all_odors_bef_ex sem_all_odors_af_ex];   %not relevant sems
bar(x,means,0.6,'k');

hold on
errorbar(x,means,general_sems,'.k','linewidth',2);

%for inhibitory, by number. sems not relevant

y3_prob=y3./total_neurons_before;
y4_prob=y4./total_neurons_before;
figure
hold on
title('inhibitory unpaired, number of responses per odor');

avg_across_odors_ex_bef=mean(y3_prob,2);
avg_across_odors_ex_af=mean(y4_prob,2);
sem_all_odors_bef_ex=std(y3_prob)./sqrt(length(y3_prob));
sem_all_odors_af_ex=std(y4_prob)./sqrt(length(y4));
x=(1:2);

means=[mean(mean(y3_prob)) mean(mean(y4_prob))];
general_sems=[sem_all_odors_bef_ex sem_all_odors_af_ex];   %not relevant sems
bar(x,means,0.6,'k');

hold on
errorbar(x,means,general_sems,'.k','linewidth',2);




%% for porportional- all odors

y1=sort(all_prop_ex_bef,2);
y1=fliplr(y1);


y2=sort(all_prop_ex_af,2);
y2=fliplr(y2)

sem_y1=std(y1)./sqrt(size(all_prop_ex_bef,1));
sem_y2=std(y2)./sqrt(size(all_prop_ex_af,1));

x=1:11;
figure
hold on
title('excitatory unpaired, prop responses per odor');
errorbar(x, squeeze(mean(y1,1)),sem_y1,'b','linewidth',3);
errorbar(x, mean(y2,1),sem_y2,'r','linewidth',3);
ylim([0 1])
xlim([0.9 11])
%% bar 
figure
hold on;
title('excitatory unpaired, prop responses per cell');
avg_across_odors_ex_bef=mean(y1,2);
avg_across_odors_ex_af=mean(y2,2);
sem_all_odors_bef_ex=std(avg_across_odors_ex_bef)./sqrt(length(avg_across_odors_ex_bef));
sem_all_odors_af_ex=std(avg_across_odors_ex_af)./sqrt(length(avg_across_odors_ex_af));
x=(1:2);

means=[mean(mean(y1)) mean(mean(y2))];
general_sems=[sem_all_odors_bef_ex sem_all_odors_af_ex];
bar(x,means,0.6,'edgecolor','k','facecolor',[0.7 0.7 0.7],'linewidth',2);
hold on
errorbar(x,means,general_sems,'.k','linewidth',2);
%%
%proportion for inhibitory responses only for all odors:

xx=1:11;
yy1=sort(all_prop_inh_bef,2);
yy1=fliplr(yy1);
yy2=sort(all_prop_inh_af,2);
yy2=fliplr(yy2);
sem_y1=std(yy1)./sqrt(size(all_prop_inh_bef,1));
sem_y2=std(yy2)./sqrt(size(all_prop_inh_af,1));


figure
hold on
title('inhibitory unpaired, prop responses per cell');
errorbar(xx, mean(yy1,1),sem_y1,'b','linewidth',3);
errorbar(xx, mean(yy2,1),sem_y2,'r','linewidth',3);
xlim([0.9 11])

%%  bars for inhibitory responses
xx=(1:2);
means=[mean(mean(yy1)) mean(mean(yy2))];
figure
hold on
title('inhibitory unpaired, prop responses per cell');
avg_across_odors_inh_bef=mean(yy1,2);
avg_across_odors_inh_af=mean(yy2,2);
sem_all_odors_bef_inh=std(avg_across_odors_inh_bef)./sqrt(length(avg_across_odors_inh_bef));
sem_all_odors_af_inh=std(avg_across_odors_inh_af)./sqrt(length(avg_across_odors_inh_af));

bar(xx,means,0.6,'edgecolor','k','facecolor',[0.7 0.7 0.7],'linewidth',2);
hold on
errorbar(xx,means,general_sems,'.k','linewidth',2);
ylim([0 0.3])

%% new- for in and exc together:
xx=1:11;
y_in_ex_bef=y1+yy1
y_in_ex_aft=y2+yy2
sem_y_in_ex_bef=std(y_in_ex_bef)./sqrt(size(y_in_ex_bef,1));
sem_y_in_ex_aft=std(y_in_ex_aft)./sqrt(size(y_in_ex_aft,1));



figure
hold on
title('both inhibitory and excitatory unpaired, prop responses per cell average for mouse');
errorbar(xx, mean(y_in_ex_bef,1),sem_y_in_ex_bef,'b','linewidth',3);
errorbar(xx, mean(y_in_ex_aft,1),sem_y_in_ex_aft,'r','linewidth',3);
xlim([0.9 11])

%%
%%  bars for both ex and inhibitory responses
xx=(1:2);
means=[mean(mean(y_in_ex_bef)) mean(mean(y_in_ex_aft))];
figure
hold on
title('inhibitory unpaired, prop responses per cell average for mouse');
avg_across_odors_inh_ex_bef=mean(y_in_ex_bef,2);
avg_across_odors_inh_ex_af=mean(y_in_ex_bef,2);
sem_all_odors_bef_inh=std(avg_across_odors_inh_ex_bef)./sqrt(length(avg_across_odors_inh_ex_bef));
sem_all_odors_af_inh=std(avg_across_odors_inh_af)./sqrt(length(avg_across_odors_inh_af));

bar(xx,means,0.6,'edgecolor','k','facecolor',[0.7 0.7 0.7],'linewidth',2);
hold on
errorbar(xx,means,general_sems,'.k','linewidth',2);
ylim([0 0.5])
%%
% plot mean of all responses over time
bef_ex_vec=zeros(length(all_excit_traces_before)*size(all_excit_traces_before,2),size(all_excit_traces_before,3));
aft_ex_vec=bef_ex_vec;
 bef_inh_vec=bef_ex_vec;
aft_inh_vec=bef_ex_vec;
i=0;
    for k=1:size(all_excit_traces_before,1);
        for j=1:size(all_excit_traces_before,2);
            i=i+1;
    bef_ex_vec(i,:)=(all_excit_traces_before(k,j,:));
    aft_ex_vec(i,:)=(all_excit_traces_after(k,j,:));
    
    bef_inh_vec(i,:)=(all_inhib_traces_before(k,j,:));
    aft_inh_vec(i,:)=(all_inhib_traces_after(k,j,:));
        end
  
    end
a=find(sum(bef_ex_vec,2)==0);
b=find(sum(aft_ex_vec,2)==0);
c=find(sum(bef_inh_vec,2)==0);
d=find(sum(aft_inh_vec,2)==0);

bef_ex_vec(a,:)=[];
aft_ex_vec(b,:)=[];
bef_inh_vec(c,:)=[];
aft_inh_vec(d,:)=[];

bef_ex_vec_mean_sorted=sort(mean(bef_ex_vec(:,49:84),2));
aft_ex_vec_mean_sorted=sort(mean(aft_ex_vec(:,49:84),2));
bef_inh_vec_mean_sorted=sort(mean(bef_inh_vec(:,49:84),2));
aft_inh_vec_mean_sorted=sort(mean(aft_inh_vec(:,49:84),2));

x1=1:length(bef_ex_vec_mean_sorted);
x2=1:length(aft_ex_vec_mean_sorted);
x3=1:length(bef_inh_vec_mean_sorted);
x4=1:length(aft_inh_vec_mean_sorted);

figure
plot(x1,bef_ex_vec_mean_sorted,'b','linewidth',1);
hold on
plot(x2,aft_ex_vec_mean_sorted,'r','linewidth',1);
title('unpaired ex responses');
ylim([0 2])
xlim([1 max(length(bef_ex_vec_mean_sorted),length(aft_ex_vec_mean_sorted))])
%%

%plot inhibitory only
figure
plot(x3,bef_inh_vec_mean_sorted,'b','linewidth',1);
hold on
plot(x4,aft_inh_vec_mean_sorted,'r','linewidth',1);
title('unpaired inh responses');

xlim([1 max(length(bef_inh_vec_mean_sorted),length(aft_inh_vec_mean_sorted))])
ylim([-0.8 0.1])
%% all traces on top of each other only excitatory unpaired

figure
title('All unpaired excitatory traces over time')
hold on
x=1:size(bef_ex_vec,2);
for i=1:size(bef_ex_vec,1);
plot(x,bef_ex_vec(i,:),'b');
end
  

xx=1:size(aft_ex_vec,2);
for i=1:size(aft_ex_vec,1);
plot(x,aft_ex_vec(i,:),'r');
end
plot(xx,mean(aft_ex_vec),'--k','linewidth',3); 
plot(x,mean(bef_ex_vec),'k','linewidth',3)  ;

%ylim([-0.5 1])


%% average unpaired on top excitatory
%new adittion: confidence interval
sem_bef_ex=std(bef_ex_vec)./sqrt(size(bef_ex_vec,1));
low_bef_ex=mean(bef_ex_vec)-sem_bef_ex;
high_bef_ex=mean(bef_ex_vec)+sem_bef_ex;

sem_aft_ex=std(aft_ex_vec)./sqrt(size(aft_ex_vec,1));
low_aft_ex=mean(aft_ex_vec)-sem_aft_ex;
high_aft_ex=mean(aft_ex_vec)+sem_aft_ex;

ax=gca;
ax.FontSize=14;
% ax.XTick= [50 63]
% ax.XTickLabels=[]



figure
title('All unpaired excitatory traces over time')
hold on
x=1:size(bef_ex_vec,2);  
xx=1:size(aft_ex_vec,2);

plot(xx,mean(aft_ex_vec),'--k','linewidth',3); 
plot(xx,low_aft_ex,'--k','linewidth',1);
plot(xx,high_aft_ex,'--k','linewidth',1);

plot(x,mean(bef_ex_vec),'k','linewidth',3) ;
plot(x,low_bef_ex,'k','linewidth',1) ;
plot(x,high_bef_ex,'k','linewidth',1)  ;


ylim([-0.25 0.1])
xlim([0 size(aft_ex_vec,2)])

figure
title('All unpaired excitatory traces over time')
hold on
xlim([0 size(aft_ex_vec,2)])
H=shadedErrorBar(x,mean(bef_ex_vec),sem_bef_ex,'b')
H=shadedErrorBar(x,mean(aft_ex_vec),sem_aft_ex,'r')
xlabel('Time')
ylabel('dF/F')
title('All unpaired excitatory responses')

ax=gca;
ax.FontSize=14;
ax.XTick=([50 63]);
 ax.XTickLabels=[]

%% all traces on top of each other only inhibitory unpaired



figure
title('All unpaired inhibitory traces over time')
hold on
x=1:size(bef_inh_vec,2);
for i=1:size(bef_inh_vec,1);
plot(x,bef_inh_vec(i,:),'b');
end
  

xx=1:size(aft_inh_vec,2);
for i=1:size(aft_inh_vec,1);
plot(x,aft_inh_vec(i,:),'r')
end
plot(xx,mean(aft_inh_vec),'--k','linewidth',3); 
plot(x,mean(bef_inh_vec),'k','linewidth',3)  ;


ylim([-0.5 0.5])
xlim([0 size(aft_inh_vec,2)])

%% average unpaired on top inhib
%new adittion: confidence interval
sem_bef_inh=std(bef_inh_vec)./sqrt(size(bef_inh_vec,1));
low_bef_inh=mean(bef_inh_vec)-sem_bef_inh;
high_bef_inh=mean(bef_inh_vec)+sem_bef_inh;

sem_aft_inh=std(aft_inh_vec)./sqrt(size(aft_inh_vec,1));
low_aft_inh=mean(aft_inh_vec)-sem_aft_inh;
high_aft_inh=mean(aft_inh_vec)+sem_aft_inh;



figure
title('All unpaired inhibitory traces over time')
hold on
x=1:size(bef_inh_vec,2);  
xx=1:size(aft_inh_vec,2);

plot(xx,mean(aft_inh_vec),'--k','linewidth',3); 
plot(xx,low_aft_inh,'--k','linewidth',1);
plot(xx,high_aft_inh,'--k','linewidth',1);

plot(x,mean(bef_inh_vec),'k','linewidth',3) ;
plot(x,low_bef_inh,'k','linewidth',1) ;
plot(x,high_bef_inh,'k','linewidth',1)  ;


ylim([-0.25 0.1])
xlim([0 size(aft_inh_vec,2)])

figure
title('All unpaired inhibitory traces over time')
hold on
xlim([0 size(aft_inh_vec,2)])
H=shadedErrorBar(x,mean(bef_inh_vec),sem_bef_inh,'b')
H=shadedErrorBar(x,mean(aft_inh_vec),sem_aft_inh,'r')
xlabel('Time')
ylabel('dF/F')
title('All unpaired inhibitory responses')


%% from here on: all paired sig traces on top of each other

all_mice_only_sig_ex_bef_paired=[];
all_mice_only_sig_inhib_bef_paired=[];
all_mice_only_sig_ex_af_paired=[];
 all_mice_only_sig_inhib_af_paired=[];

%%all traces on top for paired responses
for i=1:num_of_mice;
    full_traces_before=mice_table{1,i}.mean_traces_mat;
    full_traces_after=mice_table{2,i}.mean_traces_mat;
%     
%      full_traces_before=mice_table{1,i}.mean_traces_mat_no_smooth;    %Just change for the above lines for plot of smoothed data
%     full_traces_after=mice_table{2,i}.mean_traces_mat_no_smooth;       
    
    meta_bef=mice_table{1,i}.meta_desicion;
    meta_af=mice_table{2,i}.meta_desicion;
    meta_bef(meta_bef==-1)=-2;
    meta_af(meta_af==-1)=-2;

    total_meta_des_bef_or_af_inh_ex=meta_bef+meta_af; %allows identification of inh vs ex
    A=repmat(total_meta_des_bef_or_af_inh_ex,1,1,size( full_traces_before,3));
    
    only_sig_ex_bef_paired=full_traces_before;
    only_sig_inhib_bef_paired=full_traces_before;
    
    only_sig_ex_af_paired=full_traces_after;
    only_sig_inhib_af_paired=full_traces_after;
    
     only_sig_ex_bef_paired(A<=0)=0;  %does not include inh-ex. does include pairs in which 1 response or more was ex + significant
     only_sig_ex_af_paired(A<=0)=0;
     
      only_sig_inhib_bef_paired(A>=-1)=0; %does not include inh-ex. does include pairs in which 1 response or more was ex + significant
      only_sig_inhib_af_paired(A>=-1)=0;
     
     all_mice_only_sig_ex_bef_paired=[all_mice_only_sig_ex_bef_paired;only_sig_ex_bef_paired];
      all_mice_only_sig_inhib_bef_paired=[all_mice_only_sig_inhib_bef_paired;only_sig_inhib_bef_paired];
       all_mice_only_sig_ex_af_paired=[all_mice_only_sig_ex_af_paired;only_sig_ex_af_paired];
        all_mice_only_sig_inhib_af_paired=[all_mice_only_sig_inhib_af_paired;only_sig_inhib_af_paired];
end

%plot

%% all traces- paired responses. inhib and excit seperatly
size_1=size(all_mice_only_sig_ex_bef_paired,1);
size_2=size(all_mice_only_sig_ex_bef_paired,2);
size_3=size(all_mice_only_sig_ex_bef_paired,3);
a=reshape(all_mice_only_sig_ex_bef_paired,[size_1*size_2 size_3]);
b=reshape( all_mice_only_sig_inhib_bef_paired,[size_1*size_2 size_3]);
c=reshape(all_mice_only_sig_ex_af_paired,[size_1*size_2 size_3]);
d=reshape(all_mice_only_sig_inhib_af_paired,[size_1*size_2 size_3]);

a(sum(a,2)==0,:)=[];
b(sum(b,2)==0,:)=[];
c(sum(c,2)==0,:)=[];
d(sum(d,2)==0,:)=[];



figure
title('All paired excitatory traces over time')
hold on
x=1:size(a,2);
for i=1:size(a,1)        
plot(x,a(i,:),'b')
end


for i=1:size(c,1)        
plot(x,c(i,:),'r')
end
plot(x,mean(a),'k', 'linewidth',2)
plot(x,mean(c),'--k', 'linewidth',2)
ylim([-0.5 3.5])

ax=gca;
ax.FontSize=14;
ax.XTick=[50 63];
ax.XTickLabels=[];

%%

%new adittion: confidence interval
sem_bef_ex=std(a)./sqrt(size(a,1));
low_bef_ex=mean(a)-sem_bef_ex;
high_bef_ex=mean(a)+sem_bef_ex;

sem_aft_ex=std(c)./sqrt(size(c,1));
low_aft_ex=mean(c)-sem_aft_ex;
high_aft_ex=mean(c)+sem_aft_ex;



figure
title('All paired excitatory traces over time')
hold on
x=1:size(a,2);  
xx=1:size(c,2);

plot(xx,mean(c),'--k','linewidth',3); 
plot(xx,low_aft_ex,'--k','linewidth',1);
plot(xx,high_aft_ex,'--k','linewidth',1);

plot(x,mean(a),'k','linewidth',3) ;
plot(x,low_bef_ex,'k','linewidth',1) ;
plot(x,high_bef_ex,'k','linewidth',1)  ;


ylim([-0.05 0.4])
xlim([0 size(aft_inh_vec,2)])
ax=gca;
ax.FontSize=14;
ax.XTick= [50 63]
ax.XTickLabels=[]

%% deltas over time

deltas_over_time_ex=(a-c)
sem_deltas=std(deltas_over_time_ex)./sqrt(size(deltas_over_time_ex,1));


low_deltas_ex=mean(deltas_over_time_ex)-sem_deltas;
high_deltas_ex=mean(deltas_over_time_ex)+sem_deltas;



figure
title('BEFORE-AFTER OVER TIME EX')
hold on
x=1:size(a,2);  
xx=1:size(c,2);

ax=gca;
ax.FontSize=14;
ax.XTick= [0 50 63 84 112]
ax.XTickLabels=[-7 0 2 9]

plot(x,mean(deltas_over_time_ex),'k','linewidth',3) ;
plot(x,low_deltas_ex,'k','linewidth',1) ;
plot(x,high_deltas_ex,'k','linewidth',1)  ;


ylim([-0.05 0.15])
xlim([0 size(aft_inh_vec,2)])

deltas_over_time_in=(b-d)
sem_deltas=std(deltas_over_time_in)./sqrt(size(deltas_over_time_in,1));


low_deltas_in=mean(deltas_over_time_in)-sem_deltas;
high_deltas_in=mean(deltas_over_time_in)+sem_deltas;



figure
title('BEFORE-AFTER OVER TIME IN')
hold on
x=1:size(b,2);  
xx=1:size(d,2);

ax=gca;
ax.FontSize=14;
ax.XTick= [50 63]
ax.XTickLabels=[]

plot(x,mean(deltas_over_time_in),'k','linewidth',3) ;
plot(x,low_deltas_in,'k','linewidth',1) ;
plot(x,high_deltas_in,'k','linewidth',1)  ;


ylim([-0.05 0.15])
xlim([0 size(aft_inh_vec,2)])
%%


figure
title('All paired excitatory traces over time')
ylim([-0.4 0.4])
xlim([0 size(aft_inh_vec,2)])
hold on
H=shadedErrorBar(x,mean(a),sem_bef_ex,'b')
H=shadedErrorBar(x,mean(c),sem_aft_ex,'r')

ax=gca;
ax.FontSize=14;
ax.XTick= [50 63]
ax.XTickLabels=[]

% now normalized

a_c_norm_factor=1/max(mean(a));
a_norm=a.*a_c_norm_factor;
c_norm=c.*a_c_norm_factor;


sem_bef_ex_norm=std(a_norm)./sqrt(size(a_norm,1));
sem_aft_ex_norm=std(c_norm)./sqrt(size(c_norm,1));

figure
title('All paired normalized! excitatory traces over time')

xlim([0 size(aft_inh_vec,2)])
hold on


H=shadedErrorBar(x,mean(a_norm),sem_bef_ex_norm,'b')
H=shadedErrorBar(x,mean(c_norm),sem_aft_ex_norm,'r')

ax=gca;
ax.FontSize=14;
ax.XTick= [50 63]
ax.XTickLabels=[]

%%
figure
title('All paired inhibitory traces over time')
hold on
x=1:size(a,2);
for i=1:size(b,1)        
plot(x,b(i,:),'b')
end


for i=1:size(d,1)        
plot(x,d(i,:),'r')
end
plot(x,mean(b),'k', 'linewidth',2)
plot(x,mean(d),'--k', 'linewidth',2)


sem_bef_in=std(b)./sqrt(size(b,1));
low_bef_in=mean(b)-sem_bef_in;
high_bef_in=mean(b)+sem_bef_in;

sem_aft_in=std(d)./sqrt(size(d,1));
low_aft_in=mean(d)-sem_aft_in;
high_aft_in=mean(d)+sem_aft_in;


hold on
figure
title('All paired inhibitory traces over time')
ylim([-0.2 0.05])
xlim([0 size(aft_inh_vec,2)])
hold on
H=shadedErrorBar(x,mean(b),sem_bef_in,'b')
H=shadedErrorBar(x,mean(d),sem_aft_in,'r')

ax=gca;
ax.FontSize=14;
ax.XTick= [50 63]
ax.XTickLabels=[]

%% together paired
hold on
figure
title('All paired inhibitory and ecitatory traces over time')
ylim([-0.3 0.3])
xlim([0 size(aft_inh_vec,2)])
hold on
H=shadedErrorBar(x,mean(b),sem_bef_in,'b')
H=shadedErrorBar(x,mean(d),sem_aft_in,'r')

H=shadedErrorBar(x,mean(a),sem_bef_ex,'b')
H=shadedErrorBar(x,mean(c),sem_aft_ex,'r')

ax=gca;
ax.FontSize=14;
ax.XTick= [50 63]
ax.XTickLabels=[]


%% together normalized paired
% b_resp=b(:,49:end)
% d_resp=d(:,49:end)
% a_resp=a(:,49:end)
% c_resp=c(:,49:end)
% 
% b_min=min(b_resp')
% d_min=min(d_resp')
% 
% a_max=max(a_resp')
% c_max=max(c_resp')
% 
% 
% a_norm=a./a_max'
% c_norm=c./c_max'
% b_norm=b./abs(b_min')
% d_norm=d./abs(d_min')
% 
% 
% sem_bef_in_norm=std(b_norm)./sqrt(size(b_norm,1));
% sem_aft_in_norm=std(d_norm)./sqrt(size(d_norm,1));
% sem_bef_ex_norm=std(a_norm)./sqrt(size(a_norm,1));
% sem_aft_ex_norm=std(c_norm)./sqrt(size(c_norm,1));
% 
% 
% 
% 

a_c_norm_factor=1/max(mean(a));
a_norm=a.*a_c_norm_factor;
c_norm=c.*a_c_norm_factor;

b_d_norm_factor=abs(1/min(mean(b)));
b_norm=b.*b_d_norm_factor;
d_norm=d.*b_d_norm_factor;

sem_bef_in_norm=std(b_norm)./sqrt(size(b_norm,1));
sem_aft_in_norm=std(d_norm)./sqrt(size(d_norm,1));
sem_bef_ex_norm=std(a_norm)./sqrt(size(a_norm,1));
sem_aft_ex_norm=std(c_norm)./sqrt(size(c_norm,1));

figure
title('All paired normalized! inhibitory and ecitatory traces over time')

xlim([0 size(aft_inh_vec,2)])
hold on
H=shadedErrorBar(x,mean(b_norm),sem_bef_in_norm,'b')
H=shadedErrorBar(x,mean(d_norm),sem_aft_in_norm,'r')

H=shadedErrorBar(x,mean(a_norm),sem_bef_ex_norm,'b')
H=shadedErrorBar(x,mean(c_norm),sem_aft_ex_norm,'r')

ax=gca;
ax.FontSize=14;
ax.XTick= [50 63]
ax.XTickLabels=[]

%%'P over time

p_ex=zeros(1,size(a_norm,2));
p_in=p_ex;

d_prime_ex=zeros(1,size(a_norm,2));
d_prime_in=d_prime_ex;

t_val_ex=zeros(1,size(a_norm,2));
t_val_in=t_val_ex;

delta_size_ex=zeros(1,size(a_norm,2));
delta_size_in=delta_size_ex;

d_prime = @(x,y) abs(mean(x)-mean(y))/sqrt(0.5*((std(x))^2+(std(y))^2));


for i=1:length(p_ex);
    [h_ex, p_exit]=(ttest(a_norm(:,i),c_norm(:,i)));
    [h_in, p_inhib]=(ttest(b_norm(:,i),d_norm(:,i)));
    
    p_ex(i)=p_exit;
    p_in(i)= p_inhib;
    
    d_prime_ex(i)=d_prime(a_norm(:,i),c_norm(:,i));
    d_prime_in(i)=d_prime(b_norm(:,i),d_norm(:,i));
    
    [h,p,ci,stats] = ttest(a_norm(:,i),c_norm(:,i));
    t_val=stats.tstat;
    
     t_val_ex(i)=abs(t_val);
     
     [h,p,ci,stats] = ttest(b_norm(:,i),d_norm(:,i));
     t_val=stats.tstat;
    
     t_val_in(i)=abs(t_val);
     
     delta_size_ex(i)=abs(mean(c_norm(:,i)-a_norm(:,i)));
     delta_size_in(i)=abs(mean(d_norm(:,i)-b_norm(:,i)));
end

%% figures for difference in time, colorvec

% figure
% heatmap(p_ex)
% colorbar()
% % cmin=min([p_ex';p_in']);
% % cmax=max([p_ex';p_in']);
% cmin=min([p_ex]);
% cmax=max([p_ex]);
% % cmin=0;
% % cmax=0.001;
% caxis([cmin, cmax])
% % colormap(jet)
% 
% 
% 
% title('p values for ex traces')
% 
% figure
% heatmap(p_in)
% colorbar()
% % cmin=min([p_ex';p_in']);
% % cmax=max([p_ex';p_in']);
% 
% cmin=min([p_in]);
% cmax=max([p_in]);
% caxis([cmin, cmax])
% title('p values for in traces')
% 
% figure
% imagesc(p_ex)
% colorbar()
% cmin=min([p_ex';p_in']);
% cmax=max([p_ex';p_in']);
% caxis([cmin, cmax])
% title('p values for ex traces')
% ax=gca;
% ax.FontSize=14;
% ax.XTick= [0 49 63 112]
% ax.XTickLabels=[-7 0 2 9]
% 
% figure
% imagesc(p_in)
% colorbar()
% cmin=min([p_ex';p_in']);
% cmax=max([p_ex';p_in']);
% caxis([cmin, cmax])
% title('p values for in traces')
% ax=gca;
% ax.FontSize=14;
% ax.XTick= [0 49 63 112]
% ax.XTickLabels=[-7 0 2 9]
% 
% 
% 
% figure
% heatmap(d_prime_ex)
% % cmin=min([d_prime_ex';d_prime_in']);
% % cmax=max([d_prime_ex';d_prime_in']);
% cmin=min([d_prime_ex]);
% cmax=max([d_prime_ex]);
% 
% caxis([cmin, cmax])
% colorbar()
% 
% title('d primes for ex traces')
% 
% figure
% heatmap(d_prime_in)
% % cmin=min([d_prime_ex';d_prime_in']);
% % cmax=max([d_prime_ex';d_prime_in']);
% 
% cmin=min([d_prime_in]);
% cmax=max([d_prime_in]);
% 
% caxis([cmin, cmax])
% colorbar()
% 
% title('d primes for in traces')
% 
% figure
% imagesc(d_prime_ex)
% colorbar()
% cmin=min([d_prime_ex';d_prime_in']);
% cmax=max([d_prime_ex';d_prime_in']);
% caxis([cmin, cmax])
% title('d primes for ex traces')
% ax=gca;
% ax.FontSize=14;
% ax.XTick= [0 49 63 112]
% ax.XTickLabels=[-7 0 2 9]
% 
% 
% figure
% imagesc(d_prime_in)
% colorbar()
% cmin=min([d_prime_ex';d_prime_in']);
% cmax=max([d_prime_ex';d_prime_in']);
% caxis([cmin, cmax])
% title('d primes for in traces')
% ax=gca;
% ax.FontSize=14;
% ax.XTick= [0 49 63 112]
% ax.XTickLabels=[-7 0 2 9]
% 
% 
% figure
% heatmap(t_val_ex)
% colorbar()
% cmin=min([t_val_ex';t_val_in']);
% cmax=max([t_val_ex';t_val_in']);
% 
% cmin=min([t_val_ex]);
% cmax=max([t_val_ex]);
% 
% caxis([cmin, cmax])
% % caxis([0 12])
% colormap(jet)
% title('t values for ex traces')
% 
% figure
% heatmap(t_val_in)
% colorbar()
% % cmin=min([t_val_ex';t_val_in']);
% % cmax=max([t_val_ex';t_val_in']);
% 
% cmin=min([t_val_in;]);
% cmax=max([t_val_in;]);
% 
% caxis([cmin, cmax])
% caxis([0 12])
% title('t values for in traces')
% colormap(jet)
% 
%%
figure
hold on



ylim([-20 20])
width=10;
height=1;

imagesc(t_val_ex)
hold on



ylim([-20 20])
% cmin=min([t_val_ex';t_val_in']);
% cmax=max([t_val_ex';t_val_in']);

cmax=max([t_val_ex]);
cmax_ex=cmax;

if condition(1)=='C'
General_results_experiment=load('General_results_experiment');
cmax_ex=General_results_experiment.cmax_ex_t_val;
end

if condition(1)=='m'
General_results_awake=load('General_results_awake');
cmax_ex=General_results_awake.cmax_ex_t_val;
end

if condition(1)=='W'
General_results_awake=load('General_results_awake');
cmax_ex=General_results_awake.cmax_ex_t_val;
end


if condition(1)=='s'
General_results_awake=load('General_results_first_timelapse');
cmax_in=General_results_awake.cmax_ex_t_val;
end

if condition(1)=='t'
General_results_awake=load('General_results_first_timelapse');
cmax_in=General_results_awake.cmax_ex_t_val;
end

cmax_ex_t_val=cmax_ex;
cmax=cmax_ex;


% caxis([cmin, cmax])
caxis([0  cmax])
colorbar()
ax=gca;
ax.FontSize=14;
ax.XTick= [0 50 63 84 112]
ax.XTickLabels=[-7 0 2 9]
pos = get(gcf, 'Position');
%set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w');

title(['t values for ex traces ' 'cmax= ' num2str(cmax) ' ' condition])

%%
figure
imagesc(t_val_in)
hold on

ylim([-20 20])

cmax_in=max([t_val_in]);

if condition(1)=='C'
General_results_experiment=load('General_results_experiment');
cmax_in=General_results_experiment.cmax_in_t_val;
end

if condition(1)=='m'
General_results_awake=load('General_results_awake');
cmax_in=General_results_awake.cmax_in_t_val;
end

if condition(1)=='W'
General_results_awake=load('General_results_awake');
cmax_in=General_results_awake.cmax_in_t_val;
end

if condition(1)=='s'
General_results_awake=load('General_results_first_timelapse');
cmax_in=General_results_awake.cmax_in_t_val;
end

if condition(1)=='t'
General_results_awake=load('General_results_first_timelapse');
cmax_in=General_results_awake.cmax_in_t_val;
end



cmax_in_t_val=cmax_in;
% 
% cmin=min([t_val_ex';t_val_in']);
% cmax=max([t_val_ex';t_val_in']);

cmax=cmax_in;

% caxis([cmin, cmax])
caxis([0 cmax])
colorbar()
ax=gca;
ax.FontSize=14;
ax.XTick= [0 50 63 84 112]
ax.XTickLabels=[-7 0 2 9]
pos = get(gcf, 'Position');
%set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w');
title(['t values for in traces' 'cmax= ' num2str(cmax) ' ' condition])

%%
% 
% 
% figure
% heatmap(delta_size_ex)
% colorbar()
% cmin=min([delta_size_ex';delta_size_in']);
% cmax=max([delta_size_ex';delta_size_in']);
% %caxis([cmin, cmax])
% %caxis([0 12])
% colormap(jet)
% title('delta_size for ex traces')
% 
% figure
% heatmap(delta_size_in)
% colorbar()
% cmin=min([delta_size_ex';delta_size_in']);
% cmax=max([delta_size_ex';delta_size_in']);
% caxis([cmin, cmax])
% % caxis([0 12])
% title('delta_size for in traces')
% colormap(jet)


% width=10;
% height=1;
% 
% cmax_ex=max([delta_size_ex]);
% 
% if condition(1)=='C'
% General_results_experiment=load('General_results_experiment');
% cmax_ex=General_results_experiment.cmax_ex;
% end
% 
% if condition(1)=='m'
% General_results_awake=load('General_results_awake');
% cmax_ex=General_results_awake.cmax_ex;
% end
% 
% if condition(1)=='W'
% General_results_awake=load('General_results_awake');
% cmax_ex=General_results_awake.cmax_ex;
% end
% 
% 
% cmax_ex_deltas=cmax_ex;
%%
width=10;
height=1;

figure
imagesc(delta_size_ex)
hold on



ylim([-20 20])
% cmin=min([delta_size_ex';delta_size_in']);
% cmax=max([delta_size_ex';delta_size_in']);

% cmin_ex=min([delta_size_ex]);
cmax_ex=max([delta_size_ex]);



if condition(1)=='C'
General_results_experiment=load('General_results_experiment');
cmax_ex=General_results_experiment.cmax_ex_deltas;
end

if condition(1)=='m'
General_results_awake=load('General_results_awake');
cmax_ex=General_results_awake.cmax_ex_deltas;
end

if condition(1)=='W'
General_results_awake=load('General_results_awake');
cmax_ex=General_results_awake.cmax_ex_deltas;
end


if condition(1)=='s'
General_results_awake=load('General_results_first_timelapse');
cmax_ex=General_results_awake.cmax_ex_deltas;
end

if condition(1)=='t'
General_results_awake=load('General_results_first_timelapse');
cmax_ex=General_results_awake.cmax_ex_deltas;
end

cmax_ex_deltas=cmax_ex;

cmax=cmax_ex;

% caxis([cmin_ex, cmax_ex])
caxis([0, cmax])
% caxis([0 12])
colorbar()
ax=gca;
ax.FontSize=14;
ax.XTick= [0 49 62 112]
ax.XTickLabels=[-7 0 2 9]
pos = get(gcf, 'Position');
%set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w');



title(['delta_size for ex traces' 'cmax= ' num2str(cmax) ' ' condition])
%%

figure
imagesc(delta_size_in)
hold on



ylim([-20 20])
% cmin_in=min([delta_size_in]);
cmax_in=max([delta_size_in]);

if condition(1)=='C'
General_results_experiment=load('General_results_experiment');
cmax_in=General_results_experiment.cmax_in_deltas;
end

if condition(1)=='m'
General_results_awake=load('General_results_awake');
cmax_in=General_results_awake.cmax_in_deltas;
end

if condition(1)=='W'
General_results_awake=load('General_results_awake');
cmax_in=General_results_awake.cmax_in_deltas;
end


if condition(1)=='s'
General_results_awake=load('General_results_first_timelapse');
cmax_in=General_results_awake.cmax_in_deltas;
end

if condition(1)=='t'
General_results_awake=load('General_results_first_timelapse');
cmax_in=General_results_awake.cmax_in_deltas;
end

cmax_in_deltas=cmax_in;
cmax=cmax_in;

% cmin=min([delta_size_ex';delta_size_in']);
% cmax=max([delta_size_ex';delta_size_in']);
% caxis([cmin, cmax])
caxis([0, cmax])
% caxis([0 12])
colorbar()
ax=gca;
ax.FontSize=14;
ax.XTick= [0 49 62 112]
ax.XTickLabels=[-7 0 2 9]
pos = get(gcf, 'Position');
%set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w');


title(['delta_size for in traces' 'cmax= ' num2str(cmax) ' ' condition])


%% %% together unpaired

figure
title('All unpaired excitatory traces over time')
hold on
xlim([0 size(aft_ex_vec,2)])
H=shadedErrorBar(x,mean(bef_ex_vec),sem_bef_ex,'b')
H=shadedErrorBar(x,mean(aft_ex_vec),sem_aft_ex,'r')
xlabel('Time')
ylabel('dF/F')
title('All unpaired excitatory responses')









hold on
xlim([0 size(aft_inh_vec,2)])
H=shadedErrorBar(x,mean(bef_inh_vec),sem_bef_inh,'b')
H=shadedErrorBar(x,mean(aft_inh_vec),sem_aft_inh,'r')
xlabel('Time')
ylabel('dF/F')
title('All unpaired inhibitory responses')



title('All unpaired inhibitory and ecitatory traces over time')
ylim([-0.3 0.4])
xlim([0 size(aft_inh_vec,2)])

%%
%% average unpaired on top excitatory
%new adittion: confidence interval
sem_bef_ex=std(bef_ex_vec)./sqrt(size(bef_ex_vec,1));
low_bef_ex=mean(bef_ex_vec)-sem_bef_ex;
high_bef_ex=mean(bef_ex_vec)+sem_bef_ex;

sem_aft_ex=std(aft_ex_vec)./sqrt(size(aft_ex_vec,1));
low_aft_ex=mean(aft_ex_vec)-sem_aft_ex;
high_aft_ex=mean(aft_ex_vec)+sem_aft_ex;



figure
title('All unpaired excitatory traces over time')
hold on
x=1:size(bef_ex_vec,2);  
xx=1:size(aft_ex_vec,2);

plot(xx,mean(aft_ex_vec),'--k','linewidth',3); 
plot(xx,low_aft_ex,'--k','linewidth',1);
plot(xx,high_aft_ex,'--k','linewidth',1);

plot(x,mean(bef_ex_vec),'k','linewidth',3) ;
plot(x,low_bef_ex,'k','linewidth',1) ;
plot(x,high_bef_ex,'k','linewidth',1)  ;


ylim([-0.25 0.1])
xlim([0 size(aft_ex_vec,2)])

figure
title('All unpaired excitatory traces over time')
hold on
xlim([0 size(aft_ex_vec,2)])
H=shadedErrorBar(x,mean(bef_ex_vec),sem_bef_ex,'b')
H=shadedErrorBar(x,mean(aft_ex_vec),sem_aft_ex,'r')
xlabel('Time')
ylabel('dF/F')
title('All unpaired excitatory responses')

%% all traces on top of each other only inhibitory unpaired



figure
title('All unpaired inhibitory traces over time')
hold on
x=1:size(bef_inh_vec,2);
for i=1:size(bef_inh_vec,1);
plot(x,bef_inh_vec(i,:),'b');
end
  

xx=1:size(aft_inh_vec,2);
for i=1:size(aft_inh_vec,1);
plot(x,aft_inh_vec(i,:),'r')
end
plot(xx,mean(aft_inh_vec),'--k','linewidth',3); 
plot(x,mean(bef_inh_vec),'k','linewidth',3)  ;


ylim([-0.5 0.5])
xlim([0 size(aft_inh_vec,2)])

%% average unpaired on top inhib
%new adittion: confidence interval
sem_bef_inh=std(bef_inh_vec)./sqrt(size(bef_inh_vec,1));
low_bef_inh=mean(bef_inh_vec)-sem_bef_inh;
high_bef_inh=mean(bef_inh_vec)+sem_bef_inh;

sem_aft_inh=std(aft_inh_vec)./sqrt(size(aft_inh_vec,1));
low_aft_inh=mean(aft_inh_vec)-sem_aft_inh;
high_aft_inh=mean(aft_inh_vec)+sem_aft_inh;



figure
title('All unpaired inhibitory traces over time')
hold on
x=1:size(bef_inh_vec,2);  
xx=1:size(aft_inh_vec,2);

plot(xx,mean(aft_inh_vec),'--k','linewidth',3); 
plot(xx,low_aft_inh,'--k','linewidth',1);
plot(xx,high_aft_inh,'--k','linewidth',1);

plot(x,mean(bef_inh_vec),'k','linewidth',3) ;
plot(x,low_bef_inh,'k','linewidth',1) ;
plot(x,high_bef_inh,'k','linewidth',1)  ;


ylim([-0.25 0.1])
xlim([0 size(aft_inh_vec,2)])

figure
title('All unpaired inhibitory traces over time')
hold on
xlim([0 size(aft_inh_vec,2)])
H=shadedErrorBar(x,mean(bef_inh_vec),sem_bef_inh,'b')
H=shadedErrorBar(x,mean(aft_inh_vec),sem_aft_inh,'r')
xlabel('Time')
ylabel('dF/F')
title('All unpaired inhibitory responses')




%%



%% save for later comparisons between groups
disp('show me');
save(['Mice_table_' num2str(condition)], 'mice_table','-v7.3')

if condition(1)=='C'
save(['General_results_control'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','stats','-v7.3')%,'string')
end

if condition(1)=='E'
save(['General_results_experiment'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','cmax_in_t_val','cmax_ex_t_val','cmax_in_deltas','cmax_ex_deltas','stats','-v7.3')%,'string')
end

if condition(1)=='a'
save(['General_results_awake'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','cmax_in_t_val','cmax_ex_t_val','cmax_in_deltas','cmax_ex_deltas','stats','-v7.3')%,'string')
end

if condition(1)=='m'
save(['General_results_awake_cont'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','stats','-v7.3')%,'string')
end

if condition(1)=='f'
save(['General_results_first_timelapse'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','cmax_in_t_val','cmax_ex_t_val','cmax_in_deltas','cmax_ex_deltas','stats','-v7.3')%,'string')
end

if condition(1)=='s'
save(['General_results_sec_timelapse'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','stats','-v7.3')%,'string')
end

if condition(1)=='t'
save(['General_results_third_timelapse'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','stats','-v7.3')%,'string')
end

if condition(1)=='W'
save(['General_results_Wt_cno_cont'],'abs_deltas_for_boxplots', 'stats_for_chi', 'mice_deltas','prop_for_comparison','deltas_for_boxplots','all_total_meta_des_bef_af_inh_ex','total_prop_of_cells_before','total_prop_of_cells_after','deltas_over_time_ex','deltas_over_time_in','a','b','c','d','stats','-v7.3')%,'string')
end

%     
%     