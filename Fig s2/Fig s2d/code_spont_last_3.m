clc
clear all
close all

load('Spont_analysis_saline_cont');
saline_bef=all_bsln_bef;
saline_aft=all_bsln_aft;
deltas_bsln_saline=saline_aft-saline_bef;


load('Spont_analysis_cno_cont');
cno_cont_bef=all_bsln_bef;
cno_cont_aft=all_bsln_aft;
deltas_bsln_cno_cont=cno_cont_aft-cno_cont_bef;

deltas_bsln_general_cont=[deltas_bsln_saline;deltas_bsln_cno_cont];

load('Spont_analysis_exp_mice');
exp_mice_bef=all_bsln_bef;
exp_mice_aft=all_bsln_aft;
deltas_bsln_exp=exp_mice_aft-exp_mice_bef;

x_exp = linspace(0.7,1.3,size(deltas_bsln_exp,1))';
x_saline=linspace(1.7,2.3,size(deltas_bsln_saline,1))';
x_cno_cont=linspace(2.7,3.3,size(deltas_bsln_cno_cont,1))';

sem_exp=std(deltas_bsln_exp)/sqrt(length(deltas_bsln_exp));
sem_saline=std(deltas_bsln_saline)/sqrt(length(deltas_bsln_saline));
sem_cno_cont=std(deltas_bsln_cno_cont)/sqrt(length(deltas_bsln_cno_cont));

sems=[sem_exp sem_saline sem_cno_cont];
means=[mean(deltas_bsln_exp) mean(deltas_bsln_saline) mean(deltas_bsln_cno_cont)];
figure
hold on
plot(x_exp,deltas_bsln_exp,'ok','MarkerSize',8)
plot(x_saline,deltas_bsln_saline,'ok','MarkerSize',8)
plot(x_cno_cont,deltas_bsln_cno_cont,'ok','MarkerSize',8)

errorbar(1:3,means,sems,'.r','linewidth',2)
plot(0.5:3.5,[0 0 0 0],'--')
% ylim([-3 3])

x=[deltas_bsln_exp;deltas_bsln_saline;deltas_bsln_cno_cont];

% deltas_for_boxplots=[all_significant_responses_integrals_before-all_significant_responses_integrals_after]; %by integral
% 
% abs_deltas_for_boxplots=[abs(all_significant_responses_integrals_before)-abs(all_significant_responses_integrals_after)]; %added at 11.1.19


g = [ones(size(deltas_bsln_exp)); 2*ones(size(deltas_bsln_saline));3*ones(size(deltas_bsln_cno_cont))];
% boxplot(x,g)
% figure
hold on
boxplot(x,g,'symbol','');
%ylim([-2,2])
x=0:4;
y=[0 0 0 0 0]
plot(x,y,'--')


% 
% boxplot([deltas_exp,deltas_saline,deltas_cno_cont])
% ,
%% merged controls- only box plot
x=[deltas_bsln_exp;deltas_bsln_general_cont];

% deltas_for_boxplots=[all_significant_responses_integrals_before-all_significant_responses_integrals_after]; %by integral
% 
% abs_deltas_for_boxplots=[abs(all_significant_responses_integrals_before)-abs(all_significant_responses_integrals_after)]; %added at 11.1.19


g = [ones(size(deltas_bsln_exp)); 2*ones(size(deltas_bsln_general_cont))];
% boxplot(x,g)
figure
hold on
boxplot(x,g,'symbol','');
ylim([-2.5,2.5])
x=0:4;
y=[0 0 0 0 0]
plot(x,y,'--')


% add the dots:
x_exp = linspace(0.9,1.1,size(deltas_bsln_exp,1))';
x_general_cont=linspace(1.9,2.1,size(deltas_bsln_general_cont,1))';

plot(x_exp,deltas_bsln_exp,'ok','MarkerSize',1)
plot(x_general_cont,deltas_bsln_general_cont,'ok','MarkerSize',1)





[p,h]=signrank(deltas_bsln_exp)
[p,h]=signrank(deltas_bsln_general_cont)

ranksum(deltas_bsln_exp,deltas_bsln_general_cont)



%% compare size effect to effect in spontanous activity
cont_saline_resp=load('responses_per_cell_cont_saline')
cont_cno_resp=load('responses_per_cell_cont_cno')
exp_resp=load('responses_per_cell_exp')

resp_exp_bef=exp_resp.sort_exp_bef;
resp_exp_aft=exp_resp.sort_exp_aft;
exp_curve_sd_before=std(resp_exp_bef,0,2);
exp_curve_sd_after=std(resp_exp_aft,0,2);


merged_cont_resp_before=[cont_saline_resp.sort_cont_bef;cont_cno_resp.sort_cont_bef]
merged_cont_resp_after=[cont_saline_resp.sort_cont_aft;cont_cno_resp.sort_cont_aft]
cont_curve_sd_before=std(merged_cont_resp_before,0,2);
cont_curve_sd_after=std(merged_cont_resp_after,0,2);

merged_cont_bsln=[deltas_bsln_saline;deltas_bsln_cno_cont];
deltas_sds_cont=cont_curve_sd_before-cont_curve_sd_after;
deltas_sds_exp=exp_curve_sd_before-exp_curve_sd_after;
ratio_sds_exp=exp_curve_sd_before./exp_curve_sd_after

A=[deltas_bsln_exp';deltas_sds_exp']'
B=[merged_cont_bsln';deltas_sds_cont']'
[R,P,RL,RU] = corrcoef(A)
[R,P,RL,RU] = corrcoef(B)

[RHO,PVAL] = corr(deltas_bsln_exp,deltas_sds_exp,'Type','Spearman')
figure
plot(deltas_bsln_exp,deltas_sds_exp,'*')

figure
plot(deltas_bsln_exp,ratio_sds_exp,'*')
[RHO,PVAL] = corr(deltas_bsln_exp,ratio_sds_exp,'Type','Spearman')

max_exp_bef=max(resp_exp_bef');
max_exp_aft=max(resp_exp_aft');
deltas_max_resp_exp=max_exp_bef-max_exp_aft;

max_cont_bef=max(merged_cont_resp_before');
max_cont_aft=max(merged_cont_resp_after');

figure
plot(deltas_bsln_exp,deltas_max_resp_exp,'*')