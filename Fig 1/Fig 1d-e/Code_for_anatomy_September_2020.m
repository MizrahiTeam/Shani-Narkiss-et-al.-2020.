clc
close all
clear all

%load data from diffrent sheets into a struct
filename = 'F:\ALL DATA FOR PAPER HSN\CELL_COUNT_WITH_NO_TAM.xlsx';  %Use the Excel file provided in the folder, and change path accordingly


mice_table=struct;

no_tam={'S1','S2','S3'};
two_weeks = {'M8','S9' ,'S14', 'S15', 'S16'};
one_month={'S4','S5','S10','S11','S12','S13'};
two_months={'M1','M13','M16','M17'};
allmice={'S1','S2','S3','M8','S9' ,'S14', 'S15', 'S16','S4','S5','S10','S11','S12','S13','M1','M13','M16','M17'};
n_slice_per_mouse=3;

Mice_gcl_mat=zeros(length(allmice),n_slice_per_mouse);


for i=1:length(allmice);
datatable = readtable(filename,'Sheet',string(allmice(i)));
S= table2struct(datatable);
mice_table(i).name=(string(allmice(i)));
mice_table(i).table= S;
end

for i=1:length(allmice);
    for j=1:n_slice_per_mouse;
ClassName = class(mice_table(i).table((j*2)-1).CellsPerMm_3);
if ClassName(1)=='c'        
Mice_gcl_mat(i,j)=str2num(mice_table(i).table((j*2)-1).CellsPerMm_3);   %this is because data is located at 1 3 5 
end

if ClassName(1)=='d'        
Mice_gcl_mat(i,j)=(mice_table(i).table((j*2)-1).CellsPerMm_3);   %this is because data is located at 1 3 5 
end
    end
end

a=length(no_tam);
b=length(two_weeks);
c=length(one_month);
d=length(two_months);

no_tam_mat=Mice_gcl_mat(1:a,:);
two_weeks_mat=Mice_gcl_mat((a+1):(a+b),:);
one_month_mat=Mice_gcl_mat((a+b+1):(a+b+c),:);
two_months_mat=Mice_gcl_mat((a+b+c+1):(a+b+c+d),:);

%% plot gcl for sd by mouse
mean_gcl_per_mouse_no_tam=mean(no_tam_mat,2)
mean_gcl_per_mouse_two_weeks=mean(two_weeks_mat,2)
mean_gcl_per_mouse_one_month=mean(one_month_mat,2)
mean_gcl_per_mouse_two_months=mean(two_months_mat,2)

sem_gcl_per_mouse_no_tam=std(mean_gcl_per_mouse_no_tam)/sqrt(size(mean_gcl_per_mouse_no_tam,1))
sem_gcl_per_mouse_two_weeks=std(mean_gcl_per_mouse_two_weeks)/sqrt(size(mean_gcl_per_mouse_two_weeks,1))
sem_gcl_per_mouse_one_month=std(mean_gcl_per_mouse_one_month)/sqrt(size(mean_gcl_per_mouse_one_month,1))
sem_gcl_per_mouse_two_months=std(mean_gcl_per_mouse_two_months)/sqrt(size(mean_gcl_per_mouse_two_months,1))




means_gcl_per_mouse=[mean(mean_gcl_per_mouse_no_tam) mean(mean_gcl_per_mouse_two_weeks) mean(mean_gcl_per_mouse_one_month) mean(mean_gcl_per_mouse_two_months)];
sems_gcl_per_mouse=[sem_gcl_per_mouse_no_tam sem_gcl_per_mouse_two_weeks sem_gcl_per_mouse_one_month sem_gcl_per_mouse_two_months];

%FROM LAGACE: Because the adult mouse OB GCL has ?410,000 cells/mm3 (Parrish-Aungst et al., 2007), YFP? cells represented 0.3, 1.1, and 2.6% of the total OB GCL density at 30, 65, and 100 d after TAM, respectively


lagace0=[0];
lagace12=[410000*(0.15/100)]; %this one is not reported- estimated by eye
lagace30=[410000*(0.3/100)];   
lagace65=[410000*(1.1/100)];
lagace100=[410000*(2.6/100)]; %not in uae

x_Lagace=[0 12 30 65]
y_Lagace=[lagace0 lagace12 lagace30 lagace65]; %See estimation from Lagace et al., 2007.

x=[0 14 30 60];


all_mice_vals_times(:,1)=[0 0 0 14 14 14 14 14 30 30 30 30 30 30 60 60 60 60];

all_mice_vals_times(:,2)=[mean_gcl_per_mouse_no_tam;mean_gcl_per_mouse_two_weeks;mean_gcl_per_mouse_one_month;mean_gcl_per_mouse_two_months];

xx=all_mice_vals_times(:,1);
y=all_mice_vals_times(:,2);

b1 = xx\y;

yCalc1 = b1*xx;
scatter(xx,y)
hold on


figure
hold on
plot(xx,yCalc1,'k','linewidth',2)
errorbar(x,means_gcl_per_mouse,sems_gcl_per_mouse,'ko','linewidth',2)
title('Granular Cells Layer count. N=17 mice, 51 slices','FontSize',25)
xlabel('Weeks post Tamoxifen')
ylabel('Number of cells per mm^3')
a=gca;
axis square;
a.FontSize=20;
plot(x_Lagace,y_Lagace,'k--','linewidth',2)
xlim([-5 65])
ax=gca;
ax.FontSize=14;
ax.XTick= [0 14 30 60]
ax.XTickLabels=[0 2 4 8]

%plot mice as data points:
x1=[-1:1];
x1=[-0.5 0 0.5];
plot(x1,mean_gcl_per_mouse_no_tam,'ro')
x2=[13 13.5 14 14.5 15];
%x2=[12 13 14 15 16];
plot(x2,mean_gcl_per_mouse_two_weeks,'ro');
x3=[28.5 29 29.5 30 30.5 31];
plot(x3,mean_gcl_per_mouse_one_month,'ro');
x4=[59 59.5 60 60.5];
plot(x4,mean_gcl_per_mouse_two_months,'ro');
%% pearson correlation and lenear regression, containing all important numbers
all_mice_vals_times(:,1)=[0 0 0 14 14 14 14 14 30 30 30 30 30 30 60 60 60 60];
all_mice_vals_times(:,2)=[mean_gcl_per_mouse_no_tam;mean_gcl_per_mouse_two_weeks;mean_gcl_per_mouse_one_month;mean_gcl_per_mouse_two_months];

x=all_mice_vals_times(:,1);
y=all_mice_vals_times(:,2);

b1 = x\y;

yCalc1 = b1*x;
scatter(x,y)
hold on
plot(x,yCalc1)
xlabel('Population of state')
ylabel('Fatal traffic accidents per state')
title('Linear Regression Relation Between Accidents & Population')
grid on

Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2)   %this is r squared
b1     %this is the slope- cells added per day
[rho,pval]=corr(all_mice_vals_times) 
pval(1,2)    %this is the p_value to report