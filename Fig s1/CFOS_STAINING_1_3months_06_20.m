

%open manually all values
%they are in a folder on the desktop-right side named "measures"
%Changes the vectors of the means names to mean1...mean 6, and import as a
%column vector
clc
%close all
clear all

mice=[1 2 4 5 6 7];   %#3 was negative for genotype

for i=mice

path='F:\ALL DATA FOR PAPER HSN\ORGANIZED FOR PAPER\Fig s1\measures\';  %change path if needed~!
filename_inj=[path 'cfos' num2str(i) '_all_inj.xlsx'];
filename_non_inj=[path 'cfos' num2str(i) '_all_non_inj.xlsx'];
 data_inj = xlsread(filename_inj);
 data_non_inj = xlsread(filename_non_inj);
 INJ{i}=data_inj(:,2);
 NON_INJ{i}=data_non_inj(:,2);
end

   
thresh=2


cfos1_inj=cell2mat(INJ(1));        
cfos1_non_inj=cell2mat(NON_INJ(1));        

cfos2_inj=cell2mat(INJ(2));        
cfos2_non_inj=cell2mat(NON_INJ(2));      

cfos4_inj=cell2mat(INJ(4));        
cfos4_non_inj=cell2mat(NON_INJ(4));      

cfos5_inj=cell2mat(INJ(5));        
cfos5_non_inj=cell2mat(NON_INJ(5));      

cfos6_inj=cell2mat(INJ(6));        
cfos6_non_inj=cell2mat(NON_INJ(6));      

cfos7_inj=cell2mat(INJ(7));        
cfos7_non_inj=cell2mat(NON_INJ(7));      



%% normalize each slice by its 10 lowest measures:

%cfos1
cfos1_inj1=cfos1_inj(1:50);
cfos1_inj2=cfos1_inj(51:100);
cfos1_inj3=cfos1_inj(101:150);

cfos1_non_inj1=cfos1_non_inj(1:150);

a=(sort(cfos1_inj1))';
norm1=mean(a(1:10));

b=(sort(cfos1_inj2))';
norm2=mean(b(1:10));

c=(sort(cfos1_inj3))';
norm3=mean(c(1:10));

d=(sort(cfos1_non_inj1))';
norm4=mean(d(1:10));


cfos1_inj1_norm=cfos1_inj1./norm1;
cfos1_inj2_norm=cfos1_inj2./norm2;
cfos1_inj3_norm=cfos1_inj3./norm3;
cfos1_non_inj1_norm=cfos1_non_inj1./norm4;

all_norm_cfos1_inj=[cfos1_inj1_norm;cfos1_inj2_norm;cfos1_inj3_norm];
all_norm_cfos1_non_inj=[cfos1_non_inj1_norm];


%cfos2

cfos2_inj1=cfos2_inj(1:50);
cfos2_inj2=cfos2_inj(51:100);

cfos2_non_inj1=cfos2_non_inj(1:100);

a=(sort(cfos2_inj1))';
norm1=mean(a(1:10));

b=(sort(cfos2_inj2))';
norm2=mean(b(1:10));

c=(sort(cfos2_non_inj1))';
norm3=mean(c(1:10));


cfos2_inj1_norm=cfos2_inj1./norm1;
cfos2_inj2_norm=cfos2_inj2./norm2;
cfos2_non_inj1_norm=cfos2_non_inj1./norm3;


all_norm_cfos2_inj=[cfos2_inj1_norm;cfos2_inj2_norm];
all_norm_cfos2_non_inj=[cfos2_non_inj1_norm];


%cfos4

cfos4_inj1=cfos4_inj(1:70);
cfos4_inj2=cfos4_inj(71:120);
cfos4_inj3=cfos4_inj(121:150);

cfos4_non_inj1=cfos4_non_inj(1:50);
cfos4_non_inj2=cfos4_non_inj(51:150);

a=(sort(cfos4_inj1))';
norm1=mean(a(1:10));

b=(sort(cfos4_inj2))';
norm2=mean(b(1:10));

c=(sort(cfos4_inj3))';
norm3=mean(c(1:10));

d=(sort(cfos4_non_inj1))';
norm4=mean(d(1:10));

e=(sort(cfos4_non_inj2))';
norm5=mean(e(1:10));


cfos4_inj1_norm=cfos4_inj1./norm1;
cfos4_inj2_norm=cfos4_inj2./norm2;
cfos4_inj3_norm=cfos4_inj3./norm3;
cfos4_non_inj1_norm=cfos4_non_inj1./norm4;
cfos4_non_inj2_norm=cfos4_non_inj2./norm5;


all_norm_cfos4_inj=[cfos4_inj1_norm;cfos4_inj2_norm;cfos4_inj3_norm];
all_norm_cfos4_non_inj=[cfos4_non_inj1_norm;cfos4_non_inj2_norm];


%cfos 5- 200 cells each bulb
cfos5_inj1=cfos5_inj(1:47);            %47 cells
cfos5_inj2=cfos5_inj(48:108);          %61 cells
cfos5_inj3=cfos5_inj(109:150);         %42 cells
cfos5_inj4=cfos5_inj(151:200);         %50 cells

cfos5_non_inj1=cfos1_non_inj(1:50);      %50 cells each
cfos5_non_inj2=cfos5_non_inj(51:100);
cfos5_non_inj3=cfos5_non_inj(101:150);
cfos5_non_inj4=cfos5_non_inj(151:200);



a=(sort(cfos5_inj1))';
norm1=mean(a(1:10));

b=(sort(cfos5_inj2))';
norm2=mean(b(1:10));

c=(sort(cfos5_inj3))';
norm3=mean(c(1:10));


d=(sort(cfos5_inj4))';
norm4=mean(d(1:10));

e=(sort(cfos5_non_inj1))';
norm5=mean(e(1:10));

f=(sort(cfos5_non_inj2))';
norm6=mean(f(1:10));

g=(sort(cfos5_non_inj3))';
norm7=mean(g(1:10));

h=(sort(cfos5_non_inj4))';
norm8=mean(h(1:10));


cfos5_inj1_norm=cfos5_inj1./norm1;
cfos5_inj2_norm=cfos5_inj2./norm2;
cfos5_inj3_norm=cfos5_inj3./norm3;
cfos5_inj4_norm=cfos5_inj4./norm4;

cfos5_non_inj1_norm=cfos5_non_inj1./norm5;
cfos5_non_inj2_norm=cfos5_non_inj2./norm6;
cfos5_non_inj3_norm=cfos5_non_inj3./norm7;
cfos5_non_inj4_norm=cfos5_non_inj4./norm8;

all_norm_cfos5_inj=[cfos5_inj1_norm;cfos5_inj2_norm;cfos5_inj3_norm;cfos5_inj4_norm];
all_norm_cfos5_non_inj=[cfos5_non_inj1_norm;cfos5_non_inj2_norm;cfos5_non_inj3_norm;cfos5_non_inj4_norm];


%cfos6 150 cells
cfos6_inj1=cfos6_inj(1:28);            %28 cells
cfos6_inj2=cfos6_inj(29:70);          %42 cells
cfos6_inj3=cfos6_inj(71:120);         %50 cells
cfos6_inj4=cfos6_inj(121:150);         %30 cells


cfos6_non_inj1=cfos1_non_inj(1:60);      %60 cells 
cfos6_non_inj2=cfos6_non_inj(61:150);     %90 cells 
% cfos6_non_inj3=cfos6_non_inj(101:190);     %90 cells 



a=(sort(cfos6_inj1))';
norm1=mean(a(1:10));

b=(sort(cfos6_inj2))';
norm2=mean(b(1:10));

c=(sort(cfos6_inj3))';
norm3=mean(c(1:10));


d=(sort(cfos6_inj4))';
norm4=mean(d(1:10));

% e=(sort(cfos6_inj5))';   was double- thus not in use #2 originally
% norm5=mean(e(1:10));

f=(sort(cfos6_non_inj1))';
norm6=mean(f(1:10));

g=(sort(cfos6_non_inj2))';
norm7=mean(g(1:10));
% 
% h=(sort(cfos6_non_inj3))';
% norm8=mean(h(1:10));


cfos6_inj1_norm=cfos6_inj1./norm1;
cfos6_inj2_norm=cfos6_inj2./norm2;
cfos6_inj3_norm=cfos6_inj3./norm3;
cfos6_inj4_norm=cfos6_inj4./norm4;
% cfos6_inj5_norm=cfos6_inj5./norm5;

cfos6_non_inj1_norm=cfos6_non_inj1./norm6;
cfos6_non_inj2_norm=cfos6_non_inj2./norm7;
% cfos6_non_inj3_norm=cfos6_non_inj3./norm8;


all_norm_cfos6_inj=[cfos6_inj1_norm;cfos6_inj2_norm;cfos6_inj3_norm;cfos6_inj4_norm]; %;cfos6_inj5_norm];
all_norm_cfos6_non_inj=[cfos6_non_inj1_norm;cfos6_non_inj2_norm]; %;cfos6_non_inj3_norm];


%cfos7 180 cells

cfos7_inj1=cfos7_inj(1:60); %60 cells
cfos7_inj2=cfos7_inj(61:110); %50 cells
cfos7_inj3=cfos7_inj(111:180); %70 cells

cfos7_non_inj1=cfos7_non_inj(1:60); %60 cells
cfos7_non_inj2=cfos7_non_inj(61:120); %60 cells
cfos7_non_inj3=cfos7_non_inj(121:180); %60 cells


a=(sort(cfos7_inj1))';
norm1=mean(a(1:10));

b=(sort(cfos7_inj2))';
norm2=mean(b(1:10));

c=(sort(cfos7_inj3))';
norm3=mean(c(1:10));

d=(sort(cfos7_non_inj1))';
norm4=mean(d(1:10));

e=(sort(cfos7_non_inj2))';
norm5=mean(e(1:10));

f=(sort(cfos7_non_inj3))';
norm6=mean(f(1:10));



cfos7_inj1_norm=cfos7_inj1./norm1;
cfos7_inj2_norm=cfos7_inj2./norm2;
cfos7_inj3_norm=cfos7_inj3./norm3;
cfos7_non_inj1_norm=cfos7_non_inj1./norm4;
cfos7_non_inj2_norm=cfos7_non_inj2./norm5;
cfos7_non_inj3_norm=cfos7_non_inj3./norm6;


all_norm_cfos7_inj=[cfos7_inj1_norm;cfos7_inj2_norm;cfos7_inj3_norm];
all_norm_cfos7_non_inj=[cfos7_non_inj1_norm;cfos7_non_inj2_norm;cfos7_non_inj3_norm];

%%
figure
hold on
xlim([0 length(all_norm_cfos1_non_inj)])
bar(all_norm_cfos1_non_inj,'k')
bar(all_norm_cfos1_inj,'r')
xlabel('Cell #')
ylabel('Normalized c-Fos')
title('Mouse # 1')
legend('non injected bulb', 'injected bulb')
ylim([0 6])

x=ones(length(all_norm_cfos1_non_inj),1);
y=x.*thresh;
plot(1:length(y),y,'b--','linewidth',2)

figure
hold on
xlim([0 length(all_norm_cfos2_non_inj)])
bar(all_norm_cfos2_non_inj,'k')
bar(all_norm_cfos2_inj,'r')
xlabel('Cell #')
ylabel('Normalized c-Fos')
title('Mouse # 2')
legend('non injected bulb', 'injected bulb')
ylim([0 6])

x=ones(length(all_norm_cfos2_non_inj),1);
y=x.*thresh;
plot(1:length(y),y,'b--','linewidth',2)



figure
hold on
xlim([0 length(all_norm_cfos4_non_inj)])
bar(all_norm_cfos4_non_inj,'k')
bar(all_norm_cfos4_inj,'r')
xlabel('Cell #')
ylabel('Normalized c-Fos')
title('Mouse # 3')
legend('non injected bulb', 'injected bulb')
legend('non injected bulb', 'injected bulb')
ylim([0 6])

x=ones(length(all_norm_cfos4_non_inj),1);
y=x.*thresh;
plot(1:length(y),y,'b--','linewidth',2)

figure
hold on
xlim([0 length(all_norm_cfos5_non_inj)])
bar(all_norm_cfos5_non_inj,'k')
bar(all_norm_cfos5_inj,'r')
xlabel('Cell #')
ylabel('Normalized c-Fos')
title('Mouse # 4')                                     %this is cfos5, first of three 3 months post inj mice
legend('non injected bulb', 'injected bulb')
legend('non injected bulb', 'injected bulb')
ylim([0 6])

x=ones(length(all_norm_cfos5_non_inj),1);
y=x.*thresh;
plot(1:length(y),y,'b--','linewidth',2)

figure
hold on
xlim([0 length(all_norm_cfos6_non_inj)])
bar(all_norm_cfos6_non_inj,'k')
bar(all_norm_cfos6_inj,'r')
xlabel('Cell #')
ylabel('Normalized c-Fos')
title('Mouse # 5')
legend('non injected bulb', 'injected bulb')
legend('non injected bulb', 'injected bulb')
ylim([0 6])

x=ones(length(all_norm_cfos6_non_inj),1);
y=x.*thresh;
plot(1:length(y),y,'b--','linewidth',2)

figure
hold on
xlim([0 length(all_norm_cfos7_non_inj)])
bar(all_norm_cfos7_non_inj,'k')
bar(all_norm_cfos7_inj,'r')
xlabel('Cell #')
ylabel('Normalized c-Fos')
title('Mouse # 6')
legend('non injected bulb', 'injected bulb')
legend('non injected bulb', 'injected bulb')
ylim([0 6])

x=ones(length(all_norm_cfos7_non_inj),1);
y=x.*thresh;
plot(1:length(y),y,'b--','linewidth',2) %commented to plot graphs with the
%same length- unity in the bar widths
%xlim([0 150])
%plot(1:150,y(1:150),'b--','linewidth',2)

%% 
thresh=2

all_mice_inj_norm=[all_norm_cfos1_inj;all_norm_cfos2_inj;all_norm_cfos4_inj]
all_mice_non_inj_norm=[all_norm_cfos1_non_inj;all_norm_cfos2_non_inj;all_norm_cfos4_non_inj]

propinj1=sum(all_norm_cfos1_inj>=thresh)/length(all_norm_cfos1_inj)
propinj2=sum(all_norm_cfos2_inj>=thresh)/length(all_norm_cfos2_inj)
propinj4=sum(all_norm_cfos4_inj>=thresh)/length(all_norm_cfos4_inj)

sem_inj=std([propinj1 propinj2 propinj4])

prop_noninj1=sum(all_norm_cfos1_non_inj>=thresh)/length(all_norm_cfos1_non_inj)
prop_noninj2=sum(all_norm_cfos2_non_inj>=thresh)/length(all_norm_cfos2_non_inj)
prop_noninj4=sum(all_norm_cfos4_non_inj>=thresh)/length(all_norm_cfos4_non_inj)

sem_noninj=std([prop_noninj1 prop_noninj2 prop_noninj4])/sqrt(3)
sem_inj=std([propinj1 propinj2 propinj4])/sqrt(3)


means=[mean([prop_noninj1;prop_noninj2;prop_noninj4]) mean([propinj1;propinj2;propinj4])]
sems=[sem_noninj sem_inj]
figure
hold on
ylabel('Proportion of cells above threshold')

errorbar(means,sems,'.','linewidth', 2)
 ax=gca;
       ax.FontSize=14;
       ax.XTick=[];
      
       %% same for 3 months group:
       %% 
thresh=2

all_mice_inj_norm=[all_norm_cfos5_inj;all_norm_cfos6_inj;all_norm_cfos7_inj]
all_mice_non_inj_norm=[all_norm_cfos5_non_inj;all_norm_cfos6_non_inj;all_norm_cfos7_non_inj]

propinj5=sum(all_norm_cfos5_inj>=thresh)/length(all_norm_cfos5_inj)
propinj6=sum(all_norm_cfos6_inj>=thresh)/length(all_norm_cfos6_inj)
propinj7=sum(all_norm_cfos7_inj>=thresh)/length(all_norm_cfos7_inj)

sem_inj=std([propinj5 propinj6 propinj7])

prop_noninj5=sum(all_norm_cfos5_non_inj>=thresh)/length(all_norm_cfos5_non_inj)
prop_noninj6=sum(all_norm_cfos6_non_inj>=thresh)/length(all_norm_cfos6_non_inj)
prop_noninj7=sum(all_norm_cfos7_non_inj>=thresh)/length(all_norm_cfos7_non_inj)

sem_noninj=std([prop_noninj5 prop_noninj6 prop_noninj7])/sqrt(3)
sem_inj=std([propinj5 propinj6 propinj7])/sqrt(3)


means=[mean([prop_noninj5;prop_noninj6;prop_noninj7]) mean([propinj5;propinj6;propinj7])]
sems=[sem_noninj sem_inj]
figure
hold on
ylabel('Proportion of cells above threshold')

errorbar(means,sems,'.','linewidth', 2)
 ax=gca;
       ax.FontSize=14;
       ax.XTick=[];
       
       
       %% New figure for paper 6.2020
       
       
      
%        
       one_month_inj_prop=[propinj1 propinj2 propinj4];
       one_month_noninj_prop=[prop_noninj1 prop_noninj2 prop_noninj4];
       three_months_inj_prop= [propinj5;propinj6;propinj7];
       three_months_noninj_prop=[prop_noninj5;prop_noninj6;prop_noninj7];
       
        x=1:4;
       x_1month_inj = linspace(1.9,2.1,size(one_month_inj_prop,2))';
       x_1month_noninj = linspace(0.9,1.1,size(one_month_noninj_prop,2))';
       x_three_months_inj = linspace(3.9,4.1,size(three_months_inj_prop,1))';
       x_three_months_noninj = linspace(2.9,3.1,size(three_months_noninj_prop,1))';
       
       means_all=[mean(one_month_noninj_prop) mean(one_month_inj_prop) mean(three_months_noninj_prop) mean(three_months_inj_prop)];
       
       sem_one_month_inj_prop=std(one_month_inj_prop)/sqrt(3);
         sem_one_month_noninj_prop=std(one_month_noninj_prop)/sqrt(3);
           sem_three_months_inj_prop=std(three_months_inj_prop)/sqrt(3);
             sem_three_months_noninj_prop=std(three_months_noninj_prop)/sqrt(3);
       
        sem_all=[sem_one_month_noninj_prop sem_one_month_inj_prop sem_three_months_noninj_prop sem_three_months_inj_prop];
       
figure     
hold on
bar(x, means_all,0.3);
errorbar(x,means_all,sem_all,'.k','linewidth',2);
plot(x_1month_inj,one_month_inj_prop,'ok','MarkerFaceColor', 'k','MarkerSize',8);
plot(x_1month_noninj,one_month_noninj_prop,'ok','MarkerFaceColor', 'k','MarkerSize',8);
plot(x_three_months_inj,three_months_inj_prop,'ok','MarkerFaceColor', 'k','MarkerSize',8);
plot(x_three_months_noninj,three_months_noninj_prop,'ok','MarkerFaceColor', 'k','MarkerSize',8);

for i=1:3;
plot([x_1month_inj(i) x_1month_noninj(i)],[one_month_inj_prop(i),one_month_noninj_prop(i)],'k')
plot([x_three_months_inj(i) x_three_months_noninj(i)],[three_months_inj_prop(i),three_months_noninj_prop(i)],'k')

end

set(gca,'xTick',[1 2 3 4],'xTicklabel', {'1 Month injected','1 Month non injected','3 Month injected','3 Month non injected'})
ylabel('Proportion of c-Fos positive cells')
%% Chi test for eqaulity of porportions
   %% Statistical comparisopns based on Cells

           
    %% 1 month- inj vs non_inj
inj_total_1month=length(all_norm_cfos1_inj)+length(all_norm_cfos2_inj)+length(all_norm_cfos4_inj);          
inj_cfos_1month=sum(all_norm_cfos1_inj>=thresh)+sum(all_norm_cfos2_inj>=thresh)+sum(all_norm_cfos4_inj>=thresh);    
N1=inj_total_1month; n1=inj_cfos_1month;

non_inj_total_1month=length(all_norm_cfos1_non_inj)+length(all_norm_cfos2_non_inj)+length(all_norm_cfos4_non_inj);       
non_inj_cfos=sum(all_norm_cfos1_non_inj>=thresh)+sum(all_norm_cfos2_non_inj>=thresh)+sum(all_norm_cfos4_non_inj>=thresh);     
N2=non_inj_total_1month; n2=non_inj_cfos;
    
% Pooled estimate of proportion
       p0 = (n1+n2) / (N1+N2);
       % Expected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];
       chi2stat = sum((observed-expected).^2 ./ expected);
       p = 1 - chi2cdf(chi2stat,1)
       
     %% 3 months- inj vs non_inj
inj_total_3months=length(all_norm_cfos5_inj)+length(all_norm_cfos6_inj)+length(all_norm_cfos7_inj);          
inj_cfos_3month=sum(all_norm_cfos5_inj>=thresh)+sum(all_norm_cfos6_inj>=thresh)+sum(all_norm_cfos7_inj>=thresh);    
N1=inj_total_3months; n1=inj_cfos_3month;

non_inj_total_3months=length(all_norm_cfos5_non_inj)+length(all_norm_cfos6_non_inj)+length(all_norm_cfos7_non_inj);       
non_inj_cfos=sum(all_norm_cfos5_non_inj>=thresh)+sum(all_norm_cfos6_non_inj>=thresh)+sum(all_norm_cfos7_non_inj>=thresh);     
N2=non_inj_total_3months; n2=non_inj_cfos;
    
% Pooled estimate of proportion
       p0 = (n1+n2) / (N1+N2);
       % Expected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];
       chi2stat = sum((observed-expected).^2 ./ expected);
       p = 1 - chi2cdf(chi2stat,1);

     
 %% 1 month inj vs 3 months inj
 N1=inj_total_1month; n1=inj_cfos_1month;
 N2=inj_total_3months; n2=inj_cfos_3month;
 

 
 % Pooled estimate of proportion
       p0 = (n1+n2) / (N1+N2);
       % Expected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];
       chi2stat = sum((observed-expected).^2 ./ expected);
       p = 1 - chi2cdf(chi2stat,1);

       %% Statistical comparisopns based on mice
 mouse_group=[1 1 1 2 2 2 3 3 3 4 4 4];
 All_props=[one_month_noninj_prop';one_month_inj_prop';three_months_noninj_prop;three_months_inj_prop];
  
[p,tbl,stats] = kruskalwallis(All_props,mouse_group);
%[p,h]=signrank(one_month_noninj_prop,one_month_inj_prop);


