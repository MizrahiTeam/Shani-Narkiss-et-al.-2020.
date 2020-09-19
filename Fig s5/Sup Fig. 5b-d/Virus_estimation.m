%% THIS IS A MATLAB FILE TO ANALYSE virus specificity 9.2018
clc
clear all
close all

%%
filename = 'F:\ALL DATA FOR PAPER HSN\Virus20200601.xlsx';
datatable = readtable(filename);
datatable.norm=1000000./datatable.GCLArea_um_2_; %to move from microns to mm

% datatable.gcl_type_1_per_mm_squared=datatable.norm;
% datatable.gcl_type_2_per_mm_squared=datatable.norm;
datatable.norm_for_depth=1000./datatable.SliceUm;

filename2 = 'F:\ALL DATA FOR PAPER HSN\Virus20200601.xlsx';
datatable_pgl = readtable(filename2,'Sheet','PGL cells per mm^2'); 
datatable_pgl.norm=1000000./datatable_pgl.PGLArea_um_2_; %to move from microns to mm
datatable_pgl.norm_for_depth=1000./datatable_pgl.SliceUm;
for i=1:size(datatable,1);
datatable.gcl_type_1_per_mmsquared(i)=datatable.GCLType1(i)*datatable.norm(i);
datatable.gcl_type_2_per_mmsquared(i)=datatable.GCLType2(i)*datatable.norm(i);

datatable.gcl_type_1_per_mmcubed(i)=datatable.gcl_type_1_per_mmsquared(i)*datatable.norm_for_depth(i);
datatable.gcl_type_2_per_mmcubed(i)=datatable.gcl_type_2_per_mmsquared(i)*datatable.norm_for_depth(i);

end


% datatable_exp=datatable(1:9,:);
datatable_exp=datatable([1:9 23],:);
datatable_cont=datatable(10:19,:);
datatable_exp_2m=datatable(25:28,:);
datatable_exp_3m=datatable([20:22 29],:);


for i=1:size(datatable_pgl,1);
datatable.gcl_type_1_per_mmsquared(i)=datatable.GCLType1(i)*datatable.norm(i);
datatable.gcl_type_2_per_mmsquared(i)=datatable.GCLType2(i)*datatable.norm(i);

datatable_pgl.type_1_per_mmcubed(i)=datatable_pgl.type1PerMm_2(i)*datatable_pgl.norm_for_depth(i);
datatable_pgl.type_2_per_mmcubed(i)=datatable_pgl.type2PerMm_2(i)*datatable_pgl.norm_for_depth(i);

end


datatable_pgl_exp=datatable_pgl([1:9 23],:);
%%
exp_mean1=mean(datatable_exp.gcl_type_1_per_mmcubed);
exp_mean2=mean(datatable_exp.gcl_type_2_per_mmcubed);


exp_vec_hits=datatable_exp.gcl_type_1_per_mmcubed;
exp_vec_false=datatable_exp.gcl_type_2_per_mmcubed;

cont_mean1=mean(datatable_cont.gcl_type_1_per_mmcubed);
cont_mean2=mean(datatable_cont.gcl_type_2_per_mmcubed);

cont_vec_hits=datatable_cont.gcl_type_1_per_mmcubed;
cont_vec_false=datatable_cont.gcl_type_2_per_mmcubed;

exp_sem1=(std(datatable_exp.gcl_type_1_per_mmcubed))/sqrt(size(datatable_exp,1));
exp_sem2=(std(datatable_exp.gcl_type_2_per_mmcubed))/sqrt(size(datatable_exp,1));

cont_sem1=(std(datatable_cont.gcl_type_1_per_mmcubed))/sqrt(size(datatable_cont,1));
cont_sem2=(std(datatable_cont.gcl_type_2_per_mmcubed))/sqrt(size(datatable_cont,1));

exp_means_gcl=[exp_mean1,exp_mean2,0,0];
exp_sems_gcl=[exp_sem1,exp_sem2,0,0];

cont_means_gcl=[cont_mean1,cont_mean2];
cont_sems_gcl=[cont_sem1,cont_sem2];


exp_means_gcl_per_mm3=exp_means_gcl; 
exp_sems_gcl_per_mm3=exp_sems_gcl;
cont_means_gcl_per_mm3=cont_means_gcl;
cont_sems_gcl_per_mm3=cont_sems_gcl;

type1_means=[exp_means_gcl(1),cont_means_gcl(1)]
type2_means=[exp_means_gcl(2),cont_means_gcl(2)]

type1_sems=[exp_sems_gcl(1),cont_sems_gcl(1)]
type2_sems=[exp_sems_gcl(2),cont_sems_gcl(2)]

        

%% statistics

[p,h]=ranksum(exp_vec_false,cont_vec_false);
[p,h,stats]=ranksum(exp_vec_false,cont_vec_false); %test between the 10 mice of each group   %AS IS- used for paper
ratios=1-( exp_vec_false./(exp_vec_hits+exp_vec_false))  %what is the specificity for each exp mouse
mean(ratios)



 
%% GCL exp vs control- specificity- Fig 2C
x=1:2;
figure                          
bar(x,type1_means,0.3,'b')
hold on
errorbar(x,type1_means,type1_sems,'.k','linewidth',2)
%title('Granular Cells Layer count. N=15 mice, 45 slices','FontSize',25)

ylabel('Number of cells per mm^3')

hold on
bar(x,type2_means,0.3,'r')
hold on
errorbar(x,type2_means,type2_sems,'.k','linewidth',2)

axis square
box off

%Add individual points
x_exp = linspace(0.9,1.1,size(datatable_exp,1))
x_cont=linspace(1.9,2.1,size(datatable_cont,1))


plot(x_exp,datatable_exp.gcl_type_2_per_mmcubed,'or','MarkerSize',8)
plot(x_exp,datatable_exp.gcl_type_1_per_mmcubed,'or','MarkerFaceColor', 'b','MarkerSize',8)
plot(x_cont,datatable_cont.gcl_type_2_per_mmcubed,'or','MarkerSize',8)


%% PGL VS GCL specificity- exp. Fig S5b
exp_gcl_type1=datatable_exp.gcl_type_1_per_mmcubed;
exp_gcl_type2=datatable_exp.gcl_type_2_per_mmcubed


% 
exp_pgl_type1=datatable_pgl_exp.type_1_per_mmcubed;
exp_pgl_type2=datatable_pgl_exp.type_2_per_mmcubed;



x=1:2;

means_gcl_pgl_1=[mean(exp_gcl_type1) mean(exp_pgl_type1)];
means_gcl_pgl_2=[mean(exp_gcl_type2) mean(exp_pgl_type2)];

gcl_1_sem=(std(exp_gcl_type1))/sqrt(size(exp_gcl_type1,1));
pgl_1_sem=(std(exp_pgl_type1))/sqrt(size(exp_pgl_type1,1));

gcl_2_sem=(std(exp_gcl_type2))/sqrt(size(exp_gcl_type2,1));
pgl_2_sem=(std(exp_pgl_type2))/sqrt(size(exp_pgl_type2,1));

sem_gcl_pgl_1=[gcl_1_sem pgl_1_sem];
sem_gcl_pgl_2=[gcl_2_sem pgl_2_sem];

figure
hold on
bar(x,means_gcl_pgl_1,0.3,'m')
% bar(x,means_gcl_pgl_2,0.6,'b')
errorbar(x,means_gcl_pgl_1,sem_gcl_pgl_1,'.k','linewidth',2)


x_gcl = linspace(0.9,1.1,size(exp_gcl_type1,1))
x_pgl=linspace(1.9,2.1,size(exp_pgl_type1,1))



plot(x_gcl,exp_gcl_type1,'or','MarkerFaceColor', 'b','MarkerSize',8)
plot(x_pgl,exp_pgl_type1,'or','MarkerFaceColor', 'b','MarkerSize',8)

for i=1:size(exp_gcl_type1,1);
plot([x_gcl(i) x_pgl(i)],[exp_gcl_type1(i),exp_pgl_type1(i)],'k')
end

ylabel('Number of Myc+ & mRuby+ cells per mm^3')
set(gca,'xTick',[1 2],'xTicklabel', {'GCL','PGL'})

%% Fig S5c

exp_gcl_type1=datatable_exp.GCLType1;    %changed variables here to work with real number- not nirmalized per area
exp_gcl_type2=datatable_exp.GCLType2;



exp_pgl_type1=datatable_pgl_exp.PGLType1;
exp_pgl_type2=datatable_pgl_exp.PGLType2;

specif_gcl=(exp_gcl_type1+1)./(exp_gcl_type2+1)
specif_pgl=(exp_pgl_type1+1)./(exp_pgl_type2+1)
sem_spec_gcl=std(specif_gcl)/sqrt(size(specif_gcl,1)); 
sem_spec_pgl=std(specif_pgl)/sqrt(size(specif_gcl,1)); 

means_specif_gcl_pgl=[mean(specif_gcl) mean(specif_pgl)];
sems_specif_gcl_pgl=[sem_spec_gcl sem_spec_pgl];

x_gcl = linspace(0.9,1.1,size(exp_gcl_type1,1))
x_pgl=linspace(1.9,2.1,size(exp_pgl_type1,1))



figure
hold on
bar(x,means_specif_gcl_pgl,0.3)
errorbar(x,means_specif_gcl_pgl,sems_specif_gcl_pgl,'.k','linewidth',2)
plot(x_gcl,specif_gcl,'ok','MarkerFaceColor', 'k','MarkerSize',6)
plot(x_pgl,specif_pgl,'ok','MarkerFaceColor', 'k','MarkerSize',6)


for i=1:size(exp_gcl_type1,1);
plot([x_gcl(i) x_pgl(i)],[specif_gcl(i),specif_pgl(i)],'k')
end
%title('Granular Cells Layer count. N=15 mice, 45 slices','FontSize',25)

ylabel('Specificity Index')
set(gca,'xTick',[1 2],'xTicklabel', {'GCL','PGL'})


box off

%% GCL Myc+ mRuby+ over time Fig S5d


exp_gcl_type1_1month=datatable_exp.gcl_type_1_per_mmcubed;
exp_gcl_type1_2month=datatable_exp_2m.gcl_type_1_per_mmcubed;
exp_gcl_type1_3month=datatable_exp_3m.gcl_type_1_per_mmcubed;

x=1:3;

x_1month = linspace(0.9,1.1,size(exp_gcl_type1_1month,1))';
x_2month=linspace(1.9,2.1,size(exp_gcl_type1_2month,1))';
x_3month=linspace(2.9,3.1,size(exp_gcl_type1_3month,1))';

means_time=[mean(exp_gcl_type1_1month) mean(exp_gcl_type1_2month) mean(exp_gcl_type1_3month)];

sem_1month=std(exp_gcl_type1_1month)/sqrt(size(exp_gcl_type1_1month,1));
sem_2month=std(exp_gcl_type1_2month)/sqrt(size(exp_gcl_type1_2month,1));
sem_3month=std(exp_gcl_type1_3month)/sqrt(size(exp_gcl_type1_3month,1));

sems_time=[sem_1month sem_2month sem_3month];


figure
hold on
bar(x,means_time,0.3,'m');
errorbar(x,means_time,sems_time,'.k','linewidth',2);
plot(x_1month,exp_gcl_type1_1month,'or','MarkerFaceColor', 'b','MarkerSize',8);
plot(x_2month,exp_gcl_type1_2month,'or','MarkerFaceColor', 'b','MarkerSize',8);
plot(x_3month,exp_gcl_type1_3month,'or','MarkerFaceColor', 'b','MarkerSize',8);

set(gca,'xTick',[1 2 3],'xTicklabel', {'1 Month','2 Months','3 Months'})

mouse_group=[datatable_exp.MouseGroup;datatable_exp_2m.MouseGroup;datatable_exp_3m.MouseGroup];
GCLType1=[datatable_exp.gcl_type_1_per_mmcubed;datatable_exp_2m.gcl_type_1_per_mmcubed;datatable_exp_3m.gcl_type_1_per_mmcubed];
[p,tbl,stats] = kruskalwallis(GCLType1,mouse_group)



[h,p]=ttest2(exp_gcl_type1_1month,exp_gcl_type1_2month);
[h,p]=ttest2(exp_gcl_type1_1month,exp_gcl_type1_3month);
[h,p]=ttest2(exp_gcl_type1_2month,exp_gcl_type1_3month);

[p,tbl,stats] = anova1(GCLType1,mouse_group)
