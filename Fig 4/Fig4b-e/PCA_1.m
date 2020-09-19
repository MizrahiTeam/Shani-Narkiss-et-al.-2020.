clearvars -except  except Deltas_Dist_vec_norm_exp Deltas_Dist_vec_norm_cont
clc
%close all



for run_no=[1,2]

% CAL_PCA
% ('Valeraldehyde','Methyl propanoate','Ethyl acetate','Butyraldehyde','Blank','TMT','Female Pee', 'Male Pee', 'Peanut butter','Ethil tiglate','Propanal', 'Pups Bedding'}

clearvars -except  all_ST comp comp_SL DistMat DistMat_norm EucDis all_sds run_no Deltas_Dist_vec_norm_exp Deltas_Dist_vec_norm_cont
%clearvars -except    Deltas_Dist_vec_norm_cont
%important-CLEAR ALL between conditions but not runs
clc
%clear all
%close all


% chosen_time=56;
chosen_time=63;
%% Change only the condiotion
% 
% condition='Control';   anes %Fig s4
%  condition='Experiment'; anes %Fig s4
condition='awake'; %Fig 4
%condition='mice_awake_control'; %Fig 4
% run_no=1;    %before
% run_no=2;    %after


if condition(1)=='C'
load('Mice_table_Control.mat')
end

if condition(1)=='E'
load('Mice_table_Experiment.mat')
end

if condition(1)=='a'
load('Mice_table_awake.mat')
end

if condition(1)=='m'
load('Mice_table_mice_awake_control.mat')
end

number_of_conditions=size(mice_table,1)
all_data_cell=cell(number_of_conditions,1);

mice_num_bef=length(mice_table(1,:));
mice_num_af=length(mice_table(2,:));

num_of_fields=[mice_table{1,1}.fields_for_analysis.number_of_fields];
max_fields=max(num_of_fields)
all_data_cell{1}=cell(mice_num_bef,max_fields);
all_data_cell{2}=cell(mice_num_af,max_fields);

b=1;
c=1;
for a=1:number_of_conditions
    for b=1:mice_num_bef
        cells_per_field=mice_table{1,1}.fields_for_analysis(b).cells_per_field
        temp=1
        for c=1:num_of_fields(b)
            last_ind=temp+cells_per_field(c)-1
%         all_data_cell{a}{b,c}=mice_table{a,b}.SmoothDff(temp:last_ind,:,:,:)
        all_data_cell{a}{b,c}=mice_table{a,b}.SmoothDff_pca(temp:last_ind,:,:,:);
        all_ST{a}{b,c}=permute(all_data_cell{a}{b,c},[4 2 3 1]);          %creating the same structure as Amit's
        temp=last_ind+1
        end
    end
end


%%

 


%  n_sm = 3;
% win    = ones(n_sm,1);
% win    = win./length(win);
% % win=fir1( 2,0.2);


CellN=1;
for i=1:size(all_ST{run_no},1)
          for j=1:1:size(all_ST{run_no}(i,:),2)
              if isempty(all_ST{run_no}{i,j})==1
              continue
              end  
              
              for k=1:size(all_ST{run_no}{i,j},4) 
    
%Tria1 1
PCA_Mat_Valer1(CellN,:)=all_ST{run_no}{i,j}(:,1,1,k);    %note that val replaced the blank
PCA_Mat_Meth_prop1(CellN,:)=all_ST{run_no}{i,j}(:,1,2,k);
PCA_Mat_Eth_acet1(CellN,:)=all_ST{run_no}{i,j}(:,1,3,k);
PCA_Mat_Butyr1(CellN,:)=all_ST{run_no}{i,j}(:,1,4,k);
PCA_Mat_Blank1(CellN,:)=all_ST{run_no}{i,j}(:,1,5,k);
PCA_Mat_TMT1(CellN,:)=all_ST{run_no}{i,j}(:,1,6,k);
PCA_Mat_F_pee1(CellN,:)=all_ST{run_no}{i,j}(:,1,7,k);
PCA_Mat_M_pee1(CellN,:)=all_ST{run_no}{i,j}(:,1,8,k);
PCA_Mat_Pin_But1(CellN,:)=all_ST{run_no}{i,j}(:,1,9,k);
PCA_Mat_M_Eth_tig1(CellN,:)=all_ST{run_no}{i,j}(:,1,10,k);
PCA_Mat_Prop1(CellN,:)=all_ST{run_no}{i,j}(:,1,11,k);
PCA_Mat_M_Pups_bed1(CellN,:)=all_ST{run_no}{i,j}(:,1,12,k);

% PCA_Mat_But1(CellN,:)=filtfilt(win,1,all_ST{run_no}{i,j}{1,1}(:,1,3,k));

%Trial 2

PCA_Mat_Valer2(CellN,:)=all_ST{run_no}{i,j}(:,2,1,k);  
PCA_Mat_Meth_prop2(CellN,:)=all_ST{run_no}{i,j}(:,2,2,k);
PCA_Mat_Eth_acet2(CellN,:)=all_ST{run_no}{i,j}(:,2,3,k);
PCA_Mat_Butyr2(CellN,:)=all_ST{run_no}{i,j}(:,2,4,k);
PCA_Mat_Blank2(CellN,:)=all_ST{run_no}{i,j}(:,2,5,k);
PCA_Mat_TMT2(CellN,:)=all_ST{run_no}{i,j}(:,2,6,k);
PCA_Mat_F_pee2(CellN,:)=all_ST{run_no}{i,j}(:,2,7,k);
PCA_Mat_M_pee2(CellN,:)=all_ST{run_no}{i,j}(:,2,8,k);
PCA_Mat_Pin_But2(CellN,:)=all_ST{run_no}{i,j}(:,2,9,k);
PCA_Mat_M_Eth_tig2(CellN,:)=all_ST{run_no}{i,j}(:,2,10,k);
PCA_Mat_Prop2(CellN,:)=all_ST{run_no}{i,j}(:,2,11,k);
PCA_Mat_M_Pups_bed2(CellN,:)=all_ST{run_no}{i,j}(:,2,12,k);


%Trial 3
PCA_Mat_Valer3(CellN,:)=all_ST{run_no}{i,j}(:,3,1,k);  
PCA_Mat_Meth_prop3(CellN,:)=all_ST{run_no}{i,j}(:,3,2,k);
PCA_Mat_Eth_acet3(CellN,:)=all_ST{run_no}{i,j}(:,3,3,k);
PCA_Mat_Butyr3(CellN,:)=all_ST{run_no}{i,j}(:,3,4,k);
PCA_Mat_Blank3(CellN,:)=all_ST{run_no}{i,j}(:,3,5,k);
PCA_Mat_TMT3(CellN,:)=all_ST{run_no}{i,j}(:,3,6,k);
PCA_Mat_F_pee3(CellN,:)=all_ST{run_no}{i,j}(:,3,7,k);
PCA_Mat_M_pee3(CellN,:)=all_ST{run_no}{i,j}(:,3,8,k);
PCA_Mat_Pin_But3(CellN,:)=all_ST{run_no}{i,j}(:,3,9,k);
PCA_Mat_M_Eth_tig3(CellN,:)=all_ST{run_no}{i,j}(:,3,10,k);
PCA_Mat_Prop3(CellN,:)=all_ST{run_no}{i,j}(:,3,11,k);
PCA_Mat_M_Pups_bed3(CellN,:)=all_ST{run_no}{i,j}(:,3,12,k);

%Trial 4
PCA_Mat_Valer4(CellN,:)=all_ST{run_no}{i,j}(:,4,1,k);  
PCA_Mat_Meth_prop4(CellN,:)=all_ST{run_no}{i,j}(:,4,2,k);
PCA_Mat_Eth_acet4(CellN,:)=all_ST{run_no}{i,j}(:,4,3,k);
PCA_Mat_Butyr4(CellN,:)=all_ST{run_no}{i,j}(:,4,4,k);
PCA_Mat_Blank4(CellN,:)=all_ST{run_no}{i,j}(:,4,5,k);
PCA_Mat_TMT4(CellN,:)=all_ST{run_no}{i,j}(:,4,6,k);
PCA_Mat_F_pee4(CellN,:)=all_ST{run_no}{i,j}(:,4,7,k);
PCA_Mat_M_pee4(CellN,:)=all_ST{run_no}{i,j}(:,4,8,k);
PCA_Mat_Pin_But4(CellN,:)=all_ST{run_no}{i,j}(:,4,9,k);
PCA_Mat_M_Eth_tig4(CellN,:)=all_ST{run_no}{i,j}(:,4,10,k);
PCA_Mat_Prop4(CellN,:)=all_ST{run_no}{i,j}(:,4,11,k);
PCA_Mat_M_Pups_bed4(CellN,:)=all_ST{run_no}{i,j}(:,4,12,k);

%Trial 5
PCA_Mat_Valer5(CellN,:)=all_ST{run_no}{i,j}(:,5,1,k);  
PCA_Mat_Meth_prop5(CellN,:)=all_ST{run_no}{i,j}(:,5,2,k);
PCA_Mat_Eth_acet5(CellN,:)=all_ST{run_no}{i,j}(:,5,3,k);
PCA_Mat_Butyr5(CellN,:)=all_ST{run_no}{i,j}(:,5,4,k);
PCA_Mat_Blank5(CellN,:)=all_ST{run_no}{i,j}(:,5,5,k);
PCA_Mat_TMT5(CellN,:)=all_ST{run_no}{i,j}(:,5,6,k);
PCA_Mat_F_pee5(CellN,:)=all_ST{run_no}{i,j}(:,5,7,k);
PCA_Mat_M_pee5(CellN,:)=all_ST{run_no}{i,j}(:,5,8,k);
PCA_Mat_Pin_But5(CellN,:)=all_ST{run_no}{i,j}(:,5,9,k);
PCA_Mat_M_Eth_tig5(CellN,:)=all_ST{run_no}{i,j}(:,5,10,k);
PCA_Mat_Prop5(CellN,:)=all_ST{run_no}{i,j}(:,5,11,k);
PCA_Mat_M_Pups_bed5(CellN,:)=all_ST{run_no}{i,j}(:,5,12,k);

CellN=CellN+1;
k
              end
j
          end
i
end

PCA1=[PCA_Mat_Valer1';PCA_Mat_Meth_prop1';PCA_Mat_Eth_acet1';PCA_Mat_Butyr1';PCA_Mat_Blank1';PCA_Mat_TMT1';PCA_Mat_F_pee1';PCA_Mat_M_pee1';PCA_Mat_Pin_But1';PCA_Mat_M_Eth_tig1';PCA_Mat_Prop1';PCA_Mat_M_Pups_bed1'];
PCA2=[PCA_Mat_Valer2';PCA_Mat_Meth_prop2';PCA_Mat_Eth_acet2';PCA_Mat_Butyr2';PCA_Mat_Blank2';PCA_Mat_TMT2';PCA_Mat_F_pee2';PCA_Mat_M_pee2';PCA_Mat_Pin_But2';PCA_Mat_M_Eth_tig2';PCA_Mat_Prop2';PCA_Mat_M_Pups_bed2'];
PCA3=[PCA_Mat_Valer3';PCA_Mat_Meth_prop3';PCA_Mat_Eth_acet3';PCA_Mat_Butyr3';PCA_Mat_Blank3';PCA_Mat_TMT3';PCA_Mat_F_pee3';PCA_Mat_M_pee3';PCA_Mat_Pin_But3';PCA_Mat_M_Eth_tig3';PCA_Mat_Prop3';PCA_Mat_M_Pups_bed3'];
PCA4=[PCA_Mat_Valer4';PCA_Mat_Meth_prop4';PCA_Mat_Eth_acet4';PCA_Mat_Butyr4';PCA_Mat_Blank4';PCA_Mat_TMT4';PCA_Mat_F_pee4';PCA_Mat_M_pee4';PCA_Mat_Pin_But4';PCA_Mat_M_Eth_tig4';PCA_Mat_Prop4';PCA_Mat_M_Pups_bed4'];
PCA5=[PCA_Mat_Valer5';PCA_Mat_Meth_prop5';PCA_Mat_Eth_acet5';PCA_Mat_Butyr5';PCA_Mat_Blank5';PCA_Mat_TMT5';PCA_Mat_F_pee5';PCA_Mat_M_pee5';PCA_Mat_Pin_But5';PCA_Mat_M_Eth_tig5';PCA_Mat_Prop5';PCA_Mat_M_Pups_bed5'];

PCAall=[PCA1;PCA2;PCA3;PCA4;PCA5];
%%
x=PCAall;
%substracting the mean
for i = 1:size(x,2)
    x(:,i) = x(:,i) - mean(x(:,i));
end

% nCells = size(ResponseA,1)+size(ResponseB,1)+size(ResponseC,1);% If I
% want to normalize


c = x'*x;
[w,d] = eig(c);

w1 = w(:,end);
act1 = x*w1;

w2 = w(:,end-1);
act2 = x*w2;

w3 = w(:,end-2);
act3 = x*w3;

%%

%  EIG=eig(c);
% for i=1:length(EIG);
%     Scree(i)=EIG(i)/sum(EIG);
% end
% figure; plot(Scree)
%%
% deconcatinate back

%PCA1
WinShape=size(PCA_Mat_Valer1,2);

%Trial 1
PCA_1_1=reshape(act1(1:12*WinShape),WinShape,12);

PCA1Valer1=PCA_1_1(:,1);
PCA1Meth_prop1=PCA_1_1(:,2);
PCA1Eth_acet1=PCA_1_1(:,3);
PCA1Butyr1=PCA_1_1(:,4);
PCA1Blank1=PCA_1_1(:,5);
PCA1TMT1=PCA_1_1(:,6);
PCA1F_pee1=PCA_1_1(:,7);
PCA1M_pee1=PCA_1_1(:,8);
PCA1Pin_But1=PCA_1_1(:,9);
PCA1Eth_tig1=PCA_1_1(:,10);
PCA1Prop1=PCA_1_1(:,11);
PCA1Pups_bed1=PCA_1_1(:,12);

%Trial 2
PCA_1_2=reshape(act1(12*WinShape+1:2*12*WinShape),WinShape,12);

PCA1Valer2=PCA_1_2(:,1);
PCA1Meth_prop2=PCA_1_2(:,2);
PCA1Eth_acet2=PCA_1_2(:,3);
PCA1Butyr2=PCA_1_2(:,4);
PCA1Blank2=PCA_1_2(:,5);
PCA1TMT2=PCA_1_2(:,6);
PCA1F_pee2=PCA_1_2(:,7);
PCA1M_pee2=PCA_1_2(:,8);
PCA1Pin_But2=PCA_1_2(:,9);
PCA1Eth_tig2=PCA_1_2(:,10);
PCA1Prop2=PCA_1_2(:,11);
PCA1Pups_bed2=PCA_1_2(:,12);

%Trial 3
PCA_1_3=reshape(act1(2*12*WinShape+1:3*12*WinShape),WinShape,12);

PCA1Valer3=PCA_1_3(:,1);
PCA1Meth_prop3=PCA_1_3(:,2);
PCA1Eth_acet3=PCA_1_3(:,3);
PCA1Butyr3=PCA_1_3(:,4);
PCA1Blank3=PCA_1_3(:,5);
PCA1TMT3=PCA_1_3(:,6);
PCA1F_pee3=PCA_1_3(:,7);
PCA1M_pee3=PCA_1_3(:,8);
PCA1Pin_But3=PCA_1_3(:,9);
PCA1Eth_tig3=PCA_1_3(:,10);
PCA1Prop3=PCA_1_3(:,11);
PCA1Pups_bed3=PCA_1_3(:,12);


%Trial 4
PCA_1_4=reshape(act1(3*12*WinShape+1:4*12*WinShape),WinShape,12);

PCA1Valer4=PCA_1_4(:,1);
PCA1Meth_prop4=PCA_1_4(:,2);
PCA1Eth_acet4=PCA_1_4(:,3);
PCA1Butyr4=PCA_1_4(:,4);
PCA1Blank4=PCA_1_4(:,5);
PCA1TMT4=PCA_1_4(:,6);
PCA1F_pee4=PCA_1_4(:,7);
PCA1M_pee4=PCA_1_4(:,8);
PCA1Pin_But4=PCA_1_4(:,9);
PCA1Eth_tig4=PCA_1_4(:,10);
PCA1Prop4=PCA_1_4(:,11);
PCA1Pups_bed4=PCA_1_4(:,12);


%Trial 5
PCA_1_5=reshape(act1(4*12*WinShape+1:5*12*WinShape),WinShape,12);

PCA1Valer5=PCA_1_5(:,1);
PCA1Meth_prop5=PCA_1_5(:,2);
PCA1Eth_acet5=PCA_1_5(:,3);
PCA1Butyr5=PCA_1_5(:,4);
PCA1Blank5=PCA_1_5(:,5);
PCA1TMT5=PCA_1_5(:,6);
PCA1F_pee5=PCA_1_5(:,7);
PCA1M_pee5=PCA_1_5(:,8);
PCA1Pin_But5=PCA_1_5(:,9);
PCA1Eth_tig5=PCA_1_5(:,10);
PCA1Prop5=PCA_1_5(:,11);
PCA1Pups_bed5=PCA_1_5(:,12);


%PCA2

%Trial 1
PCA_2_1=reshape(act2(1:12*WinShape),WinShape,12);

PCA2Valer1=PCA_2_1(:,1);
PCA2Meth_prop1=PCA_2_1(:,2);
PCA2Eth_acet1=PCA_2_1(:,3);
PCA2Butyr1=PCA_2_1(:,4);
PCA2Blank1=PCA_2_1(:,5);
PCA2TMT1=PCA_2_1(:,6);
PCA2F_pee1=PCA_2_1(:,7);
PCA2M_pee1=PCA_2_1(:,8);
PCA2Pin_But1=PCA_2_1(:,9);
PCA2Eth_tig1=PCA_2_1(:,10);
PCA2Prop1=PCA_2_1(:,11);
PCA2Pups_bed1=PCA_2_1(:,12);

%Trial 2
PCA_2_2=reshape(act2(12*WinShape+1:2*12*WinShape),WinShape,12);

PCA2Valer2=PCA_2_2(:,1);
PCA2Meth_prop2=PCA_2_2(:,2);
PCA2Eth_acet2=PCA_2_2(:,3);
PCA2Butyr2=PCA_2_2(:,4);
PCA2Blank2=PCA_2_2(:,5);
PCA2TMT2=PCA_2_2(:,6);
PCA2F_pee2=PCA_2_2(:,7);
PCA2M_pee2=PCA_2_2(:,8);
PCA2Pin_But2=PCA_2_2(:,9);
PCA2Eth_tig2=PCA_2_2(:,10);
PCA2Prop2=PCA_2_2(:,11);
PCA2Pups_bed2=PCA_2_2(:,12);

%Trial 3
PCA_2_3=reshape(act2(2*12*WinShape+1:3*12*WinShape),WinShape,12);

PCA2Valer3=PCA_2_3(:,1);
PCA2Meth_prop3=PCA_2_3(:,2);
PCA2Eth_acet3=PCA_2_3(:,3);
PCA2Butyr3=PCA_2_3(:,4);
PCA2Blank3=PCA_2_3(:,5);
PCA2TMT3=PCA_2_3(:,6);
PCA2F_pee3=PCA_2_3(:,7);
PCA2M_pee3=PCA_2_3(:,8);
PCA2Pin_But3=PCA_2_3(:,9);
PCA2Eth_tig3=PCA_2_3(:,10);
PCA2Prop3=PCA_2_3(:,11);
PCA2Pups_bed3=PCA_2_3(:,12);


%Trial 4
PCA_2_4=reshape(act2(3*12*WinShape+1:4*12*WinShape),WinShape,12);

PCA2Valer4=PCA_2_4(:,1);
PCA2Meth_prop4=PCA_2_4(:,2);
PCA2Eth_acet4=PCA_2_4(:,3);
PCA2Butyr4=PCA_2_4(:,4);
PCA2Blank4=PCA_2_4(:,5);
PCA2TMT4=PCA_2_4(:,6);
PCA2F_pee4=PCA_2_4(:,7);
PCA2M_pee4=PCA_2_4(:,8);
PCA2Pin_But4=PCA_2_4(:,9);
PCA2Eth_tig4=PCA_2_4(:,10);
PCA2Prop4=PCA_2_4(:,11);
PCA2Pups_bed4=PCA_2_4(:,12);


%Trial 5
PCA_2_5=reshape(act2(4*12*WinShape+1:5*12*WinShape),WinShape,12);

PCA2Valer5=PCA_2_5(:,1);
PCA2Meth_prop5=PCA_2_5(:,2);
PCA2Eth_acet5=PCA_2_5(:,3);
PCA2Butyr5=PCA_2_5(:,4);
PCA2Blank5=PCA_2_5(:,5);
PCA2TMT5=PCA_2_5(:,6);
PCA2F_pee5=PCA_2_5(:,7);
PCA2M_pee5=PCA_2_5(:,8);
PCA2Pin_But5=PCA_2_5(:,9);
PCA2Eth_tig5=PCA_2_5(:,10);
PCA2Prop5=PCA_2_5(:,11);
PCA2Pups_bed5=PCA_2_5(:,12);

%PCA3

PCA_3_1=reshape(act3(1:12*WinShape),WinShape,12);

PCA3Valer1=PCA_3_1(:,1);
PCA3Meth_prop1=PCA_3_1(:,2);
PCA3Eth_acet1=PCA_3_1(:,3);
PCA3Butyr1=PCA_3_1(:,4);
PCA3Blank1=PCA_3_1(:,5);
PCA3TMT1=PCA_3_1(:,6);
PCA3F_pee1=PCA_3_1(:,7);
PCA3M_pee1=PCA_3_1(:,8);
PCA3Pin_But1=PCA_3_1(:,9);
PCA3Eth_tig1=PCA_3_1(:,10);
PCA3Prop1=PCA_3_1(:,11);
PCA3Pups_bed1=PCA_3_1(:,12);

%Trial 2
PCA_3_2=reshape(act3(12*WinShape+1:2*12*WinShape),WinShape,12);

PCA3Valer2=PCA_3_2(:,1);
PCA3Meth_prop2=PCA_3_2(:,2);
PCA3Eth_acet2=PCA_3_2(:,3);
PCA3Butyr2=PCA_3_2(:,4);
PCA3Blank2=PCA_3_2(:,5);
PCA3TMT2=PCA_3_2(:,6);
PCA3F_pee2=PCA_3_2(:,7);
PCA3M_pee2=PCA_3_2(:,8);
PCA3Pin_But2=PCA_3_2(:,9);
PCA3Eth_tig2=PCA_3_2(:,10);
PCA3Prop2=PCA_3_2(:,11);
PCA3Pups_bed2=PCA_3_2(:,12);

%Trial 3
PCA_3_3=reshape(act3(2*12*WinShape+1:3*12*WinShape),WinShape,12);

PCA3Valer3=PCA_3_3(:,1);
PCA3Meth_prop3=PCA_3_3(:,2);
PCA3Eth_acet3=PCA_3_3(:,3);
PCA3Butyr3=PCA_3_3(:,4);
PCA3Blank3=PCA_3_3(:,5);
PCA3TMT3=PCA_3_3(:,6);
PCA3F_pee3=PCA_3_3(:,7);
PCA3M_pee3=PCA_3_3(:,8);
PCA3Pin_But3=PCA_3_3(:,9);
PCA3Eth_tig3=PCA_3_3(:,10);
PCA3Prop3=PCA_3_3(:,11);
PCA3Pups_bed3=PCA_3_3(:,12);


%Trial 4
PCA_3_4=reshape(act3(3*12*WinShape+1:4*12*WinShape),WinShape,12);

PCA3Valer4=PCA_3_4(:,1);
PCA3Meth_prop4=PCA_3_4(:,2);
PCA3Eth_acet4=PCA_3_4(:,3);
PCA3Butyr4=PCA_3_4(:,4);
PCA3Blank4=PCA_3_4(:,5);
PCA3TMT4=PCA_3_4(:,6);
PCA3F_pee4=PCA_3_4(:,7);
PCA3M_pee4=PCA_3_4(:,8);
PCA3Pin_But4=PCA_3_4(:,9);
PCA3Eth_tig4=PCA_3_4(:,10);
PCA3Prop4=PCA_3_4(:,11);
PCA3Pups_bed4=PCA_3_4(:,12);


%Trial 5
PCA_3_5=reshape(act3(4*12*WinShape+1:5*12*WinShape),WinShape,12);

PCA3Valer5=PCA_3_5(:,1);
PCA3Meth_prop5=PCA_3_5(:,2);
PCA3Eth_acet5=PCA_3_5(:,3);
PCA3Butyr5=PCA_3_5(:,4);
PCA3Blank5=PCA_3_5(:,5);
PCA3TMT5=PCA_3_5(:,6);
PCA3F_pee5=PCA_3_5(:,7);
PCA3M_pee5=PCA_3_5(:,8);
PCA3Pin_But5=PCA_3_5(:,9);
PCA3Eth_tig5=PCA_3_5(:,10);
PCA3Prop5=PCA_3_5(:,11);
PCA3Pups_bed5=PCA_3_5(:,12);
%%
%%
% figure()
% hold all
% Num=1;
% % xlim([-0.4 0.4])
% % ylim([-0.4 0.4])
% % zlim([-0.4 0.4])
% % for i=1:size(PCA1Um1,1)%PCA1Um1
% for i=55;%:size(PCA1Um1,1)
% Igul=30;%(i)/10;
% plot3(PCA1Valer1(i,1),PCA2Valer1(i,1),PCA3Valer1(i,1),'ok','MarkerSize',Igul)
% plot3(PCA1Valer2(i,1),PCA2Valer2(i,1),PCA3Valer2(i,1),'ok','MarkerSize',Igul)
% plot3(PCA1Valer3(i,1),PCA2Valer3(i,1),PCA3Valer3(i,1),'ok','MarkerSize',Igul)
% plot3(PCA1Valer4(i,1),PCA2Valer4(i,1),PCA3Valer4(i,1),'ok','MarkerSize',Igul)
% plot3(PCA1Valer5(i,1),PCA2Valer5(i,1),PCA3Valer5(i,1),'ok','MarkerSize',Igul)
% 
% % 
% plot3(PCA1Meth_prop1(i,1),PCA2Meth_prop1(i,1),PCA3Meth_prop1(i,1),'or','MarkerSize',Igul)
% plot3(PCA1Meth_prop2(i,1),PCA2Meth_prop2(i,1),PCA3Meth_prop2(i,1),'or','MarkerSize',Igul)
% plot3(PCA1Meth_prop3(i,1),PCA2Meth_prop3(i,1),PCA3Meth_prop3(i,1),'or','MarkerSize',Igul)
% plot3(PCA1Meth_prop4(i,1),PCA2Meth_prop4(i,1),PCA3Meth_prop4(i,1),'or','MarkerSize',Igul)
% plot3(PCA1Meth_prop5(i,1),PCA2Meth_prop5(i,1),PCA3Meth_prop5(i,1),'or','MarkerSize',Igul)
% 
% % 
% plot3(PCA1Blank1(i,1),PCA2Blank1(i,1),PCA3Blank1(i,1),'.k','MarkerSize',Igul)
% plot3(PCA1Blank2(i,1),PCA2Blank2(i,1),PCA3Blank2(i,1),'.k','MarkerSize',Igul)
% plot3(PCA1Blank3(i,1),PCA2Blank3(i,1),PCA3Blank3(i,1),'.k','MarkerSize',Igul)
% plot3(PCA1Blank4(i,1),PCA2Blank4(i,1),PCA3Blank4(i,1),'.k','MarkerSize',Igul)
% plot3(PCA1Blank5(i,1),PCA2Blank5(i,1),PCA3Blank5(i,1),'.k','MarkerSize',Igul)
% 
% % 
% plot3(PCA1TMT1(i,1),PCA2TMT1(i,1),PCA3TMT1(i,1),'xb','MarkerSize',Igul)
% plot3(PCA1TMT2(i,1),PCA2TMT2(i,1),PCA3TMT2(i,1),'xb','MarkerSize',Igul)
% plot3(PCA1TMT3(i,1),PCA2TMT3(i,1),PCA3TMT3(i,1),'xb','MarkerSize',Igul)
% plot3(PCA1TMT4(i,1),PCA2TMT4(i,1),PCA3TMT4(i,1),'xb','MarkerSize',Igul)
% plot3(PCA1TMT5(i,1),PCA2TMT5(i,1),PCA3TMT5(i,1),'xb','MarkerSize',Igul)
% 
% % 
% plot3(PCA1F_pee1(i,1),PCA2F_pee1(i,1),PCA3F_pee1(i,1),'xm','MarkerSize',Igul)
% plot3(PCA1F_pee2(i,1),PCA2F_pee2(i,1),PCA3F_pee2(i,1),'xm','MarkerSize',Igul)
% plot3(PCA1F_pee3(i,1),PCA2F_pee3(i,1),PCA3F_pee3(i,1),'xm','MarkerSize',Igul)
% plot3(PCA1F_pee4(i,1),PCA2F_pee4(i,1),PCA3F_pee4(i,1),'xm','MarkerSize',Igul)
% plot3(PCA1F_pee5(i,1),PCA2F_pee5(i,1),PCA3F_pee5(i,1),'xm','MarkerSize',Igul)
% 
% %
% plot3(PCA1Eth_acet1(i,1),PCA2Eth_acet1(i,1),PCA3Eth_acet1(i,1),'ob','MarkerSize',Igul)
% plot3(PCA1Eth_acet2(i,1),PCA2Eth_acet2(i,1),PCA3Eth_acet2(i,1),'ob','MarkerSize',Igul)
% plot3(PCA1Eth_acet3(i,1),PCA2Eth_acet3(i,1),PCA3Eth_acet3(i,1),'ob','MarkerSize',Igul)
% plot3(PCA1Eth_acet4(i,1),PCA2Eth_acet4(i,1),PCA3Eth_acet4(i,1),'ob','MarkerSize',Igul)
% plot3(PCA1Eth_acet5(i,1),PCA2Eth_acet5(i,1),PCA3Eth_acet5(i,1),'ob','MarkerSize',Igul)
% 
% %
% plot3(PCA1Butyr1(i,1),PCA2Butyr1(i,1),PCA3Butyr1(i,1),'om','MarkerSize',Igul)
% plot3(PCA1Butyr2(i,1),PCA2Butyr2(i,1),PCA3Butyr2(i,1),'om','MarkerSize',Igul)
% plot3(PCA1Butyr3(i,1),PCA2Butyr3(i,1),PCA3Butyr3(i,1),'om','MarkerSize',Igul)
% plot3(PCA1Butyr4(i,1),PCA2Butyr4(i,1),PCA3Butyr4(i,1),'om','MarkerSize',Igul)
% plot3(PCA1Butyr5(i,1),PCA2Butyr5(i,1),PCA3Butyr5(i,1),'om','MarkerSize',Igul)
% 
% %
% plot3(PCA1M_pee1(i,1),PCA2M_pee1(i,1),PCA3M_pee1(i,1),'xr','MarkerSize',Igul)
% plot3(PCA1M_pee2(i,1),PCA2M_pee2(i,1),PCA3M_pee2(i,1),'xr','MarkerSize',Igul)
% plot3(PCA1M_pee3(i,1),PCA2M_pee3(i,1),PCA3M_pee3(i,1),'xr','MarkerSize',Igul)
% plot3(PCA1M_pee4(i,1),PCA2M_pee4(i,1),PCA3M_pee4(i,1),'xr','MarkerSize',Igul)
% plot3(PCA1M_pee5(i,1),PCA2M_pee5(i,1),PCA3M_pee5(i,1),'xr','MarkerSize',Igul)
% 
% %
% plot3(PCA1Pin_But1(i,1),PCA2Pin_But1(i,1),PCA3Pin_But1(i,1),'xk','MarkerSize',Igul)
% plot3(PCA1Pin_But2(i,1),PCA2Pin_But2(i,1),PCA3Pin_But2(i,1),'xk','MarkerSize',Igul)
% plot3(PCA1Pin_But3(i,1),PCA2Pin_But3(i,1),PCA3Pin_But3(i,1),'xk','MarkerSize',Igul)
% plot3(PCA1Pin_But4(i,1),PCA2Pin_But4(i,1),PCA3Pin_But4(i,1),'xk','MarkerSize',Igul)
% plot3(PCA1Pin_But5(i,1),PCA2Pin_But5(i,1),PCA3Pin_But5(i,1),'xk','MarkerSize',Igul)
% 
% %
% plot3(PCA1Eth_tig1(i,1),PCA2Eth_tig1(i,1),PCA3Eth_tig1(i,1),'og','MarkerSize',Igul)
% plot3(PCA1Eth_tig2(i,1),PCA2Eth_tig2(i,1),PCA3Eth_tig2(i,1),'og','MarkerSize',Igul)
% plot3(PCA1Eth_tig3(i,1),PCA2Eth_tig3(i,1),PCA3Eth_tig3(i,1),'og','MarkerSize',Igul)
% plot3(PCA1Eth_tig4(i,1),PCA2Eth_tig4(i,1),PCA3Eth_tig4(i,1),'og','MarkerSize',Igul)
% plot3(PCA1Eth_tig5(i,1),PCA2Eth_tig5(i,1),PCA3Eth_tig5(i,1),'og','MarkerSize',Igul)
% 
% %
% plot3(PCA1Prop1(i,1),PCA2Prop1(i,1),PCA3Prop1(i,1),'oy','MarkerSize',Igul)
% plot3(PCA1Prop2(i,1),PCA2Prop2(i,1),PCA3Prop2(i,1),'oy','MarkerSize',Igul)
% plot3(PCA1Prop3(i,1),PCA2Prop3(i,1),PCA3Prop3(i,1),'oy','MarkerSize',Igul)
% plot3(PCA1Prop4(i,1),PCA2Prop4(i,1),PCA3Prop4(i,1),'oy','MarkerSize',Igul)
% plot3(PCA1Prop5(i,1),PCA2Prop5(i,1),PCA3Prop5(i,1),'oy','MarkerSize',Igul)
% 
% %
% plot3(PCA1Pups_bed1(i,1),PCA2Pups_bed1(i,1),PCA3Pups_bed1(i,1),'xg','MarkerSize',Igul)
% plot3(PCA1Pups_bed2(i,1),PCA2Pups_bed2(i,1),PCA3Pups_bed2(i,1),'xg','MarkerSize',Igul)
% plot3(PCA1Pups_bed3(i,1),PCA2Pups_bed3(i,1),PCA3Pups_bed3(i,1),'xg','MarkerSize',Igul)
% plot3(PCA1Pups_bed4(i,1),PCA2Pups_bed4(i,1),PCA3Pups_bed4(i,1),'xg','MarkerSize',Igul)
% plot3(PCA1Pups_bed5(i,1),PCA2Pups_bed5(i,1),PCA3Pups_bed5(i,1),'xg','MarkerSize',Igul)
% 
% % str = num2str(round(i/7));
% % h = text(-0.2,-0.2, str);
% % pause
% %delete the annotation
% % delete(h)
% % % 
% end
% 
% 
% 
% 
 width = 2;     % Width in inches
 height = 2;    % Height in inches
 alw = 1;    % AxesLineWidth
 fsz = 10;      % Fontsize
 lw = 2;      % LineWidth
 msz = 8;       % MarkerSize
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
% % xlabel('PCA1')
% % ylabel('PCA2')
% % zlabel('PCA3')
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
% axis tight
% view([-20,-5,25]) 
%%



DisMapValer(:,1)=mean([PCA1Valer1(chosen_time,1),PCA1Valer2(chosen_time,1),PCA1Valer3(chosen_time,1),PCA1Valer4(chosen_time,1),PCA1Valer5(chosen_time,1)]);
DisMapValer(:,2)=mean([PCA2Valer1(chosen_time,1),PCA2Valer2(chosen_time,1),PCA2Valer3(chosen_time,1),PCA2Valer4(chosen_time,1),PCA2Valer5(chosen_time,1)]);
DisMapValer(:,3)=mean([PCA3Valer1(chosen_time,1),PCA3Valer2(chosen_time,1),PCA3Valer3(chosen_time,1),PCA3Valer4(chosen_time,1),PCA3Valer5(chosen_time,1)]);

DisMapMeth_prop(:,1)=mean([PCA1Meth_prop1(chosen_time,1),PCA1Meth_prop2(chosen_time,1),PCA1Meth_prop3(chosen_time,1),PCA1Meth_prop4(chosen_time,1),PCA1Meth_prop5(chosen_time,1)]);
DisMapMeth_prop(:,2)=mean([PCA2Meth_prop1(chosen_time,1),PCA2Meth_prop2(chosen_time,1),PCA2Meth_prop3(chosen_time,1),PCA2Meth_prop4(chosen_time,1),PCA2Meth_prop5(chosen_time,1)]);
DisMapMeth_prop(:,3)=mean([PCA3Meth_prop1(chosen_time,1),PCA3Meth_prop2(chosen_time,1),PCA3Meth_prop3(chosen_time,1),PCA3Meth_prop4(chosen_time,1),PCA3Meth_prop5(chosen_time,1)]);

DisMapEth_acet(:,1)=mean([PCA1Eth_acet1(chosen_time,1),PCA1Eth_acet2(chosen_time,1),PCA1Eth_acet3(chosen_time,1),PCA1Eth_acet4(chosen_time,1),PCA1Eth_acet5(chosen_time,1)]);
DisMapEth_acet(:,2)=mean([PCA2Eth_acet1(chosen_time,1),PCA2Eth_acet2(chosen_time,1),PCA2Eth_acet3(chosen_time,1),PCA2Eth_acet4(chosen_time,1),PCA2Eth_acet5(chosen_time,1)]);
DisMapEth_acet(:,3)=mean([PCA3Eth_acet1(chosen_time,1),PCA3Eth_acet2(chosen_time,1),PCA3Eth_acet3(chosen_time,1),PCA3Eth_acet4(chosen_time,1),PCA3Eth_acet5(chosen_time,1)]);

DisMapButyr(:,1)=mean([PCA1Butyr1(chosen_time,1),PCA1Butyr2(chosen_time,1),PCA1Butyr3(chosen_time,1),PCA1Butyr4(chosen_time,1),PCA1Butyr5(chosen_time,1)]);
DisMapButyr(:,2)=mean([PCA2Butyr1(chosen_time,1),PCA2Butyr2(chosen_time,1),PCA2Butyr3(chosen_time,1),PCA2Butyr4(chosen_time,1),PCA2Butyr5(chosen_time,1)]);
DisMapButyr(:,3)=mean([PCA3Butyr1(chosen_time,1),PCA3Butyr2(chosen_time,1),PCA3Butyr3(chosen_time,1),PCA3Butyr4(chosen_time,1),PCA3Butyr5(chosen_time,1)]);

DisMapBlank(:,1)=mean([PCA1Blank1(chosen_time,1),PCA1Blank2(chosen_time,1),PCA1Blank3(chosen_time,1),PCA1Blank4(chosen_time,1),PCA1Blank5(chosen_time,1)]);
DisMapBlank(:,2)=mean([PCA2Blank1(chosen_time,1),PCA2Blank2(chosen_time,1),PCA2Blank3(chosen_time,1),PCA2Blank4(chosen_time,1),PCA2Blank5(chosen_time,1)]);
DisMapBlank(:,3)=mean([PCA3Blank1(chosen_time,1),PCA3Blank2(chosen_time,1),PCA3Blank3(chosen_time,1),PCA3Blank4(chosen_time,1),PCA3Blank5(chosen_time,1)]);

DisMapTMT(:,1)=mean([PCA1TMT1(chosen_time,1),PCA1TMT2(chosen_time,1),PCA1TMT3(chosen_time,1),PCA1TMT4(chosen_time,1),PCA1TMT5(chosen_time,1)]);
DisMapTMT(:,2)=mean([PCA2TMT1(chosen_time,1),PCA2TMT2(chosen_time,1),PCA2TMT3(chosen_time,1),PCA2TMT4(chosen_time,1),PCA2TMT5(chosen_time,1)]);
DisMapTMT(:,3)=mean([PCA3TMT1(chosen_time,1),PCA3TMT2(chosen_time,1),PCA3TMT3(chosen_time,1),PCA3TMT4(chosen_time,1),PCA3TMT5(chosen_time,1)]);

DisMapF_pee(:,1)=mean([PCA1F_pee1(chosen_time,1),PCA1F_pee2(chosen_time,1),PCA1F_pee3(chosen_time,1),PCA1F_pee4(chosen_time,1),PCA1F_pee5(chosen_time,1)]);
DisMapF_pee(:,2)=mean([PCA2F_pee1(chosen_time,1),PCA2F_pee2(chosen_time,1),PCA2F_pee3(chosen_time,1),PCA2F_pee4(chosen_time,1),PCA2F_pee5(chosen_time,1)]);
DisMapF_pee(:,3)=mean([PCA3F_pee1(chosen_time,1),PCA3F_pee2(chosen_time,1),PCA3F_pee3(chosen_time,1),PCA3F_pee4(chosen_time,1),PCA3F_pee5(chosen_time,1)]);

DisMapM_pee(:,1)=mean([PCA1M_pee1(chosen_time,1),PCA1M_pee2(chosen_time,1),PCA1M_pee3(chosen_time,1),PCA1M_pee4(chosen_time,1),PCA1M_pee5(chosen_time,1)]);
DisMapM_pee(:,2)=mean([PCA2M_pee1(chosen_time,1),PCA2M_pee2(chosen_time,1),PCA2M_pee3(chosen_time,1),PCA2M_pee4(chosen_time,1),PCA2M_pee5(chosen_time,1)]);
DisMapM_pee(:,3)=mean([PCA3M_pee1(chosen_time,1),PCA3M_pee2(chosen_time,1),PCA3M_pee3(chosen_time,1),PCA3M_pee4(chosen_time,1),PCA3M_pee5(chosen_time,1)]);

DisMapPin_But(:,1)=mean([PCA1Pin_But1(chosen_time,1),PCA1Pin_But2(chosen_time,1),PCA1Pin_But3(chosen_time,1),PCA1Pin_But4(chosen_time,1),PCA1Pin_But5(chosen_time,1)]);
DisMapPin_But(:,2)=mean([PCA2Pin_But1(chosen_time,1),PCA2Pin_But2(chosen_time,1),PCA2Pin_But3(chosen_time,1),PCA2Pin_But4(chosen_time,1),PCA2Pin_But5(chosen_time,1)]);
DisMapPin_But(:,3)=mean([PCA3Pin_But1(chosen_time,1),PCA3Pin_But2(chosen_time,1),PCA3Pin_But3(chosen_time,1),PCA3Pin_But4(chosen_time,1),PCA3Pin_But5(chosen_time,1)]);

DisMapEth_tig(:,1)=mean([PCA1Eth_tig1(chosen_time,1),PCA1Eth_tig2(chosen_time,1),PCA1Eth_tig3(chosen_time,1),PCA1Eth_tig4(chosen_time,1),PCA1Eth_tig5(chosen_time,1)]);
DisMapEth_tig(:,2)=mean([PCA2Eth_tig1(chosen_time,1),PCA2Eth_tig2(chosen_time,1),PCA2Eth_tig3(chosen_time,1),PCA2Eth_tig4(chosen_time,1),PCA2Eth_tig5(chosen_time,1)]);
DisMapEth_tig(:,3)=mean([PCA3Eth_tig1(chosen_time,1),PCA3Eth_tig2(chosen_time,1),PCA3Eth_tig3(chosen_time,1),PCA3Eth_tig4(chosen_time,1),PCA3Eth_tig5(chosen_time,1)]);

DisMapProp(:,1)=mean([PCA1Prop1(chosen_time,1),PCA1Prop2(chosen_time,1),PCA1Prop3(chosen_time,1),PCA1Prop4(chosen_time,1),PCA1Prop5(chosen_time,1)]);
DisMapProp(:,2)=mean([PCA2Prop1(chosen_time,1),PCA2Prop2(chosen_time,1),PCA2Prop3(chosen_time,1),PCA2Prop4(chosen_time,1),PCA2Prop5(chosen_time,1)]);
DisMapProp(:,3)=mean([PCA3Prop1(chosen_time,1),PCA3Prop2(chosen_time,1),PCA3Prop3(chosen_time,1),PCA3Prop4(chosen_time,1),PCA3Prop5(chosen_time,1)]);

DisMapPups_bed(:,1)=mean([PCA1Pups_bed1(chosen_time,1),PCA1Pups_bed2(chosen_time,1),PCA1Pups_bed3(chosen_time,1),PCA1Pups_bed4(chosen_time,1),PCA1Pups_bed5(chosen_time,1)]);
DisMapPups_bed(:,2)=mean([PCA2Pups_bed1(chosen_time,1),PCA2Pups_bed2(chosen_time,1),PCA2Pups_bed3(chosen_time,1),PCA2Pups_bed4(chosen_time,1),PCA2Pups_bed5(chosen_time,1)]);
DisMapPups_bed(:,3)=mean([PCA3Pups_bed1(chosen_time,1),PCA3Pups_bed2(chosen_time,1),PCA3Pups_bed3(chosen_time,1),PCA3Pups_bed4(chosen_time,1),PCA3Pups_bed5(chosen_time,1)]);

%% adittion for d primes:
%valeraldehide
DisMapValer(:,1)=mean([PCA1Valer1(chosen_time,1),PCA1Valer2(chosen_time,1),PCA1Valer3(chosen_time,1),PCA1Valer4(chosen_time,1),PCA1Valer5(chosen_time,1)]);
DisMapValer(:,2)=mean([PCA2Valer1(chosen_time,1),PCA2Valer2(chosen_time,1),PCA2Valer3(chosen_time,1),PCA2Valer4(chosen_time,1),PCA2Valer5(chosen_time,1)]);
DisMapValer(:,3)=mean([PCA3Valer1(chosen_time,1),PCA3Valer2(chosen_time,1),PCA3Valer3(chosen_time,1),PCA3Valer4(chosen_time,1),PCA3Valer5(chosen_time,1)]);

Valer1_all_pcs=[PCA1Valer1(chosen_time,1),PCA2Valer1(chosen_time,1),PCA3Valer1(chosen_time,1)]
Valer2_all_pcs=[PCA1Valer2(chosen_time,1),PCA2Valer2(chosen_time,1),PCA3Valer2(chosen_time,1)]
Valer3_all_pcs=[PCA1Valer3(chosen_time,1),PCA2Valer3(chosen_time,1),PCA3Valer3(chosen_time,1)]
Valer4_all_pcs=[PCA1Valer4(chosen_time,1),PCA2Valer4(chosen_time,1),PCA3Valer4(chosen_time,1)]
Valer5_all_pcs=[PCA1Valer5(chosen_time,1),PCA2Valer5(chosen_time,1),PCA3Valer5(chosen_time,1)]
valer_rep_mat=[Valer1_all_pcs;Valer2_all_pcs;Valer3_all_pcs;Valer4_all_pcs;Valer5_all_pcs]

%meth_prop
DisMapMeth_prop(:,1)=mean([PCA1Meth_prop1(chosen_time,1),PCA1Meth_prop2(chosen_time,1),PCA1Meth_prop3(chosen_time,1),PCA1Meth_prop4(chosen_time,1),PCA1Meth_prop5(chosen_time,1)]);
DisMapMeth_prop(:,2)=mean([PCA2Meth_prop1(chosen_time,1),PCA2Meth_prop2(chosen_time,1),PCA2Meth_prop3(chosen_time,1),PCA2Meth_prop4(chosen_time,1),PCA2Meth_prop5(chosen_time,1)]);
DisMapMeth_prop(:,3)=mean([PCA3Meth_prop1(chosen_time,1),PCA3Meth_prop2(chosen_time,1),PCA3Meth_prop3(chosen_time,1),PCA3Meth_prop4(chosen_time,1),PCA3Meth_prop5(chosen_time,1)]);

Meth_prop1_all_pcs=[PCA1Meth_prop1(chosen_time,1),PCA2Meth_prop1(chosen_time,1),PCA3Meth_prop1(chosen_time,1)]
Meth_prop2_all_pcs=[PCA1Meth_prop2(chosen_time,1),PCA2Meth_prop2(chosen_time,1),PCA3Meth_prop2(chosen_time,1)]
Meth_prop3_all_pcs=[PCA1Meth_prop3(chosen_time,1),PCA2Meth_prop3(chosen_time,1),PCA3Meth_prop3(chosen_time,1)]
Meth_prop4_all_pcs=[PCA1Meth_prop4(chosen_time,1),PCA2Meth_prop4(chosen_time,1),PCA3Meth_prop4(chosen_time,1)]
Meth_prop5_all_pcs=[PCA1Meth_prop5(chosen_time,1),PCA2Meth_prop5(chosen_time,1),PCA3Meth_prop5(chosen_time,1)]
Meth_prop_rep_mat=[Meth_prop1_all_pcs;Meth_prop2_all_pcs;Meth_prop3_all_pcs;Meth_prop4_all_pcs;Meth_prop5_all_pcs]

%Eth acetate
DisMapEth_acet(:,1)=mean([PCA1Eth_acet1(chosen_time,1),PCA1Eth_acet2(chosen_time,1),PCA1Eth_acet3(chosen_time,1),PCA1Eth_acet4(chosen_time,1),PCA1Eth_acet5(chosen_time,1)]);
DisMapEth_acet(:,2)=mean([PCA2Eth_acet1(chosen_time,1),PCA2Eth_acet2(chosen_time,1),PCA2Eth_acet3(chosen_time,1),PCA2Eth_acet4(chosen_time,1),PCA2Eth_acet5(chosen_time,1)]);
DisMapEth_acet(:,3)=mean([PCA3Eth_acet1(chosen_time,1),PCA3Eth_acet2(chosen_time,1),PCA3Eth_acet3(chosen_time,1),PCA3Eth_acet4(chosen_time,1),PCA3Eth_acet5(chosen_time,1)]);

Eth_acet1_all_pcs=[PCA1Eth_acet1(chosen_time,1),PCA2Eth_acet1(chosen_time,1),PCA3Eth_acet1(chosen_time,1)]
Eth_acet2_all_pcs=[PCA1Eth_acet2(chosen_time,1),PCA2Eth_acet2(chosen_time,1),PCA3Eth_acet2(chosen_time,1)]
Eth_acet3_all_pcs=[PCA1Eth_acet3(chosen_time,1),PCA2Eth_acet3(chosen_time,1),PCA3Eth_acet3(chosen_time,1)]
Eth_acet4_all_pcs=[PCA1Eth_acet4(chosen_time,1),PCA2Eth_acet4(chosen_time,1),PCA3Eth_acet4(chosen_time,1)]
Eth_acet5_all_pcs=[PCA1Eth_acet5(chosen_time,1),PCA2Eth_acet5(chosen_time,1),PCA3Eth_acet5(chosen_time,1)]
Eth_acet_rep_mat=[Eth_acet1_all_pcs;Eth_acet2_all_pcs;Eth_acet3_all_pcs;Eth_acet4_all_pcs;Eth_acet5_all_pcs]

%Butyraldehyde

DisMapButyr(:,1)=mean([PCA1Butyr1(chosen_time,1),PCA1Butyr2(chosen_time,1),PCA1Butyr3(chosen_time,1),PCA1Butyr4(chosen_time,1),PCA1Butyr5(chosen_time,1)]);
DisMapButyr(:,2)=mean([PCA2Butyr1(chosen_time,1),PCA2Butyr2(chosen_time,1),PCA2Butyr3(chosen_time,1),PCA2Butyr4(chosen_time,1),PCA2Butyr5(chosen_time,1)]);
DisMapButyr(:,3)=mean([PCA3Butyr1(chosen_time,1),PCA3Butyr2(chosen_time,1),PCA3Butyr3(chosen_time,1),PCA3Butyr4(chosen_time,1),PCA3Butyr5(chosen_time,1)]);

Butyr1_all_pcs=[PCA1Butyr1(chosen_time,1),PCA2Butyr1(chosen_time,1),PCA3Butyr1(chosen_time,1)]
Butyr2_all_pcs=[PCA1Butyr2(chosen_time,1),PCA2Butyr2(chosen_time,1),PCA3Butyr2(chosen_time,1)]
Butyr3_all_pcs=[PCA1Butyr3(chosen_time,1),PCA2Butyr3(chosen_time,1),PCA3Butyr3(chosen_time,1)]
Butyr4_all_pcs=[PCA1Butyr4(chosen_time,1),PCA2Butyr4(chosen_time,1),PCA3Butyr4(chosen_time,1)]
Butyr5_all_pcs=[PCA1Butyr5(chosen_time,1),PCA2Butyr5(chosen_time,1),PCA3Butyr5(chosen_time,1)]
Butyr_rep_mat=[Butyr1_all_pcs;Butyr2_all_pcs;Butyr3_all_pcs;Butyr4_all_pcs;Butyr5_all_pcs]

%blank
DisMapBlank(:,1)=mean([PCA1Blank1(chosen_time,1),PCA1Blank2(chosen_time,1),PCA1Blank3(chosen_time,1),PCA1Blank4(chosen_time,1),PCA1Blank5(chosen_time,1)]);
DisMapBlank(:,2)=mean([PCA2Blank1(chosen_time,1),PCA2Blank2(chosen_time,1),PCA2Blank3(chosen_time,1),PCA2Blank4(chosen_time,1),PCA2Blank5(chosen_time,1)]);
DisMapBlank(:,3)=mean([PCA3Blank1(chosen_time,1),PCA3Blank2(chosen_time,1),PCA3Blank3(chosen_time,1),PCA3Blank4(chosen_time,1),PCA3Blank5(chosen_time,1)]);
Blank1_all_pcs=[PCA1Blank1(chosen_time,1),PCA2Blank1(chosen_time,1),PCA3Blank1(chosen_time,1)]
Blank2_all_pcs=[PCA1Blank2(chosen_time,1),PCA2Blank2(chosen_time,1),PCA3Blank2(chosen_time,1)]
Blank3_all_pcs=[PCA1Blank3(chosen_time,1),PCA2Blank3(chosen_time,1),PCA3Blank3(chosen_time,1)]
Blank4_all_pcs=[PCA1Blank4(chosen_time,1),PCA2Blank4(chosen_time,1),PCA3Blank4(chosen_time,1)]
Blank5_all_pcs=[PCA1Blank5(chosen_time,1),PCA2Blank5(chosen_time,1),PCA3Blank5(chosen_time,1)]
Blank_rep_mat=[Blank1_all_pcs;Blank2_all_pcs;Blank3_all_pcs;Blank4_all_pcs;Blank5_all_pcs]


%TMT
DisMapTMT(:,1)=mean([PCA1TMT1(chosen_time,1),PCA1TMT2(chosen_time,1),PCA1TMT3(chosen_time,1),PCA1TMT4(chosen_time,1),PCA1TMT5(chosen_time,1)]);
DisMapTMT(:,2)=mean([PCA2TMT1(chosen_time,1),PCA2TMT2(chosen_time,1),PCA2TMT3(chosen_time,1),PCA2TMT4(chosen_time,1),PCA2TMT5(chosen_time,1)]);
DisMapTMT(:,3)=mean([PCA3TMT1(chosen_time,1),PCA3TMT2(chosen_time,1),PCA3TMT3(chosen_time,1),PCA3TMT4(chosen_time,1),PCA3TMT5(chosen_time,1)]);
TMT1_all_pcs=[PCA1TMT1(chosen_time,1),PCA2TMT1(chosen_time,1),PCA3TMT1(chosen_time,1)]
TMT2_all_pcs=[PCA1TMT2(chosen_time,1),PCA2TMT2(chosen_time,1),PCA3TMT2(chosen_time,1)]
TMT3_all_pcs=[PCA1TMT3(chosen_time,1),PCA2TMT3(chosen_time,1),PCA3TMT3(chosen_time,1)]
TMT4_all_pcs=[PCA1TMT4(chosen_time,1),PCA2TMT4(chosen_time,1),PCA3TMT4(chosen_time,1)]
TMT5_all_pcs=[PCA1TMT5(chosen_time,1),PCA2TMT5(chosen_time,1),PCA3TMT5(chosen_time,1)]
TMT_rep_mat=[TMT1_all_pcs;TMT2_all_pcs;TMT3_all_pcs;TMT4_all_pcs;TMT5_all_pcs]

%F pee
DisMapF_pee(:,1)=mean([PCA1F_pee1(chosen_time,1),PCA1F_pee2(chosen_time,1),PCA1F_pee3(chosen_time,1),PCA1F_pee4(chosen_time,1),PCA1F_pee5(chosen_time,1)]);
DisMapF_pee(:,2)=mean([PCA2F_pee1(chosen_time,1),PCA2F_pee2(chosen_time,1),PCA2F_pee3(chosen_time,1),PCA2F_pee4(chosen_time,1),PCA2F_pee5(chosen_time,1)]);
DisMapF_pee(:,3)=mean([PCA3F_pee1(chosen_time,1),PCA3F_pee2(chosen_time,1),PCA3F_pee3(chosen_time,1),PCA3F_pee4(chosen_time,1),PCA3F_pee5(chosen_time,1)]);
F_pee1_all_pcs=[PCA1F_pee1(chosen_time,1),PCA2F_pee1(chosen_time,1),PCA3F_pee1(chosen_time,1)]
F_pee2_all_pcs=[PCA1F_pee2(chosen_time,1),PCA2F_pee2(chosen_time,1),PCA3F_pee2(chosen_time,1)]
F_pee3_all_pcs=[PCA1F_pee3(chosen_time,1),PCA2F_pee3(chosen_time,1),PCA3F_pee3(chosen_time,1)]
F_pee4_all_pcs=[PCA1F_pee4(chosen_time,1),PCA2F_pee4(chosen_time,1),PCA3F_pee4(chosen_time,1)]
F_pee5_all_pcs=[PCA1F_pee5(chosen_time,1),PCA2F_pee5(chosen_time,1),PCA3F_pee5(chosen_time,1)]
F_pee_rep_mat=[F_pee1_all_pcs;F_pee2_all_pcs;F_pee3_all_pcs;F_pee4_all_pcs;F_pee5_all_pcs]

%M_pee
DisMapM_pee(:,1)=mean([PCA1M_pee1(chosen_time,1),PCA1M_pee2(chosen_time,1),PCA1M_pee3(chosen_time,1),PCA1M_pee4(chosen_time,1),PCA1M_pee5(chosen_time,1)]);
DisMapM_pee(:,2)=mean([PCA2M_pee1(chosen_time,1),PCA2M_pee2(chosen_time,1),PCA2M_pee3(chosen_time,1),PCA2M_pee4(chosen_time,1),PCA2M_pee5(chosen_time,1)]);
DisMapM_pee(:,3)=mean([PCA3M_pee1(chosen_time,1),PCA3M_pee2(chosen_time,1),PCA3M_pee3(chosen_time,1),PCA3M_pee4(chosen_time,1),PCA3M_pee5(chosen_time,1)]);
M_pee1_all_pcs=[PCA1M_pee1(chosen_time,1),PCA2M_pee1(chosen_time,1),PCA3M_pee1(chosen_time,1)]
M_pee2_all_pcs=[PCA1M_pee2(chosen_time,1),PCA2M_pee2(chosen_time,1),PCA3M_pee2(chosen_time,1)]
M_pee3_all_pcs=[PCA1M_pee3(chosen_time,1),PCA2M_pee3(chosen_time,1),PCA3M_pee3(chosen_time,1)]
M_pee4_all_pcs=[PCA1M_pee4(chosen_time,1),PCA2M_pee4(chosen_time,1),PCA3M_pee4(chosen_time,1)]
M_pee5_all_pcs=[PCA1M_pee5(chosen_time,1),PCA2M_pee5(chosen_time,1),PCA3M_pee5(chosen_time,1)]
M_pee_rep_mat=[M_pee1_all_pcs;M_pee2_all_pcs;M_pee3_all_pcs;M_pee4_all_pcs;M_pee5_all_pcs]

%pin_but
DisMapPin_But(:,1)=mean([PCA1Pin_But1(chosen_time,1),PCA1Pin_But2(chosen_time,1),PCA1Pin_But3(chosen_time,1),PCA1Pin_But4(chosen_time,1),PCA1Pin_But5(chosen_time,1)]);
DisMapPin_But(:,2)=mean([PCA2Pin_But1(chosen_time,1),PCA2Pin_But2(chosen_time,1),PCA2Pin_But3(chosen_time,1),PCA2Pin_But4(chosen_time,1),PCA2Pin_But5(chosen_time,1)]);
DisMapPin_But(:,3)=mean([PCA3Pin_But1(chosen_time,1),PCA3Pin_But2(chosen_time,1),PCA3Pin_But3(chosen_time,1),PCA3Pin_But4(chosen_time,1),PCA3Pin_But5(chosen_time,1)]);
Pin_But1_all_pcs=[PCA1Pin_But1(chosen_time,1),PCA2Pin_But1(chosen_time,1),PCA3Pin_But1(chosen_time,1)]
Pin_But2_all_pcs=[PCA1Pin_But2(chosen_time,1),PCA2Pin_But2(chosen_time,1),PCA3Pin_But2(chosen_time,1)]
Pin_But3_all_pcs=[PCA1Pin_But3(chosen_time,1),PCA2Pin_But3(chosen_time,1),PCA3Pin_But3(chosen_time,1)]
Pin_But4_all_pcs=[PCA1Pin_But4(chosen_time,1),PCA2Pin_But4(chosen_time,1),PCA3Pin_But4(chosen_time,1)]
Pin_But5_all_pcs=[PCA1Pin_But5(chosen_time,1),PCA2Pin_But5(chosen_time,1),PCA3Pin_But5(chosen_time,1)]
Pin_But_rep_mat=[Pin_But1_all_pcs;Pin_But2_all_pcs;Pin_But3_all_pcs;Pin_But4_all_pcs;Pin_But5_all_pcs]

%Eth_tig
DisMapEth_tig(:,1)=mean([PCA1Eth_tig1(chosen_time,1),PCA1Eth_tig2(chosen_time,1),PCA1Eth_tig3(chosen_time,1),PCA1Eth_tig4(chosen_time,1),PCA1Eth_tig5(chosen_time,1)]);
DisMapEth_tig(:,2)=mean([PCA2Eth_tig1(chosen_time,1),PCA2Eth_tig2(chosen_time,1),PCA2Eth_tig3(chosen_time,1),PCA2Eth_tig4(chosen_time,1),PCA2Eth_tig5(chosen_time,1)]);
DisMapEth_tig(:,3)=mean([PCA3Eth_tig1(chosen_time,1),PCA3Eth_tig2(chosen_time,1),PCA3Eth_tig3(chosen_time,1),PCA3Eth_tig4(chosen_time,1),PCA3Eth_tig5(chosen_time,1)]);
Eth_tig1_all_pcs=[PCA1Eth_tig1(chosen_time,1),PCA2Eth_tig1(chosen_time,1),PCA3Eth_tig1(chosen_time,1)]
Eth_tig2_all_pcs=[PCA1Eth_tig2(chosen_time,1),PCA2Eth_tig2(chosen_time,1),PCA3Eth_tig2(chosen_time,1)]
Eth_tig3_all_pcs=[PCA1Eth_tig3(chosen_time,1),PCA2Eth_tig3(chosen_time,1),PCA3Eth_tig3(chosen_time,1)]
Eth_tig4_all_pcs=[PCA1Eth_tig4(chosen_time,1),PCA2Eth_tig4(chosen_time,1),PCA3Eth_tig4(chosen_time,1)]
Eth_tig5_all_pcs=[PCA1Eth_tig5(chosen_time,1),PCA2Eth_tig5(chosen_time,1),PCA3Eth_tig5(chosen_time,1)]
Eth_tig_rep_mat=[Eth_tig1_all_pcs;Eth_tig2_all_pcs;Eth_tig3_all_pcs;Eth_tig4_all_pcs;Eth_tig5_all_pcs]

%prop
DisMapProp(:,1)=mean([PCA1Prop1(chosen_time,1),PCA1Prop2(chosen_time,1),PCA1Prop3(chosen_time,1),PCA1Prop4(chosen_time,1),PCA1Prop5(chosen_time,1)]);
DisMapProp(:,2)=mean([PCA2Prop1(chosen_time,1),PCA2Prop2(chosen_time,1),PCA2Prop3(chosen_time,1),PCA2Prop4(chosen_time,1),PCA2Prop5(chosen_time,1)]);
DisMapProp(:,3)=mean([PCA3Prop1(chosen_time,1),PCA3Prop2(chosen_time,1),PCA3Prop3(chosen_time,1),PCA3Prop4(chosen_time,1),PCA3Prop5(chosen_time,1)]);
Prop1_all_pcs=[PCA1Prop1(chosen_time,1),PCA2Prop1(chosen_time,1),PCA3Prop1(chosen_time,1)]
Prop2_all_pcs=[PCA1Prop2(chosen_time,1),PCA2Prop2(chosen_time,1),PCA3Prop2(chosen_time,1)]
Prop3_all_pcs=[PCA1Prop3(chosen_time,1),PCA2Prop3(chosen_time,1),PCA3Prop3(chosen_time,1)]
Prop4_all_pcs=[PCA1Prop4(chosen_time,1),PCA2Prop4(chosen_time,1),PCA3Prop4(chosen_time,1)]
Prop5_all_pcs=[PCA1Prop5(chosen_time,1),PCA2Prop5(chosen_time,1),PCA3Prop5(chosen_time,1)]
Prop_rep_mat=[Prop1_all_pcs;Prop2_all_pcs;Prop3_all_pcs;Prop4_all_pcs;Prop5_all_pcs]

%Pups_bed
DisMapPups_bed(:,1)=mean([PCA1Pups_bed1(chosen_time,1),PCA1Pups_bed2(chosen_time,1),PCA1Pups_bed3(chosen_time,1),PCA1Pups_bed4(chosen_time,1),PCA1Pups_bed5(chosen_time,1)]);
DisMapPups_bed(:,2)=mean([PCA2Pups_bed1(chosen_time,1),PCA2Pups_bed2(chosen_time,1),PCA2Pups_bed3(chosen_time,1),PCA2Pups_bed4(chosen_time,1),PCA2Pups_bed5(chosen_time,1)]);
DisMapPups_bed(:,3)=mean([PCA3Pups_bed1(chosen_time,1),PCA3Pups_bed2(chosen_time,1),PCA3Pups_bed3(chosen_time,1),PCA3Pups_bed4(chosen_time,1),PCA3Pups_bed5(chosen_time,1)]);
Pups_bed1_all_pcs=[PCA1Pups_bed1(chosen_time,1),PCA2Pups_bed1(chosen_time,1),PCA3Pups_bed1(chosen_time,1)]
Pups_bed2_all_pcs=[PCA1Pups_bed2(chosen_time,1),PCA2Pups_bed2(chosen_time,1),PCA3Pups_bed2(chosen_time,1)]
Pups_bed3_all_pcs=[PCA1Pups_bed3(chosen_time,1),PCA2Pups_bed3(chosen_time,1),PCA3Pups_bed3(chosen_time,1)]
Pups_bed4_all_pcs=[PCA1Pups_bed4(chosen_time,1),PCA2Pups_bed4(chosen_time,1),PCA3Pups_bed4(chosen_time,1)]
Pups_bed5_all_pcs=[PCA1Pups_bed5(chosen_time,1),PCA2Pups_bed5(chosen_time,1),PCA3Pups_bed5(chosen_time,1)]
Pups_bed_rep_mat=[Pups_bed1_all_pcs;Pups_bed2_all_pcs;Pups_bed3_all_pcs;Pups_bed4_all_pcs;Pups_bed5_all_pcs]
%%
% figure
% hold on
% plot3(DisMapValer(:,1),DisMapValer(:,2),DisMapValer(:,3),'ok','MarkerSize',Igul)
% plot3(DisMapMeth_prop(:,1),DisMapMeth_prop(:,2),DisMapMeth_prop(:,3),'or','MarkerSize',Igul)
% plot3(DisMapBlank(:,1),DisMapBlank(:,2),DisMapBlank(:,3),'.k','MarkerSize',Igul)
% plot3(DisMapEth_acet(:,1),DisMapEth_acet(:,2),DisMapEth_acet(:,3),'ob','MarkerSize',Igul)
% plot3(DisMapTMT(:,1),DisMapTMT(:,2),DisMapTMT(:,3),'xb','MarkerSize',Igul)
% plot3(DisMapF_pee(:,1),DisMapF_pee(:,2),DisMapF_pee(:,3),'xm','MarkerSize',Igul)
% plot3(DisMapM_pee(:,1),DisMapM_pee(:,2),DisMapM_pee(:,3),'xr','MarkerSize',Igul)
% plot3(DisMapButyr(:,1),DisMapButyr(:,2),DisMapButyr(:,3),'om','MarkerSize',Igul)
% plot3(DisMapPin_But(:,1),DisMapPin_But(:,2),DisMapPin_But(:,3),'xk','MarkerSize',Igul)
% plot3(DisMapEth_tig(:,1),DisMapEth_tig(:,2),DisMapEth_tig(:,3),'og','MarkerSize',Igul)
% plot3(DisMapProp(:,1),DisMapProp(:,2),DisMapProp(:,3),'oy','MarkerSize',Igul)
% plot3(DisMapPups_bed(:,1),DisMapPups_bed(:,2),DisMapPups_bed(:,3),'xg','MarkerSize',Igul)

%%
DistMat{run_no}(1,1)=sqrt(sum((DisMapBlank - DisMapBlank) .^ 2));
DistMat{run_no}(1,2)=sqrt(sum((DisMapBlank - DisMapValer) .^ 2));
DistMat{run_no}(1,3)=sqrt(sum((DisMapBlank - DisMapMeth_prop) .^ 2));
DistMat{run_no}(1,4)=sqrt(sum((DisMapBlank - DisMapEth_acet) .^ 2));
DistMat{run_no}(1,5)=sqrt(sum((DisMapBlank - DisMapButyr) .^ 2));
DistMat{run_no}(1,6)=sqrt(sum((DisMapBlank - DisMapEth_tig) .^ 2));
DistMat{run_no}(1,7)=sqrt(sum((DisMapBlank - DisMapProp) .^ 2));
DistMat{run_no}(1,8)=sqrt(sum((DisMapBlank - DisMapTMT) .^ 2));
DistMat{run_no}(1,9)=sqrt(sum((DisMapBlank - DisMapF_pee) .^ 2));
DistMat{run_no}(1,10)=sqrt(sum((DisMapBlank - DisMapM_pee) .^ 2));
DistMat{run_no}(1,11)=sqrt(sum((DisMapBlank - DisMapPin_But) .^ 2));
DistMat{run_no}(1,12)=sqrt(sum((DisMapBlank - DisMapPups_bed) .^ 2));
% 2
DistMat{run_no}(2,1)=sqrt(sum((DisMapValer - DisMapBlank) .^ 2));
DistMat{run_no}(2,2)=sqrt(sum((DisMapValer - DisMapValer) .^ 2));
DistMat{run_no}(2,3)=sqrt(sum((DisMapValer - DisMapMeth_prop) .^ 2));
DistMat{run_no}(2,4)=sqrt(sum((DisMapValer - DisMapEth_acet) .^ 2));
DistMat{run_no}(2,5)=sqrt(sum((DisMapValer - DisMapButyr) .^ 2));
DistMat{run_no}(2,6)=sqrt(sum((DisMapValer - DisMapEth_tig) .^ 2));
DistMat{run_no}(2,7)=sqrt(sum((DisMapValer - DisMapProp) .^ 2));
DistMat{run_no}(2,8)=sqrt(sum((DisMapValer - DisMapTMT) .^ 2));
DistMat{run_no}(2,9)=sqrt(sum((DisMapValer - DisMapF_pee) .^ 2));
DistMat{run_no}(2,10)=sqrt(sum((DisMapValer - DisMapM_pee) .^ 2));
DistMat{run_no}(2,11)=sqrt(sum((DisMapValer - DisMapPin_But) .^ 2));
DistMat{run_no}(2,12)=sqrt(sum((DisMapValer - DisMapPups_bed) .^ 2));

% 3
DistMat{run_no}(3,1)=sqrt(sum((DisMapMeth_prop - DisMapBlank) .^ 2));
DistMat{run_no}(3,2)=sqrt(sum((DisMapMeth_prop - DisMapValer) .^ 2));
DistMat{run_no}(3,3)=sqrt(sum((DisMapMeth_prop - DisMapMeth_prop) .^ 2));
DistMat{run_no}(3,4)=sqrt(sum((DisMapMeth_prop - DisMapEth_acet) .^ 2));
DistMat{run_no}(3,5)=sqrt(sum((DisMapMeth_prop - DisMapButyr) .^ 2));
DistMat{run_no}(3,6)=sqrt(sum((DisMapMeth_prop - DisMapEth_tig) .^ 2));
DistMat{run_no}(3,7)=sqrt(sum((DisMapMeth_prop - DisMapProp) .^ 2));
DistMat{run_no}(3,8)=sqrt(sum((DisMapMeth_prop - DisMapTMT) .^ 2));
DistMat{run_no}(3,9)=sqrt(sum((DisMapMeth_prop - DisMapF_pee) .^ 2));
DistMat{run_no}(3,10)=sqrt(sum((DisMapMeth_prop - DisMapM_pee) .^ 2));
DistMat{run_no}(3,11)=sqrt(sum((DisMapMeth_prop - DisMapPin_But) .^ 2));
DistMat{run_no}(3,12)=sqrt(sum((DisMapMeth_prop - DisMapPups_bed) .^ 2));

%4
DistMat{run_no}(4,1)=sqrt(sum((DisMapEth_acet - DisMapBlank) .^ 2));
DistMat{run_no}(4,2)=sqrt(sum((DisMapEth_acet - DisMapValer) .^ 2));
DistMat{run_no}(4,3)=sqrt(sum((DisMapEth_acet - DisMapMeth_prop) .^ 2));
DistMat{run_no}(4,4)=sqrt(sum((DisMapEth_acet - DisMapEth_acet) .^ 2));
DistMat{run_no}(4,5)=sqrt(sum((DisMapEth_acet - DisMapButyr) .^ 2));
DistMat{run_no}(4,6)=sqrt(sum((DisMapEth_acet - DisMapEth_tig) .^ 2));
DistMat{run_no}(4,7)=sqrt(sum((DisMapEth_acet - DisMapProp) .^ 2));
DistMat{run_no}(4,8)=sqrt(sum((DisMapEth_acet - DisMapTMT) .^ 2));
DistMat{run_no}(4,9)=sqrt(sum((DisMapEth_acet - DisMapF_pee) .^ 2));
DistMat{run_no}(4,10)=sqrt(sum((DisMapEth_acet - DisMapM_pee) .^ 2));
DistMat{run_no}(4,11)=sqrt(sum((DisMapEth_acet - DisMapPin_But) .^ 2));
DistMat{run_no}(4,12)=sqrt(sum((DisMapEth_acet - DisMapPups_bed) .^ 2));

%5
DistMat{run_no}(5,1)=sqrt(sum((DisMapButyr - DisMapBlank) .^ 2));
DistMat{run_no}(5,2)=sqrt(sum((DisMapButyr - DisMapValer) .^ 2));
DistMat{run_no}(5,3)=sqrt(sum((DisMapButyr - DisMapMeth_prop) .^ 2));
DistMat{run_no}(5,4)=sqrt(sum((DisMapButyr - DisMapEth_acet) .^ 2));
DistMat{run_no}(5,5)=sqrt(sum((DisMapButyr - DisMapButyr) .^ 2));
DistMat{run_no}(5,6)=sqrt(sum((DisMapButyr - DisMapEth_tig) .^ 2));
DistMat{run_no}(5,7)=sqrt(sum((DisMapButyr - DisMapProp) .^ 2));
DistMat{run_no}(5,8)=sqrt(sum((DisMapButyr - DisMapTMT) .^ 2));
DistMat{run_no}(5,9)=sqrt(sum((DisMapButyr - DisMapF_pee) .^ 2));
DistMat{run_no}(5,10)=sqrt(sum((DisMapButyr - DisMapM_pee) .^ 2));
DistMat{run_no}(5,11)=sqrt(sum((DisMapButyr - DisMapPin_But) .^ 2));
DistMat{run_no}(5,12)=sqrt(sum((DisMapButyr - DisMapPups_bed) .^ 2));

%6
DistMat{run_no}(6,1)=sqrt(sum((DisMapEth_tig - DisMapBlank) .^ 2));
DistMat{run_no}(6,2)=sqrt(sum((DisMapEth_tig - DisMapValer) .^ 2));
DistMat{run_no}(6,3)=sqrt(sum((DisMapEth_tig - DisMapMeth_prop) .^ 2));
DistMat{run_no}(6,4)=sqrt(sum((DisMapEth_tig - DisMapEth_acet) .^ 2));
DistMat{run_no}(6,5)=sqrt(sum((DisMapEth_tig - DisMapButyr) .^ 2));
DistMat{run_no}(6,6)=sqrt(sum((DisMapEth_tig - DisMapEth_tig) .^ 2));
DistMat{run_no}(6,7)=sqrt(sum((DisMapEth_tig - DisMapProp) .^ 2));
DistMat{run_no}(6,8)=sqrt(sum((DisMapEth_tig - DisMapTMT) .^ 2));
DistMat{run_no}(6,9)=sqrt(sum((DisMapEth_tig - DisMapF_pee) .^ 2));
DistMat{run_no}(6,10)=sqrt(sum((DisMapEth_tig - DisMapM_pee) .^ 2));
DistMat{run_no}(6,11)=sqrt(sum((DisMapEth_tig - DisMapPin_But) .^ 2));
DistMat{run_no}(6,12)=sqrt(sum((DisMapEth_tig - DisMapPups_bed) .^ 2));

%7

DistMat{run_no}(7,1)=sqrt(sum((DisMapProp - DisMapBlank) .^ 2));
DistMat{run_no}(7,2)=sqrt(sum((DisMapProp - DisMapValer) .^ 2));
DistMat{run_no}(7,3)=sqrt(sum((DisMapProp - DisMapMeth_prop) .^ 2));
DistMat{run_no}(7,4)=sqrt(sum((DisMapProp - DisMapEth_acet) .^ 2));
DistMat{run_no}(7,5)=sqrt(sum((DisMapProp - DisMapButyr) .^ 2));
DistMat{run_no}(7,6)=sqrt(sum((DisMapProp - DisMapEth_tig) .^ 2));
DistMat{run_no}(7,7)=sqrt(sum((DisMapProp - DisMapProp) .^ 2));
DistMat{run_no}(7,8)=sqrt(sum((DisMapProp - DisMapTMT) .^ 2));
DistMat{run_no}(7,9)=sqrt(sum((DisMapProp - DisMapF_pee) .^ 2));
DistMat{run_no}(7,10)=sqrt(sum((DisMapProp - DisMapM_pee) .^ 2));
DistMat{run_no}(7,11)=sqrt(sum((DisMapProp - DisMapPin_But) .^ 2));
DistMat{run_no}(7,12)=sqrt(sum((DisMapProp - DisMapPups_bed) .^ 2));

%8

DistMat{run_no}(8,1)=sqrt(sum((DisMapTMT - DisMapBlank) .^ 2));
DistMat{run_no}(8,2)=sqrt(sum((DisMapTMT - DisMapValer) .^ 2));
DistMat{run_no}(8,3)=sqrt(sum((DisMapTMT - DisMapMeth_prop) .^ 2));
DistMat{run_no}(8,4)=sqrt(sum((DisMapTMT - DisMapEth_acet) .^ 2));
DistMat{run_no}(8,5)=sqrt(sum((DisMapTMT - DisMapButyr) .^ 2));
DistMat{run_no}(8,6)=sqrt(sum((DisMapTMT - DisMapEth_tig) .^ 2));
DistMat{run_no}(8,7)=sqrt(sum((DisMapTMT - DisMapProp) .^ 2));
DistMat{run_no}(8,8)=sqrt(sum((DisMapTMT - DisMapTMT) .^ 2));
DistMat{run_no}(8,9)=sqrt(sum((DisMapTMT - DisMapF_pee) .^ 2));
DistMat{run_no}(8,10)=sqrt(sum((DisMapTMT - DisMapM_pee) .^ 2));
DistMat{run_no}(8,11)=sqrt(sum((DisMapTMT - DisMapPin_But) .^ 2));
DistMat{run_no}(8,12)=sqrt(sum((DisMapTMT - DisMapPups_bed) .^ 2));


%9

DistMat{run_no}(9,1)=sqrt(sum((DisMapF_pee - DisMapBlank) .^ 2));
DistMat{run_no}(9,2)=sqrt(sum((DisMapF_pee - DisMapValer) .^ 2));
DistMat{run_no}(9,3)=sqrt(sum((DisMapF_pee - DisMapMeth_prop) .^ 2));
DistMat{run_no}(9,4)=sqrt(sum((DisMapF_pee - DisMapEth_acet) .^ 2));
DistMat{run_no}(9,5)=sqrt(sum((DisMapF_pee - DisMapButyr) .^ 2));
DistMat{run_no}(9,6)=sqrt(sum((DisMapF_pee - DisMapEth_tig) .^ 2));
DistMat{run_no}(9,7)=sqrt(sum((DisMapF_pee - DisMapProp) .^ 2));
DistMat{run_no}(9,8)=sqrt(sum((DisMapF_pee - DisMapTMT) .^ 2));
DistMat{run_no}(9,9)=sqrt(sum((DisMapF_pee - DisMapF_pee) .^ 2));
DistMat{run_no}(9,10)=sqrt(sum((DisMapF_pee - DisMapM_pee) .^ 2));
DistMat{run_no}(9,11)=sqrt(sum((DisMapF_pee - DisMapPin_But) .^ 2));
DistMat{run_no}(9,12)=sqrt(sum((DisMapF_pee - DisMapPups_bed) .^ 2));


%10

DistMat{run_no}(10,1)=sqrt(sum((DisMapM_pee - DisMapBlank) .^ 2));
DistMat{run_no}(10,2)=sqrt(sum((DisMapM_pee - DisMapValer) .^ 2));
DistMat{run_no}(10,3)=sqrt(sum((DisMapM_pee - DisMapMeth_prop) .^ 2));
DistMat{run_no}(10,4)=sqrt(sum((DisMapM_pee - DisMapEth_acet) .^ 2));
DistMat{run_no}(10,5)=sqrt(sum((DisMapM_pee - DisMapButyr) .^ 2));
DistMat{run_no}(10,6)=sqrt(sum((DisMapM_pee - DisMapEth_tig) .^ 2));
DistMat{run_no}(10,7)=sqrt(sum((DisMapM_pee - DisMapProp) .^ 2));
DistMat{run_no}(10,8)=sqrt(sum((DisMapM_pee - DisMapTMT) .^ 2));
DistMat{run_no}(10,9)=sqrt(sum((DisMapM_pee - DisMapF_pee) .^ 2));
DistMat{run_no}(10,10)=sqrt(sum((DisMapM_pee - DisMapM_pee) .^ 2));
DistMat{run_no}(10,11)=sqrt(sum((DisMapM_pee - DisMapPin_But) .^ 2));
DistMat{run_no}(10,12)=sqrt(sum((DisMapM_pee - DisMapPups_bed) .^ 2));

%11

DistMat{run_no}(11,1)=sqrt(sum((DisMapPin_But - DisMapBlank) .^ 2));
DistMat{run_no}(11,2)=sqrt(sum((DisMapPin_But - DisMapValer) .^ 2));
DistMat{run_no}(11,3)=sqrt(sum((DisMapPin_But - DisMapMeth_prop) .^ 2));
DistMat{run_no}(11,4)=sqrt(sum((DisMapPin_But - DisMapEth_acet) .^ 2));
DistMat{run_no}(11,5)=sqrt(sum((DisMapPin_But - DisMapButyr) .^ 2));
DistMat{run_no}(11,6)=sqrt(sum((DisMapPin_But - DisMapEth_tig) .^ 2));
DistMat{run_no}(11,7)=sqrt(sum((DisMapPin_But - DisMapProp) .^ 2));
DistMat{run_no}(11,8)=sqrt(sum((DisMapPin_But - DisMapTMT) .^ 2));
DistMat{run_no}(11,9)=sqrt(sum((DisMapPin_But - DisMapF_pee) .^ 2));
DistMat{run_no}(11,10)=sqrt(sum((DisMapPin_But - DisMapM_pee) .^ 2));
DistMat{run_no}(11,11)=sqrt(sum((DisMapPin_But - DisMapPin_But) .^ 2));
DistMat{run_no}(11,12)=sqrt(sum((DisMapPin_But - DisMapPups_bed) .^ 2));


%12

DistMat{run_no}(12,1)=sqrt(sum((DisMapPups_bed - DisMapBlank) .^ 2));
DistMat{run_no}(12,2)=sqrt(sum((DisMapPups_bed - DisMapValer) .^ 2));
DistMat{run_no}(12,3)=sqrt(sum((DisMapPups_bed - DisMapMeth_prop) .^ 2));
DistMat{run_no}(12,4)=sqrt(sum((DisMapPups_bed - DisMapEth_acet) .^ 2));
DistMat{run_no}(12,5)=sqrt(sum((DisMapPups_bed - DisMapButyr) .^ 2));
DistMat{run_no}(12,6)=sqrt(sum((DisMapPups_bed - DisMapEth_tig) .^ 2));
DistMat{run_no}(12,7)=sqrt(sum((DisMapPups_bed - DisMapProp) .^ 2));
DistMat{run_no}(12,8)=sqrt(sum((DisMapPups_bed - DisMapTMT) .^ 2));
DistMat{run_no}(12,9)=sqrt(sum((DisMapPups_bed - DisMapF_pee) .^ 2));
DistMat{run_no}(12,10)=sqrt(sum((DisMapPups_bed - DisMapM_pee) .^ 2));
DistMat{run_no}(12,11)=sqrt(sum((DisMapPups_bed - DisMapPin_But) .^ 2));
DistMat{run_no}(12,12)=sqrt(sum((DisMapPups_bed - DisMapPups_bed) .^ 2));




% %% 1
% DistMat{run_no}(1,1)=sqrt(sum((DisMapValer - DisMapValer) .^ 2));
% DistMat{run_no}(1,2)=sqrt(sum((DisMapValer - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(1,3)=sqrt(sum((DisMapValer - DisMapEth_acet) .^ 2));
% DistMat{run_no}(1,4)=sqrt(sum((DisMapValer - DisMapButyr) .^ 2));
% DistMat{run_no}(1,5)=sqrt(sum((DisMapValer - DisMapBlank) .^ 2));
% DistMat{run_no}(1,6)=sqrt(sum((DisMapValer - DisMapTMT) .^ 2));
% DistMat{run_no}(1,7)=sqrt(sum((DisMapValer - DisMapF_pee) .^ 2));
% DistMat{run_no}(1,8)=sqrt(sum((DisMapValer - DisMapM_pee) .^ 2));
% DistMat{run_no}(1,9)=sqrt(sum((DisMapValer - DisMapPin_But) .^ 2));
% DistMat{run_no}(1,10)=sqrt(sum((DisMapValer - DisMapEth_tig) .^ 2));
% DistMat{run_no}(1,11)=sqrt(sum((DisMapValer - DisMapProp) .^ 2));
% DistMat{run_no}(1,12)=sqrt(sum((DisMapValer - DisMapPups_bed) .^ 2));
% 
% % 2
% DistMat{run_no}(2,1)=sqrt(sum((DisMapMeth_prop - DisMapValer) .^ 2));
% DistMat{run_no}(2,2)=sqrt(sum((DisMapMeth_prop - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(2,3)=sqrt(sum((DisMapMeth_prop - DisMapEth_acet) .^ 2));
% DistMat{run_no}(2,4)=sqrt(sum((DisMapMeth_prop - DisMapButyr) .^ 2));
% DistMat{run_no}(2,5)=sqrt(sum((DisMapMeth_prop - DisMapBlank) .^ 2));
% DistMat{run_no}(2,6)=sqrt(sum((DisMapMeth_prop - DisMapTMT) .^ 2));
% DistMat{run_no}(2,7)=sqrt(sum((DisMapMeth_prop - DisMapF_pee) .^ 2));
% DistMat{run_no}(2,8)=sqrt(sum((DisMapMeth_prop - DisMapM_pee) .^ 2));
% DistMat{run_no}(2,9)=sqrt(sum((DisMapMeth_prop - DisMapPin_But) .^ 2));
% DistMat{run_no}(2,10)=sqrt(sum((DisMapMeth_prop - DisMapEth_tig) .^ 2));
% DistMat{run_no}(2,11)=sqrt(sum((DisMapMeth_prop - DisMapProp) .^ 2));
% DistMat{run_no}(2,12)=sqrt(sum((DisMapMeth_prop - DisMapPups_bed) .^ 2));
% 
% 
% % 3
% DistMat{run_no}(3,1)=sqrt(sum((DisMapEth_acet - DisMapValer) .^ 2));
% DistMat{run_no}(3,2)=sqrt(sum((DisMapEth_acet - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(3,3)=sqrt(sum((DisMapEth_acet - DisMapEth_acet) .^ 2));
% DistMat{run_no}(3,4)=sqrt(sum((DisMapEth_acet - DisMapButyr) .^ 2));
% DistMat{run_no}(3,5)=sqrt(sum((DisMapEth_acet - DisMapBlank) .^ 2));
% DistMat{run_no}(3,6)=sqrt(sum((DisMapEth_acet - DisMapTMT) .^ 2));
% DistMat{run_no}(3,7)=sqrt(sum((DisMapEth_acet - DisMapF_pee) .^ 2));
% DistMat{run_no}(3,8)=sqrt(sum((DisMapEth_acet - DisMapM_pee) .^ 2));
% DistMat{run_no}(3,9)=sqrt(sum((DisMapEth_acet - DisMapPin_But) .^ 2));
% DistMat{run_no}(3,10)=sqrt(sum((DisMapEth_acet - DisMapEth_tig) .^ 2));
% DistMat{run_no}(3,11)=sqrt(sum((DisMapEth_acet - DisMapProp) .^ 2));
% DistMat{run_no}(3,12)=sqrt(sum((DisMapEth_acet - DisMapPups_bed) .^ 2));
% 
% 
% %4
% DistMat{run_no}(4,1)=sqrt(sum((DisMapButyr - DisMapValer) .^ 2));
% DistMat{run_no}(4,2)=sqrt(sum((DisMapButyr - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(4,3)=sqrt(sum((DisMapButyr - DisMapEth_acet) .^ 2));
% DistMat{run_no}(4,4)=sqrt(sum((DisMapButyr - DisMapButyr) .^ 2));
% DistMat{run_no}(4,5)=sqrt(sum((DisMapButyr - DisMapBlank) .^ 2));
% DistMat{run_no}(4,6)=sqrt(sum((DisMapButyr - DisMapTMT) .^ 2));
% DistMat{run_no}(4,7)=sqrt(sum((DisMapButyr - DisMapF_pee) .^ 2));
% DistMat{run_no}(4,8)=sqrt(sum((DisMapButyr - DisMapM_pee) .^ 2));
% DistMat{run_no}(4,9)=sqrt(sum((DisMapButyr - DisMapPin_But) .^ 2));
% DistMat{run_no}(4,10)=sqrt(sum((DisMapButyr - DisMapEth_tig) .^ 2));
% DistMat{run_no}(4,11)=sqrt(sum((DisMapButyr - DisMapProp) .^ 2));
% DistMat{run_no}(4,12)=sqrt(sum((DisMapButyr - DisMapPups_bed) .^ 2));
% 
% %5
% DistMat{run_no}(5,1)=sqrt(sum((DisMapBlank - DisMapValer) .^ 2));
% DistMat{run_no}(5,2)=sqrt(sum((DisMapBlank - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(5,3)=sqrt(sum((DisMapBlank - DisMapEth_acet) .^ 2));
% DistMat{run_no}(5,4)=sqrt(sum((DisMapBlank - DisMapButyr) .^ 2));
% DistMat{run_no}(5,5)=sqrt(sum((DisMapBlank - DisMapBlank) .^ 2));
% DistMat{run_no}(5,6)=sqrt(sum((DisMapBlank - DisMapTMT) .^ 2));
% DistMat{run_no}(5,7)=sqrt(sum((DisMapBlank - DisMapF_pee) .^ 2));
% DistMat{run_no}(5,8)=sqrt(sum((DisMapBlank - DisMapM_pee) .^ 2));
% DistMat{run_no}(5,9)=sqrt(sum((DisMapBlank - DisMapPin_But) .^ 2));
% DistMat{run_no}(5,10)=sqrt(sum((DisMapBlank - DisMapEth_tig) .^ 2));
% DistMat{run_no}(5,11)=sqrt(sum((DisMapBlank - DisMapProp) .^ 2));
% DistMat{run_no}(5,12)=sqrt(sum((DisMapBlank - DisMapPups_bed) .^ 2));
% 
% %6
% DistMat{run_no}(6,1)=sqrt(sum((DisMapTMT - DisMapValer) .^ 2));
% DistMat{run_no}(6,2)=sqrt(sum((DisMapTMT - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(6,3)=sqrt(sum((DisMapTMT - DisMapEth_acet) .^ 2));
% DistMat{run_no}(6,4)=sqrt(sum((DisMapTMT - DisMapButyr) .^ 2));
% DistMat{run_no}(6,5)=sqrt(sum((DisMapTMT - DisMapBlank) .^ 2));
% DistMat{run_no}(6,6)=sqrt(sum((DisMapTMT - DisMapTMT) .^ 2));
% DistMat{run_no}(6,7)=sqrt(sum((DisMapTMT - DisMapF_pee) .^ 2));
% DistMat{run_no}(6,8)=sqrt(sum((DisMapTMT - DisMapM_pee) .^ 2));
% DistMat{run_no}(6,9)=sqrt(sum((DisMapTMT - DisMapPin_But) .^ 2));
% DistMat{run_no}(6,10)=sqrt(sum((DisMapTMT - DisMapEth_tig) .^ 2));
% DistMat{run_no}(6,11)=sqrt(sum((DisMapTMT - DisMapProp) .^ 2));
% DistMat{run_no}(6,12)=sqrt(sum((DisMapTMT - DisMapPups_bed) .^ 2));
% 
% 
% %7
% DistMat{run_no}(7,1)=sqrt(sum((DisMapF_pee - DisMapValer) .^ 2));
% DistMat{run_no}(7,2)=sqrt(sum((DisMapF_pee - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(7,3)=sqrt(sum((DisMapF_pee - DisMapEth_acet) .^ 2));
% DistMat{run_no}(7,4)=sqrt(sum((DisMapF_pee - DisMapButyr) .^ 2));
% DistMat{run_no}(7,5)=sqrt(sum((DisMapF_pee - DisMapBlank) .^ 2));
% DistMat{run_no}(7,6)=sqrt(sum((DisMapF_pee - DisMapTMT) .^ 2));
% DistMat{run_no}(7,7)=sqrt(sum((DisMapF_pee - DisMapF_pee) .^ 2));
% DistMat{run_no}(7,8)=sqrt(sum((DisMapF_pee - DisMapM_pee) .^ 2));
% DistMat{run_no}(7,9)=sqrt(sum((DisMapF_pee - DisMapPin_But) .^ 2));
% DistMat{run_no}(7,10)=sqrt(sum((DisMapF_pee - DisMapEth_tig) .^ 2));
% DistMat{run_no}(7,11)=sqrt(sum((DisMapF_pee - DisMapProp) .^ 2));
% DistMat{run_no}(7,12)=sqrt(sum((DisMapF_pee - DisMapPups_bed) .^ 2));
% 
% %8
% DistMat{run_no}(8,1)=sqrt(sum((DisMapM_pee - DisMapValer) .^ 2));
% DistMat{run_no}(8,2)=sqrt(sum((DisMapM_pee - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(8,3)=sqrt(sum((DisMapM_pee - DisMapEth_acet) .^ 2));
% DistMat{run_no}(8,4)=sqrt(sum((DisMapM_pee - DisMapButyr) .^ 2));
% DistMat{run_no}(8,5)=sqrt(sum((DisMapM_pee - DisMapBlank) .^ 2));
% DistMat{run_no}(8,6)=sqrt(sum((DisMapM_pee - DisMapTMT) .^ 2));
% DistMat{run_no}(8,7)=sqrt(sum((DisMapM_pee - DisMapF_pee) .^ 2));
% DistMat{run_no}(8,8)=sqrt(sum((DisMapM_pee - DisMapM_pee) .^ 2));
% DistMat{run_no}(8,9)=sqrt(sum((DisMapM_pee - DisMapPin_But) .^ 2));
% DistMat{run_no}(8,10)=sqrt(sum((DisMapM_pee - DisMapEth_tig) .^ 2));
% DistMat{run_no}(8,11)=sqrt(sum((DisMapM_pee - DisMapProp) .^ 2));
% DistMat{run_no}(8,12)=sqrt(sum((DisMapM_pee - DisMapPups_bed) .^ 2));
% 
% %9
% DistMat{run_no}(9,1)=sqrt(sum((DisMapPin_But - DisMapValer) .^ 2));
% DistMat{run_no}(9,2)=sqrt(sum((DisMapPin_But - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(9,3)=sqrt(sum((DisMapPin_But - DisMapEth_acet) .^ 2));
% DistMat{run_no}(9,4)=sqrt(sum((DisMapPin_But - DisMapButyr) .^ 2));
% DistMat{run_no}(9,5)=sqrt(sum((DisMapPin_But - DisMapBlank) .^ 2));
% DistMat{run_no}(9,6)=sqrt(sum((DisMapPin_But - DisMapTMT) .^ 2));
% DistMat{run_no}(9,7)=sqrt(sum((DisMapPin_But - DisMapF_pee) .^ 2));
% DistMat{run_no}(9,8)=sqrt(sum((DisMapPin_But - DisMapM_pee) .^ 2));
% DistMat{run_no}(9,9)=sqrt(sum((DisMapPin_But - DisMapPin_But) .^ 2));
% DistMat{run_no}(9,10)=sqrt(sum((DisMapPin_But - DisMapEth_tig) .^ 2));
% DistMat{run_no}(9,11)=sqrt(sum((DisMapPin_But - DisMapProp) .^ 2));
% DistMat{run_no}(9,12)=sqrt(sum((DisMapPin_But - DisMapPups_bed) .^ 2));
% 
% %10
% DistMat{run_no}(10,1)=sqrt(sum((DisMapEth_tig - DisMapValer) .^ 2));
% DistMat{run_no}(10,2)=sqrt(sum((DisMapEth_tig - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(10,3)=sqrt(sum((DisMapEth_tig - DisMapEth_acet) .^ 2));
% DistMat{run_no}(10,4)=sqrt(sum((DisMapEth_tig - DisMapButyr) .^ 2));
% DistMat{run_no}(10,5)=sqrt(sum((DisMapEth_tig - DisMapBlank) .^ 2));
% DistMat{run_no}(10,6)=sqrt(sum((DisMapEth_tig - DisMapTMT) .^ 2));
% DistMat{run_no}(10,7)=sqrt(sum((DisMapEth_tig - DisMapF_pee) .^ 2));
% DistMat{run_no}(10,8)=sqrt(sum((DisMapEth_tig - DisMapM_pee) .^ 2));
% DistMat{run_no}(10,9)=sqrt(sum((DisMapEth_tig - DisMapPin_But) .^ 2));
% DistMat{run_no}(10,10)=sqrt(sum((DisMapEth_tig - DisMapEth_tig) .^ 2));
% DistMat{run_no}(10,11)=sqrt(sum((DisMapEth_tig - DisMapProp) .^ 2));
% DistMat{run_no}(10,12)=sqrt(sum((DisMapEth_tig - DisMapPups_bed) .^ 2));
% 
% %11
% DistMat{run_no}(11,1)=sqrt(sum((DisMapProp - DisMapValer) .^ 2));
% DistMat{run_no}(11,2)=sqrt(sum((DisMapProp - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(11,3)=sqrt(sum((DisMapProp - DisMapEth_acet) .^ 2));
% DistMat{run_no}(11,4)=sqrt(sum((DisMapProp - DisMapButyr) .^ 2));
% DistMat{run_no}(11,5)=sqrt(sum((DisMapProp - DisMapBlank) .^ 2));
% DistMat{run_no}(11,6)=sqrt(sum((DisMapProp - DisMapTMT) .^ 2));
% DistMat{run_no}(11,7)=sqrt(sum((DisMapProp - DisMapF_pee) .^ 2));
% DistMat{run_no}(11,8)=sqrt(sum((DisMapProp - DisMapM_pee) .^ 2));
% DistMat{run_no}(11,9)=sqrt(sum((DisMapProp - DisMapPin_But) .^ 2));
% DistMat{run_no}(11,10)=sqrt(sum((DisMapProp - DisMapEth_tig) .^ 2));
% DistMat{run_no}(11,11)=sqrt(sum((DisMapProp - DisMapProp) .^ 2));
% DistMat{run_no}(11,12)=sqrt(sum((DisMapProp - DisMapPups_bed) .^ 2));
% 
% %12
% DistMat{run_no}(12,1)=sqrt(sum((DisMapPups_bed - DisMapValer) .^ 2));
% DistMat{run_no}(12,2)=sqrt(sum((DisMapPups_bed - DisMapMeth_prop) .^ 2));
% DistMat{run_no}(12,3)=sqrt(sum((DisMapPups_bed - DisMapEth_acet) .^ 2));
% DistMat{run_no}(12,4)=sqrt(sum((DisMapPups_bed - DisMapButyr) .^ 2));
% DistMat{run_no}(12,5)=sqrt(sum((DisMapPups_bed - DisMapBlank) .^ 2));
% DistMat{run_no}(12,6)=sqrt(sum((DisMapPups_bed - DisMapTMT) .^ 2));
% DistMat{run_no}(12,7)=sqrt(sum((DisMapPups_bed - DisMapF_pee) .^ 2));
% DistMat{run_no}(12,8)=sqrt(sum((DisMapPups_bed - DisMapM_pee) .^ 2));
% DistMat{run_no}(12,9)=sqrt(sum((DisMapPups_bed - DisMapPin_But) .^ 2));
% DistMat{run_no}(12,10)=sqrt(sum((DisMapPups_bed - DisMapEth_tig) .^ 2));
% DistMat{run_no}(12,11)=sqrt(sum((DisMapPups_bed - DisMapProp) .^ 2));
% DistMat{run_no}(12,12)=sqrt(sum((DisMapPups_bed - DisMapPups_bed) .^ 2));
%% noise estimation for normalized distmat
Valer_sd=[];
Meth_prop_sd=[];
Eth_acet_sd=[];
Butyr_sd=[];
Blank_sd=[];
F_pee_sd=[];
M_pee_sd=[];
Pin_But_sd=[];
TMT_sd=[];
Prop_sd=[];
Eth_tig_sd=[];
Pups_bed_sd=[]

i=1
for i=1:5
 Valer_sd(i)= sqrt(sum((DisMapValer - valer_rep_mat(i,:)) .^ 2))
  Meth_prop_sd(i)= sqrt(sum((DisMapMeth_prop - Meth_prop_rep_mat(i,:)) .^ 2))
   Eth_acet_sd(i)= sqrt(sum((DisMapEth_acet - Eth_acet_rep_mat(i,:)) .^ 2))
    Butyr_sd(i)= sqrt(sum((DisMapButyr - Butyr_rep_mat(i,:)) .^ 2))
     Blank_sd(i)= sqrt(sum((DisMapBlank - Blank_rep_mat(i,:)) .^ 2))
      F_pee_sd(i)= sqrt(sum((DisMapF_pee - F_pee_rep_mat(i,:)) .^ 2))
       M_pee_sd(i)= sqrt(sum((DisMapM_pee - M_pee_rep_mat(i,:)) .^ 2))
        Pin_But_sd(i)= sqrt(sum((DisMapPin_But - Pin_But_rep_mat(i,:)) .^ 2))
         TMT_sd(i)= sqrt(sum((DisMapTMT - TMT_rep_mat(i,:)) .^ 2))
          Prop_sd(i)= sqrt(sum((DisMapProp - Prop_rep_mat(i,:)) .^ 2))
           Pups_bed_sd(i)= sqrt(sum((DisMapPups_bed - Pups_bed_rep_mat(i,:)) .^ 2))
            Eth_tig_sd(i)= sqrt(sum((DisMapEth_tig - Eth_tig_rep_mat(i,:)) .^ 2))
 
 
end

all_sds{run_no}=[mean(Valer_sd),mean(Meth_prop_sd),mean(Eth_acet_sd),mean(Butyr_sd),mean(Blank_sd),mean(F_pee_sd),mean(M_pee_sd),mean(Pin_But_sd),mean(TMT_sd),mean(Prop_sd),mean(Pups_bed_sd),mean(Eth_tig_sd)]
% all_sds{run_no}=[mean(Valer_sd),mean(Meth_prop_sd),mean(Eth_acet_sd),mean(Butyr_sd),mean(Blank_sd),mean(F_pee_sd),mean(M_pee_sd),mean(Pin_But_sd),mean(TMT_sd),mean(Prop_sd),mean(Pups_bed_sd),mean(Eth_tig_sd)]


% %1
% DistMat_norm{run_no}(1,1)=sqrt(sum((DisMapValer - DisMapValer) .^ 2))/((mean(Valer_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(1,2)=sqrt(sum((DisMapValer - DisMapMeth_prop) .^ 2))/((mean(Valer_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(1,3)=sqrt(sum((DisMapValer - DisMapEth_acet) .^ 2)/((mean(Valer_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(1,4)=sqrt(sum((DisMapValer - DisMapButyr) .^ 2))/((mean(Valer_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(1,5)=sqrt(sum((DisMapValer - DisMapBlank) .^ 2))/((mean(Valer_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(1,6)=sqrt(sum((DisMapValer - DisMapTMT) .^ 2))/((mean(Valer_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(1,7)=sqrt(sum((DisMapValer - DisMapF_pee) .^ 2))/((mean(Valer_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(1,8)=sqrt(sum((DisMapValer - DisMapM_pee) .^ 2))/((mean(Valer_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(1,9)=sqrt(sum((DisMapValer - DisMapPin_But) .^ 2))/((mean(Valer_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(1,10)=sqrt(sum((DisMapValer - DisMapEth_tig) .^ 2))/((mean(Valer_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(1,11)=sqrt(sum((DisMapValer - DisMapProp) .^ 2))/((mean(Valer_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(1,12)=sqrt(sum((DisMapValer - DisMapPups_bed) .^ 2))/((mean(Valer_sd)+mean(Pups_bed_sd))/2);
% 
% % 2
% DistMat_norm{run_no}(2,1)=sqrt(sum((DisMapMeth_prop - DisMapValer) .^ 2))/((mean(Meth_prop_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(2,2)=sqrt(sum((DisMapMeth_prop - DisMapMeth_prop) .^ 2))/((mean(Meth_prop_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(2,3)=sqrt(sum((DisMapMeth_prop - DisMapEth_acet) .^ 2)/((mean(Meth_prop_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(2,4)=sqrt(sum((DisMapMeth_prop - DisMapButyr) .^ 2))/((mean(Meth_prop_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(2,5)=sqrt(sum((DisMapMeth_prop - DisMapBlank) .^ 2))/((mean(Meth_prop_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(2,6)=sqrt(sum((DisMapMeth_prop - DisMapTMT) .^ 2))/((mean(Meth_prop_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(2,7)=sqrt(sum((DisMapMeth_prop - DisMapF_pee) .^ 2))/((mean(Meth_prop_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(2,8)=sqrt(sum((DisMapMeth_prop - DisMapM_pee) .^ 2))/((mean(Meth_prop_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(2,9)=sqrt(sum((DisMapMeth_prop - DisMapPin_But) .^ 2))/((mean(Meth_prop_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(2,10)=sqrt(sum((DisMapMeth_prop - DisMapEth_tig) .^ 2))/((mean(Meth_prop_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(2,11)=sqrt(sum((DisMapMeth_prop - DisMapProp) .^ 2))/((mean(Meth_prop_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(2,12)=sqrt(sum((DisMapMeth_prop - DisMapPups_bed) .^ 2))/((mean(Meth_prop_sd)+mean(Pups_bed_sd))/2);
% 
% 
% 
% % 3
% DistMat_norm{run_no}(3,1)=sqrt(sum((DisMapEth_acet - DisMapValer) .^ 2))/((mean(Eth_acet_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(3,2)=sqrt(sum((DisMapEth_acet - DisMapMeth_prop) .^ 2))/((mean(Eth_acet_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(3,3)=sqrt(sum((DisMapEth_acet - DisMapEth_acet) .^ 2)/((mean(Eth_acet_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(3,4)=sqrt(sum((DisMapEth_acet - DisMapButyr) .^ 2))/((mean(Eth_acet_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(3,5)=sqrt(sum((DisMapEth_acet - DisMapBlank) .^ 2))/((mean(Eth_acet_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(3,6)=sqrt(sum((DisMapEth_acet - DisMapTMT) .^ 2))/((mean(Eth_acet_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(3,7)=sqrt(sum((DisMapEth_acet - DisMapF_pee) .^ 2))/((mean(Eth_acet_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(3,8)=sqrt(sum((DisMapEth_acet - DisMapM_pee) .^ 2))/((mean(Eth_acet_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(3,9)=sqrt(sum((DisMapEth_acet - DisMapPin_But) .^ 2))/((mean(Eth_acet_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(3,10)=sqrt(sum((DisMapEth_acet - DisMapEth_tig) .^ 2))/((mean(Eth_acet_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(3,11)=sqrt(sum((DisMapEth_acet - DisMapProp) .^ 2))/((mean(Eth_acet_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(3,12)=sqrt(sum((DisMapEth_acet - DisMapPups_bed) .^ 2))/((mean(Eth_acet_sd)+mean(Pups_bed_sd))/2);
% 
% 
% 
% 
% %4
% DistMat_norm{run_no}(4,1)=sqrt(sum((DisMapButyr - DisMapValer) .^ 2))/((mean(Butyr_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(4,2)=sqrt(sum((DisMapButyr - DisMapMeth_prop) .^ 2))/((mean(Butyr_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(4,3)=sqrt(sum((DisMapButyr - DisMapEth_acet) .^ 2)/((mean(Butyr_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(4,4)=sqrt(sum((DisMapButyr - DisMapButyr) .^ 2))/((mean(Butyr_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(4,5)=sqrt(sum((DisMapButyr - DisMapBlank) .^ 2))/((mean(Butyr_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(4,6)=sqrt(sum((DisMapButyr - DisMapTMT) .^ 2))/((mean(Butyr_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(4,7)=sqrt(sum((DisMapButyr - DisMapF_pee) .^ 2))/((mean(Butyr_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(4,8)=sqrt(sum((DisMapButyr - DisMapM_pee) .^ 2))/((mean(Butyr_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(4,9)=sqrt(sum((DisMapButyr - DisMapPin_But) .^ 2))/((mean(Butyr_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(4,10)=sqrt(sum((DisMapButyr - DisMapEth_tig) .^ 2))/((mean(Butyr_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(4,11)=sqrt(sum((DisMapButyr - DisMapProp) .^ 2))/((mean(Butyr_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(4,12)=sqrt(sum((DisMapButyr - DisMapPups_bed) .^ 2))/((mean(Butyr_sd)+mean(Pups_bed_sd))/2);
% 
% 
% 
% %5
% DistMat_norm{run_no}(5,1)=sqrt(sum((DisMapBlank - DisMapValer) .^ 2))/((mean(Blank_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(5,2)=sqrt(sum((DisMapBlank - DisMapMeth_prop) .^ 2))/((mean(Blank_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(5,3)=sqrt(sum((DisMapBlank - DisMapEth_acet) .^ 2)/((mean(Blank_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(5,4)=sqrt(sum((DisMapBlank - DisMapButyr) .^ 2))/((mean(Blank_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(5,5)=sqrt(sum((DisMapBlank - DisMapBlank) .^ 2))/((mean(Blank_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(5,6)=sqrt(sum((DisMapBlank - DisMapTMT) .^ 2))/((mean(Blank_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(5,7)=sqrt(sum((DisMapBlank - DisMapF_pee) .^ 2))/((mean(Blank_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(5,8)=sqrt(sum((DisMapBlank - DisMapM_pee) .^ 2))/((mean(Blank_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(5,9)=sqrt(sum((DisMapBlank - DisMapPin_But) .^ 2))/((mean(Blank_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(5,10)=sqrt(sum((DisMapBlank - DisMapEth_tig) .^ 2))/((mean(Blank_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(5,11)=sqrt(sum((DisMapBlank - DisMapProp) .^ 2))/((mean(Blank_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(5,12)=sqrt(sum((DisMapBlank - DisMapPups_bed) .^ 2))/((mean(Blank_sd)+mean(Pups_bed_sd))/2);
% 
% 
% 
% %6
% DistMat_norm{run_no}(6,1)=sqrt(sum((DisMapTMT - DisMapValer) .^ 2))/((mean(TMT_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(6,2)=sqrt(sum((DisMapTMT - DisMapMeth_prop) .^ 2))/((mean(TMT_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(6,3)=sqrt(sum((DisMapTMT - DisMapEth_acet) .^ 2)/((mean(TMT_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(6,4)=sqrt(sum((DisMapTMT - DisMapButyr) .^ 2))/((mean(TMT_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(6,5)=sqrt(sum((DisMapTMT - DisMapBlank) .^ 2))/((mean(TMT_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(6,6)=sqrt(sum((DisMapTMT - DisMapTMT) .^ 2))/((mean(TMT_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(6,7)=sqrt(sum((DisMapTMT - DisMapF_pee) .^ 2))/((mean(TMT_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(6,8)=sqrt(sum((DisMapTMT - DisMapM_pee) .^ 2))/((mean(TMT_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(6,9)=sqrt(sum((DisMapTMT - DisMapPin_But) .^ 2))/((mean(TMT_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(6,10)=sqrt(sum((DisMapTMT - DisMapEth_tig) .^ 2))/((mean(TMT_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(6,11)=sqrt(sum((DisMapTMT - DisMapProp) .^ 2))/((mean(TMT_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(6,12)=sqrt(sum((DisMapTMT - DisMapPups_bed) .^ 2))/((mean(TMT_sd)+mean(Pups_bed_sd))/2);
% 
% 
% 
% %7
% DistMat_norm{run_no}(7,1)=sqrt(sum((DisMapF_pee - DisMapValer) .^ 2))/((mean(F_pee_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(7,2)=sqrt(sum((DisMapF_pee - DisMapMeth_prop) .^ 2))/((mean(F_pee_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(7,3)=sqrt(sum((DisMapF_pee - DisMapEth_acet) .^ 2)/((mean(F_pee_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(7,4)=sqrt(sum((DisMapF_pee - DisMapButyr) .^ 2))/((mean(F_pee_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(7,5)=sqrt(sum((DisMapF_pee - DisMapBlank) .^ 2))/((mean(F_pee_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(7,6)=sqrt(sum((DisMapF_pee - DisMapTMT) .^ 2))/((mean(F_pee_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(7,7)=sqrt(sum((DisMapF_pee - DisMapF_pee) .^ 2))/((mean(F_pee_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(7,8)=sqrt(sum((DisMapF_pee - DisMapM_pee) .^ 2))/((mean(F_pee_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(7,9)=sqrt(sum((DisMapF_pee - DisMapPin_But) .^ 2))/((mean(F_pee_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(7,10)=sqrt(sum((DisMapF_pee - DisMapEth_tig) .^ 2))/((mean(F_pee_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(7,11)=sqrt(sum((DisMapF_pee - DisMapProp) .^ 2))/((mean(F_pee_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(7,12)=sqrt(sum((DisMapF_pee - DisMapPups_bed) .^ 2))/((mean(F_pee_sd)+mean(Pups_bed_sd))/2);
% 
% 
% %8
% DistMat_norm{run_no}(8,1)=sqrt(sum((DisMapM_pee - DisMapValer) .^ 2))/((mean(M_pee_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(8,2)=sqrt(sum((DisMapM_pee - DisMapMeth_prop) .^ 2))/((mean(M_pee_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(8,3)=sqrt(sum((DisMapM_pee - DisMapEth_acet) .^ 2)/((mean(M_pee_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(8,4)=sqrt(sum((DisMapM_pee - DisMapButyr) .^ 2))/((mean(M_pee_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(8,5)=sqrt(sum((DisMapM_pee - DisMapBlank) .^ 2))/((mean(M_pee_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(8,6)=sqrt(sum((DisMapM_pee - DisMapTMT) .^ 2))/((mean(M_pee_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(8,7)=sqrt(sum((DisMapM_pee - DisMapF_pee) .^ 2))/((mean(M_pee_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(8,8)=sqrt(sum((DisMapM_pee - DisMapM_pee) .^ 2))/((mean(M_pee_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(8,9)=sqrt(sum((DisMapM_pee - DisMapPin_But) .^ 2))/((mean(M_pee_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(8,10)=sqrt(sum((DisMapM_pee - DisMapEth_tig) .^ 2))/((mean(M_pee_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(8,11)=sqrt(sum((DisMapM_pee - DisMapProp) .^ 2))/((mean(M_pee_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(8,12)=sqrt(sum((DisMapM_pee - DisMapPups_bed) .^ 2))/((mean(M_pee_sd)+mean(Pups_bed_sd))/2);
% 
% %9
% DistMat_norm{run_no}(9,1)=sqrt(sum((DisMapPin_But - DisMapValer) .^ 2))/((mean(Pin_But_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(9,2)=sqrt(sum((DisMapPin_But - DisMapMeth_prop) .^ 2))/((mean(Pin_But_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(9,3)=sqrt(sum((DisMapPin_But - DisMapEth_acet) .^ 2)/((mean(Pin_But_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(9,4)=sqrt(sum((DisMapPin_But - DisMapButyr) .^ 2))/((mean(Pin_But_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(9,5)=sqrt(sum((DisMapPin_But - DisMapBlank) .^ 2))/((mean(Pin_But_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(9,6)=sqrt(sum((DisMapPin_But - DisMapTMT) .^ 2))/((mean(Pin_But_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(9,7)=sqrt(sum((DisMapPin_But - DisMapF_pee) .^ 2))/((mean(Pin_But_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(9,8)=sqrt(sum((DisMapPin_But - DisMapM_pee) .^ 2))/((mean(Pin_But_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(9,9)=sqrt(sum((DisMapPin_But - DisMapPin_But) .^ 2))/((mean(Pin_But_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(9,10)=sqrt(sum((DisMapPin_But - DisMapEth_tig) .^ 2))/((mean(Pin_But_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(9,11)=sqrt(sum((DisMapPin_But - DisMapProp) .^ 2))/((mean(Pin_But_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(9,12)=sqrt(sum((DisMapPin_But - DisMapPups_bed) .^ 2))/((mean(Pin_But_sd)+mean(Pups_bed_sd))/2);
% 
% 
% %10
% DistMat_norm{run_no}(10,1)=sqrt(sum((DisMapEth_tig - DisMapValer) .^ 2))/((mean(Eth_tig_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(10,2)=sqrt(sum((DisMapEth_tig - DisMapMeth_prop) .^ 2))/((mean(Eth_tig_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(10,3)=sqrt(sum((DisMapEth_tig - DisMapEth_acet) .^ 2)/((mean(Eth_tig_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(10,4)=sqrt(sum((DisMapEth_tig - DisMapButyr) .^ 2))/((mean(Eth_tig_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(10,5)=sqrt(sum((DisMapEth_tig - DisMapBlank) .^ 2))/((mean(Eth_tig_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(10,6)=sqrt(sum((DisMapEth_tig - DisMapTMT) .^ 2))/((mean(Eth_tig_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(10,7)=sqrt(sum((DisMapEth_tig - DisMapF_pee) .^ 2))/((mean(Eth_tig_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(10,8)=sqrt(sum((DisMapEth_tig - DisMapM_pee) .^ 2))/((mean(Eth_tig_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(10,9)=sqrt(sum((DisMapEth_tig - DisMapPin_But) .^ 2))/((mean(Eth_tig_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(10,10)=sqrt(sum((DisMapEth_tig - DisMapEth_tig) .^ 2))/((mean(Eth_tig_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(10,11)=sqrt(sum((DisMapEth_tig - DisMapProp) .^ 2))/((mean(Eth_tig_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(10,12)=sqrt(sum((DisMapEth_tig - DisMapPups_bed) .^ 2))/((mean(Eth_tig_sd)+mean(Pups_bed_sd))/2);
% 
% 
% %11
% DistMat_norm{run_no}(11,1)=sqrt(sum((DisMapProp - DisMapValer) .^ 2))/((mean(Prop_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(11,2)=sqrt(sum((DisMapProp - DisMapMeth_prop) .^ 2))/((mean(Prop_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(11,3)=sqrt(sum((DisMapProp - DisMapEth_acet) .^ 2)/((mean(Prop_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(11,4)=sqrt(sum((DisMapProp - DisMapButyr) .^ 2))/((mean(Prop_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(11,5)=sqrt(sum((DisMapProp - DisMapBlank) .^ 2))/((mean(Prop_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(11,6)=sqrt(sum((DisMapProp - DisMapTMT) .^ 2))/((mean(Prop_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(11,7)=sqrt(sum((DisMapProp - DisMapF_pee) .^ 2))/((mean(Prop_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(11,8)=sqrt(sum((DisMapProp - DisMapM_pee) .^ 2))/((mean(Prop_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(11,9)=sqrt(sum((DisMapProp - DisMapPin_But) .^ 2))/((mean(Prop_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(11,10)=sqrt(sum((DisMapProp - DisMapEth_tig) .^ 2))/((mean(Prop_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(11,11)=sqrt(sum((DisMapProp - DisMapProp) .^ 2))/((mean(Prop_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(11,12)=sqrt(sum((DisMapProp - DisMapPups_bed) .^ 2))/((mean(Prop_sd)+mean(Pups_bed_sd))/2);
% 
% 
% 
% %12
% DistMat_norm{run_no}(12,1)=sqrt(sum((DisMapPups_bed - DisMapValer) .^ 2))/((mean(Pups_bed_sd)+mean(Valer_sd))/2);
% DistMat_norm{run_no}(12,2)=sqrt(sum((DisMapPups_bed - DisMapMeth_prop) .^ 2))/((mean(Pups_bed_sd)+mean(Meth_prop_sd))/2);
% DistMat_norm{run_no}(12,3)=sqrt(sum((DisMapPups_bed - DisMapEth_acet) .^ 2)/((mean(Pups_bed_sd)+mean(Eth_acet_sd))/2));
% DistMat_norm{run_no}(12,4)=sqrt(sum((DisMapPups_bed - DisMapButyr) .^ 2))/((mean(Pups_bed_sd)+mean(Butyr_sd))/2);
% DistMat_norm{run_no}(12,5)=sqrt(sum((DisMapPups_bed - DisMapBlank) .^ 2))/((mean(Pups_bed_sd)+mean(Blank_sd))/2);
% DistMat_norm{run_no}(12,6)=sqrt(sum((DisMapPups_bed - DisMapTMT) .^ 2))/((mean(Pups_bed_sd)+mean(TMT_sd))/2);
% DistMat_norm{run_no}(12,7)=sqrt(sum((DisMapPups_bed - DisMapF_pee) .^ 2))/((mean(Pups_bed_sd)+mean(F_pee_sd))/2);
% DistMat_norm{run_no}(12,8)=sqrt(sum((DisMapPups_bed - DisMapM_pee) .^ 2))/((mean(Pups_bed_sd)+mean(M_pee_sd))/2);
% DistMat_norm{run_no}(12,9)=sqrt(sum((DisMapPups_bed - DisMapPin_But) .^ 2))/((mean(Pups_bed_sd)+mean(Pin_But_sd))/2);
% DistMat_norm{run_no}(12,10)=sqrt(sum((DisMapPups_bed - DisMapEth_tig) .^ 2))/((mean(Pups_bed_sd)+mean(Eth_tig_sd))/2);
% DistMat_norm{run_no}(12,11)=sqrt(sum((DisMapPups_bed - DisMapProp) .^ 2))/((mean(Pups_bed_sd)+mean(Prop_sd))/2);
% DistMat_norm{run_no}(12,12)=sqrt(sum((DisMapPups_bed - DisMapPups_bed) .^ 2))/((mean(Pups_bed_sd)+mean(Pups_bed_sd))/2);



%% normalized by inner variance and organized:
%1
DistMat_norm{run_no}(1,1)=sqrt(sum((DisMapBlank - DisMapBlank) .^ 2))/((mean(Blank_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(1,2)=sqrt(sum((DisMapBlank - DisMapValer) .^ 2))/((mean(Blank_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(1,3)=sqrt(sum((DisMapBlank - DisMapMeth_prop) .^ 2))/((mean(Blank_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(1,4)=sqrt(sum((DisMapBlank - DisMapEth_acet) .^ 2))/((mean(Blank_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(1,5)=sqrt(sum((DisMapBlank - DisMapButyr) .^ 2))/((mean(Blank_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(1,6)=sqrt(sum((DisMapBlank - DisMapEth_tig) .^ 2))/((mean(Blank_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(1,7)=sqrt(sum((DisMapBlank - DisMapProp) .^ 2))/((mean(Blank_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(1,8)=sqrt(sum((DisMapBlank - DisMapTMT) .^ 2))/((mean(Blank_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(1,9)=sqrt(sum((DisMapBlank - DisMapF_pee) .^ 2))/((mean(Blank_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(1,10)=sqrt(sum((DisMapBlank - DisMapM_pee) .^ 2))/((mean(Blank_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(1,11)=sqrt(sum((DisMapBlank - DisMapPin_But) .^ 2))/((mean(Blank_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(1,12)=sqrt(sum((DisMapBlank - DisMapPups_bed) .^ 2))/((mean(Blank_sd)+mean(Pups_bed_sd))/2);

% 2
DistMat_norm{run_no}(2,1)=sqrt(sum((DisMapValer - DisMapBlank) .^ 2))/((mean(Valer_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(2,2)=sqrt(sum((DisMapValer - DisMapValer) .^ 2))/((mean(Valer_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(2,3)=sqrt(sum((DisMapValer - DisMapMeth_prop) .^ 2))/((mean(Valer_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(2,4)=sqrt(sum((DisMapValer - DisMapEth_acet) .^ 2))/((mean(Valer_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(2,5)=sqrt(sum((DisMapValer - DisMapButyr) .^ 2))/((mean(Valer_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(2,6)=sqrt(sum((DisMapValer - DisMapEth_tig) .^ 2))/((mean(Valer_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(2,7)=sqrt(sum((DisMapValer - DisMapProp) .^ 2))/((mean(Valer_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(2,8)=sqrt(sum((DisMapValer - DisMapTMT) .^ 2))/((mean(Valer_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(2,9)=sqrt(sum((DisMapValer - DisMapF_pee) .^ 2))/((mean(Valer_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(2,10)=sqrt(sum((DisMapValer - DisMapM_pee) .^ 2))/((mean(Valer_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(2,11)=sqrt(sum((DisMapValer - DisMapPin_But) .^ 2))/((mean(Valer_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(2,12)=sqrt(sum((DisMapValer - DisMapPups_bed) .^ 2))/((mean(Valer_sd)+mean(Pups_bed_sd))/2);

% 3
DistMat_norm{run_no}(3,1)=sqrt(sum((DisMapMeth_prop - DisMapBlank) .^ 2))/((mean(Meth_prop_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(3,2)=sqrt(sum((DisMapMeth_prop - DisMapValer) .^ 2))/((mean(Meth_prop_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(3,3)=sqrt(sum((DisMapMeth_prop - DisMapMeth_prop) .^ 2))/((mean(Meth_prop_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(3,4)=sqrt(sum((DisMapMeth_prop - DisMapEth_acet) .^ 2))/((mean(Meth_prop_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(3,5)=sqrt(sum((DisMapMeth_prop - DisMapButyr) .^ 2))/((mean(Meth_prop_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(3,6)=sqrt(sum((DisMapMeth_prop - DisMapEth_tig) .^ 2))/((mean(Meth_prop_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(3,7)=sqrt(sum((DisMapMeth_prop - DisMapProp) .^ 2))/((mean(Meth_prop_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(3,8)=sqrt(sum((DisMapMeth_prop - DisMapTMT) .^ 2))/((mean(Meth_prop_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(3,9)=sqrt(sum((DisMapMeth_prop - DisMapF_pee) .^ 2))/((mean(Meth_prop_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(3,10)=sqrt(sum((DisMapMeth_prop - DisMapM_pee) .^ 2))/((mean(Meth_prop_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(3,11)=sqrt(sum((DisMapMeth_prop - DisMapPin_But) .^ 2))/((mean(Meth_prop_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(3,12)=sqrt(sum((DisMapMeth_prop - DisMapPups_bed) .^ 2))/((mean(Meth_prop_sd)+mean(Pups_bed_sd))/2);

%4
DistMat_norm{run_no}(4,1)=sqrt(sum((DisMapEth_acet - DisMapBlank) .^ 2))/((mean(Eth_acet_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(4,2)=sqrt(sum((DisMapEth_acet - DisMapValer) .^ 2))/((mean(Eth_acet_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(4,3)=sqrt(sum((DisMapEth_acet - DisMapMeth_prop) .^ 2))/((mean(Eth_acet_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(4,4)=sqrt(sum((DisMapEth_acet - DisMapEth_acet) .^ 2))/((mean(Eth_acet_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(4,5)=sqrt(sum((DisMapEth_acet - DisMapButyr) .^ 2))/((mean(Eth_acet_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(4,6)=sqrt(sum((DisMapEth_acet - DisMapEth_tig) .^ 2))/((mean(Eth_acet_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(4,7)=sqrt(sum((DisMapEth_acet - DisMapProp) .^ 2))/((mean(Eth_acet_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(4,8)=sqrt(sum((DisMapEth_acet - DisMapTMT) .^ 2))/((mean(Eth_acet_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(4,9)=sqrt(sum((DisMapEth_acet - DisMapF_pee) .^ 2))/((mean(Eth_acet_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(4,10)=sqrt(sum((DisMapEth_acet - DisMapM_pee) .^ 2))/((mean(Eth_acet_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(4,11)=sqrt(sum((DisMapEth_acet - DisMapPin_But) .^ 2))/((mean(Eth_acet_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(4,12)=sqrt(sum((DisMapEth_acet - DisMapPups_bed) .^ 2))/((mean(Eth_acet_sd)+mean(Pups_bed_sd))/2);

%5
DistMat_norm{run_no}(5,1)=sqrt(sum((DisMapButyr - DisMapBlank) .^ 2))/((mean(Butyr_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(5,2)=sqrt(sum((DisMapButyr - DisMapValer) .^ 2))/((mean(Butyr_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(5,3)=sqrt(sum((DisMapButyr - DisMapMeth_prop) .^ 2))/((mean(Butyr_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(5,4)=sqrt(sum((DisMapButyr - DisMapEth_acet) .^ 2))/((mean(Butyr_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(5,5)=sqrt(sum((DisMapButyr - DisMapButyr) .^ 2))/((mean(Butyr_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(5,6)=sqrt(sum((DisMapButyr - DisMapEth_tig) .^ 2))/((mean(Butyr_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(5,7)=sqrt(sum((DisMapButyr - DisMapProp) .^ 2))/((mean(Butyr_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(5,8)=sqrt(sum((DisMapButyr - DisMapTMT) .^ 2))/((mean(Butyr_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(5,9)=sqrt(sum((DisMapButyr - DisMapF_pee) .^ 2))/((mean(Butyr_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(5,10)=sqrt(sum((DisMapButyr - DisMapM_pee) .^ 2))/((mean(Butyr_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(5,11)=sqrt(sum((DisMapButyr - DisMapPin_But) .^ 2))/((mean(Butyr_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(5,12)=sqrt(sum((DisMapButyr - DisMapPups_bed) .^ 2))/((mean(Butyr_sd)+mean(Pups_bed_sd))/2);


%6
DistMat_norm{run_no}(6,1)=sqrt(sum((DisMapEth_tig - DisMapBlank) .^ 2))/((mean(Eth_tig_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(6,2)=sqrt(sum((DisMapEth_tig - DisMapValer) .^ 2))/((mean(Eth_tig_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(6,3)=sqrt(sum((DisMapEth_tig - DisMapMeth_prop) .^ 2))/((mean(Eth_tig_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(6,4)=sqrt(sum((DisMapEth_tig - DisMapEth_acet) .^ 2))/((mean(Eth_tig_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(6,5)=sqrt(sum((DisMapEth_tig - DisMapButyr) .^ 2))/((mean(Eth_tig_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(6,6)=sqrt(sum((DisMapEth_tig - DisMapEth_tig) .^ 2))/((mean(Eth_tig_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(6,7)=sqrt(sum((DisMapEth_tig - DisMapProp) .^ 2))/((mean(Eth_tig_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(6,8)=sqrt(sum((DisMapEth_tig - DisMapTMT) .^ 2))/((mean(Eth_tig_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(6,9)=sqrt(sum((DisMapEth_tig - DisMapF_pee) .^ 2))/((mean(Eth_tig_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(6,10)=sqrt(sum((DisMapEth_tig - DisMapM_pee) .^ 2))/((mean(Eth_tig_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(6,11)=sqrt(sum((DisMapEth_tig - DisMapPin_But) .^ 2))/((mean(Eth_tig_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(6,12)=sqrt(sum((DisMapEth_tig - DisMapPups_bed) .^ 2))/((mean(Eth_tig_sd)+mean(Pups_bed_sd))/2);

%7

DistMat_norm{run_no}(7,1)=sqrt(sum((DisMapProp - DisMapBlank) .^ 2))/((mean(Prop_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(7,2)=sqrt(sum((DisMapProp - DisMapValer) .^ 2))/((mean(Prop_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(7,3)=sqrt(sum((DisMapProp - DisMapMeth_prop) .^ 2))/((mean(Prop_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(7,4)=sqrt(sum((DisMapProp - DisMapEth_acet) .^ 2))/((mean(Prop_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(7,5)=sqrt(sum((DisMapProp - DisMapButyr) .^ 2))/((mean(Prop_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(7,6)=sqrt(sum((DisMapProp - DisMapEth_tig) .^ 2))/((mean(Prop_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(7,7)=sqrt(sum((DisMapProp - DisMapProp) .^ 2))/((mean(Prop_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(7,8)=sqrt(sum((DisMapProp - DisMapTMT) .^ 2))/((mean(Prop_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(7,9)=sqrt(sum((DisMapProp - DisMapF_pee) .^ 2))/((mean(Prop_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(7,10)=sqrt(sum((DisMapProp - DisMapM_pee) .^ 2))/((mean(Prop_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(7,11)=sqrt(sum((DisMapProp - DisMapPin_But) .^ 2))/((mean(Prop_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(7,12)=sqrt(sum((DisMapProp - DisMapPups_bed) .^ 2))/((mean(Prop_sd)+mean(Pups_bed_sd))/2);

%8

DistMat_norm{run_no}(8,1)=sqrt(sum((DisMapTMT - DisMapBlank) .^ 2))/((mean(TMT_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(8,2)=sqrt(sum((DisMapTMT - DisMapValer) .^ 2))/((mean(TMT_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(8,3)=sqrt(sum((DisMapTMT - DisMapMeth_prop) .^ 2))/((mean(TMT_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(8,4)=sqrt(sum((DisMapTMT - DisMapEth_acet) .^ 2))/((mean(TMT_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(8,5)=sqrt(sum((DisMapTMT - DisMapButyr) .^ 2))/((mean(TMT_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(8,6)=sqrt(sum((DisMapTMT - DisMapEth_tig) .^ 2))/((mean(TMT_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(8,7)=sqrt(sum((DisMapTMT - DisMapProp) .^ 2))/((mean(TMT_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(8,8)=sqrt(sum((DisMapTMT - DisMapTMT) .^ 2))/((mean(TMT_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(8,9)=sqrt(sum((DisMapTMT - DisMapF_pee) .^ 2))/((mean(TMT_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(8,10)=sqrt(sum((DisMapTMT - DisMapM_pee) .^ 2))/((mean(TMT_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(8,11)=sqrt(sum((DisMapTMT - DisMapPin_But) .^ 2))/((mean(TMT_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(8,12)=sqrt(sum((DisMapTMT - DisMapPups_bed) .^ 2))/((mean(TMT_sd)+mean(Pups_bed_sd))/2);


%9

DistMat_norm{run_no}(9,1)=sqrt(sum((DisMapF_pee - DisMapBlank) .^ 2))/((mean(F_pee_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(9,2)=sqrt(sum((DisMapF_pee - DisMapValer) .^ 2))/((mean(F_pee_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(9,3)=sqrt(sum((DisMapF_pee - DisMapMeth_prop) .^ 2))/((mean(F_pee_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(9,4)=sqrt(sum((DisMapF_pee - DisMapEth_acet) .^ 2))/((mean(F_pee_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(9,5)=sqrt(sum((DisMapF_pee - DisMapButyr) .^ 2))/((mean(F_pee_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(9,6)=sqrt(sum((DisMapF_pee - DisMapEth_tig) .^ 2))/((mean(F_pee_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(9,7)=sqrt(sum((DisMapF_pee - DisMapProp) .^ 2))/((mean(F_pee_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(9,8)=sqrt(sum((DisMapF_pee - DisMapTMT) .^ 2))/((mean(F_pee_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(9,9)=sqrt(sum((DisMapF_pee - DisMapF_pee) .^ 2))/((mean(F_pee_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(9,10)=sqrt(sum((DisMapF_pee - DisMapM_pee) .^ 2))/((mean(F_pee_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(9,11)=sqrt(sum((DisMapF_pee - DisMapPin_But) .^ 2))/((mean(F_pee_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(9,12)=sqrt(sum((DisMapF_pee - DisMapPups_bed) .^ 2))/((mean(F_pee_sd)+mean(Pups_bed_sd))/2);


%10

DistMat_norm{run_no}(10,1)=sqrt(sum((DisMapM_pee - DisMapBlank) .^ 2))/((mean(M_pee_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(10,2)=sqrt(sum((DisMapM_pee - DisMapValer) .^ 2))/((mean(M_pee_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(10,3)=sqrt(sum((DisMapM_pee - DisMapMeth_prop) .^ 2))/((mean(M_pee_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(10,4)=sqrt(sum((DisMapM_pee - DisMapEth_acet) .^ 2))/((mean(M_pee_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(10,5)=sqrt(sum((DisMapM_pee - DisMapButyr) .^ 2))/((mean(M_pee_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(10,6)=sqrt(sum((DisMapM_pee - DisMapEth_tig) .^ 2))/((mean(M_pee_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(10,7)=sqrt(sum((DisMapM_pee - DisMapProp) .^ 2))/((mean(M_pee_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(10,8)=sqrt(sum((DisMapM_pee - DisMapTMT) .^ 2))/((mean(M_pee_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(10,9)=sqrt(sum((DisMapM_pee - DisMapF_pee) .^ 2))/((mean(M_pee_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(10,10)=sqrt(sum((DisMapM_pee - DisMapM_pee) .^ 2))/((mean(M_pee_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(10,11)=sqrt(sum((DisMapM_pee - DisMapPin_But) .^ 2))/((mean(M_pee_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(10,12)=sqrt(sum((DisMapM_pee - DisMapPups_bed) .^ 2))/((mean(M_pee_sd)+mean(Pups_bed_sd))/2);

%11

DistMat_norm{run_no}(11,1)=sqrt(sum((DisMapPin_But - DisMapBlank) .^ 2))/((mean(Pin_But_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(11,2)=sqrt(sum((DisMapPin_But - DisMapValer) .^ 2))/((mean(Pin_But_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(11,3)=sqrt(sum((DisMapPin_But - DisMapMeth_prop) .^ 2))/((mean(Pin_But_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(11,4)=sqrt(sum((DisMapPin_But - DisMapEth_acet) .^ 2))/((mean(Pin_But_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(11,5)=sqrt(sum((DisMapPin_But - DisMapButyr) .^ 2))/((mean(Pin_But_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(11,6)=sqrt(sum((DisMapPin_But - DisMapEth_tig) .^ 2))/((mean(Pin_But_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(11,7)=sqrt(sum((DisMapPin_But - DisMapProp) .^ 2))/((mean(Pin_But_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(11,8)=sqrt(sum((DisMapPin_But - DisMapTMT) .^ 2))/((mean(Pin_But_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(11,9)=sqrt(sum((DisMapPin_But - DisMapF_pee) .^ 2))/((mean(Pin_But_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(11,10)=sqrt(sum((DisMapPin_But - DisMapM_pee) .^ 2))/((mean(Pin_But_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(11,11)=sqrt(sum((DisMapPin_But - DisMapPin_But) .^ 2))/((mean(Pin_But_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(11,12)=sqrt(sum((DisMapPin_But - DisMapPups_bed) .^ 2))/((mean(Pin_But_sd)+mean(Pups_bed_sd))/2);


%12

DistMat_norm{run_no}(12,1)=sqrt(sum((DisMapPups_bed - DisMapBlank) .^ 2))/((mean(Pups_bed_sd)+mean(Blank_sd))/2);
DistMat_norm{run_no}(12,2)=sqrt(sum((DisMapPups_bed - DisMapValer) .^ 2))/((mean(Pups_bed_sd)+mean(Valer_sd))/2);
DistMat_norm{run_no}(12,3)=sqrt(sum((DisMapPups_bed - DisMapMeth_prop) .^ 2))/((mean(Pups_bed_sd)+mean(Meth_prop_sd))/2);
DistMat_norm{run_no}(12,4)=sqrt(sum((DisMapPups_bed - DisMapEth_acet) .^ 2))/((mean(Pups_bed_sd)+mean(Eth_acet_sd))/2);
DistMat_norm{run_no}(12,5)=sqrt(sum((DisMapPups_bed - DisMapButyr) .^ 2))/((mean(Pups_bed_sd)+mean(Butyr_sd))/2);
DistMat_norm{run_no}(12,6)=sqrt(sum((DisMapPups_bed - DisMapEth_tig) .^ 2))/((mean(Pups_bed_sd)+mean(Eth_tig_sd))/2);
DistMat_norm{run_no}(12,7)=sqrt(sum((DisMapPups_bed - DisMapProp) .^ 2))/((mean(Pups_bed_sd)+mean(Prop_sd))/2);
DistMat_norm{run_no}(12,8)=sqrt(sum((DisMapPups_bed - DisMapTMT) .^ 2))/((mean(Pups_bed_sd)+mean(TMT_sd))/2);
DistMat_norm{run_no}(12,9)=sqrt(sum((DisMapPups_bed - DisMapF_pee) .^ 2))/((mean(Pups_bed_sd)+mean(F_pee_sd))/2);
DistMat_norm{run_no}(12,10)=sqrt(sum((DisMapPups_bed - DisMapM_pee) .^ 2))/((mean(Pups_bed_sd)+mean(M_pee_sd))/2);
DistMat_norm{run_no}(12,11)=sqrt(sum((DisMapPups_bed - DisMapPin_But) .^ 2))/((mean(Pups_bed_sd)+mean(Pin_But_sd))/2);
DistMat_norm{run_no}(12,12)=sqrt(sum((DisMapPups_bed - DisMapPups_bed) .^ 2))/((mean(Pups_bed_sd)+mean(Pups_bed_sd))/2);


%% plot before AND AFTER regular- not normalized by variance

if run_no==1
figure
hold on
title('Before')
clims=[0 3];
imagesc(flip((DistMat{1})),clims)
%colorbar
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
% xlabel('PCA1')
% ylabel('PCA2')
% zlabel('PCA3')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight
end

if run_no==2
figure
hold on
title('After')
clims=[0 3];
imagesc(flip(DistMat{2}),clims)
%colorbar
mean(mean(DistMat{run_no}))
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight
end
% xlabel('PCA1')
% ylabel('PCA2')
% zlabel('PCA3')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight


%% test for all distances 
% A = reshape(DistMat{1,1},[],1);
% B = reshape(DistMat{1,2},[],1);
% [~,~,EucDis{run_no}]=find(triu(DistMat{run_no}));
% [h,p]=ttest(EucDis{1},EucDis{2})
if run_no==2

A=triu(DistMat{1})
A(A==0)=[];
nans=isnan(A);
A(nans==1)=[];

B=triu(DistMat{2})
B(B==0)=[];
nans=isnan(B);
B(nans==1)=[];

%tests for all
[p,h]=signrank(A,B);
[h,p]=ttest(A,B);


%test for blank only
[h,p]=ttest(DistMat{1}(:,5),DistMat{2}(:,5));


%% create all sort of comparisons

Deltas_Dist_mat=DistMat{1}-DistMat{2};
Deltas_Dist_vec=A-B;

Deltas_Dist_mat_norm=Deltas_Dist_mat./((DistMat{1}+DistMat{2}));

nans=isnan(Deltas_Dist_mat_norm);
Deltas_Dist_mat_norm(nans==1)=0;


Deltas_Dist_vec_norm=(A-B)./((A+B));


Dist_mat_ratios=DistMat{1}./DistMat{2};
Dist_vec_ratios=A./B;

Dist_vec_ratios_log=log2(Dist_vec_ratios);

%% plot deltas
figure
hold on
title('just Deltas unnormalized distances')
clims=[-1 1];
imagesc(flip(Deltas_Dist_mat),clims)

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size

set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight
%% plot normalized deltas
figure
hold on
title('normalized Deltas unnormalized distances')
clims=[-1 1];
imagesc(flip(Deltas_Dist_mat_norm),clims)
%colorbar
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size

set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight
% CAL_MIR_PCA.m
% Displaying CAL_MIR_PCA.m.

% %% plot ratios
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

% %% plot log ratios
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


dis_vec=triu(Deltas_Dist_mat_norm);%get rid of zeros and nans
dis_vec=reshape(dis_vec,1,144);
dis_vec(dis_vec==0)=[];
nans=isnan(dis_vec);
dis_vec(nans==1)=[];

x=1:length(dis_vec);
y=zeros(size(x));
% figure
% hold on
% plot(x,dis_vec,'k','linewidth',2)
% plot(x,y,'r--','linewidth',2)
% xlim([1 length(dis_vec)])

% %%
dis_vec_ratios=triu(Dist_mat_ratios);%get rid of zeros and nans
dis_vec_ratios=reshape(dis_vec_ratios,1,144);
dis_vec_ratios(dis_vec_ratios==0)=[];
nans=isnan(dis_vec_ratios);
dis_vec_ratios(nans==1)=[];

x=1:length(dis_vec_ratios);
y=ones(size(x));
% figure
% hold on
% plot(x,dis_vec_ratios,'k','linewidth',2)
% plot(x,y,'r--','linewidth',2)
% xlim([1 length(dis_vec_ratios)])
% 
% sum(dis_vec_ratios>=1)
% length(dis_vec_ratios)
end


%% plot results for norm dist
if run_no==2
exp=Deltas_Dist_vec_norm
cont=exp%Deltas_Dist_vec_norm_cont; %chnge this to exp if no control yet

x1 = linspace(4,5,size(cont,2));
x2 = linspace(6,7,size(exp,2));
mean_vec_cont=mean(cont)*ones(size(x1))
mean_vec_exp=mean(exp)*ones(size(x2))
x3=3:8;
y=zeros(size(x3))

figure
title('Norm deltas regular distances')
xlim([3 8])
hold on
plot(x1,cont,'k*')
plot(x2,exp,'ko')
plot(x1,mean_vec_cont,'r--','linewidth',2)
plot(x2,mean_vec_exp,'r--','linewidth',2)
plot(x3,y,'k--')
XTick=[]


sem_cont=std(cont)/sqrt(length(cont))
sem_exp=std(exp)/sqrt(length(exp))

x_er=[4.5 6.5];
y_er=[mean(cont) mean(exp)];
sems=[sem_cont sem_exp];
errorbar(x_er,y_er,sems,'r.','linewidth',2)

box off

[h,p]=ttest2(exp,cont)
[h,p]=ranksum(exp,cont)

%% %% plot before AND AFTER____normalized by var
figure
hold on
title('Before norm_by_var')
clims=[0 3];
imagesc(flip(DistMat_norm{1}),clims)
colorbar
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
% xlabel('PCA1')
% ylabel('PCA2')
% zlabel('PCA3')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight

figure
hold on
title('After norm_by_var')
clims=[0 3];
imagesc(flip(DistMat_norm{2}),clims)
colorbar
mean(mean(DistMat_norm{run_no}))
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight

% xlabel('PCA1')
% ylabel('PCA2')
% zlabel('PCA3')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight
set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight


%% test for all distances 
% A = reshape(DistMat_norm{1,1},[],1);
% B = reshape(DistMat_norm{1,2},[],1);
% [~,~,EucDis{run_no}]=find(triu(DistMat_norm{run_no}));
% [h,p]=ttest(EucDis{1},EucDis{2})


A=triu(DistMat_norm{1})
A(A==0)=[];
nans=isnan(A);
A(nans==1)=[];

B=triu(DistMat_norm{2})
B(B==0)=[];
nans=isnan(B);
B(nans==1)=[];

%tests for all
[p,h]=signrank(A,B);
[h,p]=ttest(A,B);


%test for blank only
[h,p]=ttest(DistMat_norm{1}(:,1),DistMat_norm{2}(:,1));


%% create all sort of comparisons

Deltas_Dist_mat=DistMat_norm{1}-DistMat_norm{2};
Deltas_Dist_vec=A-B;

Deltas_Dist_mat_norm=Deltas_Dist_mat./((DistMat_norm{1}+DistMat_norm{2}));

nans=isnan(Deltas_Dist_mat_norm);
Deltas_Dist_mat_norm(nans==1)=0;


Deltas_Dist_vec_norm=(A-B)./((A+B));


Dist_mat_ratios=DistMat_norm{1}./DistMat_norm{2};
Dist_vec_ratios=A./B;

Dist_vec_ratios_log=log2(Dist_vec_ratios);

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
%% plot normalized deltas
figure
hold on
title('normalized Deltas')
clims=[-1 1];
imagesc(flip(Deltas_Dist_mat_norm),clims)

%colorbar                             %cancel this for figure square plot
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'color','w'); %<- Set size

set(gca, 'FontSize', fsz, 'LineWidth', alw,'FontName','times new roman'); %<- Set properties
axis tight
% CAL_MIR_PCA.m
% Displaying CAL_MIR_PCA.m.

% %% plot ratios
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


dis_vec=triu(Deltas_Dist_mat_norm);%get rid of zeros and nans
dis_vec=reshape(dis_vec,1,144);
dis_vec(dis_vec==0)=[];
nans=isnan(dis_vec);
dis_vec(nans==1)=[];

x=1:length(dis_vec);
y=zeros(size(x));
% figure
% hold on
% plot(x,dis_vec,'k','linewidth',2)
% plot(x,y,'r--','linewidth',2)
% xlim([1 length(dis_vec)])

%%
dis_vec_ratios=triu(Dist_mat_ratios);%get rid of zeros and nans
dis_vec_ratios=reshape(dis_vec_ratios,1,144);
dis_vec_ratios(dis_vec_ratios==0)=[];
nans=isnan(dis_vec_ratios);
dis_vec_ratios(nans==1)=[];

x=1:length(dis_vec_ratios);
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

exp=Deltas_Dist_vec_norm
cont=exp;%Deltas_Dist_vec_norm_cont       %need to fix this issue not to do this manually

x1 = linspace(4,5,size(cont,2));
x2 = linspace(6,7,size(exp,2));
mean_vec_cont=mean(cont)*ones(size(x1))
mean_vec_exp=mean(exp)*ones(size(x2))
x3=3:8;
y=zeros(size(x3))

figure
xlim([3 8])
hold on
plot(x1,cont,'k*')
plot(x2,exp,'ko')
plot(x1,mean_vec_cont,'r--','linewidth',2)
plot(x2,mean_vec_exp,'r--','linewidth',2)
plot(x3,y,'k--')
XTick=[]


sem_cont=std(cont)/sqrt(length(cont))
sem_exp=std(exp)/sqrt(length(exp))

x_er=[4.5 6.5];
y_er=[mean(cont) mean(exp)];
sems=[sem_cont sem_exp];
errorbar(x_er,y_er,sems,'r.','linewidth',2)

box off

[h,p]=ttest2(exp,cont)
[h,p]=ranksum(exp,cont)

end

%% plot sds_mean
if run_no==2;
x=1:12;
figure
title('Variance by odor')
hold on
plot(x,all_sds{1},'b*')
plot(x,all_sds{2},'r*')

%%saving to compare, deltas norm by var
if condition(1)=='C'
Deltas_Dist_vec_norm_cont=Deltas_Dist_vec_norm; 
end

if condition(1)=='E'
Deltas_Dist_vec_norm_exp=Deltas_Dist_vec_norm; 
end

if condition(1)=='a'
Deltas_Dist_vec_norm_exp=Deltas_Dist_vec_norm; 
end

if condition(1)=='m'
Deltas_Dist_vec_norm_cont=Deltas_Dist_vec_norm; 
end

end
%%
end