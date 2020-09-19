
   clc
   clear all
   %close all
   
   data_for_perms={};
   %%IMPORTANT PARAMETERS to choose according to the experiment
   valve_timing=49;
   
%%
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
 mice=awake_mice;       %Fig 3 EXP
% mice=mice_awake_cont;  %Fig 3 CONT
% mice=first_timelapse;  %Fig 5 8WPI
% mice=sec_timelapse;    %Fig 5 12WPI
% mice=third_timelapse   %Fig 5 16WPI
% mice=cno_cont_mice;    %Fig S3   
%Comment/uncomment the relevant groups. each run could be conducted over 1
%group only
   
%anes exp
mice=[1,2,5,9,10,39,12,13,32,36];

fields_for_analysis=struct;
fields_for_analysis(1).fields_numbers=([1 2 3]);
fields_for_analysis(2).fields_numbers=([1 2]);
fields_for_analysis(3).fields_numbers=([1]);
fields_for_analysis(4).fields_numbers=([1]);
fields_for_analysis(5).fields_numbers=([1 2]);
fields_for_analysis(6).fields_numbers=([1 2]);
fields_for_analysis(7).fields_numbers=([1]);
fields_for_analysis(8).fields_numbers=([1]);
fields_for_analysis(9).fields_numbers=([1 2]);
fields_for_analysis(10).fields_numbers=([1]);

% %anes cont
%      fields_for_analysis=struct;
% mice=[30,34,35,37,38,11,7,15,16, 40];
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1]);
% fields_for_analysis(5).fields_numbers=([1 2]);
% fields_for_analysis(6).fields_numbers=([1 2]);
% fields_for_analysis(7).fields_numbers=([1 2]);
% fields_for_analysis(8).fields_numbers=([1 2]);
% fields_for_analysis(9).fields_numbers=([1]);
% fields_for_analysis(10).fields_numbers=([1]);
% 
% %awake saline cont
% mice=[1000 1001 1002 1003 1004]
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([2 3 4]);
% fields_for_analysis(2).fields_numbers=([1 2 3 4]);
% fields_for_analysis(3).fields_numbers=([1 2 3]);
% fields_for_analysis(4).fields_numbers=([1]);
% fields_for_analysis(5).fields_numbers=([1]);
% 
% %all awake exp mice first session EXP
% mice=[17,50,200,52,54]; 
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2 3]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1 2]);
% fields_for_analysis(5).fields_numbers=([1]);
% 
% % %timelapse1 full 
% mice=[100,200,300,500,600];
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1]);
% fields_for_analysis(5).fields_numbers=([1 3]);
% 
% % % % %timelapse2 full
% mice=[101,201,301,601];
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1 3]);
% 
% % % %timelapse3 full
% mice=[102,202,302,502];
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1]);
% 
% %CNO_CONT_WT_AWAKE
% mice=[1011, 10102];
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([2]);
% fields_for_analysis(2).fields_numbers=([1 2]);


    %% 11 odors, 7 hz, 5 reps
    %blank= 5;
    HZ= 7;
    numOfOdors = 12;  %including blank
    numOfTrails = 5;
    prestim_and_odor_to_final=7; %in seconds
    stim=2;
    rec_after_stim=7; %while 3 secs are for the gap so there are 10.
    
    lengthOfOdor=(prestim_and_odor_to_final+stim+rec_after_stim)*HZ;
    lengthOfTrails=lengthOfOdor*numOfOdors;
    startOfNeurons = 2; %Don't include all uneccesary columns  
    
    
    
 
    valves_names_12 ={'Valeraldehyde','Methyl propanoate','Ethyl acetate','Butyraldehyde','Blank','TMT','Female Pee', 'Male Pee', 'Peanut butter','Ethil tiglate','Propanal', 'Pups Bedding'}
   
    valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;    %%new orfer for valves
    ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
    
    for gg=1:length( ploting)
    k_for_plot(gg)=find( valves_12==ploting(gg));
    end

    
    valves= valves_12; 
    blank=find(valves==6);    %important!!!
    

for mouse=mice;  %so fields_for_analysis is the struct that organizes all data
j=find(mice==mouse);
fields_for_analysis(j).name=(['mouse' num2str(mouse)]);
end


j=0;
data=[];
for mouse=mice;  %for external loop that runs over mice

kk=0; 
j=find(mice==mouse);
fields_for_analysis(j).cells_per_field=zeros(size(fields_for_analysis(j).fields_numbers));
fields_for_analysis(j).number_of_fields=length(fields_for_analysis(j).cells_per_field);
total_data_before=[];
total_data_after=[];
    for f=(fields_for_analysis(j).fields_numbers);
        kk=kk+1;
   % CHANGE HERE ACCORDING TO YOUR PATH:
   mouse_path='F:\ALL DATA FOR PAPER HSN\ORGANIZED FOR PAPER\General_Physiology\M';  
    
    field_path='\before\before_f';
    filename=[mouse_path num2str(mouse) field_path num2str(f) '.xlsx'];
    data = xlsread(filename);
   % data =importdata(filename); data=data.data;
    data(:,1:(startOfNeurons-1))=[];  %%%%  important!!!!! cleaning the non flourascent values
    data_size=size(data,2);
    total_data_before=[total_data_before;data'];
    fields_for_analysis(j).data_before=total_data_before;
    fields_for_analysis(j).cells_per_field(kk)=data_size;
   
    field_path='\after\after_f';
    %field_path='\after\gc_after_f'; % for gcs
    filename=[mouse_path num2str(mouse) field_path num2str(f) '.xlsx'];
    data = xlsread(filename);
    %data =importdata(filename); data=data.data;
    data(:,1:(startOfNeurons-1))=[];            %%%%  important!!!!! cleaning the non flourascent values
    total_data_after=[total_data_after;data'];
    fields_for_analysis(j).data_after=total_data_after;
   %
       
    
    end
fields_for_analysis(j).data_before=fields_for_analysis(j).data_before';
fields_for_analysis(j).data_after=fields_for_analysis(j).data_after';

   end
   

     %%
      for run=1:2;
     for mouse=mice;
         mouse
     ll=find(mice==mouse);
     numOfNeurons = size(fields_for_analysis(ll).data_before,2); %same like data after   
     mean_traces_mat=zeros(numOfNeurons,numOfOdors,lengthOfOdor);
     DataMat=zeros(numOfNeurons,numOfTrails,numOfOdors,lengthOfOdor); 
     Dff=DataMat;
     SmoothDff=DataMat;
    
    %now this is one loop that organizes all the data in a 4-dim
    %matrix (with Ami's help)+ transforming to delta_f_over_f(normalizing)+ ing
  
    ii = 0; %index for the loop to allow overcoming the difference due to the fact that first 5 lines doesn't contain Neurons
    FramesForF0start =prestim_and_odor_to_final*HZ-HZ*5;
    FramesForF0end =prestim_and_odor_to_final*HZ-HZ*2;
    
    % FramesForF0end =prestim_and_odor_to_final*HZ-2*HZ; change to check
    % predictive coding phenomena
    
    data_before=fields_for_analysis(ll).data_before;
    data_after=fields_for_analysis(ll).data_after;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HERE AND
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AT THE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                       %%%%%%%%%%%%bef after
    %data = xlsread(filename);
    %data= shuffled_data;
    
    %data=data_after;
    
  
    if run==1
    data=data_before;
    end
    
    if run==2
    data=data_after;
    end
    
%     mmm=data-data_after;
%     nnn=data-data_before;
%     if sum(sum(mmm))==0;
%         time='after';
%     end
%     
%     if sum(sum(nnn))==0;
%         time='before';
%     end
    
    temp_raw_f_mat=[];
    
   
    startOfNeurons=1;   %because it is cleaned at the begining
    ii=0;
    for neuron = startOfNeurons:numOfNeurons + startOfNeurons - 1;
        ii = ii + 1;
        
        
        for kk = 1:numOfOdors
            for jj = 1:numOfTrails;
                startInd = (jj-1) * lengthOfTrails;%lengthOfOdor ;
                tempind = (kk-1) * lengthOfOdor + 1 : kk * lengthOfOdor;
                ind = startInd + tempind;
                DataMat(ii,jj,kk,:) = data(ind,neuron);
                f0 = mean(data(ind(FramesForF0start:FramesForF0end),neuron));
                Dff(ii,jj,kk,:) = (data(ind,neuron) - f0)./f0;
                SmoothDff(ii,jj,kk,:)= smooth(Dff(ii,jj,kk,:));
                 % SmoothDff(ii,jj,kk,:)= (Dff(ii,jj,kk,:));
               
            end
            c=squeeze(SmoothDff(ii,:,kk,:));
            c=mean(c);
            mean_traces_mat(ii,kk,:)=c;
            
        end
        
    end
      data_for_perms{run,ll}.SmoothDff={SmoothDff};
      data_for_perms{run,ll}.Mean=mean_traces_mat;
     end
 
     end

%% load data   and make this a function

resp_window=49:77;
baseline_window=7:35;

% mice=[100,101,102,103,104,105,200,201,202,203,204,205,300,301,302,303,304,305,400,17];
%mice=[100,101,102,103,200,201,202,203,300,301,302,303];

for run=1:2

i=0;
for mouse=mice  
mouse
    i=i+1;
    
  
   data=(data_for_perms{run,i}.Mean);   %so the data is taken from the mean trace only. change to run later on
   
   numOfNeurons = size(data,1); 
   numOfodors = size(data,2);
   
   perm_p_vals=zeros(size(data,1),size(data,2));
    perm_effect_size=zeros(size(data,1),size(data,2));
    perm_resp_magnitude=zeros(size(data,1),size(data,2));
    perm_obs_difference=zeros(size(data,1),size(data,2));
    
   for ii = 1:numOfNeurons;  
 ii
        for kk = 1:numOfodors;
   
       trace=squeeze(data(ii,kk,:));
       
       if mean(trace(resp_window))>=0;
        a=find(trace==max(trace(resp_window))); 
        a=a(length(a));                        %this is for rare cases in which there is more than 1 point
        sample1=trace(a-3:a+3);
        
        b=find(trace==max(trace(baseline_window)));  
        b=b(length(b));                        %this is for rare cases in which there is more than 1 point
        sample2=trace(b-3:b+3);
        
       end
       
         if mean(trace(resp_window))<=0                                            %so sample 1 is for each cell_odor, and is not time locked. just 7 frames around the peack
        a=find(trace==min(trace(resp_window)));    
        a=a(length(a));
        sample1=trace(a-3:a+3);
        
        b=find(trace==min(trace(baseline_window))); 
        b=b(length(b));                        %this is for rare cases in which there is more than 1 point
        sample2=trace(b-3:b+3);
        
         end
         
         [c,d,e] = permutationTest(sample1,sample2,1,'exact',1);%,'plotresult',1,'showprogress',3432)
        
         perm_resp_magnitude(ii,kk)=mean(sample1);
        
         perm_p_vals(ii,kk)=c;
         perm_obs_difference(ii,kk)=d;
         perm_effect_size(ii,kk)=e;
       
         

        end
   end
   
   
    
      
      SmoothDff=data_for_perms{run,i}.SmoothDff;
      mean_traces_mat=data_for_perms{run,i}.Mean;
      
      
if run==1
save(['perutation results- mouse' num2str(mouse) 'results before'], 'perm_resp_magnitude','perm_p_vals','perm_obs_difference','perm_effect_size','SmoothDff','mean_traces_mat' )
disp('saved before')
end

if run==2
save(['perutation results- mouse' num2str(mouse) 'results after'], 'perm_resp_magnitude','perm_p_vals','perm_obs_difference','perm_effect_size','SmoothDff','mean_traces_mat')
disp('saved after')
end

end
       
end

%% A new loop for before vs after desicion: run code one after the first perm part and then run this



resp_window=49:77;
baseline_window=7:35;

% mice=[17,50,51,52,30,34,35,37,38,11,7,15,16,40,60,61,62,1,2,5,9,10,39,12,13,32,36];
% mice=[100,101,102,103,104,105,200,201,202,203,204,205,300,301,302,303,304,305,400,17];
% mice=[100];
%  mice=[100,101,102,103,200,201,202,203,300,301,302,303];
%  mice=[17];

i=0;

for mouse=mice  
       i=i+1;
%     mouse=100 
%     i=1   

       perm_table{1,i} = load(['perutation results- mouse' num2str(mice(i)) 'results before']); 
       perm_table{2,i} = load(['perutation results- mouse' num2str(mice(i)) 'results after']); 
            perm_p_vals_bef_vs_aft=zeros(size(perm_table{1,i}.perm_p_vals));
            perm_obs_difference_bef_vs_aft=zeros(size(perm_table{1,i}.perm_p_vals));
            perm_effect_size_bef_vs_aft=zeros(size(perm_table{1,i}.perm_p_vals));
            
            bef=  perm_table{1,i}.SmoothDff; %this is to compare before vs after for significance
            aft=  perm_table{2,i}.SmoothDff;
            bef=bef{1,1};
            aft=aft{1,1};
             numOfNeurons=size(bef,1);
             numOfodors=size(bef,3);
      for ii = 1:numOfNeurons;  
          ii
        for kk = 1:numOfodors;
         
            sample1=mean(bef(ii,:,kk,resp_window),4);
            sample2=mean(aft(ii,:,kk,resp_window),4);
           
            [c,d,e] = permutationTest(sample1,sample2,1,'exact',1);
            perm_p_vals_bef_vs_aft(ii,kk)=c;
            perm_obs_difference_bef_vs_aft(ii,kk)=d;
            perm_effect_size_bef_vs_aft(ii,kk)=e;
        end
      end
  perm_resp_delta_at_top=perm_table{1,i}.perm_resp_magnitude;- perm_table{2,i}.perm_resp_magnitude;
  save(['perutation results- mouse' num2str(mouse) ' bef_vs_after'], 'perm_resp_delta_at_top','perm_p_vals_bef_vs_aft','perm_obs_difference_bef_vs_aft','perm_effect_size_bef_vs_aft')    
end