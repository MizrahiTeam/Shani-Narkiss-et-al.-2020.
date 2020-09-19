%% load traces to estimate background code. this is code 1/3
%make sure that mice has common mice table
    clear all
    clc
    close all
    
exp_mice=[17,50,200,52,54];
saline_cont=[1000 1001 1002 1003 1004];
cno_cont=[10102 1011];

    
   
 mice=[17,50,200,52,54];  %all awake exp mice first session
BG_green_STRUCT=struct;
BG_green_STRUCT(1).fields_numbers=([1 2 3]);
BG_green_STRUCT(2).fields_numbers=([1 2]);
BG_green_STRUCT(3).fields_numbers=([1 2]);
BG_green_STRUCT(4).fields_numbers=([1 2]);
BG_green_STRUCT(5).fields_numbers=([1]);

% mice=[1000 1001 1002 1003 1004]
% BG_green_STRUCT=struct;
% BG_green_STRUCT(1).fields_numbers=([2 3 4]);
% BG_green_STRUCT(2).fields_numbers=([1 2 3 4]);
% BG_green_STRUCT(3).fields_numbers=([1 2 3]);
% BG_green_STRUCT(4).fields_numbers=([1]);
% BG_green_STRUCT(5).fields_numbers=([1]);
% % % 

% %cno control
% mice=[10102 1011]
% BG_green_STRUCT=struct;
% BG_green_STRUCT(1).fields_numbers=([1 2]);
% BG_green_STRUCT(2).fields_numbers=([2]);


  
    for run=1:2; % added for before after together
  disp(run)
   clc
   clearvars -except run BG_green_STRUCT mice exp_mice saline_cont cno_cont
    startOfNeurons=1;


    %%IMPORTANT PARAMETERS to choose according to the experiment
    %valve_timing=28;
    valve_timing=49;

    
%   %% 11 odors, 7 hz, 5 reps
    HZ= 7;
    numOfOdors = 12;  %including blank
    numOfTrails = 5;
    prestim_and_odor_to_final=7; %in seconds
    stim=2;
    rec_after_stim=7; %while 3 secs are for the gap so there are 10.
    
    lengthOfOdor=(prestim_and_odor_to_final+stim+rec_after_stim)*HZ;
    lengthOfTrails=lengthOfOdor*numOfOdors;
    
 
   
    valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;    %%new orfer for valves
    ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
    
 
%%
%DEFINE MICE AND STRUCTURE!   %manually insert fields numbers for analysis
% 


% %allmice with bg's:

     
%% first path is to open all data and save data before and after

j=0;
for mouse=mice;  %so fields_for_analysis is the struct that organizes all data
j=find(mice==mouse);
BG_green_STRUCT(j).name=(['mouse' num2str(mouse)]);
end


j=0;
data=[];
for mouse=mice;  %for external loop that runs over mice
mouse
kk=0; 
j=find(mice==mouse);
BG_green_STRUCT(j).cells_per_field=zeros(size(BG_green_STRUCT(j).fields_numbers));
BG_green_STRUCT(j).number_of_fields=length(BG_green_STRUCT(j).cells_per_field);
total_data_before=[];
total_data_after=[];
    for f=(BG_green_STRUCT(j).fields_numbers);
        kk=kk+1;
    %Change path here if needed
    
     mouse_path='F:\ALL DATA FOR PAPER HSN\ORGANIZED FOR PAPER\General_Physiology\all_physiology_data\M';
%     mouse_path='E:\2p analysis final_sun_haran\M';
   %  mouse_path='C:\2p analysis finalsham\M';
   
   
    field_path='\before\before_f';
    %field_path='\before\gc_before_f'; %for gcs
    filename=[mouse_path num2str(mouse) field_path num2str(f) '_bg.xlsx'];
    data= xlsread(filename);
 %   data =importdata(filename); %data=data.data;
  %  a=str2num(data{4,1})
% s=data(2:end,1)
% 
% sss=str2num(s)
% data(1,:)=[]
% s=str2num(cell2mat(data(:,1)'))

             %  %%%%  important!!!!! cleaning the non flourascent values
    data=data(:,size(data,2));
    data_size=size(data,2);
    total_data_before=[total_data_before;data'];
    BG_green_STRUCT(j).data_before=total_data_before;
    BG_green_STRUCT(j).cells_per_field(kk)=data_size;
   
    field_path='\after\after_f';
    %field_path='\after\gc_after_f'; % for gcs
    filename=[mouse_path num2str(mouse) field_path num2str(f) '_bg.xlsx'];
    data = xlsread(filename);
    data=data(:,size(data,2));
%     data =importdata(filename); data=data.data;
%     data(:,1:(startOfNeurons-1))=[];            %%%%  important!!!!! cleaning the non flourascent values

    total_data_after=[total_data_after;data'];
    BG_green_STRUCT(j).data_after=total_data_after;
   %
       
    
    end
BG_green_STRUCT(j).data_before=BG_green_STRUCT(j).data_before';
BG_green_STRUCT(j).data_after=BG_green_STRUCT(j).data_after';


    end


     %%
     for mouse=mice;
    ll=find(mice==mouse);
    numOfNeurons = size(BG_green_STRUCT(ll).data_before,2); %same like data after   
    mean_traces_mat=zeros(numOfNeurons,numOfOdors,lengthOfOdor);
    mean_traces_mat_no_smooth=zeros(numOfNeurons,numOfOdors,lengthOfOdor);
    DataMat=zeros(numOfNeurons,numOfTrails,numOfOdors,lengthOfOdor); 
    Dff=DataMat;
    SmoothDff=DataMat; 
     SmoothDff_pca=DataMat;
    ii = 0; %index for the loop to allow overcoming the difference due to the fact that first lines doesn't contain Neurons
    FramesForF0start =prestim_and_odor_to_final*HZ-HZ*5; % use this!!!
   %FramesForF0start =90;
   %FramesForF0end =prestim_and_odor_to_final*HZ-HZ; used that previously
    FramesForF0end =prestim_and_odor_to_final*HZ-2*HZ; %change to check
    %predictive coding phenomena use this!!!
    %FramesForF0end =105;
    
    data_before=BG_green_STRUCT(ll).data_before;
    data_after=BG_green_STRUCT(ll).data_after;

    %load the relevant data according to run
    if run==1
    data=data_before;
    end
    
    if run==2
    data=data_after;
    end
    
    mmm=data-data_after;
    nnn=data-data_before;
    if sum(sum(mmm))==0;
        time='after';
    end
    
    if sum(sum(nnn))==0;
        time='before';
    end
    
    temp_raw_f_mat=[];
    
    %% now this is one loop that organizes all the data in a 4-dim
    %matrix (with Ami's help)+ transforming to delta_f_over_f(normalizing)+ ing
    
    startOfNeurons=1;   %because it is cleaned at the begining
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
                 SmoothDff_pca(ii,jj,kk,:)=smooth(Dff(ii,jj,kk,:),28);
                
                 % SmoothDff(ii,jj,kk,:)= (Dff(ii,jj,kk,:));
               temp_raw_f_mat(ii,jj,kk)= mean(DataMat(ii,jj,kk,FramesForF0start:FramesForF0end));
               spont_est(ii,jj,kk)=std(Dff(ii,jj,kk,FramesForF0start:FramesForF0end+HZ));
            end
            c=squeeze(SmoothDff(ii,:,kk,:));
            c=mean(c);
            mean_traces_mat(ii,kk,:)=c;
           
            d=squeeze(Dff(ii,:,kk,:));            %this variable is used for plotting the responses over time
            d=mean(d);
            mean_traces_mat_no_smooth(ii,kk,:)=d;
            
        end
        
    end
        
%   % data for cleaning loop deleted
%define the response window

     startwindow=prestim_and_odor_to_final*HZ;%prestim_and_odor_to_final*HZ-HZ; %window for response so now its 42
     endwindow=startwindow+(stim*HZ)+3*HZ; %
     
   

        



    if time(1)=='b'; % run==1
    
   BG_green_STRUCT(ll).raw_data_mat_before=DataMat;
   BG_green_STRUCT(ll).SmoothDff_before=SmoothDff;
    end
    
    if time(1)=='a'; %run==2
   
      
   BG_green_STRUCT(ll).raw_data_mat_after=DataMat;
   BG_green_STRUCT(ll).SmoothDff_after=SmoothDff;
    end
 
     end
   
    end        
   
    if length(mice)==length(exp_mice);
    if sum(mice-exp_mice)==0;
        condition='exp_mice';
            save(['BG_green_STRUCT_' condition],'BG_green_STRUCT')

    end
    end

    if length(mice)==length(saline_cont);
    if sum(mice-saline_cont)==0;
        condition='saline_cont';
            save(['BG_green_STRUCT_' condition],'BG_green_STRUCT')

    end
    end

    if length(mice)==length(cno_cont);
    if sum(mice-cno_cont)==0;
        condition='cno_cont';
            save(['BG_green_STRUCT_' condition],'BG_green_STRUCT')

    end
    end


    
    % save things that are common to both before after no neeed to relate
    % to group identity yet\
        save(['BG_green_STRUCT'],'BG_green_STRUCT')

    