    %%cell to cell analysis- 2 places to change before and after NO NEED TO
    %%please read me!!!!!!
    close all
    
    
      % close all
  
    for run=1:2 % added for before after together
   
   clc
   clearvars -except run

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
    startOfNeurons = 2; %Don't include all uneccesary columns  
    
 
    valves_names_12 ={'Valeraldehyde','Methyl propanoate','Ethyl acetate','Butyraldehyde','Blank','TMT','Female Pee', 'Male Pee', 'Peanut butter','Ethil tiglate','Propanal', 'Pups Bedding'}
   
    valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;    %%new orfer for valves
    ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
    
    for gg=1:length( ploting)
    k_for_plot(gg)=find( valves_12==ploting(gg));
    end

    valves= valves_12; 
    blank=find(valves==6);    %important!!!


 
%%
%DEFINE MICE AND STRUCTURE!   %manually change and uncomment, according to
%the group you would like to observe, as described above each part
% 
 
%timelapse1 full Fig5 8WPI
% mice=[100,200,300,500,600];
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1]);
% fields_for_analysis(5).fields_numbers=([1 3]);

% 
% 
% %timelapse2 full Fig5 12WPI
% mice=[101,201,301,601];
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1 3]);
% % % 

% %timelapse3 full Fig5 16WPI
% mice=[102,202,302,502];
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1 2]);
% fields_for_analysis(4).fields_numbers=([1]);

% 
 mice=[17,50,200,52,54];  %all awake exp mice first session Fig 3 EXP
fields_for_analysis=struct;
fields_for_analysis(1).fields_numbers=([1 2 3]);
fields_for_analysis(2).fields_numbers=([1 2]);
fields_for_analysis(3).fields_numbers=([1 2]);
fields_for_analysis(4).fields_numbers=([1 2]);
fields_for_analysis(5).fields_numbers=([1]);
% 
% %all awake control before after saline first time point Fig 3 CONT
% mice=[1000 1001 1002 1003 1004]
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([2 3 4]);
% fields_for_analysis(2).fields_numbers=([1 2 3 4]);
% fields_for_analysis(3).fields_numbers=([1 2 3]);
% fields_for_analysis(4).fields_numbers=([1]);
% fields_for_analysis(5).fields_numbers=([1]);

% % % all anes  mice positive geno bef/aft CNO Fig 2 EXP
% % % 
% mice=[1,2,5,9,10,39,12,13,32,36];
% 
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2 3]);
% fields_for_analysis(2).fields_numbers=([1 2]);
% fields_for_analysis(3).fields_numbers=([1]);
% fields_for_analysis(4).fields_numbers=([1]);
% fields_for_analysis(5).fields_numbers=([1 2]);
% fields_for_analysis(6).fields_numbers=([1 2]);
% fields_for_analysis(7).fields_numbers=([1]);
% fields_for_analysis(8).fields_numbers=([1]);
% fields_for_analysis(9).fields_numbers=([1 2]);
% fields_for_analysis(10).fields_numbers=([1]);

% % % % all anes  mice negative geno bef/aft CNO Fig 2 CONT
% fields_for_analysis=struct;
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

% %cno control Fig. S3
% mice=[10102 1011]
% fields_for_analysis=struct;
% fields_for_analysis(1).fields_numbers=([1 2]);
% fields_for_analysis(2).fields_numbers=([2]);

%% added at 27.1.19 to upload the permutation test results (run without the test first, then run the test code and than run with this part) 

%make a function here to run only on the first time

num_of_mice=length(mice)
for i=1:num_of_mice;

perm_table{1,i} = load(['perutation results- mouse' num2str(mice(i)) 'results before']); %create a struct- rows for conditions, coloumns for mice
perm_table{2,i} = load(['perutation results- mouse' num2str(mice(i)) 'results after']); 
%perm_table{3,i} = load(['perutation results- mouse' num2str(mice(i)) ' bef_vs_after']); 

end


%% first path is to open all data and save data before and after

j=0;
for mouse=mice;  %so fields_for_analysis is the struct that organizes all data
j=find(mice==mouse);
fields_for_analysis(j).name=(['mouse' num2str(mouse)]);
end


j=0;
data=[];
for mouse=mice;  %for external loop that runs over mice
mouse
kk=0; 
j=find(mice==mouse);
fields_for_analysis(j).cells_per_field=zeros(size(fields_for_analysis(j).fields_numbers));
fields_for_analysis(j).number_of_fields=length(fields_for_analysis(j).cells_per_field);
total_data_before=[];
total_data_after=[];
    for f=(fields_for_analysis(j).fields_numbers);
        kk=kk+1;
    %Change path here if needed
    
    mouse_path='F:\ALL DATA FOR PAPER HSN\ORGANIZED FOR PAPER\General_Physiology\M';
   %  mouse_path='C:\2p analysis finalsham\M';
   
   
    field_path='\before\before_f';
    %field_path='\before\gc_before_f'; %for gcs
    filename=[mouse_path num2str(mouse) field_path num2str(f) '.xlsx'];
    data = xlsread(filename);
   % data =importdata(filename); data=data.data;

             %  %%%%  important!!!!! cleaning the non flourascent values
    data(:,1:(startOfNeurons-1))=[];  
    data_size=size(data,2);
    total_data_before=[total_data_before;data'];
    fields_for_analysis(j).data_before=total_data_before;
    fields_for_analysis(j).cells_per_field(kk)=data_size;
   
    field_path='\after\after_f';
    %field_path='\after\gc_after_f'; % for gcs
    filename=[mouse_path num2str(mouse) field_path num2str(f) '.xlsx'];
    %data = xlsread(filename);
    data =importdata(filename); data=data.data;
    data(:,1:(startOfNeurons-1))=[];            %%%%  important!!!!! cleaning the non flourascent values
    total_data_after=[total_data_after;data'];
    fields_for_analysis(j).data_after=total_data_after;
   %
       
    
    end
fields_for_analysis(j).data_before=fields_for_analysis(j).data_before';
fields_for_analysis(j).data_after=fields_for_analysis(j).data_after';


    end


     %%
     for mouse=mice;
    ll=find(mice==mouse);
    numOfNeurons = size(fields_for_analysis(ll).data_before,2); %same like data after   
    mean_traces_mat=zeros(numOfNeurons,numOfOdors,lengthOfOdor);
    mean_traces_mat_no_smooth=zeros(numOfNeurons,numOfOdors,lengthOfOdor);
    DataMat=zeros(numOfNeurons,numOfTrails,numOfOdors,lengthOfOdor); 
    Dff=DataMat;
    SmoothDff=DataMat; 
     SmoothDff_pca=DataMat;
    ii = 0; %index for the loop to allow overcoming the difference due to the fact that first 5 lines doesn't contain Neurons
    FramesForF0start =prestim_and_odor_to_final*HZ-HZ*5; % use this!!!
  
    FramesForF0end =prestim_and_odor_to_final*HZ-2*HZ; %change to check
    
    
    data_before=fields_for_analysis(ll).data_before;
    data_after=fields_for_analysis(ll).data_after;

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
     
    %% determine which responses are significant-peak
    
    %z scores- (a) define std as the std of the blank during the response window (b) measure the distance between peack of blank to peack of odor
    % (c) check how many times sd feets in this distance and create a matrix of
    % z scores
    %(d) set a threshold (e) add astrix
    
    %this whole thing is not in use anymore, thus marked as a comment.
 
    ii=0;
     
    Zscores=zeros(numOfNeurons,numOfTrails+1,numOfOdors);
    New_Z=zeros(numOfNeurons,numOfTrails+1,numOfOdors);
    Zscores_inhib=Zscores;
    response_magnitude=zeros(numOfNeurons,numOfOdors);
    response_magnitude_by_integral=zeros(numOfNeurons,numOfOdors);
    response_magnitude_inhib= response_magnitude;
    mean_response_magnitude_no_blankcontrol=zeros(numOfNeurons,numOfOdors);
    response_magnitude_at_peak=zeros(numOfNeurons,numOfOdors);
   
    new_t_score=zeros(numOfNeurons,numOfOdors);
    new_t_desicion= new_t_score;
    temp_exit_inhib=zeros(numOfTrails,2);   %save means for the ttest
     t_Thersh=0.05;    %changed from .05 at 25.7.18
     F_thresh=0.05;
    new_f_score=zeros(numOfNeurons,numOfTrails+1,numOfOdors);
    
    for ii=1:numOfNeurons;
        for kk=1:1:numOfOdors
        STD_MEAN= std(mean_traces_mat(ii,kk,FramesForF0start:FramesForF0end-7));  %%changed from 3 to 7 at 11.1.19
        STDMEANMAT(ii,kk,:)=STD_MEAN;
        for jj = 1:numOfTrails;
            STD= std(SmoothDff(ii,jj,kk,FramesForF0start:FramesForF0end-7));   %%changed from 3 to 7 at 11.1.19              %the statisticas are ancored to the response for the blank, each trial has its own noise
            STDMATRIX(ii,kk,jj)= STD;
            
            
            %%% important- note that the blank number is crucial here. 5 in
            %%% this case but needs to be changed
            %so now blanks are not taken into acount when determining
            %significance/magnitudes. compare to 0, since f0 is 4 secs
            %before valve opening.
         
            %   for kk = 1:numOfOdors;
                Zscores(ii,jj,kk,:) = (max(abs(SmoothDff(ii,jj,kk,startwindow:endwindow))))./STD;%-max(SmoothDff(ii,jj,blank,startwindow:endwindow)))./STD;
                New_Z(ii,jj,kk,:) = (mean(SmoothDff(ii,jj,kk,startwindow:endwindow)))/STD;
               
                
                %for inhibition:
                Zscores_inhib(ii,jj,kk,:) = (min(SmoothDff(ii,jj,kk,startwindow:endwindow)))./STD;%-max(SmoothDff(ii,jj,blank,startwindow:endwindow)))./STD;
               
              
                
               %column 1 in temp is the response. column 2 is the prestim.
               %next we'll use it for a ttest
               a=SmoothDff(ii,jj,kk,startwindow:endwindow);
               b=SmoothDff(ii,jj,kk,FramesForF0start:FramesForF0end-3);
               
                temp_exit_inhib(jj,1)=mean(a);
                temp_exit_inhib(jj,2)=mean(b);
              
                [h,p] = vartest2(a,b);
                if p<F_thresh
                new_f_score(ii,jj,kk)=1;
                else 
                new_f_score(ii,jj,kk)=0;
                end
               
        end
         [h,p] = vartest2(mean_traces_mat(ii,kk,FramesForF0start:FramesForF0end-3),mean_traces_mat(ii,kk,startwindow:endwindow));
          if p<F_thresh
                new_f_score(ii,jj+1,kk)=1;
                else 
                new_f_score(ii,jj+1,kk)=0;
                end
      
        [h,p]=ttest(temp_exit_inhib(:,1),temp_exit_inhib(:,2));
         new_t_score(ii,kk)=p;
         if new_t_score(ii,kk)<t_Thersh
           new_t_desicion (ii,kk)=1;  
         else
             new_t_desicion (ii,kk)=0;
         end
             
         Zscores_inhib(ii,numOfTrails+1,kk,:)=(min(mean_traces_mat(ii,kk,startwindow:endwindow)))./STD_MEAN;     
         Zscores(ii,numOfTrails+1,kk,:)=(max(abs(mean_traces_mat(ii,kk,startwindow:endwindow))))./STD_MEAN; 
         New_Z(ii,numOfTrails+1,kk,:)=(mean(mean_traces_mat(ii,kk,startwindow:endwindow)))/STD_MEAN;
 
         response_magnitude_inhib(ii,kk)=min(mean_traces_mat(ii,kk,startwindow:endwindow));
         response_magnitude(ii,kk)=max(mean_traces_mat(ii,kk,startwindow:endwindow));
         response_magnitude_by_integral(ii,kk)=mean(mean_traces_mat(ii,kk,startwindow:endwindow));
         response_magnitude_only_positive=response_magnitude;
         
%          %here is the condition for a response to count as inhibitory:


  
          if Zscores_inhib(ii,numOfTrails+1,kk,:)<-4;
          if abs(Zscores_inhib(ii,numOfTrails+1,kk,:))>=Zscores(ii,numOfTrails+1,kk,:)
               if  abs(response_magnitude_inhib(ii,kk))>abs(response_magnitude(ii,kk)); 
                   if response_magnitude_by_integral(ii,kk)<=0;   %
                       response_magnitude(ii,kk)=response_magnitude_inhib(ii,kk);
                       Zscores(ii,:,kk,:)=Zscores_inhib(ii,:,kk,:);
                       
                   end
                   
               end 
          end        
          end
         

                if  response_magnitude_inhib(ii,kk)>0; %in order to avoid negative responses
                    response_magnitude_inhib(ii,kk)=0;
                end
      
                
        end
        end
   % end
    


    
    %%
    %determining
%     threshold=2.7; %z score for significane
    New_Z_thresh=3;
    astrix_mat=zeros(numOfNeurons,numOfTrails+1, numOfOdors); %or size zscores
%     astrix_desicion=zeros(numOfNeurons,1, numOfOdors);
%     F_desicion= astrix_desicion;
    New_Z_mat=astrix_mat;
  
     
    for ii = 1:numOfNeurons;
        for kk=1:numOfOdors;
            for jj=1:(numOfTrails+1)
                
%                 if abs(Zscores(ii,jj,kk))>threshold;
%                     astrix_mat(ii,jj,kk)=1;
%                 end
                
                if abs(New_Z(ii,jj,kk))>New_Z_thresh
                   New_Z_mat(ii,jj,kk)=1;
                 
                         if New_Z(ii,jj,kk)<0
                        New_Z_mat(ii,jj,kk)=-1;
                         end
                   
                end
                
            end
            
%             if astrix_mat(ii,numOfTrails+1,kk)==1;      %mean condition
%                 astrix_desicion(ii,kk)=1;
%             end
%             if sum(astrix_mat(ii,1:numOfTrails,kk))<numOfTrails-2 %at least 3/5 significant responses
%                 astrix_desicion(ii,kk)=0;
%             end
%             
%             if new_f_score(ii,numOfTrails+1,kk)==1;      %mean condition
%                 F_desicion(ii,kk)=1;
%             end
%             if sum(new_f_score(ii,1:numOfTrails,kk))<numOfTrails-2 %at least 3/5 significant responses
%                 F_desicion(ii,kk)=0;
%             end 
%             
        end
        
    end
    New_Z_mean=squeeze(New_Z(:,numOfTrails+1,:));
    New_Z_mean_org=(New_Z_mean(:,k_for_plot));   % just to monitor, no use
    
    effect_size_thresh=5;
    perm_mat=perm_table{run,ll}.perm_effect_size;
    perm_in_ex=sign(perm_mat);
    perm_thresh_mat=ones(size(perm_mat))* effect_size_thresh;
    perm_des=abs(perm_mat)>=perm_thresh_mat;
    perm_des(:,blank)=0;
    perm_des=perm_des.* perm_in_ex;
    
%     New_Z_mean_des=squeeze(New_Z_mat(:,numOfTrails+1,:));
%     New_Z_mean_des(:,blank)=0;
    
    
    
     %meta_desicion=
     %new_t_desicion+squeeze(F_desicion)+squeeze(astrix_desicion)+New_Z_mean;
     %meta_desicion=New_Z_mean_des; OLD TESTS!
     
     meta_desicion=perm_des;
     meta_desicion(:,blank)=0;
     
     resp_magnitude_by_perm=perm_table{run,ll}.perm_resp_magnitude;
     
    %astrix_desicion(:,:,1)=[];  %excluding the blank
       % astrix_desicion(:,:,9)=[]; %ignoring the 9th odour
    %significant_responses_magnitudes=response_magnitude(meta_desicion>1);
    significant_responses_magnitudes=response_magnitude(meta_desicion~=0);
   %  significant_responses_magnitudes_by_integral=response_magnitude_by_integral(meta_desicion>1);
    significant_responses_magnitudes_by_integral=response_magnitude_by_integral(meta_desicion~=0);
    significant_resp_magnitude_by_perm=resp_magnitude_by_perm(meta_desicion~=0);
    %significat_responses_locations=find(meta_desicion>1);
    significat_responses_locations=find(meta_desicion~=0);

    
%     
%   %%  plotting- not crucial for the code to run
% % % %     
%     x=zeros(20,1);
%     x(20)=1;
%     size(SmoothDff);
%     % should be :   neurons     trials     oders    frames
%     NeuronsForPlot=8;
%     
%     %newk for changing the order in the plot
% %   
% % %     
%     for m=(0:NeuronsForPlot:numOfNeurons-NeuronsForPlot); %m is a vector defined in jumps of "Neurons for plot"
%         figure
% %         title(sprintf('Cells=%d-%d',m,6+m));
%         kk=-1;
%         x=ones(20,1);
%         y1=20;
%         y2=40;
%         y3=60;
%         y4=80;
%         
%         for ii=1+m:NeuronsForPlot+m; %ii for neuron
%             kk=kk+1;
%             dummy=0;
%              for k=1:numOfOdors%%for; ;% for odor k_for_plot %
%                dummy=dummy+1;
% %            subaxis(NeuronsForPlot,9, k+(9*(kk)),'SpacingHoriz', 0,'SpacingVert',0) 
% 
%              subplot(NeuronsForPlot,numOfOdors, dummy+(numOfOdors*(kk)));
% 
%                 plot(squeeze(SmoothDff(ii,:,k,:))','k');
% %                 axis off 
%                 hold on
%                 plot((squeeze(mean_traces_mat(ii,k,:))'),'r');
%                 
%                 ylim([-0.2 1.1])
%                 xlim([0 lengthOfOdor])
% %                 if astrix_desicion(ii,1,k)==1;
% %                     hold on
% %                     plot(y1,x(20),'* r')
% %                     
% %                 end
% %                  if new_t_desicion(ii,k)==1;
% %                     hold on
% %                     plot(y2,x(10),'* g')
% %                     
% %                  end
% %                  if F_desicion(ii,k)==1;
% %                     hold on
% %                     plot(y3,x(5),'* b')
% %                     
% %                  end
% %                 
% %                  if New_Z_mean_des (ii,k)==1
% %                      hold on
% %                       plot(y4,x(15),'* k')
% %                  end
% %                  
% %                   if New_Z_mean_des (ii,k)==-1
% %                      hold on
% %                       plot(y4,x(15),'o k')
% %                  end
% %                  
% 
%             end
%         end
%     end
     
%%
    
     
    
%     % plot individual trails     %uncomment this part if you'd like to
%     watch closely specific cell-odor pairs ( first neuron than odor)
%


%     c=squeeze(DataMat(1,:,6,:));
%     c=mean(c);
%     
%     figure;
%     
%     plot all trials for oder 3 for neuron number 20:
%     plot(squeeze(DataMat(2,:,1,:))')
%     ylim([-0.2 1.7]);
%     xlim([0 72]);
%     hold on
%     plot(c,'g')
%     
%     first neuron than odor
%     c=squeeze(resamp_data_mat(4,:,9,:));
%     c=mean(c);
%     
%     figure;
%     
%     plot all trials for oder 3 for neuron number 20:
%     plot(squeeze(resamp_data_mat(4,:,9,:))')
%     ylim([-0.2 1.7]);
%     xlim([0 72]);
%     hold on
%     plot(c,'g')
%     
    %% plotting proportion of cells+comperative analysis
 %% plotting proportion of cells+comperative analysis

 

    meta_desicion_cleaned=meta_desicion;
     meta_desicion_cleaned(:,blank)=[];
  
  
   
     prop_of_cells=sum(meta_desicion_cleaned~=0,2); %how many odores does a cell respond to
      prop_of_cells_inhib=sum(meta_desicion_cleaned==-1,2);
      prop_of_cells_excit=sum(meta_desicion_cleaned==1,2);
      
    num_of_cells_responding_to_odor=sum(meta_desicion_cleaned~=0,1); %how many cells responding to a given odor
   
   
     %how many cells responding to a given odor:
    x=(0:numOfOdors-1);
    figure;
    hold on
    title([num2str(run)])
    hist(prop_of_cells,x);
    A=hist(prop_of_cells,x);
    
    
    CELLS_RESP=fliplr(A);
    A=fliplr(A);
    A=A/numOfNeurons; %SO A gives the prob to belong to a group
    for ii=1:numOfOdors
        CDF(ii)=sum(A(1:ii));
    end
    
%     figure;
%     plot(fliplr(x),CDF);
%     set(gca,'XDir','reverse');
%    
    disp thanks
    
%     figure
%     hist(astrix_desicion)
%% save the magnitudes of each trial for ttest between before and after
 mat_for_diffrence_test=zeros(numOfNeurons,numOfTrails,numOfOdors);
for ii = 1:numOfNeurons;
        for kk=1:numOfTrails;
            for jj=1:numOfOdors;
             mat_for_diffrence_test(ii,kk,jj)= mean(SmoothDff(ii,kk,jj,startwindow:endwindow));  
            end
        end
end


 
%% save all unpaired trials and stats for both inhib, excit and general

    significat_responses_locations=(meta_desicion~=0);
    sig_inhib_locations=(meta_desicion==-1);
    sig_exc_locations=(meta_desicion==1);
    
    all_prop_for_odor=sum(significat_responses_locations,1);
    inhib_prop_for_odor=sum(sig_inhib_locations,1);
    exc_prop_for_odor=sum(sig_exc_locations,1);
    
    
    all_sig_trials_unpaired=zeros(size(mean_traces_mat));
    sig_inhib_trials_unpaired=zeros(size(mean_traces_mat));
    sig_excit_trials_unpaired=zeros(size(mean_traces_mat));
    
    

for ii = 1:numOfNeurons;
        for kk=1:numOfOdors
            
        if significat_responses_locations(ii,kk)~=0;
            all_sig_trials_unpaired(ii,kk,:)=mean_traces_mat(ii,kk,:);
        end
        
            if  sig_inhib_locations(ii,kk)~=0;
                
                 sig_inhib_trials_unpaired(ii,kk,:)=mean_traces_mat(ii,kk,:);
            end
            
            if sig_exc_locations(ii,kk)~=0;
                   
                  sig_excit_trials_unpaired(ii,kk,:)=mean_traces_mat(ii,kk,:);
            end
            
            
            end
end

        


%%
if time(1)=='b'
save(['cell to cell- mouse' num2str(mouse) 'results before'], 'CELLS_RESP','numOfNeurons','response_magnitude','significant_responses_magnitudes','significat_responses_locations','prop_of_cells','num_of_cells_responding_to_odor','SmoothDff','mean_traces_mat','new_t_desicion','meta_desicion','mat_for_diffrence_test','fields_for_analysis','response_magnitude_by_integral','significant_responses_magnitudes_by_integral','temp_raw_f_mat','spont_est','sig_excit_trials_unpaired','sig_inhib_trials_unpaired','all_prop_for_odor','inhib_prop_for_odor','exc_prop_for_odor','resp_magnitude_by_perm','prop_of_cells_excit','prop_of_cells_inhib','SmoothDff_pca','mean_traces_mat_no_smooth')
disp('saved before')
end

if time(1)=='a'
save(['cell to cell- mouse' num2str(mouse) 'results after'], 'CELLS_RESP','numOfNeurons','response_magnitude','significant_responses_magnitudes','significat_responses_locations','prop_of_cells','num_of_cells_responding_to_odor','SmoothDff','mean_traces_mat','new_t_desicion','meta_desicion','mat_for_diffrence_test','fields_for_analysis','response_magnitude_by_integral','significant_responses_magnitudes_by_integral','temp_raw_f_mat','spont_est','sig_excit_trials_unpaired','sig_inhib_trials_unpaired','all_prop_for_odor','inhib_prop_for_odor','exc_prop_for_odor','resp_magnitude_by_perm','prop_of_cells_excit','prop_of_cells_inhib','SmoothDff_pca','mean_traces_mat_no_smooth')
disp('saved after')
end
disp(['mouse=' num2str(mouse)])


     end
    
    end        
   
    
    % save things that are common to both before after no neeed to relate
    % to group identity yet
    