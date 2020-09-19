     %% this code is to plot responses before and after on top of each other
     
     clc
     clear all
     close all
     
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

% mice=cont_mice;         %Fig 2 CONT
  mice=exp_mice;         %Fig 2 EXP
% mice=awake_mice;       %Fig 3 EXP
% mice=mice_awake_cont;  %Fig 3 CONT
% mice=first_timelapse;  %Fig 5 8WPI
% mice=sec_timelapse;    %Fig 5 12WPI
% mice=third_timelapse   %Fig 5 16WPI
% mice=cno_cont_mice;    %Fig S3 
 %%

     blank= 5;
     HZ= 7;
     numOfOdors =12;  %including blank
     numOfTrails = 5;
     prestim_and_odor_to_final=7; %in seconds
     stim=2;
     rec_after_stim=7; %while 3 secs are for the gap so there are 10.
     lengthOfOdor=(prestim_and_odor_to_final+stim+rec_after_stim)*HZ;
     lengthOfTrails=lengthOfOdor*numOfOdors;
     startOfNeurons = 1; %incluse all uneccesary columns  
     
     
     valves_names_12 ={'Valeraldehyde','Methyl butyrate','Ethyl acetate','Butyraldehyde','Blank','TMT','Female Pee', 'Male Pee', 'Peanut butter','Ethil tiglate','Propanal', 'Pups Bedding'}; 
     valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15];    %%new order for valves plotting only
     ploting=[6 2 3 4 5 11 14 7 8 9 10 15];
     for gg=1:length( ploting)
     k_for_plot(gg)=find( valves_12==ploting(gg));
     end

    
    valves= valves_12; 
    blank=find(valves==6);    %important!!!
    

   
%% load data if using mice other than the groups
%% load data after runing first and stats
if length(mice)~=length(exp_mice);
if length(mice)~=length(cont_mice);   
for i=1:length(mice);

mice_table{1,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results before']); %create a struct- rows for conditions, coloumns for mice
mice_table{2,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results after']); 

end
end
end

if length(mice)==length(exp_mice);
if sum(mice==exp_mice)==length(exp_mice);
load Mice_table_experiment
end
end

if length(mice)==length(cont_mice);
if   sum(mice==cont_mice)==length(cont_mice);
load Mice_table_Control
end
end


if   sum(mice)~=sum(cont_mice) && sum(mice)~=sum(exp_mice) ;


for i=1:length(mice);

mice_table{1,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results before']); %create a struct- rows for conditions, coloumns for mice
mice_table{2,i} = load(['cell to cell- mouse' num2str(mice(i)) 'results after']); 

end

end


for i=1:length(mice);

perm_table{1,i} = load(['perutation results- mouse' num2str(mice(i)) 'results before']); %create a struct- rows for conditions, coloumns for mice
perm_table{2,i} = load(['perutation results- mouse' num2str(mice(i)) 'results after']); 
perm_table{3,i} = load(['perutation results- mouse' num2str(mice(i)) ' bef_vs_after']); 

end
% %% for assesing significance in the difference before/after 
%      all_mice_p_values=[]%zeros(2000,numOfOdors);
% for i=1:length(mice);
%      p_values_for_dif_mat=[];
%     group_before=mice_table{1,i}.mat_for_diffrence_test;
%      group_after=mice_table{2,i}.mat_for_diffrence_test;
%      for ii=1:size(group_before,1)
%          for kk=1:size(group_before,3)
%              [h,p]=ttest2(group_before(ii,:,kk),group_after(ii,:,kk));
%              p_values_for_dif_mat(ii,kk)=p;
%              
%          end
%      end
%      mice_table{1,i}.p_values_for_dif_mat= p_values_for_dif_mat;
%      mice_table{2,i}.p_values_for_dif_mat= p_values_for_dif_mat;
%      all_mice_p_values= [all_mice_p_values;p_values_for_dif_mat];
% end
%          thresh=0.05;
%          signif_dif_bef_af_all_mice=abs(all_mice_p_values)<thresh;
         
         %% for assesing significance in the difference before/after by permutations
     all_mice_p_values=[]%zeros(2000,numOfOdors);
     thresh=0.05;
     
for i=1:length(mice);
    
             p_values_for_dif_mat=perm_table{3,i}.perm_p_vals_bef_vs_aft;
             all_mice_p_values= [all_mice_p_values;p_values_for_dif_mat];
             mice_table{1,i}.p_values_for_dif_mat= perm_table{1,i}.perm_p_vals
             mice_table{2,i}.p_values_for_dif_mat=perm_table{2,i}.perm_p_vals
end
     
      
%      mice_table{1,i}.p_values_for_dif_mat= p_values_for_dif_mat;
%      mice_table{2,i}.p_values_for_dif_mat= p_values_for_dif_mat;
   

       
         signif_dif_bef_af_all_mice=abs(all_mice_p_values)<thresh;
         
         

         
     %%
     
  x_resp=49:63
  y_resp=ones(size(x_resp)).*(-0.45)


  x_scale=[111 111]
  y_scale=[0 1]

 
     for mouse=mice;
          
          b=find(mice==mouse);
          a{1,b}=load(['cell to cell- mouse' num2str(mouse) 'results before']);
          a{2,b}=load(['cell to cell- mouse' num2str(mouse) 'results after']);
          SmoothDff_before=(a{1,b}.SmoothDff);
          SmoothDff_after=(a{2,b}.SmoothDff);
          mean_traces_mat_before=(a{1,b}.mean_traces_mat);
          mean_traces_mat_after=(a{2,b}.mean_traces_mat);
          
%           Ps_for_dif= mice_table{1,b}.p_values_for_dif_mat;
%           signif_dif_bef_af=Ps_for_dif<thresh;

          Ps_for_dif= perm_table{3,b}.perm_p_vals_bef_vs_aft;
          signif_dif_bef_af=Ps_for_dif<thresh;



%           
          response_magnitude_before=mice_table{1,find(mice==mouse)}.response_magnitude_by_integral;
          response_magnitude_after=mice_table{2,find(mice==mouse)}.response_magnitude_by_integral;
     
       
          
          
          
          %%to bring mean traces mat in!
          
          
          
    numOfNeurons=(a{1,b}.numOfNeurons);      
  
    meta_desicion_before=(a{1,b}.meta_desicion);
    meta_desicion_after=(a{2,b}.meta_desicion);
    signif_dif_bef_af(abs(meta_desicion_before)+abs(meta_desicion_after)==0)=0;
    
    x=zeros(20,1);
    x(20)=1;
    x(10)=1;
    y1=10;
    y2=20;
    y3=100;
   % size(SmoothDff);
    % should be :   neurons     trials     oders    frames
    NeuronsForPlot=8;
    for m=(0:NeuronsForPlot:numOfNeurons-NeuronsForPlot); %m is a vector defined in jumps of "Neurons for plot"
        figure
%         title(sprintf('Cells=%d-%d',m,6+m));
        kk=-1;
        x=ones(20,1);
%     
        for ii=1+m:NeuronsForPlot+m; 
            kk=kk+1;
            dummy=0;
            for k= k_for_plot; %k for odor
               dummy=dummy+1;
%            subaxis(NeuronsForPlot,9, k+(9*(kk)),'SpacingHoriz', 0,'SpacingVert',0) 

             h= subplot_tight(NeuronsForPlot,numOfOdors, dummy+(numOfOdors*(kk)),[0.01,0.01]); % with an exterior function
%              h= subplot(NeuronsForPlot,numOfOdors,
%              dummy+(numOfOdors*(kk))); if function does not exist
             
                plot(squeeze(SmoothDff_before(ii,:,k,:))','Color', [.5 .5 1]);
%                
                hold on
                plot(squeeze(SmoothDff_after(ii,:,k,:))','Color', [1 .5 .5]);
                plot((squeeze(mean_traces_mat_before(ii,k,:))'),'LineWidth',2,'Color','b');
                plot((squeeze(mean_traces_mat_after(ii,k,:))'),'LineWidth',2,'Color','r');
         
               
                
               

%                 axis off 
                hold on
            
                
                
               % ylim([-0.3 2])
               ylim([-0.5 1.5])            %%this lim was used for all paper plots
          %      ylim([-0.25 1.25])       %this is for the cno cont, where data is presented bigger
                xlim([0 lengthOfOdor])
%                 if astrix_desicion_before(ii,1,k)==1;
%                     hold on
%                   %  plot(y,x(20),'* r')
%                     
%                 end
%                   if astrix_desicion_after(ii,1,k)==1;
%                     hold on
%                   %  plot(y,x(10),'* b')
%                     
%                   end

                   if meta_desicion_before(ii,k)~=0;
                    hold on
                    if  meta_desicion_before(ii,k)<0
                       plot(y1,x(20),'o b','markerSize',12) 
                    else
                    plot(y1,x(20),'* b','markerSize',12)
                    end
                    end
                
                  if meta_desicion_after(ii,k)~=0;
                    hold on
                     if  meta_desicion_after(ii,k)<0
                          plot(y2,x(20),'o r','markerSize',12)
                     else
                    plot(y2,x(20),'* r','markerSize',12)
                     end
                  end
if signif_dif_bef_af(ii,k)==1;
   
    hold on
                    plot(y3,x(20),'* k','markerSize',12)
end
                 set(h, 'Visible', 'off')
            end
        end
%          set(gcf,'color','w'); 
%          tight_subplot(6,9,[.1 .1],[0.1 .1],[.1 .1]);
        
    end
     %plot the rest of the neurons: 4.6.2018
last_NeuronsForPlot=mod (numOfNeurons,NeuronsForPlot);
figure

 kk=-1;
for ii=(numOfNeurons-last_NeuronsForPlot+1:numOfNeurons)      
   kk=kk+1;
            dummy=0;
            for k= k_for_plot; %k for odor
               dummy=dummy+1;

             h= subplot(last_NeuronsForPlot,numOfOdors, dummy+(numOfOdors*(kk)));

                plot(squeeze(SmoothDff_before(ii,:,k,:))','Color', [.5 .5 1]);
%                
                hold on
                plot(squeeze(SmoothDff_after(ii,:,k,:))','Color', [1 .5 .5]);
                plot((squeeze(mean_traces_mat_before(ii,k,:))'),'LineWidth',2,'Color','b');
                plot((squeeze(mean_traces_mat_after(ii,k,:))'),'LineWidth',2,'Color','r');
%         


%                 axis off 
                hold on
            
                
                
               % ylim([-0.3 2])
                ylim([-0.5 1.5])
                xlim([0 lengthOfOdor])
%                 if astrix_desicion_before(ii,1,k)==1;
%                     hold on
%                   %  plot(y,x(20),'* r')
%                     
%                 end
%                   if astrix_desicion_after(ii,1,k)==1;
%                     hold on
%                   %  plot(y,x(10),'* b')
%                     
%                   end


    
%
%
%


                if meta_desicion_before(ii,k)~=0;
                    hold on
                      if  meta_desicion_before(ii,k)<0
                       plot(y1,x(20),'o b') 
                    else
                    plot(y1,x(20),'* b')
                      end          
                end
                
                  if meta_desicion_after(ii,k)~=0;
                    hold on
                    if  meta_desicion_after(ii,k)<0
                       plot(y2,x(20),'o r') 
                    else
                    plot(y2,x(20),'* r')
                    end
                    
                  end
                  
                  
if signif_dif_bef_af(ii,k)==1;
   
    hold on
                    plot(y3,x(20),'* k')
end
                 set(h, 'Visible', 'off')
            end
end
figure
title(['So far mouse number ' num2str(mouse) 'which had' num2str(numOfNeurons) 'cells']);
%plot the rest of the neurons: 4.6.2018

     end
 
%      

        
%      bef_af_SmoothDff{1,1}= load(['cell to cell-results before']);
%      bef_af_SmoothDff{1,2}= load(['cell to cell-results after']);
%      
%      
 
    %bef_af_SmoothDff{1,1}= load(['cell to cell-results before']);
    
    %after
    %save(['cell to cell-results after'], 'SmoothDff')
   % bef_af_SmoothDff{1,2}= load(['cell to cell-results after']);