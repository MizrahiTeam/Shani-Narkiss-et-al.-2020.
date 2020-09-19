%%PCA WITH Omri
clc
%close all
clear all
%load('Mice_table_Control.mat')
%load('Mice_table_Experiment.mat')

load('Mice_table_awake.mat')

%%

% 
% %%"responses is a mat that is organized for the pca procedure size :
% %%n_cellsX(n_odors*n_trials)
% %for 1 mouse:
% responses = [] ;
% num_of_odors=12;
% for i = 1:num_of_odors;
%     responses = [responses mat_for_diffrence_test(:,:,i)] ; %here a concatination
% end
% trial_type = sort(repmat(1:size(num_of_odors),1,5));
% % responses(59:63,:) = [] ;

%%

valves_12 = [2 3 4 5 6 7 8 9 10 11 14 15] ;
%  relevant_odors=[valves_12];
% % relevant_odors=[2 3 4 5 6 7]

 relevant_odors=[2 4 6 8 10 14];
% relevant_odors=[3 5 7 9 11 15];

ntrials=5;
odor_positions=[];
for i=1:length(relevant_odors);

odor_positions(i)=find(relevant_odors(i)==valves_12);
end


%formultiple mice:
responses = [] ;
num_of_odors=length(relevant_odors); %12;
%%
for j=1:length(mice_table);
responses_before = [] ;
responses_after = [] ;
    A = [] ;
    B = [] ;
    for i =odor_positions; %1:num_of_odors    this is for all odors 
        temp_bef=mice_table{1,j}.mat_for_diffrence_test(:,:,i);
        temp_af=mice_table{2,j}.mat_for_diffrence_test(:,:,i);
        A = [A temp_bef] ; %here a concatination
        B = [B temp_af] ; %here a concatination
%        
    end
    responses_before = [responses_before ; A] ;
    responses_after = [responses_after ; B] ;

trial_type = sort(repmat(1:num_of_odors,1,ntrials));



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%determine if before or after!!!!!!!!!!!!:
responses=responses_before;
responses=responses_after;
[U,S,V] = svd(responses,'econ') ;   %check if the responses vector should be transposed for sanity check-it shouldn't
                                    % U are the Eigenvectors orderd in
                                    % colomns S is the diag matrix of the
                                    % eigen values. V are the trials,
                                    % spanned at the eigen vectors space.
                                    % namely,
                                    % the contributionss of each EV to each
                                    % observation in the original responses
                                    % matrix
% figure
% plot(diag(S))

%according to this figure, one can chose the number of relevant components
% 1 conponent
% component = 1; %the chosen component
% figure
% plot(U(:,component)) %this figure showes how much each of the cells contribute to the first component
% scatter(trial_type,V(:,component))%this is the seperation according to the chosen component alone-one can change the component at will
%%
% % 2 components
% figure
% hold on
% component1 = 1 ;
% component2 = 2 ;
% color_vec={'blue','red', 'green', 'cyan','magenta','black' ,'yellow'};
% %the loop now runs on each one of the oddors and assigns different colors
% %for each odor
% for i = 1:num_of_odors
%     Color = color_vec{i} ;
%     scatter(V(trial_type==i,component1),V(trial_type==i,component2),[],Color,'filled')
% end
% a=gca;
% a.FontSize=20;
% xlabel('component 1')
% ylabel('component 2')
% a.FontSize=20;
% title('PCA BEFORE CNO')
% %xlim([-0.4 0])
% %ylim([-0.4 0.4])
% %title('PCA AFTER CNO')
%% 3 components only for three dimensional plot when 3 components are taken into account
figure
hold on
component1 = 1 ;
component2 = 2 ;
component3 = 3 ;
color_vec={'blue','red', 'green', 'cyan','magenta','black','yellow'};
centers=zeros(length(relevant_odors),3);
surround=zeros(length(relevant_odors),3);
noise_radius=zeros(length(relevant_odors),1);
noise_j=zeros(1,ntrials)

 axis equal
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5 0.5])

for i = 1:min([length(color_vec) length(relevant_odors)])
   
    scatter3(V(trial_type==i,component1),V(trial_type==i,component2),V(trial_type==i,component3),[],char(color_vec(i)), 'filled' )
    centers(i,:)=[mean(V(trial_type==i,component1)) mean(V(trial_type==i,component2)) mean(V(trial_type==i,component3)) ];
   %surround(i,:)=[std(V(trial_type==i,component1)) std(V(trial_type==i,component2)) std(V(trial_type==i,component3)) ]; %calculated for each axis seperately
   sep_dots=[V(trial_type==i,component1) V(trial_type==i,component2) V(trial_type==i,component3)]
    
   for j=1:ntrials
   noise_j(j)= sqrt(sum((sep_dots(j,:)-centers(i,:)).^2));
   end
   noise_radius(i)= mean(noise_j);
   
 surround(i,:)=[ noise_radius(i)]; %calculated for each axis seperately


    [x,y,z]=sphere;
   if i==1 
surf(x*surround(i,1)+centers(i,1),y*surround(i,2)+centers(i,2),z*surround(i,3)+centers(i,3),'FaceColor', char(color_vec(i)),'FaceAlpha',0.2); 
   end
   
 if i==2
surf(x*surround(i,1)+centers(i,1),y*surround(i,2)+centers(i,2),z*surround(i,3)+centers(i,3),'FaceColor', char(color_vec(i)),'FaceAlpha',0.2); 
   end   

    
    
end

 if length(relevant_odors)>length(color_vec)
for i = length(color_vec)+1:length(relevant_odors)
   
    scatter3(V(trial_type==i,component1),V(trial_type==i,component2),V(trial_type==i,component3),[],char(color_vec(i-length(color_vec))) )
    centers(i,:)=[mean(V(trial_type==i,component1)) mean(V(trial_type==i,component2)) mean(V(trial_type==i,component3)) ];
    surround(i,:)=[std(V(trial_type==i,component1)) std(V(trial_type==i,component2)) std(V(trial_type==i,component3)) ];
    
    surf(x*surround(i,1)+centers(i,1),y*surround(i,2)+centers(i,2),z*surround(i,3)+centers(i,3),'FaceColor', char(color_vec(i-length(color_vec))),'FaceAlpha',0.2);

end
 end

a=gca;
a.FontSize=20;



end
    