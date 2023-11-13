close all; clear all; clc;
% addpath("Functions")

% % Data folder name. 
prompt = {'Please enter the address to the folder containing the data','Please enter the name of your data set:','Please enter the intensity threshold value:','Please enter the number of z stacks (enter 0 if each folder is representing single organoid):'};
answer = inputdlg(prompt);
address=string(answer(1));
data_name = string(answer(2));
thresh_val=str2double(answer(3));
zStackAmount = str2double(answer(4));
ref_folder = fullfile(address,data_name);
cont_folder = dir(ref_folder);
pr={'would you like to process channel 2 data? (y/n)','Do you want to apply filter to the data?(y/n)'};
ans=inputdlg(pr);
ch2_pr= string((ans(1)));
filt_active=string((ans(2)));

promp={'would you like to save the processed data in ome.tiff format? (y/n)'};
ans_o=inputdlg(promp);
ome_tif_activate= string((ans_o(1)));

output_folder = fullfile(address, data_name+ "_Analyzed");
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end


%% find the adsdresses of subfolders


% FLIM-ISS_Time-decay
% cont_folder=fullfile(current_folder,data_folder);
list_of_fold=(findEndSubfolders(ref_folder))';

for i=1:numel(list_of_fold)
folderAddresses(i) = erase(list_of_fold{i}, ref_folder);
end

folderAddresses=folderAddresses';

%% Thresholding 
G_sum = [];
S_sum = [];
cond = [];
plate=[];
name = [];
figures = [];
figuresg = [];
G_av_ch1=[];
S_av_ch1=[];
cum_G=[];
cum_S=[];
cum_int=[];
numel(list_of_fold)
head_folder={};
cond_folder={};
treatment=[];
replicate_group=[];
Sample_ID=[];


%% 
for i = 1:numel(list_of_fold)  % Looping through condition. 

    
       
slash_ind=strfind(list_of_fold{i},'\'); 
char_add=char(list_of_fold{i});
char_short_add=char(folderAddresses(i));
slash_ind_short_add=strfind(folderAddresses(i),'\');
if (strcmp(char_add(slash_ind(end)+1:end),'MetaData') && length(dir(fullfile(char_add(1:slash_ind(end)-1),"*.tif")))==0)
% if (exist(char_add(1:slash_ind(end)-1), 'dir') == 7)
  
i
   copyfile(list_of_fold{i}, fullfile((output_folder),char_short_add)); 
  figures_folder = fullfile(output_folder, folderAddresses{i}); 
  add_ind=strfind(char(figures_folder),'\');
   char_list_of_fold=char(figures_folder);
   if (i==1)
   head_folder{i}=char_list_of_fold(1:add_ind(end-1)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end-1)-1),head_folder{i-1})) 
     
    head_folder{i}=head_folder{i-1};
    
 else
    head_folder{i}=char_list_of_fold(1:add_ind(end-1)-1); 
     
 end
   end 
    
      if (i==1)
   cond_fold{i}=char_list_of_fold(1:add_ind(end-2)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end-2)-1),cond_fold{i-1})) 
     
    cond_fold{i}=cond_fold{i-1};
    
 else
    cond_fold{i}=char_list_of_fold(1:add_ind(end-2)-1); 
     
 end
      end 
   
   cond_slash_ind=strfind(cond_fold,'\'); 
   treatment_slash_ind=strfind(head_folder,'\'); 
   replicate_group_temp=cond_fold{i}(cond_slash_ind{i}(end)+1:end);  
   treatment_temp=head_folder{i}(treatment_slash_ind{i}(end)+1:end);  
   
   
   continue;
% else
    
%  fullfile(output_fold,[char(folderAddresse(jj))   
% end    
    
end  

if (strcmp(char_add(slash_ind(end)+1:end),'MetaData') && length(dir(fullfile(char_add(1:slash_ind(end)-1),"*.tif")))~=0)

i

% mkdir (fullfile(output_folder,char(folderAddresses(i))))

          figures_folder = fullfile(output_folder, folderAddresses{i});
        if ~exist(figures_folder,'dir')
          mkdir(figures_folder) 
        end
% copyfile(list_of_fold{i}, fullfile((output_folder),char_short_add));

  if (zStackAmount==0)
     
 disp(["Analyzing contents of " + string(char_add(1:slash_ind(end)-1))]);
   ref_files = dir(fullfile(char_add(1:slash_ind(end)-1),"*.tif"));  % reads the .tif file in the condition folder.  
   files=sort_nat({ref_files.name},'ascend');
   
   temp_int = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{1}),1));
   
   cur_G_Average = zeros(size(temp_int,1),size(temp_int,1));
   cur_Int_Average = zeros(size(temp_int,1),size(temp_int,1));  
   cur_S_Average=zeros(size(temp_int,1),size(temp_int,1));
 
   
    if (rem(numel(ref_files),zStackAmount)~= 0 && (zStackAmount~=0)) 
       remainders = zStackAmount - rem(numel(ref_files),zStackAmount);
       disp("You are missing a file or inputted the wrong number of z stacks. " + remainders + " files are missing." );
    end  
    

        add_ind=strfind(char(figures_folder),'\');
   char_list_of_fold=char(figures_folder);
   if (i==1)
   head_folder{i}=char_list_of_fold(1:add_ind(end-1)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end-1)-1),head_folder{i-1})) 
     
    head_folder{i}=head_folder{i-1};
    
 else
    head_folder{i}=char_list_of_fold(1:add_ind(end-1)-1); 
     
 end
   end 
    
      if (i==1)
   cond_fold{i}=char_list_of_fold(1:add_ind(end-2)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end-2)-1),cond_fold{i-1})) 
     
    cond_fold{i}=cond_fold{i-1};
    
 else
    cond_fold{i}=char_list_of_fold(1:add_ind(end-2)-1); 
     
 end
   end 
   
   cond_slash_ind=strfind(cond_fold,'\'); 
   treatment_slash_ind=strfind(head_folder,'\'); 
   treatment_temp=cond_fold{i}(cond_slash_ind{i}(end)+1:end);  
   replicate_group_temp=head_folder{i}(treatment_slash_ind{i}(end-2)+1:treatment_slash_ind{i}(end-1)-1);   
   
    
 
          
info = imfinfo(fullfile(char_add(1:slash_ind(end)-1),ref_files(1).name));
tif_length=length(info);
if ((tif_length==3) || (tif_length==6 && ch2_pr=='n'))

%        for z = 1:zStackAmount   %% Groups all the z stacks together. 
       
for j = 1:numel(ref_files)  % looping through z slices.
       
      G = []; S = []; int = [];


         cur_int = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),1));
         cur_G = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),2));
         cur_S = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),3));
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
           
 
            
       
            figures(j) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch1_Image_",files{j})]));
            saveas(figures(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch1_Image_",files{j}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(j) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch1_",files{j}(1:end-4),'.png')]));     
            
            figuresg(j) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("S_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("S_ch1_",files{j}(1:end-4),'.png')])); 
            
            
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       % z_struct = medfiltPhasor(z_struct,3);
       z_struct = wavefiltPhasor(z_struct); 

       end   
     figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch1_",files{j})]));
    saveas(gcf, fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch1_",files{j}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch1_Image_",files{j})]));
    image2=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch1_",files{j})]));
    image3=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch1_",files{j})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),['ch1_singz_',files{j}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);
if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),['ch1_singz_',files{j}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
end
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch1 = [G_av_ch1;GZ];
       S_av_ch1 = [S_av_ch1;SZ];
%        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');
       treatment = [treatment;string(treatment_temp)]            
       replicate_group=[replicate_group;string(replicate_group_temp)]
       name = [name;string(files{j})];  
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        end
       und_ind=strfind(files{j},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),string(files{j}(1:end-4)) + '_ch1.xlsx');   % Name of the excel file. 
       writetable(DataTable,filefolder,'Sheet',1)
     
            G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];
       replicate_group=[];
       Sample_ID=[];
       treatment=[];
end
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];
       replicate_group=[];
       Sample_ID=[];
       treatment=[];
%        z_struct=
%        [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
%        G_av_ch1 = [G_av_ch1;GZ];
%        S_av_ch1 = [S_av_ch1;SZ];
% %        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');
%        cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
%        plate=[plate;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        name = [name;string(files{j})];  
%     
% %        end
%        und_ind=strfind(files{j},'_');
%        DataTable=table(cond,plate,name,G_av_ch1, S_av_ch1);
% %        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
%        filefolder = fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),string(files{j}(1:end-4)) + '_ch1.xlsx');   % Name of the excel file. 
%        writetable(DataTable,filefolder,'Sheet',1)
%        G_av_ch1=[];
%        S_av_ch1=[];
%        cond = [];
%        name = [];
%        plate=[];
%        clear cond name G_av_ch1 S_av_ch1
       current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       % current_struct = medfiltPhasor(current_struct,3);
       current_struct = wavefiltPhasor(current_struct); 

       end
       plotPhasorFast(current_struct);
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S. 
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
       treatment = [treatment;string(treatment_temp)];
       name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
       replicate_group=[replicate_group;string(replicate_group_temp)]
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);
       % for ISS
       filefolder = fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_average_values_ch1.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_cumulative_phasor_ch1.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_cumulative_phasor_ch1.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all
%%
% 
%  map_res=512;
%  x=linspace(-max(cum_struct.G(:)),max(cum_struct.G(:)),512);
%     ctrs{1} = x;
%     ctrs{2} = x;
%     temp_G=cum_struct.G(cum_struct.int>0)-1.526e-05;
%     temp_S=cum_struct.S(cum_struct.int>0)-1.526e-05;
% %     temp_int=cum_int';
% [N,~] = hist3([temp_S(:),temp_G(:)],ctrs);
% figure()
% imagesc(x,x,flip(N));
% colormap(gca,jet); ax = gca; ax.Colormap(1,:)= [1 1 1]; caxis('auto');
% % colorbar; axis xy;
% 
% colorbar; axis image;
% 
% x_circle     = [0:1/512:1];
% y_circle_pos = sqrt(0.25-((x_circle-0.5).^2));
% y_circle_neg = -sqrt(0.25-((x_circle-0.5).^2));;
% hold on; plot(x_circle,[y_circle_pos;y_circle_neg],'k','LineWidth',1)
% set(gca,'ydir','reverse')
% axis([0 1 -0.625 0])
% 
% % axis xy;
% % set(gca,'XTick',[0,0.5,1], 'YTick', [0,0.25,0.5])
% % xticks([1/2:1/2^4:1]);
% % xticklabels({'0','', '','','0.5','', '','','1'});
% % 
% % yticks([0:1/2^4:1/2]);
% % yticklabels({'0','', '0.25','','0.5'});
% % xlabel('G');ylabel('S')

 %%      
       
end       
G_sum = [];
S_sum = [];
cond = [];
name = [];
figures = [];
figuresg = [];
G_av_ch1=[];
S_av_ch1=[];
cum_G=[];
cum_S=[];
int=[];
cum_int=[];
S=[];
G=[];
plate=[];
treatment=[];
replicate_group=[];
Sample_ID=[];

if ((tif_length==6 && ch2_pr=='y')) 
    
    for j = 1:numel(ref_files)  % looping through z slices.
       
      G = []; S = []; int = [];

     %   for z = 1:zStackAmount   %% Groups all the z stacks together. 
         cur_int = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),1));
         cur_G = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),2));
         cur_S = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),3));
                     
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
       
                 
            figures(j) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch1_Image_",files{j})]));
            saveas(figures(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch1_Image_",files{j}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(j) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch1_",files{j}(1:end-4),'.png')]));     
            
            figuresg(j) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("S_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("S_ch1_",files{j}(1:end-4),'.png')])); 
                       
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%        z_struct = medfiltPhasor(z_struct,3);
       z_struct = wavefiltPhasor(z_struct); 
       end   
         figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch1_",files{j})]));
    saveas(gcf, fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch1_",files{j}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch1_Image_",files{j})]));
    image2=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch1_",files{j})]));
    image3=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch1_",files{j})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),['ch1_singz_',files{j}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);
if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),['ch1_singz_',files{j}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
end
           
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch1 = [G_av_ch1;GZ];
       S_av_ch1 = [S_av_ch1;SZ];
%        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');
%        cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
%        plate=[plate;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        name = [name;string(files{j})];  

       treatment = [treatment;string(treatment_temp)]            
       replicate_group=[replicate_group;string(replicate_group_temp)]
       name = [name;string(files{j})];  
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        end
       und_ind=strfind(files{j},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);
       
%        end
%        und_ind=strfind(files{j},'_');
%        DataTable=table(cond,plate,name,G_av_ch1, S_av_ch1);
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),string(files{j}(1:end-4)) + '_ch1.xlsx');   % Name of the excel file. 
       writetable(DataTable,filefolder,'Sheet',1)
     
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
              
    end 
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];    
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
%        clear cond name G_av_ch1 S_av_ch1
       
        current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       % current_struct = medfiltPhasor(current_struct,3);
       current_struct = wavefiltPhasor(current_struct); 
       end
       plotPhasorFast(current_struct);
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S. 
%        plate=[plate;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
%        name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
%        DataTable=table(cond,name,Organoid_ave_G, Organoid_ave_S);
       
       
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%      cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
       treatment = [treatment;string(treatment_temp)];
       name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
       replicate_group=[replicate_group;string(replicate_group_temp)]
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);
       
       
       % for ISS
       filefolder = fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_average_values_ch1.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_cumulative_phasor_ch1.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_cumulative_phasor_ch1.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all
       
       
       
       
  
%Ch2 data 
G_sum = [];
S_sum = [];
cond = [];
name = [];
figures = [];
figuresg = [];
G_av_ch1=[];
S_av_ch1=[];
cum_G=[];
cum_S=[];
cum_int=[];
S=[];
G=[];
int=[];
G_av_ch2=[];
S_av_ch2=[];
plate=[];
treatment=[];
replicate_group=[];
Sample_ID=[];
clear cur_int cur_G cur_S
       
      for j = 1:numel(ref_files)
      %  for z = 1:zStackAmount   %% Groups all the z stacks together. 
           
         cur_int = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),4));
         cur_G = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),5));
         cur_S = double(imread(fullfile(char_add(1:slash_ind(end)-1),files{j}),6));
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
                  
            
            figures(j) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch2_Image_",files{j})]));
            saveas(figures(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch2_Image_",files{j}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(j) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch2_",files{j})]));
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch2_",files{j}(1:end-4),'.png')]));     
            
            figuresg(j) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("S_ch2_",files{j})]));
            saveas(figuresg(j), fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("S_ch2_",files{j}(1:end-4),'.png')]));         
            
            
       close all
       
  z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       % z_struct = medfiltPhasor(z_struct,3);
       z_struct = wavefiltPhasor(z_struct); 
       end   
         figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch2_",files{j})]));
    saveas(gcf, fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch2_",files{j}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("Int_ch2_Image_",files{j})]));
    image2=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("G_ch2_",files{j})]));
    image3=imread(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),[strcat("single_z_phasor_ch2_",files{j})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),['ch2_singz_',files{j}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);
if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),['ch2_singz_',files{j}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
end
   
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch2 = [G_av_ch2;GZ];
       S_av_ch2 = [S_av_ch2;SZ];
%        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');
%        cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
%        plate=[plate;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        name = [name;string(files{j})];  
    
       treatment = [treatment;string(treatment_temp)]            
       replicate_group=[replicate_group;string(replicate_group_temp)]
       name = [name;string(files{j})];  
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        end
%        und_ind=strfind(files{j},'_');
%        DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);   
       
       
       
%        end
       und_ind=strfind(files{j},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch2, S_av_ch2);
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),string(files{j}(1:end-4)) + '_ch2.xlsx');   % Name of the excel file. 
       writetable(DataTable,filefolder,'Sheet',1)
     
       G_av_ch2=[];
       S_av_ch2=[];
       cond = [];
       name = [];
       plate=[];
       treatment=[];
replicate_group=[];
Sample_ID=[];
       
      end
       G_av_ch2=[];
       S_av_ch2=[];
       cond = [];
       name = [];
       plate=[]; 
       treatment=[];
replicate_group=[];
Sample_ID=[];
%        clear cond name G_av_ch1 S_av_ch1
       
        current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%        current_struct = medfiltPhasor(current_struct,3);
       current_struct = wavefiltPhasor(current_struct); 
       end
       plotPhasorFast(current_struct);
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S. 
%        plate=[plate;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
%        name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
%        DataTable=table(cond,name,Organoid_ave_G, Organoid_ave_S);
       
       
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%      cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
       treatment = [treatment;string(treatment_temp)];
       name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
       replicate_group=[replicate_group;string(replicate_group_temp)]
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);
       
       
       % for ISS
       filefolder = fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_average_values_ch2.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_cumulative_phasor_ch2.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_cumulative_phasor_ch2.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all







end




S=[];
G=[];
int=[];
plate=[];
treatment=[];
replicate_group=[];
Sample_ID=[];
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       G=[];
       S=[];
       Organoid_ave_G=[];
       Organoid_ave_S=[];
% clear cond name G_av_ch1 S_av_ch1 name cond G S Organoid_ave_G Organoid_ave_S









   


 end

end

    







if (~strcmp(char_add(slash_ind(end)+1:end),'MetaData') && (zStackAmount==0))


    
    
    
   i

% mkdir (fullfile(output_folder,char(folderAddresses(i))))

          figures_folder = fullfile(output_folder, folderAddresses{i});
        if ~exist(figures_folder,'dir')
          mkdir(figures_folder) 
        end

     
 disp(["Analyzing contents of " + string(list_of_fold{i})]);
   ref_files = dir(fullfile(list_of_fold{i},"*.tif"));  % reads the .tif file in the condition folder.  
   files=sort_nat({ref_files.name},'ascend');
   
   temp_int = double(imread(fullfile(list_of_fold{i},files{1}),1));
   
   cur_G_Average = zeros(size(temp_int,1),size(temp_int,1));
   cur_Int_Average = zeros(size(temp_int,1),size(temp_int,1));  
   cur_S_Average=zeros(size(temp_int,1),size(temp_int,1));
 
   
    if (rem(numel(ref_files),zStackAmount)~= 0 && (zStackAmount~=0)) 
       remainders = zStackAmount - rem(numel(ref_files),zStackAmount);
       disp("You are missing a file or inputted the wrong number of z stacks. " + remainders + " files are missing." );
    end  
    

        add_ind=strfind(char(figures_folder),'\');
   char_list_of_fold=char(figures_folder);
   if (i==1)
   head_folder{i}=char_list_of_fold(1:add_ind(end)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end-1)-1),head_folder{i-1})) 
     
    head_folder{i}=head_folder{i-1};
    
 else
    head_folder{i}=char_list_of_fold(1:add_ind(end)-1); 
     
 end
   end 
    
      if (i==1)
   cond_fold{i}=char_list_of_fold(1:add_ind(end-1)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end-1)-1),cond_fold{i-1})) 
     
    cond_fold{i}=cond_fold{i-1};
    
 else
    cond_fold{i}=char_list_of_fold(1:add_ind(end-1)-1); 
     
 end
   end 
   
   
   cond_slash_ind=strfind(cond_fold,'\'); 
   treatment_slash_ind=strfind(head_folder,'\'); 
   treatment_temp=cond_fold{i}(cond_slash_ind{i}(end)+1:end);  
   replicate_group_temp=head_folder{i}(treatment_slash_ind{i}(end-2)+1:treatment_slash_ind{i}(end-1)-1);   
    
 
          
info = imfinfo(fullfile(list_of_fold{i},ref_files(1).name));
tif_length=length(info);
if ((tif_length==3) || (tif_length==6 && ch2_pr=='n'))

%        for z = 1:zStackAmount   %% Groups all the z stacks together. 
       
for j = 1:numel(ref_files)  % looping through z slices.
       
      G = []; S = []; int = [];


         cur_int = double(imread(fullfile(list_of_fold{i},files{j}),1));
         cur_G = double(imread(fullfile(list_of_fold{i},files{j}),2));
         cur_S = double(imread(fullfile(list_of_fold{i},files{j}),3));
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
           
 
            
       
            figures(j) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(j), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{j})]));
            saveas(figures(j), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{j}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(j) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(figures_folder,[strcat("G_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(figures_folder,[strcat("G_ch1_",files{j}(1:end-4),'.png')]));     
            
            figuresg(j) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(figures_folder,[strcat("S_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(figures_folder,[strcat("S_ch1_",files{j}(1:end-4),'.png')]));            
            
            
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%        z_struct = medfiltPhasor(z_struct,3);
        z_struct = wavefiltPhasor(z_struct); 
       end   
     figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{j})]));
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{j}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(figures_folder,[strcat("Int_ch1_Image_",files{j})]));
    image2=imread(fullfile(figures_folder,[strcat("G_ch1_",files{j})]));
    image3=imread(fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{j})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(figures_folder,['ch1_singz_',files{j}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);
if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(figures_folder,['ch1_singz_',files{j}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
end
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch1 = [G_av_ch1;GZ];
       S_av_ch1 = [S_av_ch1;SZ];
%        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');
       treatment = [treatment;string(treatment_temp)];            
       replicate_group=[replicate_group;string(replicate_group_temp)]
       name = [name;string(files{j})];  
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
       
       und_ind=strfind(files{j},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(figures_folder,string(files{j}(1:end-4)) + '_ch1.xlsx');   % Name of the excel file. 
       writetable(DataTable,filefolder,'Sheet',1)
     
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];
     treatment=[];
replicate_group=[];
Sample_ID=[];
end
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];
       treatment=[];
replicate_group=[];
Sample_ID=[];
%        z_struct=
%        [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
%        G_av_ch1 = [G_av_ch1;GZ];
%        S_av_ch1 = [S_av_ch1;SZ];
% %        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');
%        cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
%        plate=[plate;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%        name = [name;string(files{j})];  
%     
% %        end
%        und_ind=strfind(files{j},'_');
%        DataTable=table(cond,plate,name,G_av_ch1, S_av_ch1);
% %        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
%        filefolder = fullfile(output_folder,char_short_add(1:slash_ind_short_add(end)-1),string(files{j}(1:end-4)) + '_ch1.xlsx');   % Name of the excel file. 
%        writetable(DataTable,filefolder,'Sheet',1)
%        G_av_ch1=[];
%        S_av_ch1=[];
%        cond = [];
%        name = [];
%        plate=[];
%        clear cond name G_av_ch1 S_av_ch1
       current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%        current_struct = medfiltPhasor(current_struct,3);
       current_struct = wavefiltPhasor(current_struct); 
       end
       plotPhasorFast(current_struct);
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S.       
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%      cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
       treatment = [treatment;string(treatment_temp)];
       name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
       replicate_group=[replicate_group;string(replicate_group_temp)]
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);
       
       
       % for ISS
       filefolder = fullfile((head_folder{i}),[char_list_of_fold(add_ind(end)+1:end),'_average_values_ch1.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end)+1:end),'_cumulative_phasor_ch1.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end)+1:end),'_cumulative_phasor_ch1.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all
%%
% 
%  map_res=512;
%  x=linspace(-max(cum_struct.G(:)),max(cum_struct.G(:)),512);
%     ctrs{1} = x;
%     ctrs{2} = x;
%     temp_G=cum_struct.G(cum_struct.int>0)-1.526e-05;
%     temp_S=cum_struct.S(cum_struct.int>0)-1.526e-05;
% %     temp_int=cum_int';
% [N,~] = hist3([temp_S(:),temp_G(:)],ctrs);
% figure()
% imagesc(x,x,flip(N));
% colormap(gca,jet); ax = gca; ax.Colormap(1,:)= [1 1 1]; caxis('auto');
% % colorbar; axis xy;
% 
% colorbar; axis image;
% 
% x_circle     = [0:1/512:1];
% y_circle_pos = sqrt(0.25-((x_circle-0.5).^2));
% y_circle_neg = -sqrt(0.25-((x_circle-0.5).^2));;
% hold on; plot(x_circle,[y_circle_pos;y_circle_neg],'k','LineWidth',1)
% set(gca,'ydir','reverse')
% axis([0 1 -0.625 0])
% 
% % axis xy;
% % set(gca,'XTick',[0,0.5,1], 'YTick', [0,0.25,0.5])
% % xticks([1/2:1/2^4:1]);
% % xticklabels({'0','', '','','0.5','', '','','1'});
% % 
% % yticks([0:1/2^4:1/2]);
% % yticklabels({'0','', '0.25','','0.5'});
% % xlabel('G');ylabel('S')

 %%      
       
end       
G_sum = [];
S_sum = [];
cond = [];
name = [];
figures = [];
figuresg = [];
G_av_ch1=[];
S_av_ch1=[];
cum_G=[];
cum_S=[];
int=[];
cum_int=[];
S=[];
G=[];
plate=[];
treatment=[];
replicate_group=[];
Sample_ID=[];
if ((tif_length==6 && ch2_pr=='y')) 
    
    for j = 1:numel(ref_files)  % looping through z slices.
       
      G = []; S = []; int = [];

     %   for z = 1:zStackAmount   %% Groups all the z stacks together. 
         cur_int = double(imread(fullfile(list_of_fold{i},files{j}),1));
         cur_G = double(imread(fullfile(list_of_fold{i},files{j}),2));
         cur_S = double(imread(fullfile(list_of_fold{i},files{j}),3));
                     
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
       
                 
            figures(j) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(j), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{j})]));
            saveas(figures(j), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{j}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(j) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(figures_folder,[strcat("G_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(figures_folder,[strcat("G_ch1_",files{j}(1:end-4),'.png')]));     
            
            figuresg(j) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(figures_folder,[strcat("S_ch1_",files{j})]));
            saveas(figuresg(j), fullfile(figures_folder,[strcat("S_ch1_",files{j}(1:end-4),'.png')]));            
            
            
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%        z_struct = medfiltPhasor(z_struct,3);
       z_struct = wavefiltPhasor(z_struct); 

       end   
     figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{j})]));
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{j}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(figures_folder,[strcat("Int_ch1_Image_",files{j})]));
    image2=imread(fullfile(figures_folder,[strcat("G_ch1_",files{j})]));
    image3=imread(fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{j})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(figures_folder,['ch1_singz_',files{j}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);
if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(figures_folder,['ch1_singz_',files{j}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
end
           
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch1 = [G_av_ch1;GZ];
       S_av_ch1 = [S_av_ch1;SZ];
%        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');
  
       treatment = [treatment;string(treatment_temp)];            
       replicate_group=[replicate_group;string(replicate_group_temp)]
       name = [name;string(files{j})];  
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
       
       und_ind=strfind(files{j},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);   

%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(figures_folder,string(files{j}(1:end-4)) + '_ch1.xlsx');
       writetable(DataTable,filefolder,'Sheet',1)
     
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
              
    end 
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       plate=[];   
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
%        clear cond name G_av_ch1 S_av_ch1
       
        current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%      current_struct = medfiltPhasor(current_struct,3);
       current_struct = wavefiltPhasor(current_struct); 

       end
       plotPhasorFast(current_struct);
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S.     
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%      cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
       treatment = [treatment;string(treatment_temp)];
       name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
       replicate_group=[replicate_group;string(replicate_group_temp)]
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);       
       
       % for ISS
       filefolder = fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_average_values_ch1.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end)+1:end),'_cumulative_phasor_ch1.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end)+1:end),'_cumulative_phasor_ch1.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all          
            
  
%Ch2 data 
G_sum = [];
S_sum = [];
cond = [];
name = [];
figures = [];
figuresg = [];
G_av_ch1=[];
S_av_ch1=[];
cum_G=[];
cum_S=[];
cum_int=[];
S=[];
G=[];
int=[];
G_av_ch2=[];
S_av_ch2=[];
plate=[];
treatment=[];
replicate_group=[];
Sample_ID=[];
clear cur_int cur_G cur_S
       
      for j = 1:numel(ref_files)
      %  for z = 1:zStackAmount   %% Groups all the z stacks together. 
           
         cur_int = double(imread(fullfile(list_of_fold{i},files{j}),4));
         cur_G = double(imread(fullfile(list_of_fold{i},files{j}),5));
         cur_S = double(imread(fullfile(list_of_fold{i},files{j}),6));
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
                  
             figures(j) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(j), fullfile(figures_folder,[strcat("Int_ch2_Image_",files{j})]));
            saveas(figures(j), fullfile(figures_folder,[strcat("Int_ch2_Image_",files{j}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(j) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(figures_folder,[strcat("G_ch2_",files{j})]));
            saveas(figuresg(j), fullfile(figures_folder,[strcat("G_ch2_",files{j}(1:end-4),'.png')]));     
            
            figuresg(j) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(j), fullfile(figures_folder,[strcat("S_ch2_",files{j})]));
            saveas(figuresg(j), fullfile(figures_folder,[strcat("S_ch2_",files{j}(1:end-4),'.png')]));            
            
            
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%        z_struct = medfiltPhasor(z_struct,3);
       z_struct = wavefiltPhasor(z_struct); 

       end   
     figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch2_",files{j})]));
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch2_",files{j}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(figures_folder,[strcat("Int_ch2_Image_",files{j})]));
    image2=imread(fullfile(figures_folder,[strcat("G_ch2_",files{j})]));
    image3=imread(fullfile(figures_folder,[strcat("single_z_phasor_ch2_",files{j})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(figures_folder,['ch2_singz_',files{j}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);
if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(figures_folder,['ch2_singz_',files{j}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
end
      
   
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch2 = [G_av_ch2;GZ];
       S_av_ch2 = [S_av_ch2;SZ];
%        cond_find=strfind(char_add(1:slash_ind(end)-1),'_');   
       treatment = [treatment;string(treatment_temp)];            
       replicate_group=[replicate_group;string(replicate_group_temp)]
       name = [name;string(files{j})];  
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];

       und_ind=strfind(files{j},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch2, S_av_ch2);
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(figures_folder,string(files{j}(1:end-4)) + '_ch2.xlsx');
       writetable(DataTable,filefolder,'Sheet',1)
     
       G_av_ch2=[];
       S_av_ch2=[];
       cond = [];
       name = [];
       plate=[];
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
      end
       G_av_ch2=[];
       S_av_ch2=[];
       cond = [];
       name = [];
       plate=[];   
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
%        clear cond name G_av_ch1 S_av_ch1
    current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
%      current_struct = medfiltPhasor(current_struct,3);
       current_struct = wavefiltPhasor(current_struct); 
       end
       plotPhasorFast(current_struct);
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S.       
       Sample_ID=[Sample_ID;string(char_list_of_fold(add_ind(end-2)+1:add_ind(end-1)-1))];
%      cond = [cond;string(char_list_of_fold(add_ind(end-3)+1:add_ind(end-2)-1))];
       treatment = [treatment;string(treatment_temp)];
       name=string([char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_ave_values']);
       replicate_group=[replicate_group;string(replicate_group_temp)]
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);       
       
       
       
       
       
       % for ISS
       filefolder = fullfile((head_folder{i}),[char_list_of_fold(add_ind(end-1)+1:add_ind(end)-1),'_average_values_ch2.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end)+1:end),'_cumulative_phasor_ch2.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char_list_of_fold(add_ind(end)+1:end),'_cumulative_phasor_ch2.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all          






end



treatment=[];
replicate_group=[];
Sample_ID=[];
S=[];
G=[];
int=[];
plate=[];
       G_av_ch1=[];
       S_av_ch1=[];
       G_av_ch2=[];
       S_av_ch2=[];
       cond = [];
       name = [];
       G=[];
       S=[];
       Organoid_ave_G=[];
       Organoid_ave_S=[];
% clear cond name G_av_ch1 S_av_ch1 name cond G S Organoid_ave_G Organoid_ave_S

 
    
    
    
    
    
    
    
    
    
    
    
    
    
    

end












if (~strcmp(char_add(slash_ind(end)+1:end),'MetaData') && (zStackAmount~=0))


 i
 disp(["Analyzing contents of " + string(list_of_fold{i})]);
   ref_files = dir(fullfile(list_of_fold{i},"*.tif"));  % reads the ref file in the condition folder.  
   files=sort_nat({ref_files.name},'ascend');
   
   temp_int = double(imread(fullfile(list_of_fold{i},files{1}),1) );
   
   cur_G_Average = zeros(size(temp_int,1),size(temp_int,1));
   cur_Int_Average = zeros(size(temp_int,1),size(temp_int,1));  
   cur_S_Average=zeros(size(temp_int,1),size(temp_int,1));
 
   
    if rem(numel(ref_files),zStackAmount) ~= 0
       remainders = zStackAmount - rem(numel(ref_files),zStackAmount);
       disp("You are missing a file or inputted the wrong number of z stacks. " + remainders + " files are missing." );
    end  
    
          figures_folder = fullfile(output_folder, folderAddresses{i});
        if ~exist(figures_folder,'dir')
          mkdir(figures_folder) 
        end
        add_ind=strfind(char(figures_folder),'\');
   char_list_of_fold=char(figures_folder);
   if (i==1)
   head_folder{i}=char_list_of_fold(1:add_ind(end)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end)-1),head_folder{i-1})) 
     
    head_folder{i}=head_folder{i-1};
    
 else
    head_folder{i}=char_list_of_fold(1:add_ind(end)-1); 
     
 end
   end 
    
      if (i==1)
   cond_fold{i}=char_list_of_fold(1:add_ind(end-1)-1);
   else
 if (strcmp(char_list_of_fold(1:add_ind(end-1)-1),cond_fold{i-1})) 
     
    cond_fold{i}=cond_fold{i-1};
    
 else
    cond_fold{i}=char_list_of_fold(1:add_ind(end-1)-1); 
     
 end
   end 
   
   cond_slash_ind=strfind(cond_fold,'\'); 
   treatment_slash_ind=strfind(head_folder,'\'); 
   replicate_group_temp=cond_fold{i}(cond_slash_ind{i}(end-1)+1:cond_slash_ind{i}(end)-1);  
   treatment_temp=head_folder{i}(treatment_slash_ind{i}(end-1)+1:treatment_slash_ind{i}(end)-1); 
   for j = 1:numel(ref_files)/zStackAmount  % looping through z slices. 
      G = []; S = []; int = [];
    
 
          
info = imfinfo(fullfile(list_of_fold{i},ref_files(1).name));
tif_length=length(info);
if ((tif_length==3) || (tif_length==6 && ch2_pr=='n'))

       for z = 1:zStackAmount   %% Groups all the z stacks together. 
           
         cur_int = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),1));
         cur_G = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),2));
         cur_S = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),3));
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
           
 
            
       
            figures(z) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(z), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{((j-1)*3)+z})]));
            saveas(figures(z), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{((j-1)*3)+z}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(z) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(z), fullfile(figures_folder,[strcat("G_ch1_",files{((j-1)*3)+z})]));
            saveas(figuresg(z), fullfile(figures_folder,[strcat("G_ch1_",files{((j-1)*3)+z}(1:end-4),'.png')]));     
            
            figuresg(z) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(z), fullfile(figures_folder,[strcat("S_ch1_",files{((j-1)*3)+z})]));
            saveas(figuresg(z), fullfile(figures_folder,[strcat("S_ch1_",files{((j-1)*3)+z}(1:end-4),'.png')])); 
            
            
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       z_struct = medfiltPhasor(z_struct,3);
       end   
     figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{((j-1)*3)+z})]));
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{((j-1)*3)+z}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(figures_folder,[strcat("Int_ch1_Image_",files{((j-1)*3)+z})]));
    image2=imread(fullfile(figures_folder,[strcat("G_ch1_",files{((j-1)*3)+z})]));
    image3=imread(fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{((j-1)*3)+z})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(figures_folder,['ch1_singz_',files{((j-1)*3)+z}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);

    if (ome_tif_activate=='y')
    
     multiplane_ome_name=char(fullfile(figures_folder,['ch1_singz_',files{((j-1)*3)+z}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);

    end
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch1 = [G_av_ch1;GZ];
       S_av_ch1 = [S_av_ch1;SZ];
       cond_find=strfind(char_list_of_fold,'_');
       Sample_ID = [Sample_ID;string(char_list_of_fold(add_ind(end-1)+1:cond_find(end-1)-1))];
       name = [name;string(files{((j-1)*3)+z})];  
       treatment=[treatment;string(treatment_temp)];
       replicate_group=[replicate_group;string(replicate_group_temp)];
       end
       und_ind=strfind(files{((j-1)*3)+1},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(figures_folder,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '_ch1.xlsx');   % Name of the excel file. 
       writetable(DataTable,filefolder,'Sheet',1)
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
%        clear cond name G_av_ch1 S_av_ch1
       current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       current_struct = medfiltPhasor(current_struct,3);
       end
       plotPhasorFast(current_struct);
       
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S. 
       Sample_ID =string(char_list_of_fold(add_ind(end-1)+1:cond_find(end-1)-1));
       name=string([char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_ave_values_',(char_list_of_fold(add_ind(end)+1:end))]);
       treatment=[treatment;string(treatment_temp)];
       replicate_group=[replicate_group;string(replicate_group_temp)];
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);
       % for ISS
       filefolder = fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_',(char_list_of_fold(add_ind(end)+1:end)),'_ch1.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_cumulative_phasor_',(char_list_of_fold(add_ind(end)+1:end)),'_ch1.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_cumulative_phasor_',(char_list_of_fold(add_ind(end)+1:end)),'_ch1.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all
%%
% 
%  map_res=512;
%  x=linspace(-max(cum_struct.G(:)),max(cum_struct.G(:)),512);
%     ctrs{1} = x;
%     ctrs{2} = x;
%     temp_G=cum_struct.G(cum_struct.int>0)-1.526e-05;
%     temp_S=cum_struct.S(cum_struct.int>0)-1.526e-05;
% %     temp_int=cum_int';
% [N,~] = hist3([temp_S(:),temp_G(:)],ctrs);
% figure()
% imagesc(x,x,flip(N));
% colormap(gca,jet); ax = gca; ax.Colormap(1,:)= [1 1 1]; caxis('auto');
% % colorbar; axis xy;
% 
% colorbar; axis image;
% 
% x_circle     = [0:1/512:1];
% y_circle_pos = sqrt(0.25-((x_circle-0.5).^2));
% y_circle_neg = -sqrt(0.25-((x_circle-0.5).^2));;
% hold on; plot(x_circle,[y_circle_pos;y_circle_neg],'k','LineWidth',1)
% set(gca,'ydir','reverse')
% axis([0 1 -0.625 0])
% 
% % axis xy;
% % set(gca,'XTick',[0,0.5,1], 'YTick', [0,0.25,0.5])
% % xticks([1/2:1/2^4:1]);
% % xticklabels({'0','', '','','0.5','', '','','1'});
% % 
% % yticks([0:1/2^4:1/2]);
% % yticklabels({'0','', '0.25','','0.5'});
% % xlabel('G');ylabel('S')

 %%      
       
end       
G_sum = [];
S_sum = [];
cond = [];
name = [];
figures = [];
figuresg = [];
G_av_ch1=[];
S_av_ch1=[];
cum_G=[];
cum_S=[];
int=[];
cum_int=[];
S=[];
G=[];
treatment=[];
replicate_group=[];
Sample_ID=[];
if ((tif_length==6 && ch2_pr=='y')) 
    
    
        for z = 1:zStackAmount   %% Groups all the z stacks together. 
           
         cur_int = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),1));
         cur_G = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),2));
         cur_S = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),3));
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
           
 
            
       
            figures(z) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(z), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{((j-1)*3)+z})]));
            saveas(figures(z), fullfile(figures_folder,[strcat("Int_ch1_Image_",files{((j-1)*3)+z}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(z) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(z), fullfile(figures_folder,[strcat("G_ch1_",files{((j-1)*3)+z})]));
            saveas(figuresg(z), fullfile(figures_folder,[strcat("G_ch1_",files{((j-1)*3)+z}(1:end-4),'.png')]));   
            
            
            figuresg(z) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(z), fullfile(figures_folder,[strcat("S_ch1_",files{((j-1)*3)+z})]));
            saveas(figuresg(z), fullfile(figures_folder,[strcat("S_ch1_",files{((j-1)*3)+z}(1:end-4),'.png')])); 
            
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       z_struct = medfiltPhasor(z_struct,3);
       end   
     figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{((j-1)*3)+z})]));
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{((j-1)*3)+z}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(figures_folder,[strcat("Int_ch1_Image_",files{((j-1)*3)+z})]));
    image2=imread(fullfile(figures_folder,[strcat("G_ch1_",files{((j-1)*3)+z})]));
    image3=imread(fullfile(figures_folder,[strcat("single_z_phasor_ch1_",files{((j-1)*3)+z})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(figures_folder,['ch1_singz_',files{((j-1)*3)+z}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);

    if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(figures_folder,['ch1_singlez_',files{((j-1)*3)+z}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
    end
   
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch1 = [G_av_ch1;GZ];
       S_av_ch1 = [S_av_ch1;SZ];
       cond_find=strfind(char_list_of_fold,'_');    
       Sample_ID = [Sample_ID;string(char_list_of_fold(add_ind(end-1)+1:cond_find(end-1)-1))];
       name = [name;string(files{((j-1)*3)+z})];  
       treatment=[treatment;string(treatment_temp)];
       replicate_group=[replicate_group;string(replicate_group_temp)];
    
       end
       und_ind=strfind(files{((j-1)*3)+1},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(figures_folder,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '_ch1.xlsx');   % Name of the excel file. 
       writetable(DataTable,filefolder,'Sheet',1)
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
%        clear cond name G_av_ch1 S_av_ch1
       current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       current_struct = medfiltPhasor(current_struct,3);
       end
       plotPhasorFast(current_struct);
       
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S.      
       Sample_ID =string(char_list_of_fold(add_ind(end-1)+1:cond_find(end-1)-1));
       name=string([char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_ave_values_',(char_list_of_fold(add_ind(end)+1:end))]);
       treatment=[treatment;string(treatment_temp)];
       replicate_group=[replicate_group;string(replicate_group_temp)];
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);
       % for ISS
       filefolder = fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_',(char_list_of_fold(add_ind(end)+1:end)),'_ch1.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_cumulative_phasor_',(char_list_of_fold(add_ind(end)+1:end)),'_ch1.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_cumulative_phasor_',(char_list_of_fold(add_ind(end)+1:end)),'_ch1.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all
  
%Ch2 data 
G_sum = [];
S_sum = [];
cond = [];
name = [];
figures = [];
figuresg = [];
G_av_ch1=[];
S_av_ch1=[];
cum_G=[];
cum_S=[];
cum_int=[];
S=[];
G=[];
int=[];
treatment=[];
replicate_group=[];
Sample_ID=[];
clear cur_int cur_G cur_S
       
        for z = 1:zStackAmount   %% Groups all the z stacks together. 
           
         cur_int = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),4));
         cur_G = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),5));
         cur_S = double(imread(fullfile(list_of_fold{i},files{((j-1)*3)+z}),6));
                     
           G = cat(2,G,cur_G(:));
           S = cat(2,S,cur_S(:));
           int = cat(2,int,cur_int(:));
        
         
%          cum_G=[cum_G (cur_G(:))'];
%          cum_S=[cum_S (cur_S(:))'];
%          cum_int=[cum_int (cur_int(:))'];

           % Removing background below thresh value
    
           cur_int(cur_int(:,:)<thresh_val) = 0;
           cur_G(cur_int(:,:)<thresh_val)=0;
           cur_S(cur_int(:,:)<thresh_val)=0;
           
           
           cur_G_Average = cur_G_Average + cur_G;
           cur_Int_Average = cur_Int_Average + cur_int;
           cur_S_Average = cur_S_Average+cur_S;
           
 
            
       
            figures(z) = imagesc(cur_int); axis image; colormap(gca, jet); colorbar; axis image; title('Intensity'); set(gca,'XTick',[], 'YTick', []);
            saveas(figures(z), fullfile(figures_folder,[strcat("Int_ch2_Image_",files{((j-1)*3)+z})]));
            saveas(figures(z), fullfile(figures_folder,[strcat("Int_ch2_Image_",files{((j-1)*3)+z}(1:end-4),'.png')]));
            mycolormap = customcolormap(linspace(0,1,14), {'#f01114','#ff36c9','#d336ff','#8336ff','#5260ff','#6aa3f7','#76c9f5','#76f1f5','#76f5d1','#76f59e','#b1ff85','#e7ffb0','#e7ffb0','#000000'});

            figuresg(z) = imagesc(cur_G); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('G'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(z), fullfile(figures_folder,[strcat("G_ch2_",files{((j-1)*3)+z})]));
            saveas(figuresg(z), fullfile(figures_folder,[strcat("G_ch2_",files{((j-1)*3)+z}(1:end-4),'.png')])); 
            
            figuresg(z) = imagesc(cur_S); axis image; colormap(mycolormap); colorbar; caxis([0.35,0.8]); axis image; title('S'); set(gca,'XTick',[], 'YTick', []);
            saveas(figuresg(z), fullfile(figures_folder,[strcat("S_ch2_",files{((j-1)*3)+z})]));
            saveas(figuresg(z), fullfile(figures_folder,[strcat("S_ch2_",files{((j-1)*3)+z}(1:end-4),'.png')])); 
            
       close all
       
       z_struct=struct('int',cur_int(:),'G',cur_G(:),'S',cur_S(:));   
       % Applying a threshold for the phasor plot: 
       z_struct = threshPhasorStruct(z_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       z_struct = medfiltPhasor(z_struct,3);
       end   
     figure()
    plotPhasorFast(z_struct);
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch2_",files{((j-1)*3)+z})]));
    saveas(gcf, fullfile(figures_folder,[strcat("single_z_phasor_ch2_",files{((j-1)*3)+z}(1:end-4),'.png')]));
    close all 
    
    image1=imread(fullfile(figures_folder,[strcat("Int_ch2_Image_",files{((j-1)*3)+z})]));
    image2=imread(fullfile(figures_folder,[strcat("G_ch2_",files{((j-1)*3)+z})]));
    image3=imread(fullfile(figures_folder,[strcat("single_z_phasor_ch2_",files{((j-1)*3)+z})]));
    multiplaneData = cat(3, (image1), (image2),(image3));
%     multiplane_file_name=char(fullfile(figures_folder,[strcat('IIntensity_Image_',files{j+z-1})]));
    multiplane_file_name=char(fullfile(figures_folder,['ch2_singz_',files{((j-1)*3)+z}(1:end-4),'.tif']));
%     multiplane_file_name=char(fullfile(figures_folder,['mp','.tif']));

    clear options;
    options.append = true;
    saveastiff((multiplaneData),multiplane_file_name,options);

     if (ome_tif_activate=='y')
     multiplane_ome_name=char(fullfile(figures_folder,['ch2_singlez_',files{((j-1)*3)+z}(1:end-4),'.ome.tif'])); 
     bfsave(multiplaneData,multiplane_ome_name);
     end
     
       [GZ,SZ] = findCenPhasor(z_struct);  % finding the mean of G and S. 
       G_av_ch1 = [G_av_ch1;GZ];
       S_av_ch1 = [S_av_ch1;SZ];
       cond_find=strfind(char_list_of_fold,'_');
       Sample_ID = [Sample_ID;string(char_list_of_fold(add_ind(end-1)+1:cond_find(end-1)-1))];
       name = [name;string(files{((j-1)*3)+z})];  
       treatment=[treatment;string(treatment_temp)];
       replicate_group=[replicate_group;string(replicate_group_temp)];
       end
       und_ind=strfind(files{((j-1)*3)+1},'_');
       DataTable=table(replicate_group,treatment,Sample_ID,name,G_av_ch1, S_av_ch1);       
%        filefolder = fullfile((head_folder{i}),,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '.xlsx');   % Name of the excel file. 
       filefolder = fullfile(figures_folder,string(files{((j-1)*3)+1}(1:und_ind(end)-1)) + '_ch2.xlsx');   % Name of the excel file. 
       writetable(DataTable,filefolder,'Sheet',1)
       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       treatment=[];
       replicate_group=[];
       Sample_ID=[];
%        clear cond name G_av_ch1 S_av_ch1
       current_struct = struct('int',int,'G',G,'S',S);

       
         % Applying a threshold for the phasor plot: 
       current_struct = threshPhasorStruct(current_struct, thresh_val); 
       if (filt_active=='y')
       % Applying median filter for the phasor plot: 
       current_struct = medfiltPhasor(current_struct,3);
       end
       plotPhasorFast(current_struct);
       
       [Organoid_ave_G,Organoid_ave_S] = findCenPhasor(current_struct);  % finding the mean of G and S. 
       
       Sample_ID =string(char_list_of_fold(add_ind(end-1)+1:cond_find(end-1)-1));
       name=string([char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_ave_values_',(char_list_of_fold(add_ind(end)+1:end))]);
       treatment=[treatment;string(treatment_temp)];
       replicate_group=[replicate_group;string(replicate_group_temp)];
       DataTable=table(replicate_group,treatment,Sample_ID,name,Organoid_ave_G, Organoid_ave_S);
       
       
       
       % for ISS
       filefolder = fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_',(char_list_of_fold(add_ind(end)+1:end)),'_ch2.xlsx']);
       writetable(DataTable,filefolder,'Sheet',1)

       figure; 
       plotPhasorFast(current_struct); 
       hold on 
       plot((Organoid_ave_G*256)+256,(Organoid_ave_S*(-256))+256,'r*','MarkerSize',20)
       hold off
       saveas(gcf, fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_cumulative_phasor_',(char_list_of_fold(add_ind(end)+1:end)),'_ch2.tif']));
       saveas(gcf, fullfile((head_folder{i}),[char(files{((j-1)*3)+1}(1:und_ind(end)-1)),'_cumulative_phasor_',(char_list_of_fold(add_ind(end)+1:end)),'_ch2.png']));
%        plot(289,198,'r','*','MarkerSize',10)
       close all   







end



       treatment=[];
       replicate_group=[];
       Sample_ID=[];
S=[];
G=[];
int=[];

       G_av_ch1=[];
       S_av_ch1=[];
       cond = [];
       name = [];
       G=[];
       S=[];
       Organoid_ave_G=[];
       Organoid_ave_S=[];
% clear cond name G_av_ch1 S_av_ch1 name cond G S Organoid_ave_G Organoid_ave_S









  end  
end

end
head_folder=head_folder';
mergedTable_ch1 = table();
mergedTable_ch2 = table();
format long 
for kk=1:numel(list_of_fold)

  if ((kk==1 ||  ~strcmp(head_folder{kk},head_folder{kk-1})) && numel(dir(fullfile(head_folder{kk},"*.xlsx")))~=0)
     
  excell_ref_files = dir(fullfile(head_folder{kk},"*.xlsx"));  
  strings = {excell_ref_files.name};
  % Find the minimum length among the input strings
  minLen = min(cellfun(@length, strings));
  
   % Initialize the common substring
    commonPart = '';

    % Loop through characters at each position
    for i = 1:minLen
        % Get the character at position i from the first string
        currentChar = strings{1}(i);

        % Check if the character appears in all other strings at the same position
        if all(cellfun(@(str) str(i) == currentChar, strings))
            % Append the character to the common substring
            commonPart = [commonPart, currentChar];
        else
            % Break the loop if the characters don't match
            break;
        end
    end
  underline_ind=strfind(commonPart,'_');
  cum_excell_name_ch1=[commonPart(1:underline_ind(end)-1),'_organoid_stats_ch1.xlsx'];
  cum_excell_name_ch2=[commonPart(1:underline_ind(end)-1),'_organoid_stats_ch2.xlsx'];
commonPart='';
clear strings
for h=1:numel(excell_ref_files)
   if (strcmp(excell_ref_files(h).name(end-7:end-5),'ch1'))
    data = readtable(fullfile(head_folder{kk},excell_ref_files(h).name));
    mergedTable_ch1= [mergedTable_ch1; data];
    delete(fullfile(head_folder{kk},excell_ref_files(h).name));
    clear data
   end
   if(strcmp(excell_ref_files(h).name(end-7:end-5),'ch2'))
       data = readtable(fullfile(head_folder{kk},excell_ref_files(h).name));
    mergedTable_ch2= [mergedTable_ch2; data];
    delete(fullfile(head_folder{kk},excell_ref_files(h).name));
    clear data   
   end
end    
  end
 if (sum(size(mergedTable_ch1))~=0)
     
     
       filefolder = fullfile(head_folder{kk},cum_excell_name_ch1);
       writetable(mergedTable_ch1,filefolder,'Sheet',1)    
       mergedTable_ch1 = table();
 end
  
  if (sum(size(mergedTable_ch2))~=0)
     
     
       filefolder = fullfile(head_folder{kk},cum_excell_name_ch2);
       writetable(mergedTable_ch2,filefolder,'Sheet',1)    
       mergedTable_ch2 = table();
 end
  
 clear excell_ref_files 
end

cond_fold=cond_fold';
mergedTable_ch1 = table();
mergedTable_ch2 = table();
%%
for ff=1:numel(list_of_fold)


    if (ff==1 ||  ~strcmp(cond_fold{ff},cond_fold{ff-1}))
   
    cont_files = dir(cond_fold{ff});    
    dirFlags = [cont_files.isdir];    
    subFolders = cont_files(dirFlags);    
    subFolderNames = {subFolders(3:end).name}; 
    subFolderNames=subFolderNames';
    for nn=1:numel(subFolderNames)
    excell_cond_files = dir(fullfile(cond_fold{ff},subFolderNames{nn},"*.xlsx"));  
if (~isempty(excell_cond_files))
   
   for p=1:numel(excell_cond_files)
    
      if (strcmp(excell_cond_files(p).name(end-7:end-5),'ch1')) 
    data = readtable(fullfile(cond_fold{ff},subFolderNames{nn},excell_cond_files(p).name));
    mergedTable_ch1= [mergedTable_ch1; data];
    clear data   
      end
        if (strcmp(excell_cond_files(p).name(end-7:end-5),'ch2')) 
    data = readtable(fullfile(cond_fold{ff},subFolderNames{nn},excell_cond_files(p).name));
    mergedTable_ch2= [mergedTable_ch2; data];
       
    clear data
      end    
      
      
       
   end    
        
end       
        
        
    end
    
    
    
  if (sum(size(mergedTable_ch1))~=0)
     
     
       filefolder = fullfile(cond_fold{ff},'all_plates_stats_ch1.xlsx');
       writetable(mergedTable_ch1,filefolder,'Sheet',1)    
       mergedTable_ch1 = table();
  end
  if (sum(size(mergedTable_ch2))~=0)
     
     
       filefolder = fullfile(cond_fold{ff},'all_plates_stats_ch2.xlsx');
       writetable(mergedTable_ch2,filefolder,'Sheet',1)    
       mergedTable_ch2 = table();
  end  
  
  
  
  
   

  
  
  
  
  
  
  
    end 


clear cont_files excell_cond_files




end

%%
mergedTable_ch1 = table();
mergedTable_ch2 = table();
data=[];
replicate_group=[];
cond_slash_ind=strfind(cond_fold,'\'); 
for t=1:numel(cond_fold)
  replicate_group{t}=cond_fold{t}(1:cond_slash_ind{t}(end)-1);  
end   
replicate_group=replicate_group';
for q=1:numel(replicate_group)
    
   if (q==1 ||  ~strcmp(replicate_group{q},replicate_group{q-1}))  
    
    rep_group_dir = dir(replicate_group{q});    
    rep_dirFlags = [rep_group_dir.isdir];  
    
    rep_subFolders = rep_group_dir(rep_dirFlags);    
    rep_subFolderNames = {rep_subFolders(3:end).name}; 
    rep_subFolderNames=rep_subFolderNames';
     for p=1:numel(rep_subFolderNames)
    excell_rep_group_files = dir(fullfile(replicate_group{q},rep_subFolderNames{p},"*.xlsx")); 
      
     if (~isempty(excell_rep_group_files))
   
   for gg=1:numel(excell_rep_group_files)
    
      if (strcmp(excell_rep_group_files(gg).name(end-7:end-5),'ch1')) 
    data = readtable(fullfile(replicate_group{q},rep_subFolderNames{p},excell_rep_group_files(gg).name));
    mergedTable_ch1= [mergedTable_ch1; data];
    clear data   
      end
        if (strcmp(excell_rep_group_files(gg).name(end-7:end-5),'ch2')) 
    data = readtable(fullfile(replicate_group{q},rep_subFolderNames{p},excell_rep_group_files(gg).name));
    mergedTable_ch2= [mergedTable_ch2; data];
       
    clear data
        end    
      
      
       
   end    
        
    end
     
     
     end

     
     
       if (sum(size(mergedTable_ch1))~=0)
     
     
       filefolder = fullfile(replicate_group{q},'all_replicate_group_stats_ch1.xlsx');
       writetable(mergedTable_ch1,filefolder,'Sheet',1)    
       mergedTable_ch1 = table();
       end
       
       if (sum(size(mergedTable_ch2))~=0)
     
     
       filefolder = fullfile(replicate_group{q},'all_replicate_group_stats_ch2.xlsx');
       writetable(mergedTable_ch2,filefolder,'Sheet',1)    
       mergedTable_ch2 = table();
       end      
     
  end
    
end    



%%

% Define the filename of the m4a audio file
inputFilename = 'accomplish.mp3'; % Replace with the actual file path

% Read the m4a audio file
[y, Fs] = audioread(inputFilename);

% Create an audioplayer object for the audio
player = audioplayer(y, Fs);

% Play the audio file twice
play(player); % First play
pause(12); % Wait for the first play to finish
% Close the audio player object
% Stop the player
stop(player);
delete(player);
fprintf('\n');
fprintf('\n');
fprintf('All the files are analyzed successfully\n');




function folderAddresses = findEndSubfolders(parentFolder)
    % Initialize an empty cell array to store folder addresses
    folderAddresses = {};

    % Get a list of all items within the parent folder
    items = dir(parentFolder);

    % Iterate through each item
    for i = 1:numel(items)
        % Ignore current and parent directories
        if strcmp(items(i).name, '.') || strcmp(items(i).name, '..')
            continue;
        end

        % Check if the item is a folder
        if items(i).isdir
            % Get the full path of the subfolder
            subfolderPath = fullfile(parentFolder, items(i).name);

            % Recursively call the function for subfolders
            subfolderAddresses = findEndSubfolders(subfolderPath);

            % If there are no further subfolders, add the current subfolder to the addresses
            if isempty(subfolderAddresses)
                folderAddresses{end+1} = subfolderPath;
            else
                % Append the subfolder addresses to the existing addresses
                folderAddresses = [folderAddresses subfolderAddresses];
            end
        end
    end
end