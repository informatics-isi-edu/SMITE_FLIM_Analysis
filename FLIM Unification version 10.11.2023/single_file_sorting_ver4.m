clear all 
close all 
clc
format long
prompt = {'Please enter the address to the folder containing the data','Please enter the data folder name','Please enter the number of z stacks (enter 0 if each folder is representing single organoid'};
answer = inputdlg(prompt);
address=string(answer(1));

data_name=string(answer(2));
zstack=str2double(answer(3));
ref_folder = fullfile(address,data_name);

prompt={'Would you like to export data from each organoid in separate folders? (y/n) '}

fold_check_pr = inputdlg(prompt);
fold_check=string(fold_check_pr(1));

%% find the adsdresses of subfolders


% FLIM-ISS_Time-decay
% cont_folder=fullfile(current_folder,data_folder);
list_of_fold=(findEndSubfolders(ref_folder))';

for i=1:numel(list_of_fold)
folderAddresse(i) = erase(list_of_fold{i}, ref_folder);
end

folderAddresse=folderAddresse';



 %% remove Hr2 and Hr3 files 

for mm=1:length(list_of_fold) 

 tif_files=dir(fullfile(list_of_fold{mm},"*.tif"));   
 for w=1:length(tif_files)
    if (strcmp(tif_files(w).name(end-6:end-4),'Hr2') || strcmp(tif_files(w).name(end-6:end-4),'Hr3'))
        
       delete(fullfile(tif_files(w).folder,tif_files(w).name)) 
    end  
    
 end
clear tif_files
end

%%
% cont_folder = dir(ref_folder);
output_fold=fullfile(address,[answer{2}, '_unified files']);
mkdir (output_fold)
%%
for jj=1:length(folderAddresse)
    
slash_ind=strfind(list_of_fold{jj},'\'); 
char_add=char(list_of_fold{jj});
char_short_add=char(folderAddresse(jj));
slash_ind_short_add=strfind(folderAddresse(jj),'\');
if (strcmp(char_add(slash_ind(end)+1:end),'MetaData') && length(dir(fullfile(char_add(1:slash_ind(end)-1),"*.tif")))==0)
% if (exist(char_add(1:slash_ind(end)-1), 'dir') == 7)
    
   copyfile(list_of_fold{jj}, fullfile((output_fold),char_short_add)); 
   continue;
% else
    
%  fullfile(output_fold,[char(folderAddresse(jj))   
% end    
    
end    
if (strcmp(char_add(slash_ind(end)+1:end),'MetaData') && length(dir(fullfile(char_add(1:slash_ind(end)-1),"*.tif")))~=0)

mkdir (fullfile(output_fold,char(folderAddresse(jj))))
copyfile(list_of_fold{jj}, fullfile((output_fold),char_short_add));
    if (zstack==0)
    img_files=dir(fullfile(char_add(1:slash_ind(end)-1),"*.tif"));   
 
    
        
    underline_ind=strfind(img_files(1).name,'_'); 
if (img_files(1).name((underline_ind(end)+1:underline_ind(end)+2))=='ch')  
str='leica';
files=(sort_nat({img_files.name},'ascend'))';
underline_ind_plns=strfind(files,'_');
pln_detect=cell(size(files));
for i = 1:numel(files)
    pln{i} =files{i}(1:underline_ind_plns{i}(end)-1);
end

indices=strcmp(files{1}(1:underline_ind_plns{1}(end)-1),pln);
if (sum(indices)==3)
ch=1;
for ii=1:length(files)/3
% stack_indices(files,n,zstacks,ch,str) 
% stack_indices(files,n,zstacks,ch,str)
stack_ind=stack_indices(files,ii,zstack,ch,str);  
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(1)});
int=imread(file_name);
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(2)});
G=imread(file_name);
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(3)});
S=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');
single_file_name=char(fullfile(output_fold,char_short_add(1:slash_ind_short_add(end)-1),[files{stack_ind(1)}(1:ind(end)-1),'.tif']));
temp_plane(:,:,1)=single((int));
temp_plane(:,:,2)=single(standardPhase(G));
temp_plane(:,:,3)=single(standardPhase(S));
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
% imwrite(temp_plane,single_file_name,'writemode', 'append')
clear temp_plane G int S img_files stack_ind
end
else 
    ch=2;
for ii=1:(length(files)/3)/2


% stack_indices(files,n,zstacks,ch,str)
stack_ind=stack_indices(files,ii,zstack,ch,str);  
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(1)});
int_1=imread(file_name);
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(2)});
G_1=imread(file_name);
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(3)});
S_1=imread(file_name);
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(4)});
int_2=imread(file_name);
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(5)});
G_2=imread(file_name);
file_name=fullfile(char_add(1:slash_ind(end)-1),files{stack_ind(6)});
S_2=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');
single_file_name=char(fullfile(output_fold,char_short_add(1:slash_ind_short_add(end)-1),[files{stack_ind(1)}(1:ind(end)-1),'.tif']));
temp_plane(:,:,1)=single((int_1));
temp_plane(:,:,2)=single(standardPhase(G_1));
temp_plane(:,:,3)=single(standardPhase(S_1));
temp_plane(:,:,4)=single((int_2));
temp_plane(:,:,5)=single(standardPhase(G_2));
temp_plane(:,:,6)=single(standardPhase(S_2));
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
% imwrite(temp_plane,single_file_name,'writemode', 'append')
clear temp_plane G_1 int_1 S_1 G_2 int_2 S_2 img_files stack_ind
end
clear indices
end
clear files

else 
    
 str='iss';   
files=(sort_nat({img_files.name},'ascend'))';

underline_ind_plns=strfind(files,'_');

for i = 1:numel(files)
    pln{i} =[files{i}(1:underline_ind_plns{i}(end-3)-1),files{i}(underline_ind_plns{i}(end-1):end)];
end

indices=strcmp([files{1}(1:underline_ind_plns{1}(end-3)-1),files{1}(underline_ind_plns{1}(end-1):end)],pln);
 
if (sum(indices)==3)
    ch=1;
for kk=1:length(files)/3
 stack_ind=stack_indices(files,kk,zstack,ch,str);  

file_name=fullfile(list_of_fold{jj},files{stack_ind(1)});
int=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(2)});
G=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(3)});
S=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');
single_file_name=char(fullfile(output_fold,char_short_add(1:slash_ind_short_add(end)-1),[files{stack_ind(1)}(1:ind(end)-1),'.tif']));
temp_plane(:,:,1)=int;
temp_plane(:,:,2)=G;
temp_plane(:,:,3)=S;
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
clear temp_plane G int S G_2 int_2 S_2 img_files
end
else
    
    ch=2;
for kk=1:(length(files)/3)/2
       
stack_ind=stack_indices(files,kk,zstack,ch,str);  

file_name=fullfile(list_of_fold{jj},files{stack_ind(1)});
int_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(2)});
G_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(3)});
S_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(4)});
int_2=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(5)});
G_2=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(6)});
S_2=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');

single_file_name=char(fullfile(output_fold,char_short_add(1:slash_ind_short_add(end)-1),[files{stack_ind(1)}(1:ind(end)-1),'.tif']));
temp_plane(:,:,1)=int_1;
temp_plane(:,:,2)=G_1;
temp_plane(:,:,3)=S_1;
temp_plane(:,:,4)=int_2;
temp_plane(:,:,5)=G_2;
temp_plane(:,:,6)=S_2;
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
clear temp_plane G int S G_2 int_2 S_2 img_files
end


clear temp_plane G int S G_2 int_2 S_2 img_files
end
clear files   
     

end
      
    
    
    end

end


if (~strcmp(char_add(slash_ind(end)+1:end),'MetaData'))


mkdir (fullfile(output_fold,[char(folderAddresse(jj)),'_unified files']))


img_files=dir(fullfile(list_of_fold{jj},"*.tif")); 

underline_ind=strfind(img_files(1).name,'_'); 
if (img_files(1).name((underline_ind(end)+1:underline_ind(end)+2))=='ch')  
str='leica';
files=(sort_nat({img_files.name},'ascend'))';
underline_ind_plns=strfind(files,'_');
pln_detect=cell(size(files));
for i = 1:numel(files)
    pln{i} =files{i}(1:underline_ind_plns{i}(end)-1);
end

indices=strcmp(files{1}(1:underline_ind_plns{1}(end)-1),pln);
if (sum(indices)==3)
ch=1;
for ii=1:length(files)/3
% stack_indices(files,n,zstacks,ch,str) 
% stack_indices(files,n,zstacks,ch,str)
stack_ind=stack_indices(files,ii,zstack,ch,str);  
file_name=fullfile(list_of_fold{jj},files{stack_ind(1)});
int=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(2)});
G=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(3)});
S=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');
single_file_name=char(fullfile(output_fold,[char(folderAddresse(jj)),'_unified files'],[files{stack_ind(1)}(1:ind(end)-1),'.tif']));
temp_plane(:,:,1)=single((int));
temp_plane(:,:,2)=single(standardPhase(G));
temp_plane(:,:,3)=single(standardPhase(S));
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
% imwrite(temp_plane,single_file_name,'writemode', 'append')
clear temp_plane G int S img_files stack_ind
end
else 
    ch=2;
for ii=1:(length(files)/3)/2


% stack_indices(files,n,zstacks,ch,str)
stack_ind=stack_indices(files,ii,zstack,ch,str);  
file_name=fullfile(list_of_fold{jj},files{stack_ind(1)});
int_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(2)});
G_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(3)});
S_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(4)});
int_2=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(5)});
G_2=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(6)});
S_2=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');
single_file_name=char(fullfile(output_fold,[char(folderAddresse(jj)),'_unified files'],[files{stack_ind(1)}(1:ind(end)-1),'.tif']));
temp_plane(:,:,1)=single((int_1));
temp_plane(:,:,2)=single(standardPhase(G_1));
temp_plane(:,:,3)=single(standardPhase(S_1));
temp_plane(:,:,4)=single((int_2));
temp_plane(:,:,5)=single(standardPhase(G_2));
temp_plane(:,:,6)=single(standardPhase(S_2));
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
% imwrite(temp_plane,single_file_name,'writemode', 'append')
clear temp_plane G_1 int_1 S_1 G_2 int_2 S_2 img_files stack_ind
end
clear indices
end
clear files

else 
    
 str='iss';   
files=(sort_nat({img_files.name},'ascend'))';

underline_ind_plns=strfind(files,'_');

for i = 1:numel(files)
    pln{i} =[files{i}(1:underline_ind_plns{i}(end-3)-1),files{i}(underline_ind_plns{i}(end-1):end)];
end

indices=strcmp([files{1}(1:underline_ind_plns{1}(end-3)-1),files{1}(underline_ind_plns{1}(end-1):end)],pln);
 
if (sum(indices)==3)
    ch=1;
for kk=1:length(files)/3
 stack_ind=stack_indices(files,kk,zstack,ch,str);  

file_name=fullfile(list_of_fold{jj},files{stack_ind(1)});
int=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(2)});
G=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(3)});
S=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');
single_file_name=char(fullfile(output_fold,[char(folderAddresse(jj)),'_unified files'],[files{stack_ind(1)}(1:ind(end-3)),files{stack_ind(1)}(ind(end-2)+1:end)]));
temp_plane(:,:,1)=int;
temp_plane(:,:,2)=G;
temp_plane(:,:,3)=S;
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
clear temp_plane G int S G_2 int_2 S_2 img_files
end
else
    
    ch=2;
for kk=1:(length(files)/3)/2
       
stack_ind=stack_indices(files,kk,zstack,ch,str);  

file_name=fullfile(list_of_fold{jj},files{stack_ind(1)});
int_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(2)});
G_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(3)});
S_1=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(4)});
int_2=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(5)});
G_2=imread(file_name);
file_name=fullfile(list_of_fold{jj},files{stack_ind(6)});
S_2=imread(file_name);
ind=strfind(files{stack_ind(1)},'_');

single_file_name=char(fullfile(output_fold,[char(folderAddresse(jj)),'_unified files'],[files{stack_ind(1)}(1:ind(end-3)-1),files{stack_ind(1)}(ind(end-1):end)]));
temp_plane(:,:,1)=int_1;
temp_plane(:,:,2)=G_1;
temp_plane(:,:,3)=S_1;
temp_plane(:,:,4)=int_2;
temp_plane(:,:,5)=G_2;
temp_plane(:,:,6)=S_2;
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
clear temp_plane G int S G_2 int_2 S_2 img_files
end


clear temp_plane G int S G_2 int_2 S_2 img_files
end
clear files   
    
    
    
    
    

end

end

end

%% move data from each organoid into a single folder

% find all the subfolders in the output folder
if ((fold_check=='y') && (zstack~=0))
list_output_fold=(findEndSubfolders(output_fold))';


for i=1:numel(list_output_fold)
out_fold_address(i) = erase(list_output_fold{i}, output_fold);
end

out_fold_address=out_fold_address';

%%
for gg=1:numel(list_output_fold)
    
out_img_files=dir(fullfile(list_output_fold{gg},"*.tif")); 
out_files=(sort_nat({out_img_files.name},'ascend'))';
for w=1:(numel(out_files))/zstack
org_fold=[fullfile(list_output_fold{gg},['Organoid_',num2str(w)])];
mkdir ([fullfile(list_output_fold{gg},['Organoid_',num2str(w)])])
for q=1:zstack
 
source_file=fullfile(list_output_fold{gg},out_files{((w-1)*zstack)+q}); 
dest_file=fullfile(org_fold,out_files{((w-1)*zstack)+q});   
 movefile(source_file, dest_file);   
end   
    
end

img_files=dir(fullfile(list_of_fold{jj},"*.tif")); 

underline_ind=strfind(img_files(1).name,'_'); 
end


end

%%

inputFilename = 'accomplish.mp3'; 

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
fprintf('All the files are unified successfully\n');

%%

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
