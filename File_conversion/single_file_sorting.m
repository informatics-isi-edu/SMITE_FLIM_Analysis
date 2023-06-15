clear all 
close all 
clc
format long
prompt = {'Please enter the address to the folder containing the data','Please enter the data folder name','Please enter the number of z stacks'};
answer = inputdlg(prompt);
address=string(answer(1));

data_name=string(answer(2));
zstack=str2double(answer(3));
ref_folder = fullfile(address,data_name);
cont_folder = dir(ref_folder);
output_fold=fullfile(address,[answer{2}, '_concatenated files']);
mkdir (output_fold)
%%
for jj=3:length(cont_folder)
mkdir (fullfile(output_fold,[cont_folder(jj).name,'_concatenated files']))

img_files=dir(fullfile(cont_folder(jj).folder,cont_folder(jj).name,"*.tif")); 

underline_ind=strfind(img_files(1).name,'_'); 
if (img_files(1).name((underline_ind(end)+1:underline_ind(end)+2))=='ch')  

files=sort_nat({img_files.name},'ascend');
underline_ind_plns=strfind(files,'_');
pln_detect=cell(size(files));
for i = 1:numel(files)
    pln{i} =files{i}(1:underline_ind_plns{i}(end)-1);
end

indices=strcmp(files{1}(1:underline_ind_plns{1}(end)-1),pln);
if (sum(indices)==3)

for ii=1:3:length(files)
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii});
int=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii+1});
G=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii+2});
S=imread(file_name);
ind=strfind(files{ii},'_');
single_file_name=char(fullfile(output_fold,[cont_folder(jj).name,'_concatenated files'],[files{ii}(1:ind(end)-1),'.tif']));
temp_plane(:,:,1)=single((int));
temp_plane(:,:,2)=single(standardPhase(G));
temp_plane(:,:,3)=single(standardPhase(S));
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
% imwrite(temp_plane,single_file_name,'writemode', 'append')
clear temp_plane G int S img_files
end
else 
for ii=1:6:length(files)
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii});
int_1=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii+1});
G_1=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii+2});
S_1=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii+3});
int_2=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii+4});
G_2=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{ii+5});
S_2=imread(file_name);
ind=strfind(files{ii},'_');
single_file_name=char(fullfile(output_fold,[cont_folder(jj).name,'_concatenated files'],[files{ii}(1:ind(end)-1),'.tif']));
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
clear temp_plane G int S img_files
end
clear indices
end
clear files

else 
files=sort_nat({img_files.name},'ascend');
   
for kk=1:length(files)/zstack
       
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{kk});
int=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{length(files)/zstack+kk});
G=imread(file_name);
file_name=fullfile(cont_folder(jj).folder,cont_folder(jj).name,files{2*(length(files)/zstack)+kk});
S=imread(file_name);
ind=strfind(files{kk},'_');
single_file_name=char(fullfile(output_fold,[cont_folder(jj).name,'_concatenated files'],[files{kk}(1:ind(4)),files{kk}(ind(5)+1:end)]));
temp_plane(:,:,1)=int;
temp_plane(:,:,2)=G;
temp_plane(:,:,3)=S;
clear options;
options.append = true;
saveastiff(temp_plane,single_file_name,options);
clear temp_plane G int S img_files
end
clear files   
    
    
    
    
    

end
end