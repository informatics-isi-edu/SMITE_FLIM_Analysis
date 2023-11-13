function sorted_ind = stack_indices(files,n,zstack,ch,str)

if (strcmp(str,'iss'))
underline_ind_z=strfind(files,'_');
temp_ind=[];
for bb = 1:numel(files)
if (strcmp(files{bb}(underline_ind_z{bb}(end-3):underline_ind_z{bb}(end-1)),'_DC_Ch1_'))
 temp_ind=[temp_ind bb];   
end
end
sorted_ind(1)=temp_ind(n)


for bb = 1:numel(files)
if (strcmp([files{bb}(1:underline_ind_z{bb}(end-3)-1),files{bb}(underline_ind_z{bb}(end-2):end)],[files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end-3)-1),files{sorted_ind(1)}(underline_ind_z{sorted_ind(1)}(end-2):end)]) && strcmp(files{bb}(underline_ind_z{bb}(end-3):underline_ind_z{bb}(end-1)),'_G_Ch1_'))
 sorted_ind(2)=bb;   
end
end

for bb = 1:numel(files)
if (strcmp([files{bb}(1:underline_ind_z{bb}(end-3)-1),files{bb}(underline_ind_z{bb}(end-2):end)],[files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end-3)-1),files{sorted_ind(1)}(underline_ind_z{sorted_ind(1)}(end-2):end)]) && strcmp(files{bb}(underline_ind_z{bb}(end-3):underline_ind_z{bb}(end-1)),'_S_Ch1_'))
 sorted_ind(3)=bb;   
end
end

if (ch==2)
    
 for bb = 1:numel(files)
if (strcmp([files{bb}(1:underline_ind_z{bb}(end-3)-1),files{bb}(underline_ind_z{bb}(end-1):end)],[files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end-3)-1),files{sorted_ind(1)}(underline_ind_z{sorted_ind(1)}(end-1):end)]) && strcmp(files{bb}(underline_ind_z{bb}(end-3):underline_ind_z{bb}(end-1)),'_DC_Ch2_'))
 sorted_ind(4)=bb;   
end
end   
   
for bb = 1:numel(files)
if (strcmp([files{bb}(1:underline_ind_z{bb}(end-3)-1),files{bb}(underline_ind_z{bb}(end-1):end)],[files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end-3)-1),files{sorted_ind(1)}(underline_ind_z{sorted_ind(1)}(end-1):end)]) && strcmp(files{bb}(underline_ind_z{bb}(end-3):underline_ind_z{bb}(end-1)),'_G_Ch2_'))
 sorted_ind(5)=bb;   
end
end    
    
 for bb = 1:numel(files)
if (strcmp([files{bb}(1:underline_ind_z{bb}(end-3)-1),files{bb}(underline_ind_z{bb}(end-1):end)],[files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end-3)-1),files{sorted_ind(1)}(underline_ind_z{sorted_ind(1)}(end-1):end)]) && strcmp(files{bb}(underline_ind_z{bb}(end-3):underline_ind_z{bb}(end-1)),'_S_Ch2_'))
 sorted_ind(6)=bb;   
end
end         
    
end



end


if (strcmp(str,'leica'))
underline_ind_z=strfind(files,'_');
temp_ind=[];
for bb = 1:numel(files)
if (strcmp(files{bb}(underline_ind_z{bb}(end)+1:end-4),'ch0'))
 temp_ind=[temp_ind bb];   
end
end
sorted_ind(1)=temp_ind(n)

for bb = 1:numel(files)
if (strcmp(files{bb}(1:underline_ind_z{bb}(end)-1),files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end)-1)) && strcmp(files{bb}(underline_ind_z{bb}(end)+1:underline_ind_z{bb}(end)+3),'ch1'))
 sorted_ind(2)=bb;   
end
end

for bb = 1:numel(files)
if (strcmp(files{bb}(1:underline_ind_z{bb}(end)-1),files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end)-1)) && strcmp(files{bb}(underline_ind_z{bb}(end)+1:underline_ind_z{bb}(end)+3),'ch2'))
 sorted_ind(3)=bb;   
end
end

if (ch==2)
    
 for bb = 1:numel(files)
if (strcmp(files{bb}(1:underline_ind_z{bb}(end)-1),files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end)-1)) && strcmp(files{bb}(underline_ind_z{bb}(end)+1:underline_ind_z{bb}(end)+3),'ch3'))
 sorted_ind(4)=bb;   
end
end  
   
 for bb = 1:numel(files)
if (strcmp(files{bb}(1:underline_ind_z{bb}(end)-1),files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end)-1)) && strcmp(files{bb}(underline_ind_z{bb}(end)+1:underline_ind_z{bb}(end)+3),'ch4'))
 sorted_ind(5)=bb;   
end
end    
    
 for bb = 1:numel(files)
if (strcmp(files{bb}(1:underline_ind_z{bb}(end)-1),files{sorted_ind(1)}(1:underline_ind_z{sorted_ind(1)}(end)-1)) && strcmp(files{bb}(underline_ind_z{bb}(end)+1:underline_ind_z{bb}(end)+3),'ch5'))
 sorted_ind(6)=bb;   
end
end        
    
end



end 
end