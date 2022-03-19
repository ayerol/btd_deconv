function [x,state] = struct_toeplitz_slices(z,task)



%   Authors: Borbala Hunyadi (b.hunyadi@tudelft.nl) and 

%   Aybuke Erol (a.erol@tudelft.nl) 



%   This function transforms the columns of z into a Toeplitz matrix. 

%   These matrices are later placed into frontal slices of the core tensor. 

%   Note that, each column of z was already ensured to be a lagged

%   version of the preceeding column using a struct previously defined

%   (check ttoepl_cores1 in btd_deconv). This way, the frontal slices of

%   the core tensor follows a similar shifting pattern.



state = [];
K=size(z,2);
size_mat = [size(z,1) size(z,1)];


for k=1:K


    z1 = z(:,k);
    
    if ~isempty(task) && ~isempty(task.r)
        
        task1 = task;
        task1.r = squeeze(task1.r(:,:,k));
        
    end
    
    if ~isempty(task) && ~isempty(task.l)
        
        task1 = task;
        task1.l = squeeze(task1.l(:,:,k));
        z_inp = circshift(z(:,1),1-k);
        pre = flipud(z_inp(2:end));
        pre(1:k-1) = z(end-k+2:end,1);
        
        x1{k} = struct_toeplitz(z1,task1,size_mat,pre);
        
    end
    
    if isempty(task)
        
        task1 = [];
        z_inp = circshift(z(:,1),1-k);
        pre = flipud(z_inp(2:end));
        pre(1:k-1) = z(end-k+2:end,1);
        
        x1{k} = struct_toeplitz(z1,task1,size_mat,pre);
        
    end
    
end


if (size(x1{1},2))>1
    
    x=zeros(size(z,1),size(z,1),size(z,2));
    
    for k=1:size(z,2)
        
        x(:,:,k)=x1{k};
        
    end    
    
else
    
    x=cell2mat(x1);
    
end


end