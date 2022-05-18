function [sol,cost] = btd_deconv(T,M,Lacc,u,source_list)


% REMARK:

% This code is an extension of the convolutive mixtures demo offered by 

% Tensorlab (check [1]). The necessary adjustments are made by 

% Aybuke erol (a.erol@tudelft.nl)


L = length(u) - 1; K = size(T,3); R = length(source_list);
impulse = zeros(L+1,1);
impulse(floor((L+1)/2)+1) = 1;
initPnt = zeros(3,M);


% Initialize H factors

for m = 1:M
    
    for r = 1:R
        
        if source_list{r} == 't'
            
            % Initialize the HRF parameters to be estimated

            initPnt(:,m) = [2+4.75*rand rand 1.25+3.25*rand];
            model.variables.(['H' num2str(m) num2str(r)])= initPnt(:,m)';
            
        else
            
            model.variables.(['H' num2str(m) num2str(r)]) = impulse;
            
        end
        
    end
    
end


% Initialize core tensors

for r = 1:R
    
    model.variables.(['Core' num2str(r)]) = randn(L+Lacc,1);
    
end


% % Define constraints for all factors

% Apply gamma-model for HRFs
thrf = @(z,task)struct_hrf(z,task,u);

% Assign only non-negative parameters for the gamma model
tnonneg = @(z,task)struct_abs(z,task,0.001);

% Make sure that the impulsive filters for artifact sources will only
% change by the scaling
tconst = @(z,task)struct_const(z,task,impulse);

% Construct the H matrix by computing the Toeplitz matricized versions of
% the convolutive mixing filters (hence, both the HRFs and impulses)
ttoepl = @(z,task)struct_toeplitz(z,task,[Lacc L+Lacc],zeros(Lacc-1,1),...
    zeros(Lacc-1,1));

% Make sure that the adjacent columns of z are lagging one another
ttoepl_cores1 = @(z,task)struct_toeplitz(z,task,[L+Lacc,K],flipud(z(2:K)));

% Use the columns of z to form Toeplitz-structured frontal slices of the
% core tensor
ttoepl_cores2 = @(z,task)struct_toeplitz_slices(z,task); 


% Define H factors

for r = 1:R
    
    factor_H{r} = cell(M,1);
    
    for m = 1:M
        
        if source_list{r} == 't'
            
            factor_H{r}{m} = {['H' num2str(m) num2str(r)],tnonneg,thrf,...
                ttoepl};
            
        else
            
            factor_H{r}{m} = {['H' num2str(m) num2str(r)],tconst,ttoepl};
            
        end
        
    end
    
    model.factors.(['H' num2str(r)])=factor_H{r};
    
end


% Define the remaining factors (core tensors and identity)

for r = 1:R
    
    model.factors.(['Core' num2str(r)])={['Core' num2str(r)],...
        ttoepl_cores1,ttoepl_cores2};
    
end

model.factors.C = eye(K);


% Perform BTD

model.factorizations.tensor.data = T;

exbtd = cell(1,R);
for r = 1:R

    exbtd{r} = {['H' num2str(r)],['H' num2str(r)],'C',['Core' num2str(r)]};

end

model.factorizations.tensor.btd = exbtd;
options.Display = 0;
options.MaxIter = 500;
options.TolFun = 1e-8;
options.TolX = 1e-8;
[sol,output] = sdf_minf(model,options);
cost = output.fval(end);