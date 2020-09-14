function [Wb, train_err, test_err] = back_prop ( train_data, output, label,lambda, test_data, ...
    Wb, n, num_iter, batch_size, err_type )
% function [Wb, train_err, test_err] = back_prop ( train_data, test_data, ...
%   Wb, n, num_iter, batch_size, err_type );
% This function performs fine tuning of the autoencoder using
% back propagation of errors.
%
% INPUTS: train_data, test_data -- training and test data in the
%              format D x N, where D = dimension, N = number of samples. 
%         Wb -- weights and biases as described in readme.txt
%         n -- vector with autoencoder layer sizes (including original
%              data dimension).
%         num_iter -- number of iterations to use for optimization
%         batch_size -- size of batches to use (optional -- default
%              is one batch with all data)
%         err_type -- error to minimize (defaults to 1)
%                       cross_entropy = 1,
%                       reconstruction = 2.
%
% OUTPUTS: Wb -- updates weights and biases
%          train_err, test_err -- a vector with reconstruction error
%               per iteration.
%
% NOTE: Training set size must be divisible by batch size, batch size
%       is automatically determined for test set.
%


% Version 1.000
%
% Code provided by Ruslan Salakhutdinov and Geoff Hinton
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from our
% web page.
% The programs and documents are distributed without any warranty, express or
% implied.  As the programs were written for research purposes only, they have
% not been tested to the degree that would be advisable in any important
% application.  All use of these programs is entirely at the user's own risk.

% number line searchs is a permanent "compile time" value
NUM_LS = 3;

% this is the desired batch size if not specified
DES_BATCH_SIZE = 1000;

% do we split the datasets into batches?
if nargin >= 9
    
    % training set values
    num_train_cases = batch_size;
    num_train_batches = size(train_data,2)/num_train_cases;
    if num_train_batches ~= round(num_train_batches)
        error ('batch_size is not a factor of training set size!');
    end
    
    % test set division is automatic
    N = size(test_data,2);
    f = factor(N);
    rev_f = f(length(f):-1:1);
    bs_ind = min(find(cumprod(rev_f)>=DES_BATCH_SIZE));
    if ~isempty(bs_ind)
        num_test_cases = prod(rev_f(1:bs_ind));
    else
        num_test_cases = N;
    end
    num_test_batches = size(test_data,2)/num_test_cases;

else
    
    num_train_cases = size(train_data,2);
    num_train_batches = 1;
    
    num_test_cases = size(test_data,2);
    num_test_batches = 1;
    
end

if nargin <10
    err_type = 1;
end

fprintf('Fine-tuning autoencoder by minimizing cross entropy error.\n');
fprintf('There are %d batches with %d cases each.\n',num_train_batches,...
    num_train_cases);



% do specified number of iterations
for curr_iter = 1:num_iter

    % compute training reconstruction error
    AE_err = 0;
    for curr_batch = 1:num_train_batches
        
        % get current batch
        curr_data = train_data(:,(1+(curr_batch-1)*num_train_cases): ...
            (curr_batch*num_train_cases));
        %out_data=output(:,(1+(curr_batch-1)*num_train_cases): ...
          %  (curr_batch*num_train_cases));
  
        % do forward and inverse maps on batch
        for_data = AE_forward ( curr_data, Wb, n );
        recon_data = AE_inverse ( for_data, Wb, n );
        AE_err = AE_err + 1/num_train_cases*...
            sum(sum((curr_data - recon_data).^2)); 
       train_layer(:,(1+(curr_batch-1)*num_train_cases): ...
            (curr_batch*num_train_cases))=for_data;
    end
    train_err(curr_iter)=AE_err/num_train_batches;
    
    % compute test set reconstruction error
    AE_err = 0;
    for curr_batch = 1:num_test_batches
        
        % get current batch
        curr_data = test_data(:,(1+(curr_batch-1)*num_test_cases): ...
            (curr_batch*num_test_cases));
  
        % do forward and inverse maps on batch
        for_data = AE_forward ( curr_data, Wb, n );
        recon_data = AE_inverse ( for_data, Wb, n );
        AE_err = AE_err + 1/num_test_cases*...
            sum(sum((curr_data - recon_data).^2)); 
    
    end
    test_err(curr_iter)=AE_err/num_test_batches;
 
    % update user
    fprintf('Before iteration %d train squared error: %6.3f, test squared error: %6.3f.\n',...
        curr_iter,train_err(curr_iter),test_err(curr_iter));
    
    
    
     ulabel=unique(label);
 for i=1:length(label)
     [r,c]=find(label==label(i));
     index{i}=c;
     num_cc(1,i)=length(index{i});
 end
     for i=1:length(label)
     mean_class(:,i)=sum(train_layer(:, index{i}),2)/length(index{i});
     end
    mean_all=sum(train_layer,2)/size(train_layer,2); 
     
     for i=1:length(ulabel)
         [r,c]=find(label==ulabel(i));
         num_dd(i)=length(c);
         mean_dd(:,i)=sum(train_layer(:,c),2)/num_dd(i);
     end
     
     t1=(ones(size(train_layer,1),size(train_layer,2))-1./repmat(num_cc,size(train_layer,1),1));
     t2=train_layer-mean_class;
     t3=1./repmat(num_cc,size(train_layer,1),1)-1./(sum(num_cc)*ones(size(train_layer,1),size(train_layer,2)));
     t4=(mean_class-repmat(mean_all,1,length(label)));
     lda=t1.*t2-t3.*t4;
    
    % do conjugate gradient with NUM_LS linesearches in batches
    for curr_batch = 1:num_train_batches

        % get current batch
        curr_data = train_data(:,(1+(curr_batch-1)*num_train_cases): ...
            (curr_batch*num_train_cases));
         lda_data= lda(:,(1+(curr_batch-1)*num_train_cases): ...
             (curr_batch*num_train_cases));
           out_data= output(:,(1+(curr_batch-1)*num_train_cases): ...
             (curr_batch*num_train_cases));
  
        % conjugate gradient
        fprintf('Batch %d:\n',curr_batch);
        
        switch err_type
            case 1
                Wb = minimize (Wb,'cross_entropy',NUM_LS,n,curr_data,out_data,lambda, lda_data);
            case 2
                Wb = minimize (Wb,'recon_err',NUM_LS,n,curr_data,out_data,lambda, lda_data);
            otherwise
                error ('optimizing unknown error type!');
        end
            
        
    end
    
end



