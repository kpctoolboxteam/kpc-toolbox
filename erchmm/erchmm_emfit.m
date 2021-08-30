function [D0, D1, logL] = map_emfit_erchmm(trace, orders, iter_max, iter_tol, verbose)

    %% define the conditions initially that we apply to EM-algorithm. Could put these inside function. 
    if nargin< 3
        iter_max = 300;
    end
    if nargin< 4
        iter_tol = 1e-7; % stop condition on the iterations
    end
        
    %% If we want to show the progress
    if nargin<5
        verbose = true;
    end

    %% Define a function that finds all the combinations of orders of Erlang branches from value of orders.
    
    function erlangs = all_erlang(branches, sum_erlangs)
        % Define a case for when we only have 1 branch. 
        if branches==1
            erlangs={sum_erlangs};
        else
            % initiate set of erlangs
            erlangs = {};
            % create for loop to create all combinations of erlang orders
            % for each branch. 
            for k1=1:sum_erlangs-branches+1
                x = all_erlang (branches-1, sum_erlangs-k1);
                for j=1:length(x)
                    sorted_erlangs = sort([x{j} k1]);
                    % check if we have the sorted erlang branch already
                    erlang_found=false;
                    for erl=1:length(erlangs)
                        if erlangs{erl}==sorted_erlangs
                            erlang_found=true;
                            break;
                        end
                    end
                    % Define condition fro when we add the erlang order to
                    % the set. 
                    if ~erlang_found
                        erlangs{length(erlangs)+1}=sorted_erlangs;
                    end
                end
            end
        end
    end

    %% Define a functiion to print the different Erlang numbers for when we run the algorithm.
    function print_erlang_number(num)
        for index=1:length(num)
            fprintf('%g',num(index));
            if index<length(num)
                fprintf(',');
            end
        end           
    end

%% Define process that runs the algorithm with lots of different orders. We then choose the algorithm with the best logli function.  
    
    % When we do the first run through we choose this method to give the
    % different combinations of orders from the sum defined by 'orders'. 
    
    if length(orders)==1
        % Set parameters that we will populate as empty.
        X_chosen = [];
        Y_chosen = [];
        orders_best = [];
        % Set log_li_best as -inf as initial value.
        log_li_best = -inf;
        for ord_num_1=2:orders
            all_orders = all_erlang(ord_num_1, orders);
            for ord_num_2=1:length(all_orders)
                % Print the progress of choosing orders so we know.
                if verbose
                    fprintf('Calculating with orders ');
                    print_erlang_number(all_orders{ord_num_2});
                    fprintf('...\n');
                end
                % Recursive part of EM-Algorithm. Call the function again
                % to see if new choice of orders is better than the
                % previous choice. 
                [ord_X,ord_Y,l] = EM_algorithm_using_ER_CHMM(trace, all_orders{ord_num_2});
                if l > log_li_best
                    X_chosen = ord_X;
                    Y_chosen = ord_Y;
                    log_li_best = l;
                    orders_best = all_orders{ord_num_2};
                end
            end
        end
        % Define D0, D1 and print the best log-likelihood
        D0 = X_chosen;
        D1 = Y_chosen;
        logL = log_li_best;
        % print results to the terminal
        if verbose
            fprintf('Best solution: log-likelihood=%g, orders=', log_li_best);
            print_erlang_number(orders_best);
            fprintf('\n');
        end
        return;
    end
    
    %% Define a function to find the matrix D0 from the Erlang Representation
    function X = generate_D0_from_erlangs(lambda_init, orders_init)
        X = zeros(sum(orders_init));
        init_x = 1;
        for lambda_i=1:length(lambda_init)
            X(init_x:init_x+orders_init(lambda_i)-1, init_x:init_x+orders_init(lambda_i)-1) = lambda_init(lambda_i)*(diag(ones(1,orders_init(lambda_i)-1),1)-diag(ones(1,orders_init(lambda_i))));
            init_x = init_x + orders_init(lambda_i);
        end
    end
    
    % Define lengths of orders and length of trace.
    M = length(orders);
    K = length(trace);

    % initialize pi and lambda parameters such that the mean is matched
    pi_v = ones(1,M) / M;
    lambda = diag(orders) * (1:M)';
    trace_mean = sum(trace)/length(trace);
    pi_mean = sum(pi_v ./ (1:M));
    % Match the mean
    lambda = lambda * pi_mean / trace_mean;
    % initialize the transition vector
    T = ones(M,1)*pi_v;
    
    % Initialize the parameters that we will use in the main EM-algorithm
    F = zeros(M, K); % Matrix for Densities
    A_likelihoods = zeros(K, M); % Matrix for forward likelihoods
    B_likelihoods = zeros(M, K); % Matrix for backward likelihoods
    A_likelihoods_scale = zeros(1,K); % Vector to scale f. likelihoods.
    B_likelihoods_scale = zeros(1,K); % Vector to scale b. likelihoods. 
    ologli = 1;
    logL = 0;
    % initialise number of steps and the time so we know both when we run
    % the algorithm. 
    steps = 1;
    t1 = clock();
    %% Start of the Main Algorithm
    while abs((ologli-logL)/logL)>iter_tol && steps<=iter_max
        ologli = logL;
        % E-step:
        %% Compute the branch densities:
        for i=1:M
            F(i,:) = ((lambda(i)*trace).^(orders(i)-1) / factorial(orders(i)-1) * lambda(i)) .* exp(-lambda(i)*trace);
        end
        %% Compute the forward likelihood vectors:
        % update the pi value and initialize scale value
        prev_pi = pi_v;
        scaled_prev = 0;
        for k=1:K
            prev_pi = prev_pi*diag(F(:,k))*T;
            scale = log2(sum(prev_pi));
            prev_pi = prev_pi * 2^-scale;
            A_likelihoods_scale(k) = scaled_prev + scale;
            A_likelihoods(k,:) = prev_pi;
            scaled_prev = A_likelihoods_scale(k);
        end
        a_forward_likelihoods = [pi_v; A_likelihoods(1:end-1,:)];
        A_scaled_v = [0, A_likelihoods_scale(1:end-1)];
        %% compute the backward likelihood vectors:
        next_b = ones(M,1);
        scaled_prev = 0;
        for k=K:-1:1
            next_b = diag(F(:,k))*T*next_b;
            scale = log2(sum(next_b));
            next_b = next_b * 2^-scale;
            B_likelihoods_scale(k) = scaled_prev + scale;
            B_likelihoods(:,k) = next_b;
            scaled_prev = B_likelihoods_scale(k);
        end
        b_backward_likelihoods = [B_likelihoods(:,2:end), ones(M,1)];
        B_scale_v = [B_likelihoods_scale(2:end), 0];
        
        likelihood_value = pi_v*B_likelihoods(:,1);
        logL = (log(likelihood_value) + B_likelihoods_scale(1) * log(2)) / K;
        i_likelihood = 1.0 / likelihood_value;

        %% M-step:
        
        % Calculate new estimates for the parameters lambda and transition
        % matrix T.
        
        likelihoods_multiplied = a_forward_likelihoods.*B_likelihoods';
        summed_lm = sum(likelihoods_multiplied,2);
        % Multiply the a[k] and b[k] vectors for the lambda and T
        % estimations. 
        for m=1:M
            likelihoods_multiplied(:,m) = likelihoods_multiplied(:,m)./summed_lm;
        end
        % Compute the numerator and the denominator for the lambda
        % estimations as well as p_{i,j} estimations
        
        numerator_estimation = sum(likelihoods_multiplied,1);
        denominator_estimation = trace'*likelihoods_multiplied;
        pi_v = numerator_estimation/K;
        lambda = (orders.*numerator_estimation ./ denominator_estimation)';
        
        % Compute multiplication of forward likelihood and densities
        dens_mult_a_likelihood = a_forward_likelihoods.*F';
        % Find summed lambda values for estimation of T.
        summed_lm = i_likelihood*2.^(A_scaled_v+B_scale_v-B_likelihoods_scale(1))';
        for m=1:M
            dens_mult_a_likelihood(:,m) = dens_mult_a_likelihood(:,m) .* summed_lm;
        end       
        % Compute the matrix T.
        T = (dens_mult_a_likelihood'*b_backward_likelihoods').*T;
        % Normalize it
        for m=1:M
            T(m,:) = T(m,:) / sum(T(m,:));
        end        
        
        % Show increase in step if you want
        steps = steps+1;
        % Print progress report
        if verbose && etime(clock(),t1)>2
            fprintf('Num of iterations: %d, log-likelihood: %g\n', steps, logL);
            t1 = clock();
        end
    end
    
    %% Show Progress if needed
    if verbose
        fprintf('Num of iterations: %d, log-likelihood: %g\n', steps, logL);
        fprintf('EM algorithm terminated. (orders=');
        for i=1:M
            fprintf('%g',orders(i));
            if i<M
                fprintf(',');
            else
                fprintf(')\n');
            end
        end
    end
    
    %% Finalize formulation of D0 and D1 matrices
    
    % Generate the mstrix D0 from lambda and the orders
    D0 = generate_D0_from_erlangs(lambda, orders);
    % initialize D1 size
    D1 = zeros(size(D0));
    % Find indices that we will populate for D1 estimation
    indicesTo = [1 cumsum(orders(1:end-1))+1];
    indicesFrom = cumsum(orders);
    % Use diagonal of lambda and the T matrix values
    D1(indicesFrom,indicesTo) = diag(lambda)*T;
end