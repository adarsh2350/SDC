clc
clear


%% Data Accessing

% Accessing data from SDC_Data
[dist, time, t_freq, cst_C] = SDC_Data;   

% Accessing data from SDC_input
[n, K, T, D, T_avg, opr_chrgs, base_fare, reg_chrg, reg_lim, add_chrg, tax, profit_factor, dis_factor, DS, CS] = SDC_input; 



%% Determination of maximum number cabs of each type

t_total = zeros(1,K);

for k = 1:K
    for i = 1:n
        for j = 1:n
            t_total(1,k) = t_total(1,k) + (t_freq(i,j,k)*time(i,j));
        end
    end
end

max_cab = zeros(1,K);  % maximum number of cabs of paricular type

for i = 1:K
    max_cab(1,i) = 1.1*(t_total(1,i)/600);
    max_cab(1,i) = round(max_cab(1,i));
end

%%%  Cost determination

%% lower limit of pricing
 
trav_factor = zeros(1,K);
for k = 1:K
    trav_factor(k) = (opr_chrgs(k) / (30*T));  % operational cost factor
end

ll_cost = zeros(n,n,K); % lower limit

for k = 1:K
    for i = 1:n
        for j = 1:n
            if (i==j)
                continue;
            end
            temp = base_fare(1,k);
            temp = temp + (trav_factor(1,k)*time(i,j));
            
            % Amt = 
            if dist(i,j)<=reg_lim
                temp = temp + (reg_chrg(1,k)*dist(i,j));
            else
                temp = temp + (reg_chrg(1,k)*reg_lim);
                temp = temp + (add_chrg(1,k)*(dist(i,j) - 15));
            end

            ll_cost(i,j,k) = profit_factor.*temp;
            ll_cost(i,j,k) = round(ll_cost(i,j,k));
        end
    end
end

%% upper limit of pricing

ul_cost = zeros(n,n,K); % upper limit

for k = 1:K
    for i = 1:n
        for j = 1:n
            ul_cost(i,j,k) = dis_factor*cst_C(i,j,k);
        end
    end
end



%% forming horizontal matrix

temp1 = ll_cost;
h_ll = zeros(K,n*n);
for i = 1:K
  h_ll(i,:) = reshape(temp1(:,:,i),1,[]);  
end

temp2 = ul_cost;
h_ul = zeros(K,n*n);
for i = 1:K
  h_ul(i,:) = reshape(temp2(:,:,i),1,[]);  
end

temp3 = t_freq;
htf = zeros(K,n*n);
for i = 1:K
  htf(i,:) = reshape(temp3(:,:,i),1,[]);  
end



%%% Determining new cost via sTLBO

new_cost = zeros(n,n,K);  % new pricing



for q = 1:2

    Q = 10;  % number of Runs
    fun_eval = 1:Q;
    fitness_value = zeros(1,Q);
    least_val = 0;
    best_fitness_value = zeros(1,Q);


    for iteration = 1:Q
    
        %% problwm setting
        prob = @fun;   % fitness function
        lb = h_ll(q,:);     % lower bound
        ub = h_ul(q,:);   % upper bound
        %% Algorithm parameters
        Np = 10;        % Population size
        T = 50;         % No. of iterations
        %% starting of TLBO
        f = NaN(Np,1);   % vector to store fitness value
        D = length(lb);  % determining no. of decision variables
        P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);   % initial population
        
        for i = 1:Np              % evaluating fitnessfunction
           f(i) = prob(P(i,:),lb,ub,htf(q,:));
        end
        
        %% iteration loop
        
        plotx = 0:T;
        ploty = zeros(1,T+1);
        ploty(1) = min(f);
        
        for t = 1:T

            for i=1:Np
          
               % Teacher phase
               Xmean = mean(P);
              
               [~,ind] = min(f);
               Xbest = P(ind,:);
               TF = randi([1,2],1,1);
               Xnew = P(i,:) + rand(1,D).*(Xbest - TF*Xmean);
               Xnew = min(ub,Xnew);
               Xnew = max(lb,Xnew);
              
               fnew = fun(Xnew,lb,ub,htf(q,:));
               if fnew < f(i)
                   P(i,:) = Xnew;
                   f(i) = fnew;
               end
          
               % Learner phase
               p = randi([1 Np],1,1);  % selection of random partner
               while i==p
                   p = randi([1 Np],1,1);
               end
               if f(i) < f(p)
                   Xnew = P(i,:) + rand(1,D).*(P(i,:) - P(p,:));
               else
                   Xnew = P(i,:) - rand(1,D).*(P(i,:) - P(p,:));
               end
               Xnew = min(ub,Xnew);
               Xnew = max(lb,Xnew);
               fnew = fun(Xnew,lb,ub,htf(q,:));
               if fnew < f(i)
                   P(i,:) = Xnew;
                   f(i) = fnew;
               end
            end

            ploty(t+1) = min(f);
        end
        [bestfitness,ind] = min(f);
        bestsol = P(ind,:);
        fitness_value(iteration) = min(f);

        if min(f) < least_val
            least_val = min(f);
            temp = reshape(bestsol,n,[]);
            temp = transpose(temp);
            new_cost(:,:,q) = temp;
        end
        best_fitness_value = least_val;
        
        
        %% plotting best fitness vs iteration

          scatter(plotx,ploty)
          xlabel('Iterations ->')
          ylabel('Best fitness ->')
          title('best fitness vs iterations')
          ploty = [];

        

    end 

    %% plotting best fitness vs iteration
        
       scatter(fun_eval,fitness_value)
       xlabel('function evaluation ->')
       ylabel('fitness value ->')
       title('fitness value vs fitness evaluation')
       fitness_value = [];
        
    %% plotting best fitness vs iteration

      scatter(fun_eval,best_fitness_value)
      xlabel('function evaluation ->')
      ylabel('Best fitness value ->')
      title('best fitness value vs function evaluation')
      best_fitness_value = [];

end



%% removal of loss making rides

for k = 1:K
    for i = 1:n
        for j= 1:n
            if (ll_cost(i,j,k) > ul_cost(i,j,k))
                new_cost(i,j,k) = -1;
            end
            new_cost(i,j,k) = round(new_cost(i,j,k));
        end
    end
end   


%% Determination of number of cabs in a particular zone

Cabs = zeros(n,n,K);

for k = 1:K

    Temp = zeros(n,n);
    for i = 1:n
        for j = 1:n
            Temp(i,j) = (300/(time(i,j)+time(j,i)))*(new_cost(i,j)+new_cost(j,i));
        end
    end
    
    for i = 1:n
        for j = 1:n
            if i==j
                Temp(i,j) = 0;
            end
        end
    end

    Z = reshape(Temp,1,[]);
    Z = -1.*Z;


    ilb = zeros(1,n*n);
    iub = zeros(1,n*n);
    TEMP = zeros(n,n);
    for i = 1:n
        for j = 1:n
            TEMP(i,j) = (time(i,j)*t_freq(i,j,k))/600;
            TEMP(i,j) = round(TEMP(i,j));
        end
    end

    temp4 = TEMP;
    iub = reshape(TEMP,1,[]);  
    
    A = ones(1,n*n);
    b = max_cab(1,k);
    
    intcon = 1:n;

    [X, FVAL] = intlinprog(Z,intcon,A,b,[],[],ilb,iub);
    Cabs(:,:,k) = reshape(X,n,[]);
    Cabs(:,:,k) = transpose(Cabs(:,:,k));
end    
% ll_cost 
% ul_cost
% new_cost
new_cost;
Cabs;







