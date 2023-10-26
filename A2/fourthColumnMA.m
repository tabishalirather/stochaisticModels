% Tabish Ali Rather-Assignment 2
data = readmatrix("A2\F17.txt");
num_columns = width(data);
% possible types for first three columns: Rando walk, AR(1), MA(1).
first_col = data(:,1);
second_col = data(:,2);
third_col = data(:,3);
ma_data = data(:,4);
len_ma = length(ma_data);
% x_n = c + a*eta_(n-1) + eta_n , 
% Estimate the value of c. 
% Optimal value of a, use conditonal least squares, write diffrent values
% Evaluate sum of squares of residuals S(a) =  sigma(e_i)^2. for any fixed
% value of a. S(a) = sigma(e_i)^2 versus a to see if the graph is minimum.
% Use min() fxn to find minimum of graph and determine the optimal value of
% a. 
% use estimate and arima to estimate the values of c and a. Construct 95%
% predicition interval for 20 steps ahead. 
% c is just mean of the data.
mean_MA = mean(ma_data);
% x = 'fourth_column'
% e_i = x_i - xhat_i = eta_i
% e_i = x_i - mu-a_eta_(i-1)
% let eta_1 = 0.
% a = 0.2; % For example
resi_ma = zeros(len_ma,1);
resi_ma_1(1,1) = 0;
resi_ma_2(1,1) = 0;
resi_ma_3(1,1) = 0;
% resi_ma(2,1) = ma_data(2,1) - mean_MA -a*resi_ma(1,1)
% resi_ma(3,1) = ma_data(3,1) - mean_MA -a*resi_ma(2,1)
% for a=linspace(-1,1,100)
%     for idx_res=2:len_ma
%     resi_ma(idx_res,1) = ma_data(idx_res,1) - mean_MA - a* resi_ma(idx_res-1,1);
%     end
% end
%     a = 0.2
%     for idx_res=2:len_ma
%         resi_ma_1(idx_res,1) = ma_data(idx_res,1) - mean_MA - a* resi_ma_1(idx_res-1,1);
%     end
%     a = 0.3
%     for idx_res=2:len_ma
%          resi_ma_2(idx_res,1) = ma_data(idx_res,1) - mean_MA - a* resi_ma_2(idx_res-1,1);
%     end
%     a = 0.4
%     for idx_res=2:len_ma
%         resi_ma_3(idx_res,1) = ma_data(idx_res,1) - mean_MA - a* resi_ma_3(idx_res-1,1);
%     end
% sum_sqrd_res_1 = sum(resi_ma_1.^2)
% sum_sqrd_res_2 = sum(resi_ma_2.^2)
% sum_sqrd_res_3 = sum(resi_ma_3.^2)
a_values = linspace(-1,1,2800);  % example range

% Initialize matrices to store the residuals
num_a_values = length(a_values);
residuals = zeros(len_ma, num_a_values);

% Compute residuals for each a value
for idx_a = 1:num_a_values
    a = a_values(idx_a);
    for idx_res = 2:len_ma
        residuals(idx_res, idx_a) = ma_data(idx_res, 1) - mean_MA - a * residuals(idx_res-1, idx_a);
    end
end
figure

% Compute the sum of squared residuals for each a value
sum_sqrd_res = sum(residuals.^2, 1);% sum along the rows
plot(a_values, sum_sqrd_res);
xlabel('a values')
ylabel('Residuals')
[min_value_resi_sqrd, min_index_a_resi_sqrd] = min(sum_sqrd_res);
% sum_sqrd_res(k) = minr


ma_arima = arima('Constant', NaN, 'MALags', 1, 'D',0);
est_ma_arima = estimate(ma_arima, ma_data);
c_estimated = est_ma_arima.Constant
a_estimated = est_ma_arima.MA{1}
optimal_a = a_values(min_index_a_resi_sqrd)

% ttest2(a_estimated, optimal_a, "Tail","right")

% [res, ~, ~] = infer(est_ma_arima, ma_data);
% res_optimal_a = residuals(:, min_index_a_resi_sqrd);
% [h, p] = ttest2(res, res_optimal_a, 'Tail', 'both')



