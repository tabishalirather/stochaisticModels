% Tabish Ali Rather-Assignment 2
data = readmatrix("A2\F17.txt");
num_columns = width(data)
% possible types for first three columns: Rando walk, AR(1), MA(1).
first_col = data(:,1);
second_col = data(:,2);
third_col = data(:,3);
ma = data(:,4);

% plotting ACF's
plot_diagnostics(1, first_col, second_col, third_col)

% first_col
first_col_analysis(first_col)

% second_col
second_col_analysis(second_col)
third_col_analysis(third_col)




function third_col_analysis(third_col)

    % From the PACF and ACF of third column, ACF is non zero for many lags,
    % and PACF is zero after lag one. It is indicative of a AR(1) time
    % series. Let's use AR(1) modelling to forecast the data.
    % The standard model is: y_n = mu + ay_(n-1 - mu) + neeeta.
    % Let's being by estimating the paramters. 
    mu = mean(third_col)
    % Using method of moments, let's estimate a as well. as is equal to
    % first value of sample acf at lag 1 which is stored as second element
    % of autocorr vector.
    acf = autocorr(third_col);
    a = acf(2)
    % Now, let's calculate the residuals. resi = observer - fitted. 
    % Fitted_i = mu+a(x_(i-1) - mu)
    % resi = x_i - (mu + a(x_(i-1)-mu)
    len_third_col = length(third_col);
    resi_third_col = zeros(1,len_third_col);
    resi_third_col(1,1) = 0
    for idx=2:len_third_col
        resi_third_col(idx) = third_col(idx) - mu-a*(third_col(idx-1)-mu);
    end

    % Model validation: lbqtest, swtest, ttest
    lags = floor(log(len_third_col))
    [h_lbq_third_col, p_lbq_third_col] = lbqtest(resi_third_col, 'Lags', lags, 'DOF', lags-2)
    % Null hypothesis: ACF is not different from zero, there are not
    % correlations
    % Alternative hypothesis: ACF is different from zero, there are
    % correlations.
    % p < 0, Reject null hypothesis: ACF is different from zero, i.e. residuals
    % do not behave like white noise. 
    % p > o, Don't reject null hypothesis, ACF is not different from zero, i.e.
    % residuals are behaving like white noise.
    % p_scnd_lbq = 0.1837 > 0.05, ACF upto lags k is not different from
    % zero, we conclude, residuals behave like white noise.
    
    % swtest
    [h_sw_third_col, p_sw_third_col] = swtest(resi_third_col)
    %  p < 0.05 the distribution of the residuals is significantly different from normal.
    %  p > 0.05 the distribution of the residuals is NOT significantly different from normal.
    % Here p = 0.8407 > 0.05,, so residuals distribution is not significantly
    % different from normal.
    
    % ttest
    [h_third_ttest, p_third_ttest] = ttest(resi_third_col)
       % p = 0.6986,

       % Let's calulate forecast, 95% confidence interval and plot the
       % fitted values.
       fitted_vals_third_col = zeros(len_third_col,1);
       % fitted_vals_third_col(1,1) = 0

       for idx_fitted=2:len_third_col
           fitted_vals_third_col(idx_fitted) = mu+a*(third_col(idx_fitted-1)-mu);
       end
        % let's calculate the forecast values now for 20 steps.. 
        p_steps = 20;
        % forecast values = mu + (a)^p*(x_n-mu)
        for idx_forecast=1:p_steps
            forecast_values_p_steps_third(idx_forecast) = mu + ((a)^p_steps)*(third_col(idx_forecast)-mu);  
            lower_bound_third_col(idx_forecast) = mu + (a^idx_forecast)*(third_col(len_third_col) - mu) - 1.96*(std(resi_third_col))*sqrt(((1-a^(2*idx_forecast))/(1-a^2)));
            upper_bound_third_col(idx_forecast) = mu + (a^idx_forecast)*(third_col(len_third_col) - mu) + 1.96*(std(resi_third_col))*sqrt(((1-a^(2*idx_forecast))/(1-a^2)));
        end
        
       
       
       h = figure(999);
    set(h, 'Position', [300, 60, 1160, 700])
    hold on
    time_series_data = plot(third_col, 'k--', LineWidth=1)
    forecast_values_plot = plot([len_third_col+1:len_third_col+p_steps], forecast_values_p_steps_third, 'b-', LineWidth=2)
    lower_bound = plot([len_third_col+1:len_third_col+p_steps], lower_bound_third_col, 'r--');
    upper_bound = plot([len_third_col+1:len_third_col+p_steps], upper_bound_third_col, 'r--');
    plot_fittd_first_col = plot(fitted_vals_third_col, 'g--', LineWidth=1)
    legend([time_series_data forecast_values_plot lower_bound upper_bound plot_fittd_first_col], 'Observed', 'Forecast', 'Lower 95% interval', 'Upper 95% interval', 'Fitted values', 'Location','southWest')

end

function second_col_analysis(second_col)
x_test = 'calling second here'
    %% 
    % model: Naive model: x_n = x_(i-1) + n_i = 1/1-a. a = ACF of raw
% data at lag one
% acf = autocorr(second_col);
% % acf at lag 1 is stored as the second element in the acf vector 
% a = acf(2)
% % Second column. 
% ACF is gradually decreasing, non zero for many lags and PACF cuts off.
mu = mean(second_col);
std_dev= std(second_col);
len_second_col = length(second_col);
resi_second_col = zeros(len_second_col,1);
resi_second_col(1,1) = 0;
% Calculating residuals.
for j=2:len_second_col
    resi_second_col(j) = second_col(j) - second_col(j-1)
end
for idx_fit = 2:len_second_col
    fitted_scnd_col(idx_fit) =  second_col(idx_fit-1)
end



% Validating our model.
% lbqtest
lags = floor(log(len_second_col));
% dof = lags - 2 as we estimated two paramets, mu and c.
[h_scnd_lbq, p_scnd_lbq] = lbqtest(resi_second_col, 'Lags', lags, 'DOF', lags - 2)
% Null hypothesis: ACF is not different from zero, there are not
% correlations
% Alternative hypothesis: ACF is different from zero, there are
% correlations.
% p < 0, Reject null hypothesis: ACF is different from zero, i.e. residuals
% do not behave like white noise. 
% p > o, Don't reject null hypothesis, ACF is not different from zero, i.e.
% residuals are behaving like white noise.
% p_scnd_lbq = 0.1861 > 0.05, ACF upto lags k is not different from zero. 

[h_scnd_sw, p_scnd_sw] = swtest(resi_second_col)
%  p < 0.05 the distribution of the residuals is significantly different from normal.
%  p > 0.05 the distribution of the residuals is NOT significantly different from normal.
% Here p = 0.8480 > 0.05,, so residuals distribution is not significantly
% different from normal.
% test for means. 
[h_scnd_ttest, p_scnd_ttest] = ttest(resi_second_col)
% p = 0.9107,  so we fail to reject the null hypothesis, i.e there is not
% enough evidence to say the mean is statistically is different from zero.




% Now, let's get the forecast values for 20 steps ahead with 95% prediciton intervals..
p_steps = 20;
for k = 1:p_steps
    lower_bound_scnd_col(k) = second_col(k) - 1.96*sqrt(p_steps)*std_dev
    upper_bound_scnd_col(k) = second_col(k) + 1.96*sqrt(p_steps)*std_dev
    forecast_values(k) = second_col(len_second_col)
end
% x = 'this is printng'
% figure('Number', 99, 'Position', [300, 60, 1160, 700]); 
% figure('Number', 99, 'Position', [300, 60, 1160, 700], 'Name', 'Your Desired Title', 'NumberTitle', 'off');
h = figure(99); 
set(h, 'Position', [300, 60, 1160, 700])
hold on
plot(second_col, 'k--', LineWidth=1)
plot([len_second_col+1:len_second_col+p_steps], forecast_values, 'b-', LineWidth=2)
plot(fitted_scnd_col, 'g--', LineWidth=1)
plot([len_second_col+1:len_second_col+p_steps], lower_bound_scnd_col, 'r--', [len_second_col+1:len_second_col+p_steps], upper_bound_scnd_col, 'r--', LineWidth=1);

% legend('Observed','Forecast','Fitted values','95% Confidence Interval', 'Location','NorthWest');
hold off
% plot(fitted_scnd_col)
% plot(second_col, 'k--', LineWidth=1)



mean_resi_scnd = mean(resi_second_col);
std_resi_scnd= std(resi_second_col);

% figure
% histogram(resi_second_col, 'Normalization', 'pdf')
% hold on
% x = linspace(min(resi_second_col), max(resi_second_col), 201);
% y = normpdf(x, mean_resi_scnd, std_resi_scnd);
% plot(x, y, 'r-', LineWidth=2)
% xlabel('Data');
% ylabel('PDF');
% legend('Histogram with pdf of Normal Distribution');
% hold off

end



% function second_col_analysis(second_col)
% x_test = 'calling second here'
%     %% 
%     % model: y_n = mu + c(y_(n-1) - mu) + eta_n, mu = 1/1-a. a = ACF of raw
% % data at lag one
% acf = autocorr(second_col);
% % acf at lag 1 is stored as the second element in the acf vector 
% a = acf(2)
% % Second column. 
% % ACF is gradually decreasing, non zero for many lags and PACF cuts off.
% mu = mean(second_col);
% len_second_col = length(second_col);
% resi_second_col = zeros(len_second_col,1);
% for j=2:len_second_col
%     resi_second_col(j) = second_col(j) - mu - a*(second_col(j-1) - mu);
% end
% for idx_fit = 2:len_second_col
%     fitted_scnd_col(idx_fit) =  mu + a*(second_col(idx_fit-1) - mu)
% end
% % Now, let's get the forecast values for 20 steps ahead with 95% prediciton intervals..
% p_steps = 20;
% for k = 1:p_steps
%     lower_bound_scnd_col(k) = mu + (a^k)*(second_col(len_second_col) - mu) - 1.96*(std(resi_second_col))*sqrt(((1-a^(2*k))/(1-a^2)));
%     upper_bound_scnd_col(k) = mu + (a^k)*(second_col(len_second_col) - mu) + 1.96*(std(resi_second_col))*sqrt(((1-a^(2*k))/(1-a^2)));
%     forecast_values(k) = mu + a^k*(second_col(len_second_col) - mu);
% end
% x = 'this is printng'
% figure('Position', [300, 60, 1160, 700]); 
% hold on
% plot(second_col, 'ko--', LineWidth=1)
% plot([len_second_col+1:len_second_col+p_steps], forecast_values, 'b', LineWidth=2)
% plot(fitted_scnd_col, 'g--', LineWidth=1)
% plot([len_second_col+1:len_second_col+p_steps], lower_bound_scnd_col, 'r--', [len_second_col+1:len_second_col+p_steps], upper_bound_scnd_col, 'r--', LineWidth=1);
% 
% legend('Observed','Forecast','Fitted values','95% Confidence Interval', 'Location','NorthWest');
% hold off
% % plot(fitted_scnd_col)
% 
% % Validating our model.
% % lbqtest
% lags = floor(log(len_second_col));
% % dof = lags - 2 as we estimated two paramets, mu and c.
% [h_scnd_lbq, p_scnd_lbq] = lbqtest(resi_second_col, 'Lags', lags, 'DOF', lags - 2);
% % Null hypothesis: ACF is not different from zero, there are not
% % correlations
% % Alternative hypothesis: ACF is different from zero, there are
% % correlations.
% % p < 0, Reject null hypothesis: ACF is different from zero, i.e. residuals
% % do not behave like white noise. 
% % p > o, Don't reject null hypothesis, ACF is not different from zero, i.e.
% % residuals are behaving like white noise.
% % p_scnd_lbq = 0.1384 > 0.05, ACF upto lags k is not different from zero. 
% 
% [h_scnd_sw, p_scnd_sw] = swtest(resi_second_col)
% %  p < 0.05 the distribution of the residuals is significantly different from normal.
% %  p > 0.05 the distribution of the residuals is NOT significantly different from normal.
% % Here p = 0.9478 > 0.05,, so residuals distribution is not significantly
% % different from normal.
% % test for means. 
% [h_scnd_ttest, p_scnd_ttest] = ttest(resi_second_col)
% % p > 0.7712, 
% 
% mean_resi_scnd = mean(resi_second_col);
% std_resi_scnd= std(resi_second_col)
% 
% % figure
% % histogram(resi_second_col, 'Normalization', 'pdf')
% % hold on
% % x = linspace(min(resi_second_col), max(resi_second_col), 201);
% % y = normpdf(x, mean_resi_scnd, std_resi_scnd);
% % plot(x, y, 'r-', LineWidth=2)
% % xlabel('Data');
% % ylabel('PDF');
% % legend('Histogram with pdf of Normal Distribution');
% % hold off
% 
% end




function first_col_analysis(first_col)
% lbqtest
len_first_col = length(first_col);
lags = floor(log(len_first_col));
% [h_lbq, p_lbq] = lbqtest(first_col, 'lags', lags, 'DOF', lags-1)
% p_lbq = 0.3026 > 0.05, we cannot reject the null hypothesis, the first
% two sample ACF's are not different than zero. Let's proceed to model this
% data using the mean model: 
num_points = 201;
mean_first_col = mean(first_col);
resi = zeros(len_first_col,1)
for i=1:len_first_col
    resi(i) = first_col(i) - mean_first_col;
end
p_steps = 20;
forecast_value = mean_first_col;

% Now, let's create 97.5% prediction intervals. 
std_dev = std(resi);
t_score = icdf('T', 0.975, len_first_col-1);

forecast_values_p_steps = zeros(p_steps,1);
error = zeros(p_steps,1);
for i=1:p_steps
    forecast_values_p_steps(i) = mean_first_col;
    error(i) = t_score*std_dev*sqrt(1+1/len_first_col)
end
for idx = 1:length(first_col)
    fitted_vals_first_col(idx) = mean_first_col
end
% plot time series, forecast values and prediciton intervals. 
figure('Position', [300, 60, 760, 700])
time_series_data = plot(first_col, 'o-')
% x = 'this is printing'
hold on
forecast_values = plot(len_first_col+1:len_first_col+20,forecast_values_p_steps, 'b', LineWidth=1);
upper_bound = plot(len_first_col+1:len_first_col+20,forecast_values_p_steps+error, 'r:', LineWidth=2);
lower_bound = plot(len_first_col+1:len_first_col+20,forecast_values_p_steps-error, 'r:', LineWidth=1);
plot_fittd_first_col = plot(fitted_vals_first_col, 'g', LineWidth=1)
legend([time_series_data forecast_values upper_bound lower_bound plot_fittd_first_col], 'Observed', 'Forecast', 'Upper 95% interval', 'Lower 95% interval', 'Fitted values')
hold off

% Model validation:
[h_lbq, p_lbq_res] = lbqtest(resi, 'lags', lags, 'DOF', lags-1);
% p = 0.3026 > 0.05, cannot reject null hypothesis, residuals are likely
% not different from white noise process. 
[h_sw, p_sw] = swtest(resi)
% p > 0.852, cannot reject the null hypothesis. So, the distribution of
% residuals is likely not different form normal.
% Now, let's fit histogram to residuals.
%% 
%% 

mean_resi = mean(resi);
std_resi = std(resi)  
figure
histogram(resi, 'Normalization', 'pdf')
hold on
x = linspace(min(resi), max(resi), 201);
y = normpdf(x, mean_resi, std_resi);
plot(x, y, 'r-', LineWidth=2)
xlabel('Data');
ylabel('PDF');
legend('Histogram with pdf of Normal Distribution');
hold off
end

% after first lag. So it is characteristic of AR(1) process. 
% let's start modelling using AR(1). 
function plot_diagnostics(i, first_col, second_col, third_col)
series_names = {'first_`col', 'second_`col', 'third_`col'};
seriesArray = {first_col, second_col, third_col};
for i=1:length(seriesArray)
    series = seriesArray{i};
    figure('Position', [300, 60, 1160, 700]); % Set the position and size of the figure [left, bottom, width, height]
    
    % Time Series Plot
    subplot(3, 1, 1);  % 3 rows (for time series, ACF, PACF), 1 column, 1st plot
    plot(series);
    title(['Time Series Plot for Series ' series_names{i}], 'FontSize',15);

    % ACF Plot
    subplot(3, 1, 2);  % 3 rows, 1 column, 2nd plot
    autocorr(series,floor(length(series)/4));
    title(['ACF Plot for Series ' series_names{i}], 'FontSize',15);

    % PACF Plot
    subplot(3, 1, 3);  % 3 rows, 1 column, 3rd plot
    parcorr(series, floor(length(series)/4));
    title(['PACF Plot for Series ' series_names{i}], 'FontSize',15);
end
end
