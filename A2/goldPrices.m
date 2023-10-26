% % gold_price = readtable("")
% range = 'E11187:F11287';
% data = readtable('Historic-Gold-Prices.xlsx', 'Range', range)
% date = data{:,1}
% price = data{:,2}
% days = 1:101
% range_next_20_days = 'E11288:F11308';
% data_next_20_days = readtable('Historic-Gold-Prices.xlsx', 'Range', range_next_20_days);
% price_next_20_days = data_next_20_days(:,2);
% figure(1111)
% plot(days, price);
% % 
% % linear_coeff = polyfit(days, price, 1);
% % quad_fit = polyfit(days, price, 2);
% % plot_diagnostics(1, price)
% 
% plot_forecast_values_vs_actual(price_next_20_days, days,price)
% % plot_fit_vs_actual(price_next_20_days,days, price)
% arima_fxn(1,1,1,days,price)
% arima_fxn(2,1,2,days,price)
% % arima_fxn p_lbq is better for arima(1,1,1) comapred to arima(1,2,1), so
% % we assume arima(1,1,1) model.
% 
% % y_forecast_20_days = arima_fxn(2,1,3,days,price, price_next_20_days);
% % forecast_price_20_days = y_forecast_20_days(:,1);
% % actual_price_20_days = price_next_20_days.Var2;
% % % difference = forecast_price_20_days - actual_price_20_days
% %  % figure(5000)
% %  days_next_20 = (days(end)+1):(days(end)+21);
% %  % diff_data_actual = zeros(21,1)
% %  % diff_data_actual = diff(actual_price_20_days);
% %  % diff_data_actual(21,1) = 0
% % 
% %  % diff_data_forecast = diff(forecast_price_20_days)
% % 
% % figure('Position', [300, 60, 1160, 700]); % Set the position and size of the figure [left, bottom, width, height]
% % 
% % % Plot the original data
% % % plot(length(price), price, 'b-', 'LineWidth', 1);
% % plot(days_next_20, forecast_price_20_days, 'r--', 'LineWidth', 2);
% % hold on;
% % 
% % % Plot the forecasted prices for the next 20 days
% % % plot(days_next_20, forecast_price_20_days, 'r--', 'LineWidth', 2);
% % 
% % % Plot the actual prices for the next 20 days
% % plot(days_next_20, actual_price_20_days, 'g-', 'LineWidth', 2);
% % 
% % legend('Forecasted Price', 'Actual Price');
% % title('ARIMA Forecast vs Actual Prices for the Next 20 Days');
% % xlabel('Days');
% % ylabel('Gold Price');
% % 
% % hold off;
% 
% % function plot_fit_vs_actual(price_next_20_days,days, price)
% % 
% %     [y_forecast_20_days, fitted_data] = arima_fxn(1,1,1,days,price, price_next_20_days);
% %     fitted_data = fitted_data(:,1);
% %     len_price = [1:length(price)]
% %     figure('Position', [300, 60, 1160, 700]) % Set the position and size of the figure [left, bottom, width, height]    plot(len_price, price)
% %     % figure(12313)
% %         plot(len_price, price)
% %     hold on;
% %         plot (len_price, fitted_data)
% %     hold off
% % end
% 
% 
% function plot_forecast_values_vs_actual(price_next_20_days,days, price)
%     y_forecast_20_days = arima_fxn(0,1,0,days,price, price_next_20_days);
% forecast_price_20_days = y_forecast_20_days(:,1);
% actual_price_20_days = price_next_20_days.Var2;
% % difference = forecast_price_20_days - actual_price_20_days
%  % figure(5000)
%  days_next_20 = (days(end)+1):(days(end)+21);
%  % diff_data_actual = zeros(21,1)
%  % diff_data_actual = diff(actual_price_20_da,ys);
%  % diff_data_actual(21,1) = 0
% 
%  % diff_data_forecast = diff(forecast_price_20_days)
% 
% figure('Position', [300, 60, 1160, 700]); % Set the position and size of the figure [left, bottom, width, height]
% 
% % Plot the original data
% % plot(length(price), price, 'b-', 'LineWidth', 1);
% plot(days_next_20, forecast_price_20_days, 'r--', 'LineWidth', 2);
% hold on;
% 
% % Plot the forecasted prices for the next 20 days
% % plot(days_next_20, forecast_price_20_days, 'r--', 'LineWidth', 2);
% 
% % Plot the actual prices for the next 20 days
% plot(days_next_20, actual_price_20_days, 'g-', 'LineWidth', 2);
% 
% legend('Forecasted Price', 'Actual Price');
% title('ARIMA Forecast vs Actual Prices for the Next 20 Days');
% xlabel('Days');
% ylabel('Gold Price');
% 
% hold off;
% 
% end
% 
% 
% 
% function y_forecast_20_days = arima_fxn(p,d,q,days,price, price_next_20_days)
%     arima_model = arima('Constant', NaN, 'MALags', [1:p] ,'ARLags', [1:q], 'D',d)
%     estimated_model = estimate(arima_model, price);
%     resi_arima = infer(estimated_model, price);
%     fitted_data = price + resi_arima;
%     len_price = [1:length(price)]
%     figure('Position', [300, 60, 1160, 700]) % Set the position and size of the figure [left, bottom, width, height]    plot(len_price, price)
%     % figure(12313)
%         plot(len_price, price)
%     hold on;
%         plot (len_price, fitted_data)
%     hold off
%     lags = floor(log(length(resi_arima)));
%     % Model validation. lbq
%     [h_arima_lbq, p_arima_lbq] = lbqtest(resi_arima, 'lags', lags, 'DOF', lags-2)
%     % swtest
%     [h_arima_sw, p_arima_sw] = swtest(resi_arima)
%     % ttest
%     [h_arima_ttest, p_arima_ttest] = ttest(resi_arima)
% 
%     %Forecast values
%     [y_forecast_20_days, y_mean_sqr_err] = forecast(estimated_model, 21, 'Y0', price);
%     y_forecast_20_days = y_forecast_20_days;
%     length(y_forecast_20_days)
% 
%     % upper = y_forecast_20_days + 1.96*sqrt(y_mean_sqr_err);
%     % lower = y_forecast_20_days - 1.96*sqrt(y_mean_sqr_err);
%     % 
%     % figure;
%     % plot(price, 'Color', [0.7, 0.7, 0.7]);
%     % hold on;
%     % h = plot(length(price)+1:length(price)+20, y_forecast_20_days, 'r', 'LineWidth', 2);
%     % legend('Historical Data', 'Forecast');
%     % title('20-day Forecast with 95% Confidence Interval');
%     % fill([length(price)+1:length(price)+20, length(price)+20:-1:length(price)+1], [upper; flipud(lower)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     % y_forecast_20_days = y_forecast_20_days{:};
%     % price_next_20_days = price_next_20_days{:};
% 
% end
% 
% 
% % from ACF and PACF, pacf gradually cuts of, very slowly and PACF spike for
% % 0 annd 1 lag and is zero after that, the best description would be naive
% % random walk data, and using naive model.
% function plot_diagnostics(i, data)
% diff_data = diff(data);
% series_names = {'Gold Price'};
% seriesArray = {diff_data};
% for i=1:length(seriesArray)
%     series = seriesArray{i};
%     figure('Position', [300, 60, 1160, 700]); % Set the position and size of the figure [left, bottom, width, height]
% 
%     % Time Series Plot
%     subplot(3, 1, 1);  % 3 rows (for time series, ACF, PACF), 1 column, 1st plot
%     plot(series);
%     title(['Time Series Plot for Series ' series_names{i}], 'FontSize',15);
% 
%     % ACF Plot
%     subplot(3, 1, 2);  % 3 rows, 1 column, 2nd plot
%     autocorr(series,floor(length(series)/4));
%     title(['ACF Plot for Series ' series_names{i}], 'FontSize',15);
% 
%     % PACF Plot
%     subplot(3, 1, 3);  % 3 rows, 1 column, 3rd plot
%     parcorr(series, floor(length(series)/4));
%     title(['PACF Plot for Series ' series_names{i}], 'FontSize',15);
% end
% end


% Reading data
range = 'E11187:F11287';
data = readtable('Historic-Gold-Prices.xlsx', 'Range', range);
date = data{:,1};
price = data{:,2};
days = 1:101;

range_next_20_days = 'E11288:F11308';
data_next_20_days = readtable('Historic-Gold-Prices.xlsx', 'Range', range_next_20_days);
price_next_20_days = data_next_20_days(:,2);

% figure(1111)
% plot(days, price);
days_next_20 = (days(end)+1):(days(end)+21);

% [y_forecast_20_days, fitted_data] = arima_fxn(2,1,2,days,price, price_next_20_days);
for p = 0:2
    for d = 1:2
        for q = 0:2
            [y_forecast, fitted_data, upper, lower] = arima_fxn(p,d,q,days,price, price_next_20_days);

            % Plots for fitted, forecast values vs actual on the same graph
            combined_plot(p, d, q, days, price, y_forecast, fitted_data, price_next_20_days, upper, lower);
        end
    end
end



function combined_plot(p, d, q, days, price, y_forecast, fitted_data, price_next_20_days, upper, lower)
    figure('Position', [300, 60, 1160, 700]);
    len_price = [1:length(price)];
    price_next_20_days = price_next_20_days.Var2;
    days_next_20 = (days(end)+1):(days(end)+21);

    % Plotting the original data
    plot(len_price, price, 'g-', 'LineWidth', 1);
    hold on;
    
    % Plotting the fitted values
    plot(len_price, fitted_data, 'r--', 'LineWidth', 2);
    
    % Plotting the forecast values
    plot(days_next_20, y_forecast, 'bl-.', 'LineWidth', 2);
    plot(days_next_20, price_next_20_days, 'g-.', 'LineWidth', 2);
    fill([days_next_20 fliplr(days_next_20)], [upper' fliplr(lower')], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Corrected order of legend entries to match the plotting order
    legend('Actual Data', 'Fitted Data', 'Forecasted Data', 'Actual next 20 days', '95% Confidence Interval');
    title(['ARIMA Fitted, Forecasted vs Actual Data for p,d,q: ', num2str(p), ',', num2str(d), ',', num2str(q)]);
    xlabel('Days');
    ylabel('Gold Price');
    
    hold off;
end


% plot_forecast_values_vs_actual(price_next_20_days, days, price);
% plot_fitted_vs_actual(days, price, fitted_data);

function plot_fitted_vs_actual(days, price, fitted_data)
    figure('Position', [300, 60, 1160, 700]); 
    plot(days, price, 'b-', 'LineWidth', 1);
    hold on;
    plot(days, fitted_data, 'r--', 'LineWidth', 2);
    legend('Actual Data', 'Fitted Values');
    title('ARIMA Fitted Values vs Actual Data');
    xlabel('Days');
    ylabel('Gold Price');
    hold off;
end

function [y_forecast_20_days, fitted_data, upper, lower] = arima_fxn(p,d,q,days,price, price_next_20_days)
    arima_model = arima('Constant', NaN, 'MALags', [1:p] ,'ARLags', [1:q], 'D',d);
    estimated_model = estimate(arima_model, price);
    resi_arima = infer(estimated_model, price);
    fitted_data = price + resi_arima;
    
    % Model validation. lbq
    lags = floor(log(length(resi_arima)));
    [h_arima_lbq, p_arima_lbq] = lbqtest(resi_arima, 'lags', lags, 'DOF', lags-2);
    [h_arima_sw, p_arima_sw] = swtest(resi_arima);
    [h_arima_ttest, p_arima_ttest] = ttest(resi_arima);

    %Forecast values
    [y_forecast_20_days, y_mean_sqr_err,] = forecast(estimated_model, 21, 'Y0', price);
    upper = y_forecast_20_days + 1.96 * sqrt(y_mean_sqr_err);
    lower = y_forecast_20_days - 1.96 * sqrt(y_mean_sqr_err);
end

function plot_forecast_values_vs_actual(price_next_20_days,days, price)
    y_forecast_20_days = arima_fxn(0,1,0,days,price, price_next_20_days);
    forecast_price_20_days = y_forecast_20_days(:,1);
    actual_price_20_days = price_next_20_days.Var2;
    days_next_20 = (days(end)+1):(days(end)+21);
    figure('Position', [300, 60, 1160, 700]);
    plot(days_next_20, forecast_price_20_days, 'r--', 'LineWidth', 2);
    hold on;
    plot(days_next_20, actual_price_20_days, 'g-', 'LineWidth', 2);
    legend('Forecasted Price', 'Actual Price');
    title('ARIMA Forecast vs Actual Prices for the Next 20 Days');
    xlabel('Days');
    ylabel('Gold Price');
    hold off;
end


