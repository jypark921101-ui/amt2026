clear; 
close all;
load("230208_CO2_Raw.mat");
Len = length(Box_Char);

Pressure = [1014.9
1014.8
1014.8
1015
1015.5
1016.5]; 

Pressure_Time = 13:18;

filter = [12, 19];

order_previous = 'UNYHJSGCDXPKIQOBRAMT';

order_current = 'KCWFSHUJDMGRPZQXOVNA';

Smooth_Window = 78;
SW = Smooth_Window*2+1;

Smooth_Window_Measure = 60;
SW_M = Smooth_Window_Measure*2+1;

Search_Window = 60;

Time_Ref = datetime([
    2023 2 8 13 43 00;
    2023 2 8 14 42 01;
    2023 2 8 15 20 00;
    2023 2 8 16 20 01;
    2023 2 8 17 20 00;
    2023 2 8 18 19 01]);
Time_Ref = reshape(Time_Ref, [2 3])';

Cal_1.Box_Raw = zeros(seconds(Time_Ref(1, 2)-Time_Ref(1, 1)), Len);
Cal_2.Box_Raw = zeros(seconds(Time_Ref(2, 2)-Time_Ref(2, 1)), Len);
Measure = zeros(seconds(Time_Ref(3, 2)-Time_Ref(3, 1)), Len);
Measure_Diff = zeros(seconds(Time_Ref(3, 2)-Time_Ref(3, 1))-Smooth_Window_Measure*2, Len);

Cal_1.LICOR = zeros(seconds(Time_Ref(1, 2)-Time_Ref(1, 1)), 1);
Cal_2.LICOR = zeros(seconds(Time_Ref(2, 2)-Time_Ref(2, 1)), 1);

Time_Range = {
    timerange(Time_Ref(1, 1), Time_Ref(1, 2));
    timerange(Time_Ref(2, 1), Time_Ref(2, 2))};

for i = 1:Len
    Cal_1.Box_Raw(:,i) = Box{i}(Time_Range{1}, :).CO2;
    Cal_1.T_Raw(:,i) = Box{i}(Time_Range{1}, :).Temp;
    Cal_1.RH_Raw(:,i) = Box{i}(Time_Range{1}, :).RH;

    Cal_2.Box_Raw(:,i) = Box{i}(Time_Range{2}, :).CO2;
    Cal_2.T_Raw(:,i) = Box{i}(Time_Range{2}, :).Temp;
    Cal_2.RH_Raw(:,i) = Box{i}(Time_Range{2}, :).RH;

end

LICOR = retime(LICOR, unique(LICOR.Time), 'mean');

LICOR_1 = retime(LICOR, Time_Ref(1, 1)-seconds(Smooth_Window)-seconds(Search_Window):seconds(1):Time_Ref(1, 2)+seconds(Smooth_Window-1)+seconds(Search_Window), "linear");
LICOR_2 = retime(LICOR, Time_Ref(2, 1)-seconds(Smooth_Window)-seconds(Search_Window):seconds(1):Time_Ref(2, 2)+seconds(Smooth_Window-1)+seconds(Search_Window), "linear");

Cal_1.LICOR = movmean(LICOR_1.CO2, SW, 'Endpoints', 'discard');
Cal_2.LICOR = movmean(LICOR_2.CO2, SW, 'Endpoints', 'discard'); % .LICOR = smoothed vals

Cal_1.LICOR_Raw = LICOR_1.CO2;
Cal_2.LICOR_Raw = LICOR_2.CO2;

TT = NaT(2, Len);

for i = 1:Len
    TT(:, i) = [min(Box{i}.Time); max(Box{i}.Time)];
end

TT = [[min(LICOR.Time); max(LICOR.Time)], TT];

NewTimes = min(TT(1, :)):seconds(1):max(TT(2, :));

AAA = renamevars(LICOR, "CO2", "LICOR");

AAA = retime(AAA, NewTimes);

for i = 1:Len
    Working_Box = Box{i}(:,1);
    Working_Box = renamevars(Working_Box, "CO2", Box_Char(i));
    Working_Box = retime(Working_Box, NewTimes);
    AAA = [AAA, Working_Box];
end

% figure; 
% plot(Cal_1.Box);
% hold on
% plot(Cal_1.LICOR, "LineWidth", 1, "Color", [0,0,0]);
% 
% figure;
% plot(Cal_2.Box);
% hold on
% plot(Cal_2.LICOR, "LineWidth", 1, "Color", [0,0,0]);

%% De-vapor Box CO2 signals

Cal_1.AbH = zeros(size(Cal_1.RH_Raw));

for i = 1:Len
    Box_Conc = Cal_1.Box_Raw(:,i);
    Box_RH = Cal_1.T_Raw(:,i);
    Box_T = Cal_1.RH_Raw(:,i);
    P = (Pressure(1)+Pressure(2))/2;
    P_ws = zeros(size(Box_T));
    for j = 1:length(P_ws)
        T = Box_T(j);
        if T >= 0
            P_ws(j) = 6.1078 * exp(17.27 * T / (T + 237.3));
        else
            P_ws(j) = 6.1078 * exp(21.875 * T / (T + 265.5)); % Tetens' equations
        end
    end
    V_dry = 1 - (P_ws .* Box_RH/100)/P;
    Cal_1.CO2_dry(:, i) = (Box_Conc./V_dry) * (1013/P);

    Cal_1.AbH(:, i) = (P_ws .* Box_RH/100); %also in hpa
end

Cal_2.AbH = zeros(size(Cal_2.RH_Raw));

for i = 1:Len
    Box_Conc = Cal_2.Box_Raw(:,i);
    Box_RH = Cal_2.T_Raw(:,i);
    Box_T = Cal_2.RH_Raw(:,i);
    P = (Pressure(3)+Pressure(4))/2;
    P_ws = zeros(size(Box_T));
    for j = 1:length(P_ws)
        T = Box_T(j);
        if T >= 0
            P_ws(j) = 6.1078 * exp(17.27 * T / (T + 237.3));
        else
            P_ws(j) = 6.1078 * exp(21.875 * T / (T + 265.5)); % Tetens' equations, in hpa
        end
    end
    V_dry = 1 - (P_ws .* Box_RH/100)/P; % No unit
    Cal_2.CO2_dry(:, i) = (Box_Conc./V_dry) * (1013/P); % ppm

    Cal_2.AbH(:, i) = (P_ws .* Box_RH/100); %also in hpa
end

for i = 1:Len
    Box_Conc = Box{i}.CO2;
    Box_RH = Box{i}.RH;
    Box_T = Box{i}.Temp;
    P = (Pressure(5)+Pressure(6))/2;
    P_ws = zeros(size(Box_T));
    for j = 1:length(P_ws)
        T = Box_T(j);
        if T >= 0
            P_ws(j) = 6.1078 * exp(17.27 * T / (T + 237.3));
        else
            P_ws(j) = 6.1078 * exp(21.875 * T / (T + 265.5)); % Tetens' equations
        end
    end
    V_dry = 1 - (P_ws .* Box_RH/100)/P;
    Box_dry{i} = (Box_Conc./V_dry) * (1013/P);
end

%% Calculate timelag from Cal2

stdev_loop = zeros(1,2*Search_Window+1);
stdev = zeros(1, 20);
timelag = zeros(1, 20);
TR_LICOR = timerange(min(LICOR_2.Time), max(LICOR_2.Time), "closed");

for i = 1:Len
    for j = -Search_Window:Search_Window
        LL = Cal_2.LICOR(Search_Window+1-j:end-Search_Window-j);
        BB = Cal_2.CO2_dry(:,i);
        stdev_loop(j+1+Search_Window) = std(LL-BB);
    end
    timelag(i) = min(find(stdev_loop == min(stdev_loop))-Search_Window);
    stdev(i) = min(stdev_loop);
end

%% Linear regression with timelag
Regression_Consts_1 = zeros(Len, 4);
Regression_Consts_2 = zeros(Len, 4);
Regression_Consts_3 = zeros(Len, 5);
Box_Corr = cell(Len);


for i = 1:Len
    AA = Cal_1.LICOR(Search_Window+1-timelag(i):end-Search_Window-timelag(i));
    BB = [Cal_1.T_Raw(:, i), Cal_1.AbH(:, i), Cal_1.CO2_dry(:,i), ones(length(Cal_1.T_Raw(:, i)), 1)];
    C = BB\AA;
    Regression_Consts_1(i, :) = [C(1), C(2), C(3), C(4)]; % T RH (P) CO2 const
    Cal_1.CO2_cor(:, i) = BB * C;

    AA = Cal_2.LICOR(Search_Window+1-timelag(i):end-Search_Window-timelag(i));
    BB = [Cal_2.T_Raw(:, i), Cal_2.AbH(:,i), Cal_2.CO2_dry(:,i), ones(length(Cal_2.T_Raw(:, i)), 1)];
    C = BB\AA;
    Regression_Consts_2(i, :) = [C(1), C(2), C(3), C(4)];
    Cal_2.CO2_cor(:, i) = BB * C;

    % Box_Corr{i} = Box_dry{i} * C;

    AA = Cal_2.LICOR(Search_Window+1-timelag(i):end-Search_Window-timelag(i));
    BB = [Cal_2.T_Raw(:, i), Cal_2.AbH(:,i), Cal_2.CO2_dry(:,i), log(Cal_2.CO2_dry(:,i)), ones(length(Cal_2.T_Raw(:, i)), 1)];
    C = BB\AA;
    Regression_Consts_3(i, :) = [C(1), C(2), C(3), C(4), C(5)];
    Cal_2.CO2_cor_2(:, i) = BB * C;
end

Consts_avg = mean(Regression_Consts_2, 1);
for i = 1:Len
    BB = [Cal_2.T_Raw(:, i), Cal_2.AbH(:,i), Cal_2.CO2_dry(:,i), ones(length(Cal_2.T_Raw(:, i)), 1)];
    CO2_avg_consts(:, i) = BB*Consts_avg';
end

%% Calculate offset from timelag and Cal_2

LICOR_avg = mean(Cal_2.LICOR);
Box_Offset = zeros(1, Len);
Cal_2_Offset.Box = zeros(size(Cal_2.CO2_dry));
for i = 1:Len
    Box_avg = mean(Cal_2.CO2_dry(:, i));
    Box_Offset(i) = Box_avg - LICOR_avg;
end
Cal_2_Offset.LICOR = Cal_2.LICOR;
Cal_2_Offset.Box = Cal_2.CO2_dry - Box_Offset;

Rsq = zeros(20, 1);

for i = 1:Len
    Y_1 = Cal_2_Offset.Box(:,i);
    Y_orig = Cal_2.LICOR(Search_Window+1-timelag(i):end-Search_Window-timelag(i));
    Rsq(i) = 1 - sum((Y_orig-Y_1).^2)/sum((Y_orig-mean(Y_orig)).^2);
end

%% Extract box values with timelag
% Positive timelag means box lags; negative timelag means LICOR lags

index = 1:Len;
for i = 1:length(filter)
    index(index == filter(i)) = [];
end

for i = index
    Time_Range_Measure = timerange(Time_Ref(3, 1) + seconds(timelag(i)), Time_Ref(3, 2) + seconds(timelag(i)));
    Measure(:,i) = Box{i}(Time_Range_Measure, :).CO2 - Box_Offset(i);
    Measure_Diff(:,i) = movmean(Box{i}(Time_Range_Measure, :).CO2, SW_M, "Endpoints", "discard");
end

Measure_Diff = Measure_Diff(2:end, :) - Measure_Diff(1:end-1, :);

save("230208_CO2_All.mat", "Pressure", "Box", "Cal_1", "Cal_2", "timelag", "Measure_Diff", "Measure", "filter", "Box_Char", "order_current", "order_previous", "Time_Ref", "Smooth_Window_Measure", "Search_Window");
%writetimetable(AAA, "230208_SNUStation.xlsx");
save("230208_CO2_Smoothing_Window.mat", "Cal_2", "LICOR", "Time_Ref");

% %%
% figure;
% hold on
% for i = 1:20
%     plot((1:3601)+60+timelag(i), Cal_2_Offset.Box(:, i))
% end
% plot(Cal_2.LICOR, "Color", [0, 0, 0], "LineWidth", 1.5);
% 
% figure;
% hold on
% for i = 1:20
%     plot(Cal_2.LICOR(Search_Window+1-timelag(i):end-Search_Window-timelag(i)), Cal_2.CO2_cor(:, i), "LineStyle","none", "Marker",".", "MarkerEdgeColor",[0,1,1], "MarkerFaceColor",[1,0,0])
% end

%% Error per each step

RMSE_Raw = zeros(Len, 1);
RMSE_Devapor = RMSE_Raw;
RMSE_Timelag = RMSE_Raw;
RMSE_Final = RMSE_Raw;
RMSE_Cor = RMSE_Raw;

for i = 1:Len
    RMSE_Raw(i) = sqrt(mean((Cal_2.LICOR(61:end-60) - Cal_2.Box_Raw(:, i)).^2));
    RMSE_Devapor(i) = sqrt(mean((Cal_2.LICOR(61:end-60) - Cal_2.CO2_dry(:, i)).^2));
    RMSE_Timelag(i) = sqrt(mean((Cal_2.LICOR((1:3601)+60-timelag(i)) - Cal_2.CO2_dry(:, i)).^2));
    RMSE_Final(i) = sqrt(mean(((Cal_2.LICOR((1:3601)+60-timelag(i)) - Cal_2.CO2_dry(:, i)) - (mean(Cal_2.LICOR((1:3601)+60-timelag(i)) - Cal_2.CO2_dry(:, i)))).^2));
    RMSE_Cor(i) = sqrt(mean((Cal_2.LICOR((1:3601)+60-timelag(i)) - Cal_2.CO2_cor(:, i)).^2));
end

RMSE = [RMSE_Raw, RMSE_Devapor, RMSE_Timelag, RMSE_Final, RMSE_Cor];

RMSE_avg_consts = zeros(Len, 1);
for i = 1:Len
    RMSE_avg_consts(i) = sqrt(mean((Cal_2.LICOR((61:end-60)-timelag(i)) - CO2_avg_consts(:, i)).^2));
end
mean(RMSE_avg_consts)

%%
figure;
hold on
for i = 1:Len
    plot(-Cal_2.LICOR((1:3601)+60-timelag(i)) + Cal_2.CO2_cor_2(:, i));
end
%plot(Cal_2.LICOR(1:3601)-mean(Cal_2.LICOR(1:3601)));

%% Longer averages on sensors
Avg_span = 300;
Length_2 = floor(length(Cal_2.CO2_dry)/Avg_span);

LICOR_avg = zeros(Length_2, 1);
for i = 1:Length_2
    LICOR_avg(i) = mean(Cal_2.LICOR(60+(1:Avg_span)+(i-1)*Avg_span));
end

Box_avg = zeros(Length_2, Len);
T_avg = Box_avg;
AbH_avg = Box_avg;
for i = 1:Length_2
    Box_avg(i, :) = mean(Cal_2.CO2_dry((1:Avg_span)+(i-1)*Avg_span, :), 1);
    T_avg(i, :) = mean(Cal_2.T_Raw((1:Avg_span)+(i-1)*Avg_span, :), 1);
    AbH_avg(i, :) = mean(Cal_2.AbH((1:Avg_span)+(i-1)*Avg_span, :), 1);
end

for i = 1:Len
    AA = LICOR_avg;
    BB = [T_avg(:, i), AbH_avg(:,i), Box_avg(:, i), ones(Length_2, 1)];
    C = BB\AA;
    Regression_Consts_4(i, :) = [C(1), C(2), C(3), C(4)];
    Box_avg_cor(:, i) = BB * C;
end

RMSE_avg_mean = zeros(Len, 1);
for i = 1:Len
    RMSE_avg_mean(i) = sqrt(mean((LICOR_avg - Box_avg_cor(:, i)).^2));
end
mean(RMSE_avg_mean)

%% Figure 1
figure;
hold on
for i = 1:Len
    plot(Cal_1.Box_Raw(:, i));
end
plot(Cal_1.LICOR_Raw, "color", "black", "linewidth", 1.5)

figure; 
hold on
for i = 1:Len
    plot(Cal_2.Box_Raw(:, i));
end
plot(Cal_2.LICOR_Raw, "color", "black", "linewidth", 1.5)
