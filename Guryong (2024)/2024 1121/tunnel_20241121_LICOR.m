clear;

%% Preparation
Time_Ref = datetime([
    2024 11 20 20 00 00; % completed 1st point
    2024 11 20 23 00 00; % completed 2nd point
    2024 11 21 11 00 00; % Completed installation
    2024 11 22 7 00 00]); % Started removal 

Filenames = ["O", "A", "W", ...
    "D", "V", "J", "I", "Q", "G", ...
    "Z", "B", "C", "P", "M", "U", ...
    "R", "S", "X", "F", "K", "H"];

Box_locations = [0, 100, 200, ...
    300, 350, 400, 450, 500, 550, ...
    600, 650, 700, 750, 800, 850, ...
    900, 950, 1000, 1050, 1100, 1180]; % 250m omitted

Box_locations_new = 0:10:1190;

Cal_conc_1 = 398;
Cal_conc_2 = 990;

FilePath = "Raw data/";

Box_Num = length(Filenames);

Newtimes = Time_Ref(1):seconds(1):Time_Ref(4);

Box.raw = cell(Box_Num, 1);
Box.retimed = Box.raw;
Box.CO2_raw = zeros(length(Newtimes), Box_Num);
Box.Temp = Box.CO2_raw;
Box.RH = Box.CO2_raw;

%% Read xls files
for i = 1:Box_Num
    Name = FilePath + Filenames(i) + ".xls";
    temp = readtable(Name);
    time_temp = datetime(temp.Time, "InputFormat", 'MM-dd-yy HH:mm:ss');
    Box.raw{i} = timetable(time_temp, temp.CO2_PPM_, temp.Temp___, temp.Humi__RH_, 'VariableNames', ["CO2", "Temp", "RH"]);
end

%% Read LICOR file
temp = readtable("241121.txt");
time_temp = datetime(temp.Date_Y_M_D_ + temp.Time_H_M_S_, "InputFormat", "uuuu-MM-dd HH:mm:ss");
LICOR.raw =  timetable(time_temp, temp.CO2_ppm_, temp.Cell_Temperature_C_, temp.H2O_ppt_, 'VariableNames', ["CO2", "Temp", "H2O"]);

%% Retime everything to 1 second
for i = 1:Box_Num
    Box.retimed{i} = retime(Box.raw{i}, Newtimes, "linear");
    Box.CO2_raw(:, i) = Box.retimed{i}.CO2;
    Box.Temp(:, i) = Box.retimed{i}.Temp;
    Box.RH(:, i) = Box.retimed{i}.RH;
end

LICOR_times = [min(LICOR.raw.time_temp)+seconds(300), max(LICOR.raw.time_temp)];

Newtimes_LICOR = LICOR_times(1):seconds(1):LICOR_times(2);
LICOR.retimed = retime(LICOR.raw, Newtimes_LICOR, "linear");
LICOR.smooth = timetable((LICOR_times(1)+seconds(60):seconds(1):LICOR_times(2)-seconds(60))', movmean(LICOR.retimed.CO2, 121, "endpoints", "discard"), 'variablenames', "CO2"); 
LICOR.CO2_raw = LICOR.retimed.CO2;
LICOR.Temp = LICOR.retimed.Temp;
LICOR.H2O = LICOR.retimed.H2O;
Time_ind_LICOR = find(Newtimes == LICOR_times(1)+seconds(60)):find(Newtimes == LICOR_times(2)-seconds(60));
%% Generate index numbers for each sessions
Time_indices = zeros(4, 1);
for i = 1:4
    Time_indices(i) = find(Time_Ref(i) == Newtimes);
end

%% Calibrate (2pt)
Cal_ind_1 = [find(datetime([2024 11 20 21 20 00]) == Newtimes); find(datetime([2024 11 20 21 22 00]) == Newtimes)];
Cal_ind_2 = [find(datetime([2024 11 20 22 59 00]) == Newtimes); find(datetime([2024 11 20 23 00 00]) == Newtimes)];

CO2_Cal_1 = Box.CO2_raw(Cal_ind_1(1):Cal_ind_1(2), :);
CO2_Cal_1 = mean(CO2_Cal_1, 1);

CO2_Cal_2 = Box.CO2_raw(Cal_ind_2(1):Cal_ind_2(2), :);
CO2_Cal_2 = mean(CO2_Cal_2, 1);

Twopt_A = (Cal_conc_1-Cal_conc_2)./(CO2_Cal_1 - CO2_Cal_2);
Twopt_B = (Cal_conc_2*CO2_Cal_1 - Cal_conc_1*CO2_Cal_2)./(CO2_Cal_1 - CO2_Cal_2);

CO2_linear_cal = (Box.CO2_raw(Time_indices(3):Time_indices(4), :) .* Twopt_A) + Twopt_B;

%% Calibrate (LICOR, regression)
Regression_Consts_1 = zeros(2, Box_Num);
for i = 1:Box_Num
    AA = LICOR.smooth.CO2;
    BB = [Box.retimed{i}.CO2(Time_ind_LICOR), ones(length(Time_ind_LICOR), 1)];
    C = BB\AA;
    Regression_Consts_1(:, i) = [C(1); C(2)];
    Box.CO2_corr(:, i) = [Box.retimed{i}.CO2, ones(length(Box.retimed{i}.CO2), 1)] * C;
end

CO2_corr = Box.CO2_corr(Time_indices(3):Time_indices(4), :);

%% Plot (2pt)
Color = flipud(crameri("glasgow"));

Newtimes2 = Time_Ref(3):seconds(1):Time_Ref(4);

% figure;
% tiledlayout("vertical")
% ax0 = nexttile;
% hold on;
% scatter(Box_locations, zeros(size(Box_locations)), 40, Color, "filled");
% axis 'auto x'
% ax0.YAxis.Visible="off";
% set(ax0,'Color','None');
% xlabel("distance from tunnel entrance (meters)");
% 
% ax1 = nexttile;
% hold on;
% for i = 1:Box_Num
%     plot(Newtimes2, CO2_linear_cal(:, i), "Color", Color(i, :));
% end
% ylabel("corrected CO_2 concentration (ppm)");
% 
% ax2 = nexttile;
% hold on;
% plot(Traffic_newtime, Traffic_interp(:, 1)/5);
% legend(Traffic_key(1));
% ylabel("Vehicles per minute");
% 
% ax3 = nexttile;
% hold on;
% plot(Traffic_newtime, Traffic_interp(:, 2:6)/5);
% legend(Traffic_key(2:6));
% ylabel("Vehicles per minute");

% figure;
% tiledlayout("vertical")
% 
% ax1 = nexttile;
% hold on;
% imagesc(Newtimes2, [0 1], CO2_linear_cal')
% set(ax1, "Ytick", [])
% annotation('textbox',[.05 .70 .1 .2], ...
%     'String',"    Exit",'EdgeColor','none')
% annotation('textbox',[.05 .10 .1 .2], ...
%     'String',"Entrance",'EdgeColor','none')
% a = colorbar;
% a.Label.String = 'CO_2 (ppm)';

% ax2 = nexttile;
% hold on;
% plot(Traffic_newtime, Traffic_interp(:, 1)/5);
% legend(Traffic_key(1));
% ylabel("Vehicles per minute");
% 
% ax3 = nexttile;
% hold on;
% plot(Traffic_newtime, Traffic_interp(:, 2:6)/5);
% legend(Traffic_key(2:6));
% ylabel("Vehicles per minute");

% %% 
% [min(CO2_linear_cal);median(CO2_linear_cal);mean(CO2_linear_cal);max(CO2_linear_cal)]

% %% Plot (regression)
% figure;
% tiledlayout("vertical")
% 
% ax1 = nexttile;
% hold on;
% imagesc(Newtimes2, [0 1], CO2_corr')
% set(ax1, "Ytick", [])
% annotation('textbox',[.05 .70 .1 .2], ...
%     'String',"    Exit",'EdgeColor','none')
% annotation('textbox',[.05 .10 .1 .2], ...
%     'String',"Entrance",'EdgeColor','none')
% a = colorbar;
% a.Label.String = 'CO_2 (ppm)';


%% Calibration using timelag and regression
Regression_Consts_2 = zeros(2, Box_Num);
timelag = zeros(Box_Num, 1);

for i = 1:Box_Num
    
    diff_temp = zeros(241, 1);

    for t = -120:120
        diff_temp(t+121) = rmse(LICOR.smooth.CO2 - mean(LICOR.smooth.CO2), Box.retimed{i}.CO2(Time_ind_LICOR + t) - mean(Box.retimed{i}.CO2(Time_ind_LICOR + t)));
    end

    timelag(i) = find(diff_temp == min(diff_temp)) - 120;

    AA = LICOR.smooth.CO2;
    BB = [Box.retimed{i}.CO2(Time_ind_LICOR + t), ones(length(Time_ind_LICOR), 1)];
    C = BB\AA;
    Regression_Consts_2(:, i) = [C(1); C(2)];
    Box.CO2_corr_2(:, i) = [Box.retimed{i}.CO2, ones(length(Box.retimed{i}.CO2), 1)] * C;
end

for i = 1:Box_Num
    RMSE.cal2(i) = rmse(Box.CO2_corr_2(Time_ind_LICOR + t, i), LICOR.smooth.CO2);
end

for i = 1:Box_Num
    Box.CO2_corr_2_LICOR(:, i) = Box.CO2_corr_2(Time_ind_LICOR + t, i);
    CO2_corr_2(:, i) = Box.CO2_corr_2((Time_indices(3):Time_indices(4))+timelag(i) - 120, i);
end



%% Comparison of two calibration methods

% Extract sensor values belonging to LICOR

Box_comp.raw = zeros(length(Time_ind_LICOR), Box_Num);
Box_comp.cal = Box.CO2_corr(Time_ind_LICOR, :);

for i = 1:Box_Num
    Box_comp.raw(:, i) = Box.retimed{i}.CO2(Time_ind_LICOR);
end

Box_comp.tpt = (Box.CO2_raw(Time_ind_LICOR, :) .* Twopt_A) + Twopt_B;

RMSE.raw = rmse(Box_comp.raw, LICOR.smooth.CO2);
RMSE.tpt = rmse(Box_comp.tpt, LICOR.smooth.CO2);
RMSE.cal = rmse(Box_comp.cal, LICOR.smooth.CO2);

figure; 
tile1 = tiledlayout("vertical")
ax(1) = nexttile;

boxplot([RMSE.raw', RMSE.tpt', RMSE.cal', RMSE.cal2']);
xticklabels([])
ax(1).FontSize = 12;
title('2024-Nov-21', "FontSize", 16)
% text(1, median(RMSE.raw), num2str(median(RMSE.raw)))
% text(2, median(RMSE.tpt), num2str(median(RMSE.tpt)))

ax(2) = nexttile;
boxplot([RMSE.raw', RMSE.tpt', RMSE.cal', RMSE.cal2']);
xticklabels({'Raw data', 'Two point', 'Multi-point', '...w/ timelag compensation'})
ylim([0 50])
ax(2).FontSize = 12;
% text(3, median(RMSE.cal), num2str(median(RMSE.cal)))
% text(4, median(RMSE.cal2), num2str(median(RMSE.cal2)))

RMSE_median = median([RMSE.raw; RMSE.tpt; RMSE.cal; RMSE.cal2], 2);
ylabel(tile1, "RMSE(ppm)")


%% Preparing for plot
CO2_linear_cal_int = zeros([72001, 120]);
CO2_corr_int = CO2_linear_cal_int;
CO2_corr_2_int = CO2_linear_cal_int;
for i = 1:72001
    CO2_linear_cal_int(i, :) = interp1(Box_locations, CO2_linear_cal(i, :), Box_locations_new, "nearest");
    CO2_corr_int(i, :) = interp1(Box_locations, CO2_corr(i, :), Box_locations_new, "nearest");
    CO2_corr_2_int(i, :) = interp1(Box_locations, CO2_corr_2(i, :), Box_locations_new, "nearest");
end



%% Plot (both)

figure;
tile1 = tiledlayout(4, 1);
ax(1) = nexttile;
hold on;
imagesc(Newtimes2, [0 1190], Box.CO2_raw(Time_indices(3):Time_indices(4), :)')
xtickformat('HH:mm')
xsecondarylabel(Visible="off")
title("a) Raw data", "FontSize", 16)
ylim([-50 1240])
ax(1).FontSize = 12;

ax(2) = nexttile;
hold on;
imagesc(Newtimes2, [0 1190], CO2_linear_cal_int')
xtickformat('HH:mm')
xsecondarylabel(Visible="off")
title("b) 2-point calibration", "FontSize", 16)
ylim([-50 1240])
ax(2).FontSize = 12;

ax(3) = nexttile;
hold on;
imagesc(Newtimes2, [0 1190], CO2_corr_int')
xtickformat('HH:mm')
xsecondarylabel(Visible="off")
title("c) Multivalue linear regression", "FontSize", 16)
ylim([-50 1240])
ax(3).FontSize = 12;

ax(4) = nexttile;
hold on;
imagesc(Newtimes2, [0 1190], CO2_corr_2_int')
xtickformat('HH:mm')
title("d) Multivalue linear regression w/ timelag compensation", "FontSize", 16)
ylim([-50 1240])
ax(4).FontSize = 12;

cbh = colorbar(ax(end)); 
set(ax, 'Colormap', Color, 'CLim', [400 1100])
cbh.Layout.Tile = 'east'; 
cbh.Label.String = 'CO_2 (ppm)';
cbh.FontSize = 12;

xlabel(tile1, "Time (hh:mm)", "FontSize", 16)
ylabel(tile1, "Distance from the entrance(m)", "FontSize", 16)

%% Traffic

Traffic_time = datetime([2024 11 21 7 00 00]):minutes(5):datetime([2024 11 22 7 00 00]);
Traffic_speed = [14.48
31.93
27.77
17.93
14.25
14.67
17.74
15.55
14.31
47.64
46.7
58.9
55.35
55.24
58.44
61.43
63.6
65.85
64.07
46.35];

Traffic_speed_time = datetime([2024 11 21 11 30 00]):hours(1):datetime([2024 11 22 7 00 00]);

load("Traffic_Nov.mat")
Traffic_newtime = Time_Ref(3):minutes(1):Time_Ref(4); 
Traffic_interp = interp1(Traffic_time, Traffic ,Traffic_newtime);
Traffic_speed_interp = interp1(Traffic_speed_time', Traffic_speed, Traffic_newtime', "nearest", "extrap");

Car_count = sum(Traffic_interp, 2) ./ (Traffic_speed_interp * 1000 / 60) * 1190;


figure;
Tile2 = tiledlayout("vertical");

ax1 = nexttile;
hold on;
imagesc(Newtimes2, [0 1190], CO2_corr_2_int')
ylabel(["Distance from", "the entrance (m)"])
a = colorbar;
a.Label.String = 'CO_2 (ppm)';
ax1.FontSize = 12;
ylim([-50 1240])
xsecondarylabel(Visible="off")
xtickformat('HH:mm')
title("a) Corrected CO_2", "FontSize", 16)
set(ax1, 'Colormap', Color)

ax2 = nexttile;
hold on;
yyaxis right
bar(Traffic_speed_time, Traffic_speed, "FaceColor", [1, 0.8, 0.7], "EdgeAlpha", 0);
ylabel("Average traffic speed (km h^-^1)");
yyaxis left
plot(Traffic_newtime, Car_count, "LineWidth", 1.5);
ylabel("Total vehicles count in the tunnel");
title("b) Hourly average traffic speed and estimated vehicle count", "FontSize", 16)

set(gca, 'SortMethod', 'depth')
ax2.FontSize = 12;
xsecondarylabel(Visible="off")
xtickformat('HH:mm')

ax3 = nexttile;
hold on;
plot(Traffic_newtime, Traffic_interp);
legend(Traffic_key,"Location","eastoutside");
ylabel(["Number of vehicles" "entering per minute" "(Number min^-^1)"]);
ax3.FontSize = 12;
xsecondarylabel(Visible="off")
xtickformat('HH:mm')
title("c) Entering vehicles", "FontSize", 16)

xlabel(Tile2, "Time (hh:mm)")

%%

% Newtime_3_1 = datetime([2024 11 21 11 00 00]);
% Newtime_3_2 = datetime([2024 11 21 12 00 00]);
% 
% Newtimes3 = Newtime_3_1:seconds(1):Newtime_3_2;
% Traffic_newtime3 = Newtime_3_1:minutes(1):Newtime_3_2; 
% 
% Newindex_3_1 = find(Newtimes2 == Newtime_3_1);
% Newindex_3_2 = find(Newtimes2 == Newtime_3_2);
% Newindex_3 = Newindex_3_1:Newindex_3_2;
% 
% Newindex_T3_1 = find(Traffic_newtime == Newtime_3_1);
% Newindex_T3_2 = find(Traffic_newtime == Newtime_3_2);
% Newindex_T3 = Newindex_T3_1:Newindex_T3_2;
% 
% figure;
% tile1 = tiledlayout("vertical")
% 
% ax1 = nexttile;
% hold on;
% imagesc(Newtimes3, [0 1190], CO2_corr_2_int(Newindex_3, :)')
% ylabel(["Distance from", "the entrance (m)"])
% a = colorbar;
% a.Label.String = 'CO_2 (ppm)';
% ax1.FontSize = 12;
% xsecondarylabel(Visible="off")
% ylim([-50 1240])
% xlim([Newtime_3_1 Newtime_3_2])
% 
% ax2 = nexttile;
% hold on;
% plot(Traffic_newtime3, Traffic_interp(Newindex_T3, 1));
% legend(Traffic_key(1),"Location","eastoutside");
% xsecondarylabel(Visible="off")
% ylabel("Vehicles per minute");
% ax2.FontSize = 12;
% 
% ax3 = nexttile;
% hold on;
% plot(Traffic_newtime3, Traffic_interp(Newindex_T3, 2:6));
% legend(Traffic_key(2:6), "Location","eastoutside");
% ylabel("Vehicles per minute");
% ax3.FontSize = 12;
% 
% xlabel(tile1, "Time(hh:mm)", "FontSize", 16)


%%
load("July_measurements.mat")

Newtime_3_1 = datetime([2024 11 21 11 00 00]);
Newtime_3_2 = datetime([2024 11 21 12 00 00]);

Newtimes3 = Newtime_3_1:seconds(1):Newtime_3_2;
Traffic_newtime3 = Newtime_3_1:minutes(1):Newtime_3_2; 

Newindex_3_1 = find(Newtimes2 == Newtime_3_1);
Newindex_3_2 = find(Newtimes2 == Newtime_3_2);
Newindex_3 = Newindex_3_1:Newindex_3_2;

Newtimes_Jul = Time_Ref_July(3):seconds(1):Time_Ref_July(4);

figure;
tile1 = tiledlayout("vertical");

ax11(1) = nexttile;
hold on;
imagesc(Newtimes_Jul, [0 1180], CO2_linear_cal_July')
ax11(1).FontSize = 12;
ylim([-50 1240])
xlim([Time_Ref_July(3) Time_Ref_July(4)])
title(ax11(1), "a)")
ax11(1).TitleHorizontalAlignment = "left";
yticks([0 500 1000 1180]);

ax11(2) = nexttile;
hold on;
imagesc(Newtimes3, [0 1180], CO2_corr_2_int(Newindex_3, :)')
ax11(2).FontSize = 12;
ylim([-50 1240])
xlim([Newtime_3_1 Newtime_3_2])
title(ax11(2), "b)")
ax11(2).TitleHorizontalAlignment = "left";
yticks([0 500 1000 1180]);

xlabel(tile1, "Time(hh:mm)", "FontSize", 16)
ylabel(tile1, "Distance from the entrance (m)", "FontSize", 16)

cbh = colorbar(ax11(end)); 
set(ax11, 'Colormap', Color, 'CLim', [400 1100])
cbh.Layout.Tile = 'east'; 
cbh.Label.String = 'CO_2 (ppm)';
cbh.FontSize = 12;