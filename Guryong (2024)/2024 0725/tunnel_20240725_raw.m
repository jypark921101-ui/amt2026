clear;

%% Preparation
Time_Ref = datetime([
    2024 7 24 21 57 00; % inject 990; 398 injected in 22:47
    2024 7 25 0 26 00;
    2024 7 25 9 28 00; % Completed installation
    2024 7 25 11 02 00]); % Started removal 

Filenames = ["X", "B", "Q", "D", ...
    "V", "K", "W", "N", "C", "U", ...
    "S", "H", "G", "O", "F", "R", ...
    "P", "A", "Z", "J", "M", "I"];

Box_locations = [0, 100, 200, 250, ...
    300, 350, 400, 450, 500, 550, ...
    600, 650, 700, 750, 800, 850, ...
    900, 950, 1000, 1050, 1100, 1180];

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

%% Retime everything to 1 second
for i = 1:Box_Num
    Box.retimed{i} = retime(Box.raw{i}, Newtimes, "linear");
    Box.CO2_raw(:, i) = Box.retimed{i}.CO2;
    Box.Temp(:, i) = Box.retimed{i}.Temp;
    Box.RH(:, i) = Box.retimed{i}.RH;
end

%% Generate index numbers for each sessions
Time_indices = zeros(4, 1);
for i = 1:4
    Time_indices(i) = find(Time_Ref(i) == Newtimes);
end

%% Calibrate (2pt)
Cal_ind_1 = [find(datetime([2024 7 24 22 45 00]) == Newtimes); find(datetime([2024 7 24 22 46 00]) == Newtimes)];
Cal_ind_2 = [find(datetime([2024 7 25 00 26 00]) == Newtimes); find(datetime([2024 7 25 00 27 00]) == Newtimes)];

CO2_Cal_1 = Box.CO2_raw(Cal_ind_1(1):Cal_ind_1(2), :);
CO2_Cal_1 = mean(CO2_Cal_1, 1);

CO2_Cal_2 = Box.CO2_raw(Cal_ind_2(1):Cal_ind_2(2), :);
CO2_Cal_2 = mean(CO2_Cal_2, 1);

Twopt_A = (990-398)./(CO2_Cal_1 - CO2_Cal_2);
Twopt_B = (398*CO2_Cal_1 - 990*CO2_Cal_2)./(CO2_Cal_1 - CO2_Cal_2);

CO2_linear_cal = Box.CO2_raw(Time_indices(3):Time_indices(4), :) .* Twopt_A + Twopt_B;

%% Modify Traffic data
load("Traffic_data.mat")
Traffic_newtime = Time_Ref(3):minutes(1):Time_Ref(4); 
Traffic_interp = interp1(Traffic_time, Traffic ,Traffic_newtime);

%% Preparing for plotting
Box_locations_new = 0:10:1190;
CO2_linear_cal_int = zeros([5641, 120]);
for i = 1:5641
    CO2_linear_cal_int(i, :) = interp1(Box_locations, CO2_linear_cal(i, :), Box_locations_new, "nearest");
end


%% Plot
Color = parula(Box_Num);

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

figure;
tile1 = tiledlayout("vertical")

ax1 = nexttile;
hold on;
imagesc(Newtimes2, [0 1190], CO2_linear_cal')
ylabel(["Distance from", "the entrance (m)"])
a = colorbar;
a.Label.String = 'CO_2 (ppm)';
ax1.FontSize = 12;
xsecondarylabel(Visible="off")
ylim([-50 1240])

ax2 = nexttile;
hold on;
plot(Traffic_newtime, Traffic_interp(:, 1)/5);
legend(Traffic_key(1),"Location","eastoutside");
xsecondarylabel(Visible="off")
ylabel("Vehicles per minute");
ax2.FontSize = 12;

ax3 = nexttile;
hold on;
plot(Traffic_newtime, Traffic_interp(:, 2:6)/5);
legend(Traffic_key(2:6), "Location","eastoutside");
ylabel("Vehicles per minute");
ax3.FontSize = 12;

xlabel(tile1, "Time(hh:mm)", "FontSize", 16)

% %% 
% [min(CO2_linear_cal);median(CO2_linear_cal);mean(CO2_linear_cal);max(CO2_linear_cal)]
CO2_linear_cal_July = CO2_linear_cal;
Traffic_interp_July = Traffic_interp;
Time_Ref_July = Time_Ref;
save("July_measurements.mat", "CO2_linear_cal_July", "Traffic_interp_July", "Time_Ref_July")
