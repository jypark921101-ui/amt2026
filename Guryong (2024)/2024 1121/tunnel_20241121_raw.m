clear;

%% Preparation
Time_Ref = datetime([
    2024 11 20 20 00 00; % inject 990; 398 injected in 22:47
    2024 11 20 23 00 00;
    2024 11 21 11 00 00; % Completed installation
    2024 11 22 7 00 00]); % Started removal 

Filenames = ["O", "A", "W", ...
    "D", "V", "J", "I", "Q", "G", ...
    "Z", "B", "C", "P", "M", "U", ...
    "R", "S", "X", "F", "K", "H"];

Box_locations = [0, 100, 200, ...
    300, 350, 400, 450, 500, 550, ...
    600, 650, 700, 750, 800, 850, ...
    900, 950, 1000, 1050, 1100, 1190]; % 250m omitted

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
Cal_ind_1 = [find(datetime([2024 11 20 21 20 00]) == Newtimes); find(datetime([2024 11 20 21 22 00]) == Newtimes)];
Cal_ind_2 = [find(datetime([2024 11 20 22 59 00]) == Newtimes); find(datetime([2024 11 20 23 00 00]) == Newtimes)];

CO2_Cal_1 = Box.CO2_raw(Cal_ind_1(1):Cal_ind_1(2), :);
CO2_Cal_1 = mean(CO2_Cal_1, 1);

CO2_Cal_2 = Box.CO2_raw(Cal_ind_2(1):Cal_ind_2(2), :);
CO2_Cal_2 = mean(CO2_Cal_2, 1);

Twopt_A = (Cal_conc_1-Cal_conc_2)./(CO2_Cal_1 - CO2_Cal_2);
Twopt_B = (Cal_conc_2*CO2_Cal_1 - Cal_conc_1*CO2_Cal_2)./(CO2_Cal_1 - CO2_Cal_2);

CO2_linear_cal = Box.CO2_raw(Time_indices(3):Time_indices(4), :) .* Twopt_A + Twopt_B;

% %% Modify Traffic data
% load("Traffic_data.mat")
% Traffic_newtime = Time_Ref(3):minutes(1):Time_Ref(4); 
% Traffic_interp = interp1(Traffic_time, Traffic ,Traffic_newtime);

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
tiledlayout("vertical")

ax1 = nexttile;
hold on;
imagesc(Newtimes2, [0 1], CO2_linear_cal')
set(ax1, "Ytick", [])
annotation('textbox',[.05 .70 .1 .2], ...
    'String',"    Exit",'EdgeColor','none')
annotation('textbox',[.05 .10 .1 .2], ...
    'String',"Entrance",'EdgeColor','none')
a = colorbar;
a.Label.String = 'CO_2 (ppm)';

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