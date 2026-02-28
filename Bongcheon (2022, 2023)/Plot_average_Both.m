clear;

Date = ["220608"; "230208"];

load(Date(1) + '_CO2_timelag')

FN = "average.png";

Coords_Center = [37.481111, 126.952778];

Coords = [
37.480409	126.953235
37.481306	126.953625
37.481448	126.952671
37.481374	126.952945
37.482505	126.952853
37.480785	126.951668
37.483089	126.955296
37.480287	126.951115
37.48024	126.951836
37.479673	126.953716
37.482792	126.951861
37.481867	126.953668
37.481042	126.952881
37.482574	126.954414
37.480744	126.95333
37.479407	126.951006
37.479471	126.954155
37.483504	126.951667
37.481117	126.952602
37.481965   126.951885
];

Coords_Order = [
    "A", "B", "C", "D", "G", "H", "I", "J", "K", "M", "N", "O", "P", "Q",...
    "R", "S", "T", "U", "X", "Y"
];

%% Generate indices so I can match current results to previous results
% Measure(:, Index) generates order identical to Coords_Order

index = zeros(length(Box_Char), 1);
Box_Index = strings(length(Box_Char), 1);

% Step 1: relabel boxes; Find what box now corresponds to what box previously
for i = 1:length(Box_Char)
    i1 = order_current == Box_Char(i); 
    Box_Index(i) = order_previous(i1);
end

% Step 2: convert labelling to box order

% first, eliminate any missing boxes' coordinates
Empty = [];
for i = length(Coords_Order)
    if isempty(find(Box_Index == Coords_Order(i), 1))
        Empty = [Empty, i];
    end
end

Coords(Empty, :) = [];
Coords_Order(Empty) = [];

% then, find if anything in Box_Index is the same as Coords_Order
j = 1;
for i = 1:length(Coords_Order)
    index(i) = find(Box_Index == Coords_Order(i), 1);
end

% Step 3: Switch order to "correct" order
Measure = Measure(:, index);
Measure_Diff = Measure_Diff(:, index);

% Step 4: Apply Filter
filter = find(sum(index == filter, 2));

Measure(:, filter) = [];
Measure_Diff(:, filter) = [];
Coords(filter, :) = [];

%% 2D plot

Factor_Y = earthRadius * 2 * pi / 360 / 60 / 60; % How long is 1"?
Factor_X = Factor_Y * cos(Coords_Center(1) * 2*pi / 360);

Coords_X = (Coords(:,2)-Coords_Center(2))*3600 * Factor_X;
Coords_Y = (Coords(:,1)-Coords_Center(1))*3600 * Factor_Y;

Coords_X_1 = Coords(:,2);
Coords_Y_1 = Coords(:,1);

Range_X = (floor(min(Coords_X))):30:(ceil(max(Coords_X)));
Range_Y = (floor(min(Coords_Y))):30:(ceil(max(Coords_Y)));
[Mesh_X, Mesh_Y] = meshgrid(Range_X, Range_Y);

Axis_X = Range_X/3600/Factor_X + Coords_Center(2);
Axis_Y = Range_Y/3600/Factor_Y + Coords_Center(1);
LimX = [126.9503  126.9550];
LimY = [37.4790   37.4838];

AA = imread("map_en_outlier.png");
RA = imref2d(size(AA), LimX, LimY);
imAlpha = 256-rgb2gray(AA);
Black = zeros(size(AA));

Time = Time_Ref(3, 1):seconds(1):Time_Ref(3, 2);

%for t = 1:length(Measure)

F = scatteredInterpolant(Coords_X, Coords_Y, mean(Measure(1:2000, :), 1)', 'linear', 'none');

Vq = F(Mesh_X, Mesh_Y);
scAlpha=ones(size(Vq));
scAlpha(isnan(Vq))=0;

fig2D = figure;
Tile = tiledlayout(2,2,'Padding','tight');
nexttile;

hold on
title("Avg. Value from " + string(Time(1)) + " - " + string(Time(end)));
s = imagesc(Axis_X, Axis_Y, Vq, "AlphaData", scAlpha);

caxis([450, 520])
colormap(parula)
c = colorbar;
c.Label.String = "CO_2 concentration (ppm)";

h = imshow(Black, RA);
h.AlphaData = imAlpha;
axis equal
axis([RA.XWorldLimits, RA.YWorldLimits]);
set(gca,'YDir','normal')

plot(Coords_X_1, Coords_Y_1, "LineStyle", "none", "Marker", "o", "MarkerEdgeColor", "black", "MarkerFaceColor", "black")

%% Stdev

F.Values = std(Measure(1:2000, :))';

Vq = F(Mesh_X, Mesh_Y);
%scAlpha=abs(Vq);
scAlpha=ones(size(Vq));
scAlpha(isnan(Vq))=0;

nexttile;

hold on
title("Stdev of conc. from " + string(Time(1)) + " - " + string(Time(end)));
s = imagesc(Axis_X, Axis_Y, Vq, "AlphaData", scAlpha);

caxis([0, 25])
colormap(flipud(brewermap([], "RdYlBu")))
c = colorbar;
c.Label.String = "Stdev of CO2 concentration (ppm)";

h = imshow(Black, RA);
h.AlphaData = imAlpha;
axis equal
axis([RA.XWorldLimits, RA.YWorldLimits]);
set(gca,'YDir','normal')

plot(Coords_X_1, Coords_Y_1, "LineStyle", "none", "Marker", "o", "MarkerEdgeColor", "black", "MarkerFaceColor", "black")

%% second loop

load(Date(2) + '_CO2_timelag')

Coords_Center = [37.481111, 126.952778];

Coords = [
37.480409	126.953235
37.481306	126.953625
37.481448	126.952671
37.481374	126.952945
37.482505	126.952853
37.480785	126.951668
37.483089	126.955296
37.480287	126.951115
37.48024	126.951836
37.479673	126.953716
37.482792	126.951861
37.481867	126.953668
37.481042	126.952881
37.482574	126.954414
37.480744	126.95333
37.479407	126.951006
37.479471	126.954155
37.483504	126.951667
37.481117	126.952602
37.481965   126.951885
];

Coords_Order = [
    "A", "B", "C", "D", "G", "H", "I", "J", "K", "M", "N", "O", "P", "Q",...
    "R", "S", "T", "U", "X", "Y"
];

%% Generate indices so I can match current results to previous results
% Measure(:, Index) generates order identical to Coords_Order

index = zeros(length(Box_Char), 1);
Box_Index = strings(length(Box_Char), 1);

% Step 1: relabel boxes; Find what box now corresponds to what box previously
for i = 1:length(Box_Char)
    i1 = order_current == Box_Char(i); 
    Box_Index(i) = order_previous(i1);
end

% Step 2: convert labelling to box order

% first, eliminate any missing boxes' coordinates
Empty = [];
for i = length(Coords_Order)
    if isempty(find(Box_Index == Coords_Order(i), 1))
        Empty = [Empty, i];
    end
end

Coords(Empty, :) = [];
Coords_Order(Empty) = [];

% then, find if anything in Box_Index is the same as Coords_Order
j = 1;
for i = 1:length(Coords_Order)
    index(i) = find(Box_Index == Coords_Order(i), 1);
end

% Step 3: Switch order to "correct" order
Measure = Measure(:, index);
Measure_Diff = Measure_Diff(:, index);

% Step 4: Apply Filter
filter = find(sum(index == filter, 2));

Measure(:, filter) = [];
Measure_Diff(:, filter) = [];
Coords(filter, :) = [];

%% 2D plot

Factor_Y = earthRadius * 2 * pi / 360 / 60 / 60; % How long is 1"?
Factor_X = Factor_Y * cos(Coords_Center(1) * 2*pi / 360);

Coords_X = (Coords(:,2)-Coords_Center(2))*3600 * Factor_X;
Coords_Y = (Coords(:,1)-Coords_Center(1))*3600 * Factor_Y;

Coords_X_1 = Coords(:,2);
Coords_Y_1 = Coords(:,1);

Range_X = (floor(min(Coords_X))):30:(ceil(max(Coords_X)));
Range_Y = (floor(min(Coords_Y))):30:(ceil(max(Coords_Y)));
[Mesh_X, Mesh_Y] = meshgrid(Range_X, Range_Y);

Axis_X = Range_X/3600/Factor_X + Coords_Center(2);
Axis_Y = Range_Y/3600/Factor_Y + Coords_Center(1);
LimX = [126.9503  126.9550];
LimY = [37.4790   37.4838];

AA = imread("map_en_outlier.png");
RA = imref2d(size(AA), LimX, LimY);
imAlpha = 256-rgb2gray(AA);
Black = zeros(size(AA));

Time = Time_Ref(3, 1):seconds(1):Time_Ref(3, 2);

%for t = 1:length(Measure)

F = scatteredInterpolant(Coords_X, Coords_Y, mean(Measure, 1)', 'linear', 'none');

Vq = F(Mesh_X, Mesh_Y);
scAlpha=ones(size(Vq));
scAlpha(isnan(Vq))=0;

nexttile;

hold on
title("Avg. Value from " + string(Time(1)) + " - " + string(Time(end)));
s = imagesc(Axis_X, Axis_Y, Vq, "AlphaData", scAlpha);

caxis([450, 520])
colormap(parula)
c = colorbar;
c.Label.String = "CO_2 concentration (ppm)";

h = imshow(Black, RA);
h.AlphaData = imAlpha;
axis equal
axis([RA.XWorldLimits, RA.YWorldLimits]);
set(gca,'YDir','normal')

plot(Coords_X_1, Coords_Y_1, "LineStyle", "none", "Marker", "o", "MarkerEdgeColor", "black", "MarkerFaceColor", "black")

%% Stdev

F.Values = std(Measure)';

Vq = F(Mesh_X, Mesh_Y);
%scAlpha=abs(Vq);
scAlpha=ones(size(Vq));
scAlpha(isnan(Vq))=0;

nexttile;

hold on
title("Stdev of conc. from " + string(Time(1)) + " - " + string(Time(end)));
s = imagesc(Axis_X, Axis_Y, Vq, "AlphaData", scAlpha);

caxis([0, 25])
colormap(flipud(brewermap([], "RdYlBu")))
c = colorbar;
c.Label.String = "Stdev of CO2 concentration (ppm)";

h = imshow(Black, RA);
h.AlphaData = imAlpha;
axis equal
axis([RA.XWorldLimits, RA.YWorldLimits]);
set(gca,'YDir','normal')

plot(Coords_X_1, Coords_Y_1, "LineStyle", "none", "Marker", "o", "MarkerEdgeColor", "black", "MarkerFaceColor", "black")

