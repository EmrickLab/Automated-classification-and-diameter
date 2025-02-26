% Prompt the user to choose the source/location of the data
[file, path] = uigetfile('*.csv', 'Select the CSV file(s)', 'MultiSelect', 'on');
if isequal(file, 0)
    error('No file selected. Please select at least one CSV file to proceed.');
end

% Check if multiple files are selected
if iscell(file)
    full_paths = fullfile(path, file);
else
    full_paths = {fullfile(path, file)};
end

% Initialize arrays to store results
max_values = [];
ids = [];
info_values = {};
diameters = []; % Initialize array for diameters

% Loop through each selected file
for k = 1:length(full_paths)
    % Read data from the CSV file
    data_table = readtable(full_paths{k});

    % Find columns labeled 'Mean' followed by a number
    mean_columns = contains(data_table.Properties.VariableNames, 'Mean');
    mean_columns = mean_columns & ~cellfun(@isempty, regexp(data_table.Properties.VariableNames, '^Mean\d+$', 'once'));

    % Find columns labeled 'Feret' and 'MinFeret' followed by a number
    feret_columns = contains(data_table.Properties.VariableNames, 'Feret');
    feret_columns = feret_columns & ~cellfun(@isempty, regexp(data_table.Properties.VariableNames, '^Feret\d+$', 'once'));

    min_feret_columns = contains(data_table.Properties.VariableNames, 'MinFeret');
    min_feret_columns = min_feret_columns & ~cellfun(@isempty, regexp(data_table.Properties.VariableNames, '^MinFeret\d+$', 'once'));

    % Ensure that we have the same number of 'Feret' and 'MinFeret' columns
    if sum(feret_columns) ~= sum(min_feret_columns)
        error('The number of "Feret" and "MinFeret" columns do not match.');
    end

    % Loop through each mean column to find the maximum value and corresponding ID
    mean_idx = find(mean_columns);
    feret_idx = find(feret_columns);
    min_feret_idx = find(min_feret_columns);

    for i = 1:length(mean_idx)
        mean_col = data_table{:, mean_idx(i)}; % Extract mean column
        [max_value, row_idx] = max(mean_col); % Find max value and row index
        max_values(end+1) = max_value;
        ids(end+1) = data_table{row_idx, 1}; % Retrieve corresponding ID
        info_values{end+1} = data_table{row_idx, end}; % Retrieve corresponding info

        % Calculate the diameter using the formula 2 * sqrt((feret^2 + min_feret^2) / 2)
        feret = data_table{row_idx, feret_idx(i)};
        min_feret = data_table{row_idx, min_feret_idx(i)};
        diameters(end+1) = 2 * sqrt((feret^2 + min_feret^2) / 2);
    end
end

% Define custom names for each ID
id_names = containers.Map('KeyType', 'double', 'ValueType', 'char');
id_names(1) = 'background';
id_names(8) = 'trpm8 C1/C2';
id_names(2) = 'cd34 C3';
id_names(9) = 's100b C4';
id_names(10) = 'trkB C5';
id_names(11) = 's100b/calca C6';
id_names(6) = 'trpv1 C7/C9/C10';
id_names(7) = 'trpv1/trpa1 C8';
id_names(5) = 'nppb C11';
id_names(4) = 'mrgpra3 C12';
id_names(3) = 'mrgprd C13';
% Add more ID-name mappings as needed

% Extract data for 'background' ID only (id_names(1))
background_id = 1;
background_indices = ids == background_id;
background_diameters = diameters(background_indices);
background_max_values = max_values(background_indices);
background_ids = ids(background_indices);
background_info_values = info_values(background_indices);

% Exclude group 1 (background) from the data
exclude_group = 1;
valid_indices = ids ~= exclude_group;
filtered_ids = ids(valid_indices);
filtered_diameters = diameters(valid_indices);

% Create a table with the filtered diameter values
diameter_table = table(filtered_ids', filtered_diameters', ...
    'VariableNames', {'ID', 'Diameter'});

% Print the diameter table to the console
disp(diameter_table);

% Prompt the user to choose the output directory for saving the files
output_dir = uigetdir('', 'Select Output Directory');
if isequal(output_dir, 0)
    error('No directory selected. Please select a directory to save the files.');
end

% Define the output paths for the CSV file and the figures
output_csv_path = fullfile(output_dir, 'diameter_data.csv');
output_bar_figure_path = fullfile(output_dir, 'diameter_plot.png');
output_pie_figure_path = fullfile(output_dir, 'pie_chart.png');
output_hist_figure_path = fullfile(output_dir, 'histogram.png');

% Write the table to a new CSV file
writetable(diameter_table, output_csv_path);
fprintf('Diameter data has been written to %s\n', output_csv_path);

% Create a bar graph with individual points
figure('Position', [100, 100, 800, 600]); % Increase figure size
unique_ids = unique(filtered_ids);
average_diameters = zeros(size(unique_ids));

for i = 1:length(unique_ids)
    average_diameters(i) = mean(filtered_diameters(filtered_ids == unique_ids(i)));
end

% Bar graph with individual points
bar(unique_ids, average_diameters, 'FaceColor', 'flat');
hold on;
for i = 1:length(unique_ids)
    scatter(repmat(unique_ids(i), sum(filtered_ids == unique_ids(i)), 1), filtered_diameters(filtered_ids == unique_ids(i)), 'k', 'filled');
end
hold off;

xlabel('ID');
ylabel('Diameter');
title('Average Diameter with Individual Points');
xticks(unique_ids);
xticklabels(arrayfun(@(x) id_names(x), unique_ids, 'UniformOutput', false));

% Save the bar graph figure
saveas(gcf, output_bar_figure_path);
fprintf('Diameter plot has been saved to %s\n', output_bar_figure_path);

% Continue with the rest of the code as before

% Remove values labeled as 'background' (id_names(1))
valid_indices = ids ~= background_id;
max_values = max_values(valid_indices);
ids = ids(valid_indices);
info_values = info_values(valid_indices);
diameters = diameters(valid_indices);

% Aggregate values by ID
unique_ids = unique(ids);
aggregated_values = zeros(size(unique_ids));

for i = 1:length(unique_ids)
    aggregated_values(i) = sum(max_values(ids == unique_ids(i)));
end

% Generate pie chart without default labels
figure('Position', [100, 100, 800, 600]); % Increase figure size
p = pie(aggregated_values); % Create pie chart
percentValues = get(p(2:2:end), 'String'); % Get percentage values
labels = arrayfun(@(x) id_names(x), unique_ids, 'UniformOutput', false); % Get custom names for IDs
combinedLabels = strcat(labels', ': ', percentValues); % Combine labels with percentage values

% Clear the default labels
for k = 1:length(p)/2
    p(2*k).String = '';
end

title(sprintf('Distribution of Max Values (n=%d)', length(max_values)));

% Add legend to avoid label overlap
legend(combinedLabels, 'Location', 'eastoutside');

% Save the pie chart figure
saveas(gcf, output_pie_figure_path);
fprintf('Pie chart has been saved to %s\n', output_pie_figure_path);

% Ensure unique_ids is a row vector for consistent concatenation
unique_ids = unique_ids(:)'; % Ensure unique_ids is a row vector
id_counts = histcounts(ids, [unique_ids, unique_ids(end) + 1]);
figure('Position', [100, 100, 800, 600]); % Increase figure size
bar(unique_ids, id_counts);
xlabel('ID Groups');
ylabel('Count');
title('Count of Each ID Group');
xticks(unique_ids);
xticklabels(arrayfun(@(x) id_names(x), unique_ids, 'UniformOutput', false)); % Set custom names for x-axis labels

% Save the histogram figure
saveas(gcf, output_hist_figure_path);
fprintf('Histogram has been saved to %s\n', output_hist_figure_path);