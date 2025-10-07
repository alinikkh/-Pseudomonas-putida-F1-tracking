%% Extract data from plots
% Get handle to current figure
fig = gcf;

% Find all surface objects (your trajectories)
surfaces = findall(fig, 'Type', 'Surface');

% Preallocate storage
n_surfaces = numel(surfaces);
all_X = cell(n_surfaces, 1);
all_Y = cell(n_surfaces, 1);
all_C = cell(n_surfaces, 1);  % This is the straightness

for i = 1:n_surfaces
    % Get X, Y and C (straightness) data
    all_X{i} = surfaces(i).XData;
    all_Y{i} = surfaces(i).YData;
    all_C{i} = surfaces(i).CData;
end

disp('Extraction complete!');

clean_X = cellfun(@(x) x(1,:), all_X, 'UniformOutput', false);
clean_Y = cellfun(@(y) y(1,:), all_Y, 'UniformOutput', false);
clean_C = cellfun(@(c) c(1,:), all_C, 'UniformOutput', false);

%% ===== VIDEO 1: Straightness (colormap) =====
v = VideoWriter('trajectory_animation_fast_fixed.mp4','MPEG-4');
v.Quality = 100;
v.FrameRate = 25;  % high frame rate
open(v);

%% ===== FIGURE SETUP =====
fig = figure('Color','w','Position',[100 100 1920 1080]); % HD figure
hold on; axis equal; grid on;

xlim([-40 40]); ylim([-40 40]);
xticks([-40 0 40]); yticks([-40 0 40]);
set(gca,'BoxStyle','full','LineWidth',2.5, ...
    'DataAspectRatio',[1 1 1],'FontSize',24);
box on;

xlabel('x-x_0 (µm)','FontSize',24,'Interpreter','tex');
ylabel('y-y_0 (µm)','FontSize',24,'Interpreter','tex');

%% ===== COLOR MAP =====
all_straightness = cellfun(@mean, clean_C);  % optional, for info
cmap = cool(64); colormap(cmap);
caxis([0.4,0.5]);
cb = colorbar;
cb.Label.String = 'Straightness';
cb.Label.FontSize = 24;
cb.Ticks = linspace(0.4,0.5,2);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f',v), cb.Ticks,'UniformOutput',false);

%% ===== SETUP HANDLES =====
max_t = max(cellfun(@length, clean_X));
n_tracks = numel(clean_X);
s_handles = gobjects(n_tracks,1);

for i = 1:n_tracks
    s_handles(i) = plot(NaN, NaN, 'LineWidth', 2.5); % use plot for consistent thickness
end

%% ===== ANIMATION LOOP =====
for t = 1:max_t
    for i = 1:n_tracks
        x = clean_X{i}; y = clean_Y{i}; c = clean_C{i};
        if length(x) >= t
            xt = x(1:t); yt = y(1:t); ct = c(1:t);

            % Clip straightness for colormap (values <0.4 -> 0.4, >0.5 -> 0.5)
            clipped_ct = min(max(ct, 0.4), 0.5);

            % Map to colormap index
            col_idx = round((mean(clipped_ct)-0.4)/(0.5-0.4)*63) + 1;

            % Update line
            set(s_handles(i), 'XData', xt, 'YData', yt, 'Color', cmap(col_idx,:));
        end
    end
    
    % Capture frame directly from figure
    frame = getframe(fig); 
    writeVideo(v, frame);
end

close(v);
%%

%%

% I add markers so the tracks width doesnt change over time

%% ===== EXTRACT DATA FROM FIGURE =====
fig = gcf;

% Find all surface objects (your trajectories)
surfaces = findall(fig, 'Type', 'Surface');

% Preallocate storage
n_surfaces = numel(surfaces);
all_X = cell(n_surfaces, 1);
all_Y = cell(n_surfaces, 1);
all_C = cell(n_surfaces, 1);  % Straightness

for i = 1:n_surfaces
    all_X{i} = surfaces(i).XData;
    all_Y{i} = surfaces(i).YData;
    all_C{i} = surfaces(i).CData;
end

disp('Extraction complete!');

% Flatten to single row
clean_X = cellfun(@(x) reshape(x(1,:), 1, []), all_X, 'UniformOutput', false);
clean_Y = cellfun(@(y) reshape(y(1,:), 1, []), all_Y, 'UniformOutput', false);
clean_C = cellfun(@(c) reshape(c(1,:), 1, []), all_C, 'UniformOutput', false);

%% ===== SETUP VIDEO =====
v = VideoWriter('trajectory_animation_markers.mp4','MPEG-4');
v.Quality = 100;
v.FrameRate = 40; % reasonable frame rate
open(v);

%% ===== FIGURE SETUP =====
fig = figure('Color','w','Position',[100 100 1920 1080]);
hold on; axis equal; grid on;
xlim([-40 40]); ylim([-40 40]);
xticks([-40 0 40]); yticks([-40 0 40]);
set(gca,'BoxStyle','full','LineWidth',2.0,'DataAspectRatio',[1 1 1],'FontSize',24);
box on;

xlabel('x-x_0 (µm)','FontSize',24,'Interpreter','tex');
ylabel('y-y_0 (µm)','FontSize',24,'Interpreter','tex');

%% ===== COLOR MAP =====
cmap = cool(64); 
colormap(cmap);
caxis([0.0, 1]);
cb = colorbar;
cb.Label.String = 'Straightness';
cb.Label.FontSize = 24;
cb.Ticks = linspace(0,1,6);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f',v), cb.Ticks,'UniformOutput',false);

%% ===== SETUP PLOT HANDLES =====
max_t = max(cellfun(@length, clean_X));
n_tracks = numel(clean_X);
s_handles = gobjects(n_tracks,1);

% Use markers for all tracks
for i = 1:n_tracks
    s_handles(i) = plot(NaN, NaN, 'LineWidth', 2.0, 'Marker', '.', 'MarkerSize', 6);
end

%% ===== FORCE OPENGL RENDERER =====
set(gcf,'Renderer','opengl');

%% ===== ANIMATION LOOP =====
for t = 1:max_t
    for i = 1:n_tracks
        x = clean_X{i}; y = clean_Y{i}; c = clean_C{i};
        if length(x) >= t
            xt = x(1:t); 
            yt = y(1:t); 
            ct = c(1:t);

            % Clip straightness
            clipped_ct = min(max(ct, 0.0), 1);

            % Map to colormap index
            %col_idx = round((mean(clipped_ct)-0.4)/(0.5-0.4)*63) + 1;
            col_idx = round((mean(clipped_ct)-0)/(1-0)*63) + 1;
            % Update line and markers
            set(s_handles(i), 'XData', xt, 'YData', yt, 'Color', cmap(col_idx,:));
        end
    end

    drawnow; % Force redraw
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);
disp('Animation complete!');
