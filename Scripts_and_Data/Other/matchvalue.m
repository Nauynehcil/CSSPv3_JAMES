function [sim_mon_new, sta_x_new, sta_y_new] = matchvalue(sim_mon, sta_mon, sta_x, sta_y, tl, res)
% MATCHVALUE: Optimizes station-to-grid matching using KGE within a dynamic search window.
sta_num = length(sta_x);
avgvalue = mean(sim_mon, 3, 'omitnan');
% Pre-allocation
sim_mon_new = zeros(size(sta_mon));
sta_x_new = zeros(sta_num,1);
sta_y_new = zeros(sta_num,1);
% Calculate search radius R based on resolution (e.g., 0.5° -> R=1; 0.25° -> R=2)
R = round(0.5 / res); 
for grid = 1:sta_num
    % Generate coordinates for the neighborhood search window
    [sta_xc, sta_yc] = xc_yc_dynamic(sta_x(grid), sta_y(grid), R);  
    num_neighbors = length(sta_xc);
    index = ones(num_neighbors, 1) * -999; 
    % Data integrity check: Process stations with sufficient valid months
    if sum(isnan(sta_mon(grid, :)), 'all') <= (tl - 12)
        avg_obs = mean(sta_mon(grid, :), 2, 'omitnan');
        for i = 1:num_neighbors
            % Boundary protection: Ensure indices are within the grid extent
            if sta_xc(i) > 0 && sta_yc(i) > 0 && ...
               sta_xc(i) <= size(sim_mon, 1) && sta_yc(i) <= size(sim_mon, 2)
                avg_sim = avgvalue(sta_xc(i), sta_yc(i));
                ratio = avg_sim / avg_obs;
                % Quality Control: Filter out grids with obsvious misjudgment
                if avg_sim <= 0 || ratio > 5 || ratio < 1/5
                    index(i) = -999;
                else
                    sim = squeeze(sim_mon(sta_xc(i), sta_yc(i), :));
                    % Compute KGE as the similarity metric
                    index(i) = KGE(sim, sta_mon(grid, :).'); 
                end
            end
        end
    end
    % Selection: Identify the grid with the highest KGE score
    [max_val, pos_m] = max(index);
    if max_val == -999
        sim_mon_new(grid, :) = nan;
        sta_x_new(grid) = sta_x(grid);
        sta_y_new(grid) = sta_y(grid);
    else
        % Map the best grid data and store its updated coordinates
        sim_mon_new(grid, :) = sim_mon(sta_xc(pos_m), sta_yc(pos_m), :);
        sta_x_new(grid) = sta_xc(pos_m);
        sta_y_new(grid) = sta_yc(pos_m);
    end
end
end

function [xc, yc] = xc_yc_dynamic(x, y, R)
% Generate coordinate offsets based on dynamic radius R
[dx, dy] = meshgrid(-R:R, -R:R);
xc = x + dx(:);
yc = y + dy(:);
end