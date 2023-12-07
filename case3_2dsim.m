dim_x = 10;
dim_y = 10;
c_tf = 0.3;
c_cell = 0.3;
states = 5;

tissue = ones(dim_x, dim_y); 
index_x = randi([1, dim_x]);
index_y = randi([1, dim_y]);
tissue(index_x, index_y) = 2; 
tissue_cmap = [247,247,247; 0,0,0; 244,165,130; 202,0,32; 146,197,222] / 255;        
colormap(tissue_cmap);

videoFile = 'sim_video_final.avi';
videoObj = VideoWriter(videoFile);
open(videoObj);

fig = figure(1);
ax = axes('Position', [0, 0, 1, 1]);
image(tissue);
axis equal
axis tight
axis ij
gcaOpts = {'XTick', [], 'YTick', []};
set(gca, gcaOpts{:}); 
set(gcf, 'Color', 'w');

% initialize grid
tf_grid = cell(dim_x, dim_y);
for i = 1:dim_x
    for j = 1:dim_y
        tf = zeros(1, 6);
        tf_grid{i, j} = tf;
    end
end

tf_grid{index_x, index_y}(1) = 2500;
tf_grid{index_x, index_y}(6) = 1000000;

T = 10;
imageCounter = 0;

for t = 1:T
  fprintf('%d/%d (Ctrl-C to stop)\r', t, T);
  for row = 1:dim_x
    for col = 1:dim_y
        curr = tf_grid{row, col};

        %% solve ODE on [0,0.001]
        y0 = curr;
        t_range = [0 1];
        [time, y] = ode23(@case3_ode, t_range, y0);
        %% disp(y(11,:));

        %% discrete diffusion
        nei = zeros(4, 6);

        % Top
        if (row - 1) >= 1 && (row - 1) <= dim_x && col >= 1 && col <= dim_y
            for i = 5:6
                tf_grid{row - 1, col}(i) = tf_grid{row - 1, col}(i) + c_tf*(tf_grid{row,col}(i) - tf_grid{row - 1, col}(i));
            end
            for i = 1:4
                tf_grid{row - 1, col}(i) = tf_grid{row - 1, col}(i) + c_cell*(tf_grid{row,col}(i) - tf_grid{row - 1, col}(i));
            end
            nei(1,:) = tf_grid{row - 1, col};
        else
            nei(1,:) = zeros(1,6);
        end

        % Bottom
        if (row + 1) >= 1 && (row + 1) <= dim_x && col >= 1 && col <= dim_y
            for i = 5:6
                tf_grid{row + 1, col}(i) = tf_grid{row + 1, col}(i) + c_tf*(tf_grid{row,col}(i) - tf_grid{row + 1, col}(i));
            end
            for i = 1:4
                tf_grid{row + 1, col}(i) = tf_grid{row + 1, col}(i) + c_cell*(tf_grid{row,col}(i) - tf_grid{row + 1, col}(i));
            end
            nei(2,:) = tf_grid{row + 1, col};
        else
            nei(2,:) = zeros(1,6);
        end

        % Left
        if row >= 1 && row <= dim_x && (col - 1) >= 1 && (col - 1) <= dim_y
            for i = 5:6
                tf_grid{row, col - 1}(i) = tf_grid{row, col - 1}(i) + c_tf*(tf_grid{row,col}(i) - tf_grid{row, col - 1}(i));
            end
            for i = 1:4
                tf_grid{row, col - 1}(i) = tf_grid{row, col - 1}(i) + c_cell*(tf_grid{row,col}(i) - tf_grid{row, col - 1}(i));
            end
            nei(3,:) = tf_grid{row, col - 1};
        else
            nei(3,:) = zeros(1,6);
        end

        % Right
        if row >= 1 && row <= dim_x && (col + 1) >= 1 && (col + 1) <= dim_y
            for i = 5:6
                tf_grid{row, col + 1}(i) = tf_grid{row, col + 1}(i) + c_tf*(tf_grid{row,col}(i) - tf_grid{row, col + 1}(i));
            end
            for i = 1:4
                tf_grid{row, col + 1}(i) = tf_grid{row, col + 1}(i) + c_cell*(tf_grid{row,col}(i) - tf_grid{row, col + 1}(i));
            end
            nei(4,:) = tf_grid{row, col + 1};
        else
            nei(4,:) = zeros(1,6);
        end
        
        for i = 1:4
            for j = 5:6
                if nei(i,:) ~= 0
                    tf_grid{row, col}(j) = c_tf * (nei(i,j) - tf_grid{row, col}(j));
                end 
            end
            for j = 1:4
                if nei(i,:) ~= 0
                    tf_grid{row, col}(j) = c_cell * (nei(i,j) - tf_grid{row, col}(j));
                end 
            end
        end

        %% update current state
        tf_grid{row, col} = y(11,:);
        disp(tf_grid{row, col});

        %% convert to one state (fig)
        fg = [y(11,1), y(11,2), y(11,3), y(11,4)];
        disp("fg");
        disp(fg);
        [M,I] = max(fg);
        % disp("index");
        % disp(I + 1);
        tissue(row, col) = I + 1;
        disp(tissue(row, col));
    end
  end

  if mod(t, 1) == 0
      set(gcf, 'CurrentAxes', ax);
      tissue_cmap = [247,247,247; 0,0,0; 244,165,130; 202,0,32; 146,197,222] / 255;                    
      colormap(tissue_cmap);
      image(tissue);
      set(gca, gcaOpts{:});
      set(gcf, 'PaperSize', [dim_y, dim_x]/100, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, dim_y, dim_x]/100);
      frame = getframe(gcf);
      writeVideo(videoObj, frame);
      drawnow;
      pause(0.01);
   end  
end

close(videoObj);
fprintf('Saved as "%s".\n', videoFile);