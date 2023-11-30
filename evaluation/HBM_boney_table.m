%% Create HBM2024 table
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


fontsize     = 15;
boneyEvalDir = fileparts(which('HBM_boney_table.m'));

% load data
men_corr   = cat_io_csv(fullfile(boneyEvalDir,'HBM_2024','men_corr_Boney.csv'));
men_pval   = cat_io_csv(fullfile(boneyEvalDir,'HBM_2024','men_pval_Boney.csv'));
women_corr = cat_io_csv(fullfile(boneyEvalDir,'HBM_2024','women_corr_Boney.csv'));
women_pval = cat_io_csv(fullfile(boneyEvalDir,'HBM_2024','women_pval_Boney.csv'));
%both_corr  = cat_io_csv(fullfile(boneyEvalDir,'HBM_2024','both_corr_Boney.csv'));
%both_pval  = cat_io_csv(fullfile(boneyEvalDir,'HBM_2024','both_pval_Boney.csv'));
both_corr  = men_corr;
both_pval  = men_pval; 
both_corr(2:end,2:end)  = num2cell(mean(cell2mat(cat(3,men_corr(2:end,2:end),women_corr(2:end,2:end))),3));
both_pval(2:end,2:end)  = num2cell(mean(cell2mat(cat(3,men_pval(2:end,2:end),women_pval(2:end,2:end))),3));


% we want to create multipe talbes with specific focus, i.e. to evaluate 
% the bone and then fat measures separatly 
use_vals   = { % x y
  [ 9 10 11 12  ,  2:5 , 1 ]  [ 6 8 7 , 9 10 11 12  ,  2:5 ];   % bone
  [ 9 10 11 12  ,  2:5 , 1 ]  [ 6 8 7 , 9 10 11 12  ,  2:5 ];   % fat
};


for uvi = 1:size(use_vals,1)
  %% prepare smaller tables
  men_corr2   = men_corr(   [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  men_pval2   = men_pval(   [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  women_corr2 = women_corr( [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  women_pval2 = women_pval( [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  both_corr2  = both_corr(  [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  both_pval2  = both_pval(  [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 

  % matrix
  fg = 11; fgh = figure(fg); clf(fg); fgh.Color = [1 1 1]; 
  fgh.Position(3:4) = [800 400]; 
  ha = annotation('textbox',[0.005 0.92 0.99 0.08],'string','hi'); %sprintf( ...
    %'(%d) %s - site %d, N=%d (r~color,p~transparency<%0.03f)',matonly,matonlyn{matonly+1},si, numel(age), pv));
  ha.FontSize             = fontsize * 1.2;
  ha.FontWeight           = 'bold';
  ha.HorizontalAlignment  = 'center';
  ha.EdgeColor            = 'none'; 
  
  subplot('Position',[0.24 0.35 0.22 .5]); 
  im = imagesc( cell2mat(both_corr2(2:end,2:end)) ); 
  im.AlphaData = min(1,max(0,.1 * log10(1 ./ cell2mat(both_pval2(2:end,2:end) )))); 
  ax = gca; ax.TickLength = [0 0];
  ax.XTick = 1:numel( both_corr2(:,1) ); 
  ax.XTickLabel = cat_io_strrep( both_corr2(1,2:end) ,'_','\_');
  ax.XTickLabelRotation = 90;
  ax.YTick = 1:numel( both_corr2(:,1) ); 
  ax.YTickLabel = cat_io_strrep( both_corr2(2:end,1) ,'_','\_');
  ax.FontSize = fontsize; 
  clim([-1 1]);
  title('both')
  
  subplot('Position',[0.46 0.35 0.22 .5]); 
  im = imagesc( cell2mat(women_corr2(2:end,2:end)) ); 
  im.AlphaData = min(1,max(0,.1 * log10(1 ./ cell2mat(women_pval2(2:end,2:end))) )); 
  ax = gca; ax.TickLength = [0 0];
  ax.XTick = 1:numel( women_corr2(:,1) ); 
  ax.XTickLabel = cat_io_strrep( both_corr2(1,2:end) ,'_','\_');
  ax.XTickLabelRotation = 90;
  ax.YTick = []; 
  ax.FontSize = fontsize; 
  clim([-1 1]);
  title('women')
  
  subplot('Position',[0.68 0.35 0.22 .5]); 
  im = imagesc( cell2mat(men_corr2(2:end,2:end)) ); 
  im.AlphaData = min(1,max(0,.1 * log10(1./ cell2mat(men_pval2(2:end,2:end))) )); 
  ax = gca; ax.TickLength = [0 0];
  ax.XTick = 1:numel( women_corr2(:,1) ); 
  ax.XTickLabel = cat_io_strrep( both_corr2(1,2:end) ,'_','\_');
  ax.XTickLabelRotation = 90;
  ax.YTick = []; 
  ax.FontSize = fontsize; 
  clim([-1 1]);
  title('men')

  colormap(jet); 
  subplot('Position',[0.92,0.35,0.0001,0.5]); imagesc( (-1:.1:1)' ); axis equal off
  cb = colorbar; 
  
  subplot('Position',[0.96,0.35,0.0001,0.5]); axis equal off
  cb = colorbar;
  
  %subplot('Position',[0.95,0.35,0.0001,0.5]); 
  %subplot('Position',[0.98,0.35,0.0001,0.5]); imagesc( ones(1,21) ); axis equal off
   
  saveas(gcf,fullfile(boneyEvalDir,'HBM_2024',sprintf('HBM_2024_table%02d.png', uvi ))); 
  

end



%% old code 
      
      if gi==1
        ax.YTick = 1:numel(measure); ax.YTickLabel = cat_io_strrep(measure,'_','\_'); ylabel('measures'); 
      else
        ax.YTick = []; ax.YTickLabel = {};
      end
      ax.XTickLabelRotation = 90;
      switch gi
        case 2,    col = [0 0 1]; 
        case 3,    col = [1 0 0];
        otherwise, col = [0 0 0];
      end
      ax.XColor = col; ax.YColor = col; ax.LineWidth = 2;  
      title(group{gi},'Color',col);
      xlabel('parameters');
      hold on

      % grid
      nx = size(RHOm,1);
      ny = size(RHOm,2);
      edge_x = repmat((0:nx)+0.5, ny+1,1);
      edge_y = repmat((0:ny)+0.5, nx+1,1).';
      plot(edge_x , edge_y , 'Color', repmat(.9,3,1) ) % vertical lines
      plot(edge_x', edge_y', 'Color', repmat(.9,3,1) ) % horizontal lines

      % extra lines to group measures/paras
      nxi = [3 6 7  8  13 15 ];
      plot( repmat(nxi+0.5, 2,1) , repmat([0 ny+1]', 1,numel(nxi)) , 'Color', zeros(3,1) ) % vertical lines
      %nxi = [3 4 10 12 ];
      %plot( repmat(nxi+0.5, 2,1) , repmat([0 ny+1]', 1,numel(nxi)) , 'LineWidth', 2, 'Color', zeros(3,1) ) % vertical lines
      switch matonly 
        case {11,12}
          % nyi = [4 7 13 19 25   10 16 22 28]; nyim = [7 13 19 25]; % full
        case {2,3}
          nyi = [3 4 9 12]; nyim = [];
        otherwise
          nyi = []; nyim = [];
      end
      if ~isempty(nyi)
        plot( repmat([0 nx+1]', 1,numel(nyi)), repmat(nyi+0.5, 2,1)  , 'Color', zeros(3,1) ) % vertical lines
      end
      if ~isempty(nyim)
        plot( repmat([0 nx+1]', 1,numel(nyim)), repmat(nyim+0.5, 2,1)  , 'Color', col, 'LineWidth', 1.5) % vertical lines
      end
    

    colormap(jet); caxis([-1 1]); %cat_io_colormaps('BWR',64)); %graybluered
    subplot('Position',[0.95,0.27,0.0001,0.6]); imagesc( (-1:.1:1)' ); axis equal off
    colorbar

    
    