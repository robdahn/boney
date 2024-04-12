%% create figure  11
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________



    %  New measures for both sexes, males and femals
    if numel(age)>=1000, pv = 0.0001; elseif numel(age)>300, pv = 0.001; elseif numel(age)>100, pv = 0.01; else, pv = 0.05; end
    
    fg=11; fgh = figure(fg); clf(fg); fgh.Color = [1 1 1]; fgh.Position(3:4) = [1000 150 + 10*numel(measure) ]; 
    ha = annotation('textbox',[0.005 0.92 0.99 0.08],'string',sprintf( ...
      '(%d) %s%s - site %d, N=%d (r~color,p~transparency<%0.03f)',matonly,matonlyn{matonly+1},strrep(resdir,'_','\_'),si, numel(age), pv));
    ha.FontSize             = 13;
    ha.FontWeight           = 'bold';
    ha.HorizontalAlignment  = 'center';
    ha.EdgeColor            = 'none'; 

    rGMV = out.rGMV{1}; rWMV = out.rWMV{1}; rCMV = out.rCMV{1}; TIV = out.TIV{1};
      if si == 0; msk = site>0; else, msk = site==si; end % set up to process all sites for si == 0
              
    % 3 group by sex
    group = {'both sexes','males','females'};
    for gi = 1:3
      subplot('Position',[0.155  + 0.26*(gi-1)  max(0.20,min(.80, 0.60 - 0.025*numel(measure)))  ...
                          0.255                 max(0.15,min(.70, 0.05 + 0.04*numel(measure))) ]); 
      switch gi
        case 1, sexm = ' '; 
        case 2, sexm = ' & sex(:)==1 '; 
        case 3, sexm = ' & sex(:)==2 ';
      end
      for mi = 1:numel(para)
        for mj = 1:numel(measure)
          p1 = eval( sprintf('%s(        ~isnan(%s(:)) %s & msk(:) )', para{mi}   , para{mi}, sexm ));
          p2 = eval( sprintf('out.%s{1}( ~isnan(%s(:)) %s & msk(:) )', measure{mj}, para{mi}, sexm ));
          if size(p1,1)<size(p1,2), p1 = p1'; end
          if size(p2,1)<size(p2,2), p2 = p2'; end
          try
            [RHOm(mi,mj,gi),PVALm(mi,mj,gi)] = corr( p1 , p2 ,'type', corrtype); 
          catch
            RHOm(mi,mj,gi) = 0; 
            PVALm(mi,mj,gi) = 0; 
          end
        end
      end
      
      im = imagesc((RHOm(:,:,gi)')); im.AlphaData = pv ./ PVALm(:,:,gi)'; 
      ax = gca; ax.TickLength = [0 0];
      ax.XTick = 1:numel(para);    ax.XTickLabel = cat_io_strrep(para   ,'_','\_');
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
      nxi = [3 6 7  8  13 15 22];
      plot( repmat(nxi+0.5, 2,1) , repmat([0 ny+1]', 1,numel(nxi)) , 'Color', zeros(3,1) ) % vertical lines
      %nxi = [3 4 10 12 ];
      %plot( repmat(nxi+0.5, 2,1) , repmat([0 ny+1]', 1,numel(nxi)) , 'LineWidth', 2, 'Color', zeros(3,1) ) % vertical lines
      switch matonly 
        case {11,12}
          %nyi = [4 7 13 19 25   10 16 22 28]; nyim = [7 13 19 25]; % full
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
      colormap(jet);
      caxis([-1 1]); %#ok<CAXIS>
    end

     %cat_io_colormaps('BWR',64)); %graybluered
    subplot('Position',[0.95,0.27,0.0001,0.6]); imagesc( (-1:.1:1)' ); axis equal off
    colorbar;
    
    saveas(gcf,fullfile(resultdir,sprintf('mt%d_%s_site%d_n%d_%d.png',matonly,matonlyn{matonly+1},si,numel(age),plevel ) )); 
    

