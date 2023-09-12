%boney_eval_printCrossCorrelation

%% create figure 12
  %  new measures sorted by correlations
  fg=12; fgh = figure(fg); clf(fg); fgh.Position(3:4) = [800 350] + [15 5] * max(0,numel(measure)-10); fgh.Color = [1 1 1];

  for mi = 1:numel(measure)
    for mj = 1:numel(measure)
      msk = ~isnan(out.(measure{mi}){1}) & ~isnan(out.(measure{mj}){1});
      [RHOmm(mi,mj),PVALmm(mi,mj)] = corr(out.(measure{mi}){1}(msk),out.(measure{mj}){1}(msk),'type',corrtype); 
    end
  end
  
  for i=1:2
    if i==1
      subplot(1,2,i); 
      im = imagesc(RHOmm); im.AlphaData = pv ./ PVALmm; 
      ax = gca; ax.TickLength = [0 0];
      ax.XTick = 1:numel(measure); ax.XTickLabel = cat_io_strrep(measure,'_','\_');
      ax.YTick = 1:numel(measure); ax.YTickLabel = cat_io_strrep(measure,'_','\_');
      title(sprintf('correlation measures (p=%0.4f)',pv));
    else  
      subplot(1,2,i); 
      [~,RS] = sortrows(mean(RHOmm(:,1:10),2),'descend');
      %[~,RS] = sortrows(std(RHOmm)' + median(RHOmm)' + mean(RHOmm)','descend');
      im = imagesc(RHOmm(RS,RS)); im.AlphaData = pv ./ PVALmm(RS,RS); 
      ax = gca; ax.TickLength = [0 0];
      ax.XTick = 1:numel(measure); ax.XTickLabel = cat_io_strrep(measure(RS),'_','\_');
      ax.YTick = 1:numel(measure); ax.YTickLabel = cat_io_strrep(measure(RS),'_','\_');
      title(sprintf('sorted correlation measures (p=%0.4f)',pv)); 
    end
    xlabel('measures'); ylabel('measures')
    ax.XTickLabelRotation = 90;
    axis equal; hold on; xlim([0 size(RHOmm,1)] + 0.5); ylim([0 size(RHOmm,2)] + 0.5)
    colormap(jet) %cat_io_colormaps('BWR',64)); %graybluered
   
    % grid
    nx = size(RHOmm,1);
    ny = size(RHOmm,2);
    edge_x = repmat((0:nx)+0.5, ny+1,1);
    edge_y = repmat((0:ny)+0.5, nx+1,1).';
    plot(edge_x , edge_y , 'Color', repmat(.9,3,1) ) % vertical lines
    plot(edge_x', edge_y', 'Color', repmat(.9,3,1) ) % horizontal lines

    % extra lines to group measures/paras
    if i==1
      nxi = [3 4 7 9 12]; nyi = nxi; 
      plot( repmat(nxi+0.5, 2,1) , repmat([0 ny+1]', 1,numel(nxi)) , 'Color', repmat(0,3,1) ) % vertical lines
      plot( repmat([0 nx+1]', 1,numel(nyi)), repmat(nyi+0.5, 2,1)  , 'Color', repmat(0,3,1) ) % vertical lines
    end
  end

  saveas(gcf,fullfile(resultdir,sprintf('mt%d_site%d_correlation_measures.png',matonly,si)))
  