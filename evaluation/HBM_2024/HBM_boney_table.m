%% Create HBM2024 table
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


fontsize     = 13;
boneyEvalDir = fileparts(which('HBM_boney_table.m'));


% load correlation data prepared in R
men_corr   = cat_io_csv(fullfile(boneyEvalDir,'men_corr_Boney.csv'));
men_pval   = cat_io_csv(fullfile(boneyEvalDir,'men_pval_Boney.csv'));
women_corr = cat_io_csv(fullfile(boneyEvalDir,'women_corr_Boney.csv'));
women_pval = cat_io_csv(fullfile(boneyEvalDir,'women_pval_Boney.csv'));
both_corr  = cat_io_csv(fullfile(boneyEvalDir,'boney_correlation_proper.csv')); 
both_pval  = cat_io_csv(fullfile(boneyEvalDir,'boney_cor_pval.csv'));

% map upper diagnal values (corrected for multiple comparisions) to lower side
both_corrv = cell2mat( both_corr(2:end,2:end) ); 
both_corrv = both_corrv .* triu(ones(size(both_corrv)),1)  +  both_corrv' .* tril(ones(size(both_corrv)));
both_corr(2:end,2:end) = num2cell( both_corrv ); 
both_corr(2,:) = []; both_corr(:,2) = [];
both_pval(2,:) = []; both_pval(:,2) = [];
both_pvalv = cell2mat( both_pval(2:end,2:end) ); 
both_pvalv = both_pvalv .* triu(ones(size(both_pvalv)),1)  +  both_pvalv' .* tril(ones(size(both_pvalv)));
both_pval(2:end,2:end) = num2cell( both_pvalv ); 

% replace some strings in the header lines
both_corr{contains(both_corr(:,1)  ,'femur_BMD'),1} = 'BMD_femur';
both_corr{1,contains(both_corr(1,:),'femur_BMD')  } = 'BMD_femur';

% We want to create multipe talbes with specific focus, i.e. to evaluate 
% the bone and then fat measures separatly or together. 
% Define also some vertical/horicental lines to group the measures. 
%   use_values = [ x-axes entries, y-axis entries, name, x-lines, y-lines
use_vals   = { 
  [ 9 10 12  ,  2:5 , 1 13 ,      ]  [  7 8  , 9 10 12  ,  2:5 ] 'bone' [2 5 6] [3 7 8 9];   % bone [BMD, brain, age]
  [ 2:5 , 1 15:21                  ]  [ 14      ,            2:5 ] 'fat'   [1 2  ] [1 4 5   ] ;    % fat  [
  [ 9 10 12  ,  2:5 , 1 13 , 15 19 16:18 20:21 ]  [ 6 7 8 14, 9 10 12 ,  2:5 ] 'bone & fat' [3 7 8 9 14]  [3 4 7];   % bone [BMD, brain, age]
};


% loop over the cases defined per use_vals rows
for uvi = 1:size(use_vals,1)
  %% prepare smaller tables
  men_corr2   = men_corr(   [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  men_pval2   = men_pval(   [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  women_corr2 = women_corr( [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  women_pval2 = women_pval( [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  both_corr2  = both_corr(  [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 
  both_pval2  = both_pval(  [1 1+use_vals{uvi,1}] , [1 1+use_vals{uvi,2}])'; 

  % create figure
  fg = 11; fgh = figure(fg);  clf(fg);  fgh.Color = [1 1 1];  colormap(jet);
  fgh.Position(3:4) = [100 + 60*size( men_corr2 ,2 )   100 + 20*size( men_corr2 ,1 ) ]; 

  % title 
  ha = annotation('textbox',[0.005 0.92 0.99 0.08],'string',...
    sprintf(['Correlation plot of %s measures with color coded r-values ' ...
    'and transparency for p-values.'], use_vals{uvi,3}));
  ha.FontSize             = fontsize * 1.1;
  ha.FontWeight           = 'bold';
  ha.HorizontalAlignment  = 'center';
  ha.EdgeColor            = 'none'; 
  

  % plot for both sexes - print a matix a a specific position 
  sb1 = subplot('Position', [0.20  0.45 - 0.01*size( men_corr2 ,1 )  0.23   .3 + 0.02*size( men_corr2 ,1 ) ]); hold on
  im                      = imagesc( cell2mat(both_corr2(2:end,2:end)) ); 
  im.AlphaData            = min(1,max(0,1/6 * log10(10e-2 ./ cell2mat(both_pval2(2:end,2:end) )))); 
  for uvix = 1:size(men_corr2,1),  plot([uvix uvix]+.5,[1 22]-.5,'Color',[0.9 0.9 0.9]); end
  for uvix = 1:size(men_corr2,2), plot([1 22]-.5,[uvix uvix]+.5,'Color',[0.9 0.9 0.9]); end
  for uvix = use_vals{uvi,5}, plot([uvix uvix]+.5,[1 size(both_corr2,1)]-.5,'Color',[0 0 0]); end
  for uvix = use_vals{uvi,4}, plot([1 size(both_corr2,2)]-.5,[uvix uvix]+.5,'Color',[0 0 0]); end
  xlim([0.5 size(both_corr2,2)-.5]); ylim([0.5 size(both_corr2,1)-.5]); caxis([-1 1]);
  % no we have to update the axis labels 
  ax                      = gca; 
  ax.TickLength           = [0 0];
  ax.XTick                = 1:numel( both_corr2(1,2:end) ); 
  ax.XTickLabel           = cat_io_strrep( both_corr2(1,2:end) ,'_','\_');
  ax.XTickLabelRotation   = 90;
  ax.YTick                = 1:numel( both_corr2(2:end,1) ); 
  ax.YTickLabel           = cat_io_strrep( both_corr2(2:end,1) ,'_','\_');
  ax.FontSize             = fontsize; 
  ax.YDir                 = 'reverse';
  title('both sexes','FontSize',fontsize); box on
  

  % women
  sb2 = subplot('Position', sb1.Position + [0.232 0 0 0]); hold on
  im                      = imagesc( cell2mat(women_corr2(2:end,2:end)) ); 
  im.AlphaData            = min(1,max(0,1/6 * log10(10e-2 ./ cell2mat(women_pval2(2:end,2:end))) )); 
  for uvix = 1:size(men_corr2,1),  plot([uvix uvix]+.5,[1 22]-.5,'Color',[0.9 0.9 0.9]); end
  for uvix = 1:size(men_corr2,2), plot([1 22]-.5,[uvix uvix]+.5,'Color',[0.9 0.9 0.9]); end
  for uvix = use_vals{uvi,5}, plot([uvix uvix]+.5,[1 size(both_corr2,1)]-.5,'Color',[0 0 0]); end
  for uvix = use_vals{uvi,4}, plot([1 size(both_corr2,2)]-.5,[uvix uvix]+.5,'Color',[0 0 0]); end
  xlim([0.5 size(both_corr2,2)-.5]); ylim([0.5 size(both_corr2,1)-.5]); caxis([-1 1]);
  ax = gca; ax.TickLength = [0 0];
  ax.XTick                = 1:numel( both_corr2(1,2:end) ); 
  ax.XTickLabel           = cat_io_strrep( both_corr2(1,2:end) ,'_','\_');
  ax.XTickLabelRotation   = 90;
  ax.YTick                = []; 
  ax.FontSize             = fontsize; 
  ax.YDir                 = 'reverse';
  title('women'); box on
  

  % men 
  sb3 = subplot('Position', sb2.Position + [0.232 0 0 0]); hold on
  im                      = imagesc( cell2mat(men_corr2(2:end,2:end)) ); 
  im.AlphaData            = min(1,max(0,1/6 * log10(10e-2 ./ cell2mat(men_pval2(2:end,2:end))) )); 
  for uvix = 1:size(men_corr2,1),  plot([uvix uvix]+.5,[1 22]-.5,'Color',[0.9 0.9 0.9]); end
  for uvix = 1:size(men_corr2,2), plot([1 22]-.5,[uvix uvix]+.5,'Color',[0.9 0.9 0.9]); end
  for uvix = use_vals{uvi,5}, plot([uvix uvix]+.5,[1 size(both_corr2,1)]-.5,'Color',[0 0 0]); end
  for uvix = use_vals{uvi,4}, plot([1 size(both_corr2,2)]-.5,[uvix uvix]+.5,'Color',[0 0 0]); end
  xlim([0.5 size(both_corr2,2)-.5]); ylim([0.5 size(both_corr2,1)-.5]); caxis([-1 1]);
  ax                      = gca; 
  ax.TickLength           = [0 0];
  ax.XTick                = 1:numel( both_corr2(1,2:end) ); 
  ax.XTickLabel           = cat_io_strrep( both_corr2(1,2:end) ,'_','\_');
  ax.XTickLabelRotation   = 90;
  ax.YTick                = []; 
  ax.FontSize             = fontsize; 
  ax.YDir                 = 'reverse';
  title('men'); box on

 
  % Legend with colortable to code the r-values by color and p-values by
  % transparency.
  subplot('Position',[0.91 sb3.Position(2) 0.045 sb3.Position(4)]); hold on
  im                      = imagesc( repmat( flip( (-1:.1:1)' ,2 ), 1 , 4) );
  im.AlphaData            = repmat( min(1,max(0,1/6 * log10(10e-2 ./  [10e-2 10e-4 10e-6 10e-8]) )), 21 , 1);   
  ax                      = gca; 
  ax.TickLength           = [0 0];
  ax.XTick                = 1:4; 
  ax.XTickLabel           = {'10e-2','10e-4','10e-6','10e-8','10e-10'};
  ax.XTickLabelRotation   = 90;
  ax.YAxisLocation        = 'right';
  ax.YTick                = 1:5:21; 
  ax.YTickLabel           = {'-1.0','-0.5','0.0','0.5','1.0'};
  ax.FontSize             = fontsize*.8;
  for uvix = 1:4,  plot([uvix uvix]+.5,[1 22]-.5,'Color',[0.9 0.9 0.9]); end
  for uvix = 1:21, plot([1 22]-.5,[uvix uvix]+.5,'Color',[0.9 0.9 0.9]); end
  %for uvix = [1 5 6 10 11 15 16 20 21], plot([1 size(both_corr2,2)]-.5,[uvix uvix]+.5,'Color',[.2 .2 .2]); end
  ylim([.5 21.5]); xlim([.5 4.5]); box on
  xlabel('pval'); ylabel('r-val')
  title('cmap','FontSize',fontsize*0.8);
  
 
  % save image
  saveas(gcf,fullfile(boneyEvalDir,sprintf('HBM_2024_table%02d.png', uvi ))); 

end

    