%% Create HBM2024 table
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


fg = figure(3); 
fg.Position(3:4) = [800 600];
fontsize     = 13 * 2.5;
boneyEvalDir = fileparts(which('HBM_boney_scatter.m'));


% load correlation data prepared in R
ukb    = cat_io_csv( fullfile(boneyEvalDir,'boney_data.csv') );

sex    = cell2mat( ukb(2:end,3));

boneyM = [ 10];       % [9:11,  17];
evalM  = [4:8, 12:15, 16:24];
bi = 1; ei = 1;
xys = 1;

%%
scolors = { [ 0    0.5    0.75] , [ 1 .5 0]} ;
for bi = 1:numel(boneyM)
  for ei = 1:numel(evalM)
    %%
    clf(3); 

    % var data
    bmdx = cell2mat( ukb(2:end, boneyM(bi)) ); 
    val  =  ukb(2:end, evalM(ei)) ;
    if any( cellfun(@ischar,val) )
      val( cellfun(@ischar,val) ) = {nan}; 
    end
    valx = cell2mat( val ); 
    
    % var names
    bmdname = strrep( ukb{1,boneyM(bi)}, '_', '\_'); 
    valname = strrep( ukb{1,evalM(ei)} , '_', '\_'); 

    % men
    if xys
      sh = scatter( valx(sex==1), bmdx(sex==1), 's','filled'); 
    else
      sh = scatter( bmdx(sex==1), valx(sex==1), 's','filled'); 
    end
    sh.MarkerEdgeColor = scolors{1}; 
    sh.MarkerFaceColor = scolors{1};
    sh.SizeData = sh.SizeData * 6; 
    sh.MarkerEdgeAlpha = 0.1; 
    sh.MarkerFaceAlpha = 0.3; 
    hold on
    
    % women
    if xys
      sh = scatter( valx(sex==2), bmdx(sex==2), 'o','filled'); 
    else
      sh = scatter( bmdx(sex==2), valx(sex==2), 'o','filled'); 
    end
    sh.MarkerEdgeColor = scolors{2}; 
    sh.MarkerFaceColor = scolors{2}; 
    sh.SizeData = sh.SizeData * 6; 
    sh.MarkerEdgeAlpha = 0.1; 
    sh.MarkerFaceAlpha = 0.3; 
       
    % naming
   % title( sprintf( '%s vs. %s', bmdname , valname ) , 'FontSize', fontsize * 1.1);
    ax = gca; 
    ax.FontSize = fontsize * 0.9;
    ax.LineWidth = ax.LineWidth * 5;
    if xys
      xlabel( valname , 'FontSize', fontsize ); 
      ylabel( bmdname , 'FontSize', fontsize );
    else
      xlabel( bmdname , 'FontSize', fontsize ); 
      ylabel( valname , 'FontSize', fontsize );
    end
    box on; 

    xlim('padded')
    ylim('padded')

    % save image
    saveas(gcf,fullfile(boneyEvalDir, sprintf('HBM_2024_scatter_%s-%s.png', ...
      ukb{1,boneyM(bi)}, ukb{1,evalM(ei)} ))); 

  end
end

