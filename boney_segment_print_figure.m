function boney_segment_print_figure(Vo,Ym,Yc,Ybonemarrow,Si,St,Stm,seg8t,tis,tismri,Po,opt,Affine) % many things 
%print_figure. Print final report figure with volumes and surfaces.

% ToDo
% * V1: nearest interpolation with optimized data ?
% * V2: low fat issue 
% * V2: pure WM segments but jet for bone?
% * V2: countours?
% 
% * S: altas rendering (only highres) > option?
% * S: atlas nan filtering 
%

%#ok<*AGROW> 

  % in case of updates
  seg8t.Affine = Affine; 
  vx_vol = tis.res_vx_vol;
  clear -global st
  global st %#ok<GVMIS> 

  % fontsettings
  crange    = 8; 
  fontsize  = 10; 
  fontcolor = [0 0 0];
  fontname  = 'monospace';

  % colormaps
  clscol = [ lines(7) ; gray(5)];
  clscol = clscol([5,2,1,4,3,10],:); 
  grycol = gray(10);
  tiscol = cat_io_colormaps('BCGWHw',200);
  jetcol = jet(100)/1.5;
  mrkcol = cat_io_colormaps('marks',100)/1.5;
  
  % setup SPM figure
  fg = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',fg)
  spm_figure('Clear',fg);
  if isempty(fg)
    fg = spm_figure('Create','Graphics','visible','on'); 
  end
  colormap gray; 
  

  % main text report box with header (filename)
  ax = axes('Position',[0.01 0.75 0.98 0.245],'Visible','off','Parent',fg);
  text(0,0.99, ['Bone marrow extraction: ' ...
    strrep( strrep( spm_str_manip(Vo.fname,'k80d'),'\','\\'), '_','\_') '       '],...
    'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','tex','Parent',ax);




  % == table 1: full spm classes values - maybe only expert? ==
  notab1 = 1; 
  if opt.expert > 1 && ~notab1
    mgid = find(seg8t.mg > 2/numel(seg8t.lkp));
    FN  = {'lkp','mg','mn','vr'};
    FNn = {sprintf('%d/%d SPM-Classes ',numel(mgid),numel(seg8t.lkp)),'Propotion','Mean','Std'}; 
    FNt = {'%d','%0.2f','%0.2f','%0.2f'};
    % gray table rows
    for fni = 1:numel(FN)
      if fni==1
        text(0.01,0.95-            0.035,repmat(' ',1,1205), ...
          'BackgroundColor',[0.94 0.94 .94],'Parent',ax,'FontSize',fontsize*.6); 
      elseif mod(fni,2)
        text(0.01,0.95-(0.035*fni)-0.035,repmat(' ',1,1205), ...
          'BackgroundColor',[0.96 0.96 .96],'Parent',ax,'FontSize',fontsize*.6); 
      end
    end
    % add table content
    for fni = 1:numel(FN)
      seg8text(fni,1)   = text(0.01, 0.96-(0.045*fni) +0.01*(fni==1), ...
        sprintf('\\bf%s',FNn{fni}) ,'Fontname',fontname,'FontSize',...
        fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax); 
      if fni==1, seg8text(fni,1).FontWeight = 'bold'; norm = 1; else, norm = seg8t.(FN{fni})(2); end
      for lkpi = 1:numel(mgid) %numel(seg8t.lkp)
        seg8text(fni,mgid(lkpi)+1) = text(0.11 + 0.87/numel(mgid) * (lkpi), 0.95-(0.045*fni) +0.01*(fni==1), sprintf( FNt{fni}, seg8t.(FN{fni})(mgid(lkpi))/norm ),...
          'Fontname',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax,'HorizontalAlignment','right');
        if fni==1 
          seg8text(fni,mgid(lkpi)+1).FontWeight = 'bold'; seg8text(fni,mgid(lkpi)+1).Color = clscol(seg8t.lkp(mgid(lkpi)),:);  
          seg8text(fni,mgid(lkpi)+1).String = [seg8text(fni,mgid(lkpi)+1).String sprintf('(%d)',sum(seg8t.lkp == seg8t.lkp(mgid(lkpi))) )];
        end
        if fni==2, seg8text(fni,mgid(lkpi)+1).Color = grycol( round(8 - 7*seg8t.mg(mgid(lkpi))) ,:);  end
        if fni==3 
          try
            col = tiscol( max(0, min(200,round(100*seg8t.mn(mgid(lkpi))/max(seg8t.mn(2))))) ,:); 
            seg8text(fni,mgid(lkpi)+1).Color = max(0,min(1,col - max(.3,min(sum(col,2)-2)/4)));  
          end
        end
        if fni==4, seg8text(fni,mgid(lkpi)+1).Color = grycol( min(10,max(1,round(8 - 7 * seg8t.vr(mgid(lkpi)) / max(seg8t.vr(1:3))))) ,:);  end
      end
    end

    table1offset = 0.27;
  else 
    table1offset = 0.02; 
  end



  % == table 2: spm classes values ==
  tis.clsn = {'1-GM','2-WM','3-CSF','4-bone','5-head','6-BG'};
  for ci = 1:min(6,max(seg8t.lkp))
    % indicate cases were a Gaussian is used onyl for a small volume, what
    % is pointing to a to large number of Gaussian
    % ... would have to consider also the difference between Gaussians
    if sum(seg8t.lkp==ci) > max(1,sum(seg8t.mg(seg8t.lkp==ci)'<.2 | abs(gradient(seg8t.mn(seg8t.lkp==ci) ./ tis.WMth)) < .2)) 
      mark = '\color[rgb]{0,0.5,1}';
    else
      mark = ''; 
    end
    tis.clsG{ci} = sprintf('%d %s(%d)', sum(seg8t.lkp==ci), mark, ...
      max(1,sum(seg8t.mg(seg8t.lkp==ci)'<.2 | abs(gradient(seg8t.mn(seg8t.lkp==ci) ./ tis.WMth)) < .2)) );  
  end
  % #####################
  % * add QC value for tissue ? 
  % *
  % #####################
  % density ~ volume and so I use only one measure
  if seg8t.isCTseg && ~opt.normCT 
    tis.seg8x = tismri.Tth; 
    FN  = {'clsn','seg8x','vol','volr'};                                  % 'den','clsG','seg8nv',
    FNn = {'TPM class','Med.Int.','Volume mm3','Volume (%TIV)'};     % 'Volume density','n-Class(well)','Main-std',
    FNt = {'%s','%0.0f','%0.0f','%0.1f'};                                 % '%0.0f','%s','%0.3f',
  else
    FN  = {'clsn','seg8n','vol','volr'};                                  % 'den','clsG','seg8nv',
    FNn = {'TPM class','Norm.Med.Int.','Volume mm3','Volume (%TIV)'};     % 'Volume density','n-Class(well)','Main-std',
    FNt = {'%s','%0.3f','%0.0f','%0.1f'};                                 % '%0.0f','%s','%0.3f',
  end
  if isempty(tismri) || ~isfield(tismri,'vol')
    FN = FN(1:3); FNn = FNn(1:3); FNt = FNt(1:3);
  else
    tis.vol  = tismri.vol; 
    tis.den  = tismri.den;
    tis.volr = tismri.vol ./ sum(tismri.vol(1:3)) * 100; 
    % do not show BG volumes 
    tis.vol(6)  = nan; 
    tis.volr(6) = nan;
  end
  for fni = 1:numel(FN)
    if fni==1
      text(0.01,0.95-table1offset-            0.035,repmat(' ',1,1000 + 205),...
        'BackgroundColor',[0.94 0.94 .94],'Parent',ax,'FontSize',fontsize*.6); 
    elseif mod(fni,2)
      text(0.01,0.95-table1offset-(0.045*fni)-0.0,repmat(' ',1,1000 + round(205*(.6/.4))),...
        'BackgroundColor',[0.96 0.96 .96],'Parent',ax,'FontSize',fontsize*.4); 
    end
  end
  for fni = 1:numel(FN)
    tistext(fni,1)   = text(0.01              , 0.95-table1offset-(0.045*fni) +0.01*(fni==1), ...
      sprintf('\\bf%s',FNn{fni}) ,'Fontname',fontname,'FontSize',...
      fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
    if fni==1, tistext(fni,1).FontWeight = 'bold'; end
    for lkpi=numel(tis.clsn):-1:1 % down to avoid overlap
      if iscell(tis.(FN{fni})), tmpval = tis.(FN{fni}){lkpi}; else, tmpval = tis.(FN{fni})(lkpi); end
      if ~isnan(tmpval)
        tistext(fni,lkpi+1) = text(0.12 + 0.082*(lkpi), 0.95-table1offset-(0.045*fni) +0.01*(fni==1), sprintf( FNt{fni}, tmpval ),...
          'Fontname',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax,'HorizontalAlignment','right');
        if fni==1, tistext(fni,lkpi+1).FontWeight = 'bold'; tistext(fni,lkpi+1).Color = clscol(lkpi,:);  end
      end
      %if fni==3, tistext(fni,lkpi+1).Color = grycol( max(1,min(10,round(8 - tmpval))) ,:);  end
    end
    % ###################
    % * add colorrating for tissues 
    % * density vs. volume => error?
    % ###########
  end



  % == table 3: processing and qc parameter? 
  % ... processing options (method-SPM,Bone,Surf)
  % - processing [SPM,CTseg]
  % - number spm-peaks [1 1 2 3 4 2]
  % - method [mat,SPM,ref]
  % 

  % == table 4: imaging and bone parameters
  highBGs = {'low','high','high(MT)','mid'}; %twcol = [0 0 0; 0 0 0.5; 0 .5 0; .5 0 0]; 
  tis.highBGn = highBGs{tis.highBG +1}; % ##########################
  tis.res_RES = mean(vx_vol.^2 )^.5;
  if ~isempty(St)
    tis.iBone   = cat_stat_nanmean(Si.facevertexcdata(Si.facevertexcdata<100));
    tis.thBone  = cat_stat_nanmean(St.facevertexcdata(St.facevertexcdata<100));
  end
  tis.medBone0 = 0; 
  try,  tis.headthickmn = tismri.headthickmn; end
% ################## 
%  * add skull-Stripping and defacing value later 
%  * better definition of low fat intensity rating 
%  * do QC with CNR rating? 
% ################################# 
  FN{1}  = {'weightingn', 'res_RES' , 'highBGn' , 'headFatTypen', 'boneIntTypen'}; 
  FNn{1} = {'Tw'        , 'RES'     , 'BGtype'  , 'fatint',       'fatbone'     };
  FNt{1} = {'%s'        , '%0.2f'   , '%s'      , '%s',           '%s'          }; 

  if seg8t.isCTseg
    tis.minBone = tismri.iBonemn3(1);
    tis.medBone = tismri.iBonemn3(2);
    tis.maxBone = tismri.iBonemn3(3);
    FN{2}  = {'minBone'   , 'medBone' , 'maxBone' , 'thBone'  , 'headthickmn'};
    FNn{2} = {'lBone*'    , 'mBone*'  , 'hBone*'  , 'BthBone' , 'thhead'};
    FNt{2} = {'%0.0f'     , '%0.0f'   , '%0.0f'   , '%0.1f mm', '%0.1f mm'};
  else
    FN{2}  = {'minBone'   , 'medBone' , 'maxBone' , 'medBone0', 'iBone'       , 'thBone'  , 'headthickmn'};
    FNn{2} = {'lBone'     , 'mBone'   , 'hBone'   , 'medBone0', 'iBone'       , 'BthBone' , 'thhead'};
    FNt{2} = {'%0.3f'     , '%0.3f'   , '%0.3f'   , '%0.3f'   , '%0.3f'       , '%0.1f mm', '%0.1f mm'};
  end
  if isempty(St)
    FN{2}(end-1:end) = []; FNn{2}(end-1:end) = []; FNt{2}(end-1:end) = []; 
  end
  for fnj = 1:2
    text(0.01,0.64 - (fnj-1)*0.12 - table1offset,repmat(' ',1,1000 + 120),'BackgroundColor',[0.94 0.94 .94],'Parent',ax,'FontSize',fontsize*.6); 
    for fni = numel(FN{fnj}):-1:1
      segtext(fni,lkpi)   = text(0.03 + 0.075*(fni-1), 0.64-(fnj-1)*0.12-table1offset-(0.00) , ...
        sprintf('\\bf%s',FNn{fnj}{fni}) ,'Fontname',fontname,'FontSize',...
        fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax,'FontWeight','bold','HorizontalAlignment','center');
      if ~isempty( FNt{fnj}{fni} )
        segtext(fni,lkpi+1) = text(0.03 + 0.075*(fni-1), 0.64-(fnj-1)*0.12-table1offset-(0.05), sprintf( FNt{fnj}{fni},  tis.(FN{fnj}{fni}) ),...
          'Fontname',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','HorizontalAlignment','center','Parent',ax);
      end
%      if fni==1, segtext(fni,lkpi+1).Color = twcol(  max(1,min(size(twcol,1),tis.(FN{fnj}{fni}(1:end-1)))) , :); end
      if fni==2, segtext(fni,lkpi+1).Color = mrkcol( max(1,min(100,round( 100 * (tis.(FN{fnj}{fni})+1) / 4 ))) , :); end % RES
      if fnj == 1
        if fni >2 && fni < 6
          switch segtext(fni,lkpi+1).String
            case 'low',  segtext(fni,lkpi+1).Color = [ .0 .0 1.0]; 
            case 'mid',  segtext(fni,lkpi+1).Color = [ .6 .5  .0]; 
            case 'high', segtext(fni,lkpi+1).Color = [ .8 .2  .0]; 
            otherwise,   segtext(fni,lkpi+1).Color = [1.0 .0  .0]; 
          end
        end
      else
        if seg8t.isCTseg
          if fni > 0 && fni<4, segtext(fni,lkpi+1).Color = jetcol( max(1,min(100,round(0.07 * tis.(FN{fnj}{fni})))),:);  end
          if fni > 3,          segtext(fni,lkpi+1).Color = jetcol( max(1,min(100,round(5 * tis.(FN{fnj}{fni})))),:);  end
        else
          if fni > 4,          segtext(fni,lkpi+1).Color = jetcol( max(1,min(100,round(5 * tis.(FN{fnj}{fni})))),:);  end
        end
      end
    end
  end



  %% == error messages 
  % ... show only in case of no refinement or as note?
  if 0 
    wstr   = ''; % empty line
    wname  = {'note','warning','alert'}; 
    wcolor = [0 0 1; 1 0.5 0; 0.8 0 0]; 
    for wi=2:-1:0
      warn = cat_io_addwarning(wi); 
      wn   = numel(warn);
      if wn
        if cat_get_defaults('extopts.expertgui') > -wi
          msg  =  strrep(warn(1).identifier,'_','\_') ; 
          for wmi=2:wn
            msg = [msg ', ' strrep( warn(wmi).identifier,'_','\_') ];
          end
          
          if wn~=1, wnamepl='s'; else, wnamepl=''; end
          wstr = [wstr ...
            sprintf('\\color[rgb]{%0.1f %0.1f %0.1f}\\bf{%d %s%s} (%s)',...
              wcolor(wi+1,:), wn,wname{wi+1},wnamepl, msg)]; 
        end
      end
    end
    text(0.01,0.1,wstr,'Interpreter','tex'); 
  end
 

  % method 0 only read the header and we can print only the basic values/tables 
  if 0 % opt.opts.bmethod == 0
    %%    print(fg, '-djpeg', Po.report); 
    figfs10 = [ findobj(fg,'FontSize',10); findobj(fg,'FontSize',9); findobj(fg,'FontSize',11); findobj(fg,'FontSize',8)]; 
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .8; end; end %#ok<TRYNC> 
    saveas(fg,Po.report,'png')
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .8; end; end %#ok<TRYNC> 
    return; 
  end


  % histogram 
  if opt.output.report > 1 && opt.opts.bmethod > 0 
    %%
    ax2 = axes('Position',[0.61 0.855-table1offset/3 0.40 0.13],'Parent',fg,'Color',[1 1 1]); hold on;
    ax2.YColor = [1 1 1]; ax2.XColor = [1 1 1]; ax2.YTick = []; ax2.XTick = []; 
    ax2 = axes('Position',[0.62 0.855-table1offset/3 0.37 0.13],'Parent',fg,'Color',[1 1 1]); hold on;
    if seg8t.isCTseg && ~opt.normCT 
      for ci = numel(Yc)-1:-1:1
        % (tismri.Tth(2) ./ tis.WMth) correction to print the SPM 
        hhst(ci) = histogram(ax2,Ym( Yc{ci}(:)>.5 & Ym(:)>-1000 & Ym(:)<2000) ,-1000:10:2000, ...
          'LineStyle','none','FaceColor',clscol(ci,:));
        hstmax(ci) = max( hhst(ci).Values )  / 3;
      end
      xlim([-300 1800]); %title('normalized intensities');
    elseif opt.opts.bmethod > 0 
      for ci = numel(Yc)-1:-1:1
        % (tismri.Tth(2) ./ tis.WMth) correction to print the SPM 
        hhst(ci) = histogram(ax2,Ym( Yc{ci}(:)>.5 & Ym(:)>0 & Ym(:)<2) ,0:0.01:2, ...
          'LineStyle','none','FaceColor',clscol(ci,:));
        hstmax(ci) = max( hhst(ci).Values ) ;
      end
      xlim([0 2]); ax2.XTick = 0:1/6:2; % grid on;  %title('normalized intensities');
    else
      % not loaded
      hhst(4) = histogram(ax2,Ybonemarrow(Ybonemarrow(:)>0),0:0.01:2, ...
        'LineStyle','none','FaceColor',clscol(ci,:));
      hstmax(4) = max( hhst(4).Values ) ;
      xlim([0 2]); ax2.XTick = 0:1/6:2; % grid on;  %title('normalized intensities');
    end
    %%
    ylim([0 max(hstmax(1:5)*1.3 )])
    box on; 
    if ~(seg8t.isCTseg && ~opt.normCT)
      ax2.FontSize = fontsize * 0.85; 
      ax2.XTickLabelRotation = 0; 
      ax2.XTickLabel = {'0','','','0.5','','','1','','','1.5','','','2'};
      xlabel('normalized intensities with normalized SPM classes'); 
    else
      xlabel('CT intensities of SPM classes'); 
    end
    ax2.YTick = []; 
    if seg8t.isCTseg && ~opt.normCT 
      for ci = 1:3
        pl = plot([tismri.Tth(ci) tismri.Tth(ci)],[0 max(ylim)]); 
        pl.Color     = [ clscol(min(6,max(seg8t.lkp(ci))),:)];
      end
      for ci = 1:numel(tismri.iBonemn)
        pl = plot( repmat(tismri.iBonemn(ci),1,2) ,[0 max(ylim)]); 
        pl.Color = [ clscol(4,:)];
      end
      pl = plot( repmat( tis.report.fat ,1,2) ,[0 max(ylim)]); 
      pl.Color = [ clscol(5,:)];
      pl = plot( repmat( tis.report.muscle ,1,2) ,[0 max(ylim)]); 
      pl.Color = [ clscol(5,:)];
    else
      for ci = 1:numel(seg8t.mn)
        pl = plot([seg8t.mn(ci) seg8t.mn(ci)] / tis.WMth,[0 max(ylim)]); 
      end
      pl.Color     = [ clscol(min(6,max(seg8t.lkp(ci))),:) min(1,max(0,seg8t.mg(ci) * .7*sum(seg8t.lkp(ci)==seg8t.lkp).^.5))];
      pl.LineWidth = seg8t.mg(ci); %.^2 * 2; 
    end
    % backup?  plot([1 1],[0 max(ylim)],'-','Color',[.7 .7 .7]);
    if opt.opts.bmethod>0
      if numel(Yc)==6, lg = legend(flip({'GM','WM','CSF','bone','head'}),'box','off'); lg.Position(1) = .895; end 
    else
      legend(flip({'bone'}),'box','off'); 
    end
  end  
% ######################################
% horizontal boxplot for bones only with new measure as bars 
% ######################################



  %% == images ==
  if opt.output.report > 1
    spm_orthviews('Reset')
    pos = {[0.008 0.375 0.486 0.35]; [0.506 0.375 0.486 0.35]};
    % T1 + SPM segmentation
    if tis.weighting < 0 
      colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(crange+2)),1); repmat([ 0 .3 .6],round(59/(crange+2)),1); jet(59 - 2*round(59/(crange+2)))]);
      V0 = Vo; V0.dat(:,:,:) = single(.3*Yc{3} + .75*Yc{1} + 1*Yc{2} + 1.2*Yc{4} + 0.5*Yc{5} ); V0.dt(1) = 16;
      V0.pinfo  = repmat([1;0],1,size(Ym,3));
    elseif ~isempty(Ym)
    %colormap( [zeros(1,3); jet(59 - 2*round(59/(crange+2))); 0.5*ones(round(59/(crange+2)),3); ones(round(59/(crange+2)),3)])
      colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(crange+2)),1); repmat([ 0 .3 .6],round(59/(crange+2)),1); jet(59 - 2*round(59/(crange+2)))]);
      V0 = Vo; V0.dat(:,:,:) = single(Ym); V0.dt(1) = 16;
      V0.pinfo  = repmat([1;0],1,size(Ym,3));
    else
      colormap gray; 
      V0       = seg8t.image;
      Vc4      = spm_vol(spm_file(seg8t.image.fname,'prefix','c4'));
      Vc4.mat  = seg8t.Affine  * Vc4.mat; 
    end
    V0.mat    = seg8t.Affine * seg8t.image.mat; % Vo.mat;   \seg8t.tpm(1).mat
    hh0       = spm_orthviews('Image',spm_vol(fullfile(spm('dir'),'tpm','TPM.nii,1')),pos{1}); % avoid problems with high resolutions
    spm_orthviews('window',hh0,[0 10000]); % just make it black
    spm_orthviews('BB', [-85 -120 -90; 85 95 105]); % this has to be set in the low-resolution image
    hh0       = spm_orthviews('Image',V0,pos{1}); % add correct image after the other settings! 
    spm_orthviews('Reposition',[-25 0 0]);
    if tis.weighting < 0 
      spm_orthviews('Caption',hh0,sprintf('%s','CTseg'));
    else
      spm_orthviews('Caption',hh0,sprintf('%s',tis.weightingn));
    end
    if ~isempty(Ym)
      spm_orthviews('window',hh0,[0 1.2]); %
    else
      spm_orthviews('addtruecolourimage',hh0,Vc4,[0 0 0; 0.5 0.5 0.5; 1 0 0],0.4,1,0.1); 
      spm_orthviews('redraw');

      orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'),'Type','Image','Tag',''),'parent');
      orthviewlegend.YTick      = [0 .5 1];
      orthviewlegend.YTickLabel = {'BG','~Bone','Bone'};
      orthviewlegend.FontSize   =  orthviewlegend.FontSize * 0.85;
      
      try % does not work in headless mode without java
        figfs10 = [ findobj(fg,'FontSize',10); findobj(fg,'FontSize',9); findobj(fg,'FontSize',11); ...
          findobj(fg,'FontSize',8); findobj(fg,'FontSize',fontsize*.6);  findobj(fg,'FontSize',fontsize*.4); ]; 
        for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .8; end; end %#ok<TRYNC> 
        saveas(fg,Po.report,'png')
        for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .8; end; end %#ok<TRYNC> 
      catch
        cprintf('err','Cannot print report.')
      end
      
      return
    end
     %% 
    %  spm_orthviews('BB', [-90 -120 -80; 90 90 105]); %[[-78 -112 -70];[78 76 85]];
    %spm_orthviews('BB', [-90 -120 -80; 90 90 100] ./ [vx_vol;vx_vol]); % ######################## this is not corret but works ?! ... 
  
   
    
    % print cls overlay?
    if tis.weighting >= 0
      Vp0        = Vo; 
      Vp0.dt(1)  = 16;
      Vp0.mat    = seg8t.Affine * Vo.mat; 
      Vp0.pinfo  = repmat([1;0],1,size(Ym,3));
      switch 2
        case 1
          Vp0.dat(:,:,:) = single( (Yc{4}>.5) + (2*(Yc{5}>.5) + 0.5*(Yc{5}>.5)) ); 
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 1 0 0; 1 1 0],0.4,2,0.1); 
        case 2
          Vp0.dat(:,:,:) = single( (Yc{4}>.5) + (Yc{4}>.5 & Ym>.7 ) + (3*(Yc{5}>.5)) + (Yc{5}>.5 & (Ym>.9)) ); 
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 0 0.5 1; .5 0.8 1; .9 .5 0.5; 1 1 0 ],0.5,4,0); 
        case 5
          Vp0.dat(:,:,:) = single(2*(Yc{1}>.5) + 3*(Yc{2}>.5) + (Yc{3}>.5) + 4*(Yc{4}>.5) + 5*(Yc{5}>.5)); %.*(Ybone>.5) + 4*Yc{4}.*(Ybone>.5);
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 0 0 1; 0 .8 0; 1 1 1; 1 0 0; 1 1 0],0.4,5,0); 
      end
    end
  
    %% print bone marrow
    if tis.weighting < 0
      V1 = Vo; V1.dat(:,:,:) = single((1000 + Ym)/2300); V1.dt(1) = 16;
    else
      Ybonemarrow2 = Ybonemarrow; Ybonemarrow2(isnan(Ybonemarrow2)) = 0; 
      V1 = Vo; V1.dat(:,:,:) = min(crange,Ybonemarrow2) + 2*(Ybonemarrow2>0) + single(Yc{1}>.5) + single(Yc{5}>.5) + 2*single(Yc{2}>0.5);
    end
    V1.pinfo  = repmat([1;0],1,size(Ym,3));
    V1.mat    = seg8t.Affine * V1.mat; 
    V1.dt(1)  = 16;
    hh1       = spm_orthviews('Image',V1,pos{2}); 
    spm_orthviews('Interp',1);
    if tis.weighting < 0
      spm_orthviews('window',hh1,[0 1.2]);
      spm_orthviews('Caption',hh1,sprintf('%s',tis.weightingn));
    else
      spm_orthviews('window',hh1,[0 crange+2]);
      spm_orthviews('Caption',hh1,'bonemarrow + GM/HD + 2*WM');
    end
    spm_orthviews('Reposition',[-25 0 0]);
    spm_orthviews('redraw'); 
   
    
    % this is replaced by spm_orthviews('redraw'); 
    if exist('orthviewlegend','var')
      if opt.opts.bmethod > 0
        orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'),'Type','Image','Tag',''),'parent');
        orthviewlegend.YTickLabel = {'BG','lBone','hBone','muscle','fat'};
        orthviewlegend.FontSize   =  orthviewlegend.FontSize * 0.85;
      else
        orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'),'Type','Image','Tag',''),'parent');
        orthviewlegend.YTick      = {0 1};
        orthviewlegend.YTickLabel = {'BG','Bone'};
        orthviewlegend.FontSize   =  orthviewlegend.FontSize * 0.85;
      end
    end

    %% bone surface 
    if isfield(Po,'central') && ~isempty(St)
      for idi = 2:3 %1:2
        spm_orthviews('AddContext',idi);  
        spm_ov_mesh('display',idi,Po.central );
        % apply affine transformation
        V = (seg8t.Affine * ([st.vols{idi}.mesh.meshes(end).vertices,...
             ones(size(st.vols{idi}.mesh.meshes(end).vertices,1),1)])' )';
        V(:,4) = [];
        st.vols{idi}.mesh.meshes = subsasgn(st.vols{idi}.mesh.meshes,struct('subs','vertices','type','.'),single(V));
        % change line style
        hM = findobj(st.vols{idi}.ax{1}.cm,'Label','Mesh');
        UD = get(hM,'UserData');  
        UD.width = 1;
        if tis.weighting == -1
          UD.style = {'w-'};  
        else
          if idi==2, UD.style = {'r-'}; else, if tis.boneIntType>1, UD.style = {'w-'}; else, UD.style = {'r-'}; end; end
          %UD.style = {'r-'};  
        end
        set(hM,'UserData',UD); clear UD hM
        warning('off','MATLAB:subscripting:noSubscriptsSpecified');
        spm_ov_mesh('redraw',idi);
      end
    end
  end
 




  %% == surface ==
  clear h; 
  if ~isempty(Stm) && opt.output.report > 1
    % thickness map
    hCS{1} = subplot('Position',[0.015 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{1} = 'l';
    hCS{2} = subplot('Position',[0.255 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{2} = 'r';
    hCS{3} = subplot('Position',[0.015 0.01 0.23 0.19],'Parent',fg,'visible','off');  sview{3} = 't';
    hCS{4} = subplot('Position',[0.255 0.06 0.23 0.14],'Parent',fg,'visible','off');  sview{4} = 'p';
    hCS{5} = subplot('Position',[0.255 0.03 0.23 0.02],'Parent',fg,'visible','off'); 
    
    imat = spm_imatrix(seg8t.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
    V = (Rigid * ([Stm.vertices, ones(size(Stm.vertices,1),1)])' )'; V(:,4) = []; Stm.vertices = V;

    for ci = 1:4
      h{ci} = cat_surf_render2(Stm,'parent',hCS{ci}); 
      cat_surf_render2('Clim',h{ci},[0 20]); 
      switch lower(sview{ci})
        case {'r'},  cat_surf_render('view',h{ci},[  90   0]); 
        case {'l'},  cat_surf_render('view',h{ci},[ -90   0]);  
        case {'t'},  cat_surf_render('view',h{ci},[   0  90]); 
        case {'p'},  cat_surf_render('view',h{ci},[   0   0]); 
      end
    end
    
    cb = cat_surf_render2('Colorbar',h{4});
    cb.colourbar.Location = 'South';
    cb.colourbar.Position = [0.275 0.015 0.2 0.005];
    text(0.5,0.015,'Bone thickness (tBone)','fontsize',fontsize-1,'FontName',fontname,...
        'HorizontalAlignment','center','Parent',hCS{5}); 

  end

  if ~isempty(Si) && opt.output.report > 1
    hCS{1} = subplot('Position',[0.5 + 0.015 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{1} = 'l';
    hCS{2} = subplot('Position',[0.5 + 0.255 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{2} = 'r';
    hCS{3} = subplot('Position',[0.5 + 0.015 0.01 0.23 0.19],'Parent',fg,'visible','off');  sview{3} = 't';
    hCS{4} = subplot('Position',[0.5 + 0.255 0.06 0.23 0.14],'Parent',fg,'visible','off');  sview{4} = 'p';
    hCS{5} = subplot('Position',[0.5 + 0.255 0.03 0.23 0.02],'Parent',fg,'visible','off'); 

    imat = spm_imatrix(seg8t.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
    V = (Rigid * ([Si.vertices, ones(size(Si.vertices,1),1)])' )'; V(:,4) = []; Si.vertices = V;

    for ci = 1:4
      h{ci} = cat_surf_render2(Si,'parent',hCS{ci}); 
      if tis.weighting < 0 
        cat_surf_render2('Clim',h{ci},[500 1500]); 
      else
        cat_surf_render2('Clim',h{ci},[0 crange]); 
      end
      switch lower(sview{ci})
        case {'r'},  cat_surf_render('view',h{ci},[  90   0]); 
        case {'l'},  cat_surf_render('view',h{ci},[ -90   0]);  
        case {'t'},  cat_surf_render('view',h{ci},[   0  90]); 
        case {'p'},  cat_surf_render('view',h{ci},[   0   0]); 
      end
    end

    if 1
      cb = cat_surf_render2('Colorbar',h{4});
      cb.colourbar.Location = 'South';
      cb.colourbar.Position = [0.5 + 0.275 0.015 0.2 0.005];
      text(0.5,0.015,'Bone intensity (iBone)','fontsize',fontsize-1,'FontName',fontname,...
        'HorizontalAlignment','center','Parent',hCS{5}); 

    else
      %% colormap
      if strcmpi(spm_check_version,'octave')
        axes('Position',[0.965 0.03 0.01 0.28],'Parent',fg); image(flip(121:1:120+surfcolors)','Parent',cc{4});
        cc{4} = gca; 
      else
        cc{4} = axes('Position',[0.965 0.03 0.01 0.28],'Parent',fg); image(flipud(121:1:120+surfcolors)','Parent',cc{4});
      end
      
      %% histogram line
      if strcmpi(spm_check_version,'octave')
        axes('Position',[0.936 0.03 0.03 0.28],'Parent',fg,'Visible', 'off','tag', 'cat_surf_results_hist', ...
        'xcolor',fontcolor,'ycolor',fontcolor);
        cc{5} = gca; 
      else
        cc{5} = axes('Position',[0.936 0.03 0.03 0.28],'Parent',fg,'Visible', 'off','tag', 'cat_surf_results_hist', ...
          'xcolor',fontcolor,'ycolor',fontcolor);
      end
      side  = hSD{1}.cdata;
      [d,h] = hist( side(~isinf(side(:)) & ~isnan(side(:)) &  side(:)<6 & side(:)>0) ,  hrange); %#ok<HIST> 
      d = d./numel(side);
      d = d./max(d);
      text(cc{5},'test')
      
      % print histogram
      hold(cc{5},'on');  
      for bi = 1:numel(d)
        b(bi) = barh(cc{5},h(bi),-d(bi),boxwidth); 
        set(b(bi),'Facecolor',cmap3(bi,:),'Edgecolor',fontcolor); 
      end
      ylim([0,20]); xlim([-1 0]);
    end
  end


  %
  
  colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(crange+2)),1); repmat([ 0 .3 .6],round(59/(crange+2)),1); jet(59 - 2*round(59/(crange+2)))]);
    

  try % does not work in headless mode without java
    figfs10 = [ findobj(fg,'FontSize',fontsize+1); findobj(fg,'FontSize',fontsize); findobj(fg,'FontSize',fontsize-1);  ...
      findobj(fg,'FontSize',fontsize*0.85); findobj(fg,'FontSize',fontsize*.6);  findobj(fg,'FontSize',fontsize*.4); ]; 
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .75; end; end %#ok<TRYNC> 
    saveas(fg,Po.report,'png')
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .75; end; end %#ok<TRYNC> 
  end

end
%=======================================================================
function [Theader,Tline,Tavg, Cheader, MAfn, matm, mmatm]   = prepare_print(P,opt)
%prepare_print. Create table elements for the in-line report.
  if 0
    MAfn    = {'Tw','Thbg','Tfat','Tbone','Tres','Tcnr',      'BG',  'CSF','GM',   'muscle','fat',   'bone','marrow',   'MBI' , 'BT' , 'HDT' }; % 'Pll', Tw, Tfat, Thg, Tcnr
    MAffn   = {'s' ,'s'   ,'s'   ,'s'    ,'f'   ,'f'   ,      'f' ,  'f'  ,'f' ,   'f'     ,'f'  ,   'f'   ,'f'     ,   'f'   , 'f'  , 'f'   }; % 'e'  ,
    MAfsep  = [0   ,0     ,0     ,0      ,0     ,0     ,      1   ,  0    ,0   ,   1       ,0    ,   1     ,0       ,   0     , 0    , 0     ]; % 0    ,
  else
    MAfn    = {'Tw','Thbg','Tfat','Tbone','Tres','Tcnr',      'BG',  'CSF','GM',   'Tskull','Thead',  'MED' , 'MBI'  , 'BT' , 'HDT' }; % 'Pll', Tw, Tfat, Thg, Tcnr
    MAffn   = {'s' ,'s'   ,'s'   ,'s'    ,'f'   ,'f'   ,      'f' ,  'f'  ,'f' ,   'f'     ,'f'    ,  'f'   , 'f'    , 'f'  , 'f'   }; % 'e'  ,
    MAfsep  = [0   ,0     ,0     ,0      ,0     ,0     ,      1   ,  0    ,0   ,   1       ,0      ,  1     ,  0     , 0    , 0 ]; % 0    ,
  end
  
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',opt.snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',opt.snspace(1)-8);
  Tavg    = sprintf('%%%ds:',opt.snspace(1)-1);
  for fi = 1:numel(MAfn)
    Cheader = [Cheader MAfn{fi}]; 
    AMfni   = strrep( MAfn{fi} ,'_','');
    if MAfsep(fi), sep = ' |'; else, sep = ''; end
    Theader = sprintf(sprintf('%%s%%s%%%ds' ,opt.snspace(2)),Theader, sep, AMfni ); 
    switch MAffn{fi}
      case {'s','d'}
        Tline   = sprintf('%s%s%%%d%s', Tline , sep, opt.snspace(2), MAffn{fi});
        Tavg    = sprintf('%s%s%%%d%s', Tavg  , sep, opt.snspace(2), MAffn{fi});
      otherwise
        Tline   = sprintf('%s%s%%%d.%d%s',Tline , sep, opt.snspace(2), opt.snspace(3) .* (MAffn{fi}=='f'),MAffn{fi});
        Tavg    = sprintf('%s%s%%%d.%d%s',Tavg  , sep, opt.snspace(2), opt.snspace(3) .* (MAffn{fi}=='f'),MAffn{fi});
    end
  end
  % add time 
  Tline   = sprintf('%s%%4.0fs',Tline);
  
  % main result table
  matm    = num2cell(nan(numel(P),numel(MAfn)));
  mmatm   = 10.5*ones(numel(P),numel(MAfn));
  
  % print title
  if opt.opts.verb
    fprintf('\nBone Preprocessing:\n');
    if opt.opts.verb 
      methodstr = {'SPMmat8','MedBoneSurf','rMedBoneSurf'};
      reportstr = {'Table','Table + Volumes','Table + Volumes + Surfaces'};
      fprintf('  Method:   %d (%s)\n', opt.opts.bmethod, methodstr{opt.opts.bmethod+1});
      fprintf('  Report:   %d (%s)\n', opt.output.report, reportstr{opt.output.report});
      if opt.opts.bmethod>0
        fprintf('  Reduce:   %d (%s)\n', opt.reduce, reportstr{opt.reduce});
      end
      % parameter ######################
      %  - method id + name it 
      %  - report option + name it (only tables, imags, full)
      %  - warnings/information for atypical settings
      % ################################
    end
    fprintf('\n%s\n%s\n',  Theader,repmat('-',size(Theader)));  
  end
end
