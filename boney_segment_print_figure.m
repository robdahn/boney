function boney_segment_print_figure(Vo,Ym,Yc,Ybonemarrow, Si,St, out,job,Affine) % many things 
%print_figure. Print final report figure with volumes and surfaces.
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

% ToDo
% * V1: nearest interpolation with optimized data ?
% * V2: low fat tissue 
% * V2: pure WM segments but jet for bone?
% * V2: contours?
% 
% * S: atlas rendering (only highres) > option?
% * S: atlas nan filtering 
%

  %#ok<*AGROW,*TRYNC> 
  %#ok<*HIST>   % have to use hist function in case of octave
  %#ok<*GVMIS>  % this is to get the spm_orthviews "st" variable to control the volume print


  % in case of updates
  out.spm8.Affine = Affine; 
  

  % == setup of SPM figure and print options == 
  popts = prepareReport(out.P.org,out.P.mrirdir);
  
  
  % == Table 1: full spm classes values ==
  %  - this includes a lot of details and is only useful for crazy developers 
  popts.table1offset = printCLStable(job,out,popts,1); % last value = don't print 


  % == Table 2: spm classes values ==
  %  - this one is quite nice
  %  - ###### add tissue QC paraemters? #####
  printTPMtable(job,out,popts);


  % == Table 3: print processing parameter ==
  %  - processing       [SPM,CAT,CTseg]
  %  - method           [mat,SPM,ref]
  %  - number spm-peaks [1 1 2 3 4 2]
  

  % == Table 4: imaging and bone parameters ==
  printMainTable(job,out,popts,Si,St);


  % == error messages ==
  %  - show only in case of no refinement or as note?
  % printErrors
 

 
  % == histogram ==
  %  - not working for Octave so far because I used the newer matlab histogram 
  %    function that is not available 
  printHistogram(Ym,Yc,Ybonemarrow,job,out,popts);


% ######################################
% horizontal boxplot for bones only with new measure as bars 
% ######################################


  %% == images ==
  %  - plot volumes (and draw surface)
  printVolumes(out.P,Vo,Ym,Yc,Ybonemarrow,St,job,out,popts)
  

  % == surface ==
  %  - plot surfaces with bone thickness (with atlas lines) and bone marrow intensity 
  printSurfaces(St,Si,job,out,popts)

end
function popts = prepareReport(fname,resdir)
%prepareReport. Setup of SPM figure and print options.

  % fontsettings
  popts.crange    = 8; 
  popts.fontsize  = 10; 
  popts.fontcolor = [0 0 0];
  popts.fontname  = 'monospace';
  
  % colormaps
  popts.clscol    = [ lines(7) ; gray(8)];
  popts.clscol    = popts.clscol([5,2,1,4,3,13],:); 
  popts.grycol    = gray(10);
  popts.tiscol    = cat_io_colormaps('BCGWHw',200);
  popts.jetcol    = jet(100)/1.5;
  popts.mrkcol    = cat_io_colormaps('marks',100)/1.5;
  
  % setup SPM figure
  popts.fig = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',popts.fig)
  spm_figure('Clear',popts.fig);
  if isempty(popts.fig)
    popts.fig = spm_figure('Create','Graphics','visible','on'); 
  end
  colormap gray; 

  % main text report box with header (filename)
  popts.ax = axes('Position',[0.01 0.75 0.98 0.245],'Visible','off','Parent',popts.fig);
  [pp,ff,ee] = spm_fileparts( spm_str_manip( fname ,sprintf('k%d',80 - numel(resdir))) ); 
  pp = strrep( strrep( pp ,'\','\\'), '_','\_'); ff = strrep( strrep( ff ,'\','\\'), '_','\_');
  resdir2 = strrep( strrep( spm_str_manip( [filesep resdir filesep]),'\','\\'), '_','\_');
  fname2 = [pp '\color[rgb]{0,0.5,1}' resdir2 '\color[rgb]{0,0,0}' ff ee]; 
  text(0,0.99, ['Bone marrow extraction: ' fname2 '       '],...
    'FontSize',popts.fontsize+1,'FontWeight','Bold','Interpreter','tex','Parent',popts.ax);

end
function table1offset = printCLStable(job,out,popts,notab1)
% printCLStable. Values of each SPM cls variable (only useful for developers)

  if job.opts.expert > 1 && ~notab1
    mgid = find(out.spm8.mg > 2/numel(out.spm8.lkp));
    FN  = {'lkp','mg','mn','vr'};
    FNn = {sprintf('%d/%d SPM-Classes ',numel(mgid),numel(out.spm8.lkp)),'Propotion','Mean','Std'}; 
    FNt = {'%d','%0.2f','%0.2f','%0.2f'};
    % gray table rows
    for fni = 1:numel(FN)
      if fni==1
        text(0.01,0.95-            0.035,repmat(' ',1,1205), ...
          'BackgroundColor',[0.94 0.94 .94],'Parent',popts.ax,'FontSize',popts.fontsize*.6); 
      elseif mod(fni,2)
        text(0.01,0.95-(0.035*fni)-0.035,repmat(' ',1,1205), ...
          'BackgroundColor',[0.96 0.96 .96],'Parent',popts.ax,'FontSize',popts.fontsize*.6); 
      end
    end
    % add table content
    for fni = 1:numel(FN)
      seg8text(fni,1)   = text(0.01, 0.96-(0.045*fni) +0.01*(fni==1), ...
        sprintf('\\bf%s',FNn{fni}) ,'Fontname',popts.fontname,'FontSize',...
        popts.fontsize,'color',popts.fontcolor,'Interpreter','tex','Parent',popts.ax); 
      if fni==1, seg8text(fni,1).FontWeight = 'bold'; norm = 1; else, norm = out.spm8.(FN{fni})(2); end
      for lkpi = 1:numel(mgid) 
        seg8text(fni,mgid(lkpi)+1) = text(0.11 + 0.87/numel(mgid) * (lkpi), 0.95-(0.045*fni) + ...
          0.01*(fni==1), sprintf( FNt{fni}, out.spm8.(FN{fni})(mgid(lkpi))/norm ),...
          'Fontname',popts.fontname,'FontSize',popts.fontsize,'color',popts.fontcolor, ...
          'Interpreter','tex','Parent',popts.ax,'HorizontalAlignment','right');
        if fni==1 
          seg8text(fni,mgid(lkpi)+1).FontWeight = 'bold'; 
          seg8text(fni,mgid(lkpi)+1).Color      = popts.clscol(out.spm8.lkp(mgid(lkpi)),:);  
          seg8text(fni,mgid(lkpi)+1).String     = [seg8text(fni,mgid(lkpi)+1).String ...
            sprintf('(%d)',sum(out.spm8.lkp == out.spm8.lkp(mgid(lkpi))) )];
        end
        if fni==2, seg8text(fni,mgid(lkpi)+1).Color = popts.grycol( round(8 - 7*out.spm8.mg(mgid(lkpi))) ,:);  end
        if fni==3 
          try
            col = popts.tiscol( max(0, min(200,round(100*out.spm8.mn(mgid(lkpi))/max(out.spm8.mn(2))))) ,:); 
            seg8text(fni,mgid(lkpi)+1).Color = max(0,min(1,col - max(.3,min(sum(col,2)-2)/4)));  
          end
        end
        if fni==4
          seg8text(fni,mgid(lkpi)+1).Color = ...
            popts.grycol( min(10,max(1,round(8 - 7 * out.spm8.vr(mgid(lkpi)) / max(out.spm8.vr(1:3))))) ,:);  
        end
      end
    end
    table1offset = 0.27;
  else 
    table1offset = 0.02; 
  end
end
function printTPMtable(job,out,popts)
%printTPMtable. Print averaged values for each TPM class.

  out.tis.clsn = {'1-GM','2-WM','3-CSF','4-bone','5-head','6-BG'};
  for ci = 1:min(6,max(out.spm8.lkp))
    % indicate cases where a Gaussian is used only for a small volume, that
    % is pointing to a too large number of Gaussians
    % ... would have to consider also the difference between Gaussians
    if sum(out.spm8.lkp==ci) > max(1,sum(out.spm8.mg(out.spm8.lkp==ci)'<.2 | ...
        abs(gradient(out.spm8.mn(out.spm8.lkp==ci) ./ out.tis.WMth)) < .2)) 
      mark = '\color[rgb]{0,0.5,1}';
    else
      mark = ''; 
    end
    out.tis.clsG{ci} = sprintf('%d %s(%d)', sum(out.spm8.lkp==ci), mark, ...
      max(1,sum(out.spm8.mg(out.spm8.lkp==ci)'<.2 | abs(gradient(out.spm8.mn(out.spm8.lkp==ci) ./ out.tis.WMth)) < .2)) );  
  end
  % #####################
  % * add QC value for tissue ? 
  %   - you would need some colorcoding for it and right now the second 
  %     volume plot gives already some overview 
  % #####################
  if out.spm8.isCTseg && ~job.opts.normCT 
    out.tis.seg8x = out.tismri.Tth;  % ... see comments for further measures but better to keep it simple       
    FN  = {'clsn','seg8x','vol','volr'};                                  % 'den','clsG','seg8nv',
    FNn = {'TPM class','Med.Int.','Volume mm3','Volume (%TIV)'};          % 'Volume density','n-Class(well)','Main-std',
    FNt = {'%s','%0.0f','%0.0f','%0.1f'};                                 % '%0.0f','%s','%0.3f',
  else
    FN  = {'clsn','seg8n','vol','volr'};                                  % 'den','clsG','seg8nv',
    FNn = {'TPM class','Norm.Med.Int.','Volume mm3','Volume (%TIV)'};     % 'Volume density','n-Class(well)','Main-std',
    FNt = {'%s','%0.3f','%0.0f','%0.1f'};                                 % '%0.0f','%s','%0.3f',
  end
  if isempty(out.tismri) || ~isfield(out.tismri,'vol')
    FN = FN(1:3); FNn = FNn(1:3); FNt = FNt(1:3);
  else
    out.tis.vol     = out.tismri.vol; 
    out.tis.volr    = out.tismri.vol ./ sum(out.tismri.vol(1:3)) * 100; 
    % do not show BG volumes 
    out.tis.vol(6)  = nan; 
    out.tis.volr(6) = nan;
  end
  for fni = 1:numel(FN)
    if fni==1
      text(0.01,0.95-popts.table1offset-            0.035,repmat(' ',1,1000 + 205),...
        'BackgroundColor',[0.94 0.94 .94],'Parent',popts.ax,'FontSize',popts.fontsize*.6); 
    elseif mod(fni,2)
      text(0.01,0.95-popts.table1offset-(0.045*fni)-0.0,repmat(' ',1,1000 + round(205*(.6/.4))),...
        'BackgroundColor',[0.96 0.96 .96],'Parent',popts.ax,'FontSize',popts.fontsize*.4); 
    end
  end
  for fni = 1:numel(FN)
    tistext(fni,1)   = text(0.01              , 0.95-popts.table1offset-(0.045*fni) +0.01*(fni==1), ...
      sprintf('\\bf%s',FNn{fni}) ,'Fontname',popts.fontname,'FontSize',...
      popts.fontsize,'color',popts.fontcolor,'Interpreter','tex','Parent',popts.ax);
    if fni==1, set(tistext(fni,1),'FontWeight','bold'); end
    for lkpi=numel(out.tis.clsn):-1:1 % down to avoid overlap
      if iscell(out.tis.(FN{fni})), tmpval = out.tis.(FN{fni}){lkpi}; else, tmpval = out.tis.(FN{fni})(lkpi); end
      if ~isnan(tmpval)
        tistext(fni,lkpi+1) = text(0.12 + 0.082*(lkpi), ...
          0.95-popts.table1offset-(0.045*fni) +0.01*(fni==1), sprintf( FNt{fni}, tmpval ),...
          'Fontname',popts.fontname,'FontSize',popts.fontsize,'color',popts.fontcolor, ...
          'Interpreter','tex','Parent',popts.ax,'HorizontalAlignment','right');
        if fni==1, set(tistext(fni,lkpi+1),'FontWeight','bold','Color',popts.clscol(lkpi,:));  end
      end
      %if fni==3, tistext(fni,lkpi+1).Color = popts.grycol( max(1,min(10,round(8 - tmpval))) ,:);  end
    end
    % ###################
    % * add colorrating for tissues 
    % * density vs. volume => error?
    % ###########
  end
end
function printMainTable(job,out,popts,Si,St)
%printMainTable. Plot table with major bone and head measures. 

  % short variables and names
  highBGs     = {'low','high','high(MT)','mid'};  
  out.tis.highBGn = highBGs{out.tis.highBG + 1}; 
  out.tis.res_RES = mean(out.tis.res_vx_vol.^2 )^.5;
  if ~isempty(St)
    out.tis.iBone   = cat_stat_nanmean(Si.facevertexcdata(Si.facevertexcdata<100));
    out.tis.thBone  = cat_stat_nanmean(St.facevertexcdata(St.facevertexcdata<100 & St.facevertexcdata>0));
  end
  try out.tis.headthickmn = out.tismri.headthickmn; end 
% ################## 
%  * add skull-stripping and defacing value later 
%  * better definition of low fat intensity rating 
%  * do QC with CNR rating? 
% ################################# 
  FN{1}  = {'weightingn', 'res_RES' , 'highBGn' , 'headFatTypen', 'headBoneTypen', 'bonecortex', 'bonemarrow' , 'bonedensity' }; 
  FNf{1} = {'tis'       , 'tis'     , 'tis'     , 'tis'         , 'tis'          , 'tis'       , 'tis'        , 'tis'         };
  FNn{1} = {'Tw'        , 'Tres'    , 'Tbg'     , 'Tfat'        , 'Tbone'        , 'Tbcor'     , 'Tbmar'      , 'Tbdns'       };
  FNt{1} = {'%s'        , '%0.2f'   , '%s'      , '%s'          , '%s'           , '%0.3f'     , '%0.3f'      , '%0.3f'       }; 
  
  
  mgid = find(out.spm8.mg > 2/numel(out.spm8.lkp));
  lkpi = numel(mgid);
  % - the SPM measures minbone and maxBone are not good but medbone is ok > (1) medBoneS 
  % - the classic > (2) medBoneC ( correction option only internal test ) and (3) fst_vol4 (masked?)  
  % - the new     > (3) OccBone
  if 0 %out.spm8.isCTseg 
    out.tis.minBone = out.tismri.iBonemn3(1);
    out.tis.medBone = out.tismri.iBonemn3(2);
    out.tis.maxBone = out.tismri.iBonemn3(3);
    FN{2}  = {'minBone'   , 'medBone' , 'maxBone' , 'thBone'  , 'headthickness'};
    FNf{2} = {'tismri'    , 'tismri'  , 'tismri'  , 'tismri'  , 'tismri'};
    FNn{2} = {'lBone*'    , 'mBone*'  , 'hBone*'  , 'BthBone' , 'thhead'};
    FNt{2} = {'%0.0f'     , '%0.0f'   , '%0.0f'   , '%0.1f mm', '%0.1f mm'};
    FNi{2} = {1           , 1         , 1         , 1         , 1         };
  else
    %%
    if job.opts.bmethod==1, xr = 'v'; else, xr = 's'; end
    if job.opts.bmethod==1, bmf = 'classic'; else; bmf = 'tismri'; end
    FN{2}  = {'bone_med' , 'bonecortex' , 'bonemarrow' , 'bonethickness' , 'headthickness' , 'head'      };
    FNf{2} = {bmf        , [xr 'ROI']   , [xr 'ROI']   , [xr 'ROI']      , [xr 'ROI']      , [xr 'ROI']  };
    FNn{2} = {'boneMed'  , [xr 'Bcor']  , [xr 'Bmar']  , [xr 'Bth']      , [xr 'Hth']      , [xr 'Hmed'] };
    FNt{2} = {'%0.3f'    , '%0.3f'      , '%0.3f'      , '%0.3f'         , '%0.3f'         , '%0.3f'     };
    FNi{2} = {1          , 3            , 4            , 1               , 1               , 1           };
  end
  %%
  for fnj = 1:2
    text(0.01,0.64 - (fnj-1)*0.12 - popts.table1offset,repmat(' ',1,1000 + 120),... 
      'BackgroundColor',[0.94 0.94 .94],'Parent',popts.ax,'FontSize',popts.fontsize*.6); 
    for fni = numel(FN{fnj}):-1:1
      if isfield(out, FNf{fnj}{fni}) && isfield(out.(FNf{fnj}{fni}), FN{fnj}{fni})
        segtext(fni,lkpi)   = text(0.04 + 0.075*(fni-1), 0.64-(fnj-1)*0.12-popts.table1offset-(0.00) , ...
          sprintf('\\bf%s',FNn{fnj}{fni}) ,'Fontname',popts.fontname,'FontSize',...
          popts.fontsize,'color',popts.fontcolor,'Interpreter','tex','Parent',popts.ax, ...
          'FontWeight','bold','HorizontalAlignment','center');
        if ~isempty( FNt{fnj}{fni} )
          if isnumeric( out.(FNf{fnj}{fni}).(FN{fnj}{fni}) ) && numel(  out.(FNf{fnj}{fni}).(FN{fnj}{fni}) )>1
            val =  out.(FNf{fnj}{fni}).(FN{fnj}{fni})(FNi{fnj}{fni}); 
          else
            val =  out.(FNf{fnj}{fni}).(FN{fnj}{fni}); 
          end
          segtext(fni,lkpi+1) = text(0.04 + 0.075*(fni-1), ...
            0.64-(fnj-1)*0.12-popts.table1offset-(0.05), sprintf( FNt{fnj}{fni}, val ),...
            'Fontname',popts.fontname,'FontSize',popts.fontsize,'color',popts.fontcolor, ...
            'Interpreter','tex','HorizontalAlignment','center','Parent',popts.ax);
        else
          val =  out.(FNf{fnj}{fni}).(FN{fnj}{fni});
        end
        if fni==2
          try
            set( segtext(fni,lkpi+1), 'Color', ...
              popts.mrkcol( max(1,min(100,round( 100 * (val+1) / 4 ))) , :)); 
          end
        end % RES
        if fnj == 1
          if fni >2 && fni < 6
            switch get(segtext(fni,lkpi+1),'String')
              case 'low',  set(segtext(fni,lkpi+1),'Color',[ .0 .0 1.0]); 
              case 'mid',  set(segtext(fni,lkpi+1),'Color',[ .6 .3  .0]); 
              case 'high', set(segtext(fni,lkpi+1),'Color',[ .8 .2  .0]); 
              otherwise,   set(segtext(fni,lkpi+1),'Color',[1.0 .0  .0]); 
            end
          end
        else
          if 0 %out.spm8.isCTseg
            if fni > 0 && fni<4
              set(segtext(fni,lkpi+1),'Color',popts.jetcol( max(1,min(100,round(0.07 * out.(FNf{fnj}{fni}).(FN{fnj}{fni})))),:));  
            end
            if fni > 3        
              set(segtext(fni,lkpi+1),'Color',popts.jetcol( max(1,min(100,round(5 * out.(FNf{fnj}{fni}).(FN{fnj}{fni})))),:));  
            end
          else
            if fni > 4
              try
                set(segtext(fni,lkpi+1),'Color',popts.jetcol( max(1,min(100,round(5 * out.(FNf{fnj}{fni}).(FN{fnj}{fni})))),:)); 
              end
            end
          end
        end
      end
    end
  end

end
function printHistogram(Ym,Yc,Ybonemarrow,job,out,popts)
%printHistogram.

  if job.output.report > 1  && ~strcmpi(spm_check_version,'octave') % && job.opts.bmethod > 0
    %%
    % plot a white box over the too long table 
    axes('Position',[0.61 0.855 - popts.table1offset/3 - 0.05 0.40 0.18],'Parent',popts.fig,'Color',[1 1 1]); 
    ax2 = gca; set(ax2,'YColor',[1 1 1],'XColor',[1 1 1],'YTick',[],'XTick',[]); 
    
    % real plot for histogram
    axes('Position',[0.62 0.855-popts.table1offset/3 0.37 0.13],'Parent',popts.fig,'Color',[1 1 1]); 
    ax2 = gca; hold on;
    
    % print histogram
    %plot([1 1],[0 max(ylim)],'--','Color',[.7 .7 .7],'LineWidth',1.5);
    clsi = numel(Yc):-1:1; clsi( flip(cellfun('isempty',Yc)) ) = []; 
    if out.spm8.isCTseg && ~job.opts.normCT 
      % CT case
      for ci = clsi
        if ~strcmpi(spm_check_version,'octave') 
          hhst(ci) = histogram(ax2,Ym( Yc{ci}(:)>.5 & Ym(:)>-1000 & Ym(:)<2000) ,-1000:10:2000, ...
            'LineStyle','none','FaceColor',popts.clscol(ci,:));
        end
        hstmax(ci) = max( hhst(ci).Values )  / 3;
      end
      xlim([-300 1800]); %title('normalized intensities');
    else
      % default MRI case with all volumes 
      for ci = clsi % inverse order to have relevant things in the foreground
        if ~strcmpi(spm_check_version,'octave') 
          hhst(ci) = histogram(ax2,Ym( Yc{ci}(:)>.5 & Ym(:)>0 & Ym(:)<2) , 0:0.01:2, ...
            'LineStyle','none','FaceColor',popts.clscol(ci,:));
          hstmax(ci) = max( hhst(ci).Values ) ;
        else
          %hhst(ci) = 
          hist(ax2,Ym( Yc{ci}(:)>.5 & Ym(:)>-1000 & Ym(:)<2000) ,-1000:10:2000); 
          %, ...
          %  'LineStyle','none','FaceColor',popts.clscol(ci,:));
        end
      end
      xlim([0 2]); set(ax2,'XTick',0:1/6:2); % ax2.XTick = 0:1/6:2; % grid on;  %title('normalized intensities');
    end

    %% setup axes of the histogram box
    if numel(clsi)>1, clsi(1) = []; end
    ylim([0 max([1,hstmax(clsi)*1.3] )]); box on; 
    if ~(out.spm8.isCTseg && ~job.opts.normCT)
      set(ax2,'FontSize',popts.fontsize * 0.85, 'XTickLabelRotation', 0, ... 
        'XTickLabel', {'0','','','0.5','','','1','','','1.5','','','2'});
      xlabel('normalized intensities with normalized SPM classes'); 
    else
      xlabel('CT intensities of SPM classes'); 
    end
    ax2.YTick = []; 


    % add tissue peaks to the histogram
    lstyle = {':','-'};
    if 0 % out.spm8.isCTseg && ~job.opts.normCT 
% #####################      
      for ci = 1:3
        pl = plot([out.tismri.Tth(ci) out.tismri.Tth(ci)],[0 max(ylim)]); 
        set(pl,'Color',[ popts.clscol(min(6,max(out.spm8.lkp(ci))),:)]);
      end
      for ci = 1:numel(out.tismri.iBonemn)
        pl = plot( repmat(out.tismri.iBonemn(ci),1,2) ,[0 max(ylim)]); 
        set(pl,'Color',[ popts.clscol(4,:)]);
      end
      pl = plot( repmat( out.tis.report.fat ,1,2) ,[0 max(ylim)]); 
      set(pl,'Color',[ popts.clscol(5,:)]);
      pl = plot( repmat( out.tis.report.muscle ,1,2) ,[0 max(ylim)]); 
      set(pl,'Color',[ popts.clscol(5,:)]);
    else
      for ci = 1:numel(out.spm8.mn)
        pl = plot([out.spm8.mn(ci) out.spm8.mn(ci)] / out.tis.WMth,[0 max(ylim)]); 
        set(pl,'LineWidth',out.spm8.mg(ci) * 2,'LineStyle', ... 
          lstyle{ 1 + ( out.spm8.mg(ci) == max(out.spm8.mg( out.spm8.lkp==out.spm8.lkp(ci) ) )) } );
        set(pl,'Color',[ popts.clscol( min(6,max(out.spm8.lkp(ci))),:) ...
          min(1,max(0,out.spm8.mg(ci) * .7*sum(out.spm8.lkp(ci)==out.spm8.lkp).^.5))]);
      end
    end

    if numel(Yc)==6 
      lg = legend(flip({'GM','WM','CSF','bone','head','BG'}),'box','off'); 
    else 
      lg = legend({'bone','GM','WM','CSF'},'box','off'); % ######### this need refinement but it is ok
    end 
    lgpos = get(lg,'Position'); lgpos(1) = .89; set(lg,'Position',lgpos); 
  end  
end
function printErrors
%printErrors. Print errors ... in development

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
function printVolumes(Po,Vo,Ym,Yc,Ybonemarrow,St,job,out,popts)
%printVolumes. Call SPM orthview and print the MRI image and its segments.

  if job.output.report > 1 % report 1 is only the table 

    clear -global st % this is the spm_orthview variable for the slice images
    global st 

    spm_orthviews('Reset')
    pos = {[0.008 0.375 0.486 0.35]; [0.506 0.375 0.486 0.35]};
    % T1 + SPM segmentation
    if out.tis.weighting < 0 
      colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(popts.crange+2)),1); ...
        repmat([ 0 .3 .6],round(59/(popts.crange+2)),1); jet(59 - 2*round(59/(popts.crange+2)))]);
      V0 = Vo; V0.dat(:,:,:) = single(.3*Yc{3} + .75*Yc{1} + 1*Yc{2} + 1.2*Yc{4} + 0.5*Yc{5} ); V0.dt(1) = 16;
      V0.pinfo  = repmat([1;0],1,size(Ym,3));
    elseif job.opts.bmethod == 0
      colormap gray; 
      V0       = out.spm8.image;
      Vc4      = spm_vol(spm_file(out.spm8.image.fname,'prefix','c4'));
      Vc4.mat  = out.spm8.Affine  * Vc4.mat; 
    else
      colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(popts.crange+2)),1); ...
        repmat([ 0 .3 .6],round(59/(popts.crange+2)),1); jet(59 - 2*round(59/(popts.crange+2)))]);
      V0 = Vo; V0.dat(:,:,:) = single(Ym); V0.dt(1) = 16;
      V0.pinfo  = repmat([1;0],1,size(Ym,3));
    end
    V0.mat    = out.spm8.Affine * out.spm8.image.mat; 
    % load first to avoid problems with high resolutions
    hh0       = spm_orthviews('Image',spm_vol(fullfile(spm('dir'),'tpm','TPM.nii,1')),pos{1}); 
    spm_orthviews('window',hh0,[0 10000]); % just make it black
    spm_orthviews('BB', [-85 -120 -90; 85 95 105]); % this has to be set in the low-resolution image
    hh0       = spm_orthviews('Image',V0,pos{1}); % add correct image after the other settings! 
    spm_orthviews('Reposition',[-25 0 0]);
    if out.tis.weighting < 0 
      spm_orthviews('Caption',hh0,sprintf('%s','CTseg'));
    else
      spm_orthviews('Caption',hh0,sprintf('%s',out.tis.weightingn));
    end
    if job.opts.bmethod == 0
      spm_orthviews('addtruecolourimage',hh0,Vc4,[0 0 0; 0.5 0.5 0.5; 1 0 0],0.4,1,0.1); 
      spm_orthviews('redraw');

      orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'), ...
        'Type','Image','Tag',''),'parent');
      set(orthviewlegend,'YTick',[0 .5 1],'YTickLabel',{'BG','~Bone','Bone'}, ...
        'FontSize',get(orthviewlegend,'FontSize') * 0.85,'YAxisLocation','right', ...
        'Position',get(orthviewlegend,'Position') - [0.04 0 0.04 0]); 
      
      try % does not work in headless mode without java
        figfs10 = [ 
          findobj(popts.fig,'FontSize',11); 
          findobj(popts.fig,'FontSize',10); 
          findobj(popts.fig,'FontSize',9); 
          findobj(popts.fig,'FontSize',8); 
          findobj(popts.fig,'FontSize',popts.fontsize*.6);  
          findobj(popts.fig,'FontSize',popts.fontsize*.4); ]; 

        % print figure with reduced fontsizes
        for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .8; end; end 
        saveas(popts.fig,out.P.report,'png')
        for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .8; end; end  
      catch
        cprintf('err','Cannot print report.')
      end
      
      return
    else
      spm_orthviews('window',hh0,[0 1.2]); 
    end
    
   
    
    % print cls overlay?
    if out.tis.weighting >= 0
      Vp0        = Vo; 
      Vp0.dt(1)  = 16;
      Vp0.mat    = out.spm8.Affine * Vo.mat; 
      Vp0.pinfo  = repmat([1;0],1,size(Ym,3));
      switch 2
        case 1
          Vp0.dat(:,:,:) = single( (Yc{4}>.5) + (2*(Yc{5}>.5) + 0.5*(Yc{5}>.5)) ); 
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 1 0 0; 1 1 0],0.4,2,0.1); 
        case 2
          Vp0.dat(:,:,:) = single( (Yc{4}>.5) + (Yc{4}>.5 & Ym>.7 ) + (3*(Yc{5}>.5)) + (Yc{5}>.5 & (Ym>.9)) ); 
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 0 0.5 1; .5 0.8 1; .9 .5 0.5; 1 1 0 ],0.5,4,0); 
        case 5
          Vp0.dat(:,:,:) = single(2*(Yc{1}>.5) + 3*(Yc{2}>.5) + (Yc{3}>.5) + 4*(Yc{4}>.5) + 5*(Yc{5}>.5));
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 0 0 1; 0 .8 0; 1 1 1; 1 0 0; 1 1 0],0.4,5,0); 
      end
   end
  
    % print bone marrow
    if out.tis.weighting < 0
      V1 = Vo; V1.dat(:,:,:) = single((1000 + Ym)/2300); V1.dt(1) = 16;
    else
      Ybonemarrow2 = Ybonemarrow; Ybonemarrow2(isnan(Ybonemarrow2)) = 0; 
      V1 = Vo; V1.dat(:,:,:) = min(popts.crange,Ybonemarrow2) + 2*(Ybonemarrow2>0) + ...
        single(Yc{1}>.5) + single(Yc{5}>.5) + 2*single(Yc{2}>0.5);
    end
    V1.pinfo  = repmat([1;0],1,size(Ym,3));
    V1.mat    = out.spm8.Affine * V1.mat; 
    V1.dt(1)  = 16;
    hh1       = spm_orthviews('Image',V1,pos{2}); 
    spm_orthviews('Interp',1);
    if out.tis.weighting < 0
      spm_orthviews('window',hh1,[0 1.2]);
      spm_orthviews('Caption',hh1,sprintf('%s',out.tis.weightingn));
    else
      spm_orthviews('window',hh1,[0 popts.crange+2]);
      spm_orthviews('Caption',hh1,'bonemarrow + GM/HD + 2*WM');
    end
    spm_orthviews('Reposition',[-25 0 0]); pause(0.01)
    
    % this is replaced by spm_orthviews('redraw');
    spm_orthviews('redraw');
    orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'), ...
      'Type','Image','Tag',''),'parent');
    if exist('orthviewlegend','var') && numel(orthviewlegend)>0 % no CTseg so far
      if numel(orthviewlegend)>1, orthviewlegend = orthviewlegend{1}; end
      if job.opts.bmethod > 0
        set(orthviewlegend(1),'YTick',[0 1 2 3 4], ...
          'YTickLabel', {' background',' bonecortex',' bonemarrow',' muscle',' fat'}, ...
          'FontSize',get(orthviewlegend,'FontSize') * 0.85,'YAxisLocation','right', ...
          'Position',get(orthviewlegend,'Position') - [0.04 0 0.04 0]); 
      else
        set(orthviewlegend(1),'YTick',[0 1],'YTickLabel',{' background',' bone'}, ...
          'FontSize',get(orthviewlegend,'FontSize') * 0.85,'YAxisLocation','right', ...
          'Position',get(orthviewlegend,'Position') - [0.04 0 0.04 0]); 
      end
    end

    
    %% bone surface 
    if isfield(Po,'central') && ~isempty(St)
      for idi = 2:3 
        spm_orthviews('AddContext',idi);  
        spm_ov_mesh('display',idi,out.P.central );
        % apply affine transformation
        V = (out.spm8.Affine * ([st.vols{idi}.mesh.meshes(end).vertices,...
             ones(size(st.vols{idi}.mesh.meshes(end).vertices,1),1)])' )';
        V(:,4) = [];
        st.vols{idi}.mesh.meshes = subsasgn(st.vols{idi}.mesh.meshes, ...
          struct('subs','vertices','type','.'),single(V));
        % change line style
        hM = findobj(st.vols{idi}.ax{1}.cm,'Label','Mesh');
        UD = get(hM,'UserData');  
        UD.width = 1;
        if out.tis.weighting == -1
          UD.style = {'w-'};  
        else
          if idi==2
            UD.style = {'r-'}; 
          elseif out.tis.headBoneType>1
            UD.style = {'w-'}; 
          else
            UD.style = {'r-'};  
          end
        end
        set(hM,'UserData',UD); clear UD hM
        warning('off','MATLAB:subscripting:noSubscriptsSpecified');
        spm_ov_mesh('redraw',idi);
        if exist('orthviewlegend','var')
          orthviewlegend.Position(3) = 0.01;  
        end
      end
    end
  end
end
function printSurfaces(St,Si,job,out,popts)
%printSurfaces.

  if job.opts.bmethod > 1
    % print thickness surface map 
    if ~isempty(St) && job.output.report > 1 
      % thickness map
      hCS{1} = subplot('Position',[0.015 0.20 0.23 0.15],'Parent',popts.fig,'visible','off');  sview{1} = 'l';
      hCS{2} = subplot('Position',[0.255 0.20 0.23 0.15],'Parent',popts.fig,'visible','off');  sview{2} = 'r';
      hCS{3} = subplot('Position',[0.015 0.01 0.23 0.19],'Parent',popts.fig,'visible','off');  sview{3} = 't';
      hCS{4} = subplot('Position',[0.255 0.06 0.23 0.14],'Parent',popts.fig,'visible','off');  sview{4} = 'p';
      hCS{5} = subplot('Position',[0.255 0.03 0.23 0.02],'Parent',popts.fig,'visible','off'); 
      
      imat = spm_imatrix(out.spm8.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
      V = (Rigid * ([St.vertices, ones(size(St.vertices,1),1)])' )'; V(:,4) = []; St.vertices = V;
  
      % in case of octave the rendering works different
      if strcmpi(spm_check_version,'octave') 
        Si.vertices = [Si.vertices(:,2) Si.vertices(:,1) Si.vertices(:,3)];
        St.vertices = [St.vertices(:,2) St.vertices(:,1) St.vertices(:,3)];
      end
  
      for ci = 1:4
        %if ~strcmpi(spm_check_version,'octave') 
          h{ci} = cat_surf_render2(St,'parent',hCS{ci}); 
          cat_surf_render2('Clim',h{ci},[0 20]); 
          switch lower(sview{ci})
            case {'r'},  cat_surf_render('view',h{ci},[  90   0]); 
            case {'l'},  cat_surf_render('view',h{ci},[ -90   0]);  
            case {'t'},  cat_surf_render('view',h{ci},[   0  90]); 
            case {'p'},  cat_surf_render('view',h{ci},[   0   0]); 
          end
        %else
        % in case of octave we want to use the volume render but it is not working yet 
        %  cat_surf_renderv(Stm,[],struct('view',sview{ci},'mat',spm_imatrix(eye(4)),'h',hCS{2},'interp',interp*0.9));
        %end
      end
      
      cb1 = cat_surf_render2('Colorbar',h{4});
      set( cb1.colourbar, 'Location','South',  'Position',[0.275 0.015 0.2 0.005] );
      
      text(0.5,0.015, 'Bone thickness (mm)', 'fontsize',popts.fontsize-1, ...
        'FontName',popts.fontname, 'HorizontalAlignment','center', 'Parent',hCS{5}); 
  
    end
  
  
  
    if ~isempty(Si) && job.output.report > 1
      hCS{1} = subplot('Position',[0.5 + 0.015 0.20 0.23 0.15],'Parent',popts.fig,'visible','off');  sview{1} = 'l';
      hCS{2} = subplot('Position',[0.5 + 0.255 0.20 0.23 0.15],'Parent',popts.fig,'visible','off');  sview{2} = 'r';
      hCS{3} = subplot('Position',[0.5 + 0.015 0.01 0.23 0.19],'Parent',popts.fig,'visible','off');  sview{3} = 't';
      hCS{4} = subplot('Position',[0.5 + 0.255 0.06 0.23 0.14],'Parent',popts.fig,'visible','off');  sview{4} = 'p';
      hCS{5} = subplot('Position',[0.5 + 0.255 0.03 0.23 0.02],'Parent',popts.fig,'visible','off'); 
  
      imat   = spm_imatrix(out.spm8.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
      V      = (Rigid * ([Si.vertices, ones(size(Si.vertices,1),1)])' )'; V(:,4) = []; Si.vertices = V;
  
      for ci = 1:4
        try
          h{ci} = cat_surf_render2(Si,'parent',hCS{ci}); 
        catch 
          cat_io_cprintf('err','cat_surf_render2 error\n');
        end
        if out.tis.weighting < 0 
          cat_surf_render2('Clim',h{ci},[500 1500]); 
        else
          cat_surf_render2('Clim',h{ci},[0 popts.crange]); 
        end
        switch lower(sview{ci})
          case {'r'},  cat_surf_render('view',h{ci},[  90   0]); 
          case {'l'},  cat_surf_render('view',h{ci},[ -90   0]);  
          case {'t'},  cat_surf_render('view',h{ci},[   0  90]); 
          case {'p'},  cat_surf_render('view',h{ci},[   0   0]); 
        end
      end
  
      if 1
        % no bone histogram
        cb2 = cat_surf_render2('Colorbar',h{4});
        if isfield(cb2,'colourbar')
          set(cb2.colourbar,'Location','South','Position',[0.5 + 0.275 0.015 0.2 0.005]);
          text(0.5,0.015,'Bonemarrow intensity','fontsize',popts.fontsize-1,'FontName',popts.fontname,...
            'HorizontalAlignment','center','Parent',hCS{5}); 
        end
      else
        %% colormap
        axes('Position',[0.965 0.03 0.01 0.28],'Parent',popts.fig); 
        image(flip(121:1:120+surfcolors)','Parent',cc{4});
        cc{4} = gca; % for octave!
        
        %% histogram line
        axes('Position',[0.936 0.03 0.03 0.28], 'Parent', popts.fig, ...
          'Visible', 'off','tag', 'cat_surf_results_hist', ...
          'xcolor', popts.fontcolor, 'ycolor', popts.fontcolor);
        cc{5} = gca; % for octave!
      
        side  = hSD{1}.cdata;
        [d,h] = hist( side(~isinf(side(:)) & ~isnan(side(:)) &  side(:)<6 & side(:)>0) ,  hrange);
        d = d./numel(side);
        d = d./max(d);
        text(cc{5},'test')
        
        % print histogram
        hold(cc{5},'on');  
        for bi = 1:numel(d)
          b(bi) = barh(cc{5}, h(bi), -d(bi), boxwidth); 
          set(b(bi), 'Facecolor', cmap3(bi,:), 'Edgecolor', popts.fontcolor); 
        end
        ylim([0,20]); xlim([-1 0]);
      end
    end
  
  
    % restore colormap?
    colormap( [ [0 0.02 0.07]; repmat([0.05 0.15 .35], round(59/(popts.crange+2)),1); ...
      repmat([ 0 .3 .6],  round(59/(popts.crange+2)), 1); jet(59 - 2*round(59/(popts.crange+2)))] );
    cb1.colourbar.Limits      = [5 20];
    cb1.colourbar.Ticks       = 5:15/4:20;
    cb1.colourbar.TickLabels  = 0:5:20;
    cb2.colourbar.Limits      = [2 8];
    cb2.colourbar.Ticks       = 2:6/2:8;
    cb2.colourbar.TickLabels  = {'bone (0)' '~WM (4) ' 'fat (8)'};

  
    % final print with try to avoid crashing for unknown reasons
    try 
      %% does not work in headless mode without java
      figfs10 = [ 
        findobj(popts.fig, 'FontSize', popts.fontsize + 1    ); 
        findobj(popts.fig, 'FontSize', popts.fontsize        );
        findobj(popts.fig, 'FontSize', popts.fontsize - 1    );  
        findobj(popts.fig, 'FontSize', popts.fontsize * 0.85 ); 
        findobj(popts.fig, 'FontSize', popts.fontsize * 0.6  ); 
        findobj(popts.fig, 'FontSize', popts.fontsize * 0.4  ); ]; 
  
      % print with reduced fontsize
      for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .75; end; end 
      saveas(popts.fig,out.P.report,'png')
      for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .75; end; end 
    end
  end
end



