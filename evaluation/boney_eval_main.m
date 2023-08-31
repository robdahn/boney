% Bone aging test script.
% * processing of different datasets such as IXI, CAMCAN and NKI datasets 
%   of healthy (adults) subjects with full age range
% * Evaluation of 
%   * BMD measures in the UKB (non public by UKB)
%   * IXI / NKI / CamCAN >> aging effects
%   * OASIS >> CT vs. MRI , longitudinal, test-retest
%   * Buchert, Simon >> large-scale test-retest 
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


%
% ToDo: 
%  * restructure these batches and add parts of the data and description 
%    how to use (eg. download IXI and set path) 
%  * include and organize different test conditions and subdirs with
%    detailed reports 
%  * how to include the version or date?
%  * extend help 
%  * include setting for our data paths

clear groups

% data - groups (bias corrected m*.nii)
% * fasttest includes a smaller two group test subset of young and old IXI
%   subjects.
% * matonly==1 uses only the SPM preprocessing parameter of the *_seg8.mat 
%   (named as SPM) whereas matonly==0 call the so called MED pipeline that 
%   loads and refines the SPM tissue segments to estimate own volume and
%   intensity parameters.
% 0-MED, 1-SPMo;   2-SPMn, 3-MED2, 4-SPMo0, 5-MED0;    10-final0, 11-final1, 12-final2 
matonlyn  = {'MED','SPM',   'NEWSPM','NEWMED','SPM0','MED0','','','','','BoneySPM','BoneyV','BoneyS'}; 
if 1
  maindirs  = {
    '/Volumes/SG5TB/MRData/202301_MRIbones/'; 
    '/Volumes/SG5TB/MRData/202301_MRIbones/'; 
    '/Volumes/SG5TB/MRData/202301_MRIbones/'; 
    };
else
  maindirs  = {
    '/Users/polona/Downloads/202301_MRIbones';
    '/Volumes/OneTouch/MRData/202301_MRIbones';
    '/Users/polona/Downloads/IXI'; 
    }; 
end
resultdir = fullfile(maindirs{1},'results_Boney_20230828'); 
if ~exist(resultdir,'dir'), mkdir(resultdir); end
dname = 'aging';

% get the image data
%testset = % IXI, UKB
groups{1} = [ 
    ...cat_vol_findfiles(maindirs{3},'mIXI*.nii');
    ...NA - cat_vol_findfiles(maindirs{2},'camcan','BIDSsep','anat','msub*T1w.nii');
    ...NA - cat_vol_findfiles(maindirs{2},'NKI_RAW','NKI_T1','msub*T1w.nii');
    ...NA - cat_vol_findfiles(maindirs{2},'ukb_bones','msub*T1w.nii');
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_testset')      ,'mG*.nii'  ,struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_testsetExt')   ,'mG*.nii'  ,struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_subsample100a'),'msub*.nii',struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_subsample100b'),'msub*.nii',struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_subsample100c'),'msub*.nii',struct('maxdepth',1)));
    %cellstr(cat_vol_findfiles(fullfile(maindirs{1},'OASIS3_CT')    ,'c04_*CTseg.nii',struct('maxdepth',4)));
    ... setdiff( ...
    ...    cat_vol_findfiles(maindirs{2},'OASIS3','BIDS','msub*T1w.nii'), ...
    ...    cat_vol_findfiles(maindirs{2},'OASIS3','BIDS','msub*echo*T1w.nii'));
    ...setdiff( ...
      ...cat_vol_findfiles(maindirs{2},'NKI_RS','T1w','m*T1w.nii'), ...
      ...[cat_vol_findfiles(maindirs{2},'NKI_RS','T1w','mw*T1w.nii');
       ...cat_vol_findfiles(maindirs{2},'NKI_RS','T1w','ME*T1w.nii')]);
    ]; 

% get the phenotypic data
Pcsv{1} = {
    ...fullfile(maindirs{3},'IXI.csv');
    ...fullfile(maindirs{2},'camcan','BIDSsep','cc700-scored','participant_data.csv'); 
    ...fullfile(maindirs{2},'NKI_RAW','participants.tsv'); 
    ...fullfile(maindirs{2},'NKI_RAW',NKI_RS/201506/NKI.1-39.phenotypic.csv';
    %fullfile(maindirs{1},'OASIS3_CT','+tables/OASIS3_data_files/scans/demo-demographics/resources/csv/files/OASIS3_demographics.csv');
    fullfile(maindirs{1},'ukb_bones_testset','ukb_3000.csv');
    };

%
boney_eval_readCSV




% single tests
if 0
%% Processing
  matonly = 12;   
  out1 = boney_getBoneMeasures({groups{1}(1)},struct('matonly',matonly,'rerun',0));
  %out = getBoneMeasures(groups,struct('matonly',matonly,'rerun',0));
end

%%
%for matonly = [12, 11, 10, 0:5] % 0-MED, 1-SPMo;   2-SPMn, 3-MED2, 4-SPMo0, 5-MED0;    10-final0, 11-final1, 12-final2 
matonly = 12;   
out = boney_getBoneMeasures(groups,struct('matonly',matonly,'rerun',0));
% end
%%
plevel = 2; % 1-basic, 2-3-full
boney_eval_prepareMeasures

% scatter plot for each test paramter with linear fits
% print overview    
printfg = 0; 
boney_eval_printDetails % (this also involve some other steps >> seperate!)
boney_eval_printMain

%% print for each (bone) measure of one bone processing pipeline a detailed
printfg = 1;
boney_eval_printDetails
   







    %% create figure 14
    %  med surface-measures for both sexes, male and femal
    fg=14; fgh = figure(fg); clf(fg); fgh.Color = [1 1 1]; fgh.Position(3:4) = [900 50 + 20*numel(measure) ]; 
    ha = annotation('textbox',[0.005 0.92 0.99 0.08],'string',sprintf( ...
      'measures vs. parameter (r~color,p~transparency) - site %d, N=%d ', si, numel(age)));
    ha.FontSize             = 13;
    ha.FontWeight           = 'bold';
    ha.HorizontalAlignment  = 'center';
    ha.EdgeColor            = 'none'; 

    rGMV = out.rGMV{1}; rWMV = out.rWMV{1}; rCMV = out.rCMV{1}; TIV = out.TIV{1};
    clear RHOm PVALm
    % 3 group by sex
    group = {'both','males','females'}; 
    switch matonly
      case 2
        measure2 = [measure; ...
          {'bonemarrowmin_med' 'bonemarrowmin_std' 'bonemarrowmax_med' 'bonemarrowmax_std' ...
           'bonemarrow_med' 'bonemarrow_std' 'bonethickness_med' 'bonethickness_std' ...
           'head_med' 'head_std' 'headthickness_med' 'headthickness_std'}'];
      case {11,12}
        measure2 = measure; %
         [measure; ...
          {
          'bone_med'; 'tis_bone';  'tis_head';  'tismri_volfatr'; ... volume with high intensity (fatvolume)
          'vROI_bonecortex3'; 'vROI_bonethickness3'; 'vROI_headthickness3'; 
          'sROI_bonecortex3'; 'sROI_bonethickness3'; 'sROI_headthickness3'; }];
    end
    for gi = 1:3
      subplot('Position',[0.155  + 0.26*(gi-1)  0.20 + 8/(numel(measure)+1).^2  ...
                          0.255                 0.78 - 2/(numel(measure)+1)]); 
      switch gi
        case 1, sexm = ' '; 
        case 2, sexm = ' & sex(:)==1 '; 
        case 3, sexm = ' & sex(:)==2 ';
      end
      for mi = 1:numel(para)
        mjj = 0;
        for mj = 1:numel(measure2)
          p1 = eval( sprintf('%s(        ~isnan(%s(:)) %s & msk(:))', para{mi}   , para{mi}, sexm ));
          if size(p1,1)<size(p1,2), p1 = p1'; end
          if size( out.(measure2{mj}){1} , 2) == 1
            mjj = mjj + 1;
            p2 = eval( sprintf('out.%s{1}( ~isnan(%s(:)) %s & msk(:))', measure2{mj}, para{mi}, sexm ));
            if size(p2,1)<size(p2,2), p2 = p2'; end
            [RHOm(mi,mjj,gi),PVALm(mi,mjj,gi)] = corr( p1 , p2 ,'type', corrtype); 
            measure3{mjj} = measure2{mj};
          else
            % (1)=global, (2) 1=f, (3) 3=o, (4) 8=p, (5) 9=xx, (6) 10=x (7) 11=xx, (8) 12=p, (9) 13=xx
            for fi = [1 2 4 3 ] % 8 5 9 6 7] %1: size( out.(measure2{mj}){1} , 2) 
              mjj = mjj + 1;
              p2c{fi} = eval( sprintf('out.%s{1}( ~isnan(%s(:)) %s & msk(:), %d)', measure2{mj}, para{mi}, sexm , fi ));
              if size(p2c{fi},1)<size(p2c{fi},2), p2c{fi} = p2c{fi}'; end
              [RHOm(mi,mjj,gi),PVALm(mi,mjj,gi)] = corr( p1 , p2c{fi} ,'type', corrtype);
              measure3{mjj} = sprintf('%s(%d)',measure2{mj},fi);
            end
          end
        end
      end
      
      im = imagesc( (RHOm(:,:,gi)')); im.AlphaData = .001 ./ max(eps,PVALm(:,:,gi)'); 
      ax = gca; ax.TickLength = [0 0]; 
      ax.XTick = 1:numel(para);    ax.XTickLabel = cat_io_strrep(para   ,'_','\_');
      if gi==1
        ax.YTick = 1:numel(measure3); ax.YTickLabel = cat_io_strrep(measure3,'_','\_'); ylabel('measures'); 
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
      ax.FontSize = 9;
      title(group{gi})
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
      nxi = [3 4 10 12 ];
      plot( repmat(nxi+0.5, 2,1) , repmat([0 ny+1]', 1,numel(nxi)) , 'Color', zeros(3,1) ) % vertical lines
      %nxi = [3 4 10 12 ];
      %plot( repmat(nxi+0.5, 2,1) , repmat([0 ny+1]', 1,numel(nxi)) , 'LineWidth', 2, 'Color', zeros(3,1) ) % vertical lines
      switch matonly 
        case {11,12}
          % nyi = [4 7 13 19 25   10 16 22 28]; nyim = [7 13 19 25]; % full
          nyi = [4 ]; nyim = [6 9]; % short
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

      caxis([-1 1] * round( max(abs(RHO(RHO<1)))*10)/10 )
    
    end
    colormap(jet) %cat_io_colormaps('BWR',64)); %graybluered
    subplot('Position',[0.95,0.27,0.0001,0.6]); imagesc( (-1:.1:1)' ); axis equal off
      caxis([-1 1] * round( max(abs(RHO(RHO<1)))*10)/10 )
    colorbar
    




    %% create figure 12
    %  new measures sorted by correlations
    fg=12; fgh = figure(fg); clf(fg); fgh.Position(3:4) = [800 350]; fgh.Color = [1 1 1];
    for mi = 1:numel(measure)
      for mj = 1:numel(measure)
        msk = ~isnan(out.(measure{mi}){1}) & ~isnan(out.(measure{mj}){1});
        [RHOmm(mi,mj),PVALmm(mi,mj)] = corr(out.(measure{mi}){1}(msk),out.(measure{mj}){1}(msk),'type',corrtype); 
      end
    end
    title('correlation measures')
    for i=1:2
      if i==1
        subplot(1,2,i); 
        im = imagesc(RHOmm); im.AlphaData = .05 ./ PVALmm; 
        ax = gca; ax.TickLength = [0 0];
        ax.XTick = 1:numel(measure); ax.XTickLabel = cat_io_strrep(measure,'_','\_');
        ax.YTick = 1:numel(measure); ax.YTickLabel = cat_io_strrep(measure,'_','\_');
        title('correlation measures'); xlabel('measures'); ylabel('measures')
      else  
        subplot(1,2,i); 
        [~,RS] = sortrows(std(RHOmm)' + median(RHOmm)' + mean(RHOmm)','descend');
        im = imagesc(RHOmm(RS,RS)); im.AlphaData = .05 ./ PVALmm(RS,RS); 
        ax = gca; ax.TickLength = [0 0];
        ax.XTick = 1:numel(measure); ax.XTickLabel = cat_io_strrep(measure(RS),'_','\_');
        ax.YTick = 1:numel(measure); ax.YTickLabel = cat_io_strrep(measure(RS),'_','\_');
        title('sorted correlation measures'); xlabel('measures'); ylabel('measures')
      end
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
    
   


%%

figure(90), clf
msk = site>0; %msk = site==4;
f = scatter(out.rGMV{1}(sex==1 & msk), out.bone_med{1}(sex==1 & msk)','filled'); f.MarkerEdgeColor = [ 0 0 1]; f.MarkerFaceAlpha = 0.2; f.MarkerEdgeAlpha = 0.2;hold on
f = scatter(out.rGMV{1}(sex==2 & msk), out.bone_med{1}(sex==2 & msk)','filled'); f.MarkerEdgeColor = [ 1 0 0]; f.MarkerFaceAlpha = 0.2; f.MarkerEdgeAlpha = 0.2;
[f,r(1),p] = fit(out.rGMV{1}(sex==1 & msk), out.bone_med{1}(sex==1 & msk)','poly2','Normalize','on','Robust','Bisquare'); ph=plot(f); ph.Color = [0 0 1]; 
[f,r(2),p] = fit(out.rGMV{1}(sex==2 & msk), out.bone_med{1}(sex==2 & msk)','poly2','Normalize','on','Robust','Bisquare'); ph=plot(f);  
legend( ...
                    sprintf('Male (R^2=%0.3f)',r(1).adjrsquare), ...
                    sprintf('Female (R^2=%0.3f)',r(2).adjrsquare), ...
                    'Location','NorthWest');                    
xlabel('rGMV'); ylabel('median bone density'), ylim([0 2]);
saveas(gcf,'Camcan.png'); 
%%


%% create table
tab = cell(1,numel(Pcsv));
for csvi = 1:numel(Pcsv)
   tab{csvi}{1,1}       = 'Filename'; 
   tab{csvi}(2:numel(groups{csvi})+1,1) = spm_str_manip(groups{csvi},'t');

   [pp,ff,ee] = spm_fileparts(groups{1}{si}); 
   if contains(pp,'OASIS3')
       %%
       subj = groups{csvi}; 
       
       tab{csvi}(1,2:4) = {'OASISID','DAY','RUN'}; 
       for si = 1:numel(subj)
           tab{csvi}{si+1,2} = spm_str_manip(subj{si},'hhht'); tab{csvi}{si+1,2} = tab{csvi}{si+1,2}(5:end);
           tab{csvi}{si+1,3} = spm_str_manip(subj{si},'hht');  tab{csvi}{si+1,3} = tab{csvi}{si+1,3}(6:end); 
           if contains(spm_str_manip(subj{si},'t'),'run')
             tab{csvi}{si+1,4} = spm_str_manip(subj{si},'t'); tab{csvi}{si+1,4} = tab{csvi}{si+1,4}(strfind(tab{csvi}{si+1,4},'run-')+[4 5]); 
           else 
             tab{csvi}{si+1,4} = '01';
           end
       end
           
       
   end
   
   
   FN = {'bone_num','bone_med','rGMV','rWMV','rCSF'}; 
   for fni = 1:numel(FN)
       tab{csvi}{1,end+1} = FN{fni};
       if size( out.(FN{fni}){csvi} , 2) > 1  
          tab{csvi}(2:numel(groups{csvi})+1,end) = num2cell(out.(FN{fni}){csvi}');
       else
          tab{csvi}(2:numel(groups{csvi})+1,end) = num2cell(out.(FN{fni}){csvi});
       end
   end
   
   [pp,ff,ee] = fileparts(Pcsv{csvi}); 
   resultcsv{csvi} = fullfile(pp,sprintf('%s_bonemeasures.csv',ff));
   cat_io_csv(resultcsv{csvi},tab{csvi});
end


