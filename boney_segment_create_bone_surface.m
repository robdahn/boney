function [Si, St, Stm, Sth, sROI] = boney_segment_create_bone_surface(Vo, ...
  Ym, Yc, Ybonepp, Ybonethick, Ybonemarrow, Yheadthick, Ya, Ymsk, out, job)
%create_bone_surface. Surface-bases processing pipeline.
% Final bone processing function that create the bone surfaces to extract
% values from the surrounding bone tissue. It uses the percentage possition
% map Ybonepp that runs in the middle of the bone and map the thickness
% map Ybonethick on it. 
% In additon the thickness map of the head is also mapped.
%
% 
%  [Si, St, Stm, Sth, sROI] = boney_segment_create_bone_surface(Vo, Ym, ...
%   Yc, Ybonepp, Ybonethick, Ybonemarrow, Yheadthick, Ya, Ymsk, out, job)
%
%  Si          .. bone intensity surface
%  St          .. bone thickness surface
%  Stm         .. head
%  Sth         .. head thickness surface
%  sROI        .. 
%
%  Vo          .. original file header
%  Ym          .. intensity normalized input image
%  Ybonepp     .. percentage map of the bone (0-head to 1-brain)
%  Ybonethick  .. bone thickness map
%  Ybonemarrow .. bone marrow map (masked normalized image)
%  Yheadthick  .. head thickness map
%  Ya          .. bone atlas map (major regions)
%  Ymsk        .. bone mask map (to avoid bad regions)
%  out         .. main results
%  job         .. main parameters
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

  vx_vol = sqrt(sum(Vo.mat(1:3,1:3).^2));
  
  % surface coordinate transformation matrix
  matlab_mm  = Vo.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % CAT internal space
 

  %% create surface 
  Ypp  = Ybonepp + 0; spm_smooth(Ypp,Ypp,4);  
  [Yboneppr,res] = cat_vol_resize(smooth3(Ybonepp),'reduceV',vx_vol,job.opts.reduce,6,'meanm'); %#ok<ASGLU>        ./min(vx_vol)
  txt = evalc(sprintf('[Yppc,CBS.faces,CBS.vertices] = cat_vol_genus0(Yboneppr,.5,0);')); %#ok<NASGU>  % try to use lower values to avoid running into the vene
  CBS.vertices = CBS.vertices .* repmat(res.vx_red,size(CBS.vertices,1),1) - repmat((res.vx_red-1)/2,size(CBS.vertices,1),1); %#ok<NODEF> 
  CBS = cat_surf_fun('smat',CBS,matlab_mm); % transform to mm 
  CBS.EC = size(CBS.vertices,1) + size(CBS.faces,1) - size(spm_mesh_edges(CBS),1);  
  saveSurf(CBS,out.P.central);

  % optimize surface for midbone position 
  % simple blurring 
  Vpp  = cat_io_writenii(Vo, Ypp , '', sprintf('%s.pp','skull') ,...
      'percentage bone position map', 'uint8', [0,1/255],[1 0 0],struct());
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d', out.P.central ,out.P.central, 16 );  % simple smoothing to remove stair artifacts - 8-16 iterations in red 2
  cat_system(cmd,0);
  if 0 % not working?
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,out.P.central ,out.P.central ,-0.00005, 0.00005, 100, 0.01); %#ok<UNRCH> % -0.05, 0.05
    cat_system(cmd,0);
  end
  CBS = loadSurf(out.P.central);
  if exist(Vpp.fname,'file'), delete(Vpp.fname); end


  %% extract values from the bone 
  % create a smoothed thickness map for the mapping
  Ybonethick2 = cat_vol_approx(Ybonethick .* (Ybonethick>1 & Ybonethick<100),'nh',1,4);
  Ybonethick(Ybonethick>100) = Ybonethick2(Ybonethick>100);
  spm_smooth(Ybonethick,Ybonethick2,2./vx_vol);

  % map smoothed thickness map to surface
  Si = CBS; Si.facevertexcdata = cat_surf_fun('isocolors',max(3,Ybonethick), CBS, matlab_mm); 
  M  = spm_mesh_smooth(Si);
  % add min
  Si.facevertexcdata = single(spm_mesh_smooth(M,double(Si.facevertexcdata),10));
  Si.facevertexcdata(Si.facevertexcdata<=1) = nan; 
  if ~strcmpi(spm_check_version,'octave') 
    Si = cat_surf_fun('approxnans',Si); % ################  currently not working under octave ...
  end

  cat_io_FreeSurfer('write_surf_data',out.P.thick,Si.facevertexcdata);

  %% map values
  % write bone marrow map 
  Vpp  = cat_io_writenii(Vo, Ybonemarrow , '', 'skull.marrow' ,'bone marrow', 'single', [0,1],[1 0 0],struct());
  estminmax = 1; 
  if estminmax
    %% estimate the minimum value > real bone
    Ybonemarrow3 = Ybonemarrow; Ybonemarrow3(Ybonepp==0 | Ybonepp==1) = cat_stat_nanmean(Ybonemarrow3(Ybonepp>0 & Ybonepp<1)); 
    Vppmin  = cat_io_writenii(Vo, Ybonemarrow3 , '', 'skull.bone' ,'bone', 'single', [0,1],[1 0 0],struct());
    mappingstr = sprintf('-linear -min -steps "10" -start "-.5" -end ".5" -thickness "%s" ', out.P.thick); % weighted_avg
    cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vppmin.fname , out.P.marrowmin );
    cat_system(cmd,0); delete(Vppmin.fname);
    %% estimate the maximum value > bone marrow
    mappingstr = sprintf('-linear -max -steps "10" -start "-0.25" -end "0.25" -thickness "%s" ', out.P.thick); % weighted_avg
    cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrowmax );
    cat_system(cmd,0);
    marrowmin = cat_io_FreeSurfer('read_surf_data',out.P.marrowmin);
    marrowmax = cat_io_FreeSurfer('read_surf_data',out.P.marrowmax);
  end
  %%
  mappingstr = sprintf('-linear -weighted_avg -steps "10" -start "-.25" -end ".25" -thickness "%s" ', out.P.thick); % weighted_avg
  cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrow );
  cat_system(cmd,0);
  %%
  Si.facevertexcdata = cat_io_FreeSurfer('read_surf_data',out.P.marrow);
  % get atlas information
  Satlas = cat_surf_fun('isocolors',Ya, CBS, matlab_mm,'nearest');
  Si.facevertexcdata = Si.facevertexcdata .* (Satlas>0);
  cat_io_FreeSurfer('write_surf_data',out.P.thick,Si.facevertexcdata);


  if 0
    %% just for tests
    Si.facevertexcdata = cat_io_FreeSurfer('read_surf_data',out.P.marrowmax); %#ok<UNRCH> 
    %Si.facevertexcdata = single(spm_mesh_smooth(M,double(Si.facevertexcdata),2));
    Si.facevertexcdata = Si.facevertexcdata .* double(cat_io_FreeSurfer('read_surf_data',out.P.marrowmin));
    Sht = cat_surf_render2(Si); title(sprintf('bone thickness %s/%s',spm_str_manip(pp,'t'),ff));
    cat_surf_render2('Clim',Sht,[0 8]); cat_surf_render2('Colorbar',Sht);
  end
  %%
  if 0
    %% fat
    Yppi = max(0,min(1, smooth3(cat_vol_morph(Ybone<.5,'e') + min(0,-Ydiv*4)) + smooth3(Ybonepp))); %#ok<UNRCH> 
    Vpp  = cat_io_writenii(Vo, Yppi , '', sprintf('%s.pp','skull') ,...
        'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values', 'uint8', [0,1/255],...
        min([1 1 2],[1 0 0]),struct());
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.1"  "0.1"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,Pouter,Pouter,-.1, .1, 50, 0.01); 
    cat_system(cmd,job.opts.opts.verb-3);
    CBSo = loadSurf(Pouter);
  end




  %% get peak bone threshold
  % -use of nan for unwanted rois ?
  Si.vertices = CBS.vertices; Si.faces = single(CBS.faces); St = Si; Stm = Si; Sth = Si; Smsk = Si; 
  S.atlas     = cat_surf_fun('isocolors',Ya         , CBS, matlab_mm, 'nearest');
  S.thick     = cat_surf_fun('isocolors',Ybonethick , CBS, matlab_mm);
  S.hdthick   = cat_surf_fun('isocolors',Yheadthick , CBS, matlab_mm);
  St.facevertexcdata = S.thick;
  Sth.facevertexcdata = S.hdthick; 
  if ~strcmpi(spm_check_version,'octave') 
    St.facevertexcdata(St.facevertexcdata<=2) = nan; 
    Sth.facevertexcdata(Sth.facevertexcdata<=2) = nan; 
    St  = cat_surf_fun('approxnans',St);
    Si  = cat_surf_fun('approxnans',Si);
  end
  Stm.facevertexcdata = S.thick .* cat_surf_fun('isocolors',max(.1,1 - (cat_vol_grad(Ya*1000)>0.1) * .9 ),CBS, matlab_mm); 
  Smsk.facevertexcdata( cat_surf_fun('isocolors', single( Ymsk ) ,CBS, matlab_mm) > 1.5) = nan; 
  Sth.facevertexcdata( isnan( Smsk.facevertexcdata) ) = nan; 
  if ~isempty( job.opts.Pmask{1} )
    mskROIs = 0; % define unwanted ROIs
    for mski = mskROIs, Si.facevertexcdata(S.atlas==mski) = nan; end
  end
  cat_io_FreeSurfer('write_surf_data',out.P.thick,St.facevertexcdata);
  cat_io_FreeSurfer('write_surf_data',out.P.headthick,St.facevertexcdata);
 
  
  if 0
    %% just to test surfaces and values here
    Sx.vertices = CBS.vertices; Sx.faces = single(CBS.faces); 
    switch 6 %#ok<UNRCH> 
      case 0, Sx.facevertexcdata = cat_surf_fun('isocolors',Ybonemarrow, CBS, matlab_mm); 
      case 1, Sx.facevertexcdata = cat_surf_fun('isocolors',Ybonethick,CBS, matlab_mm); 
      case 2, Sx.facevertexcdata = cat_surf_fun('isocolors',Ym .* max(0, 1 + (cat_vol_grad(Ya*1000)>0.1) * 10 ),CBS, matlab_mm); 
      case 3, Sx.facevertexcdata = cat_surf_fun('isocolors',Ya,CBS, matlab_mm,'nearest'); 
      case 4, Sx.facevertexcdata = cat_surf_fun('isocolors',Ybgdist - Yheaddist,CBS, matlab_mm,'nearest'); 
      case 5, Sx.facevertexcdata = cat_surf_fun('isocolors',Yheadthick,CBS, matlab_mm); 
      case 6, Sx.facevertexcdata = Sth.facevertexcdata; 
    end  

    Sht = cat_surf_render2(Sx); title(sprintf('bone thickness %s/%s',strrep(spm_str_manip(out.P.pp,'t'),'_','\_'),out.P.orgff));
    cat_surf_render2('Clim',Sht,[0 4]); cat_surf_render2('Colorbar',Sht);
  end

  
  %% global + regional measures as column elements
  %  ----------------------------------------------------------------------
  rii = 1;
  sROI.help = 'ROI=0 is defined masked global values excluding the lower parts of the skull, whereas all other ROIs are without masking';
  for ri = 0:max(Ya(Ya(:)<intmax('uint16'))) % region ri==0 is the global value
    if ri == 0 || isnan(ri)
      ri = 0; %#ok<FXSET> % case of failed atlas mapping 
      sROI.boneatlas_id(1,rii)              = inf;
      sROI.nonnanvol(1,rii)                 = sum(S.atlas>intmax('uint16')) ./ numel(S.atlas);
      if ~isempty( job.opts.Pmask{1} ), sROI.boneatlas_name{1,rii} = 'full-masked'; 
      else,                             sROI.boneatlas_name{1,rii} = 'full-unmasked'; 
      end
      % bone marrow intensity 
      sROI.bonemarrow_mean(1,rii)           = cat_stat_nanmean(Si.facevertexcdata .* Smsk.facevertexcdata); 
      sROI.bonemarrow_std(1,rii)            = cat_stat_nanstd(Si.facevertexcdata .* Smsk.facevertexcdata); 
      sROI.bonemarrow_med(1,rii)            = cat_stat_nanmedian(Si.facevertexcdata .* Smsk.facevertexcdata); 
      sROI.bonemarrow_iqr(1,rii)            = iqr(Si.facevertexcdata .* Smsk.facevertexcdata); 
      % bone thickness
      sROI.bonethickness_mean(1,rii)        = cat_stat_nanmean(S.thick); 
      sROI.bonethickness_std(1,rii)         = cat_stat_nanstd(S.thick); 
      sROI.bonethickness_med(1,rii)         = cat_stat_nanmedian(S.thick); 
      sROI.bonethickness_iqr(1,rii)         = iqr(S.thick); 
      % head thickness
      sROI.headthickness_mean(1,rii)        = cat_stat_nanmean(S.hdthick); 
      sROI.headthickness_std(1,rii)         = cat_stat_nanstd(S.hdthick); 
      sROI.headthickness_med(1,rii)         = cat_stat_nanmedian(S.hdthick); 
      sROI.headthickness_iqr(1,rii)         = iqr(S.thick); 
      if estminmax
        sROI.bonemarrowmin_mean(1,rii)      = cat_stat_nanmean(marrowmin); 
        sROI.bonemarrowmin_std(1,rii)       = cat_stat_nanstd(marrowmin); 
        sROI.bonemarrowmin_med(1,rii)       = cat_stat_nanmedian(marrowmin); 
        sROI.bonemarrowmin_iqr(1,rii)       = iqr(marrowmin); 
        sROI.bonemarrowmax_mean(1,rii)      = cat_stat_nanmean(marrowmax); 
        sROI.bonemarrowmax_std(1,rii)       = cat_stat_nanstd(marrowmax); 
        sROI.bonemarrowmax_med(1,rii)       = cat_stat_nanmedian(marrowmax); 
        sROI.bonemarrowmax_iqr(1,rii)       = iqr(marrowmax); 
      end
      rii = rii + 1;
    else
      if sum(S.atlas==ri)>0
        sROI.boneatlas_id(1,rii)            = ri;  
        sROI.boneatlas_name{1,rii}          = sprintf('ROI%d',ri); 
        sROI.nonnanvol(1,rii)               = sum(S.atlas==ri) ./ numel(S.atlas);
        % bone marrow intensity 
        sROI.bonemarrow_mean(1,rii)         = cat_stat_nanmean(Si.facevertexcdata(S.atlas==ri)); 
        sROI.bonemarrow_std(1,rii)          = cat_stat_nanstd(Si.facevertexcdata(S.atlas==ri)); 
        sROI.bonemarrow_med(1,rii)          = cat_stat_nanmedian(Si.facevertexcdata(S.atlas==ri)); 
        sROI.bonemarrow_iqr(1,rii)          = iqr(Si.facevertexcdata(S.atlas==ri)); 
        % thickness
        sROI.bonethickness_mean(1,rii)      = cat_stat_nanmean(S.thick(S.atlas==ri)); 
        sROI.bonethickness_std(1,rii)       = cat_stat_nanstd(S.thick(S.atlas==ri)); 
        sROI.bonethickness_med(1,rii)       = cat_stat_nanmedian(S.thick(S.atlas==ri)); 
        sROI.bonethickness_iqr(1,rii)       = iqr(S.thick(S.atlas==ri)); 
        % head thickness
        sROI.headthickness_mean(1,rii)      = cat_stat_nanmean(S.hdthick(S.atlas==ri)); 
        sROI.headthickness_std(1,rii)       = cat_stat_nanstd(S.hdthick(S.atlas==ri)); 
        sROI.headthickness_med(1,rii)       = cat_stat_nanmedian(S.hdthick(S.atlas==ri)); 
        sROI.headthickness_iqr(1,rii)       = iqr(S.thick(S.atlas==ri)); 
        if estminmax
          sROI.bonemarrowmin_mean(1,rii)    = cat_stat_nanmean(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmin_std(1,rii)     = cat_stat_nanstd(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmin_med(1,rii)     = cat_stat_nanmedian(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmin_iqr(1,rii)     = iqr(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmax_mean(1,rii)    = cat_stat_nanmean(marrowmax(S.atlas==ri)); 
          sROI.bonemarrowmax_std(1,rii)     = cat_stat_nanstd(marrowmax(S.atlas==ri)); 
          sROI.bonemarrowmax_med(1,rii)     = cat_stat_nanmedian(marrowmax(S.atlas==ri)); 
          sROI.bonemarrowmax_iqr(1,rii)     = iqr(marrowmax(S.atlas==ri)); 
        end
        rii = rii + 1;
      end
    end
  end
end
% helping functions to load/write gifti surfaces
%=======================================================================
function saveSurf(CS,P)
  save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),P,'Base64Binary'); %#ok<USENS> 
end
%=======================================================================
function CS1 = loadSurf(P)
  CS = gifti(P);
  CS1.vertices = CS.vertices; CS1.faces = CS.faces; 
  if isfield(CS,'cdata'), CS1.cdata = CS.cdata; end
end
%=======================================================================
