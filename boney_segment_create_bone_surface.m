function [Si, Stm, sROI] = boney_segment_create_bone_surface ...
  (Vo, Ybonepp, Ybonemarrow, Ybonethick, Yheadthick, Ya, Ymsk, YaROIname, out, job)
%create_bone_surface. Surface-bases processing pipeline.
% Final bone processing function that create the bone surfaces to extract
% values from the surrounding bone tissue. It uses the percentage possition
% map Ybonepp that runs in the middle of the bone and map the thickness
% map Ybonethick on it. 
% In additon the thickness map of the head is also mapped.
%
% 
%  [Si, Stm, sROI] = boney_segment_create_bone_surface(Vo, ...
%    Ybonepp, Ybonemarrow, Ybonethick, Yheadthick, Ya, Ymsk, out, job)
%
%  Si          .. bone intensity surface
%  Stm         .. bone thickness surface with regional boundaries
%  sROI        .. 
%
%  Vo          .. original file header
%  Ybonepp     .. percentage map of the bone (0-head to 1-brain)
%  Ybonemarrow .. bone marrow map (masked normalized image)
%  Ybonethick  .. bone thickness map
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

  
  % surface coordinate transformation matrix
  matlab_mm  = Vo.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % CAT internal space
  vx_vol     = sqrt(sum(Vo.mat(1:3,1:3).^2));
 

  %% == create surface ==
  %  Worked best for the already smooth Ybonepp map and additional strong 
  %  surface smoothing. Surface deformation was not helpful. 
  [Yboneppr,res] = cat_vol_resize(smooth3(Ybonepp),'reduceV',vx_vol,job.opts.reduce,6,'meanm'); %#ok<ASGLU>       
  txt = evalc(sprintf('[Yppc,CBS.faces,CBS.vertices] = cat_vol_genus0(Yboneppr,.5,0);')); %#ok<NASGU>  
  CBS.vertices = CBS.vertices .* repmat(res.vx_red,size(CBS.vertices,1),1) - repmat((res.vx_red-1)/2,size(CBS.vertices,1),1); %#ok<NODEF> 
  CBS = cat_surf_fun('smat',CBS,matlab_mm); % transform to mm 
  CBS.EC = size(CBS.vertices,1) + size(CBS.faces,1) - size(spm_mesh_edges(CBS),1);  
  saveSurf(CBS,out.P.central);

  % optimize surface for midbone position by simple blurring 
  % simple smoothing to remove stair artifacts - 8-16 iterations in red 2
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d', out.P.central ,out.P.central, 16 );  
  cat_system(cmd,0);
  CBS = loadSurf(out.P.central); 
   
  % create a (smoothed) thickness map for the mapping extract values from the bone
  % .. however, the smoothing was not improving the mapping
  Ybonethick2 = cat_vol_approx(Ybonethick .* (Ybonethick>1 & Ybonethick<100),'nh',1,3);
  Si = CBS; Si.facevertexcdata = cat_surf_fun('isocolors',max(3,Ybonethick2), CBS, matlab_mm); 
  cat_io_FreeSurfer('write_surf_data',out.P.thick,Si.facevertexcdata);

  

  % == map values ==
  % estimate the local minimum to get the hard bone
  bonemed      = cat_stat_nanmedian(Ybonemarrow(Ybonepp>0 & Ybonepp<1)); 
  Ybonemarrow3 = Ybonemarrow; Ybonemarrow3(Ybonepp==0 | Ybonepp==1) = bonemed; % limit by median
  Vppmin       = cat_io_writenii(Vo, Ybonemarrow3 , '', 'skull.bone' ,'bone', 'single', [0,1],[1 0 0],struct());
  mappingstr   = sprintf('-linear -min -steps "9" -start "-.5" -end ".5" -thickness "%s" ', out.P.thick); % min - larger range to assure minimum
  cmd          = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vppmin.fname , out.P.cortex );
  cat_system(cmd,0); delete(Vppmin.fname);
  cortex       = min( cat_stat_nanmedian(Ybonemarrow3(Ybonepp>0 & Ybonepp<1)) , cat_io_FreeSurfer('read_surf_data',out.P.cortex));
  
  % estimate maximum for bone marrow with full range
  % ######## RD20230831: just for tests.
  Vpp          = cat_io_writenii(Vo, Ybonemarrow , '', 'skull.marrow' ,'bone marrow', 'single', [0,1],[1 0 0],struct());
  mappingstr   = sprintf('-linear -max -steps "9" -start "-.5" -end ".5" -thickness "%s" ', out.P.thick); % weighted_avg
  cmd          = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrow );
  cat_system(cmd,0);
  marrowmax    = cat_io_FreeSurfer('read_surf_data',out.P.marrow);

  % estimate average for bone marrow with even limited range
  Vpp          = cat_io_writenii(Vo, Ybonemarrow , '', 'skull.marrow' ,'bone marrow', 'single', [0,1],[1 0 0],struct());
  mappingstr   = sprintf('-linear -weighted_avg -steps "5" -start "-.1" -end ".1" -thickness "%s" ', out.P.thick); % weighted_avg
  cmd          = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrow );
  cat_system(cmd,0);
  Si.facevertexcdata = cat_io_FreeSurfer('read_surf_data',out.P.marrow);

  % get atlas information
  Satlas = cat_surf_fun('isocolors',Ya, CBS, matlab_mm,'nearest');
  Si.facevertexcdata = Si.facevertexcdata .* (Satlas>0);
  cat_io_FreeSurfer('write_surf_data',out.P.thick,Si.facevertexcdata);

  if 0 %#ok<*UNRCH> 
    %% only for debugging and tests!
    Sh  = cat_surf_render2(Si);                                          cat_surf_render2('Clim',Sh,[0 6]); cat_surf_render2('Colorbar',Sh); title('soft') 
    Stx = Si; Stx.facevertexcdata = cortex; Sh = cat_surf_render2(Stx);  cat_surf_render2('Clim',Sh,[0 3]); cat_surf_render2('Colorbar',Sh);  title('hard')
  end



  %% get peak bone threshold
  % -use of nan for unwanted rois ?
  Si.vertices = CBS.vertices; Si.faces = single(CBS.faces); St = Si; Stm = Si; Sth = Si; 
  S.atlas     = cat_surf_fun('isocolors',Ya         , CBS, matlab_mm, 'nearest');
  S.mask      = cat_surf_fun('isocolors',Ymsk       , CBS, matlab_mm, 'nearest') == 1;
  S.thick     = cat_surf_fun('isocolors',Ybonethick , CBS, matlab_mm);
  S.hdthick   = cat_surf_fun('isocolors',Yheadthick , CBS, matlab_mm);
  St.facevertexcdata  = S.thick;
  Sth.facevertexcdata = S.hdthick; 

  % thickess with atlas borders
  Stm.facevertexcdata = S.thick .* cat_surf_fun('isocolors',max(.1,1 - (cat_vol_grad(Ya*1000)>0.1) * .9 ),CBS, matlab_mm);
  

  
  %% global + regional measures as column elements
  %  ----------------------------------------------------------------------
  %  - similar to the volume measures we also focus here on the global/
  %    regional mean values and ignore median, std, and iqr
  rii = 1;  
  sROI.help = 'ROI=0 is defined masked global values excluding the lower parts of the skull, whereas all other ROIs are without masking';
  for ri = 0:max(Ya(Ya(:)<intmax('uint16'))) % region ri==0 is the global value
    if ri == 0 || isnan(ri)
      ri = 0; %#ok<FXSET> % case of failed atlas mapping 
      sROI.boneatlas_id(1,rii)       = inf;
      sROI.nonnanvol(1,rii)          = sum(S.atlas>intmax('uint16')) ./ numel(S.atlas);
      if ~isempty( job.opts.Pmask{1} ), sROI.boneatlas_name{1,rii} = 'full-masked'; 
      else,                             sROI.boneatlas_name{1,rii} = 'full-unmasked'; 
      end
      sROI.bonemarrow(1,rii)         = cat_stat_nanmean(Si.facevertexcdata(S.mask)); 
      sROI.bonemarrowmax(1,rii)      = cat_stat_nanmean(marrowmax(S.mask)); 
      sROI.bonecortex(1,rii)         = cat_stat_nanmean(cortex(S.mask)); 
      sROI.bonethickness(1,rii)      = cat_stat_nanmean(S.thick(S.mask)); 
      sROI.headthickness(1,rii)      = cat_stat_nanmean(S.hdthick(S.mask)); 
      rii = rii + 1;
    else
      if sum(S.atlas==ri)>0
        sROI.boneatlas_id(1,rii)     = ri;  
        if isempty(YaROIname) %|| numel(YaROIname)>max(Ya(Ya(:)<intmax('uint16')))
          sROI.boneatlas_name{1,rii} = sprintf('ROI%d',ri); 
        else
          sROI.boneatlas_name{1,rii} = YaROIname{rii};
        end
        sROI.nonnanvol(1,rii)        = sum(S.atlas==ri) ./ numel(S.atlas);
        sROI.bonemarrow(1,rii)       = cat_stat_nanmean(Si.facevertexcdata(S.atlas==ri)); 
        sROI.bonemarrowmax(1,rii)    = cat_stat_nanmean(marrowmax(S.atlas==ri)); 
        sROI.bonecortex(1,rii)       = cat_stat_nanmean(cortex(S.atlas==ri)); 
        sROI.bonethickness(1,rii)    = cat_stat_nanmean(S.thick(S.atlas==ri)); 
        sROI.headthickness(1,rii)    = cat_stat_nanmean(S.hdthick(S.atlas==ri)); 
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
