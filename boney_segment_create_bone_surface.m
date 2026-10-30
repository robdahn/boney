function [Si, Stm, sROI] = boney_segment_create_bone_surface ...
  (Vo, Ybonepp, Ybonemarrow, Ybonethick, Yheadthick, Ya, Ymsk, YaROIname, YaROIid, out, job)
%create_bone_surface. Surface-based processing pipeline.
% Final bone processing function that creates the bone surfaces to extract
% values from the surrounding bone tissue. It uses the percentage position
% map Ybonepp that runs in the middle of the bone and maps the thickness
% map Ybonethick on it.
% In additon, the thickness map of the head is also mapped.
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
  vx_vol     = sqrt(sum(Vo.mat(1:3,1:3).^2)) ;


  %% == create surface ==
  %  Worked best for the already smooth Ybonepp map and additional strong
  %  surface smoothing. Surface deformation was not helpful.
  [Yboneppr,res] = cat_vol_resize(smooth3(Ybonepp),'reduceV',vx_vol,job.opts.reduce,6,'meanm'); %#ok<ASGLU>
  txt = evalc(sprintf('[Yppc,CBS.faces,CBS.vertices] = cat_vol_genus0(Yboneppr,.5,0);')); %#ok<NASGU>
  if isempty(CBS.vertices) %#ok<NODEF>
    error('boney:createBoneSurface:empty','Empty bone surface (failed bone segmentation).');
  end
  CBS.vertices = CBS.vertices .* repmat(res.vx_red * [0 1 0; 1 0 0; 0 0 1],size(CBS.vertices,1),1) ... 
    - repmat((res.vx_red * [0 1 0; 1 0 0; 0 0 1]-1)/2,size(CBS.vertices,1),1); %#ok<NODEF> 
  CBS = cat_surf_fun('smat',CBS,matlab_mm); % transform to mm
  CBS.EC = size(CBS.vertices,1) + size(CBS.faces,1) - size(spm_mesh_edges(CBS),1);
  saveSurf(CBS,out.P.central);
  
  % optimize surface for midbone position by simple blurring
  % simple smoothing to remove stair artifacts - 8-16 iterations in red 2
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d', out.P.central ,out.P.central, 16 );
  cat_system(cmd,0);

  CBS = loadSurf(out.P.central);

  % create a (smoothed) thickness map for the mapping extracted values from the bone
  % .. however, the smoothing did not improve the mapping
  Ybonethick2 = cat_vol_approx(Ybonethick .* (Ybonethick>1 & Ybonethick<100),'nh',1,3);
  Si = CBS; Si.facevertexcdata = cat_surf_fun('isocolors',max(3,Ybonethick2), CBS, matlab_mm);
  cat_io_FreeSurfer('write_surf_data',out.P.thick,Si.facevertexcdata);



  % == map values ==
  % estimate the local minimum to get the hard bone (cortex)
  bonemed      = cat_stat_nanmedian(Ybonemarrow(Ybonepp>0 & Ybonepp<1));
  Ybonemarrow3 = Ybonemarrow; Ybonemarrow3(Ybonepp==0 | Ybonepp==1) = bonemed; % limit by median
  Vppmin       = cat_io_writenii(Vo, Ybonemarrow3 , out.P.surfpath, 'skull.bone' ,'bone', 'single', [0,1],[1 0 0],struct());
  mappingstr   = sprintf('-min -steps "9" -start "-.5" -end ".5" -thickness "%s" ', out.P.thick); % min - larger range to assure minimum
  cmd          = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vppmin.fname , out.P.cortex );
  cat_system(cmd,0); delete(Vppmin.fname);
  cortex       = min( cat_stat_nanmedian(Ybonemarrow3(Ybonepp>0 & Ybonepp<1)) , cat_io_FreeSurfer('read_surf_data',out.P.cortex));

  %% estimate maximum for bone marrow with full range (marrowmax)
  %  and the average for bone marrow with even limited range 
  Vpp          = cat_io_writenii(Vo, Ybonemarrow , out.P.surfpath, 'skull.marrow' ,'bone marrow', 'single', [0,1],[1 0 0],struct());
  mappingstr   = sprintf('-max -steps "9" -start "-.5" -end ".5" -thickness "%s" ', out.P.thick); 
  cmd          = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrow );
  cat_system(cmd,0); 
  marrowmax    = cat_io_FreeSurfer('read_surf_data',out.P.marrow);
  mappingstr   = sprintf('-weighted_avg -steps "5" -start "-.25" -end ".25" -thickness "%s" ', out.P.thick); 
  cmd          = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrow );
  cat_system(cmd,0); delete(Vpp.fname);
  Si.facevertexcdata = cat_io_FreeSurfer('read_surf_data',out.P.marrow);

  % mask the bone marrow by the atlas (the thickness file is not changed)
  Satlas = cat_surf_fun('isocolors',Ya{1}, CBS, matlab_mm,'nearest');
  Si.facevertexcdata = single(Si.facevertexcdata .* (Satlas>0));
  cat_io_FreeSurfer('write_surf_data',out.P.marrow,Si.facevertexcdata);



  %% get peak bone threshold
  % -use of nan for unwanted rois ?
  Si.vertices = CBS.vertices; Si.faces = single(CBS.faces); St = Si; Stm = Si; Sth = Si;
  for ai = 1:numel(Ya)
    S.atlas{ai} = cat_surf_fun('isocolors',Ya{ai}      , CBS, matlab_mm, 'nearest');
  end
  S.mask      = cat_surf_fun('isocolors',single(Ymsk), CBS, matlab_mm) > .5;
  S.thick     = cat_surf_fun('isocolors',Ybonethick  , CBS, matlab_mm);
  S.hdthick   = cat_surf_fun('isocolors',Yheadthick  , CBS, matlab_mm);
  St.facevertexcdata  = S.thick;
  Sth.facevertexcdata = S.hdthick;

  % thickess with atlas borders
  Stm.facevertexcdata = S.thick .* cat_surf_fun('isocolors',max(.1,1 - (cat_vol_grad(Ya{1}*1000)>0.1) * .9 ),CBS, matlab_mm);



  %% global + regional measures as column elements
  %  ----------------------------------------------------------------------
  %  - similar to the volume measures we also focus here on the global/
  %    regional mean values and ignore median, std, and iqr
  clear sROI
  for ai = 1:numel(Ya)
    sROI(ai).file = '';
    sROI(ai).help = [
        'A masked image is used for global values to extract only the upper part of the skull, ' ...
        'whereas no masking is used in case of atlas regions. '];
    for ri = 0:numel( YaROIid{ai} )  % region ri==0 is the global value
      sROI(ai).boneatlas_nr(1,ri+1) = ri;  
      if ri == 0 
        sROI(ai).boneatlas_id(1,ri+1)       = inf;
        sROI(ai).nonnanvol(1,ri+1)          = sum(S.atlas{ai}>intmax('uint16')) ./ numel(S.atlas{ai});
        if ~isempty( job.opts.Pmask ) && ~isempty( job.opts.Pmask{1} )
          sROI(ai).boneatlas_name{1,ri+1}   = 'full-masked';
        else
          sROI(ai).boneatlas_name{1,ri+1}   = 'full-unmasked';
        end
        sROI(ai).bonemarrow(1,ri+1)         = cat_stat_nanmean(Si.facevertexcdata(S.mask));
        sROI(ai).bonemarrowmax(1,ri+1)      = cat_stat_nanmean(marrowmax(S.mask));
        sROI(ai).bonecortex(1,ri+1)         = cat_stat_nanmean(cortex(S.mask));
        sROI(ai).bonethickness(1,ri+1)      = cat_stat_nanmean(S.thick(S.mask));
        sROI(ai).headthickness(1,ri+1)      = cat_stat_nanmean(S.hdthick(S.mask));
      else
        sROI(ai).boneatlas_id(1,ri+1)     = YaROIid{ai}(ri);
        if iscell(YaROIname{ai}) && ri <= numel(YaROIname{ai}) && ~isempty(YaROIname{ai}{ri}) 
          sROI(ai).boneatlas_name{1,ri+1} = YaROIname{ai}{ri};
        else
          sROI(ai).boneatlas_name{1,ri+1} = sprintf('ROI%d',ri);
        end
        sROI(ai).bonemarrow(1,ri+1)       = cat_stat_nanmean(Si.facevertexcdata(S.atlas{ai}==YaROIid{ai}(ri)));
        sROI(ai).bonemarrowmax(1,ri+1)    = cat_stat_nanmean(marrowmax(S.atlas{ai}==YaROIid{ai}(ri)));
        sROI(ai).bonecortex(1,ri+1)       = cat_stat_nanmean(cortex(S.atlas{ai}==YaROIid{ai}(ri)));
        sROI(ai).bonethickness(1,ri+1)    = cat_stat_nanmean(S.thick(S.atlas{ai}==YaROIid{ai}(ri)));
        sROI(ai).headthickness(1,ri+1)    = cat_stat_nanmean(S.hdthick(S.atlas{ai}==YaROIid{ai}(ri)));
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
