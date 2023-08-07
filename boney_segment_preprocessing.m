function P = boney_segment_preprocessing(P,out,method,bias,rerun)
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________
  
  % check for already processed files
  PC = P; PC{1} = out(1).P.org; for i=2:numel(PC), PC{i}  = out(i).P.org; end
  if rerun
    PC = P; 
  else
    % extract prcoessed filenames
    Ppc = P; Ppc{1} = out(1).P.cls{1}; for i=2:numel(P), Ppc{i} = out(i).P.cls{1}; end
    
    % have to use CAT developer GUI for rerun function 
    oldexpertgui =  cat_get_defaults('extopts.expertgui');
    cat_get_defaults('extopts.expertgui',2);
    rpc = cat_io_rerun(PC,Ppc,0) > 0;
    cat_get_defaults('extopts.expertgui',oldexpertgui);
    PC  = PC(rpc); 
  end

  %% processing
  if ~isempty(PC)
    switch method
      case {1,'spm'}
        matlabbatch = SPM_preprocessing(PC, bias);
      case {2,'cat'} % not working
        matlabbatch = CAT_preprocessing(PC, 0, bias); %nproc
      otherwise
        error('unkown preprocessing method')
    end
  
    % data optimization 
    

    % run SPM batch
    spm_jobman('run',matlabbatch); 
  end


end
% subfunction with the main CAT matlabbatch
function matlabbatch = SPM_preprocessing(P,bias)
  matlabbatch{1}.spm.spatial.preproc.channel.vols       = P;
  matlabbatch{1}.spm.spatial.preproc.channel.biasreg    = 0.001;
  matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm   = 60 - 30*bias;
  matlabbatch{1}.spm.spatial.preproc.channel.write      = [0 1]; % write bias corrected image
  ngaus = [1 1 2 4 3 2]; 
  for ti = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm    = ...
      {fullfile(spm('dir'),'tpm',sprintf('TPM.nii,%d',ti))};
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus  = ngaus(ti);
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [ti<6 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
  end
  matlabbatch{1}.spm.spatial.preproc.warp.mrf           = 1;
  matlabbatch{1}.spm.spatial.preproc.warp.cleanup       = 1;
  matlabbatch{1}.spm.spatial.preproc.warp.reg           = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.spatial.preproc.warp.affreg        = 'subj';
  matlabbatch{1}.spm.spatial.preproc.warp.fwhm          = 0;
  matlabbatch{1}.spm.spatial.preproc.warp.samp          = 5;
  matlabbatch{1}.spm.spatial.preproc.warp.write         = [0 0];
  matlabbatch{1}.spm.spatial.preproc.warp.vox           = NaN;
  matlabbatch{1}.spm.spatial.preproc.warp.bb            = [NaN NaN NaN; NaN NaN NaN];
end
function matlabbatch = CAT_preprocessing(P, nproc, bias)
% subfunction with the main CAT matlabbatch
  matlabbatch{1}.spm.tools.cat.estwrite.data                              = P;
  matlabbatch{1}.spm.tools.cat.estwrite.data_wmh                          = {''};
  matlabbatch{1}.spm.tools.cat.estwrite.nproc                             = nproc;
  matlabbatch{1}.spm.tools.cat.estwrite.useprior                          = '';
  matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm                          = {'/Users/dahnke/Documents/MATLAB/spm12g/tpm/TPM.nii'};
  matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg                       = 'sub';
  matlabbatch{1}.spm.tools.cat.estwrite.opts.ngaus                        = [1 1 2 3 4 2];
  matlabbatch{1}.spm.tools.cat.estwrite.opts.warpreg                      = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.tools.cat.estwrite.opts.bias.biasstr                 = 0.5 + 0.5*bias;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.acc.accstr                   = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.redspmres                    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.3];
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.setCOM       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP          = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.affmod       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr        = 0; % faster ? 
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr       = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASmyostr    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr      = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr   = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr       = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC         = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC          = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf          = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.T1           = {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/T1.nii')};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.brainmask    = {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/brainmask.nii')};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.cat12atlas   = {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/cat.nii')};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm    = {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/Template_1_Dartel.nii')};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm  = {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/Template_0_GS.nii')};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr       = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.bb           = 12;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox          = 1.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres            = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod         = 'pbt2x';
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.SRP               = 30;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist             = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex      = 0.7;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp      = 0.1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp    = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental        = 0; 
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy                = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors        = 1; %
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb                = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno                = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.surface                    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures              = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics                       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40                                   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra                                    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers                                  = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus                                 = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamic_nuclei                          = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3                                 = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain                              = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Tian_Subcortex_S4_7T                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas                                 = {''};
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 1; % GM
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 1; % WM
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native     = 1; % CSF
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native    = 1; % head
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.native   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel    = 1; % image (m*)
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.native     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [0 0];
  matlabbatch{1}.spm.tools.cat.estwrite.output.rmat           = 0;
end