function [Vo, Yo, Yc, Ya, Ymsk, Ym, Affine, YaROIname, RES, BB] = ...
  boney_segment_loadMRI(P, job, seg8t, tis, cls, bd)
% loadMRI. Load MRI and segmentation maps with general size limit.
%
%  [Vo,Yo,Yc,Ya,Ym, Affine, YaROIname, RES, BB] = ...
%    boney_segment_loadMRI(P, job, seg8t, tis, cls, bd)
%
%  P             .. structure with prepared filenames of this subject
%  job           .. main SPM job structure (for affreg options)
%   .opts.normCT .. how normalize CT data
%   .opts.affreg .. do affine registration
%   .opts.Patlas .. atlas map path
%   .opts.Pmask  .. mask map path
%   .opts.reslim .. general resolution limit
%  seg8t         .. SPM mat structure (with Affine matrix and intensities)
%   .isCTseg     .. use of CTseg for segmentation
%   .tpm         .. for affine registration
%   .Affine      .. Affine registration matrix
%   .image       .. main image header
%  tis           .. our tissue intensity structure for intensity normalization
%   .res_vx_vol  .. voxel properties
%   .seg8o       .. tissue intensity values
%  cls           .. used classes (fast approach only load class 4; default=1:5)
%  bd            .. brain distance (need to limit the extraction of values
%                   in general; default=25)
%
%  Vo            .. original image header
%  Yo            .. original image
%  Yc            .. segment class images (cell)
%  Ya            .. atlas image
%  Ymsk          .. head mask (to avoid facial bones in global estimation)
%  Ym            .. intensity normalized image
%  Affine        .. (reprocessed) affine transformation to MNI
%  YaROIname     .. names of ROIs of the atlas
%  RES           .. resolution structure to avoid ultra-high resolutions
%  BB            .. boundary box to temporary remove the background
%
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


% TODO:
%  (1) tissue-based intensity scaling (currently only BG-WM based)
%  (2) use/eval affreg? (current affreg focus on all tissues, i.e. it's ok)


  if ~exist('cls'   ,'var'), cls    = 1:5; end
  if ~exist('bd'    ,'var'), bd     = 25;  end

  if exist(P.boneyPPmat,'file')
    try %#ok<TRYNC> % just standard load if not exist
      load(P.boneyPPmat); 
    end
  end
  if ~exist('Yc','var') || ~exist('Yo','var') || ~exist('Vc','var') || ~exist('Vo','var')
    % get bias corrected original image
    Vo = spm_vol(P.bc);
    Yo = single(spm_read_vols(Vo));

    % load segmentation
    Pc = cell(1,6); Yc = cell(1,6); Yc{6} = ones(Vo.dim,'single');
    for ci = cls
      if seg8t.isCTseg % CTseg
        Pc{ci}  = fullfile(P.orgpp,sprintf('c%02d%s%s',ci,P.ppff(4:end),P.ee));
      else
        Pc{ci}  = P.cls{ci}; %fullfile(P.orgpp,sprintf('c%d%s%s',ci,P.orgff,P.ee));
      end
      Vc(ci)  = spm_vol(Pc{ci}); %#ok<AGROW>
      Yc{ci}  = single(spm_read_vols(Vc(ci)));
      Yc{6}   = Yc{6} - Yc{ci};
    end
  
    % SPM/CAT segmentation
    if job.output.writeseg == 2
      save(P.boneyPPmat,'Yc','Yo','Vc','Vo')
    end
  end


  % create a linear intensity normalized image
  % .. unclear side effects ... and the histogram is not looking nice ?
  % .. maybe use some log scaling based approach later
% ###############################
  if 0 %~isempty(Pa)
    minimg = min( Yo(:) );
    maximg = max( Yo(:) );
    mintis = min( tis.seg8o );
    maxtis = max( tis.seg8o(2)*1.5 , max( seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.01)) ); % fat=maxhead
    switch tis.weighting
      case 1 % T1
        isc   = 1;
        T3th  = [minimg mintis tis.seg8o(3) tis.seg8o(1) tis.seg8o(2) tis.seg8o(2)+diff(tis.seg8o(1:2)) maxtis maximg];
        T3thx = [0 0.05 1 2 3 4 5 6];
        T3th  = interp1(T3th ,1:1/isc:numel(T3th )*isc,'pchip'); T3th  = smooth(T3th ,16*isc); %spm_smooth(T3th ,T3th ,.2*isc)
        T3thx = interp1(T3thx,1:1/isc:numel(T3thx)*isc,'pchip'); T3thx = smooth(T3thx,16*isc); %spm_smooth(T3thx,T3thx,.2*isc)
        Ym    = cat_main_gintnorm(Yo,struct('T3th',T3th,'T3thx',T3thx));
      case 2 % T2
        Ym = cat_main_gintnorm(Yo,struct('T3th',[minimg mintis tis.seg8o(2) ...
          tis.seg8o(1) tis.seg8o(3) maxtis maximg],'T3thx',[0 0.05 1 2 3 5 6]));
      case 3
        Ym = (Yo - min([ max([-.5 minimg]) tis.seg8o ])) / (tis.seg8o(2) - min([ max([-.5 minimg]) tis.seg8o ]));
      otherwise
        Ym = (Yo - min([ 0 tis.seg8o(3),tis.seg8o(end)])) / (tis.seg8o(2) - min([ 0 tis.seg8o(3),tis.seg8o(end)]));
    end
  elseif job.opts.fmethod==2
  % Not required in case of CAT preprocessing
    Ym = Yo;
  else
  % just a simple BG/WM based normalization
    if tis.weighting == 2 % MT
      Ym = (Yo - min([-.5 tis.intnorm(1) ])) / ( tis.intnorm(2) - min([-.5 tis.intnorm(1) ]));
    elseif tis.weighting == -1 % CT
      if job.opts.normCT %|| ~seg8t.isCTseg
        Ym = tis.fnormCT(Yo); 
        Yo = Ym; 
      else
        Ym = Yo;
      end
    else
      Ym = (Yo - tis.intnorm(1) ) / (tis.intnorm(2)  - tis.intnorm(1) );
    end
  end



  % == do affine registration ==
  % ##################### & job.opts.refine ????
  if job.opts.affreg > 0
    VG            = seg8t.tpm(1);
    if ~exist('Ytpmbrain','var')
      Ytpmbrain = spm_read_vols(seg8t.tpm(1)) +  spm_read_vols(seg8t.tpm(2)) +  spm_read_vols(seg8t.tpm(3));
    end
    VG.dat(:,:,:) = single(Ytpmbrain);
    VG.dt         = 16;
    VG.pinfo      = repmat([1;0],1,size(VG,3));
    VG            = cat_spm_smoothto8bit(VG,6);

    VF            = spm_vol(seg8t.image(1));
    VF.dat(:,:,:) = single(Yc{1} + Yc{2} + Yc{3});
    VF.dt         = 16;
    VF.pinfo      = repmat([1;0],1,size(VF,3));
    VF            = cat_spm_smoothto8bit(VF,6);

    if job.opts.affreg == 2
      evalc('Affine_com  = cat_vol_set_com(VF);'); % avoid output
      Affine_com(1:3,4) = -Affine_com(1:3,4); %#ok<NODEF>
    elseif job.opts.affreg == 3
      Affine_com = eye(4);
    else
      Affine_com = seg8t.Affine;
    end

    % prepare affine parameters
    aflags = struct('sep',12, ... max(6,max(sqrt(sum(VG(1).mat(1:3,1:3).^2)))), ...
      'regtype','subj','WG',[],'WF',[],'globnorm',1); % job.job.opts.opts.affreg
    warning off
    Affine  = spm_affreg(VG, VF, aflags, Affine_com);
    warning on
  elseif job.opts.affreg<0 || isempty(seg8t.Affine) || all(all(seg8t.Affine==eye(4)))
    %%
    VF            = spm_vol(seg8t.image(1));
    VF.dat(:,:,:) = single(Yc{1} + Yc{2} + Yc{3});
    VF.dt         = 16;
    VF.pinfo      = repmat([1;0],1,size(VF,3));
    VF            = cat_spm_smoothto8bit(VF,6); %#ok<NASGU>

    if abs(job.opts.affreg) == 2 || isempty(seg8t.Affine) || all(all(seg8t.Affine==eye(4)))
      evalc('Affine_com  = cat_vol_set_com(VF);'); % avoid output
      Affine_com(1:3,4) = -Affine_com(1:3,4); %#ok<NODEF>
    elseif abs(job.opts.affreg) == 3
      Affine_com = eye(4);
    else
      Affine_com = seg8t.Affine;
    end
    warning off
    Affine = spm_maff8(Vo,4,16,seg8t.tpmA,Affine_com,'subj',80);
    Affine = spm_maff8(Vo,4,4 ,seg8t.tpmA,Affine    ,'subj',40);
    warning on
  else
    Affine = seg8t.Affine;
  end
% ##############################
% quantify qc of affine registration (and update)



% #########
% - add non-linear mapping
% ##########

  % load atlas in individual space by applying the affine transformation
  if ~isempty(job.opts.Patlas)
    for ai = 1:numel( job.opts.Patlas )
      if ~isempty(job.opts.Patlas{1})
        Va     = spm_vol(job.opts.Patlas{ai});
        Ya{ai} = zeros(size(Ym),'single');
        for zi = 1:size(Ym,3)
          Ya{ai}(:,:,zi) = single(spm_slice_vol( Va , ...
            (Va.mat \ Affine * Vo.mat) * spm_matrix([0 0 zi]), ... % apply affine transformation
            [size(Ym,1), size(Ym,2)],[0,NaN])); % nearest neighbor interpolation
        end
        if max(Ya{ai}(:))==1, Ya{ai} = Ya{ai} + 1; Yareset(ai) = 1; else, Yareset(ai) = 0; end
        Ya{ai}(isnan(Ya{ai}(:))) = 0; 
        clear Va;
        [~,YD] = cat_vbdist(single(Ya{ai}>0.5),smooth3(Yc{6})<.5); Ya{ai} = Ya{ai}(YD);
        Pacsv{ai} = spm_file(job.opts.Patlas{ai},'ext','.csv');
        if exist(Pacsv{ai},'file')
          csv = cat_io_csv(Pacsv{ai});
          YaROIname{ai} = ['background';csv(2:end,2)];
        else
          YaROIname{ai} = unique(Ya{ai}(:));
        end
      else
        Ya{ai}        = zeros(size(Ym),'single');
        YaROIname{ai} = 0;
      end
    end
  else
    Ya{ai}        = zeros(size(Ym),'single');
    YaROIname{ai} = 0;
  end

  % load mask in individual space by applying the affine transformation
  if ~isempty(job.opts.Pmask{1})
    Vmsk = spm_vol(job.opts.Pmask{1});
    Ymsk = zeros(size(Ym),'single');
    for zi = 1:size(Ym,3)
      Ymsk(:,:,zi) = single(spm_slice_vol( Vmsk , ...
        (Vmsk.mat \ Affine * Vo.mat) * spm_matrix([0 0 zi]), ... % apply affine transformation
        [size(Ym,1),size(Ym,2)],[0,NaN])); % nearest neighbor interpolation
    end
    clear Vmsk;
    % fill zeros and select the upper part (region==1)
    [~,YD] = cat_vbdist(single(Ymsk>0),smooth3(Yc{6})<.5); Ymsk = Ymsk(YD);
    Ymsk = Ymsk<1.5; % this mask is limited ... ################## prepare masks ones for faster processing ! #############
  else
    Ymsk = false(size(Ym));
  end

  % extend atlas to all voxels
  if ~isempty(job.opts.Patlas{1}) || ~isempty(job.opts.Pmask)
    [~,YI] = cat_vbdist(single(Ya{ai}>0)); Ya{ai} = Ya{ai}(YI);
  end
  if Yareset(ai), Ya{ai} = Ya{ai} - 1; end


  % limit boundary box
  Yb  = ( Yc{1} + Yc{2} + Yc{3} ) >.5;
  [Yo,Ym,Ymsk,BB] = cat_vol_resize({Yo,Ym,Ymsk} ,'reduceBrain',tis.res_vx_vol,bd,Yb);
  for ai = 1:numel(Ya)
    Ya{ai} = cat_vol_resize(Ya{ai} ,'reduceBrain',tis.res_vx_vol,bd,Yb);
  end
  for ci = 1:numel(Yc)
    Yc{ci} = cat_vol_resize(Yc{ci} ,'reduceBrain',tis.res_vx_vol,bd,Yb);
  end

  % limit resolution
  [Yo,Ym,RES] = cat_vol_resize({Yo,Ym} ,'reduceV' ,tis.res_vx_vol,job.opts.reslim,16,'meanm');
  for ai = 1:numel(Ya)
    Ya{ai}    = cat_vol_resize(Ya{ai}  ,'reduceV' ,tis.res_vx_vol,job.opts.reslim,16,'nearest');
  end
  Ymsk        = cat_vol_resize(Ymsk    ,'reduceV' ,tis.res_vx_vol,job.opts.reslim,16,'meanm') > 0.5;
  Ysum        = zeros(size(Ym),'single'); 
  for ci = 1:numel(Yc)
    Yc{ci} = cat_vol_resize(Yc{ci},'reduceV' ,tis.res_vx_vol,job.opts.reslim,16,'meanm');
    Ysum   = Ysum + Yc{ci}; 
  end
  for ci = 1:numel(Yc)
    Yc{ci} = Yc{ci} ./ Ysum; 
  end

end
