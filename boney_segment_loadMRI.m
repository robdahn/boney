function [Vo, Yo, Yc, Ya, Ymsk, Ym, Affine, RES, BB] = boney_segment_loadMRI(P,opt,seg8t,tis,cls,bd) % TODO: (1) tissue-base intensity scaling; (2) use/eval affreg?
% loadMRI. Load MRI and segmentation.
%
%  [Vo,Yo,Yc,Ya,Ym, Affine] = loadMRI(P,Pa)
%
%  Po,Vo,Yo .. original image
%  Pc,Vc,Yo .. segment class images (cell)
%  Pa,Va,Ya .. atlas image
%

%opt.affreg, 1:5, opt.reslim, 25 ,opt); 
  
  if ~exist('cls'   ,'var'), cls    = 1:5; end 
  if ~exist('bd'    ,'var'), bd     = 25;  end 

  % get bias corrected original image
  Vo = spm_vol(P.bc);
  Yo = single(spm_read_vols(Vo));
  
  % load segmentation 
  Pc = cell(1,6); Yc = cell(1,6); Yc{6} = ones(Vo.dim,'single');
  for ci = cls % 1:5
    if seg8t.isCTseg % CTseg 
      Pc{ci}  = fullfile(P.orgpp,sprintf('c%02d%s%s',ci,P.ppff(4:end),P.ee));
    else
      Pc{ci}  = fullfile(P.orgpp,sprintf('c%d%s%s',ci,P.orgff,P.ee));
    end
    Vc(ci)  = spm_vol(Pc{ci}); 
    Yc{ci}  = single(spm_read_vols(Vc(ci)));
    Yc{6}   = Yc{6} - Yc{ci}; 
  end

  % create a linear intensity normalized image
  % .. unclear side effects ... and the histogram is not looking nice ?
  if 0 %~isempty(Pa)
    minimg = min( Yo(:) ); 
    maximg = max( Yo(:) );
    mintis = min( tis.seg8o );
    maxtis = max( tis.seg8o(2)*1.5 , tis.maxHead*tis.seg8o(2) ); 
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
  else
    if tis.weighting == 2 % MT
      Ym = (Yo - min([-.5 tis.seg8o ])) / (tis.seg8o(2) - min([-.5 tis.seg8o ]));
    elseif tis.weighting == -1 % CT
      if opt.normCT 
        Ym = (Yo - min( tis.seg8o )) / max(tis.seg8o(:) - min(tis.seg8o(:)));
      else
        Ym = Yo; 
      end
    else
      Ym = (Yo - min( tis.seg8o )) / (tis.seg8o(2) - min([tis.seg8o(3),tis.seg8o(end)]));
    end
  end
  


  % == do affine registration ==
  % ##################### & opt.refine ????
  if opt.affreg>0  
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
    
    if opt.affreg == 2  
      evalc('Affine_com  = cat_vol_set_com(VF);'); % avoid output
      Affine_com(1:3,4) = -Affine_com(1:3,4); %#ok<NODEF> 
    elseif opt.affreg == 3
      Affine_com = eye(4);
    else
      Affine_com = seg8t.Affine; 
    end
  
    % prepare affine parameter 
    aflags = struct('sep',12, ... max(6,max(sqrt(sum(VG(1).mat(1:3,1:3).^2)))), ...
      'regtype','subj','WG',[],'WF',[],'globnorm',1); % job.opt.opts.affreg
    warning off
    Affine  = spm_affreg(VG, VF, aflags, Affine_com); 
    warning on
  elseif opt.affreg<0 || isempty(seg8t.Affine) || all(all(seg8t.Affine==eye(4)))
    %%
    VF            = spm_vol(seg8t.image(1));
    VF.dat(:,:,:) = single(Yc{1} + Yc{2} + Yc{3});
    VF.dt         = 16; 
    VF.pinfo      = repmat([1;0],1,size(VF,3));
    VF            = cat_spm_smoothto8bit(VF,6); %#ok<NASGU> 

    if abs(opt.affreg) == 2 || isempty(seg8t.Affine) || all(all(seg8t.Affine==eye(4))) 
      evalc('Affine_com  = cat_vol_set_com(VF);'); % avoid output
      Affine_com(1:3,4) = -Affine_com(1:3,4); %#ok<NODEF> 
    elseif abs(opt.affreg) == 3
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




  % load atlas in individual space by appling the affine transformation
  if ~isempty(opt.output.Patlas{1})
    Va = spm_vol(opt.output.Patlas{1});
    Ya = zeros(size(Ym),'single'); 
    for zi = 1:size(Ym,3)
      Ya(:,:,zi) = single(spm_slice_vol( Va , ...
        (Va.mat \ Affine * Vo.mat) * spm_matrix([0 0 zi]), ... % apply affine tranformation
        size(Ym,1:2),[0,NaN])); % nearest neighbor interpolation 
    end
    clear Va;
    [~,YD] = cat_vbdist(single(Ya>0),smooth3(Yc{6})<.5); Ya = Ya(YD);
  else
    Ya   = zeros(size(Ym),'single'); 
  end
  
  % load mask in individual space by appling the affine transformation
  if ~isempty(opt.output.Pmask{1})
    Vmsk = spm_vol(opt.output.Pmask{1});
    Ymsk = zeros(size(Ym),'single'); 
    for zi = 1:size(Ym,3)
      Ymsk(:,:,zi) = single(spm_slice_vol( Vmsk , ...
        (Vmsk.mat \ Affine * Vo.mat) * spm_matrix([0 0 zi]), ... % apply affine tranformation
        size(Ym,1:2),[0,NaN])); % nearest neighbor interpolation 
    end
    clear Vmsk; 
    [~,YD] = cat_vbdist(single(Ymsk>0),smooth3(Yc{6})<.5); Ymsk = Ymsk(YD);
   % Ymsk = Ymsk>1.5; % this mask is limited ... ################## preapare masks ones for faster processing ! #############
  else
    Ymsk = false(size(Ym));
  end
  
  % extend atlas to all voxels
  if ~isempty(opt.output.Patlas{1}) || ~isempty(opt.output.Pmask)
    [~,YI] = cat_vbdist(single(Ya>0)); Ya = Ya(YI);   
  end  
  

  % limit boundary box
  Yb  = ( Yc{1} + Yc{2} + Yc{3} ) >.5; 
  [Yo,Ym,Ya,Ymsk,BB] = cat_vol_resize({Yo,Ym,Ya,Ymsk} ,'reduceBrain',tis.res_vx_vol,bd,Yb); 
  for ci = 1:numel(Yc)
    Yc{ci} = cat_vol_resize(Yc{ci} ,'reduceBrain',tis.res_vx_vol,bd,Yb);      
  end

  % limit resolution 
  [Yo,Ym,RES] = cat_vol_resize({Yo,Ym}  ,'reduceV' ,tis.res_vx_vol,opt.reslim,16,'meanm');
  [Ya,Ymsk]   = cat_vol_resize({Ya,Ymsk},'reduceV' ,tis.res_vx_vol,opt.reslim,16,'nearest');
  for ci = 1:numel(Yc)
    Yc{ci} = cat_vol_resize(Yc{ci},'reduceV' ,tis.res_vx_vol,opt.reslim,16,'meanm');
  end
end