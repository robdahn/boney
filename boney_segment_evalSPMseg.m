function [tismri, Ybraindist0] = boney_segment_evalSPMseg(Yo,Ym,Yc,Ymsk,vx_vol, fast, job, seg8t, tis) 
%evalSPMseg. Evaluation of the SPM segmentation. 
% This function extracts tissue intensities by the segmentation Yc from the
% intensity scaled image Ym. 
%
% [tismri, Ybraindist0] = boney_segment_evalSPMseg(Yo,Ym,Yc,Ymsk, ...
%   vx_vol,fast, opt, seg8t, tis)
%
%  Yo           .. original image
%  Ym           .. intensity normalized image
%  Yc           .. segment class images (cell)
%  Ymsk         .. head mask (to avoid facial bones in global estimation)
%  vx_vol       .. voxel-size
%  fast         .. use further tissue specific refinement or just extract  
%                  the median (default=0)
%  job          .. main SPM job structure (for affreg options) 
%   .opts.verb  .. be verbose
%  seg8t        .. spm8 structure without larger fields 
%   .isCTseg    .. special case in case of CT data
%  tis          .. measures based on SPM tissue thresholds and 
%                  basic information based on image and tissue properties
%   .seg8n      .. intensity peaks in Ym 
%
%  tismri       .. structure with MRI based measures
%  Ybraindist0  .. mask to limit the distance to the brain
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


% TODO: (1) test for WMHs; (2) 

  % == evaluate SPM segmentation ==
  % - estimate a brain distance map to further evaluate head/bone/background classes 
  Ybrain         = Yc{1} + Yc{2} + Yc{3}.^2; % underweight low CSF values 
  [Ybrainr,resR] = cat_vol_resize(Ybrain ,'reduceV' ,vx_vol,3,32,'meanm');  
  Ybraindist0r   = cat_vbdist( single( Ybrainr>.5 ) , true(size(Ybrainr)), resR.vx_volr); 
  Ybraindist0    = cat_vol_resize(Ybraindist0r ,'dereduceV' ,resR);         
  

  % test if the head-class has a lot of very low intensity values
  % the correction is done later if the refinement is selected
  if ~seg8t.isCTseg && tis.weighting < 0
    tismri.warning.clsUpdate = (sum( Yc{5}(:)>.5 & Ym(:)<.3 & Ybraindist0(:)<50 ) / sum( Yc{5}(:)>.5 & Ybraindist0(:)<50 ) ) > .2; 
    if tismri.warning.clsUpdate && 0 % ###### INACTIVE ######
      % create a (silent) warning
      try 
        cat_io_addwarning( [mfilename ':needSegUpdate'], ...
          sprintf('Bad SPM tissue class 5 - probably overestimated (%0.2f)', ...
          ( sum(Yc{5}(:)>.5 & Ym(:)<.1) ./ sum(Yc{5}(:)>.5 ) ) ) , 1, [1 0], [], 0, job.opts.verb>1);
      catch
        if job.opts.verb > 1 % print only for expert
          fprintf('\n');  
          cat_io_addwarning(  [mfilename ':needSegUpdate'], ...
            sprintf('Bad SPM tissue class 5 - probably overestimated (%0.2f)', ...
            ( sum(Yc{5}(:)>.5 & Ym(:)<.1) ./ sum(Yc{5}(:)>.5 ) ) ) , 1, [1 0], [], 0, job.opts.verb>1);
        end
      end
    end
  end
% ################### 
% early correction of fatal errors that are easy to correct for ?
% >> no, this is refined segmentation 
% ###################

% ###################
% * detect skull-stripping 
% * detect defacing 
% * detect WMHs? (as a second peak in WM class)
% * detect high intensity blood vessels
% ###################


  % separate fat and muscles
  os = round( mean(.9 ./ vx_vol)); 
  if seg8t.isCTseg %|| tis.weighting < 0% CT images
    Yfat = cat_vol_morph( Yc{5}>0.5 & Ym>-500 & Ym<0   & Ymsk>.5 ,'o',os); 
    Ymus = cat_vol_morph( Yc{5}>0.5 & Ym>0    & Ym<110 & Ymsk>.5 ,'o',os);
    Yhd  = cat_vol_morph( Yc{5}>0.5 & Ym>-500 & Ym<110 & Ymsk>.5 ,'o',os);
  elseif tis.weighting < 0% CT images
    %%
    Yhd  = cat_vol_morph((Yc{1} + Yc{2} + Yc{3} + Yc{4} + Yc{5})>.5 & Ym>-500,'ldo',2); 
    Yfat = Yc{5}>0.5 & Ym>-500 & Ym<-10 & Yhd; 
    Ymus = Yc{5}>0.5 & Ym>-10  & Ym<100  & Yhd;
  elseif tis.headFatType==2 % high bone intensity
    Yfat = cat_vol_morph( Yc{5}>0.5 & Ym>mean(tis.seg8n(2))                               & Ymsk>.5,'o',os); 
    Ymus = cat_vol_morph( Yc{5}>0.5 & Ym>mean(tis.seg8n([1,3])) & Ym<mean(tis.seg8n(1:2)) & Ymsk>.5,'o',os);
    Yhd  = cat_vol_morph( Yc{5}>0.5 & Ym>mean(tis.seg8n([1,3]))                           & Ymsk>.5,'o',os); 
  else % low fat intensity fat suppression
    Yfat = cat_vol_morph( Yc{5}>0.5 & (Ym>mean(tis.seg8n(3))     & Ym<mean(tis.seg8n(1))   |  Ym>mean(tis.seg8n(2))) & Ymsk>.5,'o',os); 
    Ymus = cat_vol_morph( Yc{5}>0.5 &  Ym>mean(tis.seg8n([1,3])) & Ym<mean(tis.seg8n(1:2))                           & Ymsk>.5,'o',os);
    Yhd  = cat_vol_morph( Yc{5}>0.5 & Ym>mean(tis.seg8n([1,3]))                                                      & Ymsk>.5,'o',os); 
  end
  


  % == evaluate/test the SPM segmentation == 
  tismri.help.TIV     = 'Total intracranial volume (GM + WM + CSF).';
  tismri.help.vol     = 'Volume of the SPM tissues classes in mm (probability >.5).';
  tismri.help.volr    = 'Relative volume of the SPM tissues classes (probability >.5) normalized by TIV.';
  tismri.help.volfat  = 'Volume of fat tissue in the masked upper head (simple threshold to separate the head tissues, in mm).'; 
  tismri.help.volfatr = 'Relative volume of fat tissue in the masked upper head (simple threshold to separate the head tissues).'; 
  tismri.help.volmus  = 'Volume of muscle-like tissue in the masked upper head (simple threshold to separate the head tissues, in mm).'; 
  tismri.help.volmusr = 'Relative volume of muscle-like tissue in the masked upper head (simple threshold to separate the head tissues).'; 
  tismri.help.clsQC   = 'Relation of voxel with high vs. low density within 30 mm distance. ';
  %tismri.help.Tth     = 'Median intensity of the (optimized) tissue class (~peak intensity).'; 
  %tismri.help.Tiqr    = 'IQR of the intensity of the (optimized) tissue class (~peak intensity).'; 
  tismri.help.int     = ['Intensity based evaluated tissue classes: ' ...
    '(1) GM:   median (=tismri.Tth(1)); ' ...
    '(2) WM:   median (=tismir.Tth(2)); ' ...
    '(3) CSF:  median (=tismri.Th(3)); ' ...
    '(4) bone: kmeans with 3 classes for bone cortex (1), low-intensity (2) and high-intensity bone marrow (3); ' ...
    '(5) head: kmeans with 3 classes head (1), muscle (2), and fat (3); ' ...
    '(6) bg:   median (=tismri.Th(6)). ' ...
    ];



  tismri.TIV     = sum( (Yc{1}(:) + Yc{2}(:) + Yc{3}(:)) > 0.5) .* prod(vx_vol) / 1000; 
  tismri.volfat  = cat_stat_nansum( Yfat(:) ) .* prod(vx_vol) / 1000; % tissue volume 
  tismri.volmus  = cat_stat_nansum( Ymus(:) ) .* prod(vx_vol) / 1000; % tissue volume 
  tismri.THV     = cat_stat_nansum( Yhd(:) )  .* prod(vx_vol) / 1000; % tissue volume 
 % tismri.volfatr = tismri.volfat ./ tismri.TIV;
 % tismri.volmusr = tismri.volmus ./ tismri.TIV;
  tismri.volfatr = tismri.volfat ./ tismri.THV;
  tismri.volmusr = tismri.volmus ./ tismri.THV;
  tismri.vol     = nan(1,6); tismri.volr = nan(1,6); 
  for ci = 1:6
    % estimate tissue volumes and density values
    tismri.vol(ci)  = cat_stat_nansum(Yc{ci}(:)>0.5) .* prod(vx_vol) / 1000; % tissue volume 
    tismri.volr(ci) = tismri.vol(ci) ./ tismri.TIV;
  
    % test if any class has more low probability values than high
    tismri.clsQC(ci) = sum(Yc{ci}(:)>.5) ./ sum(Yc{ci}(:)>eps & Yc{ci}(:)<.5 & Ybraindist0(:)<30); 
    if tismri.clsQC(ci)<.5
      % we add only a note as this one is not so relevant so far
      try
        cat_io_addwarning( sprintf('%s:badSPMcls%d',mfilename,ci) , ...
          sprintf('Bad SPM tissue class %d - probably underestimated (%0.2f)', ci, tismri.clsQC(ci)), ...
          0,[1 0],0,0,job.opts.verb>1);
      catch
        % for the old function only if highly verbose
        if job.opts.verb > 1
          fprintf('\n');  
          cat_io_addwarning( sprintf('%s:badSPMcls%d',mfilename,ci) , ...
            sprintf('Bad SPM tissue class %d - probably underestimated (%0.2f)', ci, tismri.clsQC(ci)), ...
            0,[1 0],0,0);
        end
      end
    end
 
    % refined evaluation with class optimization 
    % otherwise just the median
% ###################   
    switch ci .* (1-fast)
      case 1
        tismri.int.GM   = cat_stat_nanmedian(Yo(Yc{ci}>0.5));
        tismri.Tth(ci)  = tismri.int.GM;
      case 2
        % refinements to avoid PVE, WMHs and PVS  
        Yw                = cat_vol_morph(Yc{ci}>0.9 & cat_vol_morph(Ym>.95,'c'),'e'); 
        if sum(Yw(:)) < 1000, Yw = cat_vol_morph(Yc{ci}>0.5 & cat_vol_morph(Ym>.95,'c'),'e'); end 
        if sum(Yw(:)) < 1000, Yw = cat_vol_morph(Yc{ci}>0.5 & Ym>.9,'e'); end 
        if sum(Yw(:)) < 1000, Yw = Yc{ci}>0.5 & Ym>.9; end 
        tismri.int.WM     = cat_stat_nanmedian(Yo( Yw )); % avoid WMHs and PVE
        tismri.Tth(ci)    = tismri.int.WM;
        % ########################
        % WMHs ? 
        % ########################
        clear Yw; 
      case 3
        % refinements to avoid PVE 
        Yw                = cat_vol_morph(Yc{ci}>0.9 & Ym<.5,'e'); 
        if sum(Yw(:)) < 1000, Yw = Yc{ci}>0.9 | cat_vol_morph(Yc{ci}>0.5 & Ym<.5,'e'); end 
        if sum(Yw(:)) < 1000, Yw = Yc{ci}>0.9 & Ym<.5; end 
        tismri.int.CSF    = cat_stat_nanmedian( Yo( Yw )); % avoid PVE
        tismri.Tth(ci)    = tismri.int.CSF;
        clear Yw;
      case 4
        % ################### needs refinement depending on number ###########
        if seg8t.isCTseg || tis.weighting < 0
          tismri.int.bone         = cat_stat_kmeans(Yo(cat_vol_morph(Yc{ci}>0.5,'e')),2); 
          tismri.int.bone_cortex  = tismri.int.bone(2);
          tismri.int.bone_marrow  = tismri.int.bone(1);

          tismri.Tth(ci)          = mean(tismri.int.bone);
        else
          %%
          Ymx = Yc{4}>.25 & Ym<.5;
          Ybx = cat_vol_localstat(Ym ,Ymx,1,2,1);  
          Ybx = cat_vol_localstat(Ybx,Ymx,1,1,1);
          tismri.int.bone_cortex = median(Ybx(Ymx(:))); 
          
          Ymx = smooth3( Yc{4}>.5 & Ym>tismri.int.bone_cortex * 2 ) > .5;
          Ybx = cat_vol_localstat(Ym ,Ymx,1,3,2);  
          Ybx = cat_vol_localstat(Ybx,Ymx,1,1,2);
          tismri.int.bone_marrow = median(Ybx(Ymx(:) & Ym(:)>1));

          Ymx = smooth3( Yc{4}>.5 & Ym>tismri.int.bone_cortex * 2) > .5 & Ym<1.5;
          Ybx = cat_vol_localstat(Ybx,Ymx,2,1,2);
          mn  = cat_stat_kmeans(Ybx(Ymx(:) & Ybx(:)>0),3); 
          tismri.int.bone_marrow_red = mn(2);
          
          tismri.int.bone = cat_stat_kmeans(Yo(cat_vol_morph(Yc{ci}>0.5,'e')),3); 
          tismri.Tth(ci)  = tismri.int.bone(3);
        end
       
      case 5
        %% ################### needs refinement depending on number ###########
        Ybrainskulldist   = cat_vbdist( Yc{1} + Yc{2} + Yc{3} + Yc{4} ) / mean(vx_vol);
       
        Ymx = smooth3( Yc{5}>.5 & Ym > tismri.int.bone_cortex * 4 ) > .5 & Ybrainskulldist<60;
        Ybx = cat_vol_localstat(Ym ,Ymx,1,3,3);  
        Ybx = cat_vol_localstat(Ybx,Ymx,1,1,2);
        mn  = cat_stat_kmeans(Ybx(Ymx(:) & Ybx(:)>0),2); 
        tismri.int.head_fat = mn(2);
        tismri.head_fat = mn(2);
%%   
        Yg  = cat_vol_grad(Ym); 
        Ymx = smooth3( Yc{5}>.5 & Ym > tis.bonecortex * 2 & Ym<1 & Yg<.5 ) > .5 & Ybrainskulldist<60;
        Ybx = cat_vol_localstat(Ym ,Ymx,1,1,1);  
        mn  = cat_stat_kmeans(Ybx(Ymx(:)>.5),2); 
        if isnan(mn(2))
          tismri.int.head_muscle = tis.seg8mns(1);
          tismri.head_muscle     = tis.seg8mns(1);
        else
          tismri.int.head_muscle = mn(2);
          tismri.head_muscle     = mn(2);
        end
        %{
        tismri.int.head   = cat_stat_kmeans(Yo(cat_vol_morph(Yc{ci}>0.9,'e')),3); 
        tismri.Tth(ci)    = tismri.int.head(3);
        
          
        Ymus              = Yc{5}>.5 & Ym>mean([tis.CSF,tis.GM]) & Ym<.9; 
        Ymus              = cat_vol_morph(smooth3(Ymus)>.5,'o'); 
        if tis.headFatType == 0
          Yfat            = Yc{5}>.5 & Ym>( cat_stat_nanmedian(Ym(Ymus(:)))  +  2*cat_stat_nanmedian(Ym(Ymus(:))) ); 
        else
          if tis.highBG 
            Yfat          = Yc{5}>.5 & ~Ymus;
          else
            Yfat          = Yc{5}>.5 & ~Ymus & Ym>tis.background*2;
          end
        end
        tismri.head_muscle = cat_stat_nanmedian(Ym(Ymus(:) & Ybrainskulldist(:) < 60));
        tismri.head_fat    = cat_stat_nanmedian(Ym(Yfat(:) & Ybrainskulldist(:) < 60));
        %}

      otherwise
        [Yir,Ymr]         = cat_vol_resize( {Yo .* (Yc{ci}>.5),(Yc{ci}>.9)} ,'reduceV' ,1,2,32,'meanm'); 
        tismri.Tth(ci)    = cat_stat_nanmedian(Yir(Ymr>.5));
    end
    % SD in tissue class
    %{
      if ~exist('Yir','var')
        [Yir,Ymr]           = cat_vol_resize( {Yo .* (Yc{ci}>.5),(Yc{ci}>.9)} ,'reduceV' ,1,2,32,'meanm'); 
      end
      tismri.Tsd(ci)        = cat_stat_nanstd(Yir(Ymr>.5));
      tismri.Tiqr(ci)       = iqr(Yir(Ymr(:)>.5 & ~isnan(Yir(:))));
      clear Yir Ymr;
    %}
  end
  
  
end
