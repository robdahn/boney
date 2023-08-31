function out = boney_getBoneMeasures(groups,opt)
% The function returnes cranial bone measures by using the output of the SPM 
% segmentation of the native tissue output of classes 1 to 5 (c#*.nii) and  
% bias corrected images (m*.nii).
%
% out = getBoneMeasures(groups,opt)
% 
% groups .. cell of cell with m*.ni
%
%Input:  modulated nifti images 
%Output: bone_num: bone volume
%        bone_med: median intensity of bone tissue
%        bone_std: SD of bone intensity
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


    if ~exist('opt','var'), opt = struct(); end
    def.matonly  = 0; 
    def.printimg = 1;
    def.printvol = 1;
    def.rerun    = 0; 
    opt = cat_io_checkinopt(opt,def); 
    
    % prepare output parameter
    out.bone_num = {}; out.bone_med = {}; out.bone_std = {}; 
    out.int_ci = {}; out.vol_ci = {}; 
    out.TIV = {}; out.con = {}; 
    for gi = 1:numel(groups)
        fprintf('Group %d/%d\n',gi,numel(groups)); 
        if opt.matonly > 9
        %% new approach
          clear opt2;
          opt2.files             = groups{gi}; 
          opt2.opts.bmethod      = min(2,opt.matonly - 10); 
          opt2.output.writevol   = opt.printvol * 0;
          opt2.output.writeseg   = opt.printvol;
          opt2.output.writesurf  = opt.printvol * 0; 
          opt2.output.report     = 3 * opt.printimg; 
          opt2.opts.subdirs      = 1; 
          opt2.opts.affreg       = 0; % 0-
          opt2.opts.reduce       = 2; %
          opt2.opts.rerun        = opt.rerun; 
          opt2.opts.verb         = 1; % 1-default, 2-details

          %%
          [~,nout] = boney_segment(opt2);

          %%
          if opt.matonly == 10 %  classic SPM values for all processings
            for si = 1:numel(groups{gi})
              out.bone_num{gi}(si,1)  = nout(si).tismri.volr(4);                    % bone volume normalised by TIV
              out.bone_med{gi}(si,1)  = nout(si).tismri.bone_med; % median bone intensity value
              out.bone_std{gi}(si,1)  = nan; 
              out.TIV{gi}(si,1)       = nout(si).tismri.TIV; 
              out.int_ci{gi}(si,:)    = nan(1,6); 
              out.vol_ci{gi}(si,:)    = nan(1,6); 
              out.con{gi}(si,1)       = nan;
            end
          elseif opt.matonly == 11 || opt.matonly == 12
            for si = 1:numel(groups{gi})
              out.bone_num{gi}(si,1)  = nout(si).classic.bone_num;                 % bone volume normalised by TIV
              out.bone_med{gi}(si,1)  = nout(si).classic.bone_med; % median bone intensity value
              out.bone_std{gi}(si,1)  = nout(si).classic.bone_std;
              out.TIV{gi}(si,1)       = nout(si).tismri.TIV; 
              out.int_ci{gi}(si,:)    = nout(si).tismri.Tth; 
              out.vol_ci{gi}(si,:)    = nout(si).tismri.volr;    % relative volume by TIV 
              out.con{gi}(si,1)       = nout(si).tis.seg8conr;     % GM/WM contrast
              out.tis_bonecortex{gi}(si,1)  = nout(si).tis.bonecortex; 
              out.tis_bonemarrow{gi}(si,1)  = nout(si).tis.bonemarrow; 
              out.tis_bonedensity{gi}(si,1) = nout(si).tis.bonedensity; 
              out.tis_bone{gi}(si,1)        = nout(si).tis.bone; 
              out.tis_head{gi}(si,1)        = nout(si).tis.head; 
              out.tismri_volfatr{gi}(si,1)  = nout(si).tismri.volfatr; 
              out.tismri_volmusr{gi}(si,1)  = nout(si).tismri.volmusr; 
              
              % vol
              out.vROI_bonecortex1{gi}(si,1)    = nout(si).vROI.bonecortex(1);
              out.vROI_bonecortex2{gi}(si,1)    = nout(si).vROI.bonecortex(2);
              out.vROI_bonecortex3{gi}(si,1)    = nout(si).vROI.bonecortex(3);
              out.vROI_bonecortex4{gi}(si,1)    = nout(si).vROI.bonecortex(4) + nout(si).vROI.bonecortex(5);

              out.vROI_bonemarrow1{gi}(si,1)    = nout(si).vROI.bonemarrow(1);
              out.vROI_bonemarrow2{gi}(si,1)    = nout(si).vROI.bonemarrow(2);
              out.vROI_bonemarrow3{gi}(si,1)    = nout(si).vROI.bonemarrow(3);
              out.vROI_bonemarrow4{gi}(si,1)    = nout(si).vROI.bonemarrow(4) + nout(si).vROI.bonemarrow(5);
              out.vROI_bonemarrow0{gi}(si,1)    = nout(si).vROI.bonemarrow(1) + nout(si).vROI.bonemarrow(2) ...
                                               -2*nout(si).vROI.bonemarrow(2) - nout(si).vROI.bonethickness(1) ...
                                               +  nout(si).vROI.bonethickness(2) ;

              out.vROI_bonethickness1{gi}(si,1) = nout(si).vROI.bonethickness(1);
              out.vROI_bonethickness2{gi}(si,1) = nout(si).vROI.bonethickness(2);
              out.vROI_bonethickness3{gi}(si,1) = nout(si).vROI.bonethickness(3);
              out.vROI_bonethickness4{gi}(si,1) = nout(si).vROI.bonethickness(4) + nout(si).vROI.bonethickness(5);

              out.vROI_headthickness2{gi}(si,1) = nout(si).vROI.headthickness(2);
              out.vROI_headthickness3{gi}(si,1) = nout(si).vROI.headthickness(3);
              out.vROI_headthickness4{gi}(si,1) = nout(si).vROI.headthickness(4) + nout(si).vROI.headthickness(5);
              %out.vROI_head3{gi}(si,1)          = nout(si).vROI.head(3);
              %out.vROI_head4{gi}(si,1)          = nout(si).vROI.head(4);
              % surf
              if opt.matonly == 12
                out.sROI_bonecortex1{gi}(si,1)    = nout(si).sROI.bonecortex(1);
                out.sROI_bonecortex2{gi}(si,1)    = nout(si).sROI.bonecortex(2);
                out.sROI_bonecortex3{gi}(si,1)    = nout(si).sROI.bonecortex(3);
                out.sROI_bonecortex4{gi}(si,1)    = nout(si).sROI.bonecortex(4) + nout(si).sROI.bonecortex(5);
                out.sROI_bonemarrow1{gi}(si,1)    = nout(si).sROI.bonemarrow(1);
                out.sROI_bonemarrow2{gi}(si,1)    = nout(si).sROI.bonemarrow(2);
                out.sROI_bonemarrow3{gi}(si,1)    = nout(si).sROI.bonemarrow(3);
                out.sROI_bonemarrow4{gi}(si,1)    = nout(si).sROI.bonemarrow(4) + nout(si).sROI.bonecortex(6);
                out.sROI_bonethickness1{gi}(si,1) = nout(si).sROI.bonethickness(1);
                out.sROI_bonethickness2{gi}(si,1) = nout(si).sROI.bonethickness(2);
                out.sROI_bonethickness3{gi}(si,1) = nout(si).sROI.bonethickness(3);
                out.sROI_bonethickness4{gi}(si,1) = nout(si).sROI.bonethickness(4) + nout(si).sROI.bonethickness(5);
                out.sROI_headthickness1{gi}(si,1) = nout(si).sROI.headthickness(1);
                out.sROI_headthickness2{gi}(si,1) = nout(si).sROI.headthickness(2);
                out.sROI_headthickness3{gi}(si,1) = nout(si).sROI.headthickness(3);
                out.sROI_headthickness4{gi}(si,1) = nout(si).sROI.headthickness(4) + nout(si).sROI.headthickness(5);
              end

              % tried to normalize this but it is not working
              %TIVnorm = ((nout(si).tismri.TIV / (pi*4/3)).^(1/3)) / ((1350 / (pi*4/3)).^(1/3)); 
              TIVnorm = (nout(si).tismri.TIV / 1600 ) .^ (1/3); 
              out.vROI_headthickness3t{gi}(si,1) = nout(si).vROI.headthickness(3) / TIVnorm;
              out.sROI_headthickness3t{gi}(si,1) = nout(si).sROI.headthickness(3) / TIVnorm;

              % combined head measures - not working
              out.vhdt1{gi}(si,1) =  nout(si).tis.head + nout(si).tismri.volfatr;
              out.shdt3{gi}(si,1) = (nout(si).tis.head + nout(si).tismri.volfatr) .* nout(si).sROI.headthickness(3); 

              % combined bone measures - working!
              out.vROI_BMDH{gi}(si,1)    = -nout(si).vROI.bonecortex(3) + nout(si).vROI.bonethickness(1) ...
                                           -nout(si).vROI.bonecortex(3) + nout(si).vROI.bonethickness(3);
              out.vROI_BMDH2{gi}(si,1)   = -nout(si).vROI.bonecortex(3) + nout(si).vROI.bonethickness(1) ...
                                           -nout(si).vROI.bonecortex(3) + nout(si).vROI.bonethickness(1);
              out.vROI_BMDH4{gi}(si,1)   = -nout(si).vROI.bonecortex(3) + nout(si).vROI.bonethickness(1) ...
                                           -nout(si).vROI.bonemarrow(3);
              out.vROI_BMDH1{gi}(si,1)   = -nout(si).vROI.bonecortex(1) + nout(si).vROI.bonethickness(3);
              out.vROI_BMDH3{gi}(si,1)   = -nout(si).vROI.bonecortex(3) + nout(si).vROI.bonethickness(1);
              if opt.matonly == 12
                out.sROI_BMDH{gi}(si,1)  = -nout(si).sROI.bonecortex(3) + nout(si).sROI.bonethickness(1) ...
                                           -nout(si).classic.bone_med   + nout(si).sROI.bonethickness(3);
                out.sROI_BMDH2{gi}(si,1) = -nout(si).sROI.bonecortex(3) + nout(si).sROI.bonethickness(3) ...
                                           -nout(si).classic.bone_med   + nout(si).sROI.bonemarrow(3);
              end       
            end
          end
        else

          if opt.matonly > 1
          % new approach
  
            opt2.files      = groups{gi}; 
            opt2.method     = min(2,opt.matonly - 1); 
            opt2.report     = 3 * opt.printimg; 
            opt2.writevol   = opt.printvol;
            opt2.writesurf  = opt.printvol; 
            opt2.subdirs    = 1; 
            opt2.affreg     = 0; % 0-
            opt2.reduce     = 2; %
            opt2.rerun      = opt.rerun; 
            opt2.Patlas     = fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-regions.nii');
            opt2.Pmask      = fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-mask.nii');
            %opt2.Patlas     = fullfile(spm('dir'),'tpm','labels_Neuromorphometrics_bones.nii');
            opt2.verb       = 1; % 1-default, 2-details
  
            nout = cat_run_boneseg(opt2);
  
            if opt.matonly == 4 % classic SPM values for all processings
              for si = 1:numel(groups{gi})
                out.bone_num{gi}(si,1)  = nan;                    % bone volume normalised by TIV
                out.bone_med{gi}(si,1)  = nout(si).fst.spm_bone_med1; % median bone intensity value
                out.bone_std{gi}(si,1)  = nan; 
                out.TIV{gi}(si,1)       = nan; 
                out.int_ci{gi}(si,:)    = nan(1,6); 
                out.vol_ci{gi}(si,:)    = nan(1,6); 
                out.con{gi}(si,1)       = nan;
              end
            elseif opt.matonly == 5 % classic MED values (simple corrections)
              for si = 1:numel(groups{gi})
                out.bone_num{gi}(si,1)  = nout(si).fst.bone_num; % bone volume normalised by TIV
                out.bone_med{gi}(si,1)  = nout(si).fst.bone_med; % median bone intensity value
                out.bone_std{gi}(si,1)  = nout(si).fst.bone_std; 
                out.TIV{gi}(si,1)       = nout(si).fst.TIV; 
                out.int_ci{gi}(si,:)    = nout(si).fst.int_ci;
                out.vol_ci{gi}(si,:)    = nout(si).fst.vol_ci; 
                out.con{gi}(si,1)       = nout(si).fst.con;
              end    
            else
              vFN = setdiff(fieldnames(nout(1).vROI),{'help','boneatlas_id','boneatlas_name'});
              sFN = setdiff(fieldnames(nout(1).sROI),{'help','boneatlas_id','boneatlas_name'});
                
              for si = 1:numel(groups{gi})
                % extract volume ROI measures
                for vFNi = 1:numel(vFN)
                  for vFNij = 1:numel(nout(si).vROI.(vFN{vFNi}))
                    out.(vFN{vFNi}){gi}(si,vFNij) = nout(si).vROI.(vFN{vFNi})(vFNij);
                  end
                end
                % extract surface ROI measures
                for vFNi = 1:numel(sFN)
                  for vFNij = 1:numel(nout(si).sROI.(sFN{vFNi}))
                    out.(sFN{vFNi}){gi}(si,vFNij) = nout(si).sROI.(sFN{vFNi})(vFNij);
                  end
                end
  
                % other 
                out.bone_med{gi}(si,1)  = nout(si).fst.bone_med; % median bone intensity value
                out.bone_std{gi}(si,1)  = nout(si).fst.bone_std; 
                out.bone_num{gi}(si,1)  = nout(si).fst.bone_num; % bone volume normalised by TIV
                out.bone_vol{gi}(si,1)  = nout(si).tismri.volr(4); % bone volume normalised by TIV
                out.TIV{gi}(si,1)       = nout(si).tismri.TIV;
                out.int_ci{gi}(si,:)    = nout(si).tismri.Tth;       
                out.vol_ci{gi}(si,:)    = nout(si).tismri.volr;    % relative volume by TIV 
                out.con{gi}(si,1)       = nout(si).tis.seg8conr;     % GM/WM contrast
                
                if ~isnan(nout(si).tismri.Tth(1) )
    % this part is still under construction and map processed parameters to 
    % the variables used in our analysis script
                  out.mnBone{gi}(si,1)   = nout(si).tismri.iBonemn2(1);   % [.72 |.84*][.70 |.76*][.24|.22]
                  out.mnMarrow{gi}(si,1) = nout(si).tismri.iBonemn2(2);   % [.72 |.84*][.70 |.76*][.24|.22]
                end
                
                % bone/head thickness 
                out.bdt{gi}(si,1)      = nout(si).tismri.tBone;     
                out.bdt31{gi}(si,1)    = nout(si).tismri.tBonemn(1);
                out.bdt32{gi}(si,1)    = nout(si).tismri.tBonemn(2);
                out.bdt33{gi}(si,1)    = nout(si).tismri.tBonemn(3);
                out.hdt{gi}(si,1)      = nout(si).tismri.headthickmn;     
  
                % bone+head thickness               
                out.bvht{gi}(si,1)       = nout(si).tismri.tBone / (nout(si).tismri.tBone + nout(si).tismri.headthickmn);     
                out.b2ht{gi}(si,1)       = nout(si).tismri.tBone \ nout(si).tismri.headthickmn;     
                out.b2htx{gi}(si,1)      = nout(si).tismri.tBone \ nout(si).tismri.headthickmn * log(out.TIV{gi}(si,1)); %^(1/3);     
                out.bht{gi}(si,1)        = nout(si).tismri.tBone + nout(si).tismri.headthickmn;     
                out.bhttiv{gi}(si,1)     = nout(si).tismri.tBone + nout(si).tismri.headthickmn * out.TIV{gi}(si,1)^(1/3);     
                out.bdtx{gi}(si,1)       = nout(si).tismri.tBonemn(1) + nout(si).tismri.tBonemn(2) - + nout(si).tismri.tBonemn(3);     
                out.bhttivx{gi}(si,1)    = nout(si).tismri.tBonemn(1) + nout(si).tismri.tBonemn(2) - + nout(si).tismri.tBonemn(3) ...
                                           + nout(si).tismri.headthickmn * out.TIV{gi}(si,1)^(1/3);     
                
  
                %% global threhholds
                out.bonessmn{gi}(si,:) = nout(si).tismri.surfmn;  
                out.bonessmd{gi}(si,:) = nout(si).tismri.surfmd;  
                out.bonesssd{gi}(si,:) = nout(si).tismri.surfsd;  
                out.bonessiq{gi}(si,:) = nout(si).tismri.surfiqr;  
                out.bonesvmn{gi}(si,:) = nout(si).tismri.volmn;  
                out.bonesvmd{gi}(si,:) = nout(si).tismri.volmd;  
                out.bonesvsd{gi}(si,:) = nout(si).tismri.volsd;  
                out.bonesviq{gi}(si,:) = nout(si).tismri.voliqr;  
                
                % surface/volume based kmeans measures ... 
                % ... but kmeans worse than global MN or MD
                out.boneskmeans21{gi}(si,1)   = nout(si).tismri.iBonemn2(1);  
                out.boneskmeans31{gi}(si,1)   = nout(si).tismri.iBonemn3(1); 
                out.boneskmeans22{gi}(si,1)   = nout(si).tismri.iBonemn2(2);  
                out.boneskmeans32{gi}(si,1)   = nout(si).tismri.iBonemn3(2);   % winner surf 
                out.boneskmeans33{gi}(si,1)   = nout(si).tismri.iBonemn3(3);   
                out.bonevkmeans21{gi}(si,1)   = nout(si).tismri.volmn2(1);  
                out.bonevkmeans31{gi}(si,1)   = nout(si).tismri.volmn3(1);     % winner vol
                out.bonevkmeans22{gi}(si,1)   = nout(si).tismri.volmn2(2);  
                out.bonevkmeans32{gi}(si,1)   = nout(si).tismri.volmn3(2); 
                out.bonevkmeans33{gi}(si,1)   = nout(si).tismri.volmn3(3); 
  % > use proposion value to replace intensities ################           
  
                out.bonehdt1{gi}(si,1) = -nout(si).tismri.iBonemn2(1) + out.hdt{gi}(si,1) ...
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone; 
                out.bonehdt2{gi}(si,1) = -nout(si).tismri.iBonemn2(2) + out.hdt{gi}(si,1) ...
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone;
                out.bonehdt3{gi}(si,1) = -mean(nout(si).tismri.iBonemn2) + ...
                                         -out.bone_med{gi}(si,1) + out.hdt{gi}(si,1) * out.TIV{gi}(si,1)^(1/3) ...
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone;
  
                out.bonehdt4{gi}(si,1) =                           out.hdt{gi}(si,1) * out.TIV{gi}(si,1)^(1/3) ...
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone;
                out.bonehdt5{gi}(si,1) = -out.bone_med{gi}(si,1) + out.hdt{gi}(si,1)  ... 
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone; 
                out.bonehdt6{gi}(si,1) = -out.bone_med{gi}(si,1) + out.hdt{gi}(si,1) ... 
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone; % BEST SO FAR
                out.bonehdt7{gi}(si,1) = -out.bone_med{gi}(si,1) + out.hdt{gi}(si,1) * out.TIV{gi}(si,1)^(1/3) ... 
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone * out.TIV{gi}(si,1)^(1/3); % BEST SO FAR
                out.bonehdt8{gi}(si,1) = -nout(si).fst.bone_med  + out.hdt{gi}(si,1) * out.TIV{gi}(si,1)^(1/3) ... 
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone; % BEST SO FAR
                out.bonehdt9{gi}(si,1) = -nout(si).tismri.surfmd - nout(si).tismri.volmd  + out.hdt{gi}(si,1) * out.TIV{gi}(si,1)^(1/3) ... 
                                         -nout(si).fst.bone_med  + nout(si).tismri.tBone; 
                out.bonehdt10{gi}(si,1) = -nout(si).fst.bone_med  + out.hdt{gi}(si,1) * log(out.TIV{gi}(si,1)) ... 
                                          -nout(si).fst.bone_med  + nout(si).tismri.tBone; % BEST SO FAR
                out.bonehdt11{gi}(si,1) = -nout(si).fst.bone_med  + out.hdt{gi}(si,1) * log(out.TIV{gi}(si,1)) ... 
                                          -nout(si).fst.bone_med  + nout(si).tismri.tBone * log(out.TIV{gi}(si,1)); 
                out.bonehdt12{gi}(si,1) = -log(nout(si).fst.bone_med)  + log(out.hdt{gi}(si,1) / (out.TIV{gi}(si,1))) ... 
                                          -log(nout(si).fst.bone_med)  + log(nout(si).tismri.tBone / (out.TIV{gi}(si,1)));
                out.bonehdt13{gi}(si,1) = ( -log(nout(si).fst.bone_med) -nout(si).tismri.surfmd - nout(si).tismri.volmd) .* ( (out.hdt{gi}(si,1)) + (nout(si).tismri.tBone )); 
                  %           
              end
            end
          else
            % old approach
            for si = 1:numel(groups{gi})
           
              fprintf('  Subject %d/%d\n',si,numel(groups{gi})); 
  
              %% load data
              Vm = spm_vol(groups{gi}{si}); 
              Ym = spm_read_vols(Vm); 
              
              % get SPM TPMs for GM, WM, CSF, and skull   
              [pp,ff,ee] = spm_fileparts(groups{gi}{si});
              tmpmat = load(fullfile(pp,sprintf('%s%s',ff(2:end),'_seg8.mat'))); 
              out.mat{gi}(si) = rmfield(tmpmat,{'Twarp','Tbias'});
   
              if 0 %max(out.mat{gi}(si).lkp) > 6
                  error('class_error',sprintf('%d-TPM not preparted yet.', max(out.mat{gi}(si).lkp))); 
              end
              
              % get high prob value of CSF tissue
              out.spm_con{gi}(si) = sum(out.mat{gi}(si).mn( out.mat{gi}(si).lkp == 3 )' .* ...
                                     ( out.mat{gi}(si).mg( out.mat{gi}(si).lkp == 3 ) ==  ...
                                       max(out.mat{gi}(si).mg( out.mat{gi}(si).lkp == 3 ) ) ));
              if max(out.mat{gi}(si).lkp) > 6
                  out.spm_bone_med1{gi}(si) = out.mat{gi}(si).mn( out.mat{gi}(si).lkp == 4 ) / out.spm_con{gi}(si);
                  out.spm_bone_med2{gi}(si) = out.mat{gi}(si).mn( out.mat{gi}(si).lkp == 5 ) / out.spm_con{gi}(si);
              else
                  % default case for 6 TPM classes
                  bone = sort(out.mat{gi}(si).mn( out.mat{gi}(si).lkp == 4 ) / out.spm_con{gi}(si)); % order intensities within class 4 (normalized by CSF)
                  out.spm_bone_med1{gi}(si) = bone(3); % highest value
                  out.spm_bone_med2{gi}(si) = bone(1); 
                  out.spm_bone_med3{gi}(si) = bone(2); 
              end
             
              if opt.matonly
                  out.bone_num{gi}(si) = nan;                      % bone volume normalised by TIV
                  out.bone_med{gi}(si) = out.spm_bone_med1{gi}(si); % median bone intensity value
                  out.bone_std{gi}(si) = nan; 
                  out.TIV{gi}(si,1)    = nan; 
                  out.int_ci{gi}(si,:) = nan(1,6);
                  out.vol_ci{gi}(si,:) = nan(1,6); 
                  out.con{gi}(si)      = nan;
              else
                  outnamenii = fullfile(pp,sprintf('MED_%s%s',ff,ee)); 
                  outnamemat = fullfile(pp,sprintf('MED_%s.mat',ff)); 
                  
                  clear S; 
                  if ~exist(outnamemat,'file') || opt.rerun
                      % load SPM tissue segments
                      for ci = 1:5 
                         Vp{ci} = spm_vol(fullfile(pp,sprintf('c%d%s%s',ci,ff(2:end),ee))); 
                         Yp{ci} = spm_read_vols(Vp{ci}); 
  
                         S.int_ci(1,ci) = median(Ym(Yp{ci}>0.9)); %intensity
                         S.vol_ci(1,ci) = sum(Yp{ci}(:)>0.5) / 1000; %volume
                      end
                      % estimate TIV
                      S.TIV = sum(S.vol_ci(1,1:3));
  
                      % create bone mask
                      if max(out.mat{gi}(si).lkp) > 6
                          Ymsk = smooth3(Yp{4} + Yp{5})>0.5;     % remove small noise
                      else
                          Ymsk = smooth3(Yp{4})>0.25;             % remove small noise
                      end
                      Ymsk = cat_vol_morph(Ymsk,'l');         % get the largest bones element
                      Ymsk = cat_vol_morph(Ymsk,'dc',6);      % remove holes in the bone
                  
                      % extract values (hist)
                      S.con       = S.int_ci(1,3); % CSF contrast
                      S.bone_num  = (sum(Ymsk(:)) / 1000) / S.TIV; % bone volume normalised by TIV
                      S.bone_med  = median(Ym(Ymsk(:))) / S.con; % median bone intensity value
                      S.bone_std  = std(Ym(Ymsk(:)))    / S.con;                     
                      
                      % write masked image
                      if opt.printvol
                          Vout = Vm; Vout.fname = outnamenii;
                          spm_write_vol(Vout,( Ym / S.con ) .* Ymsk);
                      end
  
                      save(outnamemat,'S'); 
                  else
                      load(outnamemat,'S'); 
                  end
                  
                  out.bone_num{gi}(si) = S.bone_num; % bone volume normalised by TIV
                  out.bone_med{gi}(si) = S.bone_med; % median bone intensity value
                  out.bone_std{gi}(si) = S.bone_std; 
                  out.TIV{gi}(si,1)    = S.TIV; 
                  out.int_ci{gi}(si,:) = S.int_ci;
                  out.vol_ci{gi}(si,:) = S.vol_ci; 
                  out.con{gi}(si)      = S.con;
                 
              end
            end
          end
        
        end
        out.rGMV{gi} =  out.vol_ci{gi}(:,1); % ./ out.TIV{gi};
        out.rWMV{gi} =  out.vol_ci{gi}(:,2); % ./ out.TIV{gi};
        out.rCMV{gi} =  out.vol_ci{gi}(:,3); % ./ out.TIV{gi};
        
    end
   
    %%
    FN = setdiff( fieldnames(out) , {'int_ci','vol_ci'} ) ; 
    for fni = 1:numel(FN)
        out.(FN{fni}) = reshape( out.(FN{fni}) , size(groups) ); 
    end
end
