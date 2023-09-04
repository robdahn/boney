function [Yc,Ye,Ya,clscor] = boney_segment_refineSPM(Yo,Ym,Yc,Ya,Ybraindist0,tis,tismri)
% == REFINE SPM DATA ==
% * Create meninges class? 
%   - This is not realy working cause the intensities are changing too much.
%     Hence, it is better to evaluate the whole mid structure. 
%   - Would thresholding by thickness work? 
%     Expected, but it is better to use the thickness as correction factor.
% * What about the blood vessels?
%   - Local filtering of outliers might help. 
%     Hence, we use a median filter of bone volume.
% * Split-up skull class by bone and bone marrow?
%   - This was not really working as the contrast changes too much over aging
%     and is not robust enough for both structures in various protocols, with 
%     shift artifacts of, e.g., fat. 
% * Analyse head by muscles and fat (and skin or other things?)
%   - Not yet.
% * New affine registration for bone only?
%   - should be similar in adults?


% #######################
% Test for errors in SPM segmentation and suggest strategies.  
% * E.g. when the cls4 is also in the real background 
% Strategies: 
% * run SPM with higher separation (sep)
%
% * create warning in case of too small BB or strong or strange defacing? 

% #######################
% * get FAT peak from SPM to quantify if fat suppression was used
% * get FAT peak from closes head regions 
% >>>  
% * use FAT peak for general scaling 

Yco    = Yc; 
vx_vol = tis.res_vx_vol; 
vxmm3  = prod(vx_vol) * 1000; 
Ybraindist0s = cat_vol_smooth3X(Ybraindist0,6);  
clscor.help = 'Boney refinement of SPM tissue classes described by the transferred tissue volume transfer between classes.';

  %%
  Yc = Yco; 
  
  % defacing >> background 
  Ydeface = cat_vol_morph(Yo==0,'l',[10,0.1])>0; % need multiple objects for face and ears 
  if sum(Ydeface(:)) > 10000
    for ci = 1:5, Ycn = Yc{ci} .* Ydeface; Yc{6} = Yc{6} + Ycn; Yc{ci} = Yc{ci} - Ycn; end
  end

  % background >> bone | head
  % background should be only outside the head ... not optimal for air
  Yc6    = Yc{6} .* cat_vol_morph(cat_vol_morph(Yc{6}<.5,'lc'),'e'); 
  clscor.BG2BN = sum(Yc6(:) .* ~(Yo(:)<=tismri.Tth(3) & (Ybraindist0s(:)<10))) / vxmm3 / tismri.TIV; 
  clscor.BG2HD = sum(Yc6(:) .*  (Yo(:)<=tismri.Tth(3) & (Ybraindist0s(:)<10))) / vxmm3 / tismri.TIV; 
  Yc{6}  = Yc{6} - Yc6;
  Yc{5}  = Yc{5} + Yc6 .* ~(Yo<=tismri.Tth(3) & (Ybraindist0s<10)); 
  Yc{4}  = Yc{4} + Yc6 .*  (Yo<=tismri.Tth(3) & (Ybraindist0s<10));
  % head >> background 
  % In children the head/bone class is too large and needs correction. 
  if 1 % seg8t.BGtype == ... % low intensity background
    Yc6a = Yc{5} .* smooth3(cat_vol_morph( ((Yc{6} + Yc{5})>.5 & ~Ydeface & Yo<tismri.Tth(3)), 'lo'));
  else
    % noisy background of MP2RAGE/MT sequences need some other (gradient-based) definition 
  end
  clscor.HD2BN = sum(Yc6a(:)) / vxmm3  / tismri.TIV; 
  Yc{5}  = Yc{5} - Yc6a;
  Yc{6}  = Yc{6} + Yc6a; 
  % final background closing 
  Yc6    = cat_vol_morph(cat_vol_morph(Yc{6}<.5,'lc')<.5,'e'); 
  for ci = 1:5, Ycn = Yc{ci} .* Yc6; Yc{6} = Yc{6} + Ycn; Yc{ci} = Yc{ci} - Ycn; end


  %% head (+ background) ~ bonemarrow>> bone
  if tis.headBoneType
    Yhead  = min(1,single(Yc{5} + Yc{6} + 0.1 * max(0,Ybraindist0s-15))); % no smoothing here!
    Yhead  = smooth3(cat_vol_morph(cat_vol_morph(smooth3(Yhead)>.7,'lo',4),'d')) .* Yhead .* (Ybraindist0s>0);
    cn = 0; for ci = 5:6, Ycn = Yc{ci} .* Yhead; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
    Ye{1} = zeros(size(Yc{1}),'single'); 
  else
    dmn = cat_stat_kmeans( Ybraindist0s(Yc{4}>.5 & ~cat_vol_morph(Yc{4}>.5,'e') & Ybraindist0s<(median(Ybraindist0s(Yc{4}>.5))*4) ) , 2 );
    Ybone = cat_vol_morph(cat_vol_morph(Yc{4}>.5 & Ybraindist0s>(dmn(1)*.5) & Ybraindist0s<(dmn(2)*1.2),'lc'),'o',1);
    cn = 0; for ci = [1:3,5:6], Ycn = Yc{ci} .* ~Ybone; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
   
    %% get the veins
    Yhead  = min(1,single(Yc{5} + Yc{6} + 0.1 * max(0,Ybraindist0s-15))); % no smoothing here!
    Yhead2 = smooth3(cat_vol_morph(cat_vol_morph(smooth3(Yhead)>.7,'lo',4),'d')) .* Yhead .* (Ybraindist0s>0);
    Yhead  = 1 - (Yhead - Yhead2); 
    % Yhead  = cat_vol_morph(Yhead + (Yc{4}>.5),'ldc',3); 
    % Yhead  = cat_vol_smooth3X(Yhead,4)>.5; 
    % Yhead  = Yhead .* cat_vol_morph(cat_vol_morph((Yc{1}+Yc{2}+Yc{3}+Yc{4})>.5,'lc',2),'do',7); 
    %%
    cn = 0; Ye{1} = zeros(size(Yc{1}),'single'); for ci = 5:6, Ycn = Yc{ci} .* Yhead; cn = cn + sum(Ycn(:)); Ye{1} = Ye{1} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
  end

  %% brain > bone
  Ybrain = single(Yc{1} + Yc{2} + Yc{3}); 
  Ybrain = (cat_vol_morph( Ybrain > .5,'ldo',3)); % & (Ybraindist0>10)) | (Ybrain & (Ybraindist0<10));           % remove skull
  Ybrain = smooth3(cat_vol_morph( Ybrain > 0.5,'dc',3));  % this closing includes meninges!
  Ybrain = max( Ybrain , cat_vol_smooth3X(Ybrain,2)>.6); % remove vein
  cn = 0; for ci = 1:3, Ycn = Yc{ci} .* Ybrain; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
  clscor.BR2BN = cn / vxmm3 / tismri.TIV; 
  
  % Yc4 outside the TPM definition >> head | background
  Yc4    = cat_vol_morph(Yc{4}>0,'l') & (Ybraindist0s>10); 
  clscor.BN2HD = sum(Yc4(:)) / vxmm3 / tismri.TIV; 
  Yc{5}  = Yc{5} + Yc{4} .* Yc4; 
  Yc{4}  = Yc{4} - Yc{4} .* Yc4; 

% These corrections are not enough in children with inoptimal TPM overlay (e.g. the NIH templates). 
% In the head worst-case (but good brain case) the bone may fully include the head. 
% Separation maybe by region-growing of (high) head and (low) bone regions. 
% 

  %% head > bone
  %{
  Ybrain = single(Yc{1} + Yc{2} + Yc{3}); 
  Ybrain = (cat_vol_morph( Ybrain > .5,'lo',2) & (Ybraindist0>10)) | (Ybrain & (Ybraindist0<10));           % remove skull
  Ybrain = smooth3(cat_vol_morph( Ybrain > 0.5,'c',3));  % this closing includes meninges!
  Ybrain = max( Ybrain , cat_vol_smooth3X(Ybrain,2)>.6); % remove vein
  for ci = 1:3, Ycn = Yc{ci} .* Ybrain; Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
  

  
  %% head >> BG 
  Yc5lth = min(.5,tismri.Tth(5)/tismri.Tth(2));
  Yc5    = Yc{5} .* (Ym<min(0.2,Yc5lth)) .* cat_vol_morph(Yc{5}+Yc{6} > .5,'e'); 
  Yc{6}  = Yc{6} + Yc5; 
  Yc{5}  = Yc{5} - Yc5; 
  %Yc6    = Yc{5} .* (Ym>min(0.2,Yc5lth)) .* cat_vol_morph(Yc{5}+Yc{6},'e'); 
  
  clear Yc5 Yc5lth; 
  %}
  
  %% close bone (by getting from all classes)
  %{
  Yc4    = 1 - (cat_vol_morph(Yc{4}>.25 & ~Ybrain,'lc',5) - Yc{4});
  for ci = [1:3 5:6], Ycn = Yc{ci} .* Yc4; Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
  %}
% TPM layber model?
% - distance from WM and from BG
%%
  if 0
    %% meninges detection in low int ? 
      Yn  = cat_vol_localstat(single(Ym),(Yc{4}+Yc{3})>.5 & ~cat_vol_morph(Yc{1}>.5,'d'),1,4);
      Yn2 = cat_vol_localstat(single(Ym),(Yc{4}+Yc{3})>.5 & ~cat_vol_morph(Yc{1}>.5,'d'),2,4);
      Ym .* (Yc{4} + Yc{3}),Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}) .* Ym .* 3; 
      Ym .* (1- smooth3(Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}).*Ym * 3))
      Ym .* (1- smooth3(Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}).*Ym * 3)) % ... this is quite nice to remove mengs and BVs (but also affects the GM/WM boundary !
  end


  %% refinement of atlas map
  if 0 
    Ymg = cat_vol_div( Ym .* cat_vol_smooth3X(Yc{4}>.5,1));
    Ymg = cat_vol_smooth3X( Ymg ./ mean(Ymg(Ymg(:)<-0.1)) , 4) ;
    Yag = cat_vol_grad(Ya);
    Ya2 = Ya .* (cat_vbdist(Yag,Yhead<1,vx_vol)>20 & Yc{4}>.5); 
    Ya2(Ya2==0 & Yc{4}<.1) = nan; 
    %
    Yx = cat_vol_smooth3X(Yc{4},2) .* Ymg; 
    Ya3 = cat_vol_downcut(Ya2,Yx,-.05);
    Ya3(Ya3==0 & Yc{4}<.1) = nan; 
    Ya3 = cat_vol_downcut(Ya3,Yx,.2);
    [Yd,Yi] = cat_vbdist(Ya3); Ya3f = Ya3(Yi);
    Ya3g = cat_vol_grad(Ya3f);
  end

  % ########
  % The problem is to get all of the bone marrow, also the misaligned things that are
  % now part of the head or the brain (aligned as meninges-like thing to
  % CSF/GM/WM). 
  %% magic operation that is now working
  test = 0;
  if test
    Ybg = smooth3(Yc{6})>.5 | Ybraindist0s>15 ; 
    Ybe = Ybrain>.5; % & smooth3(Yc{1} + Yc{2} + (Yc{3} .* (Ym<.5)))>.8; % ############## T1! 
    Ybe = cat_vol_morph(Ybe,'lc');
    %Yhd = single( Ybg ); Yhd( Ybe ) = nan; [~,Yhd] = cat_vol_downcut( Yhd , max(eps,min(10,  3 - Ym*3)), .5);
    %Yhd = single( Ybg ); Yhd( Ybe ) = nan; [~,Yhd] = cat_vol_downcut( Yhd , max(eps,min(10,  Ym*3)), .5);
    Ybd = single( Ybe ); Ybd( Ybg ) = nan; [~,Ybd] = cat_vol_downcut( Ybd , min(10,  Ym*3) , .5);
    %Yboneh = ( (min(Yhd,Ybd*1e10) ./ (Yhd+Ybd))>.5);
    Yboneh = Ybd./max(1,Ybraindist0s - median(Ybraindist0s(Yc{4}(:)>.5))) / 100; 
    Yboneh(Ybg | isnan(Yboneh) | Yboneh>2) = 0;  
    Ych = zeros(Yboneh);
    for ci = [1:3 5:6], Ych = Ych + Yc{5} .* Yboneh; end
    Yc5   = Yc{5} .* Yboneh;
    Yc{4} = Yc{4} + Yc5;
    Yc{5} = Yc{5} - Yc5; 
    clear Yc5; 
  end

  
end
