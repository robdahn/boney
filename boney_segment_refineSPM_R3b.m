function [Yc,Ye,Ya,clscor] = boney_segment_refineSPM_R3(Yo,Ym,Yc,Ya,Ybraindist0,tis,tismri,refine)
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
  Yc    = Yco; 
  Ybox  = true(size(Ym)); Ybox(2:end-1,2:end-1,2:end-1) = false;  
  Ybox4 = true(size(Ym)); Ybox4(5:end-4,5:end-4,5:end-4) = false;  
  sx    = min(1,1/mean(vx_vol)); 
   
  %% RD202402
  if refine > 1
  % new refinement block 
    Yc = Yco;
    
    % 1) Rought correction of usual low-intensity backgrounds 
    %    Important if SPM describes noisy backgounds, e.g., 
    %    BuchertHTP_c056p065r00_scan495_SiemensAvanto1.5T_t191s.nii
    if tis.highBG == 0
      % use lower for speedup but 3 or 4 mm can have issues with to close image boundaries
      [Ymr,Yc6r,RES] = cat_vol_resize({Ym,Yc{end}} ,'reduceV'  ,tis.res_vx_vol,2,16,'meanm');
      Ybgr = (cat_vol_smooth3X(Yc6r,sx) > 0.5 & (Ymr<.3)) | cat_vol_morph( Ymr < .2,'ldo',2);
      Ybgr = cat_vol_morph( Ybgr , 'l'  , [10 0.1]); 
      Ybgr = cat_vol_morph( Ybgr , 'dc' , 2); 
      Yboxr = true(size(Ymr)); Yboxr(2:end-1,2:end-1,2:end-1) = false; Ybgr(Yboxr) = max(0,1 - Ymr(Yboxr)); clear Yboxr;
      Ybg   = cat_vol_resize( cat_vol_smooth3X( Ybgr , sx), 'dereduceV', RES);
      % reduce image background boundary effects
      Yboxx = cat_vol_smooth3X( cat_vol_morph( Ym<0.8 & Ybox4 ,'lo',1,vx_vol),2)*.6; 
      Ybg  = max(Yboxx,Ybg); %clear Yboxx
      Ybg  = min(1,max(0, 4*tan( Ybg - 0.5) + .5)); % higher value = sharper
      clear Yc6r Ybgr RES; 
      for ci = 1:numel(Yc) - 1
        Yc2b    = Yc{ci}  .* Ybg; 
        Yc{end} = Yc{end} + Yc2b; 
        Yc{ci}  = Yc{ci}  - Yc2b; 
      end
      Yb2c      = Yc{end} .* (1-Ybg);
      Yc{end-1} = Yc{end-1} + Yb2c;
      Yc{end}   = Yc{end}   - Yb2c; 
      clear Yc2b Yb2c;
    else
      Ybg = cat_vol_morph( Yc{end} > .5 , 'lo', 1,vx_vol); 
    end


    % 2) Rought correct skull-stripping
    %    Important if SPM head is in the brain, e.g., 
    %    BuchertHTP_c062p074r00_scan076_SiemensEspree1.5T_t191s.nii
    %    BuchertHTP_c093p115r00_scan316_SiemensSkyra3T_t192s
    Yb = cat_vol_smooth3X(cat_vol_morph( (Yc{1} + Yc{2} + Yc{3}) > 0.5 , 'ldo' , 5), sx); 
    if 1
    % brain tissue threshold may cause side effects (uncovered in a fast low-res testcase)  
      Ycsf = Yb .* Yc{3};  
    elseif tis.weighting == 2 % T2
      Ycsf = Yb .* cat_vol_smooth3X(Yc{3}>.5 & Ym>0.8*max(tis.seg8n(1:2))*1.2, sx); 
    else  
      Ycsf = Yb .* cat_vol_smooth3X(Yc{3}>.5 & Ym<1.2*mean([ min(tis.seg8n(1:2)), tis.seg8n(3)]),sx); 
    end
    [Ybr,RES] = cat_vol_resize( Yc{1} + Yc{2}*2 + Ycsf,'reduceV'  ,tis.res_vx_vol,2,16,'meanm'); % weight WM a bit more
    Ybr = smooth3(Ybr) > .5;
    Ybr = cat_vol_morph( Ybr , 'ldo', 1); 
    Ybr = cat_vol_morph( Ybr , 'dc' , 2); 
    Yb  = min(1,max(0, 10*tan(cat_vol_resize( smooth3( Ybr ), 'dereduceV', RES) - 0.5) + .5));
    clear Ybr RES; 
    for ci = 4:5
      Yh2b   = Yc{ci} .* Yb; 
      Yc{3}  = Yc{3}  + Yh2b;
      Yc{ci} = Yc{ci} - Yh2b;
    end
    clear Yh2b
    for ci = 1:3
      Yb2h   = Yc{ci} .* (1-Yb); 
      Yc{4}  = Yc{4}  + Yb2h;
      Yc{ci} = Yc{ci} - Yb2h;
    end
    clear Yb2h;


    %% Position-based bone/skull corrections
    % In many cases the bone-head relation is about 50:50 or at least 
    % between 20:80 or 80:20.
    Ybad = cat_vbdist( single( Ybg ) , Yb  < 1 , vx_vol); 
    Ybrd = cat_vbdist( single( Yb )  , Ybg < 1 , vx_vol); 
    Ypt  = Ybad + Ybrd;
    Ypp  = min(1,Ybad ./ Ypt);
    % to allow higher values close to the center
    Yppc = cat_vol_smooth3X( min(1, min(Ybad*1.5,Ybrd*8) ./ (Ybrd+Ybad)) .^ max(1,min(2,Ybrd)),sx);
    

    % Shell tissue model:
    % use a simple label map to apply some "shell" filter
    % with intensities from low to high going into the object 
    tth = 0.8; 
    % start with certain voxels 
    if tis.highBG == 0
      ith  = mean([tismri.int.bone_cortex, tismri.head_muscle]); 
    else
      ith  = 0; % don't do anything 
    end
    Ypx  = single( -3 + (Ybg>.5) );
    Ypx  = Ypx + (Ypx==-3 & ~Ybox) .* ( ...
      + 2*cat_vol_morph(Yc{5}>tth & Ym>ith & Ypp<.5,'lo',1,vx_vol) ...
      + 3*(Yc{4}>.5 & Ypp>.33 & Yppc>.33) ...
      + 4*(Yc{3}>tth) + 5*(Yc{1}>tth) + 6*(Yc{2}>tth) );
    % further refinements
    Ypx((Ybox4>0 & Ym<ith) | Ym==0) = -2; % defaced regions or close to image bouundaries 
    if numel(Ya)
      Ypx(Ypx==-3 & Ypt<20 & cat_vol_smooth3X(Yc{4} .* Ypp,sx)>0.01 & Ypp>.33 & Yb<.5 & (Ym<tismri.head_muscle.*0.5.*Ypp.*2.*Yppc & Ym<2*Yppc .* (tis.CSF*.5+.5*tis.GM)) & (Ya{1}<9 | Ya{1}==12)) = 0;  
      Ypx(Ypx==-3 & Ypt<20 & cat_vol_smooth3X(Yc{4} .* Ypp,sx)>0.01 & Ypp>.66 & Yb<.5 & Ym<(3.*Yppc) & (Ym<tismri.head_muscle.*Yppc*.4 | Ym>tismri.head_muscle*1.2 ) & (Ya{1}<9 | Ya{1}==12)) = 0;  
      Ypx(Ypx==-3 & Ypt>20 & Ybrd>1 & Ym>tismri.head_muscle*1.2) = -1;
      Ypx(Ypx==-3 & Ypp<0.5 & Ybg<.5) = -1;
      Ypx(Ypx==-3 & Ypt<20 & Yb<.5 & Ypp>.5 & Yppc>.8 & (Ym<ith | Ym>tismri.head_muscle.*1.2.*(1+Ypp)) & (Ya{1}<9 | Ya{1}==12)) = 0;
    else
      %Ypx(Ypx==-3 & Ypt<20 & cat_vol_smooth3X(Yc{4} .* Ypp,sx)>0.01 & Ypp>.33 & Yb<.5 & (Ym<tismri.head_muscle.*0.5.*Ypp.*2.*Yppc & Ym<2*Yppc .* (tis.CSF*.5+.5*tis.GM)) ) = 0;  
      %Ypx(Ypx==-3 & Ypt<20 & cat_vol_smooth3X(Yc{4} .* Ypp,sx)>0.01 & Ypp>.66 & Yb<.5 & Ym<(3.*Yppc) & (Ym<tismri.head_muscle.*Yppc*.4 | Ym>tismri.head_muscle*1.2 )) = 0;   
      Ypx(Ypx==-3 & Ypt>20 & Ybrd>1 & Ym>tismri.head_muscle*1.2) = -1;
      Ypx(Ypx==-3 & Ypp<0.5 & Ybg<.5) = -1;
      %Ypx(Ypx==-3 & Ypt<20 & Yb<.5 & Ypp>.5 & Yppc>.8 & (Ym<ith | Ym>tismri.head_muscle.*1.2.*(1+Ypp))) = 0;
    end
    %%
   
  %  Ypx(Ybrd > 10  & Ym > 0.4)   = -1; % bone should have low intensity if it is far from the brain (close to the brain it could be marrow) > head
  %  Ypx(Ypp < 0.33 & Ypx >- 0.5) = -1; % close to background > more likely head
  %  Ypx(Ypp > 0.66 & Ypx <- 0.5  & ~Yb & Ybrd<5) =  0; % close to brain > more likely bone
  %  Ypx(Ypp > 0.33 & Ypp <  0.65 ) = -3; % undefined region inbetween 
    
    %% first approximation between brain and background as skull-head layer
    Ypxa0 = min(3,Ypx); %Ypxa0(Ypxa0>=-1 & Ypxa0<=0) = -3;
    Ypxa0 = cat_vol_approx(Ypxa0 + 3,'rec',2) - 3;
    %
    Ypx(Ypx==-3) = Ypxa0(Ypx==-3);
    Ypxa = round( cat_vol_smooth3X(Ypx,sx) ); 
    for i=3:-1:-3
      [Ymskr,RES] = cat_vol_resize(Ypxa,'reduceV'  ,tis.res_vx_vol,4,16,'meanm');
      Ymskr  = cat_vol_morph( cat_vol_morph( round(Ymskr) >= i-.5 ,'ldo',1) , 'ldc', 4); 
      Ymsk   = min(1,max(0,4*tan(cat_vol_resize( cat_vol_smooth3X( Ymskr , 2 ), 'dereduceV', RES) - 0.25) + .25));
      Ypxa   = max(Ypxa, Ymsk * (i+2) - 3);  
    end

    
    % further morphometric corrections
    Ypx(Ypx==-1 & ~cat_vol_morph( Ypxa<-0.00 | Ypx==-1 ,'ldo',3,vx_vol) & Ypxa0>-0.5 & ~Ybox4 & Ypxa0>-.5) = 0; % remove head fragments
    

    % optimize segments
    %%{
    Ypx(Ypx>0   & ~cat_vol_morph( cat_vol_smooth3X(Yc{1}+Yc{2}+Yc{3},sx)>tth  ,'ldo',2,vx_vol))   = -3;   % remove brain fragments 
    Ypx(Ypx==-1 & ~cat_vol_morph( Ypxa0<-0.00 | Ypx==-1 ,'ldo',3,vx_vol) & Ypxa0>-0.5 & ~Ybox4 & Ypxa0>-.5) = 0; % remove head fragments
    Ypx(Ypx==-2 & ~cat_vol_morph(Ypxa0<-0.50 | Ypx==-2,'ldo',2,vx_vol) & Ypxa0>0) = -3; % remove background fragments
    Ypx(Ypx==-3 & Ybox & Ym<0.2  & Ypxa0<0) = -3;                   % extend background 
    %%}

    %% approximation as shell filter
    Ypxa = cat_vol_approx(Ypx + 3,'rec',2) - 3;
    Ypx(Ypx==-3 & Ypxa>-2 & Ypxa<-1 & Ym<.1) = -2; % background
    Ypx(Ypx==-3 & Ypxa>-1 & Ypxa<1)  = 0;
    Ypx(Ypx==-3 & Ypxa>-1 & Ypxa<1)  = 0;
    Ypx(Ypx>-3 & Ypx<1 & cat_vol_smooth3X( cat_vol_morph( smooth3(Yc{1}+Yc{2}+Yc{3})>tth  ,'l'),sx)) = 1; % remove brain fragments 
    [~,I] = cat_vbdist(single(Ypx>-3)); Ypx = Ypx(I); 

    %% adopts other tissues
    p2c = @(Y,x,y) 1 - min(1,max(0, abs(Y-x) * y )); 
    % head > background
    Yh2b  = Yc{5} .* p2c(Ypx,-2,10); % (Ypx>-2.1 & Ypx<-1.9); 
    Yc{6} = Yc{6} + Yh2b; 
    Yc{5} = Yc{5} - Yh2b; 
    clear Yh2b
    % background > head
    Yb2h  = Yc{6} .* p2c(Ypx,-1,10) .* (1-Yppc.^2); % (Ypx>-1.1 & Ypx<-0.9 & ~Ybox4); 
    Yc{5} = Yc{5} + Yb2h; 
    Yc{6} = Yc{6} - Yb2h; 
    clear Yb2h
    % head > skull 
    Yh2s  = Yc{5} .* p2c(Ypx,0,10) .* Yppc.^2 .*  min(1,max(0, 4*tan( cat_vol_smooth3X( (~Ybox4 & Ypxa>-.5 & Ypxa<0.5) | (Ybrd<3 & Ypp>0.8) , 1/mean(vx_vol) ) - 0.5) + .5));
    Yc{4} = Yc{4} + Yh2s; 
    Yc{5} = Yc{5} - Yh2s; 
    clear Yh2s
    % skull > head
    Yss   = cat_vol_smooth3X(Ypxa > -0.5 ,4); 
    Ys2h  = Yc{4} .* min(1,max(0, 4*tan( cat_vol_smooth3X( ( (Yss < 0 | Ypp < .5 ) & Ym > 0.3 & Ybrd < 10 ) | Ym==0 | Ypp < 0.2 | (Ybrd > 10  & Ym > 0.4) , 1/mean(vx_vol) ) - 0.5) + .5)); 
    Yc{5} = Yc{5} + Ys2h; 
    Yc{4} = Yc{4} - Ys2h; 
    clear Ys2h;
  else
    Yss =  0;
  end
  % debugging
  Ypxo = 1/5*Yco{5} + 2/5*Yco{4} + 3/5*Yco{3} + 4/5*Yco{1} + 5/5*Yco{2};   
  Ypxu = 1/5*Yc{5}  + 2/5*Yc{4}  + 3/5*Yc{3}  + 4/5*Yc{1}  + 5/5*Yc{2}; 



  %%
  % defacing >> background 
  Ydeface = cat_vol_morph(Yo==0,'l',[10,0.1])>0; % need multiple objects for face and ears 
  if sum(Ydeface(:)) > 10000
    for ci = 1:5, Ycn = Yc{ci} .* Ydeface; Yc{6} = Yc{6} + Ycn; Yc{ci} = Yc{ci} - Ycn; end
  end

  % background >> bone | head
  % background should be only outside the head ... not optimal for air
  Yc6    = Yc{6} .* cat_vol_morph(cat_vol_morph(Yc{6}<.5,'lc',1,vx_vol),'e',1,vx_vol); 
  clscor.BG2BN = sum(Yc6(:) .* ~(Yo(:)<=tismri.Tth(3) & (Ybraindist0s(:)<10))) / vxmm3 / tismri.TIV; 
  clscor.BG2HD = sum(Yc6(:) .*  (Yo(:)<=tismri.Tth(3) & (Ybraindist0s(:)<10))) / vxmm3 / tismri.TIV; 
  Yc{6}  = Yc{6} - Yc6;
  Yc{5}  = Yc{5} + Yc6 .* ~(Yo<=tismri.Tth(3) & (Ybraindist0s<10)); 
  Yc{4}  = Yc{4} + Yc6 .*  (Yo<=tismri.Tth(3) & (Ybraindist0s<10) & Yss>.5);
  % head >> background 
  % In children the head/bone class is too large and needs correction. 
  if tis.highBG == 0 % low intensity background
    if tis.weighting == 1 % T1w
      Yc6a = Yc{5} .* cat_vol_smooth3X(cat_vol_morph( ((Yc{6} + Yc{5})>.5 & ~Ydeface & Yo<tismri.Tth(3)), 'lo',1,vx_vol),sx);
    else % PDw
      Yc6a = Yc{5} .* cat_vol_smooth3X(cat_vol_morph( ((Yc{6} + Yc{5})>.5 & ~Ydeface & Yo<min(tismri.Tth(1:3))/4 ), 'lo',1,vx_vol),sx);
    end
  else
    % noisy background of MP2RAGE/MT sequences need some other (gradient-based) definition 
    % nothing that we can do here, right?
    Yc6a = 0; 
  end
  clscor.HD2BN = sum(Yc6a(:)) / vxmm3  / tismri.TIV; 
  Yc{5}  = Yc{5} - Yc6a;
  Yc{6}  = Yc{6} + Yc6a; 
  % final background closing 
  Yc6    = cat_vol_morph(cat_vol_morph(Yc{6}<.5,'lc')<.5,'e',1,vx_vol); 
  for ci = 1:5, Ycn = Yc{ci} .* Yc6; Yc{6} = Yc{6} + Ycn; Yc{ci} = Yc{ci} - Ycn; end


  %% head (+ background) ~ bonemarrow>> bone
  if ~exist('Ypp','var')
    % Position-based bone/skull corrections
    % In many cases the bone-head relation is about 50:50 or at least 
    % between 20:80 or 80:20.
    Ybg  = Yc{6}; 
    Yb   = Yc{1} + Yc{2} + Yc{3}; 
    Ybad = cat_vbdist( single( Ybg ) , Yb  < 1 , vx_vol); 
    Ybrd = cat_vbdist( single( Yb )  , Ybg < 1 , vx_vol); 
    Ypt  = Ybad + Ybrd;
    Ypp  = min(1,Ybad ./ Ypt);
    clear Ybg Yb Ybad Ybrd Ypt; 
  end

  if tis.headBoneType
    Yhead  = min(1,single(Yc{5} + Yc{6} + 0.1 * max(0,Ybraindist0s-15))); % no smoothing here!
    Yhead  = smooth3(cat_vol_morph(cat_vol_morph(smooth3(Yhead)>.7,'lo',4,vx_vol),'d',1,vx_vol)) .* Yhead .* (Ybraindist0s>0);
    cn = 0; for ci = 5:6, Ycn = Yc{ci} .* Yhead; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
    Ye{1} = zeros(size(Yc{1}),'single'); 
  else
    dmn = cat_stat_kmeans( Ybraindist0s(Yc{4}>.5 & ~cat_vol_morph(Yc{4}>.5,'e',1,vx_vol) & Ybraindist0s<(median(Ybraindist0s(Yc{4}>.5))*4) ) , 2 );
    Ybone = cat_vol_morph(cat_vol_morph(Yc{4}>.5 & Ybraindist0s>(dmn(1)*.5) & Ybraindist0s<(dmn(2)*1.2),'lc',1,vx_vol),'o',1,vx_vol);
    Ybone = Ybone .* min(1,Ypp*5) .* min(1,(1 - Ypp)*20);
    cn = 0; for ci = [1:3,5:6], Ycn = Yc{ci} .* ~Ybone; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
   
    %% get the veins
    Yhead  = min(1,single(Yc{5} + Yc{6} + 0.1 * max(0,Ybraindist0s-15))); % no smoothing here!
    Yhead2 = smooth3(cat_vol_morph(cat_vol_morph(smooth3(Yhead)>.7,'lo',4,vx_vol),'d',vx_vol)) .* Yhead .* (Ybraindist0s>0);
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
  Ybrain = (cat_vol_morph( Ybrain > .5,'ldo',3,vx_vol)); % & (Ybraindist0>10)) | (Ybrain & (Ybraindist0<10));           % remove skull
  Ybrain = smooth3(cat_vol_morph( Ybrain > 0.5,'dc',3,vx_vol));  % this closing includes meninges!
  Ybrain = max( Ybrain , cat_vol_smooth3X(Ybrain,2)>.6); % remove vein
  cn = 0; for ci = 1:3, Ycn = Yc{ci} .* Ybrain; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
  clscor.BR2BN = cn / vxmm3 / tismri.TIV; 
  
  % Yc4 outside the TPM definition >> head | background
  if all( tis.res_vx_vol < 1.5 )
    Yc4 = cat_vol_morph(Yc{4}<0.5,'l') & Yc{4}>0 & (Ybraindist0s>10) & Yss>.5; 
    clscor.BN2HD = sum(Yc4(:)) / vxmm3 / tismri.TIV; 
    Yc{5}  = Yc{5} + Yc{4} .* Yc4; 
    Yc{4}  = Yc{4} - Yc{4} .* Yc4; 
  end
  
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
  %}

  Ybrain = smooth3(cat_vol_smooth3X(Yc{3}+Yc{2}+Yc{1}>.05,4)>.4).^4;
  Yc4nc  = max(0,Yc{4} - (1-Ybrain));  
  Yc{3}  = Yc{3} + Yc4nc; 
  Yc{4}  = Yc{4} - Yc4nc;
  
  %% head >> BG 
  Yc5lth = min(.5,tismri.Tth(5)/tismri.Tth(2));
  Yc5    = Yc{5} .* ( Ym<min(0.2,Yc5lth)  &  abs(Ym)>0.001  &  cat_vol_morph(Yc{5}+Yc{6} > .5,'e',1,vx_vol)); 
  Yc{6}  = Yc{6} + Yc5; 
  Yc{5}  = Yc{5} - Yc5; 
  %Yc6    = Yc{5} .* (Ym>min(0.2,Yc5lth)) .* cat_vol_morph(Yc{5}+Yc{6},'e'); 
 
  %% avoid wholes in head
  Yc5c  = max(Yc{5},smooth3(cat_vol_morph(Yc{5}>.5,'dc',1.5,vx_vol)).^4);
  Yc{4} = Yc{4} - max(0,Yc5c - Yc{5}); 
  Yc{5} = Yc5c; 


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
      Yn  = cat_vol_localstat(single(Ym),(Yc{4}+Yc{3})>.5 & ~cat_vol_morph(Yc{1}>.5,'d',1,vx_vol),1,4);
      Yn2 = cat_vol_localstat(single(Ym),(Yc{4}+Yc{3})>.5 & ~cat_vol_morph(Yc{1}>.5,'d',1,vx_vol),2,4);
      Ym .* (Yc{4} + Yc{3}),Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}) .* Ym .* 3; 
      Ym .* (1- smooth3(Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}).*Ym * 3))
      Ym .* (1- smooth3(Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}).*Ym * 3)) % ... this is quite nice to remove mengs and BVs (but also affects the GM/WM boundary !
  end


  %% refinement of atlas map
  if 0 
    Ymg = cat_vol_div( Ym .* cat_vol_smooth3X(Yc{4}>.5,1));
    Ymg = cat_vol_smooth3X( Ymg ./ mean(Ymg(Ymg(:)<-0.1)) , 4) ;
    Yag = cat_vol_grad(Ya{1});
    Ya2 = Ya{1} .* (cat_vbdist(Yag,Yhead<1,vx_vol)>20 & Yc{4}>.5); 
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
    Ybe = cat_vol_morph(Ybe,'lc',1,vx_vol);
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

  %%
  Ypxo = 1/5*Yco{5} + 2/5*Yco{4} + 3/5*Yco{3} + 4/5*Yco{1} + 5/5*Yco{2}; 
  Ypxu = 1/5*Yc{5}  + 2/5*Yc{4}  + 3/5*Yc{3}  + 4/5*Yc{1}  + 5/5*Yc{2}; 

end
