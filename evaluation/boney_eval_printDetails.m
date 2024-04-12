%boney_extractBoneMeasures_printDetails. Batch to create a figure
%
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

    fname    = measure; 
    corrtype = 'Pearson'; %'Spearman'; 'Pearson'
    marksize = max(20,(out.TIV{1} - 800) / 1000 * 30); % bmi;
    markname = 'TIV';
    tblue = '\color[rgb]{0 0 1}'; tblack = '\color[rgb]{0 0 0}'; 

 
    fprintf('Print details: ')
    clear RHO PVAL
    for mi = 1:numel(measure)
      %%
        for si = [unique(site)] % loop over sites
            if  ~all(isnan(out.(measure{mi}){1}(:))) && contains(Pcsv{1}{csvi},'ukb')
              % create figure
              if printfg
                fg=10; fgh = figure(fg); clf(fg); fgh.Position(3:4) = [1400 1100]; 
                fgh.Name = sprintf('mt%02d_%s_%s_site%d_n%d',matonly,matonlyn{matonly+1},measure{mi},si,numel(age)); 
                ha = annotation('textbox',[0.005 0.92 0.99 0.08],'string',sprintf( ...
                  '%d-%s%s: %s - site %d, N=%d, radius~%s ',matonly,matonlyn{matonly+1},strrep(resdir,'_','\_'), ...
                  strrep(measure{mi},'_','\_'), si, numel(age), markname));
                ha.FontSize             = 14;
                ha.FontWeight           = 'bold';
                ha.HorizontalAlignment  = 'center';
                ha.EdgeColor            = 'none'; 
              end

              %% print subplot names
              for bmdi = 1:numel(para)
                % measures/definitions of each subplot
                clear bmdx
                switch bmdi
                  % bone measures
                  case 1,  bmdx = bmdh;         bmdname = 'BMDhead';
                  case 2,  bmdx = bmdt;         bmdname = 'BMDtotal';
                  case 3,  bmdx = bmdf;         bmdname = 'BMDfemur';
                  case 4,  bmdx = age;          bmdname = 'age';
                  % fat measures
                  case 5,  bmdx = [];           bmdname = 'Histogram'; 
                  case 6,  bmdx = VAT;          bmdname = 'VAT';
                  case 7,  bmdx = ASAT;         bmdname = 'ASAT';
                  case 8,  bmdx = fatp;         bmdname = 'rFAT';
                  case 9,  bmdx = bmi;          bmdname = 'BMI';
                  case 10, bmdx = waist;        bmdname = 'waist';
                  % brain measures
                  case 11, bmdx = out.rGMV{1};  bmdname = 'rGMV';
                  case 12, bmdx = out.rWMV{1};  bmdname = 'rWMV';
                  case 13, bmdx = out.rCMV{1};  bmdname = 'rCMV';
                  case 14, bmdx = out.TIV{1};   bmdname = 'TIV';
                  case 15, bmdx = sittingheight; bmdname = 'sittingheight';
                  % other
                  case 16, bmdx = birthweight;  bmdname = 'birthweight'; 
                  case 17, bmdx = timesummer;   bmdname = 'time out summer'; 
                  case 18, bmdx = walk;         bmdname = 'walk';
                  case 19, bmdx = smoking;      bmdname = 'smoking'; 
                  case 20, bmdx = alc;          bmdname = 'alcohol'; 
                  otherwise, continue
                end
                if all(isnan(bmdx)), bmdx = zeros(size(bmdx)); end
                if any(size(bmdx) ~= size(sex)), bmdx = bmdx'; end
                fx = 1 + (bmdi>3);
         
                if si == 0; msk = site>0; else, msk = site==si; end % set up to process all sites for si == 0
               % msk = msk & age>60 & age<70; 
                

                if printfg
                  % set up subplot position
                  px=5; py=4; [xx,yy] = ind2sub([px,py],bmdi);
                  sp = subplot('Position',[(xx-1)/px+.03, .95*(py-yy)/py+.05, 1/px*.78, 1/py*.72]);
                  if bmdi==5
%                    histogram(out.(measure{mi}){1}, max(6,min(100,round(log(numel(bmdh)) * 6 )))); 
                    hh(1)=histogram(out.(measure{mi}){1}(sex==2), max(6,min(100,round(log(numel(bmdh)) * 3 )))); hold on
                    hh(2)=histogram(out.(measure{mi}){1}(sex==1), max(6,min(100,round(log(numel(bmdh)) * 3 )))); 
                    hh(1).FaceColor = [ 1    0         0      ]; hh(1).FaceAlpha = .6; 
                    hh(2).FaceColor = [ 0    0.4470    0.7410 ]; hh(2).FaceAlpha = .6; 
                    hh(2).BinEdges  = hh(1).BinEdges; 
                    title(sprintf('Histogram %s',strrep(fname{mi},'_','\_'))); 
                    xlabel(strrep(fname{mi},'_','\_')); 
                    ylim(ylim.*[1 1.2]); grid on
                    continue; 
                  end

                  % plot values
                  sh = scatter(bmdx(sex==1 & msk), out.(measure{mi}){1}(sex==1 & msk),(marksize(sex==1 & msk)-10)*2,'filled'); 
                  sh.MarkerEdgeColor = [ 0    0.4470    0.7410]; sh.MarkerFaceColor = [ 0    0.4470    0.7410];sh.MarkerEdgeAlpha = 0.1; sh.MarkerFaceAlpha = 0.3; 
                  hold on
                  sh = scatter(bmdx(sex==2 & msk), out.(measure{mi}){1}(sex==2 & msk),(marksize(sex==2 & msk)-10)*2,'filled'); 
                  sh.MarkerEdgeColor = [ 1 0 0]; sh.MarkerFaceColor = [ 1 0 0]; sh.MarkerEdgeAlpha = 0.1; sh.MarkerFaceAlpha = 0.3; 
                else
                  if bmdi==5, continue; end
                end

                fitval2 = 'poly1'; % fitval2 = fitval{mi};    
                msk  = msk(:) & ~isinf(out.(measure{mi}){1}(:)) & ~isnan(out.(measure{mi}){1}(:)); 
                msk2 = ~isinf(out.(measure{mi}){1}(:)); if any(size(msk) ~= size(msk2)), msk2 = msk2'; end
                msk3 = ~isnan(out.(measure{mi}){1}(:)); if any(size(msk) ~= size(msk3)), msk3 = msk3'; end
                msk  = msk & msk2 & msk3; 
                if any(size(msk) ~= size(sex)), msk = msk'; end

                % fits
                if any( size( out.(measure{mi}){1} ) ~= size( bmdx' ) ),  out.(measure{mi}){1} = out.(measure{mi}){1}'; end
                [RHOa,PVALa] = corr(bmdx(~isnan(bmdx) &          msk)',out.(measure{mi}){1}(~isnan(bmdx) &          msk),'type',corrtype); 
                [RHOm,PVALm] = corr(bmdx(~isnan(bmdx) & sex==1 & msk)',out.(measure{mi}){1}(~isnan(bmdx) & sex==1 & msk),'type',corrtype); 
                [RHOf,PVALf] = corr(bmdx(~isnan(bmdx) & sex==2 & msk)',out.(measure{mi}){1}(~isnan(bmdx) & sex==2 & msk),'type',corrtype);
                if printfg
                  warning off
                  [f,r(3)]     = fit(bmdx( ~isnan(bmdx) &          msk)',double(out.(measure{mi}){1}(~isnan(bmdx) &          msk)),...
                    fitval2,'Normalize','on','Robust','Bisquare'); ph=plot(f); ph.Color = [0.5    0.5    0.5]; 
                  [f,r(1)]     = fit(bmdx( ~isnan(bmdx) & sex==1 & msk)',double(out.(measure{mi}){1}(~isnan(bmdx) & sex==1 & msk)),...
                    fitval2,'Normalize','on','Robust','Bisquare'); ph=plot(f); ph.Color = [0    0.4470    0.7410]; 
                  [f,r(2)]     = fit(bmdx( ~isnan(bmdx) & sex==2 & msk)',double(out.(measure{mi}){1}(~isnan(bmdx) & sex==2 & msk)),...
                    fitval2,'Normalize','on','Robust','Bisquare'); ph=plot(f); ph.Color = [1    0.0000    0.0000];  
                  warning on
                end

                % save for other figure
                RHO(mi,bmdi,1) = RHOa; PVAL(mi,bmdi,1) = PVALa; 
                RHO(mi,bmdi,2) = RHOm; PVAL(mi,bmdi,2) = PVALm; 
                RHO(mi,bmdi,3) = RHOf; PVAL(mi,bmdi,3) = PVALf; 

                if printfg
                  % colors PVAL
                  if     PVALa<1e-10, signia='\bf\color[rgb]{1 0 0}';  elseif PVALa<1e-5, signia='\color[rgb]{.5 0 .5}'; 
                  elseif PVALa<0.01,  signia='\color[rgb]{0 0 1}';     elseif PVALa<0.05, signia='\color[rgb]{0  0 .5}'; 
                  else,  signia='\color[rgb]{0.4 0.4 0.5}'; 
                  end 
                  if     PVALm<1e-10, signim='\bf\color[rgb]{1 0 0}';  elseif PVALm<1e-5, signim='\color[rgb]{.5 0 .5}'; 
                  elseif PVALm<0.01,  signim='\color[rgb]{0 0 1}';     elseif PVALm<0.05, signim='\color[rgb]{0  0 .5}'; 
                  else, signim='\color[rgb]{0.4 0.4 0.5}'; 
                  end 
                  if     PVALf<1e-10, signif='\bf\color[rgb]{1 0 0}';  elseif PVALf<1e-5, signif='\color[rgb]{.5 0 .5}'; 
                  elseif PVALf<0.01,  signif='\color[rgb]{0 0 1}';     elseif PVALf<0.05, signim='\color[rgb]{0  0 .5}'; 
                  else,  signif='\color[rgb]{0.4 0.4 0.5}'; 
                  end 
                  % colors RHO if signficant
                  tblacka = tblack;  tblackm = tblack;  tblackf = tblack; 
                  if PVALa<0.05
                    if     abs(RHOa)>0.9/fx, signra='\color[rgb]{1 0 0}'; elseif abs(RHOa)>0.8/fx, signra='\color[rgb]{.5 0 .5}'; 
                    elseif abs(RHOa)>0.5/fx, signra='\color[rgb]{0 0 1}'; elseif abs(RHOa)>0.3/fx, signra='\color[rgb]{0 0 .5}'; 
                    elseif abs(RHOa)>0.1/fx, signra=''; else, signra='\color[rgb]{0.3 0.3 .4}'; tblacka = signra;
                    end 
                  else
                    signra='\color[rgb]{0.4 0.4 0.5}'; tblacka = signra;
                  end
                  if PVALm<0.05
                    if     abs(RHOm)>0.9/fx, signrm='\color[rgb]{1 0 0}'; elseif abs(RHOm)>0.8/fx, signrm='\color[rgb]{.5 0 .5}'; 
                    elseif abs(RHOm)>0.5/fx, signrm='\color[rgb]{0 0 1}'; elseif abs(RHOm)>0.3/fx, signrm='\color[rgb]{0 0 .5}'; 
                    elseif abs(RHOm)>0.1/fx, signrm=''; else, signrm='\color[rgb]{0.3 0.3 .4}'; tblackm = signrm;
                    end
                  else
                     signrm='\color[rgb]{0.4 0.4 0.5}'; tblackm = signrm;
                  end
                  if PVALf<0.05 
                    if     abs(RHOf)>0.9/fx, signrf='\color[rgb]{1 0 0}'; elseif abs(RHOf)>0.8/fx, signrf='\color[rgb]{.5 0 .5}'; 
                    elseif abs(RHOf)>0.5/fx, signrf='\color[rgb]{0 0 1}'; elseif abs(RHOf)>0.3/fx, signrf='\color[rgb]{0 0 .5}'; 
                    elseif abs(RHOf)>0.1/fx, signrf=''; else, signrf='\color[rgb]{0.3 0.3 .4}'; tblackf = signrf;
                    end 
                  else
                    signrf='\color[rgb]{0.4 0.4 0.5}'; tblackf = signrf;
                  end
                  % color main RHO in title
                  if PVALa<0.05
                    if     abs(RHOa)>0.6 || abs(RHOm)>0.8 || abs(RHOf)>0.8, signrt='\color[rgb]{0.0 0.0 0.1}'; 
                    elseif abs(RHOa)>0.4 || abs(RHOm)>0.6 || abs(RHOf)>0.6, signrt='\color[rgb]{0.1 0.1 0.2}'; 
                    elseif abs(RHOa)>0.2 || abs(RHOm)>0.4 || abs(RHOf)>0.4, signrt='\color[rgb]{0.2 0.2 0.3}'; 
                    elseif abs(RHOa)>0.1 || abs(RHOm)>0.2 || abs(RHOf)>0.2, signrt='\color[rgb]{0.3 0.3 0.4}'; 
                    else,                                                   signrt='\color[rgb]{0.4 0.4 0.5}'; 
                    end
                  else
                    signrt='\color[rgb]{0.4 0.4 0.5}'; 
                  end
  
                  % title, legend, etc. 
                  title(sprintf('%s %s (R^2=%0.3f, r=%s%0.3f%s, p=%s%0.0e%s)',signrt,strrep(bmdname,'_','\_'), r(3).adjrsquare, signra,RHOa,signrt, signia,PVALa,signrt)); 
                  xlabel(strrep(bmdname,'_','\_')); ylabel(strrep(fname{mi},'_','\_')); box on; grid on; ylim( ylim + [0 0.3*diff(ylim)]);
                  lg = legend( ...
                      sprintf('%sMale     (R^2=%0.3f, r=%s%0.3f%s, p=%s%0.0e%s)',tblackm,r(1).adjrsquare,signrm,RHOm,tblackm,signim,PVALm,tblackm), ...
                      sprintf('%sFemale (R^2=%0.3f, r=%s%0.3f%s, p=%s%0.0e%s)'  ,tblackf,r(2).adjrsquare,signrf,RHOf,tblackf,signif,PVALf,tblackf), ...
                      'Location','NorthWest','Interpreter','tex');    
                  lg.EdgeColor = min([1 1 1],1.5 * mean(cell2mat(cellfun(@eval,cat_io_strrep({tblackm;tblackf},{'\color[rgb]{','}'},{'[',']'}),'UniformOutput',false))));
                  ax = gca; 
                  ax.XColor = min([1 1 1],1.5 * min(cell2mat(cellfun(@eval,cat_io_strrep({signrt;tblackm;tblackf},{'\color[rgb]{','}'},{'[',']'}),'UniformOutput',false))));
                  ax.YColor = ax.XColor; 
                end

              end
            end
            if printfg
              resultsubdir = fullfile(resultdir,sprintf('mt%02d_%s_details_site%d_n%d',matonly,matonlyn{matonly+1},si,numel(age)));
              if ~exist(resultsubdir,'dir'), mkdir(resultsubdir); end
              saveas(gcf,fullfile(resultsubdir,sprintf('mt%02d_%s_site%d_n%d_%s.png',matonly,matonlyn{matonly+1},si,numel(age),measure{mi})));
            end
        end
    end


    fprintf('done.\n')