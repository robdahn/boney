function [Theader,Tline,Tavg, Cheader, MA, matm, mmatm] = ...
  boney_segment_prepare_print(P,job)
%prepare_print. Create table elements for the in-line report and fst print.
%
% [Theader,Tline,Tavg, Cheader, MAfn0, MAfn, matm, mmatm] = ...
%   boney_segment_prepare_print(P,job)
%  
%  P          .. input files
%  job        .. SPM batch structure 
%  Theader    .. header line of the table with all variables names
%  Tline      .. empty subject column 
%  Tavg       .. empty average column
%  Cheader    .. cell header 
%  MAfn       .. cell table with header
%  matm       .. main result table for values
%  mmatm      .. main result table for marks/ratings
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


% TODO: 
% * final update of table entries and variable names
% * add further parameter in overview
% * T1 values 
%   + (brain values - expert/developer)
%   + bone values
%      MAT ( mBC , mBM , mBMN)
%      VOL EXOR SURF (bone-cortex + bone-marrow + bone-thickness + head-thickness) 
%          ([v/s]BC , ~BM , ~BMN , vBT , v

  % Cell structure to define the table header fieldnames and the datatype
  % of the variable and if a seperator (|) should be used. 
% ############# RD202308: We have update the entries and naminges with the finally selected most relevant measures.   
  
  global boned %#ok<GVMIS> 
  ex  = boned.expertgui>0;  % show expert values
  exc = job.opts.bmethod>0 & job.opts.classic;   % show classic 
  exm = job.opts.bmethod;   % show vol and/or surface results
  MA  = {
  ... name    field0    field1          i format  sep dips
    'Tw'      'tis'     'weightingn'    1    's'  0   1
    'Tbg'     'tis'     'highBGn'       1    's'  0   1
    'Tfat'    'tis'     'headFatTypen'  1    's'  0   1
    'Tbone'   'tis'     'headBoneTypen' 1    's'  0   1
    'Tres'    'tis'     'res_RES'       1    'f'  0   1   
    'Tcnr'    'tis'     'seg8CNR'       1    'f'  0   1
    'Tbg'     'tis'     'background'    1    'f'  1   ex>1
    'Tcsf'    'tis'     'CSF'           1    'f'  0   ex>1
    'Tgm'     'tis'     'GM'            1    'f'  0   ex>1
    'Thead'   'tis'     'head'          1    'f'  0   ex>1
    'vBmar'   'vROI'    'bonemarrow'    3    'f'  1   exm>=1
    'vBcor'   'vROI'    'bonecortex'    4    'f'  0   exm>=1 
    'vBth'    'vROI'    'bonethickness' 3    'f'  0   exm>=1
    'vHth'    'vROI'    'headthickness' 4    'f'  0   exm>=1 
    'vHmed'   'vROI'    'head'          1    'f'  0   exm>=1 
    'sBmar'   'sROI'    'bonemarrow'    3    'f'  1   exm==2
    'sBcor'   'sROI'    'bonecortex'    4    'f'  0   exm==2
    'sBth'    'sROI'    'bonethickness' 3    'f'  0   exm==2
    'sHth'    'sROI'    'headthickness' 4    'f'  0   exm==2 
    'sHmed'   'vROI'    'head'          1    'f'  0   exm==2 
    'Tbcor'   'tis'     'bonecortex'    1    'f'  1   1
    'Tbmar'   'tis'     'bonemarrow'    1    'f'  0   1
    'Tbdns'   'tis'     'bonedensity'   1    'f'  0   1
    'Bmed'    'tismri'  'bone_med'      1    'f'  1   exm==0
    'Bmedc'   'classic' 'bone_med'      1    'f'  0   exc
    };
  MA( [MA{:,7}]<1 , : ) = [];

  % Creat a string that can be used with the sprintf function to format the
  % different table variables in a nice way.
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',job.opts.snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',job.opts.snspace(1)-8);
  Tavg    = sprintf('%%%ds:',job.opts.snspace(1)-1);
  for fi = 1:size(MA,1)
      Cheader = [Cheader MA{fi,1}]; %#ok<AGROW> 
      MAfni   = strrep( MA{fi,1} ,'_','');
    if MA{fi,6}, sep = ' |'; else, sep = ''; end
    Theader = sprintf(sprintf('%%s%%s%%%ds' ,job.opts.snspace(2)),Theader, sep, MAfni ); 
    switch MA{fi,5}
      case {'s','d'}
        Tline   = sprintf('%s%s%%%d%s', Tline , sep, job.opts.snspace(2), MA{fi,5});
        Tavg    = sprintf('%s%s%%%d%s', Tavg  , sep, job.opts.snspace(2), MA{fi,5});
      otherwise
        Tline   = sprintf('%s%s%%%d.%d%s',Tline , sep, job.opts.snspace(2), job.opts.snspace(3) .* (MA{fi,5}=='f'),MA{fi,5});
        Tavg    = sprintf('%s%s%%%d.%d%s',Tavg  , sep, job.opts.snspace(2), job.opts.snspace(3) .* (MA{fi,5}=='f'),MA{fi,5});
    end
  end
  % add field for the subject processing time 
  Tline   = sprintf('%s%%4.0fs',Tline);
  

  % main result table
  matm    = num2cell(nan(numel(P),size(MA,1)));
  mmatm   = 10.5*ones(numel(P),size(MA,1));
 

  % print title
  if job.opts.verb
    fprintf('\nBone Preprocessing:\n');
    if job.opts.verb 
      verbosestr = {'No','Yes','Yes - Details'};
      pmethodstr = {'previous','SPM','CAT'};
      tmpstr     = {'Adult TPM','Children TPM'};
      bmethodstr = {'SPMmat8','Volume-based','Surface-based', ...
        'Volume-based (old version without refinement)','Volume-based (old version with refinement)'};
      reportstr  = {'No','Yes - Basic','Yes - Details','Yes - Details'};
      %reportstr = {'Table','Table + Volumes','Table + Volumes + Surfaces'};
      fprintf('  Tissue Probability Map:   %d (%s)\n', job.opts.ctpm,     tmpstr{job.opts.ctpm});
      fprintf('  Preprocessing   method:   %d (%s)\n', job.opts.pmethod,  pmethodstr{job.opts.pmethod+1});
      fprintf('  Bone processing method:   %d (%s)\n', job.opts.bmethod,  bmethodstr{job.opts.bmethod+1});
      fprintf('  Verbose processing:       %d (%s)\n', job.opts.verb,     verbosestr{job.opts.verb + 1});
      if ex
        fprintf('  Rerun preprocessing:      %d (%s)\n', job.opts.prerun,  verbosestr{job.opts.prerun  + 1});
        fprintf('  Rerun bone processing:    %d (%s)\n', job.opts.rerun,   verbosestr{job.opts.rerun   + 1});
        fprintf('  Affine registration:      %d (%s)\n', job.opts.affreg,  verbosestr{job.opts.affreg  + 1});
        fprintf('  Strong bias correction:   %d (%s)\n', job.opts.bias,    verbosestr{job.opts.bias    + 1});
        fprintf('  Refine segmentation:      %d (%s)\n', job.opts.refine,  verbosestr{job.opts.refine  + 1});
        fprintf('  Use sub directories:      %d (%s)\n', job.opts.subdirs, verbosestr{job.opts.subdirs + 1});
        fprintf('  Use sub directories:      %d (%s)\n', job.opts.subdirs, verbosestr{job.opts.subdirs + 1});
        fprintf('  Classic measures:         %d (%s)\n', job.opts.classic, verbosestr{job.opts.classic + 1});
        fprintf('  Resolution limit:         %0.2f  \n', job.opts.reslim);
        if job.opts.bmethod>1
          fprintf('  Reduce-factor:            %dx (1-org,2-half,...) \n', job.opts.reduce);
        end
        fprintf('  Bone atlas:               %s\n'     , job.opts.Patlas{1} );
        fprintf('  Bone mask:                %s\n'     , job.opts.Pmask{1} );
        fprintf('  Write ouput: \n');
        fprintf('    Write volumes:          %d (%s)\n', job.output.writevol, reportstr{job.output.writevol + 1});
        fprintf('    Write surfaces:         %d (%s)\n', job.output.writevol, reportstr{job.output.writevol + 1});
        fprintf('    Report:                 %d (%s)\n', job.output.report, reportstr{job.output.report});
      else
        fprintf('  Report:                   %d (%s)\n', job.output.report, reportstr{job.output.report});
      end
    end
    fprintf('\n%s\n%s\n',  Theader,repmat('-',size(Theader)));  
  end
end
