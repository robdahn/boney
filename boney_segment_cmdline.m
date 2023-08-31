function boney_segment_cmdline(job,out,i,stime2,rerunstr)
% boney_segment_cmdline. Just print to the command line
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


  if job.opts.verb
    %% Todo
    %  -  Link volume with SEG overlay 
    % [-] Volume with surf overlay + "[MBI3D,TH3D]" Links
    %  *  create a table header / structure!
    %  *  save processing output mat/xml
    %
    %  -  link processed file (vol + surf?) wherefore?
    %  -  colorcoded values by CAT intensity map?
    %  *  add warning Letter!
    %  *  tiss error with class 
    %  *  thin bone thick head error?
    %
  
    job.opts.verb = 0;
    [~, Tline, ~, ~ ,MA,matm,mmatm] = boney_segment_prepare_print(job.files,job);

    %% create one line for the output table depending on the MAfn defintion
    for fni = 1:size(MA,1)
      if isfield( out(i) , MA{fni,2}) && isfield( out(i).(MA{fni,2}) , MA{fni,3} ) 
        if ischar( out(i).(MA{fni,2}).(MA{fni,3}) )
          matm{i,fni}  = out(i).(MA{fni,2}).(MA{fni,3});
          mmatm(i,fni) = 1; %mmatm(i,fni) = out(i).(MA{fni,2}).(MA{fni,3});
        elseif isnumeric( out(i).(MA{fni,2}).(MA{fni,3})(MA{fni,4}) )
          matm{i,fni}  = out(i).(MA{fni,2}).(MA{fni,3})(MA{fni,4});
          mmatm(i,fni) = out(i).(MA{fni,2}).(MA{fni,3})(MA{fni,4});
        else
          mmatm(i,fni) = 1; 
        end
      end
    end


    if out(i).CTseg
      %%
      matm{i,end-2} = matm{i,end-2} / 1000;
      matm{i,end-1} = matm{i,end-1} / 1000;
      
      %%
      cat_io_cprintf( job.opts.MarkColor( min(size(job.opts.MarkColor,1),max(1,floor( (matm{i,end-2} - 1) * ...
        size(job.opts.MarkColor,1)))),:),sprintf( Tline,i,...
          sprintf( sprintf('%%%ds',job.opts.snspace(1)-8) , spm_str_manip(job.files{i},['a' num2str(job.opts.snspace(1) - 14)])), ...
          ... spm_file( sprintf( sprintf('%%%ds',job.snspace(1)-8) , ...
          ... spm_str_manip(P{i},['a' num2str(job.snspace(1) - 14)])),'link',sprintf('spm_display(''%s'');',P{i})), ...
          matm{i,:},etime(clock,stime2)));
    else    
      cat_io_cprintf( job.opts.MarkColor( min(size(job.opts.MarkColor,1),max(1,floor( matm{i,end-2} * 3 / 9.5 * ...
        size(job.opts.MarkColor,1)))),:),sprintf( Tline,i,...
          sprintf( sprintf('%%%ds',job.opts.snspace(1)-8) , spm_str_manip(job.files{i},['a' num2str(job.opts.snspace(1) - 14)])), ...
          ... spm_file( sprintf( sprintf('%%%ds',job.snspace(1)-8) , ....
          ... spm_str_manip(P{i},['a' num2str(job.snspace(1) - 14)])),'link',sprintf('spm_display(''%s'');',P{i})), ...
          matm{i,:},etime(clock,stime2)));
    end
  
    fprintf(' %s',rerunstr);
  
    warn = cat_io_addwarning(1); if numel(warn)>0, cat_io_cprintf('warn',sprintf(' %dW',numel(warn))); end
    errn = cat_io_addwarning(2); if numel(errn)>0, cat_io_cprintf('err' ,sprintf(' %dE',numel(errn))); end
    fprintf('\n');

  end
end