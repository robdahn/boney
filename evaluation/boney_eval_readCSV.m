% boney_eval_readCSV. Read and evaluate CSV files to setup parameters.
%
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________



clear sex age height weight bmi site;
for csvi = 1:numel(Pcsv{1})
    if contains(Pcsv{1}{csvi},'IXI')
        csv  = cat_io_csv(Pcsv{1}{csvi},'','',struct('delimiter',',','komma','.'));
        siten{1} = 'IXI-Guys';
        siten{2} = 'IXI-HH';
        siten{3} = 'IXI-IO';
        for si = 1:numel(groups{1})
            %extracting demographics
            [pp,ff,ee] = spm_fileparts(groups{1}{si}); 
            if ~contains(pp,'IXI'), continue; end
            id = find( contains(csv(:,1),ff)==1 );
            try 
                sex(si)     = csv{id,3}; % 1 - male, 2 - female
                age(si)     = csv{id,13};
                height(si)  = csv{id,4}/100; 
                weight(si)  = csv{id,5}; 
                bmi(si)     = weight(si) / height(si).^2;  
            catch
                sex(si)     = nan; 
                age(si)     = nan; 
                height(si)  = nan; 
                weight(si)  = nan; 
                bmi(si)     = nan; 
            end
            if     ~isempty( strfind(ff,'Guys') ), site(si) = 1; 
            elseif ~isempty( strfind(ff,'HH') ),   site(si) = 2; 
            else,                                  site(si) = 3; 
            end
        end
    elseif contains(Pcsv{1}{csvi},'camcan')
        siten{4} = 'CamCan';
        csv = cat_io_csv(Pcsv{1}{csvi},'','',struct('delimiter',',','komma','.'));
        for si = 1:numel(groups{1})
            %extracting demographics
            [pp,ff,ee] = spm_fileparts(groups{1}{si}); 
            if ~contains(pp,'camcan'), continue; end
            id = find( contains(csv(:,1),ff(6:13))==1 );
            sex(si)     = csv{id,5}; % 1 - male, 2 - female
            age(si)     = csv{id,2};
            height(si)  = 1.8; 
            weight(si)  = 80; 
            bmi(si)     = 25;  
            site(si)    = 4;
        end
    elseif contains(Pcsv{1}{csvi},'NKI_RAW')
        siten{5} = 'NKI';
        csv = cat_io_csv(Pcsv{1}{csvi},'','',struct('delimiter','\t','komma','.'));
        for si = 1:numel(groups{1})
            %extracting demographics
            [pp,ff,ee] = spm_fileparts(groups{1}{si}); 
            if ~contains(pp,'NKI_RAW'), continue; end
            id = find( contains(csv(:,1),ff(6:14))==1 ); 
            sex(si)     = contains(csv{id,3},'MALE') + contains(csv{id,3},'FEMALE'); 
            age(si)     = csv{id,2};
            height(si)  = 1.8; 
            weight(si)  = 80; 
            bmi(si)     = 25;  
            site(si)    = 5; 
        end
     elseif contains(Pcsv{1}{csvi},'NKI_RS')
        siten{6} = 'NKI-RS';
        csv = cat_io_csv(Pcsv{1}{csvi},'','',struct('delimiter',',','komma','.'));
        for si = 1:numel(groups{1})
            %extracting demographics
            [pp,ff,ee] = spm_fileparts(groups{1}{si}); 
            if ~contains(pp,'NKI_RS'), continue; end
            id = find( [csv{2:end,1}] == str2double(ff(2:8)) , 1) + 1; 
            if ~isempty(id)
                sex(si)     = contains(csv{id,3},'male') + contains(csv{id,3},'female'); 
                age(si)     = csv{id,9};
                try
                    height(si)  = csv{id,5}/100; 
                    weight(si)  = csv{id,6}; 
                    bmi(si)     = weight(si) / height(si).^2;   
                catch
                    height(si)  = 1.8; 
                    weight(si)  = 80; 
                    bmi(si)     = 25;  
                end                
                site(si)    = 6;     
            end
        end 
    elseif contains(Pcsv{1}{csvi},'OASIS3')
        siten{7} = 'OASIS3';
        csv = cat_io_csv(Pcsv{1}{csvi},'','',struct('delimiter',',','komma','.'));
        for si = 1:numel(groups{1})
            %% extracting demographics
            [pp,ff,ee]  = spm_fileparts(groups{1}{si}); isCT = contains(ff,'CT_CTseg');

         %   if ~contains(pp,'OASIS3'), continue; end
            pps         = spm_str_manip(pp,'hht');
            pp          = spm_str_manip(pp,'ht');
            id          = find( contains( csv(:,1), pps(5:end) )==1 );
            sex(si)     = csv{id,5}; % 1 - male, 2 - female
            age(si)     = csv{id,3} + (str2double(pp(end-3:end)))./365;
            sub{si}     = pps; %ff(6:13);
            height(si)  = 1.8; 
            weight(si)  = 80; 
            bmi(si)     = 25;  
            site(si)    = 7; 
        end
        
      elseif contains(Pcsv{1}{csvi},'ukb')
        siten{8} = 'UKB';
        csv = cat_io_csv(Pcsv{1}{csvi},'','',struct('delimiter',',','komma','.'));
        for si = 1:numel(groups{1})
            %% extracting demographics
            [pp,ff,ee] = spm_fileparts(groups{1}{si}); 
          %  if ~contains(pp,'ukb'), continue; end
            pp2 = spm_str_manip(pp,'ht');
            if strcmp(pp2,'202301_MRIbones') && ~contains(pp,'testsetExt') && ~contains(pp,'subsample') 
              seps = find(ff=='_'); 
              id = find(contains( cellfun( @(x) sprintf('%d',x), csv(:,1), 'UniformOutput', false) ,ff(seps(4)+1+4:seps(5)-1))==1,1,'first');
            elseif strcmp(pp2,'202301_MRIbones') && strcmp( ff(1:4) , 'msub' )
              seps = find(ff=='_'); 
              id = find(contains( cellfun( @(x) sprintf('%d',x), csv(:,1), 'UniformOutput', false) ,ff(2+4:seps(1)-1))==1,1,'first');
            else
              id = find(contains( cellfun( @(x) sprintf('%d',x), csv(:,1), 'UniformOutput', false) ,ff(6:12))==1,1,'first');
            end
            if ~isempty(id)
                sex(si)           = csv{id,3}; % 1 - male, 2 - female
                age(si)           = csv{id,2};
                bmi(si)           = csv{id,4}; 
                bmdh(si)          = csv{id,5}; 
                bmdt(si)          = csv{id,6}; 
                age_meno(si)      = csv{id,7};
                HRT(si)           = csv{id,8};
                fIQ(si)           = csv{id,9};
                birthweight(si)   = csv{id,10};
                sunburnchild(si)  = csv{id,11};
                bmdf(si)          = csv{id,12};
                VAT(si)           = csv{id,13};
                ASAT(si)          = csv{id,14};
                fatp(si)          = csv{id,15};
                sittingheight(si) = csv{id,16};
                waist(si)         = csv{id,17};
                smoking(si)       = csv{id,18};
                timesummer(si)    = csv{id,19};
                timewinter(si)    = csv{id,20};
                happyness(si)     = csv{id,21};
                healthsat(si)     = csv{id,22};
                height(si)        = 1.8;
                weight(si)        = 80; 
                site(si)          = 8; 
            else
                sex(si)           = nan;
                age(si)           = nan; 
                bmi(si)           = nan; 
                bmdh(si)          = nan; 
                bmdt(si)          = nan; 
                age_meno(si)      = nan;
                HRT(si)           = nan;
                fIQ(si)           = nan;
                birthweight(si)   = nan;
                sunburnchild(si)  = nan;
                bmdf(si)          = nan;
                VAT(si)           = nan;
                ASAT(si)          = nan;
                fatp(si)          = nan;
                sittingheight(si) = nan;
                waist(si)         = nan;
                smoking(si)       = nan;
                timesummer(si)    = nan;
                timewinter(si)    = nan;
                happyness(si)     = nan;
                healthsat(si)     = nan;
                height(si)        = 1.8;
                weight(si)        = 80; 
                site(si)          = 8; 
            end
        end
    end
    % correct some outliers
    bmi(bmi<=0 | isinf(bmi)) = 1; 
    height(height<=0 | isinf(height)) = 1; 
    weight(weight<=0 | isinf(weight)) = 1;    
end
csvtab{csvi} = csv; 