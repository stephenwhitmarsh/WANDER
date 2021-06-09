
function [artefact] = WANDER_artdef_window(isubject,force,rootpath)

% PARAMETERS
MEGwindow = 1;
EGGwindow = 1;

if isempty(rootpath)
    rootpath = 'd:\analysis\WANDER\data\';
end

fname_artdef_window = [rootpath filesep 'artefacts\s' num2str(isubject) '_windowed.mat'];

if exist(fname_artdef_window,'file') && force~=1
    fprintf('Returning windowed artefact definition\n');
    artefact = load(fname_artdef_window);
else
    fprintf('Windowed artefact definition not found, creating it now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    WANDER_subjectinfo;
       
    % read artifact definitions
    artdef_EOG = WANDER_blink_detection(isubject,0);       
    artdef_MEG = WANDER_artefact_detection_MEG(isubject,0);
    artdef_EGG = WANDER_artefact_detection_EGG(isubject,0);
    
    % read trial data for sampleinfo & timeaxis
    data_EGG   = WANDER_redefine_EGG_to_probe(isubject,0);

    % extend artifact intervals according to predefined filtering window
    for ipart = 1 : 4
        for itrial = 1:size(artdef_EGG{ipart}.EGG.artifact,1)
            artdef_ext_EGG{ipart}(itrial,1)     = artdef_EGG{ipart}.EGG.artifact(itrial,1)      - data_EGG{ipart}.fsample*EGGwindow/2;
            artdef_ext_EGG{ipart}(itrial,2)     = artdef_EGG{ipart}.EGG.artifact(itrial,2)      + data_EGG{ipart}.fsample*EGGwindow/2;
        end
        for itrial = 1:size(artdef_MEG{ipart}.MEG.artifact,1)
            artdef_ext_MEG{ipart}(itrial,1)     = artdef_MEG{ipart}.MEG.artifact(itrial,1)      - data_EGG{ipart}.fsample*MEGwindow/2;
            artdef_ext_MEG{ipart}(itrial,2)     = artdef_MEG{ipart}.MEG.artifact(itrial,2)      + data_EGG{ipart}.fsample*MEGwindow/2;
        end
        for itrial = 1:size(artdef_MEG{ipart}.muscle.artifact,1)
            artdef_ext_muscle{ipart}(itrial,1)  = artdef_MEG{ipart}.muscle.artifact(itrial,1)   - data_EGG{ipart}.fsample*MEGwindow/2;
            artdef_ext_muscle{ipart}(itrial,2)  = artdef_MEG{ipart}.muscle.artifact(itrial,2)   + data_EGG{ipart}.fsample*MEGwindow/2;
        end
        for itrial = 1:size(artdef_EOG{ipart}.zvalue.artifact,1)
            artdef_ext_EOG{ipart}(itrial,1)     = artdef_EOG{ipart}.zvalue.artifact(itrial,1)   - data_EGG{ipart}.fsample*MEGwindow/2;
            artdef_ext_EOG{ipart}(itrial,2)     = artdef_EOG{ipart}.zvalue.artifact(itrial,2)   + data_EGG{ipart}.fsample*MEGwindow/2;
        end
    end
    
    % create long array of artefact samples
    clear *_samples *_samples_ext
    for ipart = 1 : 4
        EGG_samples{ipart}          = [];
        MEG_samples{ipart}          = [];
        EOG_samples{ipart}          = [];
        muscle_samples{ipart}       = [];
        EGG_samples_ext{ipart}      = [];
        MEG_samples_ext{ipart}      = [];
        EOG_samples_ext{ipart}      = [];
        muscle_samples_ext{ipart}   = [];
        
        try
            for itrial = 1:size(artdef_ext_EGG{ipart},1)
                EGG_samples{ipart}      = [EGG_samples{ipart} artdef_EGG{ipart}.EGG.artifact(itrial,1):artdef_EGG{ipart}.EGG.artifact(itrial,2)];
                EGG_samples_ext{ipart}  = [EGG_samples_ext{ipart} artdef_ext_EGG{ipart}(itrial,1):artdef_ext_EGG{ipart}(itrial,2)];
            end
        catch
            fprintf('No EGG artifacts found in block %d \n',ipart);
        end
        try
            for itrial = 1:size(artdef_ext_MEG{ipart},1)
                MEG_samples{ipart}      = [MEG_samples{ipart} artdef_MEG{ipart}.MEG.artifact(itrial,1):artdef_MEG{ipart}.MEG.artifact(itrial,2)];
                MEG_samples_ext{ipart}  = [MEG_samples_ext{ipart} artdef_ext_MEG{ipart}(itrial,1):artdef_ext_MEG{ipart}(itrial,2)];
            end
        catch
            fprintf('No MEG artifacts found in block %d \n',ipart);
        end
        try
            for itrial = 1:size(artdef_ext_EOG{ipart},1)
                EOG_samples{ipart}      = [EOG_samples{ipart} artdef_EOG{ipart}.zvalue.artifact(itrial,1):artdef_EOG{ipart}.zvalue.artifact(itrial,2)];
                EOG_samples_ext{ipart}  = [EOG_samples_ext{ipart} artdef_ext_EOG{ipart}(itrial,1):artdef_ext_EOG{ipart}(itrial,2)];
            end
        catch
            fprintf('No EOG artifacts found in block %d \n',ipart);
        end
        try
            for itrial = 1:size(artdef_ext_muscle{ipart},1)
                muscle_samples{ipart}       = [muscle_samples{ipart} artdef_MEG{ipart}.muscle.artifact(itrial,1):artdef_MEG{ipart}.muscle.artifact(itrial,2)];
                muscle_samples_ext{ipart}   = [muscle_samples_ext{ipart} artdef_ext_muscle{ipart}(itrial,1):artdef_ext_muscle{ipart}(itrial,2)];
            end
        catch
            fprintf('No muscle artifacts found in block %d \n',ipart);
        end
    end
    
    % create overlap between trial samples overlap and artefact samples
    for ipart = 1 : 4
        for itrial = 1:size(data_EGG{ipart}.time,2)
            fprintf('Calculating overlap with EGG artifacts, block%d, trial%d\n',ipart,itrial);
            artefact.overlap_EGG{ipart}{itrial}          = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),EGG_samples{ipart});
            artefact.overlap_EGG_ext{ipart}{itrial}      = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),EGG_samples_ext{ipart});
        end
        for itrial = 1:size(data_EGG{ipart}.time,2)
            fprintf('Calculating overlap with MEG artifacts, block%d, trial%d\n',ipart,itrial);
            artefact.overlap_MEG{ipart}{itrial}          = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),MEG_samples{ipart});
            artefact.overlap_MEG_ext{ipart}{itrial}      = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),MEG_samples_ext{ipart});
        end
        for itrial = 1:size(data_EGG{ipart}.time,2)
            fprintf('Calculating overlap with EOG artifacts, block%d, trial%d\n',ipart,itrial);
            artefact.overlap_EOG{ipart}{itrial}          = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),EOG_samples{ipart});
            artefact.overlap_EOG_ext{ipart}{itrial}      = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),EOG_samples_ext{ipart});
        end      
        for itrial = 1:size(data_EGG{ipart}.time,2)
            fprintf('Calculating overlap with muscle artifacts, block%d, trial%d\n',ipart,itrial);
            artefact.overlap_muscle{ipart}{itrial}       = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),muscle_samples{ipart});
            artefact.overlap_muscle_ext{ipart}{itrial}   = ismember(data_EGG{ipart}.sampleinfo(itrial,1):data_EGG{ipart}.sampleinfo(itrial,2),muscle_samples_ext{ipart});
            
            % mask trials that have a combined rejection of samples above threshold ratio
            if mean(artefact.overlap_muscle_ext{ipart}{itrial}) > 0.2
                artefact.overlap_muscle_ext{ipart}{itrial} = true(size(artefact.overlap_muscle_ext{ipart}{itrial}));
            end
        end
    end
    
    % report on overlap
    for ipart = 1 : 4
        
        maxlength(ipart) = 0;
        for itrial = 1:size(data_EGG{ipart}.time,2)
            maxlength(ipart) = max(maxlength(ipart),size(artefact.overlap_EGG{ipart}{itrial},2));
        end
        
        matrix_EGG{ipart}       = zeros(size(data_EGG{ipart}.time,2),maxlength(ipart));
        matrix_MEG{ipart}       = zeros(size(data_EGG{ipart}.time,2),maxlength(ipart));
        matrix_EOG{ipart}       = zeros(size(data_EGG{ipart}.time,2),maxlength(ipart));
        matrix_jump{ipart}      = zeros(size(data_EGG{ipart}.time,2),maxlength(ipart));
        matrix_muscle{ipart}    = zeros(size(data_EGG{ipart}.time,2),maxlength(ipart));
        
        for itrial = 1:size(data_EGG{ipart}.time,2)
            matrix_EGG{ipart}(itrial,~artefact.overlap_EGG{ipart}{itrial})           = 1;
            matrix_EGG{ipart}(itrial,artefact.overlap_EGG_ext{ipart}{itrial})        = 3;       
            matrix_EGG{ipart}(itrial,artefact.overlap_EGG{ipart}{itrial})            = 2;          
            matrix_MEG{ipart}(itrial,~artefact.overlap_MEG{ipart}{itrial})           = 1;
            matrix_MEG{ipart}(itrial,artefact.overlap_MEG_ext{ipart}{itrial})        = 3;
            matrix_MEG{ipart}(itrial,artefact.overlap_MEG{ipart}{itrial})            = 2;       
            matrix_EOG{ipart}(itrial,~artefact.overlap_EOG{ipart}{itrial})           = 1;
            matrix_EOG{ipart}(itrial,artefact.overlap_EOG_ext{ipart}{itrial})        = 3;       
            matrix_EOG{ipart}(itrial,artefact.overlap_EOG{ipart}{itrial})            = 2; 
            matrix_muscle{ipart}(itrial,~artefact.overlap_muscle{ipart}{itrial})     = 1;
            matrix_muscle{ipart}(itrial,artefact.overlap_muscle_ext{ipart}{itrial})  = 3; 
            matrix_muscle{ipart}(itrial,artefact.overlap_muscle{ipart}{itrial})      = 2;
        end
    end
    
    fig = figure;   
    cmap = [0 0 0; 0 1 0; 0 0 1; 1 0 0];
    for ipart = 1 : 4
        h = subplot(4,4,(ipart-1)*4+1);
        image(matrix_EGG{ipart}+1);     title('EGG artefacts');       set(h,'XTick',[]);  set(h,'YTick',[]); colormap(cmap);     
        h = subplot(4,4,(ipart-1)*4+2);
        image(matrix_MEG{ipart}+1);     title('Visual inspection');   set(h,'XTick',[]);  set(h,'YTick',[]);
        h = subplot(4,4,(ipart-1)*4+3);
        image(matrix_EOG{ipart}+1);     title('Blinks');              set(h,'XTick',[]);  set(h,'YTick',[]);
        h = subplot(4,4,(ipart-1)*4+4);
        image(matrix_muscle{ipart}+1);  title('Muscle artefacts');    set(h,'XTick',[]);  set(h,'YTick',[]);
    end
    set(fig,'PaperSize',[11*2 8.5*2 ]);
    set(fig,'PaperPosition',[0 0 11*2 8.5*2]);
    set(fig,'PaperOrientation','landscape');
    print(fig,'-dpdf',['d:\analysis\WANDER\images\artefacts\s' num2str(isubject) '_artefacts.pdf']);
    print(fig,'-dpng',['d:\analysis\WANDER\images\artefacts\s' num2str(isubject) '_artefacts.png']);
    
    % save artefact definitions
    save(fname_artdef_window,'artefact','-v7.3');
end
