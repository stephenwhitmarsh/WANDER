function [dataset_task dataset_rs] = WANDER_subjectinfo(rootpath)

i = 0;
if rootpath == 1
    rootpath_raw = 'D:\WANDER\data\raw';
else
    rootpath_raw = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/WANDER/data/raw';
end

for isubject = 1 : 26
    subject_dir = sprintf('wander%02d_s%02d',isubject,isubject);
    date_dir    = dir([rootpath_raw,filesep,subject_dir,filesep '16*']);
    date_dir    = date_dir.name;    
    for irun = 1 : 4
        if isubject == 21 && irun == 4
            run_fname = sprintf('s%02drun%01d_trans_sss.fif',isubject,irun);
        else
            run_fname = sprintf('s%02drun%01d_trans_tsss.fif',isubject,irun);
        end
        %         run_fname = sprintf('s%02drun%01d.fif',isubject,irun);
        dataset_task{isubject,irun} = [rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep,run_fname];
    end
    try
        EGG_placement{isubject} = dir([rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep '*.jpg']);
        EGG_placement{isubject} = [rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep,EGG_placement{isubject}.name];
    catch end
    nrofsets(isubject)      = 4;
    bad_fnames{isubject}    = dir([rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep '*.bad']);
    for ifile = 1 : size(bad_fnames{isubject},1)
        i = i + 1;
        flist{i} = [rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep,bad_fnames{isubject}(ifile).name];   
    end
end

i = 0;
for isubject = 1 : 26
    subject_dir = sprintf('wander%02d_s%02d',isubject,isubject);
    date_dir    = dir([rootpath_raw,filesep,subject_dir,filesep '16*']);
    date_dir    = date_dir.name;
    run_fname   = sprintf('s%02dresting_state_tsss.fif',isubject);
    dataset_rs{isubject} = [rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep,run_fname];
    try
        EGG_placement{isubject} = dir([rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep '*.jpg']);
        EGG_placement{isubject} = [rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep,EGG_placement{isubject}.name];
    catch
    end
    nrofsets(isubject)      = 4;
    bad_fnames{isubject}    = dir([rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep '*.bad']);
    for ifile = 1 : size(bad_fnames{isubject},1)
        i = i + 1;
        flist{i} = [rootpath_raw,filesep,subject_dir,filesep,date_dir,filesep,bad_fnames{isubject}(ifile).name];   
    end
end
