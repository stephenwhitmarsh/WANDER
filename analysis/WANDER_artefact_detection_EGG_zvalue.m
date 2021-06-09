function [artdef] = WANDER_artefact_detection_EGG_zvalue(isubject,force,rootpath,restingstate)

% EGG artefact detection based on ztransform of second derivative
if restingstate == 1
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_EGG_rs.mat'];
    else
        fname_artdef = ['/home/swhitmarsh/WANDER/data/artefacts/s' num2str(isubject) '_artdef_EGG_rs.mat'];
    end
else
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_EGG.mat'];
    else
        fname_artdef = ['/home/swhitmarsh/WANDER/data/artefacts/s' num2str(isubject) '_artdef_EGG.mat'];
    end
end

if exist(fname_artdef,'file') && force ~= 1
    fprintf('Returning EGG artefact definitions \n');
    load(fname_artdef);
else
    fprintf('EGG artefact definition not found, creating it now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    [dataset_task, dataset_rs] = WANDER_subjectinfo;
    
    EGG = WANDER_filter_EGG_broadband(isubject,0,rootpath,restingstate);
    FFT = WANDER_FFT_EGG(isubject,0,rootpath,restingstate);
    
    if exist(fname_artdef,'file')
        load(fname_artdef);
    end
    
    if restingstate == 1
        
        EGG_bb = WANDER_filter_EGG_broadband(isubject,0,rootpath,restingstate);
        EGG_bb.trial{1} = zscore(EGG.trial{1},0,2);
        EGG = WANDER_filter_EGG(isubject,0,rootpath,restingstate);
        EGG.trial{1} = zscore(EGG.trial{1},0,2);
        
        %smooth data
        EGG.trial{1}(3,:) = smooth(EGG.trial{1}(3,:),5);
        
        %differentials
        EGG.trial{1}(5,:) = [0 diff(EGG.trial{1}(3,:),1,2)];
        EGG.trial{1}(6,:) = [0 0 diff(EGG.trial{1}(3,:),2,2)];
        artperiod         = artedge:size(EGG.trial{1},2)-artedge;
        
        %ztransform without edges
        EGG.trial{1}(5,:) = EGG.trial{1}(5,:) - mean(EGG.trial{1}(5,artperiod));
        EGG.trial{1}(5,:) = EGG.trial{1}(5,:) ./ std(EGG.trial{1}(5,artperiod));
        EGG.trial{1}(6,:) = EGG.trial{1}(6,:) - mean(EGG.trial{1}(6,artperiod));
        EGG.trial{1}(6,:) = EGG.trial{1}(6,:) ./ std(EGG.trial{1}(6,artperiod));
        
        %smooth data
        EGG.trial{1}(5,:) = smooth(EGG.trial{1}(5,:),5);
        EGG.trial{1}(6,:) = smooth(EGG.trial{1}(6,:),5);
        
        %define artifacts
        artval =zeros(1,size(EGG.trial{1},2));
        artval(artperiod) = EGG.trial{1}(6,artperiod) > threshold;
        
        left = find(artval==1)-1;
        left = left(left>0);
        artval(left) = 1;
        right = find(artval==1)+1;
        right = right(right<length(artval));
        artval(right) = 1;
        
        artbeg = find(diff([0 artval])== 1);
        artend = find(diff([artval 0])==-1);
        artifact = [artbeg(:) artend(:)];
        artifact_ext = [artbeg(:)- artpad artend(:) + artpad];
        
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(7,1,1);
        plot(EGG.trial{1}(1,:));    axis tight; title([num2str(isubject) ' orig']);
        ax = axis;
        axis([ax(1),ax(2),-5,5]);
        line([ax(1) ax(2)],[4 4],'color','r');
        line([ax(1) ax(2)],[-4 -4],'color','r');
        
        subplot(7,1,2);
        plot(EGG.trial{1}(2,:));    axis tight; title([num2str(isubject) ' phase']);
        ax = axis;
        axis([ax(1),ax(2),-5,5]);
        line([ax(1) ax(2)],[4 4],'color','r');
        line([ax(1) ax(2)],[-4 -4],'color','r');
        
        if ~isempty(artifact)
            for iart = 1:size(artifact,1)
                patch([artifact_ext(iart,1) artifact_ext(iart,1) artifact_ext(iart,2) artifact_ext(iart,2)],[-7 7 7 -7],'g','FaceAlpha',0.3,'Edgecolor','none');
                patch([artifact(iart,1) artifact(iart,1) artifact(iart,2) artifact(iart,2)],[-7 7 7 -7],'r','FaceAlpha',0.3,'Edgecolor','none');
            end
        end
        
        subplot(7,1,3);
        plot(EGG.trial{1}(3,:));    axis tight; title([num2str(isubject) ' amplitude']);
        ax = axis;
        axis([ax(1),ax(2),-5,5]);
        line([ax(1) ax(2)],[4 4],'color','r');
        line([ax(1) ax(2)],[-4 -4],'color','r');
        
        subplot(7,1,4);
        plot(EGG.trial{1}(5,:));    axis tight; title([num2str(isubject) ' diff(amplitude)']);
        ax = axis;
        axis([ax(1),ax(2),-5,5]);
        line([ax(1) ax(2)],[4 4],'color','r');
        line([ax(1) ax(2)],[-4 -4],'color','r');
        
        subplot(7,1,5);
        plot(EGG.trial{1}(6,:));    axis tight; title([num2str(isubject) ' diff2(amplitude)']);
        ax = axis;
        axis([ax(1),ax(2),-threshold*1.5,threshold*1.5]);
        line([ax(1) ax(2)],[threshold threshold],'color','r','linestyle',':');
        line([ax(1) ax(2)],[-threshold -threshold],'color','r','linestyle',':');
        line([artedge artedge],[-threshold*1.5 threshold*1.5],'color','k');
        line([ax(2)-artedge ax(2)-artedge],[-threshold*1.5 threshold*1.5],'color','k');
        
        if ~isempty(artifact)
            for iart = 1:size(artifact,1)
                patch([artifact_ext(iart,1) artifact_ext(iart,1) artifact_ext(iart,2) artifact_ext(iart,2)],[-threshold*1.5 threshold*1.5 threshold*1.5 -threshold*1.5],'g','FaceAlpha',0.3,'Edgecolor','none');
                patch([artifact(iart,1) artifact(iart,1) artifact(iart,2) artifact(iart,2)],[-threshold*1.5 threshold*1.5 threshold*1.5 -threshold*1.5],'r','FaceAlpha',0.3,'Edgecolor','none');
            end
        end
        
        subplot(7,1,6);
        plot(EGG.trial{1}(4,:));    axis tight; title([num2str(isubject) ' narrow filter']);
        ax = axis;
        axis([ax(1),ax(2),-5,5]);
        line([ax(1) ax(2)],[4 4],'color','r');
        line([ax(1) ax(2)],[-4 -4],'color','r');
        
        subplot(7,1,7);
        plot(EGG_bb.trial{1}(4,:)); axis tight; title([num2str(isubject) ' broad filter']);
        ax = axis; hold;
        axis([ax(1),ax(2),-5,5]);
        line([ax(1) ax(2)],[4 4],'color','r');
        line([ax(1) ax(2)],[-4 -4],'color','r');
        
        artdef = artifact_ext;
        
        print(fig,'-painters','-dpng',['d:\analysis\WANDER\images\EGG\overview_s' num2str(isubject) '.png']);
        
    else
        
        for ipart = 1:4
            EGG{ipart}.trial{1}         = zscore(EGG{ipart}.trial{1},0,2);
            
            cfg                         = [];
            cfg.channel                 = {FFT.max_chan,[FFT.max_chan '_phase']};
            cfg.preproc.demean          = 'yes';
            cfg.viewmode                = 'vertical';
            cfg.blocksize               = 60*5;
            
            if exist(fname_artdef,'file')
                cfg.artfctdef.EGG.artifact  = artdef{ipart}.EGG.artifact;
            else
                cfg.artfctdef.EGG.artifact  = [];
            end
            
            cfg.selectfeature           = 'EGG';
            cfg                         = ft_databrowser(cfg,EGG{ipart});
            artdef{ipart}               = cfg.artfctdef;
        end
    end
    save(fname_artdef,'artdef');
end

