cd w:\WANDER\images\article
f = dir('*.fig');
for ifile = 1 : size(f,1)
    fig = openfig(f(ifile).name);
    [PATHSTR,fname,EXT] = fileparts(f(ifile).name);
    set(gcf, 'renderer', 'painters', 'paperpositionmode', 'auto');
    print('-dpdf',[fname '.pdf']);
    print('-deps',[fname '.eps']);

end

close all


