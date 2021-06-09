function seed = WANDER_plantseeds(isubject,force)

fname_seeds = ['d:\analysis\WANDER\data\ICA\seeds.mat'];

if exist(fname_seeds,'file') && force ~= 1
    fprintf('Returning seed\n');
    load(fname_seeds);
else
    fprintf('Seeds not found, planting them now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    
    for iseed = 1 : 26
        seeds(iseed) = rand * 2^32;
    end
    save(fname_seeds,'seeds');
end
seed = seeds(isubject);
