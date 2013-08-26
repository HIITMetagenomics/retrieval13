% This script generates test samples for metagenomics using metaSim
% The code modifies a single taxon profile to generate samples with different abundance
% This code requires metaSim and Jellyfish, install them and provide their paths to appropriate variables
% The generated files are stored in the storeDir with names storeName0, storeName1, ...
% Copyright Sohan Seth

% Error profile for empirical sequencing error in metaSim
errorfile = '/triton/ics/scratch/mi/metaSim/errorModel/errormodel-80bp.mconf';

% Samples for evaluating performance of string mining
flagJellyfish = 1; noSamples = 200; noReads = 5000000;
if flagJellyfish
    jellyfishPath = ['/triton/ics/scratch/mi/sohan/ProcessedData/Matlab'' code''/jellyfish-1.1.5/bin/jellyfish'];%'/triton/ics/scratch/mi/jellyfish-1.1.6/bin/jellyfish';
    storeDir = ['/triton/ics/scratch/mi/metaSim/fasta/stringMining-',date]; mkdir(storeDir);
    % storeDir = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r'; mkdir(storeDir);
    mprffile = '/triton/ics/scratch/mi/metaSim/project/HWFAs.mprf';
    storeName = 'HWFA-';
    flagExistingAbundanceProfile = 0; % Load abundance values from existing profile
    if flagExistingAbundanceProfile
        load fasta/stringMining-pHWFA-8g-46s.mat
    end
end

% Count number of lines in the taxon profile to decide number of species
[status noLines] = system(['wc -l < ',mprffile]);
noLines = str2num(noLines); noSpecies = floor(noLines/2); clear noLines

% Arbirary division of abundance profile in two classes
if flagSequenceError | flagJellyfish
    if ~flagExistingAbundanceProfile
        label = rand(noSamples,1) > 0.75; abundanceLevelList = zeros(noSamples,noSpecies);
    end
end

if flagJellyfish
    classA = rand(noSpecies,1); classB = classA; % Dirichlet prior
    % temp = randperm(noSpecies); temp2 = randperm(floor(noSpecies/2));
    temp = randperm(noSpecies); temp2 = randperm(floor(1*noSpecies/1));
    classB(temp(temp2)) = classA(temp(1:length(temp2)));
end

for countSample = 1:noSamples
    % Generate random abundance level $put special generation procedure here$
    
    if flagJellyfish
        if ~flagExistingAbundanceProfile
            if label(countSample)
                abundanceLevel = gamrnd(classA,1,size(classA));
            else
                abundanceLevel = gamrnd(classB,1,size(classB));
            end
            abundanceLevel = abundanceLevel / sum(abundanceLevel);
            abundanceLevel = floor(1000*abundanceLevel);
            abundanceLevelList(countSample,:) = abundanceLevel;
        else
            abundanceLevel = floor(1000*abundanceLevelList(countSample,:));
        end
    end
    
    status = changeAbundanceScript(mprffile,abundanceLevel);
    
    tic
    % Generate fna file
    [~, result] = system(['/home/seths1/metasim/MetaSim cmd -mg ',errorfile,...
        ' -r',num2str(noReads),' -d ',storeDir,' ',mprffile]);
    savedFile = regexp(result,'Read data saved as `(.+)''.','tokens'); savedFile = char(savedFile{1,1});
    
    movefile([storeDir,'/',savedFile],[storeDir,'/',storeName,num2str(countSample),'.fna'])
    display = regexp(result,'(Generated \d+ Reads)','tokens');
    fprintf('[%d sample generated. %s.]\n',countSample,char(display{1,1}));
    fprintf('[fna file generated %0.6f s]\n',toc)
    
    % Translate file in fastq format
    status = system(['./fasta2fastq ',storeDir,'/',storeName,num2str(countSample),'.fna ',storeDir,'/',storeName,num2str(countSample),'.fastq']);
    fprintf('[fastq file generated %0.6f s]\n',toc)
    
    % Compute k-mer count and save
    if flagSequenceError
        kmerCount = computeKmer(3,[storeDir,'/',storeName,num2str(countSample),'.fastq'],'all',[]);
        save([storeDir,'/',storeName,num2str(countSample),'.mat'],'kmerCount');
    end
    
    tic
    % Computer k-mer and save
    if flagJellyfish
        for countK = [3,8]
            status = computekmerjf(countK,[storeDir,'/',storeName,num2str(countSample),'.fastq'],jellyfishPath);
        end
    end
    fprintf('[k-mer computed %0.6f s]\n',toc)
    
    % Store status of the simulation
    fid = fopen([storeDir,'/status.txt'],'a');
    fprintf(fid,'%s',result);
    fclose(fid);
end
save([storeDir,'/abundance.mat'],'abundanceLevelList','label')
toc