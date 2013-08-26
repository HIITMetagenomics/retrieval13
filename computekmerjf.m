function status = computekmerjf(k,fastqPath,jellyfishPath,matlab)
% This function computes the k-mer counts of a fastq file using jellyfish
% Input: k denotes which k-mer to be computed, put a number between 2-31
%        fastqPath is the path where the fastq file is stored
%         output is stored in the same directory with same name and _k.txt
%           e.g. /home/ss.fastq and k=3 would produce in /home/ss_3.txt
%        jellyfishPath is the the path of the jellyfish shell file
%        matlab is an optional flag, if true then store result in a .m file
% Output: status is either true (if execution went fine) or error message
% Copyright Sohan Seth

if nargin == 3
    matlab = false;
end
try
    cell101 = @(x)(char(x{1,1}));
    dirName = cell101(regexp(fastqPath,'(.+/).+\.fastq','tokens'));
    fileName = cell101(regexp(fastqPath,'.+/([^/])+\.fastq','tokens'));
    outputFileName = [fileName,'-',num2str(k)];
    if k > 9
        [status result] = system([jellyfishPath,' count -m ',num2str(k),' -o ',dirName,'temp -c 3 -s 1000000 -t 32 ',fastqPath]);
        [status result] = system([jellyfishPath,' merge -o ',dirName,outputFileName,' ',dirName,'temp*']);
        [status result] = system([jellyfishPath,' dump -o ',[dirName,outputFileName,'.txt'],' ',dirName,outputFileName]);
        [status result] = system(['rm ',dirName,'temp*']);
        [status result] = system(['rm ',dirName,outputFileName]);
    else
        [status result] = system([jellyfishPath,' count -m ',num2str(k),' -o ',dirName,'temp -c 3 -s 1100000 -t 32 ',fastqPath]);
        [status result] = system([jellyfishPath,' dump -o ',[dirName,outputFileName,'.txt'],' ',dirName,'temp_0']);
        [status result] = system(['rm ',dirName,'temp_0']);
    end
catch err
    status = err; return
end

if matlab
    file = [dirName,outputFileName,'.txt'];
    try
        fid = fopen(file);
        [status noLines] = system(['wc -l < ',file]);
        noLines = str2num(noLines);
        kmerId = cast(zeros(noLines/2,1),'uint64');
        kmerCount = cast(zeros(noLines/2,1),'uint64');
        for countLine = 1:noLines/2
            C = fgetl(fid);
            kmerCount(countLine,1) = str2num(C(2:end));
            C = fgetl(fid);
            C = double(C) - 64;
            C(C == 7) = 0; %G
            C(C == 20) = 2; %T
            kmerId(countLine,1) = cast(C * 4.^(k-1:-1:0)'+1,'uint64');
        end
        fclose(fid);
        save([dirName,outputFileName,'.mat'],'kmerId','kmerCount');
    catch err
        status = err;
    end
end