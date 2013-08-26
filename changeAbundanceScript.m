function status = changeAbundanceScript(origFile,abundanceLevel)
% Supporting script to change existing .mprf file with new abundance values
% Input:    name of the fna file
%           new abundance values in a vector
% Output:   status is 0 if execution was fine, error message if error
% Copyright  Sohan Seth

try
    status = 0;
    path = regexp(origFile,'(.+/).+\.','tokens');
    name = regexp(origFile,'/([^/]+)\.','tokens');
    newFile =  [char(path{1,1}),char(name{1,1}),'2.mprf'];
    
    [status noLines] = system(['wc -l < ',origFile]);
    noLines = str2num(noLines);
    
    % Open original newFile and a new newFile to write the new abundance values
    fid = fopen(origFile); fidW = fopen(newFile,'w'); 
    countLine = 0;
    while 1
        C = fgetl(fid);
        countLine = countLine + 1;
        if C == -1
            break
        end
        % First line stays intact
        if countLine == 1
            fprintf(fidW,'%s',C);
            continue;
        end
        % Replace abundance values in the next lines
        if ~mod(countLine,2)
            fprintf(fidW,'\n%s',regexprep(C,'\d+ name "',[num2str(abundanceLevel(floor(countLine/2))),' name "']));
        end
        if mod(countLine,2)
            fprintf(fidW,'\n%s',regexprep(C,'\d+        taxid',[num2str(abundanceLevel(floor(countLine/2))),'        taxid']));
        end
    end
    % Close files
    fclose(fid); fclose(fidW);
    % Delete original and rename the new taxon profile
    delete(origFile); movefile(newFile,origFile);
catch err
    status = err;
end