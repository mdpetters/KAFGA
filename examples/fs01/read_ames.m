function [r] = read_ames(file)
    fileID = fopen(file);
    tline = fgetl(fileID);
    c = strsplit(tline);
    nx = str2num(char(c(1)));
    for i = 1:nx-1; tline = fgetl(fileID); end
    
    j = 1;
    while ~feof(fileID)
        tline = fgetl(fileID);
        c = str2num(char(strsplit(tline)));
        r.x1(j) = c(1);
        r.x2(j) = c(2);
        r.g1(j) = c(3);
        r.g2(j) = c(4);
        j = j+1;
    end
    fclose(fileID);
end
