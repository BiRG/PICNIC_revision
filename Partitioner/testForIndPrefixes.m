function testForIndPrefixes
%%
load('..\hmm\config\config.mat','SNP_ind');

% Because of a spare tab at the end of the header line, we need to read one
% additional empty column off the header line before reading the rest of
% the file.
rawASCIIDataFID = fopen('..\hmm\data\cancer\raw\20110805_A3_NT12.feature_intensity');
header = textscan(rawASCIIDataFID,'%s%s%s%s%s%s%s%s%s%s%s%s%s',1,'Delimiter','\t\n');
rawData = textscan(rawASCIIDataFID,'%s%d%d%d%d%d%d%d%d%d%d%d','Delimiter','\t\n');
fclose(rawASCIIDataFID);

probeIDs = rawData{1};
clear rawData;

true_SNP_inds = probeIDs(SNP_ind==1);
false_SNP_inds = probeIDs(SNP_ind==0);
clear probeIDs SNP_ind;

%%
for i = 1:length(true_SNP_inds)
    if (~strcmp(true_SNP_inds{i}(1:3), 'SNP'))
        disp '';
        disp ''; 'Found a non-SNP prefix''d ID...';
        disp '';
        break;
    end
end

%%
for i = 1:length(false_SNP_inds)
    if (~strcmp(false_SNP_inds{i}(1:2), 'CN'))
        disp '';
        disp ''; 'Found a non-CN prefix''d ID...';
        disp '';
        break;
    end
end

%%
for i = 1:length(true_SNP_inds)
    if (strcmp(true_SNP_inds{i}(1:2), 'CN'))
        disp '';
        disp ''; 'Found a CN prefix''d ID...';
        disp '';
        break;
    end
end

%%
for i = 1:length(false_SNP_inds)
    if (strcmp(false_SNP_inds{i}(1:3), 'SNP'))
        disp '';
        disp ''; 'Found a SNP prefix''d ID...';
        disp '';
        break;
    end
end

end
