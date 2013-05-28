function Features_2CnvGeno

%inputFile=[inputDir '/' inputFileName]
currDir=pwd;
t1=[currDir '/data/cancer/raw'];
t2='*.feature_intensity';

% do a check that the directory exists..
if (exist(t1,'dir'))
    cell_names = Utils_ListNames(t1,t2);
    addpath(t1);
else
    disp(['Directory ',t1, ' does not exist.']);
end

load('config.mat','no_features');
load('config.mat','SNP_ind');
load('config.mat','chr_info');


chr_start = chr_info(:,1);
chr_end = chr_info(:,2);
no_SNPs_both = chr_info(22,2);
no_SNPs_male = chr_info(24,2);
no_SNPs_female = chr_info(23,2);
no_SNPs = size(no_features,1);

load('config.mat','T_coeffs');
load('config.mat','A_coeffs');
load('config.mat','B_coeffs');

A0 = T_coeffs(:,1);
A1 = T_coeffs(:,2);
B0 = T_coeffs(:,3);
B1 = T_coeffs(:,4);

% Process each file in the working directory
num_cell_names = size(cell_names,2);
for fileNo=1:num_cell_names
    inputFile=[t1 filesep cell_names{fileNo}];
    if exist(inputFile,'file')~=2
       % somehow quit the program
        disp(['File ' inputFile ' does not exist. Please check the location and filename']);
        return;
    else
        SNP_id = dlmread(inputFile,'\t',[1,1,no_SNPs,1]);
        raw_data = dlmread(inputFile,'\t',[1,4,no_SNPs,11]);

    end
    
    raw_data = raw_data.*(raw_data>0);
    norm_consts = sum(raw_data,2)./(2*sum(no_features,2)-1+SNP_ind);
    raw_data = raw_data/sum(norm_consts)*no_SNPs_female;
    A = sum(raw_data(:,1:4),2)./no_features;
    B = sum(raw_data(:,5:8),2)./no_features;

    % Obtain Cnvs
    cnv = A_coeffs.*A + B_coeffs.*B;

    % Obtain Genos (Allelic Ratios)
    B_rescaled = (B-B0+B1)./B1;
    A_rescaled = (A-A0+A1)./A1;
    geno = 2/pi*atan(B_rescaled./A_rescaled);
    geno = abs(geno);
    geno = max(eps,geno);

    CNVGeno=[SNP_id,cnv,geno];
    try
        outdir=[currDir '/data/cancer/normalised/'];
        if (exist(outdir,'dir'))
            dlmwrite(['data/cancer/normalised/' cell_names{fileNo}],CNVGeno,'delimiter',',','precision',8);
        else
            disp(['File does not exist: ' outdir]);
        end
    catch ME
        disp(['Problems writing to data/cancer/normalised/' cell_names{file_no}]);
        return;
    end
end
return
