
% Copyright (c) 2008, Genome Research Ltd (GRL).
% All rights reserved.
% Author: Chris Greenman <cdg@sanger.ac.uk>
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above
%       copyright notice, this list of conditions and the following
%       disclaimer in the documentation and/or other materials
%       provided with the distribution.
%     * Neither the name of the <organization> nor the
%       names of its contributors may be used to endorse or promote
%       products derived from this software without specific prior
%       written permission.
%
% THIS SOFTWARE IS PROVIDED BY GRL ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED.  IN NO EVENT SHALL GRL BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGE.


function HMM_RunB(name_no,init_val,delta_val,varargin)

% function HMMRunB : runs the programs for batch processing the HMM.
% 	Inputs: name_no: number of files to process - default is 1.
%		init_val: Copy number intensity of deletions.
%		delta_val: Change in intensity between successive copy numbers.
%		Optional arguments: no_levels - the maximum segmented copy number.   

optargin = size(varargin,2);

if ischar(name_no),
    name_no = str2num(name_no);
end
if ischar(init_val),
    init_val = str2num(init_val);
end
if ischar(delta_val),
    delta_val = str2num(delta_val);
end

% Set Global Variables
global G
global H
global b
global SNP_ind
global SNP_pos
global alpha
global beta
global c
global seg_info
global p
global params
global states
global genotypes
global stats
global genes
global segments

currDir=pwd;
t1=[currDir '/data/cancer/normalised'];
t2='*.feature_intensity';



% Load reference files
cell_names = Utils_ListNames(t1,t2);
no_names = size(cell_names,2);

load('config.mat','chr_info');
load('config.mat','p_coeffs');
p_raw= 1-p_coeffs;
p(:,1,1) = p_raw.^2;
p(:,1,2) = (1-p_raw).^2;
p(:,2,1) = p_raw.*(1-p_raw);
p(:,2,2) = p(:,2,1);

% Construct levels to use
load([currDir '/output/',cell_names{name_no},'/params_',cell_names{name_no},'.mat']);

if(optargin==0)
	seg_info.no_levels = 15;
else
	seg_info.no_levels=str2num(varargin{1});
end

seg_info.levels = init_val+[0:seg_info.no_levels-1]*delta_val;
seg_info.no_genos_ind = [];

for level_no = 1:seg_info.no_levels,
    seg_info.no_genos_ind(level_no) = floor((level_no-1+0.1)/2)+1;
end
seg_info.no_genos_max = max(seg_info.no_genos_ind);



for chr_no = 1:seg_info.no_chr,
    SNP_ind{chr_no} = seg_info.SNP_ind(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),1);
    SNP_pos{chr_no} = seg_info.SNP_pos(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),1);
end
    
 
% Load data
vals = dlmread([currDir '/data/cancer/normalised/',cell_names{name_no}],',');

for chr_no = 1:seg_info.no_chr,
    SNP_id{chr_no} = vals(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),1);
    H{chr_no} = max(eps,vals(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),2));
    G{chr_no} = max(eps,vals(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),3));
end
    

disp(['Starting HMM']);
% Run HM
HMM_InitParamSet;
disp(['Starting BaumWelch']);
HMM_BaumWelch;
disp('starting Viterbi');
HMM_Viterbi;
HMM_GenoType;
HMM_GenomeStats;
HMM_GeneAnalysis;

disp(['Done HMM']); 
 
% Write params
mkdir([currDir '/output2/',cell_names{name_no}]);
save([currDir '/output2/',cell_names{name_no},'/params_',cell_names{name_no},'.mat'],'seg_info','params','states');
  
% Define Public Output Information Matrix
output_mat(1:seg_info.no_SNPs,1:17) = 0;
for chr_no = 1:seg_info.no_chr,
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),1) = SNP_id{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),2) = H{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),3) = G{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),4) = (H{chr_no}-params.init)/params.delta;
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),5) = states.cnv{chr_no}-1;
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),6) = states.geno{chr_no}-1;
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),8) = 2/pi*atan(1+(states.cnv{chr_no}-1)*params.eta)-0.5;
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),7) = (1/2-max(0,(states.geno{chr_no}-1)./(states.cnv{chr_no}-1))).*(4/pi*atan(1+(states.cnv{chr_no}-1)*params.eta)-1);
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),9) = 1*(states.geno{chr_no}==1);
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),10:11) = genotypes.geno_class{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),12) = genotypes.geno_state_change{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),13) = genotypes.geno_prob{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),14) = genotypes.geno_cond_prob{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),15) = genotypes.geno_HET_prob{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),16) = genotypes.geno_A_LOH_prob{chr_no};
    output_mat(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),17) = genotypes.geno_B_LOH_prob{chr_no};
end
mat2min = [abs(0.5+output_mat(:,8)-output_mat(:,3)),abs(0.5-output_mat(:,8)-output_mat(:,3))];
mat2min = [mat2min,abs(0.5+output_mat(:,7)-output_mat(:,3)),abs(0.5-output_mat(:,7)-output_mat(:,3))];
qual = min(mat2min,[],2);
qualval(seg_info.no_chr+1) = sum(qual)/(seg_info.chr_end(seg_info.no_chr));
qualval(seg_info.no_chr+2) = params.gamma;
qualval(seg_info.no_chr+3) = params.sd1;
qualval(seg_info.no_chr+4) = params.sd2;
for chr_no = 1:seg_info.no_chr,
    qualval(chr_no) = sum(qual(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no)))/(seg_info.chr_end(chr_no)-seg_info.chr_start(chr_no)+1);
end

dlmwrite([currDir '/output2/','/output_',cell_names{name_no},'.csv'],output_mat,'delimiter',',','precision',7);
csvwrite([currDir '/output2/',cell_names{name_no},'/complexity_',cell_names{name_no},'.csv'],stats.complexity);
csvwrite([currDir '/output2/',cell_names{name_no},'/LOH_',cell_names{name_no},'.csv'],stats.LOH);
csvwrite([currDir '/output2/',cell_names{name_no},'/genes_',cell_names{name_no},'.csv'],[genes.raw,genes.mean,genes.del,genes.change]);
dlmwrite([currDir '/output2/',cell_names{name_no},'/qual_',cell_names{name_no},'.csv'],qualval','delimiter',',','precision',7);
dlmwrite([currDir '/output2/',cell_names{name_no},'/segments_',cell_names{name_no},'.csv'],segments,'delimiter',',','precision',7);

 
% Plot figures
upos = [1 1 1600 700];
close all;
for chr_no = 1:seg_info.no_chr,
    Plot_Chrom_HMM(chr_no);
    u = gcf;
    set(u,'Position',upos);
    saveas(u,[currDir '/output2/',cell_names{name_no},'/fig_',num2str(chr_no),'_',cell_names{name_no},'.fig']);
    close all;
end
Plot_Genome;
u = gcf;
set(u,'Position',upos);
saveas(u,[currDir '/output2/',cell_names{name_no},'/genome_fig_',cell_names{name_no},'.fig']);
close all;


return
