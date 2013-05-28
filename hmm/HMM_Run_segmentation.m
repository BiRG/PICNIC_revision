
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


function HMM_Run_segmentation


% function HMM_Run_segmentation.
% Input: reads in data from /data/cancer/raw to perform normalisation
% Output:  Segmentation data produced in /output/sample_name directory

% Set Global Variables
global G
global H
global b

global alpha
global beta
global c
global seg_info
global p
global params

% Load reference files
currDir=pwd;
t1=[currDir '/data/cancer/normalised'];
t2='*.feature_intensity';
cell_names = Utils_ListNames(t1,t2);
no_names = size(cell_names,2);

% load in config.mat which contains all the parameters necessary
load('config.mat','chr_info');
seg_info.chr_start = chr_info(:,1);
seg_info.chr_end = chr_info(:,2);
seg_info.Tc = seg_info.chr_end - seg_info.chr_start + 1;
seg_info.no_chr = length(seg_info.chr_start);
seg_info.no_SNPs = seg_info.chr_end(seg_info.no_chr);


load('config.mat','SNP_ind');
load('config.mat','SNP_pos');
load('config.mat','p_coeffs');
seg_info.p = 1-p_coeffs;
snp2=SNP_ind;
pos2=SNP_pos;
clear SNP_ind;
clear SNP_pos;
global SNP_ind;
global SNP_pos;
seg_info.SNP_ind = snp2;
seg_info.SNP_pos = pos2;

p(:,1,1) = seg_info.p.^2;
p(:,2,1) = seg_info.p.*(1-seg_info.p);
p(:,2,2) = seg_info.p.*(1-seg_info.p);
p(:,1,2) = (1-seg_info.p).^2;

for chr_no = 1:seg_info.no_chr,
    SNP_ind{chr_no} = seg_info.SNP_ind(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),1);
    SNP_pos{chr_no} = seg_info.SNP_pos(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),1);
end

for name_no=1:size(cell_names,2)

    % Load data
    inputFile=[t1 '/' cell_names{name_no}];
    vals = dlmread(inputFile,',');
    for chr_no = 1:seg_info.no_chr,
        SNP_id{chr_no} = vals(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),1);
        H{chr_no} = vals(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),2);
        G{chr_no} = vals(seg_info.chr_start(chr_no):seg_info.chr_end(chr_no),3);
    end

    % Run HMM
     HMM_InitSeg;

    % Write params
    t1=[currDir '/output/' cell_names{name_no} ];
    if (exist(t1,'dir'))==0
        mkdir(t1);
    end

    save([t1,'/params_',cell_names{name_no},'.mat'],'seg_info');

    % Plot figures
    upos = [1 1 1600 700];
    close all;
    for chr_no = 1:seg_info.no_chr,
        Plot_Chrom_Raw(chr_no);
        u = gcf;
        set(u,'Position',upos);
        saveas(u,[t1,'/fig_',num2str(chr_no),'_',cell_names{name_no},'.fig']);
        close all;
    end
    f = figure('visible','off');
    Plot_Genome(f);
    saveas(f,[t1,'/genome_fig_',cell_names{name_no},'.fig']);
    close all;
	
	% DCW -- below repairs a bug on the second iteration due to t1
    %        getting clobbered; needs a better solution
    t1=[currDir '/data/cancer/normalised'];
end

return
