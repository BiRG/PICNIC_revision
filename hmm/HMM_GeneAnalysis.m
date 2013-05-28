
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

function HMM_GeneAnalysis

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
global params
global states
global genotypes
global stats
global genes

% Initialisations
del_prob_min = 0.01;
no_sims = 100;
eps = 1e-300;
chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;
Tc = seg_info.Tc;
no_chr = seg_info.no_chr;
no_SNPs = seg_info.no_SNPs;
no_levels = seg_info.no_levels;
no_genos_ind = seg_info.no_genos_ind;
vals = csvread('config/GeneFootPrint2Read.csv');
%vals = csvread('6k_Info/GeneFootPrint2Read.csv');
no_genes = size(vals,1);
gene_start = vals(:,2);
gene_end = vals(:,3);
gene_chr = vals(:,1);
genes.mean(1:no_genes,1) = 0;
genes.raw(1:no_genes,1) = 0;
genes.change(1:no_genes,1) = 0;
genes.del(1:no_genes,1) = 0;

% Need transition  transpose tensor
for level_no1 = 1:no_levels
    for geno_no1 = 1:no_genos_ind(level_no1),
        level_mat(level_no1,geno_no1) = level_no1-1;
        a_diag(level_no1,geno_no1) = params.chi(level_no1,geno_no1,level_no1,geno_no1);
        for level_no2 = 1:no_levels
            for geno_no2 = 1:no_genos_ind(level_no1),
                a_trans(level_no2,geno_no2,level_no1,geno_no1) = params.chi(level_no1,geno_no1,level_no2,geno_no2);
            end
        end
    end
end

% Calc mean ampl, state change prob and prob of deletion per gene
for gene_no = 1:no_genes,
    gene_no;
    gene_SNP1 = sum(gene_start(gene_no,1)>seg_info.SNP_pos(seg_info.chr_start(gene_chr(gene_no,1)):seg_info.chr_end(gene_chr(gene_no,1)),1));
    gene_SNP2 = 1+sum(gene_end(gene_no,1)>=seg_info.SNP_pos(seg_info.chr_start(gene_chr(gene_no,1)):seg_info.chr_end(gene_chr(gene_no,1)),1));
    chr_no = gene_chr(gene_no);
    t1 = max(1,gene_SNP1);
    t2 = min(gene_SNP2,Tc(chr_no)-1);
    if t2<t1,
        t1=t2;
    end
    
    % Mean Ampl
    denom1 = t2-t1+1;
    val2add = sum(squeeze(sum(sum(repmat(level_mat,[1,1,denom1]).*alpha{chr_no}(:,:,t1+1:t2+1).*beta{chr_no}(:,:,t1+1:t2+1),1),2))./c{chr_no}(t1+1:t2+1)');
    genes.mean(gene_no,1) = val2add/denom1;
    genes.raw(gene_no,1) = mean(H{chr_no}(t1:t2));
        
    % Deletion Prop
    if sum(squeeze(alpha{chr_no}(1,1,t1+1:t2+1).*beta{chr_no}(1,1,t1+1:t2+1))./c{chr_no}(t1+1:t2+1)') > del_prob_min,
        del_sum = 0;
        for sim_no = 1:no_sims,
            prob_mat = alpha{chr_no}(:,:,t1+1).*beta{chr_no}(:,:,t1+1)/c{chr_no}(t1+1);
            h = randsample([1:no_levels],1,true,sum(prob_mat,2));
            g = randsample([1:no_genos_ind(h)],1,true,prob_mat(h,1:no_genos_ind(h)));
            exit_flag = 1*(h==1);
            t = t1;
            while exit_flag == 0 & t<=t2-1,
                t=t+1;
                prob_mat = squeeze(a_trans(:,:,h,g).*b{chr_no}(:,:,t+1).*beta{chr_no}(:,:,t+1)/c{chr_no}(t)/beta{chr_no}(h,g,t));
                h_new = randsample([1:no_levels],1,true,sum(prob_mat,2));
                g_new = randsample([1:no_genos_ind(h_new)],1,true,prob_mat(h_new,1:no_genos_ind(h_new)));
                exit_flag = 1*(h_new==1);
                h = h_new;
                g = g_new;
            end
            del_sum = del_sum + exit_flag;
        end
        genes.del(gene_no,1) = del_sum/no_sims;
    else
        genes.del(gene_no,1) = 0;
    end
    
    % State change prob
    if t1==t2,
        t1 = t1 - 1;
    end
    denom2 = sum(log(c{chr_no}(t1+2:t2)));
    mat2add = log(alpha{chr_no}(:,:,t1+1)+eps)+log(eps+beta{chr_no}(:,:,t2+1));
    for u = 0:t2-t1-1,
        mat2add = mat2add+log(eps+a_diag)+log(eps+b{chr_no}(:,:,t1+u+1));
    end
    genes.change(gene_no,1) = 1 - sum(sum(squeeze(exp(mat2add+denom2))));

    
end

% Now do paired gene analyses


% Calculate probability of LOH indiced fusion


% Calculate probability of paired state changes


return
