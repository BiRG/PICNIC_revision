
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


function HMM_GenoType


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

% Initialisations
chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;
no_chr = seg_info.no_chr;
a0 = params.a0;
eta = params.eta;
q = params.q;
sd1 = params.sd1;
sd2 = params.sd2;
Tc = seg_info.Tc;
no_levels = seg_info.no_levels;
no_genos_ind = seg_info.no_genos_ind;
for level_no = 1:no_levels,
    for geno_no = 1:no_genos_ind(level_no),
        a_diag(level_no,geno_no) = params.chi(level_no,geno_no,level_no,geno_no);
    end
end
for chr_no = 1:no_chr,
    chr_no;
    mix_full(1:Tc(chr_no),1:no_levels,1:no_genos_ind(no_levels)) = 0;
    mix_IK(1:Tc(chr_no),1:no_levels,1:no_genos_ind(no_levels),1:2,1:2) = 0;
    geno_like(1:Tc(chr_no),1:no_levels,no_genos_ind(no_levels),1:2,1:2) = 0;
    geno_like2max(1:Tc(chr_no),1:no_levels,no_genos_ind(no_levels),1:2,1:2) = 0;
    geno_class{chr_no}(1:Tc(chr_no),1:2) = 0;
    geno_prob{chr_no}(1:Tc(chr_no),1) = 0;
    geno_cond_prob{chr_no}(1:Tc(chr_no),1) = 0;
    geno_HET_prob{chr_no}(1:Tc(chr_no),1) = 0;
    geno_A_LOH_prob{chr_no}(1:Tc(chr_no),1) = 0;
    geno_B_LOH_prob{chr_no}(1:Tc(chr_no),1) = 0;
    geno_state_change{chr_no}(1:Tc(chr_no),1) = 0;
    for level_no = 1:no_levels,
        for geno_no = 1:no_genos_ind(level_no),
            mix_full(1:Tc(chr_no),level_no,geno_no) = MixtureFull(G{chr_no},p(chr_start(chr_no):chr_end(chr_no),:,:),q,eta,sd1,sd2,level_no,geno_no)';
            mix_full(1:Tc(chr_no),level_no,geno_no) = mix_full(1:Tc(chr_no),level_no,geno_no).*SNP_ind{chr_no} + (1 - SNP_ind{chr_no});
            for ii = 1:2,
                for kk = 1:2,
                    mix_IK(1:Tc(chr_no),level_no,geno_no,ii,kk) = MixtureIK(G{chr_no},p(chr_start(chr_no):chr_end(chr_no),:,:),q,eta,sd1,sd2,level_no,geno_no,ii,kk)';
                    mix_IK(1:Tc(chr_no),level_no,geno_no,ii,kk) = mix_IK(1:Tc(chr_no),level_no,geno_no,ii,kk).*SNP_ind{chr_no} + 1*(1 - SNP_ind{chr_no});
                    geno_like(1:Tc(chr_no),level_no,geno_no,ii,kk) = mix_IK(1:Tc(chr_no),level_no,geno_no,ii,kk)./...
                                                                     mix_full(1:Tc(chr_no),level_no,geno_no).*...
                                                                     squeeze(alpha{chr_no}(level_no,geno_no,2:1+Tc(chr_no))).*...
                                                                     squeeze(beta{chr_no}(level_no,geno_no,2:1+Tc(chr_no)))./...
                                                                     c{chr_no}(2:Tc(chr_no)+1)';
                    geno_like2max(1:Tc(chr_no),level_no,geno_no,ii,kk) = geno_like(1:Tc(chr_no),level_no,geno_no,ii,kk);
                end
            end
            % Combine for less than four clusters
            if level_no == 1 & geno_no == 1,                                                % Deletions
                xx = sum(sum(geno_like(1:Tc(chr_no),level_no,geno_no,:,:),4),5);
                for ii = 1:2,
                    for kk = 1:2,
                        geno_like2max(1:Tc(chr_no),level_no,geno_no,ii,kk) = xx;
                    end
                end
            end
            if level_no > 1 & geno_no == 1,                                                 % LOH
                for kk = 1:2,
                    xx = sum(geno_like(1:Tc(chr_no),level_no,geno_no,:,kk),4);
                    for ii = 1:2,
                        geno_like2max(1:Tc(chr_no),level_no,geno_no,ii,kk) = xx;
                    end
                end
            end
            if level_no > 1 & geno_no == (level_no+1)/2,                                    % Even
                xx = sum(geno_like(1:Tc(chr_no),level_no,geno_no,2,:),5); 
                for kk = 1:2,       
                    geno_like2max(1:Tc(chr_no),level_no,geno_no,2,kk) = xx;
                end
            end
        end
    end
    [m4,i4] = max(geno_like2max(1:Tc(chr_no),:,:,:,:),[],5);
    [m3,i3] = max(m4,[],4);
    [m2,i2] = max(m3,[],3);
    [m1,i1] = max(m2,[],2);
    level_ind = i1;
    geno_ind = m1*0;
    pk_ind = geno_ind;
    all_ind = geno_ind;
    for t = 1:Tc(chr_no),
        geno_ind(t) = i2(t,i1(t));
        pk_ind(t) = i3(t,level_ind(t),geno_ind(t));
        all_ind(t) = i4(t,level_ind(t),geno_ind(t),pk_ind(t));
        if geno_ind(t)==1,
            pk_vec = [1:2];
        else
            pk_vec = pk_ind(t);
        end
        if 2*geno_ind(t)==level_ind(t)+1 & pk_ind(t)==2,
            all_vec = [1:2];
        elseif level_ind(t)==1,
            all_vec = [1:2];
        else
            all_vec = all_ind(t);
        end
        geno_class{chr_no}(t,1) = 1*(pk_ind(t)==1)*(all_ind(t)==1)*(level_ind(t)-1)+...
                                  1*(pk_ind(t)==2)*(all_ind(t)==1)*(level_ind(t)-geno_ind(t))+...
                                  1*(pk_ind(t)==2)*(all_ind(t)==2)*(geno_ind(t)-1)+...
                                  1*(pk_ind(t)==1)*(all_ind(t)==2)*(0);
        geno_class{chr_no}(t,2) = level_ind(t)-1-geno_class{chr_no}(t,1);
        geno_prob{chr_no}(t,1) = sum(sum(squeeze(geno_like(t,level_ind(t),geno_ind(t),pk_vec,all_vec))));
        geno_cond_prob{chr_no}(t,1) = geno_prob{chr_no}(t,1)/squeeze(sum(sum(geno_like(t,level_ind(t),geno_ind(t),:,:))));
        geno_HET_prob{chr_no}(t,1) = sum(sum(sum(squeeze(geno_like(t,2:no_levels,2:no_genos_ind(no_levels),2,1:2)))));
        geno_A_LOH_prob{chr_no}(t,1) = sum(sum(squeeze(geno_like(t,2:no_levels,1,1:2,1))))+sum(sum(squeeze(geno_like(t,2:no_levels,2:no_genos_ind(no_levels),1,1))));
        geno_B_LOH_prob{chr_no}(t,1) = sum(sum(squeeze(geno_like(t,2:no_levels,1,1:2,2))))+sum(sum(squeeze(geno_like(t,2:no_levels,2:no_genos_ind(no_levels),1,2))));
        geno_state_change{chr_no}(t,1) = 1 - sum(sum(alpha{chr_no}(:,:,t).*a_diag.*b{chr_no}(:,:,t).*beta{chr_no}(:,:,t+1)));
    end
end

% Assign outputs
genotypes.geno_class = geno_class;
genotypes.geno_prob = geno_prob;
genotypes.geno_cond_prob = geno_cond_prob;
genotypes.geno_HET_prob = geno_HET_prob;
genotypes.geno_A_LOH_prob = geno_A_LOH_prob;
genotypes.geno_B_LOH_prob = geno_B_LOH_prob;
genotypes.geno_state_change = geno_state_change;

return

function val = MixtureFull(vec,p,q,eta,sd1,sd2,level_no,geno_no)


mu1 = 2/pi*atan(1+(level_no-1)*eta)-0.5;
mu2 = (1/2-max(0,(geno_no-1)/(level_no-1)))*(4/pi*atan(1+(level_no-1)*eta)-1);
val = p(:,1,1)*q.*normpdf(vec,0.5+mu1,sd1)+p(:,1,2)*q.*normpdf(vec,0.5-mu1,sd1)+...
      p(:,1,1)*(1-q).*normpdf(vec,0.5+mu1,sd2)+p(:,1,2)*(1-q).*normpdf(vec,0.5-mu1,sd2)+...
      p(:,2,1)*q.*normpdf(vec,0.5+mu2,sd1)+p(:,2,2)*q.*normpdf(vec,0.5-mu2,sd1)+...
      p(:,2,1)*(1-q).*normpdf(vec,0.5+mu2,sd2)+p(:,2,2)*(1-q).*normpdf(vec,0.5-mu2,sd2);

return
  
function val = MixtureIK(vec,p,q,eta,sd1,sd2,level_no,geno_no,i,k)


mu = (i==1)*(2/pi*atan(1+(level_no-1)*eta)-0.5)+(i==2)*(1/2-max(0,(geno_no-1)/(level_no-1)))*(4/pi*atan(1+(level_no-1)*eta)-1);
sgn = (k==1)*1-1*(k==2);
val = p(:,i,k)*q.*normpdf(vec,0.5+sgn*mu,sd1)+p(:,i,k)*(1-q).*normpdf(vec,0.5+sgn*mu,sd2);

return
