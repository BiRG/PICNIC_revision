
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


function HMM_BaumWelch(varargin)

optargin = size(varargin,2);



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

% Structure conversions
chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;
no_levels = seg_info.no_levels;
no_genos_ind = seg_info.no_genos_ind;
a0 = params.a0;
chi = params.chi;
init = params.init;
delta = params.delta;
gamma = params.gamma;
zeta = params.zeta;
q = params.q;
r = params.r;
eta = params.eta;
sd1 = params.sd1;
sd2 = params.sd2;

% Initialisations
warning off
if(optargin==0)
	no_loops_max = 0;
else
	no_loops_max=str2num(varargin{1});
end
eps = 1e-300;
a(no_levels,no_genos_ind(no_levels),no_levels,no_genos_ind(no_levels)) = 0;

% Prior parameters
for level_no1 = 1:no_levels,
    for geno_no1 = 1:no_genos_ind(level_no1),
        ar0(level_no1,geno_no1) = 1;
        ind0(level_no1,geno_no1) = 1;
    end
end
ep1 = 0.74*1000;
ep2 = 0.26*1000;
eq1 = 0.7*10;
eq2 = 0.3*10;
er1 = 0.995*1000;
er2 = 0.005*1000;
echi_change = 1e-18;
echi_change2del = 1e-14;
echi_change4del = 1e-8;
for level_no1 = 1:no_levels,
    for geno_no1 = 1:no_genos_ind(level_no1),
        for level_no2 = 1:no_levels,
            for geno_no2 = 1:no_genos_ind(level_no2),
                if level_no1==level_no2 & geno_no1==geno_no2,
                    echi(level_no1,geno_no1,level_no2,geno_no2) = 1;
                else
                    if level_no1==1,
                        echi(level_no1,geno_no1,level_no2,geno_no2) = echi_change4del;
                    else
                        if level_no2==1,
                            echi(level_no1,geno_no1,level_no2,geno_no2) = echi_change2del;
                        else
                            echi(level_no1,geno_no1,level_no2,geno_no2) = echi_change;
                        end
                    end
                end
            end
        end
        echi(level_no1,geno_no1,:,:) = echi(level_no1,geno_no1,:,:)/sum(sum(echi(level_no1,geno_no1,:,:)));
    end
end
echi = echi*10*1e18;
B_init = 1.5*0.1;          % 1.5 is COV=sd/mn, other value is sd - widen for less stringency
B_delta = 1.5*0.1;
B_gamma = 1.5*0.1;
B_zeta = 1.5*1;
B_eta = 1.5*0.20;
B_sd1 = 1.5*0.05;
B_sd2 = 1.5*0.05;
A_init = 0.5/B_init;       % numerator is mean
A_delta = 0.2/B_delta;
A_gamma = 0.2/B_gamma;
A_zeta = 5/B_zeta;
A_eta = 0.60/B_eta;
A_sd1 = 0.08/B_sd1;
A_sd2 = 0.12/B_sd2;

% Initialisations
like_iter_lim = 1e-6;
like_delta = inf;
like_0 = -inf;
loop_no = 0;
no_chr = length(chr_start);
no_states = sum(no_genos_ind);
Tc = chr_end-chr_start+1;
for chr_no = 1:no_chr,
    b{chr_no}(1:no_levels,1:no_genos_ind(no_levels),1:Tc(chr_no)) = 0;
    alpha{chr_no}(1:no_levels,1:no_genos_ind(no_levels),1:Tc(chr_no)) = 0;
    beta{chr_no}(1:no_levels,1:no_genos_ind(no_levels),1:Tc(chr_no)) = 0;
    c{chr_no}(1:Tc(chr_no)) = 0;
end

% loop through BaumWelch iterations
while loop_no < no_loops_max & abs(like_delta) > like_iter_lim,
    
    % Initialisations
    loop_no = loop_no + 1
    a_trans = permute(a,[3,4,1,2]);
    qq(1) = q;
    qq(2) = 1-q;
    rr(1) = r;
    rr(2) = 1-r;
    sd(1) = sd1;
    sd(2) = sd2;

    % Calculate emissions
    for level_no  = 1:no_levels,
        for geno_no = 1:no_genos_ind(level_no),
            for chr_no = 1:no_chr,
                NM = NormMix(log(H{chr_no}(1:Tc(chr_no))),log(init+(level_no-1)*delta),gamma,r,zeta);
                MF = MixtureFull(G{chr_no}(1:Tc(chr_no)),p(chr_start(chr_no):chr_end(chr_no),:,:),q,eta,sd1,sd2,level_no,geno_no);
                MF = MF.*SNP_ind{chr_no}+(1-SNP_ind{chr_no});
                b{chr_no}(level_no,geno_no,1:Tc(chr_no)) = NM.*MF;
            end
        end
    end
    
    % Construct Transition matrix
    for level_no1 = 1:no_levels,
        for geno_no1 = 1:no_genos_ind(level_no1),
            for level_no2 = 1:no_levels,
                for geno_no2 = 1:no_genos_ind(level_no1),
                    a(level_no1,geno_no1,level_no2,geno_no2) = chi(level_no1,geno_no1,level_no2,geno_no2);
                    a_trans(level_no2,geno_no2,level_no1,geno_no1) = chi(level_no1,geno_no1,level_no2,geno_no2);
                end
            end
        end
    end
    
    % Calculate forward and backward equations for alpha and beta
    for chr_no = 1:no_chr,
        alpha{chr_no}(:,:,1) = a0.*b{chr_no}(:,:,1);
        c{chr_no}(1) = 1/sum(sum(alpha{chr_no}(:,:,1)));
        alpha{chr_no}(:,:,1) = alpha{chr_no}(:,:,1)*c{chr_no}(1);
        for t = 1:Tc(chr_no),
            alpha{chr_no}(:,:,t+1) = squeeze(sum(sum(repmat(squeeze(alpha{chr_no}(:,:,t)),[1,1,no_levels,no_genos_ind(no_levels)]).*a,1),2));
            alpha{chr_no}(:,:,t+1) = b{chr_no}(:,:,t).*alpha{chr_no}(:,:,t+1);
            c{chr_no}(t+1) = 1/sum(sum(alpha{chr_no}(:,:,t+1)));
            alpha{chr_no}(:,:,t+1) = alpha{chr_no}(:,:,t+1)*c{chr_no}(t+1);
        end
        beta{chr_no}(:,:,Tc(chr_no)+1) = ind0*c{chr_no}(Tc(chr_no)+1);
        for t = Tc(chr_no):-1:1,
            beta_tens = repmat(squeeze(beta{chr_no}(:,:,t+1)),[1,1,no_levels,no_genos_ind(no_levels)]);
            emiss_out = repmat(b{chr_no}(:,:,t),[1,1,no_levels,no_genos_ind(no_levels)]);
            beta{chr_no}(:,:,t) = c{chr_no}(t)*squeeze(sum(sum(beta_tens.*a_trans.*emiss_out,1),2));
        end
    end
    
    % Update initial distribution params
    init_mat = alpha{1}(:,:,1).*beta{1}(:,:,1);
    for chr_no = 2:no_chr,
        init_mat = init_mat + alpha{chr_no}(:,:,1).*beta{chr_no}(:,:,1);
    end
    a0_new = ind0.*(init_mat + (ar0 - 1))/(no_chr+sum(sum(ar0))-no_states);
    
    % Update transition matrix parameter
    chi_new(1:no_levels,1:no_genos_ind(no_levels),1:no_levels,1:no_genos_ind(no_levels)) = 0;
    for level_no1 = 1:no_levels,
        for geno_no1 = 1:no_genos_ind(level_no1),
            for level_no2 = 1:no_levels,
                for geno_no2 = 1:no_genos_ind(level_no2),
                    for chr_no = 1:no_chr,
                        chi_new(level_no1,geno_no1,level_no2,geno_no2) = chi_new(level_no1,geno_no1,level_no2,geno_no2) +...
                                                                         sum(alpha{chr_no}(level_no1,geno_no1,1:Tc(chr_no)).*...
                                                                         a(level_no1,geno_no1,level_no1,geno_no1).*...
                                                                         b{chr_no}(level_no1,geno_no1,1:Tc(chr_no)).*...
                                                                         beta{chr_no}(level_no1,geno_no1,2:Tc(chr_no)+1),3);
                    end
                end
            end
            chi_new(level_no1,geno_no1,:,:) = chi(level_no1,geno_no1,:,:).*chi_new(level_no1,geno_no1,:,:) + echi(level_no1,geno_no1,:,:) - 1;
            chi_new(level_no1,geno_no1,:,:) = chi_new(level_no1,geno_no1,:,:)/sum(sum(chi_new(level_no1,geno_no1,:,:)));
        end
    end
    
    % Calculate emissions params updates
    q1 = eq1 - 1;
    q2 = eq2 - 1;
    r1 = er1 - 1;
    r2 = er2 - 1;
    X(1:no_levels,1:no_genos_ind(no_levels)) = 0;
    Y(1:no_levels,1:no_genos_ind(no_levels)) = 0;
    Z(1:no_levels,1:no_genos_ind(no_levels)) = 0;
    XX(1:no_levels,1:no_genos_ind(no_levels),1:2,1:2,1:2) = 0;
    YY(1:no_levels,1:no_genos_ind(no_levels),1:2,1:2,1:2) = 0;
    ZZ(1:no_levels,1:no_genos_ind(no_levels),1:2,1:2,1:2) = 0;
    X_zeta(1:no_levels,1:no_genos_ind(no_levels)) = 0;
    Y_zeta(1:no_levels,1:no_genos_ind(no_levels)) = 0;
    for level_no  = 1:no_levels,
        level_no;
        mu(1) = 2/pi*atan(1+(level_no-1)*eta)-0.5;
        for geno_no = 1:no_genos_ind(level_no),
            geno_no;
            mu(2) = (1/2-max(0,(geno_no-1)/(level_no-1)))*(4/pi*atan(1+(level_no-1)*eta)-1);
            for chr_no = 1:no_chr,
                chr_no;
                LH = log(H{chr_no}(1:Tc(chr_no)));
                NM = NormMix(LH,log(init+(level_no-1)*delta),gamma,r,zeta);
                MF = MixtureFull(G{chr_no}(1:Tc(chr_no)),p(chr_start(chr_no):chr_end(chr_no),:,:),q,eta,sd1,sd2,level_no,geno_no);
                MF = MF.*SNP_ind{chr_no} + (1 - SNP_ind{chr_no});
                MJ1 = MixtureJ(G{chr_no}(1:Tc(chr_no)),p(chr_start(chr_no):chr_end(chr_no),:,:),q,eta,sd1,sd2,level_no,geno_no,1);
                MJ1 = MJ1.*SNP_ind{chr_no} + q*(1 - SNP_ind{chr_no});
                MJ2 = MF - MJ1;
                NM2 = (1-r)*normpdf(LH,0,zeta);
                NM1 = NM - NM2;
                q1 = q1 + sum(SNP_ind{chr_no}.*...
                              squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                              squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                              c{chr_no}(2:Tc(chr_no)+1)'./...
                              MF.*MJ1);
                q2 = q2 + sum(SNP_ind{chr_no}.*...
                              squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                              squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                              c{chr_no}(2:Tc(chr_no)+1)'./...
                              MF.*MJ2);
                r1 = r1 + sum(squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                              squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                              c{chr_no}(2:Tc(chr_no)+1)'./...
                              NM.*NM1);
                r2 = r2 + sum(squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                              squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                              c{chr_no}(2:Tc(chr_no)+1)'./...
                              NM.*NM2);
                X(level_no,geno_no) = X(level_no,geno_no) + sum(squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                                                                c{chr_no}(2:Tc(chr_no)+1)'.*...
                                                                NM1./NM.*...
                                                                (LH).^2);
                Y(level_no,geno_no) = Y(level_no,geno_no) + sum(squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                                                                c{chr_no}(2:Tc(chr_no)+1)'.*...
                                                                NM1./NM.*...
                                                                LH);
                Z(level_no,geno_no) = Z(level_no,geno_no) + sum(squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                                                                c{chr_no}(2:Tc(chr_no)+1)'.*...
                                                                NM1./NM);
                X_zeta(level_no,geno_no) = X_zeta(level_no,geno_no) + sum(squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                          squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                                                                          c{chr_no}(2:Tc(chr_no)+1)'.*...
                                                                          NM2./NM.*...
                                                                          (LH).^2);
                Y_zeta(level_no,geno_no) = Y_zeta(level_no,geno_no) + sum(squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                          squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1))./...
                                                                          c{chr_no}(2:Tc(chr_no)+1)'.*...
                                                                          NM2./NM);
                for ii = 1:2,
                    for jj = 1:2,
                        for kk = 1:2,
                            sgn_tens(level_no,geno_no,ii,jj,kk) = (-1)^(kk+1);
                            NP = normpdf(G{chr_no}(1:Tc(chr_no)),(0.5+(-1)^(kk+1)*mu(ii)),sd(jj));
                            ZZ(level_no,geno_no,ii,jj,kk) = ZZ(level_no,geno_no,ii,jj,kk) + sum(SNP_ind{chr_no}.*...
                                                                              squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                              squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                              (p(chr_start(chr_no):chr_end(chr_no),ii,kk)*qq(jj).*NP)./...
                                                                              MF./...
                                                                              c{chr_no}(2:Tc(chr_no)+1)'.*...
                                                                              (G{chr_no}(1:Tc(chr_no))).^2);
                            YY(level_no,geno_no,ii,jj,kk) = YY(level_no,geno_no,ii,jj,kk) + sum(SNP_ind{chr_no}.*...
                                                                              squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                              squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                              (p(chr_start(chr_no):chr_end(chr_no),ii,kk)*qq(jj).*NP)./...
                                                                              MF./...
                                                                              c{chr_no}(2:Tc(chr_no)+1)'.*...
                                                                              G{chr_no}(1:Tc(chr_no)));
                            XX(level_no,geno_no,ii,jj,kk) = XX(level_no,geno_no,ii,jj,kk) + sum(SNP_ind{chr_no}.*...
                                                                              squeeze(alpha{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                              squeeze(beta{chr_no}(level_no,geno_no,2:Tc(chr_no)+1)).*...
                                                                              (p(chr_start(chr_no):chr_end(chr_no),ii,kk)*qq(jj).*NP)./...
                                                                              MF./...
                                                                              c{chr_no}(2:Tc(chr_no)+1)');
                        end
                    end
                end
            end
        end
    end
    q_new = q1/(q1+q2);
    r_new = r1/(r1+r2);
    cnvs1 = [init,delta,gamma];
    option_set = optimset('MaxIter',50,'MaxFunEvals',1000);
    params_est = fsolve(@cnv_fn,cnvs1,option_set,X,Y,Z,ind0,no_levels,no_genos_ind,A_init,B_init,A_delta,B_delta,A_gamma,B_gamma);
    init_new = params_est(1);
    delta_new = params_est(2);
    gamma_new = params_est(3);
    cnvs2 = [zeta];
    option_set = optimset('MaxIter',50,'MaxFunEvals',1000);
    params_est = fsolve(@cnv_fn2,cnvs2,option_set,X_zeta,Y_zeta,ind0,no_levels,no_genos_ind,A_zeta,B_zeta);
    zeta_new = params_est(1);
    genos1 = [eta,sd1,sd2];
    option_set = optimset('MaxIter',50,'MaxFunEvals',1000);
    params_est = fsolve(@geno_fn,genos1,option_set,XX,YY,ZZ,no_levels,no_genos_ind,A_eta,B_eta,A_sd1,B_sd1,A_sd2,B_sd2,sgn_tens);
    eta_new = params_est(1);
    sd1_new = params_est(2);
    sd2_new = params_est(3);
    %'Calculated Emission Params'
    
    disp('update likelihood');
    % Calculate updated likelihood and difference
    LL = 0;
    for chr_no = 1:no_chr,
        LL = LL - sum(log(c{chr_no}));
    end
    LL = LL + (A_init-1)*log(init)-init/B_init;
    LL = LL + (A_delta-1)*log(delta)-delta/B_delta;
    LL = LL + (A_gamma-1)*log(gamma)-gamma/B_gamma;
    LL = LL + (A_zeta-1)*log(zeta)-zeta/B_zeta;
    LL = LL + (A_eta-1)*log(eta)-eta/B_eta;
    LL = LL + (A_sd1-1)*log(sd1)-sd1/B_sd1;
    LL = LL + (A_sd2-1)*log(sd2)-sd2/B_sd2;
    LL = LL + sum(sum(sum(sum((echi-1).*log(chi)))));
    LL = LL + (ep1-1)*log(p) + (ep2-1)*log(1-p);
    LL = LL + (eq1-1)*log(q) + (eq2-1)*log(1-q);
    LL = LL + (er1-1)*log(r) + (er2-1)*log(1-r);
    for level_no1 = 1:no_levels,
        for geno_no1 = 1:no_genos_ind(level_no1),
            LL = LL + (ar0(level_no1,geno_no1)-1).*log(eps+a0(level_no1,geno_no1));
        end
    end
    like_delta = LL - like_0;
    %'Likelihood Done'
    
    % Renew parameters
    like_0 = LL;
    a0 = a0_new;
    chi = chi_new;
    q = q_new;
    r = r_new;
    init = init_new;
    delta = delta_new;
    gamma = gamma_new;
    zeta = zeta_new;
    eta = eta_new;
    sd1 = sd1_new;
    sd2 = sd2_new;
    %'Params renewed'
    
end

% Do emission probs calcs for last time
for level_no  = 1:no_levels,
    for geno_no = 1:no_genos_ind(level_no),
        for chr_no = 1:no_chr,
            NM = NormMix(log(H{chr_no}(1:Tc(chr_no))),log(init+(level_no-1)*delta),gamma,r,zeta);
            MF = MixtureFull(G{chr_no}(1:Tc(chr_no)),p(chr_start(chr_no):chr_end(chr_no),:,:),q,eta,sd1,sd2,level_no,geno_no);
            MF = MF.*SNP_ind{chr_no}+(1-SNP_ind{chr_no});
            b{chr_no}(level_no,geno_no,1:Tc(chr_no)) = NM.*MF;
        end
    end
end

% Construct Transition matrix
for level_no1 = 1:no_levels
    for geno_no1 = 1:no_genos_ind(level_no1),
        for level_no2 = 1:no_levels
            for geno_no2 = 1:no_genos_ind(level_no1),
                a(level_no1,geno_no1,level_no2,geno_no2) = chi(level_no1,geno_no1,level_no2,geno_no2);
                a_trans(level_no2,geno_no2,level_no1,geno_no1) = chi(level_no1,geno_no1,level_no2,geno_no2);
            end
        end
    end
end

% Do alphas & Betas for last time
for chr_no = 1:no_chr,
    chr_no
    alpha{chr_no}(:,:,1) = a0.*b{chr_no}(:,:,1);
    c{chr_no}(1) = 1/sum(sum(alpha{chr_no}(:,:,1)));
    alpha{chr_no}(:,:,1) = alpha{chr_no}(:,:,1)*c{chr_no}(1);
    for t = 1:Tc(chr_no),
        [chr_no,t];
        alpha{chr_no}(:,:,t+1) = squeeze(sum(sum(repmat(squeeze(alpha{chr_no}(:,:,t)),[1,1,no_levels,no_genos_ind(no_levels)]).*a,1),2));
        alpha{chr_no}(:,:,t+1) = b{chr_no}(:,:,t).*alpha{chr_no}(:,:,t+1);
        c{chr_no}(t+1) = 1/sum(sum(alpha{chr_no}(:,:,t+1)));
        alpha{chr_no}(:,:,t+1) = alpha{chr_no}(:,:,t+1)*c{chr_no}(t+1);
    end
    beta{chr_no}(:,:,Tc(chr_no)+1) = ind0*c{chr_no}(Tc(chr_no)+1);
    for t = Tc(chr_no):-1:1,
        [chr_no,t];
        beta_tens = repmat(squeeze(beta{chr_no}(:,:,t+1)),[1,1,no_levels,no_genos_ind(no_levels)]);
        emiss_out = repmat(b{chr_no}(:,:,t),[1,1,no_levels,no_genos_ind(no_levels)]);
        beta{chr_no}(:,:,t) = c{chr_no}(t)*squeeze(sum(sum(beta_tens.*a_trans.*emiss_out,1),2));
    end
end
%'Alpha, Betas finalised'

% Finalise parameters
params.a0 = a0;
params.init = init;
params.delta = delta;
params.gamma = gamma;
params.zeta = zeta;
params.chi = chi;
params.q = q;
params.r = r;
params.eta = eta;
params.sd1 = sd1;
params.sd2 = sd2;

return

% ----------------------------------------------------------------
% Kernel functions
% ----------------------------------------------------------------

function val = NormMix(vec,mn,sd,r,zeta)


val = r*normpdf(vec,mn,sd)+(1-r)*normpdf(vec,0,zeta);

return

function val = MixtureFull(vec,p,q,eta,sd1,sd2,level_no,geno_no)


mu1 = 2/pi*atan(1+(level_no-1)*eta)-0.5;
mu2 = (1/2-max(0,(geno_no-1)/(level_no-1)))*(4/pi*atan(1+(level_no-1)*eta)-1);
val = p(:,1,1)*q.*normpdf(vec,0.5+mu1,sd1)+p(:,1,2)*q.*normpdf(vec,0.5-mu1,sd1)+...
      p(:,1,1)*(1-q).*normpdf(vec,0.5+mu1,sd2)+p(:,1,2)*(1-q).*normpdf(vec,0.5-mu1,sd2)+...
      p(:,2,1)*q.*normpdf(vec,0.5+mu2,sd1)+p(:,2,2)*q.*normpdf(vec,0.5-mu2,sd1)+...
      p(:,2,1)*(1-q).*normpdf(vec,0.5+mu2,sd2)+p(:,2,2)*(1-q).*normpdf(vec,0.5-mu2,sd2);

return

function val = MixtureJ(vec,p,q,eta,sd1,sd2,level_no,geno_no,j)


mu1 = 2/pi*atan(1+(level_no-1)*eta)-0.5;
mu2 = (1/2-max(0,(geno_no-1)/(level_no-1)))*(4/pi*atan(1+(level_no-1)*eta)-1);
q_j = (j==1)*q+(j==2)*(1-q);
sd = (j==1)*sd1+(j==2)*sd2;
val = p(:,1,1)*q_j.*normpdf(vec,0.5+mu1,sd)+p(:,1,2)*q_j.*normpdf(vec,0.5-mu1,sd)+...
      p(:,2,1)*q_j.*normpdf(vec,0.5+mu2,sd)+p(:,2,2)*q_j.*normpdf(vec,0.5-mu2,sd);
      
return

function val = cnv_fn(params,X,Y,Z,ind0,no_levels,no_genos_ind,A_init,B_init,A_delta,B_delta,A_gamma,B_gamma)


init = params(1);
delta = params(2);
gamma = params(3);

mat = repmat([0:no_levels-1]',[1,no_genos_ind(no_levels)]);

val1 = sum(sum(ind0.*(Y./(init+mat*delta)-Z.*log(init+mat*delta)./(init+mat*delta))))+gamma^2*((A_init-1)/init-1/B_init);
val2 = sum(sum(ind0.*(Y.*(mat)./(init+mat*delta)-Z.*mat.*log(init+mat*delta)./(init+mat*delta))))+gamma^2*((A_delta-1)/delta-1/B_delta);
val3 = sum(sum(ind0.*(X-2*Y.*log(init+mat*delta)+Z.*(-gamma^2+log(init+mat*delta).^2))))+gamma^2*((A_gamma-1)-gamma/B_gamma);

val = [val1;val2;val3];

return

function val = cnv_fn2(params,X_zeta,Y_zeta,ind0,no_levels,no_genos_ind,A_zeta,B_zeta)


zeta = params(1);

val = zeta^3+zeta^2*B_zeta*(sum(sum(Y_zeta))-(A_zeta-1))-B_zeta*sum(sum(X_zeta));

return

function val = geno_fn(params,XX,YY,ZZ,no_levels,no_genos_ind,A_eta,B_eta,A_sd1,B_sd1,A_sd2,B_sd2,sgn)


eta = params(1);
sd1 = params(2);
sd2 = params(3);

mu = sgn*0;
mu_diff = sgn*0;

for level_no = 1:no_levels,
    for geno_no = 1:no_genos_ind(level_no),
        mu0(1) = 2/pi*atan(1+(level_no-1)*eta)-0.5;
        mu0(2) = (1/2-max(0,(geno_no-1)/(level_no-1)))*(4/pi*atan(1+(level_no-1)*eta)-1);
        mu_diff0(1) = 2/pi*(level_no-1)/(1+(1+(level_no-1)*eta)^2);
        mu_diff0(2) = 1/pi*(2*(level_no-1)-4*(geno_no-1))/(1+(1+(level_no-1)*eta)^2);
        for ii = 1:2,
            for jj = 1:2,
                for kk = 1:2,
                    mu(level_no,geno_no,ii,jj,kk) = mu0(ii);
                    mu_diff(level_no,geno_no,ii,jj,kk) = mu_diff0(ii);
                end
            end
        end
    end
end

sd = mu*0;
sd(:,:,:,1,:) = sd1;
sd(:,:,:,2,:) = sd2;

val1 = sum(sum(sum(sum(sum(sgn./sd.^2.*mu_diff.*(YY-(0.5+sgn.*mu).*XX))))))+((A_eta-1)/eta-1/B_eta);
val2 = sum(sum(sum(sum(-XX(:,:,:,1,:).*sd(:,:,:,1,:).^2+ZZ(:,:,:,1,:)-2*YY(:,:,:,1,:).*(0.5+sgn(:,:,:,1,:).*mu(:,:,:,1,:))+...
                       XX(:,:,:,1,:).*(0.5+sgn(:,:,:,1,:).*mu(:,:,:,1,:)).^2))))+...
                       sd1^3*((A_sd1-1)/sd1-1/B_sd1);
val3 = sum(sum(sum(sum(-XX(:,:,:,2,:).*sd(:,:,:,2,:).^2+ZZ(:,:,:,2,:)-2*YY(:,:,:,2,:).*(0.5+sgn(:,:,:,2,:).*mu(:,:,:,2,:))+...
                       XX(:,:,:,2,:).*(0.5+sgn(:,:,:,2,:).*mu(:,:,:,2,:)).^2))))+...
                       sd2^3*((A_sd2-1)/sd2-1/B_sd2);

val = [val1;val2;val3];

return
