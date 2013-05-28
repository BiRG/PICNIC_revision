
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


function [p,q,mu1,mu2,sd1,sd2,fit_score,cat] = Utils_ModalQuadFit(data_vec);


% Initialisations
warning off

p_min = 0.1;
mu_cut = 0.07;
mu_diff_min = 1.5;
fit_cut = 1.0;
block_len_min = 5000;

% Estimate parameters
p_init = 0.74;
q_init = 0.5;
mu1_init = 0.3;
mu2_init = 0.1;
sd1_init = 0.01;
sd2_init = 0.04;

params_init = [p_init,q_init,mu1_init,mu2_init,sd1_init,sd2_init]';
lb = [0.01,0.01,0.01,0.01,0.01,0.01]';
ub = [0.99,0.99,0.99,0.99,0.06,0.10]';

option_set = optimset('MaxIter',50,'MaxFunEvals',1000,'Display','off');
%option_set = optimset('MaxIter',50,'MaxFunEvals',1000,'Display','iter');
params_est = fmincon(@norm_mix,params_init,[],[],[],[],lb,ub,[],option_set,data_vec);

p = params_est(1);
q = params_est(2);
mu1 = params_est(3);
mu2 = params_est(4);
sd1 = params_est(5);
sd2 = params_est(6);

if mu2 > mu1,
    dum = mu1;
    mu1 = mu2;
    mu2 = dum;
    p = 1-p;
end
if sd2 < sd1,
    dum = sd1;
    sd1 = sd2;
    sd2 = dum;
    q = 1-q;
end

% Calculate fit
x = [0.0125:0.0250:0.9875];
y_fit = p/2*q*normpdf(x,0.5+mu1,sd1)+p/2*q*normpdf(x,0.5-mu1,sd1)+(1-p)/2*q*normpdf(x,0.5+mu2,sd1)+(1-p)/2*q*normpdf(x,0.5-mu2,sd1)+...
        p/2*(1-q)*normpdf(x,0.5+mu1,sd2)+p/2*(1-q)*normpdf(x,0.5-mu1,sd2)+(1-p)/2*(1-q)*normpdf(x,0.5+mu2,sd2)+(1-p)/2*(1-q)*normpdf(x,0.5-mu2,sd2);
bin_heights = hist(data_vec,x);
bin_heights = bin_heights/sum(bin_heights)*40;
fit_score = sum(max(0,y_fit-bin_heights).^2);

% Calculate genotype class
if p*(1-p)>p_min*(1-p_min),
    if mu1<=mu_cut & mu2<=mu_cut,
        cat = 1;
    elseif mu1>mu_cut & mu2<=mu_cut & abs(mu1-mu2)/(sd1+sd2)>mu_diff_min,
        cat = 3;
    elseif mu1>mu_cut & mu2>mu_cut & abs(mu1-mu2)/(sd1+sd2)>mu_diff_min,
        cat = 4;
    elseif mu1>mu_cut & mu2>mu_cut & abs(mu1-mu2)/(sd1+sd2)<=mu_diff_min,
        cat = 2;
    else
        cat = -1;
    end
else
    if mu1<=mu_cut,
        cat = 1;
    elseif mu1>mu_cut,
        cat = 2;
    else
        cat = -1;
    end
end

% Check for bad fit
if fit_score > fit_cut | length(data_vec)<block_len_min,
    cat = -1;
end

return

function neg_ll_val = norm_mix(params,x);


p = params(1);
q = params(2);
mu1 = params(3);
mu2 = params(4);
sd1 = params(5);
sd2 = params(6);

neg_ll_val = -sum(log(p/2*q*normpdf(x,0.5+mu1,sd1)+p/2*q*normpdf(x,0.5-mu1,sd1)+(1-p)/2*q*normpdf(x,0.5+mu2,sd1)+(1-p)/2*q*normpdf(x,0.5-mu2,sd1)+...
    p/2*(1-q)*normpdf(x,0.5+mu1,sd2)+p/2*(1-q)*normpdf(x,0.5-mu1,sd2)+(1-p)/2*(1-q)*normpdf(x,0.5+mu2,sd2)+(1-p)/2*(1-q)*normpdf(x,0.5-mu2,sd2)));

return
