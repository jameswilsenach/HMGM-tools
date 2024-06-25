function IR = getIR(trans)
%% getIR.m
% Calculates the information (or entropy) rate of the HMM state space model.
%
%% Input:
%
% trans : nxn stochastic matrix (rows will be normalised if not stochastic)
%
%% Output:
%
% IR : information gained (in bits) by sampling from the HMM at equilibtrium
%%
if ~all(sum(trans,2)==1)
    trans = trans./sum(trans,2);
end

if ~all(sum(trans,2)==1)
    error('The transition matrix cannot contain a row with all zero entries or any negative entries')
end
mc = dtmc(trans);
if ~isergodic(mc)
    error('This Markov process is not ergodic and so has no unique stationary distribution')
end

[~,~,W] = eig(trans); %Get the stationary distribution from the eigenvectors
stat = real(W(:,1))';
stat = stat/sum(stat); %Normalise the eigenvector so that it sums to one
IR = -sum(diag(stat)*trans.*log(trans),'all');%calculate the IR
end