

function good_or_bad = test_function(W_standard,W_to_be_compared)

% The following code is prepared by The Laboratory of Neuroscience for
% Education (NfE lab) at The University of Hong Kong for the use in PAC2018
% competition only.
% 
% Please note that due to the time limit, the code has not been cleaned up,
% optimized and structrued. Therefore, the comprehensibility of the code
% can be low.
%
% NfE lab owns the copyright of this program. The program is not allowed to be redistributed.


    Wi = [];
    Wi(:,:,1) = inv(W_standard);
    Wi(:,:,2) = inv(W_to_be_compared);
    Wi_original = Wi;
    

%%%%%%%clustering (matching) part%%%
ELE  = size(Wi,1);
ch = size(Wi,2);n_set = size(Wi,3);
chains = [1:ch]';
chains_ave_corr = [1:ch]';


%iterate
chains = repmat(1:ELE,n_set,1)';
chains_ave_corr_sub1 = [];
for iter = 1:5
for s1 = 1:n_set
%     disp([s1,iter]);
    new_ch = [];
    for ch1 = 1:ch
        corr_s1 = [];
        corr_s2 = [];
        left_elec = setdiff([1:ch],new_ch(1:ch1-1));
    for ch2 = 1:length(left_elec)
        factor1 = [];factor2 = [];
        for chain_i = 1:size(chains,2)
            
            map1 = Wi(:,chains(ch1,chain_i),chain_i);
            map2 = Wi(:,left_elec(ch2),s1);
            %covariance
            tmp = mean(mean(map1.*map2));
            factor1(chain_i) = abs(tmp);
        end
        factor1(s1) = [];
        corr_s1(ch2) = mean(factor1);
    end
    
    
    
    [a b] = max(corr_s1);new_ch(ch1) = left_elec(b(1));
    chains_ave_corr_sub1(ch1,s1) = corr_s1(b(1));
    end
    chains(:,s1) = new_ch(:);
    
end
chains_ave_corr = mean(chains_ave_corr_sub1,2);
[tem6,tem7] = sort(-chains_ave_corr);
chains = chains(tem7,:);
end


%%%correaltion

corrs = [];
for c = 1:ELE
    topoi1 = Wi_original(:,chains(c,1),1);
    topoi2 = Wi_original(:,chains(c,2),2);
    corrs(c) = abs(corr(topoi1,topoi2));
end

good_or_bad = sum(corrs>0.95)>80;
disp(good_or_bad);
