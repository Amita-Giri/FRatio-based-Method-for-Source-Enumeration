function [S_AP,S_AP_2] = alternating_projections(Y,ngrids, gainMatrix, nsources,max_iters,mode)
% AP localization of 1st source
Y=(Y);
gainMatrix=(gainMatrix);
S_AP=[];
switch mode
    case 'AP'
        mu=0;%(1e-3);
        C=Y*Y'+mu*trace(Y*Y')*eye(size(Y,1)); % Array Covariance matrix
    case 'AP-w-MUSIC'
        [u,s,~]=svd(Y*Y'); % Array Covariance matrix
        C=u(:,1:nsources)*s(1:nsources,1:nsources)*u(:,1:nsources)';
    case 'AP-MUSIC'
        [u,~,~]=svd(Y*Y'); % Array Covariance matrix
        C=u(:,1:nsources)*u(:,1:nsources)';
end
ap_val1=(zeros(ngrids,1));
C=(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st Phase
% a) Initialization: search the 1st source location over the entire
% dipoles topographies space (ndipoles=15,002 topographies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ap_val1 = (sum(gainMatrix.*(C*gainMatrix))./sum(gainMatrix.*gainMatrix))'; %vectorized, faster
% for p=1:ngrids
%     l_p=gainMatrix(:,p);
%     ap_val1(p)=(l_p'*C*l_p)/(l_p'*l_p);
% end
[~,s1_idx]=max(ap_val1); % obtain the 1st source location
S_AP=[S_AP s1_idx];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b) Now, add one source at a time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q = 2:nsources
    % AP localization of q-th source
    ap_val2=(zeros(ngrids,1));
    A=gainMatrix(:,S_AP);
    P_A=A*pinv(A'*A)*A';
    Q=eye(size(P_A,1))-P_A;
    ap_val2 = (sum(gainMatrix.*(Q*C*(Q*gainMatrix)))./sum(gainMatrix.*(Q*gainMatrix)))'; %vectorized, faster
    %     for k=1:length(S_AP)
    %         p = S_AP(k);
    %         l_p=gainMatrix(:,p);
    %         ap_val2(p)=(l_p'*Q*C*Q*l_p)/(l_p'*Q*l_p);
    %     end
    %         for p=1:ngrids
    %             l_p=gainMatrix(:,p);
    %             ap_val2(p)=(l_p'*Q*C*Q*l_p)/(l_p'*Q*l_p);
    %         end
    [~,s2_idx]=max(ap_val2);
    S_AP=[S_AP, s2_idx];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_AP_2 = S_AP;
for iter = 1:max_iters
    S_AP_2_Prev = S_AP_2;
    for q = 1:nsources
        % AP localization of q-th source
        ap_val2=(zeros(ngrids,1));
        S_AP_TMP = S_AP_2;
        S_AP_TMP(:,q)=[];
        A=gainMatrix(:,S_AP_TMP);
        P_A=A*pinv(A'*A)*A';
        Q=eye(size(P_A,1))-P_A;
        ap_val2 = (sum(gainMatrix.*(Q*C*(Q*gainMatrix)))./sum(gainMatrix.*(Q*gainMatrix)))'; %vectorized, faster
        %         for k=1:length(S_AP_2)
        %             p = S_AP_2(k);
        %             l_p=gainMatrix(:,p);
        %             ap_val2(p)=(l_p'*Q*C*Q*l_p)/(l_p'*Q*l_p);
        %         end
        %                 for p=1:ngrids
        %                     l_p=gainMatrix(:,p);
        %                     ap_val2(p)=(l_p'*Q*C*Q*l_p)/(l_p'*Q*l_p);
        %                 end
        [~,sq_idx]=max(ap_val2);
        S_AP_2(q)=sq_idx;
    end
    if iter>1 && isequal(S_AP_2,S_AP_2_Prev)
        % No improvement vs. previous iteration
        break
    end
end

if nsources ==0
    S_AP_2 = [];
end