function [S_APP,Loc_APP,q_APP,iter] = AP_Free_Orient_Amita(Y, Gain, nsources,max_iters, mode, GridLoc)
% AP localization of 1st source
if nsources == 0
    S_APP = []; Loc_APP = []; q_APP = []; iter=0;
else
    S_AP=[];
    switch mode
        case 'AP'
            C=Y*Y'+(1e-3)*trace(Y*Y')*eye(size(Y,1)); % Array Covariance matrix
        case 'AP-ISW'
            [u,s,~]=svd(Y*Y'); % Array Covariance matrix
            C=u(:,1:nsources)*s(1:nsources,1:nsources)*u(:,1:nsources)';
        case 'AP-MUSIC'
            [u,~,~]=svd(Y*Y'); % Array Covariance matrix
            C=u(:,1:nsources)*u(:,1:nsources)';
    end

    [~,Cols]=size(Gain);
    ndipoles=Cols/3;
    ap_val1=zeros(ndipoles,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1st Phase
    % a) Initialization: search the 1st source location over the entire
    % dipoles topographies space (ndipoles=15,002 topographies)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for p=1:ndipoles
        gain_idx=(p-1)*3+1:p*3;
        L_p=Gain(:,gain_idx);
        F=L_p'*C*L_p;
        G=L_p'*L_p+(1e-3)*trace(L_p'*L_p)*eye(3);
        [~,D] = eig(F,G);
        D = diag(D);
        ap_val1(p)=max(D);
    end
    % obtain the 1st source location
    [~,s1_idx]=max(ap_val1);
    % obtain the 1st source orientation
    gain_idx=(s1_idx-1)*3+1:s1_idx*3;
    L_p=Gain(:,gain_idx);
    F=L_p'*C*L_p;
    G=L_p'*L_p+(1e-3)*trace(L_p'*L_p)*eye(3);
    [V,D] = eig(F,G);
    D=diag(D);
    [~,max_ev_idx]=max(D);
    q=V(:,max_ev_idx);
    %q=q/norm(q);
    S_AP=[S_AP s1_idx];
    A=L_p*q;

    %store first source loc and orientation
    Loc_AP{1} = GridLoc(s1_idx,:)';
    q_AP{1} = q;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (b) Now, add one source at a time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for q = 2:nsources
        % AP localization of q-th source
        ap_val2=zeros(ndipoles,1);
        P_A=A*pinv(A'*A)*A';
        Q=eye(size(P_A,1))-P_A;
        for p=1:ndipoles
            gain_idx=(p-1)*3+1:p*3;
            L_p=Gain(:,gain_idx);
            F=L_p'*Q*C*Q*L_p;
            G=L_p'*Q*Q*L_p+(1e-3)*trace(L_p'*L_p)*eye(3);
            [~,D] = eig(F,G);
            D=diag(D);
            ap_val2(p)=max(D);
        end
        [~,s2_idx]=max(ap_val2);
        S_AP=[S_AP, s2_idx];

        % obtain the 2nd source orientationReduced
        gain_idx=(s2_idx-1)*3+1:s2_idx*3;
        L_p=Gain(:,gain_idx);
        F=L_p'*Q*C*Q*L_p;
        G=L_p'*Q*Q*L_p+(1e-3)*trace(L_p'*L_p)*eye(3);
        [V,D] = eig(F,G);
        D=diag(D);
        [~,max_ev_idx]=max(D);
        t=V(:,max_ev_idx);
        A=[A L_p*t];

        Loc_AP{q} = GridLoc(s2_idx,:)';
        q_AP{q} = t;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2nd Phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_AP_2=S_AP;
    for iter = 1:max_iters
        for q = 1:nsources
            % AP localization of q-th source
            ap_val2=zeros(ndipoles,1);
            A_TMP=A;
            A_TMP(:,q)=[];
            %GainIndexes = MapDipoles2GainIndexes(S_AP_TMP);
            %A=Gain(:,GainIndexes);
            P_A=A_TMP*pinv(A_TMP'*A_TMP)*A_TMP';
            Q=eye(size(P_A,1))-P_A;
            for p=1:ndipoles
                gain_idx=(p-1)*3+1:p*3;
                L_p=Gain(:,gain_idx);
                F=L_p'*Q*C*Q*L_p;
                G=L_p'*Q*Q*L_p+(1e-3)*trace(L_p'*Q*Q*L_p)*eye(3);
                [~,D] = eig(F,G);
                D=diag(D);
                ap_val2(p)=max(D);
            end
            [~,sq_idx]=max(ap_val2);
            S_AP_2(iter,q)=sq_idx;
            gain_idx=(sq_idx-1)*3+1:sq_idx*3;
            %S_AP(:,q)=sq_idx;
            L_p=Gain(:,gain_idx);
            F=L_p'*Q*C*Q*L_p;
            G=L_p'*Q*Q*L_p+(1e-3)*trace(L_p'*Q*Q*L_p)*eye(3);
            [V,D] = eig(F,G);
            D=diag(D);
            [~,max_ev_idx]=max(D);
            t=V(:,max_ev_idx);
            A(:,q)=L_p*t;


            q_AP2{iter,q} = t;
            Loc_AP2{iter,q} = GridLoc(sq_idx,:)';

        end
        %S_AP_2
        if iter>1 && isequal(S_AP_2(iter,:),S_AP_2(iter-1,:))
            % No improvement vs. previous iteration
            break
        end
    end

    %disp(['Iter = ' num2str(iter)]);

    S_APP=S_AP_2(iter,:);
    Loc_APP = Loc_AP2(iter,:);
    q_APP = q_AP2(iter,:);
end


