function [Fratiosupek_Combined,S_Reduced,S_Full,Q_FU,Q_RE] = Fstatistics_Free_OrientSimrand4(Y,Gain,AP_max_iters,mode,GridLoc,nmax,Multfactor)
T = size(Y,2);
N_sensors = size(Gain,1);
Y_White = Y;
Gain1 = Multfactor*Gain;
for nsources = 0:nmax
    if nsources == 0
        [S_Reduced{nsources+1,1},~,q_Red,~] = AP_Free_Orient_Amita(Y_White, Gain1, nsources, AP_max_iters,mode,GridLoc);
    else
        S_Reduced{nsources+1,1} = S_Full{nsources,1};  q_Red = q_Fu;
    end
    [S_Full{nsources+1,1},~,q_Fu,~] = AP_Free_Orient_Amita(Y_White, Gain1, nsources+1, AP_max_iters,mode,GridLoc);

    gain_idxR = []; gain_idxF = [];
    q_Reduced = []; q_Full =[];
    Y_Reduced = zeros(N_sensors,T); Y_Full = zeros(N_sensors,T);
    for i = 1:length(S_Reduced{nsources+1,1})
        gain_idxR = [gain_idxR (S_Reduced{nsources+1,1}(1,i)-1)*3+1:S_Reduced{nsources+1,1}(1,i)*3];
        q_Reduced = blkdiag(q_Reduced, (q_Red{1,i}/norm(q_Red{1,i})));
    end
    A_Reduced = Gain1(:,gain_idxR)*q_Reduced; S_Est_Reduced = (A_Reduced'*A_Reduced)\(A_Reduced'*Y_White);
    Y_Reduced = A_Reduced*S_Est_Reduced;
    for i = 1:length(S_Full{nsources+1,1})
        gain_idxF= [gain_idxF (S_Full{nsources+1,1}(1,i)-1)*3+1:S_Full{nsources+1,1}(1,i)*3];
        q_Full = blkdiag(q_Full, (q_Fu{1,i}/norm(q_Fu{1,i})));
    end
    A_Full = Gain1(:,gain_idxF)*q_Full; S_Est_Full = (A_Full'*A_Full)\(A_Full'*Y_White);
    Y_Full = A_Full*S_Est_Full;
    Grand_Reduced = sum((Y_White-Y_Reduced).^2);
    Grand_Full = sum((Y_White-Y_Full).^2);
    SSE_Reduced{nsources+1,1} = Grand_Reduced;
    SSE_Full{nsources+1,1} = Grand_Full;
    Q_FU{nsources+1,1} = q_Full;
    Q_RE{nsources+1,1} = q_Reduced;

    %% F-Ratio ANlysis
    req = 4+T;
    DOF_Reduced = N_sensors*T-req*(nsources);
    DOF_Full = N_sensors*T-req*(nsources+1);
    Fratiosupek_Combined(nsources+1,1) = ((sum(SSE_Reduced{nsources+1,1}))/(DOF_Reduced))/((sum(SSE_Full{nsources+1,1})/DOF_Full));
end




