function [Fratio, ASFULL] = Fstatistics(Y,ngrids,gainMatrix,AP_max_iters,T,N_sensors,nmax,flag)
%%
for nsources = 0:nmax
    if strcmp('RAP',flag)
        [S_Reduced] = RAP_MUSIC(Y, gainMatrix,nsources);
        [S_Full] = RAP_MUSIC(Y, gainMatrix,nsources+1);
    elseif strcmp('AP',flag)
        [~,S_Reduced] = alternating_projections(Y, ngrids, gainMatrix, nsources , AP_max_iters,'AP');
        [~,S_Full] = alternating_projections(Y, ngrids, gainMatrix, nsources + 1, AP_max_iters,'AP');
    end
    A_Reduced = gainMatrix(:,S_Reduced); S_Est_Reduced = (A_Reduced'*A_Reduced)\(A_Reduced'*Y);
    A_Full = gainMatrix(:,S_Full); S_Est_Full = (A_Full'*A_Full)\(A_Full'*Y);
    Y_Reduced = A_Reduced*S_Est_Reduced;
    Y_Full = A_Full*S_Est_Full;
    Grand_Reduced = (sum((Y-Y_Reduced).^2));
    Grand_Full = (sum((Y-Y_Full).^2));
    %%
    ASFULL{1,nsources+1} = S_Full;

    SSE_Reduced = Grand_Reduced;
    SSE_Full = Grand_Full;
    req = 3+T;
    DOF_Reduced = N_sensors*T-req*(nsources);
    DOF_Full = N_sensors*T-req*(nsources+1);
    Fratio(nsources+1,1) = ((sum(SSE_Reduced))/(DOF_Reduced))/((sum(SSE_Full)/DOF_Full));
end
