function [u] = RAP_MUSIC(Y,A,nsources)
%let's say we have p1,p2,...pn locations in region of interest (ROI)
%lj = l(pj) j = 1.... n. n is the number of locations in the ROI
%n is number of dipoles as well.
%lj is the topograph of location n, independent of the orientations dipole
%m = number of sensors, N time points
%A = [l1,l2,...,ln] size of mxn
%S = spirce nxN time cource metrix
%epsilon = Noise mxN
% Y = AS + epsilon -> measurement data metrix that we get from MEG
% sometimes works on noisy data as well! because we are aware to the number
% of sources
X = size(A);
m = X(1);
n = X(2);
[U,~,~] = svd(Y*Y');
u = [];

%by page 160, first step is to use music localizer
all_val1=zeros(n,1);
%u(p) localizer is ||Psg*l(p)||^2/||l(p)||^2 by page 155
%for fixed dipoles : ||Psg*l(p)||^2 = l(p)'*Psg*l(p) by page 158
% ||l(p)||^2 = l(p)'*l(p)
Psg = U(:,1:nsources); % 7.2 page 154
PsgtPsg = Psg*Psg';

%% Vectorization
all_val1 = (sum(A.*((PsgtPsg)*A))./sum(A.*A))';
% for dipole = 1:n
%     l_p = A(:,dipole);
%     all_val1(dipole) = (l_p'*PsgtPsg*l_p)/(l_p'*l_p);
% end

% Global maximum
[~,u]=max(all_val1);
% figure
% plot(all_val1)

%next step:

for j = 2:nsources
    B = A(:,u);
    BPB = B*pinv(B);
    Q_k = eye(size(BPB)) - BPB;
    QkU = Q_k*U(:,1:nsources);
    %P_k = (QkU)*(QkU)';
    %by page 161, u(p) = ||P_k*Q_K*l(p)||^2/||Q_K*l(p)||^2
    %by page 173 u(p) = (l(p)'*Q_K*P_K*Q_K*l(p))/(l(p)'*Q_k*l(p))
    QkPktPkQk = Q_k'*(QkU)*(QkU)'*Q_k;
    QktQk = Q_k'*Q_k;
    %% Vectorization
    all_val1 = (sum(A.*((QkPktPkQk)*A))./sum(A.*(QktQk*A)))';
    %     for dipole = 1:n
    %         l_p = A(:,dipole);
    %         all_val1(dipole) = (l_p'*QkPktPkQk*l_p)/(l_p'*QktQk*l_p);
    %     end

    % ignore previously found dipoles (NaN them)
    for i = 1:size(u,2)
        all_val1(u(i)) = NaN;
    end

    [~,m]=max(all_val1);
    u = [u,m];
end

if nsources ==0
    u = [];
end


