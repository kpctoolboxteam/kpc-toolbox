function [p, Q, nConnComp, connComp]=ctmc_solve(Q,options)
% PROB=ctmc_solve(P) - Equilibrium distribution of the continuous-time
% Markov chain
%
%  Input:
%  P: infinitesimal generator of the continuous-time Markov chain
%
%  Output:
%  PROB: equilibrium distribution of the continuous-time Markov chain
%
%  Examples:
%  - ctmc_solve([-0.5,0.5;0.2,-0.2])


if length(Q) > 6000 && (nargin==1 || ~options.force)
    fprintf(1,'ctmc_solve: the order of Q is greater than 6000, i.e., %d elements. Press key to continue.',length(Q));
    pause;
end

if size(Q)==1
    p = 1;
    return
end

Q = ctmc_makeinfgen(Q); % so that spurious diagonal elements are set to 0
n = length(Q);
if nargout > 2
    [nConnComp, connComp] = weaklyconncomp(Q+Q');
    if nConnComp > 1
        % reducible generator - fail fast
        p = nan(1,n);
        return
    end
end
if all(Q==0)
    p = ones(1,n)/n;
    return
end
p = zeros(1,n);
b = zeros(n,1);

nnzel = 1:n;
Qnnz = Q; bnnz = b;
Qnnz_1 = Qnnz; bnnz_1 = bnnz;

isReducible = false;
goon = true;
while goon
    
    %
    %     zerorow=find(sum(abs(Qnnz),2)==0);
    %         if length(zerorow)>=1
    %             if nargin==1 || options.verbose
    %                 %warning('ctmc_solve: the infinitesimal generator is reducible (zero row)');
    %                 fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
    %                 isReducible = true;
    %             end
    %         end
    %     nnzrow = setdiff(nnzel, zerorow);
    %
    %     zerocol=find(sum(abs(Qnnz),1)==0);
    %     nnzcol = setdiff(nnzel, zerocol);
    %         if length(zerocol)>=1
    %             if ~isReducible && (nargin==1 || options.verbose)
    %                 %warning('ctmc_solve: the infinitesimal generator is reducible (zero column)');
    %                 fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
    %             end
    %         end
    %     nnzel = intersect(nnzrow, nnzcol);
    
    nnzel = find(sum(abs(Qnnz),1)~=0 & sum(abs(Qnnz),2)'~=0);
    if length(nnzel) < n && ~isReducible
        isReducible = true;
        if (nargin > 1 && options.verbose == 2) % debug
            fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
        end
    end
    Qnnz = Qnnz(nnzel, nnzel);
    bnnz = bnnz(nnzel);
    Qnnz = ctmc_makeinfgen(Qnnz);
    if all(size(Qnnz_1(:)) == size(Qnnz(:))) && all(size(bnnz_1(:)) == size(bnnz(:)))
        goon = false;
    else
        Qnnz_1 = Qnnz; bnnz_1 = bnnz; nnzel = 1:length(Qnnz);
    end
end
if isempty(Qnnz)
    p = ones(1,n)/n;
    return
end
Qnnz(:,end) = 1;
bnnz(end) = 1;

if ~isdeployed
    if issym(Q)
        p = sym(p);
    end
end

if nargin == 1
    p(nnzel)=Qnnz'\ bnnz;
else
    switch options.method
        case 'gpu'
            try
                gQnnz = gpuArray(Qnnz');
                gbnnz = gpuArray(bnnz);
                pGPU = gQnnz \ gbnnz;
                gathered_pGPU = gather(pGPU);
                p(nnzel) = gathered_pGPU; % transfer from GPU to local env
            catch
                warning('ctmc_solve: GPU either not available or execution failed. Switching to default method.');
                p(nnzel) = Qnnz'\ bnnz;
            end
        otherwise
            p(nnzel)=Qnnz'\ bnnz;
    end
end

end