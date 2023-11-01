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
    line_warning(mfilename,'The order of Q is greater than 6000, i.e., %d elements. Press key to continue.\n',length(Q));
    pause;
end

if size(Q)==1
    p = 1;
    nConnComp = 1;
    connComp = 1:length(Q);
    return
end

Q = ctmc_makeinfgen(Q); % so that spurious diagonal elements are set to 0
n = length(Q);

if issym(Q)
    symvariables = symvar(Q); % find all symbolic variables
    B = double(subs(Q+Q',symvariables,ones(size(symvariables)))); % replace all symbolic variables with 1.0
else
    B = abs(Q+Q')>0;
end
[nConnComp, connComp] = weaklyconncomp(B);
if nConnComp > 1
    % reducible generator - solve each component recursively
    line_warning(mfilename,'Reducible generator. No initial vector available, decomposing and solving each component recursively.\n');
    if issym(Q)
        p = sym(zeros(1,n));
    else
        p = zeros(1,n);
    end

    for c=1:nConnComp
        Qc = Q(connComp==c,connComp==c);
        Qc = ctmc_makeinfgen(Qc);
        p(connComp==c) = ctmc_solve(Qc);
    end
    p = p /sum(p);
    return
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
            line_warning(mfilename,'The infinitesimal generator is reducible.\n');
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
Qnnz_1 = Qnnz;
Qnnz(:,end) = 1;
bnnz_1 = Qnnz;
bnnz(end) = 1;

if ~isdeployed
    if issym(Q)
        p = sym(p);
    end
end

warning('off','MATLAB:singularMatrix');
if nargin == 1
    p(nnzel)=Qnnz'\ bnnz;
    if any(isnan(p))
        % verify if this has become reducible        
        if issym(Qnnz)
            symvariables = symvar(Qnnz); % find all symbolic variables
            B = double(subs(Qnnz+Qnnz',symvariables,ones(size(symvariables)))); % replace all symbolic variables with 1.0
        else
            B = abs(Qnnz+Qnnz')>0;
        end
        [nConnComp, connComp] = weaklyconncomp(B);
        if nConnComp > 1
            % reducible generator - solve each component recursively
            if issym(Qnnz)
                p(nnzel) = sym(zeros(1,n));
            else
                p(nnzel) = zeros(1,n);
            end

            for c=1:nConnComp
                Qc = Q(connComp==c,connComp==c);
                Qc = ctmc_makeinfgen(Qc);
                p(intersect(find(connComp==c),nnzel)) = ctmc_solve(Qc);
            end
            p = p /sum(p);            
            return
        end
    end
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

if issym(Q)
    Q=simplify(Q);
end
warning('on','MATLAB:singularMatrix');
end
