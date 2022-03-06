function [IDIk,support]=trace_idi(S,kset,OPTION,n)
% TRACE-IDI  Computes the index of dispersion for intervals of a trace
% [IDIk,J]=trace_idi(S,k) - computes IDI(k)- J is the number of points used
%
% [IDIk,J]=trace_idi(S,k,'aggregate',n) - computes IDI(k) using samples
% of S=X1+...X_n
%
% [IDIk,J]=trace_idi(S,k,'aggregate-mix',n) - computes IDI(k) using samples
% of S_t=X1+...X_n(t), i.e., the cardinality of the t-th sum sample is n(t)
%
% see Sriram, Whitt, "Characterizing Superposition Arrival Processes...", JSAC 6, 1986
IDIk=[];
S=S(:);
if nargin<2
    [IDIk,support]=trace_idi(S,min([1000,ceil(length(S)/30)]));
    return
end
for k=kset
    if nargin<3
        support=length(S)-k-1;
        Sk=zeros(1,length(S)-k);
        for t=1:(length(S)-k-1)
            Sk(t)=sum(S(t:(t+k-1))); % for all subsets of length k
        end
        IDIk(end+1)=k*var(Sk)/mean(Sk)^2;
    elseif strcmpi(OPTION,'aggregate')
        % data is already aggregated
        % S is a sample of S=X1+...X_n
        keff=floor(k/n);
        Sk=zeros(1,length(S)-keff);
        support=length(S)/(keff);
        for t=1:(length(S)-keff-1)
            Sk(t)=sum(S(t:(t+keff-1))); % for all subsets of length k
        end
        IDIk(end+1)=k*var(Sk)/mean(Sk)^2;
    elseif strcmpi(OPTION,'aggregate-mix')
        % data is already aggregated
        % S is a sample of S=X1+...X_n(i)
        Sk=[];
        % find subsets that sum the closest possible to k
        y=[];
        for t0=1:(length(n)-ceil(sum(n)/k))
            acc=0;
            indexes=[1,0];
            for t=t0:length(n)
                acc=acc+n(t);
                if acc > k
                    indexes(end,2)=t-1;
                    indexes(end+1,1:2)=[t,0];
                    y(end+1)=acc-n(t);
                    acc=n(t);
                end
            end
            indexes(end,:)=[];
            for i=1:size(indexes,1)
                Sk(i)=sum(S(indexes(i,1):indexes(i,2)));
            end
        end
        IDIk(end+1)=k*var(Sk)/mean(Sk)^2;
    end
end
end
