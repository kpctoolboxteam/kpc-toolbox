function JM=map_joint(MAP,a,i)
% map_joint computes the joint moments of a MAP
%
% JM=map_joint(MAP,a,i) returns the joint moment of a MAP
% 
%   Input:      
%           MAP = a map in the form of {D0,D1}
%           a   = a vector (a1, a2, a3, ... ) which are the subscripts 
%                 of each term in E[(X_a1)^i1*(X_a2)^i2*...]
%           i   = a vector (i1, i2, i3, ... ) which specifying the 
%                 power of each element in the joint moments 
%                 E((X_a1)^i1*(X_a2)^i2*(X_a3)^i3*(X_a4)^i4*...]
%
%   Output: 
%           JM = E[(X_a1)^i1*(X_{a1+a2})^i2*(X_{a1+a2+a3})^i3*... ]
%
%
%  
%
a=cumsum(a);                 % the cumulative vector of vector a
P=inv(-MAP{1})*MAP{2};
invD0=inv(-MAP{1});
JM=1;
K=length(a);
for k=1:(K-1)
%    a(k+1)-a(k)
    JM=JM*factorial(i(k))*invD0^(i(k))*(P^(a(k+1)-a(k)));
end
%a(K)
%JM=map_pie(MAP)*JM*factorial(i(K))*invD0^(i(K))*(P^a(K))*e(length(P));
JM=map_pie(MAP)*JM*factorial(i(K))*invD0^(i(K))*e(length(P));
                            % e is a column vector of ones

    
        
    
        function p=map_pie(MAP)
        % function p=map_pie(MAP)
        %
        %   Input: a map in the form of {D0,D1}
        %
        %   Output: p = the probability vector of states in a MAP, where 
        %               p*(-D0)^(-1)*D1 = p
        %   
        %   Comments: in this probability vector p, p(i) means that after one
        %             arrival has been generated, the probability the next arrival 
        %             will enter the markov chain through state i 
        %             
            P  = inv(-MAP{1})*MAP{2}; 
            Qt=(P-eye(size(P)));
            p=ctmc_solve(Qt);

                  

                            function p=ctmc_solve(Q)
                            % FUNCTION p=ctmc_solve(Q) - solve ctmc by global balance
                            % 
                            %   Input: parameter Q is the infinitesimal generator matrix Q
                            %   Output: p is the static probability vector 
                            %   
                            %    Note:  this solver didn't use any efficient algorithm, 
                            %           just for simple small Q matrices
                                zerocol=find(sum(abs(Q))==0); 
                                if length(zerocol)>1 
                                    warning('two cols of Q are composed by zeros only.');        
                                    b=zeros(size(Q,1),1);
                                elseif length(zerocol)==0
                                    zerocol=size(Q,1);
                                end
                                b=zeros(size(Q,1),1); b(zerocol)=1;
                                Q(:,zerocol)=1; % add normalization condition    
                                p=(Q')\b; % compute steady state probability vector
                                p=p'; % p is a row vector
                            end
                 

        end
    
        function n=e(t)
            if nargin<1
                n=[1;1]
            else
                n=ones(t,1);
            end
        end

    


end