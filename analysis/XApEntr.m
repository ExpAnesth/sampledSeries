%calculates Cross Approximate Entropy of time series X Y
%r=percentage of the signal standard-deviation used as noise filter (give the value in %, usually 20)
%M=length of the subsequences to be compared, usually 1 or 2
%k -->we compare the prevalence of repetitive patterns of length M with those of length M+k, usually k=1


%Requirements: X,Y are of same size and normalized, i.e.: x^(i) = (x(i) - mean X)/SD(X); same for Y (see Pincus 1996: PNAS 1996;93;14100-14105)

% Author: Matthias Kreuzer (?)

function[XApEn]=XApEntr(X,Y,M,r,k)


N=length(X);
Cm=[];
r=std(X)*r/100
for m=M:k:M+k 
	C=[];
	for i=1:(N-m+1)
		V=[X(i:m+i-1)];
		count=0;
		for j=1:(N-m+1)
			Z=[Y(j:m+j-1)];
			dif=(abs(V-Z)<r);	%two subsequences are similar if the difference between any pair 
                 				%of corrisponding measurements in the pattern is less than r
			A=all(dif);
			count=count+A;
		end
	C=[C count/(N-m+1)];	%vector containing the 
                    		%Cim=(number of patterns similar to the one beginning at interval i)/total number  
                   		%of pattern with the same length M
	end
	Cm=[Cm sum(C)/(N-m+1)];%vector containing the means of the Cim for subsequences of length M and of length M+k%
end
XApEn=log(Cm(1)/Cm(2))