function FRSAF=FRSAD(data,delta,sigma)
%%%input:
% Data is data matrix without decisions, where rows for samples and columns for attributes.
% Delta is used to adjust the adaptive fuzzy radius.
%%%output
% Fuzzy rough anomaly subspace factor (FRSAF).

% Handle input parameters and set default values
if nargin < 1  % No input parameters provided
    error('At least the data matrix "data" must be provided as an input parameter');
elseif nargin < 2  % Only "data" is provided (1 input parameter)
    delta = 1;       % Default value for delta
    sigma = 0.01;    % Default value for sigma
elseif nargin < 3  % "data" and "delta" are provided (2 input parameters)
    sigma = 0.01;    % Default value for sigma
end

subattributeset=reduct(data,delta,sigma);

% disp(subattributeset); 

subdata=data(:,subattributeset);

[n,m]=size(subdata);
% disp(m);

LA=1:m;
weight1=zeros(n,m);
weight3=zeros(n,m);
Acc_A_a=zeros(n, m, m);
for col=1:m
    lA_d=setdiff(LA,col);
    Acc_A_a_tem=zeros(n, m);
    
    frm=hamming_matrix(data(:, col),delta);
    [frm_temp,~,ic]=unique(frm,'rows');
    for i=1:size(frm_temp,1)
        RM_temp=repmat(frm_temp(i,:),n,1);
        i_tem=find(ic==i);
        for k=0:m-1
            if k==0
                lA_de=lA_d;
            else
                lA_de=setdiff(lA_d,lA_d(k));
            end
            frm_tmp = hamming_matrix(data(:, lA_de),delta);
            Low_A=sum(min(max(1- frm_tmp,RM_temp),[],2));
            Upp_A=sum(max(min(frm_tmp,RM_temp),[],2));
            Acc_A_a_tem(i_tem,k+1)=Low_A/Upp_A;
        end
        Acc_A_a(:, :, col)=Acc_A_a_tem;
        weight3(i_tem,col)=1-(sum(frm_temp(i,:))/n)^(1/3);
        weight1(i_tem,col)=(sum(frm_temp(i,:))/n);
    end
end
%%
FSDOG=zeros(n,m);
for col = 1:m
    Acc_A_a_tem = Acc_A_a(:,:,col);
    FSDOG(:,col) = 1 - (mean(Acc_A_a_tem,2)).*weight1(:,col);
end
FRSAF=mean(FSDOG.*weight3,2);
end