clear,close all,clc
% 0: nominal , 1: interval
% X(1)= Age
% X(2)= Sex ( 0:male ; 1:female)
% X(3)= FHxCVD (Family History for CVD)
% X(4)= iDM2  (Diabetes Mellitus (FPG, 2hPP, Drug))
% X(5)= iHWHR  (High waist to hip ratio)
% X(6)= sms2cat (smoking)
% X(7)= sbpcat2
% X(8)= sbpcat3
% X(9)= sbpcat4
% X(10)= Tchcat2
% X(11)= Tchcat3
% X(12)= Tchcat4
% X(13)= Tchcat4

% data
load('E:\arshad\master_project\master project\beta.mat')
load('E:\arshad\master_project\master project\maindata.mat')        % calculated from cox regression

% cox regression
A=B(:,[1:6,9:15]);
C=B(:,7);
D=not(B(:,8));
[b,logL,H,stats] = coxphfit(A,C,'censoring',D);



%type of inputs
for k=1:13
R=A(:,k);
RR=R(1,1);
r=false;
for i=2:length(R)
    for j=1:length(RR)
        if R(i,1)==RR(j,1)
            r=true;
            break;
        end
    end
    if r==false
         RR(j+1,1)=R(i,1);
end
r=false;

end

if length(RR)==2
    type(k,1)=0;
else
    type(k,1)=1;
end
end

% calculate Mi for each Risk Factor
for i=1:13 
if type(i,1)==1
    M(i,1)=mean(A(:,i));
else type(i,1)==0
    M(i,1)=length(find(A(:,i)==1))/length(A(:,i));
end
end

for m=1:length(A)
   
% Inputs
% X=[66;1;0;0;1;0;0;1;0;1;0;0;0];
X=A(m,:);
X=X';

% calculate Risk of CVD with PARS function
for j=1:13
    F(j,1)=beta(j,1)*(X(j,1)-M(j,1));
    j=j+1;
end
F=sum(F(1:13,1));
S=0.97292179;
% S=0.97577599;
P(m,1)= 1-S^exp(F);
end

PP=[];
for y=1:length(P)
   
    if P(y,1) >= 0.15555
        PP(y,1)=1;
    else
        PP(y,1)=0;
    end
    
end

%% AUC calcuation
GOLD=B(:,8);
[x,y,T,auc]=perfcurve(GOLD,P,1);