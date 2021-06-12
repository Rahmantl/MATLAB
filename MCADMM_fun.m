function [ Z,relerr,i,X,Y,timead ] = MCADMM_fun( M,E,tol,Y,Z,maxiter,Lamda,gr,ro1,M_E )
% relerr=norm((M-X*Y),'fro')/(max(1,norm(M,'fro')));
%reserr=norm((M-X*Y).*E,'fro')/(max(1,norm(M.*E,'fro')));
% M_E=M.*E;ro1=1/ro;gr=gama*ro;
[m,n]=size(E);
ro=1e8;%penalty coefficient
Lamda=zeros(m,n);
F=(E<1);%Omega complement
tic
for i=1:maxiter 
    X=(Z*Y')/(Y*Y');
    Y=(X'*X)\(X'*Z);
    XY=X*Y;M_XY=M-XY;
    %Z=XY-ro1*Lamda+(M_XY).*E;
    % Z new (28nov)
    Zomega=((1/(1+ro))*XY-(1/(1+ro))*Lamda+(ro/1+ro)*M).*E;
    ZomegaC=XY.*F;
    Z=Zomega+ZomegaC;
    Lamda=Lamda.*E+gr*(Z-M).*E;
    %
    %    Lamda=Lamda+gr*(Z-M).*E;
    relerr=norm(M_XY,'fro')/(max(1,norm(M,'fro')))
    if relerr<tol
        break
    end
end
timead=toc;
end

