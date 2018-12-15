function [Iout,intensity,fitness,time]=new_segment(I,level)
tic;

I = imread(I);

if size(I,3)==3 %RGB image
    I = rgb2gray(I);
end

[n_countR,x_valueR] = imhist(I);

Nt=size(I,1)*size(I,2);
Lmax=256;   %256 different maximum levels are considered in an image (i.e., 0 to 255)

for i=1:Lmax
    probR(i)=n_countR(i)/Nt;
end

N = 150; %predefined PSO population for multi-segmentation

N_PAR = level-1;  %number of thresholds (number of levels-1)

N_GER = 150; %number of iterations of the PSO algorithm

PHI1 = 0.8;  %individual weight of particles
PHI2 = 0.8;  %social weight of particles
W = 1.2;   %inertial factor

vmin=-5;
vmax=5;

vR=zeros(N,N_PAR);  %velocities of particles
X_MAXR = Lmax*ones(1,N_PAR);
X_MINR = ones(1,N_PAR);
gBestR = zeros(1,N_PAR);
gbestvalueR = -10000;
gauxR = ones(N,1);
xBestR=zeros(N,N_PAR);
fitBestR=zeros(N,1);
fitR = zeros(N,1);
xR = zeros(N,N_PAR);
for i = 1: N
    for j = 1: N_PAR
        xR(i,j) = fix(rand(1,1) * ( X_MAXR(j)-X_MINR(j) ) + X_MINR(j));
    end
end

% for si=1:length(xR)
%     xR(si,:)=sort(xR(si,:));
% end

nger=1;

for j=1:N
    fitR(j)=sum(probR(1:xR(j,1)))*(sum((1:xR(j,1)).*probR(1:xR(j,1))/sum(probR(1:xR(j,1)))) - sum((1:Lmax).*probR(1:Lmax)) )^2;
    for jlevel=2:level-1
        fitR(j)=fitR(j)+sum(probR(xR(j,jlevel-1)+1:xR(j,jlevel)))*(sum((xR(j,jlevel-1)+1:xR(j,jlevel)).*probR(xR(j,jlevel-1)+1:xR(j,jlevel))/sum(probR(xR(j,jlevel-1)+1:xR(j,jlevel))))- sum((1:Lmax).*probR(1:Lmax)))^2;
    end
    fitR(j)=fitR(j)+sum(probR(xR(j,level-1)+1:Lmax))*(sum((xR(j,level-1)+1:Lmax).*probR(xR(j,level-1)+1:Lmax)/sum(probR(xR(j,level-1)+1:Lmax)))- sum((1:Lmax).*probR(1:Lmax)))^2;
    fitBestR(j)=fitR(j);
end
[aR,bR]=max(fitR);
gBestR=xR(bR,:);
gbestvalueR = fitR(bR);
xBestR = xR;

while(nger<=N_GER)
    i=1;
    
    randnum1 = rand ([N, N_PAR]);
    randnum2 = rand ([N, N_PAR]);
    
    vR = fix(W.*vR + randnum1.*(PHI1.*(xBestR-xR)) + randnum2.*(PHI2.*(gauxR*gBestR-xR)));
    vR = ( (vR <= vmin).*vmin ) + ( (vR > vmin).*vR );
    vR = ( (vR >= vmax).*vmax ) + ( (vR < vmax).*vR );
    xR = xR+vR;
    xR = ( (xR <= X_MINR(1)).*X_MINR(1) ) + ( (xR > X_MINR(1)).*xR );
    xR = ( (xR >= X_MAXR(1)).*X_MAXR(1) ) + ( (xR < X_MAXR(1)).*xR );
    
    
    for j = 1:N
        for k = 1:N_PAR
            if (k==1)&&(k~=N_PAR)
                if xR(j,k) < X_MINR(k)
                    xR(j,k) = X_MINR(k);
                elseif xR(j,k) > xR(j,k+1)
                    xR(j,k) = xR(j,k+1);
                end
            end
            if ((k>1)&&(k<N_PAR))
                if xR(j,k) < xR(j,k-1)
                    xR(j,k) = xR(j,k-1);
                elseif xR(j,k) > xR(j,k+1)
                    xR(j,k) = xR(j,k+1);
                end
            end
            if (k==N_PAR)&&(k~=1)
                if xR(j,k) < xR(j,k-1)
                    xR(j,k) = xR(j,k-1);
                elseif xR(j,k) > X_MAXR(k)
                    xR(j,k) = X_MAXR(k);
                end
            end
            if (k==1)&&(k==N_PAR)
                if xR(j,k) < X_MINR(k)
                    xR(j,k) = X_MINR(k);
                elseif xR(j,k) > X_MAXR(k)
                    xR(j,k) = X_MAXR(k);
                end
            end
        end
    end
    
    while(i<=N)
        if(i==N)
            for j=1:N
                fitR(j)=sum(probR(1:xR(j,1)))*(sum((1:xR(j,1)).*probR(1:xR(j,1))/sum(probR(1:xR(j,1)))) - sum((1:Lmax).*probR(1:Lmax)) )^2;
                for jlevel=2:level-1
                    fitR(j)=fitR(j)+sum(probR(xR(j,jlevel-1)+1:xR(j,jlevel)))*(sum((xR(j,jlevel-1)+1:xR(j,jlevel)).*probR(xR(j,jlevel-1)+1:xR(j,jlevel))/sum(probR(xR(j,jlevel-1)+1:xR(j,jlevel))))- sum((1:Lmax).*probR(1:Lmax)))^2;
                end
                fitR(j)=fitR(j)+sum(probR(xR(j,level-1)+1:Lmax))*(sum((xR(j,level-1)+1:Lmax).*probR(xR(j,level-1)+1:Lmax)/sum(probR(xR(j,level-1)+1:Lmax)))- sum((1:Lmax).*probR(1:Lmax)))^2;
                if fitR(j) > fitBestR(j)
                    fitBestR(j) = fitR(j);
                    xBestR(j,:) = xR(j,:);
                end
            end
            [aR,bR] = max (fitR);
            if (fitR(bR) > gbestvalueR)
                gBestR=xR(bR,:)-1;
                gbestvalueR = fitR(bR);
            end
            nger=nger+1;
        end
        i=i+1;
    end
end
%     gbestvalueR

    gBestR=sort(gBestR);
    Iout=imageGRAY(I,gBestR);

if nargout>1
        gBestR=sort(gBestR);
        intensity=gBestR;     %return optimal intensity
        if nargout>2
            fitness=gbestvalueR;    %return fitness value
            if nargout>3
                time=toc;   %return CPU time
            end
        end
end

function imgOut=imageGRAY(img,Rvec)
limites=[0 Rvec 255];
tamanho=size(img);
imgOut(:,:)=img*0;
        
cores=colormap(lines)*255;
close all;
%tic
k=1;
    for i= 1:tamanho(1,1)
        for j=1:tamanho(1,2)
            while(k<size(limites,2))
                if(img(i,j)>=limites(1,k) && img(i,j)<=limites(1,k+1))
                    imgOut(i,j,1)=limites(1,k);
                end
                k=k+1;
            end
            k=1;
        end
    end