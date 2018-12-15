function [Iout,time,threshold_mat,fitness,peaksnr]=my_abc1(I,level)
tic;

%% Reading Image

try %#ok<TRYNC>
    disp('Reading Image');
    I = imread(I);                  % Trying to read the image if source is passed
end

if size(I,3)==3 %RGB image
    disp('Converting RGB image to GRAYSCALE');
    I = rgb2gray(I);                % Comverting the RGB image to grayscale for analysis
end

n_countR = imhist(I);             % Calculating the features of the image I

Nt = size(I,1)*size(I,2);
Lmax = 255;                         % 256 different maximum levels are considered in an image (i.e., 0 to 255)

probR = (n_countR/Nt);            % calculating the probability of each pixel in the image
disp(size(probR));
%% Problem Definition (STEP 1)

N_PAR = level-1;                    % Number of threshold values 

%% ABC Settings

MaxIt = 100;                         % Maximum Number of Iterations

N = 30;                            % Population Size (Colony Size)

nOnlooker = N;                      % Number of Onlooker Bees

L = round(0.6*N_PAR*N);             % Abandonment Limit Parameter (Trial Limit)

%% Intialization of parameters for Population

X_MAXR = Lmax*ones(1,N_PAR);        % Adjusting max value x can go upto
X_MINR = zeros(1,N_PAR);             % Adjusting min value x can go upto
fitBestR = zeros(N,1);              % Used for taking the best fit values out of the fitR vector
xR = zeros(N,N_PAR);                % Used for storing the random positions for optimizations

count = 1;                          % Setting the cycle value to 1

C = zeros(N,1);                     % Setting the trail for each solution equal to 0

%% Initializing the Positions

for i = 1:N
    for j = 1:N_PAR
        check = X_MINR(j) + fix(unifrnd(0,1)*(X_MAXR(j)-X_MINR(j)));      % Initially setting up the positions
        while(j>1 && check<xR(i,j-1))
            check = X_MINR(j) + fix(unifrnd(0,1)*(X_MAXR(j)-X_MINR(j)));
        end
        xR(i,j) = check;
    end
end

fitR = Kapur(N,N_PAR,xR,probR);     % Finding the fitness of all points in xR

[~,bR] = max(fitR);                 % bR is storing the index at which fitR gets its maximum value
gBestR = xR(bR,:);                  % Storing the best position till now
gbestvalueR = fitR(bR);             % Storing the best fitness value till now
xBestR = xR;                        % Stores the till now best solution

while(count <= MaxIt)
    % Employed Bees
    %% Initializing the Positions (STEP 2)
    for i = 1:N
        
        K = [1:i-1 i+1:N];          % Maintaining a list for random search such that i is not present
        k = K(randi([1 numel(K)])); % Choose k randomly, not equal to i 
        
        phi = unifrnd(-1,+1); % Choose a random value between the range [-1,1] and size [1,N_PAR]
        
        % Define New Positions
%         check_update = xR(i) + phi.*(xR(i)-xR(k)); % Storing new position for checking
        for j = 1:N_PAR
            check_update = xR(i,j) + fix(phi.*(xR(i,j)-xR(k,j)));
            if check_update>=0 && check_update<=255
                xR(i,j) = check_update;
            end
        end
        
        fitR = Kapur(N,N_PAR,xR,probR);
        for j = 1:N
            if fitR(j) > fitBestR(j)
                fitBestR(j) = fitR(j);
                xBestR(j) = xR(j);
            else
                C(i) = C(i)+1;
            end
        end
    end
        
    % Calculate Fitness Value and selecting Probability
    new_fit = Kapur(N,N_PAR,xR,probR);
    P = new_fit/sum(new_fit);
    
    % Onlooker Bees
    %% Sending onlooker bees to food sources
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:N];
        k=K(randi([1 numel(K)]));
        
        phi=unifrnd(-1,+1);
        
        for j = 1:N_PAR
            check_update = fix(xR(i,j)+(phi*(xR(i,j)-xR(k,j))));
            if check_update>0 && check_update<255
                xR(i,j) = check_update;
            end
        end
        
        fitR = Kapur(N,N_PAR,xR,probR);
        for j = 1:N
            if fitR(j) > fitBestR(j)
                fitBestR(j) = fitR(j);
                xBestR(j,:) = xR(j,:);
            else
                C(i) = C(i)+1;
            end
        end
    end
    
    % Scout Bees
    %% Sending Scouts to search area
    for i = 1:N
        if C(i) >= L
            for j = 1: N_PAR
                xR(i,j) = X_MINR(j)+fix(unifrnd(0,1)*(X_MAXR(j)-X_MINR(j)));
            end
            PI0 = probR(1:xR(i,1));
            ind = PI0 == 0;
            ind = ind .* eps;
            PI0 = PI0 + ind;
            
            w0 =  sum(PI0);
            H0 = -sum((PI0/w0).*(log2(PI0/w0)));
            fitR(i) = fitR(i) + H0;
            
            for jl = 2: N_PAR
                PI0 = probR(xR(i,jl-1)+1:xR(i,jl));
                ind = PI0 == 0;
                ind = ind .* eps;
                PI0 = PI0 + ind;
                
                w0 =  sum(PI0);
                H0 = -sum((PI0/w0).*(log2(PI0/w0)));
                fitR(i) = fitR(i) + H0;
            end
            C(i) = 0;
        end
    end
    
    [~,bR] = max(fitR);
    if (fitR(bR) > gbestvalueR)
        gBestR=sort(xR(bR,:));
        gbestvalueR = fitR(bR);
    end
    count = count+1;
end


%% Results

Iout = imageGRAY(I,gBestR);
threshold_mat = gBestR;
fitness = gbestvalueR;
peaksnr = psnr(I,Iout);
time = toc;

function imgOut=imageGRAY(img,Rvec)
limites=[0 Rvec 255];
img_size=size(img);
imgOut(:,:)=img*0;

k=1;
    for i= 1:img_size(1,1)
        for j=1:img_size(1,2)
            while(k<size(limites,2))
                if(img(i,j)>=limites(1,k) && img(i,j)<=limites(1,k+1))
                    imgOut(i,j,1)=limites(1,k);
                end
                k=k+1;
            end
            k=1;
        end
    end
