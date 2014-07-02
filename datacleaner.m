function [X y] = datacleaner(X_raw, y_raw, ratedTons)
%% Abstract
% datacleaner() performs basic data cleaning and generates the following
% matrices to be used in the data visualization:
% X = [Status Tons PercentCapacity Fevap Tei Teo Fcond Tco Tci Pei Peo Pco Pci A_kW B_kW Date_Time]
% y = [KWperTon]

%% Create useful variables
% specific heat of water, BTU/lb*degF
c_p = 500; 
% convert BTU/hr to Tons
BTUperHr2Tons = 12000; 

%% Add additional parameters
% add empty columns for Status, Tons, and percent capacity variables
m = size(X_raw, 1);
X_raw = [zeros(m, 1) zeros(m,1) zeros(m,1) X_raw zeros(m,1) zeros(m,1)];
% create Status, Tons, and percent capacity variables
X_raw(:,1) = gt(X_raw(:,4), 300); 
X_raw(:,2) = (X_raw(:,4).*(X_raw(:,5) - X_raw(:,6)).*c_p)/BTUperHr2Tons;
% account for situations where only 1 compressor (of 2) is operational 
for n = 1:m
    if X_raw(n,14) < (X_raw(n,15) - 200)
        X_raw(n,17) = 1;
    end
    if X_raw(n,15) < (X_raw(n,14) - 200)
        X_raw(n,17) = 1;
    end
    if X_raw(n,17) == 1
    X_raw(n,3) = X_raw(n,2)./(ratedTons/2);
    else
    X_raw(n,3) = X_raw(n,2)./(ratedTons);
    end
end

%% Clean data for outliers
% combine X and y for cleaning
Xy = [y_raw X_raw];
% Status OFF filter
l1 = size(Xy,1);
condition1=Xy(:,2)==0;
Xy(condition1,:)=[];
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('Status rows removed %d\n', rmvd); 
% kW/Ton out of range filter
l1 = size(Xy,1);
Xy=Xy(Xy(:,1)>0.15 & Xy(:,1)<0.9, :);
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('kW/Ton rows removed %d\n', rmvd); 
% Percent capacity filter
l1 = size(Xy,1);
Xy=Xy(Xy(:,4)>0.2 & Xy(:,4)<0.95, :);
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('PER rows removed %d\n', rmvd); 
% Fevap filter
l1 = size(Xy,1);
Xy=Xy(Xy(:,5) > 1000 & Xy(:,5)< 7000, :);
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('CHWF rows removed %d\n', rmvd); 
% Teo filter
l1 = size(Xy,1);
Xy=Xy(Xy(:,7) > 37 & Xy(:,7)< 41, :);
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('CHWST rows removed %d\n', rmvd); 
% Tei filter
l1 = size(Xy,1);
Xy=Xy(Xy(:,6) > 45 & Xy(:,6)< 60, :);
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('CHWRT rows removed %d\n', rmvd); 
% Fcond filter
l1 = size(Xy,1);
Xy=Xy(Xy(:,8) > 3000 & Xy(:,8)< 11000, :);
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('CWF rows removed %d\n', rmvd); 
% Repeat data filter
l1 = size(Xy,1);
condition = diff(Xy(:,3));
condition = [1 ; condition];
condition1 = condition == 0;
Xy(condition1,:)=[];
l2 = size(Xy,1);
rmvd = l1 - l2;
fprintf('Repeat rows removed %d\n', rmvd); 


%% Divide data into X (predictors) and y (response)
X = Xy(1:end,2:17);
y = Xy(1:end,1);


end










