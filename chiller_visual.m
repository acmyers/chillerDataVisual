%% Chiller Data Visualization
% @author: Andrew Myers
    
%% Findings
% - The visual maps out conditions where the chiller typically operates over course of a year
% - Tei tends to increase as Tci increases
% - kW/Ton tends to increase with Tci
% - Greater uncertainty (i.e. more outliers) under warmer conditions
% - Visualizaiton helps identify additional outliers
% - Potential useful method for visualizing other operating characteristics of other turbomachinery

%% Inititialization
clear ; close all; clc

%% Load Data
% import raw data from csv file
data = csvread('chiller_data.csv',1,0);

% define parameters
Date_Time = data(:,1); % data and time in excel serial format
KWperTon = data(:,2); % measured kW per Ton of cooling
Teo = data(:,3); % temperature of water at evaporator outlet (degrees F)
Tei = data(:,4); % temperature of water at evaporator inlet (degress F)
Fevap = data(:,5); % flow rate of water through evaporator (gpm)
Tci = data(:,6); % temperature of water at condenser inlet (degrees F)
Tco = data(:,7); % temperature of water at condenser outlet (degress F)
Fcond = data(:,8); % flow rate of water through condenser (gpm)
Pei = data(:,9); % pressure reading at evaporator inlet (psi)
Peo = data(:,10); % pressure reading at evaporator outlet (psi)
Pco = data(:,11); % pressure reading at condenser outlet (psi)
Pci = data(:,12); % pressure reading at condenser inlet (psi)
A_kW = data(:,13); % power consumption for compressor A (kW)
B_kW = data(:,14); % power consumption for compressor B (kW)

% combine parameters
X_raw = [Fevap Tei Teo Fcond Tco Tci Pei Peo Pco Pci A_kW B_kW Date_Time];
y_raw = KWperTon;

% add chiller properties
ratedTons = 3700;

% clean the data
[X y] = datacleaner(X_raw, y_raw, ratedTons);

%% Organize data into subsets
% This section will take the data and create data subsets based on specified
% predictors. Note, we are trying to visualize the change in response 
% variable given a change in these (and another) predictor variables. 
% In the example below, Tei is binned by intervals of 1 degF across its 
% range of possible values and Tci is binned by intervals of 2 degF 
% across its range of possible values.  When the indices from both Tei and 
% Tco bins match it constitutes a data subset i.e. matrix entry for the visual.
% In order to experiment with other predictors simply change the variable
% that is binned.  In order to add additional predictors (i.e. dimensions)
% create a new binVariable and add it to where subsets and subsetLabels are
% created.  This will require adding an additional for loop for each new 
% variable and a nested intersect() function. Note: changing/adding predictors 
% will also require some manipulation of the output labels throughout.

% create intervals to bin by Tei
binTei = num2cell((46:1:55)');
% retrieve indices for each bin
for i = 1 : length(binTei)-1
    binTei{i,2} = find(X(:,5) <= cell2mat(binTei(i+1,1)) & X(:,5) > cell2mat(binTei(i,1))); %creates subsets based on Tei %: means all xdir
end

% create intervals to bin by Tci
binTci = num2cell((54:2:80)');
% retrieve indices for each bin
for i = 1 : length(binTci)-1
    binTci{i,2} = find(X(:,9) <= cell2mat(binTci(i+1,1)) & X(:,9) > cell2mat(binTci(i,1)));
end

% combine Tei and Tci bins to create data subsets
subsets = cell(length(binTei)-1,length(binTci)-1);
subsetLabels = cell(length(binTei)-1,length(binTci)-1);
for i = 1 : length(binTei)-1
    for n = 1 : length(binTci)-1
        subsets{i,n} = intersect( binTei{i,2}(:,1),binTci{n,2}(:,1));
        subsetLabels{i,n} = ['Tei >', num2str(binTei{i,1}), ' and Tci >', num2str(binTci{n,1})];
    end
end


%% Create Chiller Data Visualization (using subsets as matrix entries)
% To create the visual I used the subplot_pos() function available here: 
% {http://p-martineau.com/perfect-subplot-in-matlab/}
% Given the subset indices, Ihave visualized kW/Ton (response variable) against
% percent capacity (additonal predictor variable mentioned above).  The 
% result is a matrix with each entry displaying kW/Ton vs. Percent_Capacity, 
% while Tci changes the matrix x-axis and Tei changes on the matrix y-axis.
% NOTE: I've included a polynomial fit for visualization only, meaning it 
% has not been checked for accuracy on validation and/or test data.  
% However, the code could still be useful for those interested.

% parameters for figure and panel size
xdir = size(subsets,2); % setting Tci change in x direction of subplot
ydir = size(subsets,1); % setting Tei change in y direction of subplot
plotheight=30;
plotwidth=30;
subplotsx=xdir;
subplotsy=ydir;   
leftedge=2.5;
rightedge=0.4;   
topedge=1;
bottomedge=2.5;
spacex=0.2;
spacey=0.2;
fontsize=11;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% set the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

% define axes
perrange = 0.2:0.01:1; % range of percent capacity values used for prediction interval
kptaxis = 0.2:0.2:0.8; % kW/Ton axis labels
peraxis = 0.25:0.25:0.75; % percent capacity axis labels

% preallocate space
modelSummary = cell(xdir*ydir,10);
p_ints = cell(xdir*ydir,10);

j = 1;
for i = 1 : xdir;
    for ii = 1 : ydir;
        
        dataSubset = subsets{ii,i}(:,:); 
        modelSummary{j,1} = subsetLabels{ii,i}; 
        
        if length(dataSubset) > 25
            
        % create polynomial fit (power N = 2)
        [quadmod S] = polyfit(X(dataSubset,3), y(dataSubset), 2);
        quadmodcurve = polyval(quadmod, X(dataSubset,3)); 
        [r2 rmse] = rsquare(y(dataSubset), quadmodcurve, true);
        res = y(dataSubset) - quadmodcurve;
        
        % output summary of subset data and polyfits
        modelSummary{j,2} = mean(y(dataSubset));
        modelSummary{j,3} = mean(X(dataSubset,2));
        modelSummary{j,4} = length(dataSubset);
        modelSummary{j,5} = r2;
        modelSummary{j,6} = rmse;
        
        % add prediction intervals to the plots
        [Y,DELTA] = polyconf(quadmod,perrange,S);
        
        % plot the subplots
        quadmodcurve = polyval(quadmod,perrange); 
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top'); %#ok<LAXES>
        plot(X(dataSubset,3), y(dataSubset), '.', perrange, quadmodcurve, 'g-', perrange,Y+DELTA,'r--', perrange,Y-DELTA,'r--');
        axis([0 1 0 1]);  
        set(gca,'XTick', peraxis, 'YTick', kptaxis,'TickLength',[0.025 0.025]);
        
        else
        quadmodcurve = 0; 
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top'); %#ok<LAXES>
        plot(0, 0);
        axis([0 1 0 1]);
        set(gca,'XTick', peraxis, 'YTick', kptaxis,'TickLength',[0.025 0.025]);
        end
        
        if ii>1
        set(ax,'xticklabel',[])
        end
        
        if i>1
        set(ax,'yticklabel',[])
        end
        
        if i==1
        ylab = {['T_{EI}',blanks(1),'>',blanks(1), num2str(binTei{ii,1})]; 'kW/Ton'};
        ylabel(ylab);
        end
 
        if ii==1
        xlab = {'%capacity',['T_{CI}',blanks(1),'>',blanks(1),num2str(binTci{i,1})]};
        xlabel(xlab)
        end
 
        j = j + 1;
    end
end



