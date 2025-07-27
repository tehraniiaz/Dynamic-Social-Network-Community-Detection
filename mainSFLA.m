%% SFLA Community Detection

clc;
clear;
close all;

%% Problem Definition

fitnessFuncChoice = 2;  % 1: for Modularity 2: for NMI
fitnessFuncStr = {'fitnessModularity','fitnessNMI'};

% load data
load dolphins.mat
data = Problem.A;
data = full(data);
community = community; %#ok

N = size(data , 1);

% Adjency Cell
for i=N :-1: 1
    M{i} = find(data(i,:)==1);
end

% Adjency Matrix
maxLen = max(cellfun(@(x)length(x),M));
adjM = cell2mat(cellfun(@(x)[x(:);zeros(maxLen - length(x),1)],M,'UniformOutput',false));

% Objective Function
fitnessFuncs = {@(x)fitnessModularity(x,adjM),@(x)fitnessNMI(x,community,adjM)};
CostFunction = fitnessFuncs{fitnessFuncChoice};

nVar = N;                                   % Number of Unknown Variables
VarSize = [1 nVar];                         % Unknown Variables Matrix Size

VarMax = cellfun(@length,M);                % Upper Bound of Unknown Variables
% VarMin = ones(size(VarMax)).*(VarMax>0);          % Lower Bound of Unknown Variables
VarMin = zeros(1,N);                        % Lower Bound of Unknown Variables

%% SFLA Parameters

MaxIt = 100;        % Maximum Number of Iterations
maxStillIter = 25;
stillResultsTol = 1e-3;

nPopMemeplex = 5;                          % Memeplex Size
nPopMemeplex = max(nPopMemeplex, nVar+1);   % Nelder-Mead Standard

nMemeplex = 20;                  % Number of Memeplexes
nPop = nMemeplex*nPopMemeplex;	% Population Size

I = reshape(1:nPop, nMemeplex, []);

% FLA Parameters
fla_params.q = max(round(0.2*nPopMemeplex),2);   % Number of Parents
fla_params.alpha = 5;   % Number of Offsprings
fla_params.beta = 5;    % Maximum Number of Iterations
fla_params.sigma = 2;   % Step Size
fla_params.CostFunction = CostFunction;
fla_params.VarMin = VarMin;
fla_params.VarMax = VarMax;

%% Initialization

% Empty Individual Template
empty_individual.Position = [];
empty_individual.Cost = [];

% Initialize Population Array
pop = repmat(empty_individual, nPop, 1);

% Initialize Population Members
for i=1:nPop
    pop(i).Position = arrayfun(@(ind)randi([VarMin(ind), VarMax(ind)],1),1:nVar);
    pop(i).Cost = CostFunction(pop(i).Position);
end

% Sort Population
pop = SortPopulation(pop);

% Update Best Solution Ever Found
BestSol = pop(1);

% Initialize Best Costs Record Array
BestCosts = nan(MaxIt, 1);

%% SFLA Main Loop

for it = 1:MaxIt

    fla_params.BestSol = BestSol;

    % Initialize Memeplexes Array
    Memeplex = cell(nMemeplex, 1);

    % Form Memeplexes and Run FLA
    for j = 1:nMemeplex
        % Memeplex Formation
        Memeplex{j} = pop(I(j,:));

        % Run FLA
        Memeplex{j} = RunFLA(Memeplex{j}, fla_params);

        % Insert Updated Memeplex into Population
        pop(I(j,:)) = Memeplex{j};
    end

    % Sort Population
    pop = SortPopulation(pop);

    % Update Best Solution Ever Found
    BestSol = pop(1);

    % Store Best Cost Ever Found
    BestCosts(it) = BestSol.Cost;

    % Show Iteration Information
    fprintf('Iteration %03d: --- Modularity = %10.6f --- NMI = %10.6f \n', it,-fitnessFuncs{1}(BestSol.Position),-fitnessFuncs{2}(BestSol.Position));

    % Stopping the Loop after getting almost still results for some iterations
    if it > maxStillIter
        minVal = min(abs(BestCosts(it+(-maxStillIter+1:0))));
        maxVal = max(abs(BestCosts(it+(-maxStillIter+1:0))));
        if (maxVal/minVal) <= (1 + stillResultsTol)
            break
        end
    end

end

%% Results Extraction

bestPos = BestSol.Position;
colVec = 1:length(bestPos);
colVec(bestPos==0) = [];
rowVec = bestPos(bestPos~=0);
bestL=zeros(size(bestPos));
bestL(bestPos~=0) = adjM(sub2ind(size(adjM),rowVec,colVec));

%% Results

figure;
semilogy(1:it,-BestCosts(1:it),'LineWidth', 2), grid on, grid minor
xlabel('Iteration'), ylabel('Best Cost');
title([ fitnessFuncStr{fitnessFuncChoice}(8:end),' Values'])

%% Graph Formation & Plot

startGraph = graph(data);

commSepMat = zeros(length(bestPos));
colVec = 1:length(bestPos);
colVec(bestPos==0) = [];
indVec = sub2ind(size(data),bestL(bestPos~=0),colVec);
commSepMat(indVec) = 1;
commSepGraph = digraph(commSepMat);

bins = conncomp(commSepGraph,'Type','weak');
colorOrd = colormap(hsv(max(max(bins),max(community))));
colorMat = colorOrd(bins,:);

figure
subplot(1,3,1),plot(startGraph,'EdgeColor','#0072BD','NodeColor',"#7E2F8E",'LineStyle','--','MarkerSize',7),title('Starting Graph')

subplot(1,3,2),plot(commSepGraph,'EdgeColor','#0072BD','NodeColor',colorMat,'LineStyle','--','MarkerSize',7)
title(['After Community Detection by ', fitnessFuncStr{fitnessFuncChoice}(8:end),' Function'])
dim = [.1 .95 .08 .04];

XR = xlim;
X0 = dim(1)*diff(XR) + XR(1);
YR = ylim;
Y0 = dim(2)*diff(YR) + YR(1);
text(X0,Y0,['Modularity = ',num2str(-fitnessFuncs{1}(BestSol.Position)),newline,' NMI = ', num2str(-fitnessFuncs{2}(BestSol.Position))],'FontSize',12)

trueCommMat = trimCommunityMat(data,community);
trueCommMat = digraph(trueCommMat);

subplot(1,3,3),plot(trueCommMat,'EdgeColor','#0072BD','NodeColor',colorOrd(community,:),'LineStyle','--','MarkerSize',7)
title('True Communities')


%% 
timeStr = char(string(datetime('now','Format','MMMdd_HH_mm_ss')));

fitnessStr = {'Q','NMI'};
hndl = gcf;
while(~isempty(hndl))
    hndl.Units = 'normalized';
    hndl.Position = [0 0 1 1];
    print(hndl,['fig',fitnessStr{fitnessFuncChoice},'-',num2str(hndl.Number),'_',timeStr],'-dpng');
    % exportgraphics(hndl,['fig',fitnessStr{fitnessFuncChoice},'-',num2str(hndl.Number),'_',timeStr,'.png']);
    % close(hndl)
    hndl = findobj('Type','figure','Number',hndl.Number - 1);
end
