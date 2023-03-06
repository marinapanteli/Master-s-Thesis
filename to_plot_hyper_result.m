%% plot viable points for OutV1 after hyperspace
%load('outv1.mat')

% For out chosen topology

% 1) Before Experiment
load('U:\THESIS\Thesis\topodesign-2.0-master\examples\Marina\Save\OutV_7.mat', 'OutV')
paramSpecs_before = 'Copy_of_parametersTHESIS_after_arabinose2.txt'; 
newTable = readtable(paramSpecs_before);
newTable([1,5,23:24,31,39,42,45],:)=[];
bmin = newTable.bmin;
bmax = newTable.bmax;

newTable.islog = zeros(length(bmin),1);
data = OutV; %OutV
data.rowmat = data.V;
paramtoplot = newTable.names;
n = length(paramtoplot);
indices = zeros(1,n);
x = log10(data.rowmat);



% 2) After Experiment
%load('U:\THESIS\Thesis\topodesign-2.0-master\examples\Marina\24_pos_in\Saved_viablePoints_sampled_2104153_00000000000000001100000000000000000000000000000111.mat')
load('resulting_robustnesses_after_module4.mat', 'res')
paramSpecs_after = 'Copy_2_of_parametersTHESIS_finaleval.txt'; 
newTable = readtable(paramSpecs_after);
newTable(5,:)=[];
bmin = newTable.bmin;
bmax = newTable.bmax;

newTable.islog = zeros(length(bmin),1);
data = res(3); %OutV
data.rowmat = data.V;
paramtoplot = newTable.names;
n = length(paramtoplot);
indices = zeros(1,n);
x = log10(data.rowmat);


%apodw

for i = 1:length(bmax)
if bmin(i)==bmax(i)    
        bmin(i)=bmin(i)/2;
        bmax(i)=bmax(i)*2;
end
end

xmin = log10(bmin);
xmax = log10(bmax);

axmin=xmin;
axmax=xmax;


factor = 1.25;
map = copper();

for k=1:(n-1),
    subplot(n,n,(k-1)*n+1)
    [~,density,X,Y] = kde2d([x(:,1),x(:,k)],100);
    contour(X,Y,density,'LineWidth',1);
    hold on;
%         plot(log10(paraoptestimated(1)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
    pos = get(gca, 'Position');
    pos(3) = factor*pos(3);
    pos(4) = factor*pos(4);
    set(gca, 'Position', pos, 'XTick',[]);    
    


   
    axis([axmin(1) axmax(1) axmin(k) axmax(k)]);
    for l =2:(k-1),
        subplot(n,n,(k-1)*n+l)
        [~,density,X,Y] = kde2d([x(:,l),x(:,k)],100);
        contour(X,Y,density,'LineWidth',1);
        hold on;
%             plot(log10(paraoptestimated(l)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
        pos = get(gca, 'Position');
        pos(3) = factor*pos(3);
        pos(4) = factor*pos(4);
        set(gca, 'Position', pos, 'XTick',[],'YTick',[]);       
        axis([axmin(l) axmax(l) axmin(k) axmax(k)]);   
    end
    subplot(n,n,(k-1)*n+k)
    pd = fitdist(x(:,k),'Kernel');
    X = xmin(k):0.1:xmax(k);
    Y = pdf(pd,X);
    plot(X,Y,'Color','black','LineWidth',1);
    title(paramtoplot{k});
    pos = get(gca, 'Position');
    pos(3) = factor*pos(3);
    pos(4) = factor*pos(4);
    set(gca, 'Position', pos, 'XTick',[],'YTick',[]); 

    if max(Y)~=0
        axis([[axmin(k) axmax(k)] 0 max(Y)]);    
    else
        axis([[axmin(k) axmax(k)] 0 0.5]);    
    end
end

k = n;
subplot(n,n,(k-1)*n+1)
[~,density,X,Y] = kde2d([x(:,1),x(:,k)],100);
contour(X,Y,density,'LineWidth',1);
hold on;
%     plot(log10(paraoptestimated(1)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
pos = get(gca, 'Position');
pos(3) = factor*pos(3);
pos(4) = factor*pos(4);
set(gca, 'Position', pos); 
axis([axmin(1) axmax(1) axmin(k) axmax(k)]);

for l =2:(k-1),
    subplot(n,n,(k-1)*n+l)
    [~,density,X,Y] = kde2d([x(:,l),x(:,k)],100);
    contour(X,Y,density,'LineWidth',1);
    hold on;
%         plot(log10(paraoptestimated(l)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
    pos = get(gca, 'Position');
    pos(3) = factor*pos(3);
    pos(4) = factor*pos(4);
%     set(gca, 'Position', pos, 'YTick',[]);     
    set(gca, 'Position', pos, 'XTick',[],'YTick',[]);  
    axis([axmin(l) axmax(l) axmin(k) axmax(k)]);    
end
l = k;
subplot(n,n,(k-1)*n+l)
pd = fitdist(x(:,k),'Kernel');
X = xmin(k):0.1:xmax(k);
Y = pdf(pd,X);
plot(X,Y,'Color','black','LineWidth',1);
title(paramtoplot{k});
pos = get(gca, 'Position');
pos(3) = factor*pos(3);
pos(4) = factor*pos(4);
% set(gca, 'Position', pos, 'YTick',[]);      
set(gca, 'Position', pos, 'XTick',[],'YTick',[]);  
axis([[axmin(k) axmax(k)] 0 max(Y)]); 


colormap(map);

set(gcf,'PaperOrientation','landscape','PaperUnits','Normalized','PaperPosition',[0 0 1 1]);
    print('before_cor', '-dpdf');
