clear; clc;

meanredaway = [];
meanredtoward = [];
meanblueaway = [];
meanbluetoward = [];

for zz=1:1
    gridlenx = 100;
    gridleny = 100;
    M_AMT_1 = 200;
    M_AMT_2 = 200;
    D1 = 5*10^(-10);
    D2 = 5*10^(-10);
    S1 = 0.1;
    S2 = 0.1;
    ThresholdforU1 = 1000;
    ThresholdforU2 = 1000;
    UptakeOfM1 = 200;
    UptakeOfM2 = 200;
    dt = 600;
    dl = 3*10^(-3);
    niter = 1500;
    divisions=zeros(niter,gridlenx,gridleny);
    GRID = zeros(niter,gridlenx,gridleny);
    U1 = zeros(niter,gridlenx,gridleny);
    U2 = zeros(niter,gridlenx,gridleny);
    keepcells=zeros(niter,gridlenx,gridleny);
    
    cells=cellclass.empty;
    for i=45:54
        cells(end+1)=cellclass(1,i,50);
        cells(end+1)=cellclass(2,i,60);
    end
    
    for time=2:niter
        %Inject molecules from cell
        for c=1:length(cells)
            row=cells(c).locationx; col=cells(c).locationy;
            if cells(c).ID==1
                U1(time-1,row,col)=U1(time-1,row,col)+M_AMT_1;
            end
            if cells(c).ID==2
                U2(time-1,row,col)=U2(time-1,row,col)+M_AMT_2;
            end
        end
        
        %diffusion
        for row=1:gridlenx
            for col=1:gridleny
                if row-1>0 && col-1>0 && row+1<gridlenx && col+1<gridleny
                    U1(time,row,col) = U1(time-1,row,col) + D1*dt/(dl^2)*(U1(time-1,row-1,col) + U1(time-1,row+1,col) + U1(time-1,row,col-1) + U1(time-1,row,col+1) - 4*U1(time-1,row,col));
                    U2(time,row,col) = U2(time-1,row,col) + D2*dt/(dl^2)*(U2(time-1,row-1,col) + U2(time-1,row+1,col) + U2(time-1,row,col-1) + U2(time-1,row,col+1) - 4*U2(time-1,row,col));
                end
            end
        end
        
        %divisions
        for c=1:length(cells)
            row=cells(c).locationx; col=cells(c).locationy;
            newfitness=cells(c).fitness;
            if cells(c).ID==1
                if U2(time,row,col)>ThresholdforU2
                    U2(time,row,col) = U2(time,row,col) - UptakeOfM2;
                    newfitness = cells(c).fitness + S1;
                end
            else
                if U1(time,row,col)>ThresholdforU1
                    U1(time,row,col) = U1(time,row,col) - UptakeOfM1;
                    newfitness = cells(c).fitness + S2;
                end
            end
            newcell=celldivision(cells(c),cells,newfitness);
            if isa(newcell,'cellclass')
                cells(end+1)=newcell;
                divisions(time,row,col)=1;
            end
        end
        
        for c=1:length(cells)
            if cells(c).ID==1
                keepcells(time,cells(c).locationx,cells(c).locationy)=1;
            end
            if cells(c).ID==2
                keepcells(time,cells(c).locationx,cells(c).locationy)=2;
            end
        end
        %     disp(time)
    end
    % figure
    % for i=1:niter
    %     drawnow
    %     pcolor(reshape(U2(i,:,:),[gridlenx,gridleny]))
    % end
    %%
    
    map=[0,0,0; 1,0,0; 0,0,1];
    figure
    for i=1:niter
        drawnow
        pcolor(reshape(keepcells(i,:,:),[gridlenx,gridleny]))
        colormap(map)
        title(['Iteration = ',num2str(i),' | ','Number of cells = ',num2str(nnz(reshape(keepcells(i,:,:),[gridlenx,gridleny])))])
    end
    
    %%
    redaway=0;
    redtoward=0;
    blueaway=0;
    bluetoward=0;
    for c=1:length(cells)
        if cells(c).ID==1
            if cells(c).locationy<50 && cells(c).locationy>=46
                redaway=redaway+1;
            elseif cells(c).locationy>50
                redtoward=redtoward+1;
            end
        end
        if cells(c).ID==2
            if cells(c).locationy>60 && cells(c).locationy<=64
                blueaway=blueaway+1;
            elseif cells(c).locationy<60
                bluetoward=bluetoward+1;
            end
        end
    end
    meanredaway(end+1) = redaway;
    meanredtoward(end+1) = redtoward;
    meanblueaway(end+1)= blueaway;
    meanbluetoward(end+1)=bluetoward;
    
end
%%
load data.mat
figure()
data=[meanbluetoward';meanblueaway';meanredtoward';meanredaway'];
g1 = repmat({'Strain 1 towards' },length(meanbluetoward),1);
g2 = repmat({'Strain 1 away'},length(meanblueaway),1);
g3 = repmat({'Strain 2 towards'},length(meanredtoward),1);
g4 = repmat({'Strain 2 away'},length(meanredaway),1);
g=[g1;g2;g3;g4];
boxplot(data,g)
% xlabel(['Strain 1 towards Strain 2','Strain 1 away from Strain 2',' Strain 2 towards Strain 1  ',' Strain 2 away from Strain 1 '])
ylabel('Mean number of cell-blocks');
title('Metabolite-Mediated Model');
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
set(gca, 'Xtick', 1:4);
xt = get(gca, 'XTick');
[h1,p1] = ttest(meanbluetoward,meanblueaway)
[h2,p2] = ttest(meanredtoward,meanredaway)
