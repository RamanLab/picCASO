clear; clc;

%% parameters
gridlenx = 100;
gridleny = 100;
QS_AMT_1 = 1.6*10^(-1);
QS_AMT_2 = 1.6*10^(-1);
D = 8*10^(-10);
dgAHL = 0.0046;

%The fitness values are selected after FBA-based gene selection is
%performed
if isfile('Fitnesses.mat')
    load Fitnesses.mat;
else
    BaseFitness1 = 0.082;
    FitnessGain1 = 0.082;
    BaseFitness2 = 0.082;
    FitnessGain2 = 0.082;
end
% The threshold values are set based on the device-AHL pairs chosen
if isfile('Thresholds.mat')
    load Thresholds.mat;
    ThresholdforU2 = Threshold1*10^9;
    ThresholdforU1 = Threshold2*10^9;
else
    ThresholdforU1 = 1;
    ThresholdforU2 = 1;
end
dt = 600;
dl = 3*10^(-3);
niter = 100;
divisions=zeros(niter,gridlenx,gridleny);
GRID = zeros(niter,gridlenx,gridleny);
U1 = zeros(niter,gridlenx,gridleny);
U2 = zeros(niter,gridlenx,gridleny);
keepcells=zeros(niter,gridlenx,gridleny);

%% seeding cells
cells=cellclass.empty;
cells(end+1)=cellclass(1,50,50);
cells(end+1)=cellclass(1,55,55);
cells(end+1)=cellclass(2,50,55);
cells(end+1)=cellclass(2,55,50);

%% model run
for time=2:niter
    
    
    %Inject molecules from cell
    for c=1:length(cells)
        row=cells(c).locationx; col=cells(c).locationy;
        if cells(c).ID==1
            U1(time-1,row,col)=U1(time-1,row,col)+QS_AMT_1;
        end
        if cells(c).ID==2
            U2(time-1,row,col)=U2(time-1,row,col)+QS_AMT_2;
        end
    end
    
    %diffusion
    for row=1:gridlenx
        for col=1:gridleny
            if row-1>0 && col-1>0 && row+1<gridlenx && col+1<gridleny
                U1(time,row,col) = U1(time-1,row,col) + D*dt/(dl^2)*(U1(time-1,row-1,col) + U1(time-1,row+1,col) + U1(time-1,row,col-1) + U1(time-1,row,col+1) - 4*U1(time-1,row,col));
                U2(time,row,col) = U2(time-1,row,col) + D*dt/(dl^2)*(U2(time-1,row-1,col) + U2(time-1,row+1,col) + U2(time-1,row,col-1) + U2(time-1,row,col+1) - 4*U2(time-1,row,col));
            end
        end
    end
    
    %AHL Degradation
    for row=1:gridlenx
        for col=1:gridleny
            U1(time,row,col) = U1(time,row,col) - dgAHL*U1(time,row,col);
            U2(time,row,col) = U2(time,row,col) - dgAHL*U2(time,row,col);
        end
    end
    
    
    %divisions
    for c=1:length(cells)
        row=cells(c).locationx; col=cells(c).locationy;
        if cells(c).ID==1
            newfitness = BaseFitness1;
            if U2(time,row,col)>ThresholdforU2
                newfitness = newfitness + FitnessGain1;
            end
        else
            newfitness = BaseFitness2;
            if U1(time,row,col)>ThresholdforU1
                newfitness = newfitness + FitnessGain2;
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
    disp(time)
end

%% plot generation
% map=[0,0,0; 1,0,0; 0,0,1];
% figure
% for i=1:niter
%     drawnow
%     pcolor(reshape(keepcells(i,:,:),[gridlenx,gridleny]))
%     colormap(map)
%     title(['Number of cells = ',num2str(nnz(reshape(keepcells(i,:,:),[gridlenx,gridleny])))])
% end

%% video generation
v = VideoWriter('Square.avi');
open(v)
map=[0,0,0; 1,0,0; 0,0,1];
figure
for i=1:niter
    pcolor(reshape(keepcells(i,:,:),[gridlenx,gridleny]))
    colormap(map)
    title(['Iteration = ',num2str(i),' | No. of cell blocks = ',num2str(nnz(reshape(keepcells(i,:,:),[gridlenx,gridleny])))])
    frame=getframe(gcf);
    writeVideo(v,frame);
end
close(v)