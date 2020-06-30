clear; clc;
% 1 red 2 blue
AllNcells1=[];
AllNcells2=[];

for zzz=1:10
    niter = 100;
    gridlenx = 100;
    gridleny = 100;
    % GRID=zeros(100,niter,gridlenx,gridleny);
    Ncells1=[];
    Ncells2=[];
    
    for zz=1:100
        QS_AMT_1 = 1.6*10^(-1);
        QS_AMT_2 = 1.6*10^(-1);
        D = 8*10^(-10);
        S1 = zzz*0.01;
        S2 = 0.05;
        ThresholdforU2 = 1;
        ThresholdforU1 = 1;
        dt = 600;
        dl = 3*10^(-3);
        divisions=zeros(niter,gridlenx,gridleny);
        U1 = zeros(niter,gridlenx,gridleny);
        U2 = zeros(niter,gridlenx,gridleny);
        keepcells=zeros(niter,gridlenx,gridleny);
        
        cells=cellclass.empty;
        cellsoftype1=4;
        cellsoftype2=4;
        rows1=randsample(40:60,cellsoftype1+cellsoftype2);
        cols1=randsample(40:60,cellsoftype1+cellsoftype2);
        for i=1:cellsoftype1
            cells(end+1)=cellclass(1,rows1(i),cols1(i));
        end
        for j=cellsoftype1+1:cellsoftype1+cellsoftype2
            cells(end+1)=cellclass(2,rows1(j),cols1(j));
        end
        % c1=cellclass(1,45,45);
        % c2=cellclass(1,50,50);
        % c3=cellclass(2,45,50);
        % c4=cellclass(2,50,45);
        % cells=[c1,c2,c3,c4];
        
        
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
            dgAHL = 0.0046;
            for row=1:gridlenx
                for col=1:gridleny
                    U1(time,row,col) = U1(time,row,col) - dgAHL*U1(time,row,col);
                    U2(time,row,col) = U2(time,row,col) - dgAHL*U2(time,row,col);
                end
            end
            
            %divisions
            for c=1:length(cells)
                row=cells(c).locationx; col=cells(c).locationy;
                newfitness=cells(c).fitness;
                if cells(c).ID==1
                    if U2(time,row,col)>ThresholdforU2
                        newfitness = cells(c).fitness + S1;
                    end
                else
                    if U1(time,row,col)>ThresholdforU1
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
        end
        disp(['Iteration number ',num2str(zz)]);
        type1=0;
        type2=0;
        for c=1:length(cells)
            if cells(c).ID==1
                type1=type1+1;
            else
                type2=type2+1;
            end
        end
        Ncells1(end+1)=type1;
        Ncells2(end+1)=type2;
    end
    AllNcells1=[AllNcells1;Ncells1];
    AllNcells2=[AllNcells2;Ncells2];

end

%%
% figure
% for i=1:niter
%     drawnow
%     pcolor(reshape(U1(i,:,:),[gridlenx,gridleny]))
% end
%%
% end
% map=[0,0,0; 1,0,0; 0,0,1];
% figure
% for i=1:niter
%     drawnow
%     pcolor(reshape(keepcells(i,:,:),[gridlenx,gridleny]))
%     colormap(map)
%     title(i)
% end

% figure()
% data=[Ncells1';Ncells2'];
% g1 = repmat({'Red'},length(Ncells1),1);
% g2 = repmat({'Blue'},length(Ncells2),1);
% g=[g1;g2];
% boxplot(data,g)
save('somedata2.mat')

%%
load somedata2
x = 0.01:0.01:0.1;
y1 = mean(AllNcells1');
err1 = std(AllNcells1');
y2 = mean(AllNcells2');
err2 = std(AllNcells2');

figure
xlim([0.005, 0.105]);
c = [0 0 1];
d = [1 0 0];
sz = 2;
for i=1:length(x)
scatter(repmat(x(i)-0.0005,[1,100]),AllNcells1(i,:),sz,d,'filled');
hold on
scatter(repmat(x(i),[1,100]),AllNcells2(i,:),sz,c,'filled');
end

errorbar(x-0.0005,y1,err1,'--xr','LineWidth',1.2);
hold on
errorbar(x,y2,err2,'--xb','LineWidth',1.2);
legend("Strain 1","Strain 2");
xlabel("Gain of fitness of strain 1");
ylabel("Number of cell-blocks");
title("Two strains");
xlim([0.005, 0.105]);

