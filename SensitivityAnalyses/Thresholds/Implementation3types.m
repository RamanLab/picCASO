clear; clc;
%red 1, green 2, blue 3
% QS_AMT_2 -> U2 -> more red cells
% QS_AMT_3 -> U3 -> more green cells
% QS_AMT_1 -> U1 -> more blue cells
AllNcells1=[];
AllNcells2=[];
AllNcells3 =[];

for zzz=1:10
    gridlenx = 100;
    gridleny = 100;
    niter = 100;
    Ncells1=[];
    Ncells2=[];
    Ncells3=[];
    
    for zz=1:100
        QS_AMT_1 = 1.6*10^(-1);
        QS_AMT_2 = 1.6*10^(-1);
        QS_AMT_3 = 1.6*10^(-1);
        D = 8*10^(-10);
        S1 = 0.08;
        S2 = 0.08;
        S3 = 0.08;
        Thresholds=linspace(0.1,10,10);
        Threshold1 = Thresholds(zzz);
        Threshold2 = 5;
        Threshold3 = 5;
        dt = 600;
        dl = 3*10^(-3);
        keepcells = zeros(niter,gridlenx,gridleny);
        U1 = zeros(niter,gridlenx,gridleny);
        U2 = zeros(niter,gridlenx,gridleny);
        U3 = zeros(niter,gridlenx,gridleny);
        
        % c1=cellclass(2,50,45);
        % c2=cellclass(2,50,50);
        % c6=cellclass(2,50,55);
        % c3=cellclass(1,45,50);
        % c4=cellclass(1,45,55);
        % c5=cellclass(1,45,45);
        % cells=[c1,c2,c3,c4,c5,c6];
        
        %     %red
        %     c1=cellclass(1,45,45);
        %     c2=cellclass(1,50,50);
        %     %green
        %     c3=cellclass(2,45,50);
        %     c4=cellclass(2,50,45);
        %     %blue
        %     c5=cellclass(3,47,47);
        %
        %     cells=[c1,c2,c3,c4,c5];
        %     keepcells=zeros(niter,gridlenx,gridleny);
        cells=cellclass.empty;
        cellsoftype1=4;
        cellsoftype2=4;
        cellsoftype3=4;
        rows1=randsample(40:60,cellsoftype1+cellsoftype2+cellsoftype3);
        cols1=randsample(40:60,cellsoftype1+cellsoftype2+cellsoftype3);
        for i=1:cellsoftype1
            cells(end+1)=cellclass(1,rows1(i),cols1(i));
        end
        for j=cellsoftype1+1:cellsoftype1+cellsoftype2
            cells(end+1)=cellclass(2,rows1(j),cols1(j));
        end
        for k=cellsoftype1+cellsoftype2+1:cellsoftype1+cellsoftype2+cellsoftype3
            cells(end+1)=cellclass(3,rows1(k),cols1(k));
        end
        
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
                if cells(c).ID==3
                    U3(time-1,row,col)=U3(time-1,row,col)+QS_AMT_3;
                end
            end
            
            %diffusion
            for row=1:gridlenx
                for col=1:gridleny
                    if row-1>0 && col-1>0 && row+1<gridlenx && col+1<gridleny
                        U1(time,row,col) = U1(time-1,row,col) + D*dt/(dl^2)*(U1(time-1,row-1,col) + U1(time-1,row+1,col) + U1(time-1,row,col-1) + U1(time-1,row,col+1) - 4*U1(time-1,row,col));
                        U2(time,row,col) = U2(time-1,row,col) + D*dt/(dl^2)*(U2(time-1,row-1,col) + U2(time-1,row+1,col) + U2(time-1,row,col-1) + U2(time-1,row,col+1) - 4*U2(time-1,row,col));
                        U3(time,row,col) = U3(time-1,row,col) + D*dt/(dl^2)*(U3(time-1,row-1,col) + U3(time-1,row+1,col) + U3(time-1,row,col-1) + U3(time-1,row,col+1) - 4*U3(time-1,row,col));
                    end
                end
            end
            
            %AHL Degradation
            dgAHL = 0.0046;
            for row=1:gridlenx
                for col=1:gridleny
                    U1(time,row,col) = U1(time,row,col) - dgAHL*U1(time,row,col);
                    U2(time,row,col) = U2(time,row,col) - dgAHL*U2(time,row,col);
                    U3(time,row,col) = U3(time,row,col) - dgAHL*U3(time,row,col);
                end
            end
            
            %divisions
            for c=1:length(cells)
                row=cells(c).locationx; col=cells(c).locationy;
                newfitness=cells(c).fitness;
                if cells(c).ID==1
                    if U2(time,row,col)>Threshold2
                        newfitness = cells(c).fitness + S1;
                    end
                end
                if cells(c).ID==2
                    if U3(time,row,col)>Threshold3
                        newfitness = cells(c).fitness + S2;
                    end
                end
                if cells(c).ID==3
                    if U1(time,row,col)>Threshold1
                        newfitness = cells(c).fitness + S3;
                    end
                end
                newcell=celldivision(cells(c),cells,newfitness);
                if isa(newcell,'cellclass')
                    cells(end+1)=newcell;
                end
            end
            
            
            for c=1:length(cells)
                if cells(c).ID==1
                    keepcells(time,cells(c).locationx,cells(c).locationy)=1;
                end
                if cells(c).ID==2
                    keepcells(time,cells(c).locationx,cells(c).locationy)=2;
                end
                if cells(c).ID==3
                    keepcells(time,cells(c).locationx,cells(c).locationy)=3;
                end
            end
        end
        disp(['Iteration number ',num2str(zz)]);
        type1=0;
        type2=0;
        type3=0;
        for c=1:length(cells)
            if cells(c).ID==1
                type1=type1+1;
            elseif cells(c).ID==2
                type2=type2+1;
            elseif cells(c).ID==3
                type3=type3+1;
            end
        end
        Ncells1(end+1)=type1;
        Ncells2(end+1)=type2;
        Ncells3(end+1)=type3;
    end
    AllNcells1=[AllNcells1;Ncells1];
    AllNcells2=[AllNcells2;Ncells2];
    ALlNcells3=[AllNcells3;Ncells3];
end

% figure
% for i=1:niter
%     drawnow
%     pcolor(reshape(U3(i,:,:),[gridlenx,gridleny]))
% end
%%
% map=[0,0,0; 1,0,0; 0,0,1; 0 1 0];
% figure
% for i=1:niter
%     drawnow
%     pcolor(reshape(keepcells(i,:,:),[gridlenx,gridleny]))
%     colormap(map)
%     title(i)
% end

%%
% figure()
% data=[Ncells1';Ncells2';Ncells3'];
% g1 = repmat({'Red'},length(Ncells1),1);
% g2 = repmat({'Green'},length(Ncells2),1);
% g3 = repmat({'Blue'},length(Ncells3),1);
% g=[g1;g2;g3];
% boxplot(data,g)
save('somedata3.mat')

%%
load somedata3

y1 = mean(AllNcells1');
err1 = std(AllNcells1');
y2 = mean(AllNcells2');
err2 = std(AllNcells2');
y3 = mean(AllNcells3');
err3 = std(AllNcells3');

figure
b = [0 1 0];
c = [0 0 1];
d = [1 0 0];
sz = 2;
for i=1:length(Thresholds)
scatter(repmat(Thresholds(i)-0.05,[1,100]),AllNcells1(i,:),sz,d,'filled');
hold on
scatter(repmat(Thresholds(i),[1,100]),AllNcells2(i,:),sz,c,'filled');
scatter(repmat(Thresholds(i)+0.05,[1,100]),AllNcells3(i,:),sz,b,'filled');
end

errorbar(Thresholds-0.05,y1,err1,'--xr','LineWidth',1.2);
hold on
errorbar(Thresholds,y2,err2,'--xb','LineWidth',1.2);
errorbar(Thresholds+0.05,y3,err3,'--xg','LineWidth',1.2);
legend("Strain 1","Strain 2","Strain 3");
xlabel("Threshold of QS molecule amount for gene expression in strain 1 (nmol/min)");
ylabel("Number of cell-blocks");
title("Three strains");
xlim([0 10.1]);


