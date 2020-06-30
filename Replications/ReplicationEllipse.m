clear; clc;
% 1 red 2 blue
load gridellipse.mat
gridlenx = 100;
gridleny = 100;
QS_AMT_1 = 3*10^(-1);
% QS_AMT_2 = 1.6*10^(-1);
D = 8*10^(-10);
% S1 = 0.08;
% S2 = 0.08;
% Threshold = 1;
dt = 600;
dl = 3*10^(-3);
niter = 450;
dgAHL = 0.0046;
divisions=zeros(niter,gridlenx,gridleny);
GRID = zeros(niter,gridlenx,gridleny);
U1 = zeros(niter,gridlenx,gridleny);
U2 = zeros(niter,gridlenx,gridleny);
reds=zeros(niter,gridlenx,gridleny);

% cells=cellclass.empty;
% c1=cellclass(1,50,45);
% c2=cellclass(1,50,55);
% cells=[c1,c2];
% xs=[];
% ys=[];
% for c=1:length(cells)
%     xs(end+1)=cells(c).locationx;
%     ys(end+1)=cells(c).locationy;
% end
% for i=10:90
%     for j=10:90
%         skip=0;
%         for z=1:length(xs)
%             if i==xs(z) && j==ys(z)
%                 skip=1;
%             end
%         end
%         if skip==1
%             continue
%         end
%         cells=[cells,cellclass(2,i,j)];
%     end
%     disp(['Populating row ',num2str(i)])
% end
% keepcells=zeros(niter,gridlenx,gridleny);
% reds=zeros(niter,gridlenx,gridleny);
% disp('The grid has been populated')

for time=2:niter
    
    
    %Inject molecules from cell
    for c=1:length(cells)
        row=cells(c).locationx; col=cells(c).locationy;
        if cells(c).ID==1
            U1(time-1,row,col)=U1(time-1,row,col)+QS_AMT_1;
        end
%         if cells(c).ID==2
%             U2(time-1,row,col)=U2(time-1,row,col)+QS_AMT_2;
%         end
    end
    
    %diffusion
    for row=1:gridlenx
        for col=1:gridleny
            if row-1>0 && col-1>0 && row+1<gridlenx && col+1<gridleny
                U1(time,row,col) = U1(time-1,row,col) + D*dt/(dl^2)*(U1(time-1,row-1,col) + U1(time-1,row+1,col) + U1(time-1,row,col-1) + U1(time-1,row,col+1) - 4*U1(time-1,row,col));
%                 U2(time,row,col) = U2(time-1,row,col) + D*dt/(dl^2)*(U2(time-1,row-1,col) + U2(time-1,row+1,col) + U2(time-1,row,col-1) + U2(time-1,row,col+1) - 4*U2(time-1,row,col));
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
    
    %fluorescence
    for c=1:length(cells)
        row=cells(c).locationx; col=cells(c).locationy;
        if cells(c).ID==1
            reds(time,cells(c).locationx,cells(c).locationy)=1;
        else
            if U1(time,row,col)>10^(-2) && U1(time,row,col)<2*10^(-1)
                reds(time,cells(c).locationx,cells(c).locationy)=2;
            end
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
    disp(['Iteration ',num2str(time)])
end



% figure
% for i=1:niter
%     drawnow
%     pcolor(reshape(U1(i,:,:),[gridlenx,gridleny]))
% end

% map=[0,0,0; 1,0,0; 0,0,1];
% figure
% for i=1:niter
%     drawnow
%     pcolor(reshape(keepcells(i,:,:),[gridlenx,gridleny]))
%     colormap(map)
%     title(i)
% end

map=[0,0,0; 1,0,0;0 1 0];
figure
for i=1:niter
    drawnow
    pcolor(reshape(reds(i,:,:),[gridlenx,gridleny]))
    colormap(map)
    title(i)
end


