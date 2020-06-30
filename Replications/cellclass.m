classdef cellclass
    
    properties
        ID
        locationx
        locationy
        fitness =  1/12;
    end
    
    methods
        
        function c = cellclass(ID,locationx,locationy)
            c.ID=ID;
            c.locationx=locationx;
            c.locationy=locationy;
        end
        
        function newcell = celldivision(c,cells,newfitness)
            p=rand;
            spaces=emptyneighbours(c,cells);
            [options,~]=size(spaces);
            if newfitness>p && ~isempty(spaces)
                newcelllocation=spaces(randi(options),:);
                newcell=cellclass(c.ID,newcelllocation(1),newcelllocation(2));
            else
                newcell=0;
            end
        end
        function en=emptyneighbours(c,cells)
            en=[c.locationx+1,c.locationy;c.locationx,c.locationy-1;c.locationx-1,c.locationy;c.locationx,c.locationy+1];
            todelete=[];
            for ii = 1:length(cells)
                for jj=1:length(en)
                    if cells(ii).locationx==en(jj,1) && cells(ii).locationy==en(jj,2)
                        todelete(end+1)=jj;
                    end
                end
            end
            en(todelete,:)=[];
        end
        
    end
end