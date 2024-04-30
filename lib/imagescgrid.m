function [hi,hl,hlr,hld]=imagescgrid(hi,colwidth_grid,colwidth_rect,colwidth_diag)
if nargin<2, colwidth_grid=[]; end
if nargin<3, colwidth_rect=[]; end
if nargin<4, colwidth_diag=[]; end

if isempty(colwidth_grid), colwidth_grid=[1 1 1]*0.3; end
if isempty(colwidth_rect), colwidth_rect=[1 1 1]*0; end
if isempty(colwidth_diag), colwidth_diag=colwidth_grid; end

if isnumeric(hi) 
    hi=imagesc(hi);
    axis image;axis off; 
end
ncg=size(get(hi,'CData'),1);
hold on;
c1=1; hl=[]; hlr=[];
for i = 1:ncg
    hl(c1)=plot([.5,ncg+0.5],[i-.5,i-.5],'k-');c1=c1+1;
    hl(c1)=plot([i-.5,i-.5],[.5,ncg+0.5],'k-');c1=c1+1;
end
set(hl,'color',colwidth_grid(1:3));
if length(colwidth_grid)>3, set(hlr,'linewidth',colwidth_grid(4)); end

nc=ncg;
xs=[0.5 0.5 nc+0.5 0.5;
    nc+0.5 0.5 nc+0.5 nc+0.5;
    nc+0.5 nc+0.5 0.5 nc+0.5;
    0.5 nc+0.5 0.5 0.5]; c1=1;
for i=1:4 
    hlr(c1)=plot(xs(i,1:2),xs(i,3:4),'k-');c1=c1+1;
end
set(hlr,'color',colwidth_rect(1:3));
if length(colwidth_rect)>3, set(hlr,'linewidth',colwidth_rect(4)); end

hld=[];
if ~all(colwidth_grid(1:3)==colwidth_rect(1:3))
    c1=1;nc=1;       
    xs=[0.5 0.5 nc+0.5 0.5;
    nc+0.5 0.5 nc+0.5 nc+0.5;
    nc+0.5 nc+0.5 0.5 nc+0.5;
    0.5 nc+0.5 0.5 0.5]; 
    for j = 1:ncg        
        for i=1:4 
            hlr(c1)=plot(xs(i,1:2)+j-1,xs(i,3:4)+j-1,'k-');c1=c1+1;
        end
    end    
    set(hlr,'color',colwidth_diag(1:3));
    if length(colwidth_diag)>3, set(hlr,'linewidth',colwidth_diag(4)); end
end

if 0
    set(gca,'LineWidth',2,'XColor',[0 0 0],'YColor',[0 0 0]);         
    set(gca,'XTickLabel',[],'YTickLabel',[]); axis on;                
end

hold off;
