% % extract best-point and ranges
function boption_new = extract_parm(results,boption,cntnum,margin)
if nargin < 4, margin=0.01; end
if nargin < 3, cntnum=1; end
%cntnum=1; margin=0.01;
% update range for parameter serarch
boption_new = get_ranges(results, cntnum, margin,boption);
end

function boption_new = get_ranges(results, cutnum, margin,boption)

[val,tmpid]=sort(results.ObjectiveTrace);
%subplot(1,2,1); plot(val,'LineWidth',2); ylabel('-F'); title('all'); 
%set(gca,'FontSize',24,'LineWidth',1.5,'FontName','arial narrow');
%subplot(1,2,2); plot(val(1:cutnum),'LineWidth',2); ylabel('-F'); title('best10')
%set(gca,'FontSize',24,'LineWidth',1.5,'FontName','arial narrow');


tmp2 = results.XTrace(tmpid(1:cutnum),:);
parm_names = fieldnames(tmp2);

for ii=1:(length(parm_names)-3)
    tmpchk = table2array(tmp2(:,ii));
    min1=min(tmpchk);
    max1=max(tmpchk);
    
    if( min1 > tmpchk(1) - margin ), min1 = tmpchk(1) - margin; end
    if( max1 < tmpchk(1) + margin ), max1 = tmpchk(1) + margin; end

    cmd=sprintf('name= %s, min= %f max= %f best= %f', parm_names{ii}, min1, max1,tmpchk(1));
    disp(cmd)
    try
        boption_new(ii).pEname =  boption(ii).pEname;
    catch
    end    
    boption_new(ii).name =  parm_names{ii};
    boption_new(ii).range = [min1 max1];
    boption_new(ii).estBest = results.bestPoint.(parm_names{ii});
    
end
end