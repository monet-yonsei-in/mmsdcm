function [pE,pC] = med_remove_prior(pE,pC,pF)
if ~isempty(pF)
    fn_pE = fieldnames(pE);
    fn_pF = fieldnames(pF);
    for nField_pE = 1:numel(fn_pE)
        for nField_pF = 1:numel(fn_pF)
            if strcmp(fn_pE{nField_pE},fn_pF{nField_pF})
                pE = rmfield(pE,fn_pF{nField_pF});
                pC = rmfield(pC,fn_pF{nField_pF});
            end
        end
    end
end
return;
