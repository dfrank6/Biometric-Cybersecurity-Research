function istrash = trashcheck(ERP)
istrash = true;
for i = 1:length(ERP)
    if ERP(i) > 400
        istrash = false;
        break;
    end
end
end