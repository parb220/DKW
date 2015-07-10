function pall=pallGet(pvar,pvec,pidx);
pall=pvec;
for i=1:length(pvar)
     pall(pidx(i))=pvar(i);   % revised pall 
end
