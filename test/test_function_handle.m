defF = @findDefFtw;
defF1 = @(x,y,pa) defF(x,y,pa,1);
defF2 = @(x,y,pa) defF(x,y,pa,2);

f1 = defF1(1,2,3)
f2 = defF2(1,2,3)

function valF = findDefFtw(xx,yy,pa,sub)
    if sub==1
        valF=1;
    else
        valF=2;
    end
end

