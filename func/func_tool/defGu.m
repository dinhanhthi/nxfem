%% Function handle of g(u)
% List of all function handles of g(u) used in getGM* and getL*

function gu = defGu
    gu.change = @defChange; % changeable
    gu.dchange = @defDiffChange; % derivative (depends on chage)
    gu.G1 = @defG1; % use in main_sys_linda
end


%% Definition

function val = defChange(u) % changeable g(u)
    val = u.^2;
end

function val = defDiffChange(u) % derivative of g(u) changeable
    val = 2*u;
end

function val = defG1(u) % *ug'(u)-g(u)
    gu = defChange(u);
    dgu = defDiffChange(u);
    val = u.*dgu - gu;
end