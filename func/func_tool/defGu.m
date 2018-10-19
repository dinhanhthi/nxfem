%% Function handle of g(u)
% List of all function handles of g(u) used in getGM* and getL*

function gu = defGu
    gu.change = @defChange; % changeable
    gu.dchange = @defDiffChange; % derivative (depends on chage)
    gu.G1 = @defG1; % used in main_sys_linda
end


%% Definition

function val = defChange(u,pa) % changeable g(u)
%     val = u.^2;
    val = u./(pa.K0 + u); % chopp06combine
end

function val = defDiffChange(u,pa) % derivative of g(u) changeable
%     val = 2*u;
    val = pa.K0/(pa.K0+u).^2; % chopp06combine
end

function val = defG1(u) % *ug'(u)-g(u)
    gu = defChange(u);
    dgu = defDiffChange(u);
    val = u.*dgu - gu;
end