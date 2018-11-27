%% Function handle of g(u)
% List of all function handles of g(u) used in getGM* and getL*
% NOTE: if change gu.change, you must change gu.dchange and gu.changeFu!!!

function gu = defGu
    gu.change = @defChange; % changeable
    gu.dchange = @defDiffChange; % derivative (depends on change)
    gu.changeFu = @defChangeFu; % derivative (depends on change)
    gu.G1 = @defG1; % used in main_sys_linda
end


%% Definition

function val = defChange(u,pa) % changeable g(u)
%     val = u.^2; % main_sys_linda, main_chap4
    val = u./(pa.K0 + u); % chopp06combine
end

function val = defDiffChange(u,pa) % derivative of g(u) changeable
%     val = 2*u; % main_sys_linda, main_chap4
    val = pa.K0/(pa.K0+u).^2; % chopp06combine
end

function val = defChangeFu(u,pa)
% Only used for finding Fdel = Au(uold)*uold
% We use 1/(K0+uold) instead of uold/(K0+uold)
% defChange = defChangeFu * u
val = 1/(pa.K0 + u); % chopp06combine
end

function val = defG1(u,pa) % *ug'(u)-g(u)
    % main_sys_linda
    gu = defChange(u,pa);
    dgu = defDiffChange(u,pa);
    val = u.*dgu - gu;
end