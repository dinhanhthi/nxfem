function vtest = test_inputParser(a,b,c,varargin)
% Wanna know variable input and reuired input are the same or not?
% Result: 
%   - must put required input in both line of function and "addRequired"
%   - must use p.Results for those created by addParameter
%   - can use both p.Results.a and a if a is required
%   - a in input and 'a' in addRequired are the same


p = inputParser;
addRequired(p,'a');
addRequired(p,'b');
addRequired(p,'c');
addParameter(p,'d',1);
parse(p,a,b,c,varargin{:});

vtest = a+2*b+3*c-100*p.Results.d;
% vtest = p.Results.a - a;

end