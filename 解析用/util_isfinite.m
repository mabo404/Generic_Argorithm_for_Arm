function tf = util_isfinite(x, type)
if nargin<2, type='auto'; end
xv = x(:);
switch lower(type)
    case 'vec3', tf = (numel(xv)==3) && all(isfinite(xv));
    case 'quat', tf = (numel(xv)==4) && all(isfinite(xv));
    case {'mat','matrix'}, tf = all(isfinite(xv));
    otherwise, tf = all(isfinite(xv));
end
end
