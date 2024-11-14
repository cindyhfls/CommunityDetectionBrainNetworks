function ll = vline(XX,varargin)
%In the newer MATLAB versions you can use vline(y0).
YY=get(gca,'ylim');
ll=line([XX XX],[YY(1) YY(2)]);
ll.LineStyle = '--';
set(ll,'color','k')

if length(varargin) > 0
  set(ll, varargin{:});
end
