function ll = vline(XX,varargin)
YY=get(gca,'ylim');
ll=line([XX XX],[YY(1) YY(2)]);
ll.LineStyle = '--';
set(ll,'color','k')

if length(varargin) > 0
  set(ll, varargin{:});
end
