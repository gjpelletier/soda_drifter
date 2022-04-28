function [ck] = jet2(nk)

% Matlab function to obtain nk rows of the jet colormap 
%
% Written by:
%
% Greg Pelletier 
% gregp@sccwrp.org
%

% r,g,b=zeros(256),zeros(256),zeros(256)
% for i = 0:255
  % n=4*i/256
  % r[i+1]=255*min(max(min(n-1.5,-n+4.5),0),1);
  % g[i+1]=255*min(max(min(n-0.5,-n+3.5),0),1);
  % b[i+1]=255*min(max(min(n+0.5,-n+2.5),0),1);
% end

inc=256/(nk-1);
for ik = 1:nk
	if ik==1
		i=1;
	elseif ik==2
		i=inc;
	else
		i=i+inc;
	end
	if i>256; i=256; end
	n=4*(i-1)/256;
    ck(ik,1)=round(255*min(max(min(n-1.5,-n+4.5),0),1))/255;
    ck(ik,2)=round(255*min(max(min(n-0.5,-n+3.5),0),1))/255;
    ck(ik,3)=round(255*min(max(min(n+0.5,-n+2.5),0),1))/255;
end		% for ik=1:nk

end    % function