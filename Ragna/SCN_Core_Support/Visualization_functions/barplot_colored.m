function [h,s]=barplot_colored(data,varargin)

% haven't had time to document this yet
% 
% [h,s]=barplot_colored(data,varargin)
% 
% h is the axes handle, s is the handle to its children
% 
%     if strcmp(varargin{k},'colormap')
%         eval(['colormapfun=@' varargin{k+1} ';'])
%     elseif strcmp(varargin{k},'title')
%         title=varargin{k+1};
%     elseif strcmp(varargin{k},'XTickLabels')
%         XTickLabels=varargin{k+1};
%     elseif strcmp(varargin{k},'ylabel')
%         ylabel=varargin{k+1};
%     elseif strcmp(varargin{k},'xlabel')
%         xlabel=varargin{k+1};
%     end

if iscell(data)
    for k=1:length(data)
        means(k)=mean(data{k});
        stderr(k)=std(data{k})/sqrt(length(data{k}));
    end
else
    means=mean(data);
    stderr=std(data)/sqrt(size(data,1));
end

for k=1:2:length(varargin)
    if strcmp(varargin{k},'colormap')
        eval(['colormapfun=@' varargin{k+1} ';'])
    elseif strcmp(varargin{k},'title')
        title=varargin{k+1};
    elseif strcmp(varargin{k},'XTickLabel')
        XTickLabel=varargin{k+1};
    elseif strcmp(varargin{k},'Ylabel')
        Ylabel=varargin{k+1};
    elseif strcmp(varargin{k},'Xlabel')
        Xlabel=varargin{k+1};
    end
end


if ~exist('colormapfun','var')
    colormapfun=@hsv;
end

h=bar(means);
s=get(h,'Children');

colormap(colormapfun(length(means)));
set(s,'CData',1:length(means));

hold on

errorbar(means,stderr,'k','LineWidth',2,'LineStyle','none')
set(gca,'Xlim',[0 length(means)+1])
if exist('XTickLabel','var')
    set(gca,'XTickLabel',XTickLabel)
else
    set(gca,'XTickLabel',[])
end
if exist('title','var')
    title(title)
end
if exist('Ylabel','var')
    ylabel(Ylabel)
end
if exist('Xlabel','var')
    xlabel(Xlabel)
end

hold off