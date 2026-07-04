function fontset(fig,size)
if nargin==0
    fig=gcf;
    size=24;
end
set(findall(fig,'-property','MarkerSize'),'MarkerSize',16)
set(findall(fig,'-property','FontSize'),'FontSize',size)
set(findall(fig,'-property','LineWidth'),'LineWidth',1.6)
set(fig,'Position',[100 100 600 400]);
set(findall(fig,'-property','ColorOrder'),'ColorOrder',colororder(["blue";"red";"black";"#77AC30";"#EDB120";"#7E2F8E"]))
legend('Location','Best')
end