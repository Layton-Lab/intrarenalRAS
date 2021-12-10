function save_pdf(title)

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
saveas(gcf,strcat(title,'.pdf'));
%print -dpdf -painters nonAngII_variable_fig

end