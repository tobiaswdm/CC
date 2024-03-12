function DrawDispersion(sys,color,savepath)

k_continuous =linspace(0,floor(sys.N_s/2),100*sys.N_s);
[r_continuous,~] = DispersionRelation(k_continuous,sys);

figure;
plot(k_continuous,r_continuous,'LineWidth',1.5,'Color',color.reference)
hold on;
plot(sys.k,sys.r_k,'o','Color',color.reference,'MarkerSize',8,'LineWidth',1.5)
box on;
title('Dispersion Diagram')
xlabel('Wave number - $k$')
ylabel('Eigenfrequency ratio - $r_k$')
xlim([0 floor(sys.N_s/2)])
ylim([0.99 1.01*sys.r_k(end)])
if ~isempty(savepath)
    savefig([savepath 'dispersion_diagram.fig'])
end


end