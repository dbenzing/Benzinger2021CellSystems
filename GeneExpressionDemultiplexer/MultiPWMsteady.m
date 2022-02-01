function [RFP, YFP] = MultiPWMsteady(p,p2,TFtot,TFtot2,Reptot,I, period, width,kdegProt,tspan,initial)

	pulsenumber = tspan(end) / period;

	for numb = 1 : pulsenumber
		tspan = [0 width];
		[T,Y] = ode23s(@(t,y) detDemulti(t,y,p,p2,TFtot,TFtot2,Reptot,I,kdegProt), tspan, initial);
		initial = Y(end,:);
			
		tspan = [0 period-width];
		[T,Y] = ode23s(@(t,y) detDemulti(t,y,p,p2,TFtot,TFtot2,Reptot,0,kdegProt), tspan, initial);
		initial = Y(end,:);
    end
	
    RFP = Y(end,4);
	YFP = Y(end,7);
	
end


